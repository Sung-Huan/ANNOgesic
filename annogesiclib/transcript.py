import os
import sys
import shutil
from subprocess import call
from annogesiclib.helper import Helper
from annogesiclib.multiparser import Multiparser
from annogesiclib.converter import Converter
from annogesiclib.combine_frag_tex import combine
from annogesiclib.stat_TA_comparison import stat_ta_tss, stat_ta_gff
from annogesiclib.transcript_assembly import assembly
from annogesiclib.fill_gap import fill_gap, longer_ta


class TranscriptAssembly(object):

    def __init__(self, out_folder):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.converter = Converter()
        self.gff_outfolder = os.path.join(out_folder, "gffs")
        self.tran_path = os.path.join(self.gff_outfolder, "tmp")
        self.stat_path = os.path.join(out_folder, "statistics")
        self.tmps = {"gff": "tmp.gff", "merge": "tmp_merge",
                     "tran": os.path.join(out_folder, "tmp_tran"),
                     "tss_ta": os.path.join(self.gff_outfolder, "tmp_tss_ta"),
                     "ta_tss": os.path.join(self.gff_outfolder, "tmp_ta_tss"),
                     "ta_gff": os.path.join(self.gff_outfolder, "tmp_ta_gff"),
                     "gff_ta": os.path.join(self.gff_outfolder, "tmp_gff_ta"),
                     "uni": os.path.join(self.gff_outfolder, "tmp_uni"),
                     "overlap": os.path.join(self.gff_outfolder, "tmp_overlap")}
        self.frag = "transcript_assembly_fragment.gff"
        self.tex = "transcript_assembly_tex.gff"
        self.endfix_tran = "transcript.gff"

    def _compute_transcript(self, wig_f, wig_r, height, width, wig_folder,
                            out_folder, tolerance, replicates, wig_type, strain,
                            libs, tex, low_cutoff):
        print("Computing transcript assembly for {0}...".format(strain))
        out = os.path.join(out_folder, "_".join([strain, wig_type]))
        assembly(wig_f, wig_r, height, width, tolerance, low_cutoff,
                 wig_folder, tex, libs, replicates, out)

    def _combine_wigs(self, wig_folder, merges, out_wig, direct, strain):
        for wig in os.listdir(wig_folder):
            if (direct in wig):
                filename = wig.split("_STRAIN_")
                if filename[1][:-4] == strain:
                    print("Merge {0} now....".format(wig))
                    self.helper.merge_file(os.path.join(wig_folder, wig),
                                           os.path.join(merges, out_wig))

    def _compute(self, merges, wig_folder, tex, height, width, tolerance,
                 replicates, out_folder, wig_type, libs, low_cutoff):
        strains = []
        wigs = os.path.join(wig_folder, "tmp")
        for wig in os.listdir(wigs):
            filename = wig.split("_STRAIN_")
            if filename[1][:-4] not in strains:
                strains.append(filename[1][:-4])
        for strain in strains:
            self._combine_wigs(wigs, merges, "_".join([strain, "forward.wig"]),
                               "forward", strain)
            print("Merge {0}_forward.wig complete...".format(strain))
            self._combine_wigs(wigs, merges, "_".join([strain, "reverse.wig"]),
                               "reverse", strain)
            print("Merge {0}_reverse.wig complete...".format(strain))
        for strain in strains:
            f_file = os.path.join(merges, "_".join([strain, "forward.wig"]))
            r_file = os.path.join(merges, "_".join([strain, "reverse.wig"]))
            self._compute_transcript(f_file, r_file, height, width, wig_folder,
                                     out_folder, tolerance, replicates,
                                     wig_type, strain, libs, tex, low_cutoff)
        return strains

    def _compare_tss(self, tas, gff_outfolder, tss_path, stat_path, fuzzy):
        self.multiparser.parser_gff(tss_path, "TSS")
        self.multiparser.combine_gff(gff_outfolder,
                          os.path.join(tss_path, "tmp"), "transcript", "TSS")
        print("Comaring of Transcript assembly and TSS file...")
        tss_folder = os.path.join(tss_path, "tmp")
        for ta in tas:
            ta_file = os.path.join(gff_outfolder,
                                   "_".join([ta, self.endfix_tran]))
            stat_tss_out = os.path.join(stat_path,
                           "".join(["stat_compare_Transcriptome_assembly_TSS_",
                                    ta, ".csv"]))
            for tss in os.listdir(tss_folder):
                filename = tss.split("_TSS")
                if (filename[0] == ta) and (tss.endswith(".gff")):
                    stat_ta_tss(ta_file, os.path.join(tss_folder, tss),
                                stat_tss_out, self.tmps["ta_tss"],
                                self.tmps["tss_ta"], fuzzy)
                    os.remove(ta_file)
                    os.remove(os.path.join(tss_folder, tss))
                    self.helper.sort_gff(self.tmps["ta_tss"], ta_file)
                    self.helper.sort_gff(self.tmps["tss_ta"],
                                         os.path.join(tss_path, tss))
                    os.remove(self.tmps["tss_ta"])
                    os.remove(self.tmps["ta_tss"])

    def _compare_cds(self, tas, gff_outfolder, cds_path, stat_path):
        self.multiparser.parser_gff(cds_path, None)
        self.multiparser.combine_gff(gff_outfolder,
                         os.path.join(cds_path, "tmp"), "transcript", None)
        print("Comaring of Transcript assembly and CDS file...")
        cds_folder = os.path.join(cds_path, "tmp")
        for ta in tas:
            ta_file = os.path.join(gff_outfolder,
                                   "_".join([ta, self.endfix_tran]))
            stat_gff_out = os.path.join(stat_path,
                           "".join(["stat_compare_Transcriptome_assembly_CDS_",
                           ta, ".csv"]))
            for gff in os.listdir(cds_folder):
                if (gff[:-4] == ta) and (gff.endswith(".gff")):
                    cds_file = os.path.join(cds_folder, gff)
                    stat_ta_gff(ta_file, cds_file, stat_gff_out,
                                self.tmps["ta_gff"], self.tmps["gff_ta"])
                    os.remove(ta_file)
                    os.remove(os.path.join(cds_path, gff))
                    self.helper.sort_gff(self.tmps["ta_gff"], ta_file)
                    self.helper.sort_gff(self.tmps["gff_ta"],
                                         os.path.join(cds_path, gff))
                    os.remove(self.tmps["ta_gff"])
                    os.remove(self.tmps["gff_ta"])

    def _compare_tss_cds(self, compare_tss, compare_cds, gff_outfolder,
                         stat_path, fuzzy, tas):
        if compare_tss is not None and compare_cds is not None:
            self.multiparser.parser_gff(gff_outfolder, "transcript")
            self._compare_cds(tas, gff_outfolder, compare_cds, stat_path)
            self._compare_tss(tas, gff_outfolder, compare_tss, stat_path, fuzzy)
        elif compare_cds is not None and compare_tss is None:
            self.multiparser.parser_gff(gff_outfolder, "transcript")
            self._compare_cds(tas, gff_outfolder, compare_cds, stat_path)
        elif compare_cds is None and compare_tss is not None:
            self.multiparser.parser_gff(gff_outfolder, "transcript")
            self._compare_tss(tas, gff_outfolder, compare_tss, stat_path, fuzzy)

    def _for_one_wig(self, type_, wigs, tex, height, width, tolerance,
                     replicates, out_folder, libs, gff_outfolder, low_cutoff):
        print("Computing {0} wig files....".format(type_))
        folder = wigs.split("/")
        folder = "/".join(folder[:-1])
        merges = os.path.join(folder, "merge_wigs")
        self.helper.check_make_folder(merges)
        wig_path = os.path.join(wigs, "tmp")
        self.multiparser.parser_wig(wigs)
        strains = self._compute(merges, wigs, tex, height, width, tolerance,
                        replicates, out_folder, type_, libs, low_cutoff)
        for strain in strains:
            out = os.path.join(gff_outfolder,
                  "_".join([strain, "transcript_assembly", type_ + ".gff"]))
            self.helper.sort_gff(os.path.join(out_folder,
                                 "_".join([strain, type_])), out)
            os.remove(os.path.join(out_folder, "_".join([strain, type_])))
        return strains

    def _for_two_wigs(self, frag_wigs, tex_wigs, strains,
                      gff_outfolder, tolerance):
        if (frag_wigs is not None) and (tex_wigs is not None):
            print("merge fragment and tex treat one ....")
            for strain in strains:
                frag_gff = os.path.join(gff_outfolder,
                           "_".join([strain, self.frag]))
                tex_gff = os.path.join(gff_outfolder,
                          "_".join([strain, self.tex]))
                final_gff = os.path.join(gff_outfolder,
                            "_".join([strain, self.endfix_tran]))
                for gff in os.listdir(gff_outfolder):
                    if "transcript_assembly" in gff:
                        filename = gff.split("_transcript_assembly_")
                        if (strain == filename[0]) and (
                            "tex.gff" == filename[1]):
                            tex_file = gff
                        elif (strain == filename[0]) and (
                              "fragment.gff" == filename[1]):
                            frag_file = gff
                combine(os.path.join(gff_outfolder, frag_file),
                        os.path.join(gff_outfolder, tex_file), tolerance,
                        os.path.join(gff_outfolder,
                                     "_".join([strain, self.endfix_tran])))
                os.remove(frag_gff)
                os.remove(tex_gff)
        else:
            if frag_wigs is not None:
                for strain in strains:
                    frag_gff = os.path.join(gff_outfolder,
                               "_".join([strain, self.frag]))
                    final_gff = os.path.join(gff_outfolder,
                                "_".join([strain, self.endfix_tran]))
                    shutil.move(frag_gff, final_gff)
            elif tex_wigs is not None:
                for strain in strains:
                    tex_gff = os.path.join(gff_outfolder,
                              "_".join([strain, self.tex]))
                    final_gff = os.path.join(gff_outfolder,
                                "_".join([strain, self.endfix_tran]))
                    shutil.move(tex_gff, final_gff)

    def _post_modify(self, tas, gffs, tran_path, length,
                     gff_outfolder, out_folder):
        for ta in tas:
            for gff in os.listdir(gffs):
                if (".gff" in gff) and (gff[:-4] == ta):
                    break
            print("Modifying {0} refering to {1}...".format(ta, gff))
            fill_gap(os.path.join(gffs, gff),
                     os.path.join(tran_path, "_".join([ta, self.endfix_tran])),
                     "overlap", self.tmps["overlap"])
            fill_gap(os.path.join(gffs, gff),
                     os.path.join(tran_path, "_".join([ta, self.endfix_tran])),
                     "uni", self.tmps["uni"])
            tmp_merge = os.path.join(gff_outfolder, self.tmps["merge"])
            if self.tmps["merge"] in gff_outfolder:
                os.remove(tmp_merge)
            self.helper.merge_file(self.tmps["overlap"], tmp_merge)
            self.helper.merge_file(self.tmps["uni"], tmp_merge)
            tmp_out = os.path.join(gff_outfolder, "_".join(["tmp", ta]))
            self.helper.sort_gff(tmp_merge, tmp_out)
            os.remove(self.tmps["overlap"])
            os.remove(self.tmps["uni"])
            os.remove(tmp_merge)
            final_out = os.path.join(gff_outfolder, "_".join(["final", ta]))
            longer_ta(tmp_out, length, final_out)
            shutil.move(final_out, os.path.join(self.tmps["tran"],
                      "_".join([ta, self.endfix_tran])))
            os.remove(tmp_out)
        shutil.rmtree(gff_outfolder)
        shutil.move(self.tmps["tran"], gff_outfolder)

    def _remove_file(self, frag_wigs, tex_wigs, gffs,
                     compare_cds, compare_tss, out_folder):
        if frag_wigs is not None:
            self.helper.remove_wigs(frag_wigs)
        if tex_wigs is not None:
            self.helper.remove_wigs(tex_wigs)
        if gffs is not None:
            self.helper.remove_tmp(gffs)
        if compare_cds is not None:
            self.helper.remove_tmp(compare_cds)
        if compare_tss is not None:
            self.helper.remove_tmp(compare_tss)
        self.helper.remove_tmp(os.path.join(out_folder, "gffs"))
        self.helper.remove_tmp(self.gff_outfolder)

    def run_transcript_assembly(self, frag_wigs, tex_wigs, sort, tex, length,
            gffs, height, width, tolerance, low_cutoff, replicates_tex,
            replicates_frag, out_folder, compare_tss, compare_cds, fuzzy,
            tlibs, flibs):
        if (replicates_tex is not None) and (
            replicates_frag is not None):
            replicates = {"tex": int(replicates_tex),
                          "frag": int(replicates_frag)}
        elif replicates_tex is not None:
            replicates = {"tex": int(replicates_tex), "frag": -1}
        elif replicates_frag is not None:
            replicates = {"tex": -1, "frag": int(replicates_frag)}
        else:
            print("Error:No replicates number assign!!!")
            sys.exit()
        if (frag_wigs is None) and (tex_wigs is None):
            print("Error: there is no wigs files!!!!\n")
            sys.exit()
        if frag_wigs is not None:
            strains = self._for_one_wig("fragment", frag_wigs, tex, height,
                            width, tolerance, replicates, out_folder, flibs,
                            self.gff_outfolder, low_cutoff)
        if tex_wigs is not None:
            strains = self._for_one_wig("tex", tex_wigs, tex, height, width,
                            tolerance, replicates, out_folder, tlibs,
                             self.gff_outfolder, low_cutoff)
        self._for_two_wigs(frag_wigs, tex_wigs, strains,
                           self.gff_outfolder, tolerance)
        tas = []
        if gffs is not None:
            if sort:
                for gff in os.listdir(gffs):
                    if gff.endswith(".gff"):
                        self.helper.sort_gff(os.path.join(gffs, gff), self.tmps["gff"])
                        shutil.move(self.tmps["gff"], os.path.join(gffs, gff))
            self.multiparser.parser_gff(gffs, None)
            self.multiparser.combine_gff(gffs, os.path.join(gffs, "tmp"), None, None)
            self.multiparser.parser_gff(self.gff_outfolder, "transcript")
            self.multiparser.combine_gff(gffs, self.tran_path, None, "transcript")
            self.helper.check_make_folder(self.tmps["tran"])
            for ta in os.listdir(self.tran_path):
                if ta.endswith(".gff"):
                    if os.path.getsize(os.path.join(self.tran_path, ta)) != 0:
                        tas.append(ta.replace("_" + self.endfix_tran, ""))
            self._post_modify(tas, gffs, self.tran_path, length,
                              self.gff_outfolder, out_folder)
        self._compare_tss_cds(compare_tss, compare_cds, self.gff_outfolder,
                              self.stat_path, fuzzy, tas)
        self._remove_file(frag_wigs, tex_wigs, gffs, compare_cds,
                          compare_tss, out_folder)
