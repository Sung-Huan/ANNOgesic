import os
import sys
import shutil
from subprocess import call
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.sORF_intergenic import get_intergenic
from annogesiclib.sORF_detection import sorf_detection
from annogesiclib.stat_sorf import stat


class sORFDetection(object):

    def __init__(self, tsss, srnas, out_folder, trans, fastas):
        self.multiparser = Multiparser()
        self.helper = Helper()
        if tsss is not None:
            self.tss_path = os.path.join(tsss, "tmp")
        else:
            self.tss_path = None
        if srnas is not None:
            self.srna_path = os.path.join(srnas, "tmp")
        else:
            self.srna_path = None
        self.gff_output = os.path.join(out_folder, "gffs")
        self.table_output = os.path.join(out_folder, "tables")
        self.tran_path = os.path.join(trans, "tmp")
        self.fasta_path = os.path.join(fastas, "tmp")
        self.all_cand = "all_candidates"
        self.best = "best"

    def _check_gff(self, gffs):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(gffs, gff))

    def _check_necessary_files(self, gffs, trans, tex_wigs, frag_wigs,
                               utr_detect, tsss, srnas):
        if (gffs is None) or \
           (trans is None) or \
           ((tex_wigs is None) and (frag_wigs is None)):
            print("Error: lack required files!!!!")
            sys.exit()
        if utr_detect:
            if (tsss is None):
                print("Error: lack required files for UTR derived sORF detection!!!!")
                sys.exit()
        self._check_gff(gffs)
        self.multiparser.parser_gff(gffs, None)
        if tsss is not None:
            self._check_gff(tsss)
            self.multiparser.parser_gff(tsss, "TSS")
            self.multiparser.combine_gff(gffs, self.tss_path, None, "TSS")
        self._check_gff(trans)
        if srnas is not None:
            self._check_gff(srnas)
            self.multiparser.parser_gff(srnas, "sRNA")
            self.multiparser.combine_gff(gffs, self.srna_path, None, "sRNA")

    def _merge_wigs(self, tex_wigs, frag_wigs):
        if (tex_wigs is not None) and (frag_wigs is not None):
            folder = tex_wigs.split("/")
            folder = "/".join(folder[:-1])
            merge_wigs = os.path.join(folder, "merge_wigs")
            self.helper.check_make_folder(merge_wigs)
            for wig in os.listdir(tex_wigs):
                if os.path.isdir(os.path.join(tex_wigs, wig)):
                    pass
                else:
                    shutil.copy(os.path.join(tex_wigs, wig), merge_wigs)
            for wig in os.listdir(frag_wigs):
                if os.path.isdir(os.path.join(frag_wigs, wig)):
                    pass
                else:
                    shutil.copy(os.path.join(frag_wigs, wig), merge_wigs)
        elif (tex_wigs is not None):
            merge_wigs = tex_wigs
        elif (frag_wigs is not None):
            merge_wigs = frag_wigs
        return merge_wigs

    def _combine_libs(self, tlibs, flibs):
        if (tlibs is None) and (flibs is None):
            print("Error: please input proper libraries!!")
        if (tlibs is not None) and (flibs is not None):
            libs = tlibs + flibs
        elif (tlibs is not None):
            libs = tlibs
        elif (flibs is not None):
            libs = flibs
        return libs

    def _start_stop_codon(self, prefixs, out_folder, utr_length, libs,
                          tex_notex, replicate, cutoff_inter, cutoff_3utr,
                          cutoff_5utr, cutoff_intercds, wig_path, merge_wigs,
                          start_codon, stop_codon, max_len, min_len, table_best,
                          utr_detect, background, fuzzy_rbs,
                          print_all, no_srna, noafter_tss, no_tss):
        for prefix in prefixs:
            if self.srna_path is not None:
                srna_file = os.path.join(self.srna_path,
                                         "_".join([prefix, "sRNA.gff"]))
            else:
                srna_file = None
            if self.tss_path is not None:
                tss_file = os.path.join(self.tss_path,
                                        "_".join([prefix, "TSS.gff"]))
            else:
                tss_file = None
            sorf_detection(os.path.join(self.fasta_path, prefix + ".fa"),
                           srna_file, os.path.join(out_folder,
                                      "_".join([prefix, "inter.gff"])),
                           tss_file, utr_length, utr_detect, libs, tex_notex,
                           replicate, cutoff_inter, cutoff_3utr, cutoff_5utr,
                           cutoff_intercds, os.path.join(wig_path,
                                            "_".join([prefix, "forward.wig"])),
                           os.path.join(wig_path,
                           "_".join([prefix, "reverse.wig"])), merge_wigs,
                           start_codon, stop_codon, table_best, max_len, min_len,
                           os.path.join(self.gff_output, self.all_cand,
                           "_".join([prefix, "sORF"])), background,
                           fuzzy_rbs, print_all, no_srna, noafter_tss, no_tss)
            if "_".join([prefix, "sORF_all.gff"]) in os.listdir(
                         os.path.join(self.gff_output, self.all_cand)):
                shutil.move(os.path.join(self.gff_output, self.all_cand,
                            "_".join([prefix, "sORF_all.gff"])),
                            os.path.join(self.gff_output, self.all_cand,
                            "_".join([prefix, "sORF.gff"])))
                shutil.move(os.path.join(self.gff_output, self.all_cand,
                            "_".join([prefix, "sORF_best.gff"])),
                            os.path.join(self.gff_output, self.best,
                            "_".join([prefix, "sORF.gff"])))
                shutil.move(os.path.join(self.gff_output, self.all_cand,
                            "_".join([prefix, "sORF_all.csv"])),
                            os.path.join(self.table_output, self.all_cand,
                            "_".join([prefix, "sORF.csv"])))
                shutil.move(os.path.join(self.gff_output, self.all_cand,
                            "_".join([prefix, "sORF_best.csv"])),
                            os.path.join(self.table_output, self.best,
                            "_".join([prefix, "sORF.csv"])))

    def _remove_tmp(self, out_folder, fastas, gffs, tsss,
                    trans, srnas, tex_wigs, frag_wigs):
        self.helper.remove_all_content(out_folder, ".gff", "file")
        self.helper.remove_tmp(fastas)
        self.helper.remove_tmp(gffs)
        self.helper.remove_tmp(tsss)
        self.helper.remove_tmp(trans)
        self.helper.remove_tmp(srnas)
        self.helper.remove_wigs(tex_wigs)
        self.helper.remove_wigs(frag_wigs)

    def _compare_tran_cds(self, gffs, out_folder, utr_detect, hypo):
        prefixs = []
        for gff in os.listdir(gffs): ## compare transcript and CDS
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                prefixs.append(prefix)
                print("comparing transcript and CDS of {0}".format(prefix))
                get_intergenic(os.path.join(gffs, gff),
                               os.path.join(self.tran_path,
                               "_".join([prefix, "transcript.gff"])),
                               os.path.join(out_folder,
                               "_".join([prefix, "inter.gff"])),
                               utr_detect, hypo)
        return prefixs

    def _merge_libs(self, flibs, tlibs):
        if (flibs is not None) and (tlibs is not None):
            input_libs = flibs + tlibs
        elif (flibs is not None):
            input_libs = flibs
        elif (tlibs is not None):
            input_libs = tlibs
        else:
            print("Error: No libs be assigned!!!")
            sys.exit()
        return input_libs

    def run_sorf_detection(self, out_folder, utr_detect, trans, gffs, tsss,
            utr_length, min_len, max_len, tex_wigs, frag_wigs, cutoff_inter,
            cutoff_5utr, cutoff_3utr, cutoff_intercds, fastas, tlibs, flibs,
            tex_notex, replicates_tex, replicates_frag, table_best, srnas,
            start_codon, stop_codon, background, fuzzy_rbs, noafter_tss,
            print_all, no_srna, no_tss, hypo):
        if fuzzy_rbs > 6:
            print("Error: --fuzzy_rbs should be equal or less than 6!!")
            sys.exit()
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
        self._check_necessary_files(gffs, trans, tex_wigs,
                                    frag_wigs, utr_detect, tsss, srnas)
        merge_wigs = self._merge_wigs(tex_wigs, frag_wigs)
        input_libs = self._merge_libs(flibs, tlibs)
        wig_path = os.path.join(merge_wigs, "tmp")
        self.multiparser.parser_wig(merge_wigs)
        self.multiparser.combine_wig(gffs, wig_path, None, input_libs)
        self.multiparser.parser_gff(trans, "transcript")
        self.multiparser.combine_gff(gffs, self.tran_path,
                                      None, "transcript")
        self.multiparser.parser_fasta(fastas)
        self.multiparser.combine_fasta(gffs, self.fasta_path, None)
        libs = self._combine_libs(tlibs, flibs)
        prefixs = self._compare_tran_cds(gffs, out_folder, utr_detect, hypo)
        self._start_stop_codon(prefixs, out_folder, utr_length, libs,
             tex_notex, replicates, cutoff_inter, cutoff_3utr, cutoff_5utr,
             cutoff_intercds, wig_path, merge_wigs, start_codon, stop_codon,
             max_len, min_len, table_best, utr_detect, background,
             fuzzy_rbs, print_all, no_srna, noafter_tss, no_tss)
        for sorf in os.listdir(os.path.join(self.gff_output, self.all_cand)):
            print("statistics of {0}".format(sorf))
            if sorf.endswith("_sORF.gff"):
                stat(os.path.join(self.gff_output, self.all_cand, sorf),
                     os.path.join(self.gff_output, self.best, sorf),
                     os.path.join(out_folder, "statistics",
                     "_".join(["stat", sorf.replace(".gff", ".csv")])),
                     utr_detect)
        self._remove_tmp(out_folder, fastas, gffs, tsss, trans,
                         srnas, tex_wigs, frag_wigs)
