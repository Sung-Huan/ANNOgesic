#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
from subprocess import call
from transaplib.helper import Helper
from transaplib.multiparser import Multiparser
from transaplib.converter import Converter
from transaplib.combine_frag_tex import Combine
from transaplib.stat_TA_comparison import Stat_TA_TSS, Stat_TA_GFF
from transaplib.transcript_assembly import Assembly
from transaplib.fill_gap import Fill_gap, Longer_TA


class Transcript_Assembly(object):

    def __init__(self):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.converter = Converter()    

    def _compute_transcript(self, wig_f, wig_r, height, width, wig_folder, out_folder,
                            tolerance, support, wig_type, strain, libs, tex):
        print("Computing transcript assembly for " + strain)
        out = out_folder + "/" + strain + "_" + wig_type
        Assembly(wig_f, wig_r, height, width, tolerance, wig_folder, tex, libs, support, out)

    def _combine_wigs(self, wig_folder, merge, out, direct, strain):
        for wig in os.listdir(wig_folder):
            if (direct in wig):
                filename = wig.split("_STRAIN_")
                if filename[1][:-4] == strain:
                    print("Merge " + wig + " now....")
                    call(["cat", wig_folder + "/" + wig], stdout=out)

    def _compute(self, merges, wig_folder, tex, height, width, tolerance, 
                 support, out_folder, wig_type, libs):
        strains = []
        wigs = wig_folder + "/tmp/"
        for wig in os.listdir(wigs):
            filename = wig.split("_STRAIN_")
            if filename[1][:-4] not in strains:
                strains.append(filename[1][:-4])
        for strain in strains:
            f_file = merges + "/" + strain + "_forward.wig"
            r_file = merges + "/" + strain + "_reverse.wig"
            out_f = open(f_file, "w")
            out_r = open(r_file, "w")
            self._combine_wigs(wigs, merges, out_f, "forward", strain)
            print("Merge " + strain + "_forward.wig complete...")
            self._combine_wigs(wigs, merges, out_r, "reverse", strain)
            print("Merge " + strain + "_reverse.wig complete...")
            out_f.close()
            out_r.close()
        #######################################################
        # compute_transcript is seperated from previous loop  #
        # because compute_transcript is based on the number   #
        #######################################################
        for strain in strains:
            f_file = merges + "/" + strain + "_forward.wig"
            r_file = merges + "/" + strain + "_reverse.wig"
            self._compute_transcript(f_file, r_file, height, width, 
                                     wig_folder, out_folder, tolerance, support,
                                     wig_type, strain, libs, tex)
        return strains

    def _compare_TSS(self, tas, gff_outfolder, tss_path, stat_path, fuzzy):
        self.multiparser._parser_gff(tss_path, "TSS")
        self.multiparser._combine_gff(gff_outfolder, tss_path + "/tmp", "transcript", "TSS")
        print("Comaring of Transcript assembly and TSS file...")
        tss_folder = tss_path + "/tmp/"
        for ta in tas:
            ta_file = gff_outfolder + ta + "_transcript.gff"
            stat_tss_out = stat_path + "stat_compare_Transcriptome_assembly_TSS_" + ta + ".csv"
            for tss in os.listdir(tss_folder):
                filename = tss.split("_TSS")
                if (filename[0] == ta) and (tss.endswith(".gff")):
                    Stat_TA_TSS(ta_file, tss_folder + tss, stat_tss_out, "tmp_ta_tss", "tmp_tss", fuzzy)
                    call(["rm", ta_file])
                    call(["rm", tss_folder + "/" + tss])
                    self.helper.sort_gff("tmp_ta_tss", ta_file)
                    self.helper.sort_gff("tmp_tss", tss_path + "/" + tss)
                    call(["rm", "tmp_tss"])
                    call(["rm", "tmp_ta_tss"])
    def _compare_CDS(self, tas, gff_outfolder, cds_path, stat_path):
        self.multiparser._parser_gff(cds_path, None)
        self.multiparser._combine_gff(gff_outfolder, cds_path + "/tmp", "transcript", None)
        print("Comaring of Transcript assembly and CDS file...")
        cds_folder = cds_path + "/tmp/"
        for ta in tas:
            ta_file = gff_outfolder + ta + "_transcript.gff"
            stat_gff_out = stat_path + "stat_compare_Transcriptome_assembly_CDS_" + ta + ".csv"
            for gff in os.listdir(cds_folder):
                if (gff[:-4] == ta) and (gff.endswith(".gff")):
                    cds_file = cds_folder + gff
                    Stat_TA_GFF(ta_file, cds_file, stat_gff_out, "tmp_ta_gff", "tmp_gff_ta")
                    call(["rm", ta_file])
                    call(["rm", cds_path + "/" + gff])
                    self.helper.sort_gff("tmp_ta_gff", ta_file)
                    self.helper.sort_gff("tmp_gff_ta", cds_path + "/" + gff)
                    self.helper.sort_gff("tmp_ta_gff", "tmp1" + ta)
                    self.helper.sort_gff("tmp_gff_ta", "tmp2" + ta)
                    call(["rm", "tmp_ta_gff"])
                    call(["rm", "tmp_gff_ta"])

    def _compare_TSS_CDS(self, compare_TSS, compare_CDS, gff_outfolder, stat_path, fuzzy, tas):
        if (compare_TSS is not False) and \
           (compare_CDS is not False):
            self.multiparser._parser_gff(gff_outfolder, "transcript")
            self._compare_CDS(tas, gff_outfolder, compare_CDS, stat_path)
            self._compare_TSS(tas, gff_outfolder, compare_TSS, stat_path, fuzzy)
        elif (compare_CDS is not False) and \
             (compare_TSS is False):
            self.multiparser._parser_gff(gff_outfolder, "transcript")
            self._compare_CDS(tas, gff_outfolder, compare_CDS, stat_path)
        elif (compare_CDS is False) and \
             (compare_TSS is not False):
            self.multiparser._parser_gff(gff_outfolder, "transcript")
            self._compare_TSS(tas, gff_outfolder, compare_TSS, stat_path, fuzzy)

    def _remove_wigs(self, wigs):
        folder = wigs.split("/")
        folder = "/".join(folder[:-1])
        call(["rm", "-rf", folder + "/merge_wigs"])
        call(["rm", "-rf", wigs + "/tmp"])
        os.system("rm -rf " + wigs + "/*_folder")

    def _for_one_wig(self, type_, wigs, tex, height, width, tolerance, support,
                         out_folder, libs, gff_outfolder):
        print("Computing " + type_ + " wig files....")
        folder = wigs.split("/")
        folder = "/".join(folder[:-1])
        merges = folder + "/merge_wigs"
        if "merge_wigs" in os.listdir(folder):
            call(["rm", "-rf", merges])
        call(["mkdir", merges])
        wig_path = wigs + "/tmp/"
        self.multiparser._parser_wig(wigs)
        strains = self._compute(merges, wigs, tex, height, width, tolerance, 
                                support, out_folder, type_, libs)
        for strain in strains:
            out = gff_outfolder + strain + "_transcript_assembly_" + type_ + ".gff"
            self.helper.sort_gff(out_folder + "/" + strain + "_" + type_, out)
            call(["rm", out_folder + "/" + strain + "_" + type_])
        return strains

    def _for_two_wigs(self, frag_wigs, tex_wigs, strains, gff_outfolder, tolerance):
        if (frag_wigs is not False) and (tex_wigs is not False):
            print("merge fragment and tex treat one ....")
            for strain in strains:
                print(strain)
                for gff in os.listdir(gff_outfolder):
                    if "transcript_assembly" in gff:
                        filename = gff.split("_transcript_assembly_")
                        if (strain == filename[0]) and \
                           ("tex.gff" == filename[1]):
                            tex_file = gff
                        elif (strain == filename[0]) and \
                             ("fragment.gff" == filename[1]):
                            frag_file = gff
                Combine(gff_outfolder + frag_file, 
                        gff_outfolder + tex_file, tolerance,
                        gff_outfolder + strain + "_transcript.gff")
                call(["rm", gff_outfolder + strain + "_transcript_assembly_fragment.gff"])
                call(["rm", gff_outfolder + strain + "_transcript_assembly_tex.gff"])
        else:
            if frag_wigs is not False:
                for strain in strains:
                    call(["mv", gff_outfolder + strain + "_transcript_assembly_fragment.gff",
                          gff_outfolder + strain + "_transcript.gff"])
            elif tex_wigs is not False:
                for strain in strains:
                    call(["mv", gff_outfolder + strain + "_transcript_assembly_tex.gff",
                          gff_outfolder + strain + "_transcript.gff"])

    def _post_modify(self, tas, gffs, tran_path, length, gff_outfolder):
        for ta in tas:
            overlap = open("tmp_overlap", "w")
            uni = open("tmp_uni", "w")
            final = open("tmp_merge", "w")
            for gff in os.listdir(gffs):
                if (".gff" in gff) and (gff[:-4] == ta):
                    break
            print("Modifying " + ta + " refering to " + gff + "...")
            Fill_gap(gffs + "/" + gff, tran_path + ta + "_transcript.gff",
                     "overlap", "tmp_overlap")
            Fill_gap(gffs + "/" + gff, tran_path + ta + "_transcript.gff",
                     "uni", "tmp_uni")
            call(["cat", "tmp_overlap"], stdout=final)
            call(["cat", "tmp_uni"], stdout=final)
            overlap.close()
            uni.close()
            final.close()
            output = "tmp_" + ta
            self.helper.sort_gff("tmp_merge", output)
            call(["rm", "tmp_overlap"])
            call(["rm", "tmp_uni"])
            call(["rm", "tmp_merge"])
            Longer_TA("tmp_" + ta, length, "final_" + ta)
            print(gff_outfolder + ta + "_transcript.gff")
            call(["mv", "final_" + ta, "tmp_tran/" + ta + "_transcript.gff"])
            call(["rm", "tmp_" + ta])
        os.system("rm -rf " + gff_outfolder)
        call(["mv", "tmp_tran", gff_outfolder])

    def run_Transcript_Assembly(self, project_path, bin_path, frag_wigs, tex_wigs, sort, 
                                tex, length, gffs, height, width, tolerance, support, out_folder, 
                                compare_TSS, compare_CDS, fuzzy, replicates, tlibs, flibs):
        gff_outfolder = out_folder + "/gffs/"
        stat_path = out_folder + "/statistics/"
        if (frag_wigs is False) and (tex_wigs is False):
            print("Error: there is no wigs files!!!!\n")
            sys.exit()
        if frag_wigs is not False:
            strains = self._for_one_wig("fragment", frag_wigs, tex, height, width, tolerance, 
                                        support, out_folder, flibs, gff_outfolder)
        if tex_wigs is not False:
            strains = self._for_one_wig("tex", tex_wigs, tex, height, width, tolerance, 
                                        support, out_folder, tlibs, gff_outfolder)
        self._for_two_wigs(frag_wigs, tex_wigs, strains, gff_outfolder, tolerance)
        if gffs is not False:
            if sort is True:
                for gff in os.listdir(gffs):
                    if gff.endswith(".gff"):
                        self.helper.sort_gff(gffs + "/" + gff, "tmp_sort_gff")
                        call(["mv", "tmp_sort_gff", gffs + "/" + gff])
            self.multiparser._parser_gff(gffs, None)
            self.multiparser._combine_gff(gffs, gffs + "/tmp", None, None)
            self.multiparser._parser_gff(gff_outfolder, "transcript")
            self.multiparser._combine_gff(gffs, gff_outfolder + "/tmp", None, "transcript")
            tran_path = gff_outfolder + "/tmp/"
            tas = []
            self.helper.check_make_folder(os.getcwd() + "/", "tmp_tran")
            for ta in os.listdir(tran_path):
                if ta.endswith(".gff"):
                    if os.path.getsize(tran_path + ta) != 0:
                        tas.append(ta.replace("_transcript.gff", ""))
            self._post_modify(tas, gffs, tran_path, length, gff_outfolder)
        self._compare_TSS_CDS(compare_TSS, compare_CDS, gff_outfolder, stat_path, fuzzy, tas)
        ####################################################
        # Remove temperary files and folders               #
        ####################################################
        if frag_wigs is not False:
            self.helper.remove_wigs(frag_wigs)
        if tex_wigs is not False:
            self.helper.remove_wigs(tex_wigs)
        if gffs is not False:
            self.helper.remove_tmp(gffs)
        if compare_CDS is not False:
            self.helper.remove_tmp(compare_CDS)
        if compare_TSS is not False:
            self.helper.remove_tmp(compare_TSS)
        self.helper.remove_tmp(out_folder + "/gffs")
