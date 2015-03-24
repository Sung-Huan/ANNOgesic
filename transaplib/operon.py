#!/usr/bin/python
from Bio import SeqIO
import os
import sys
from subprocess import call
sys.path.append(os.environ["Transap_BIN"])
from sort_gff import Sort_GFF
import multiparser
import check_gff_attributes
class Operon_detection(object):

    def _get_correct_file(self, folder, feature, prefix):
        detect = False
        for file_ in os.listdir(folder):
            if (file_.replace(feature, "") == prefix):
                detect = True
                return folder + file_
        if detect is False:
            print("Error: there is no proper file - " + prefix + feature)
            sys.exit()

    def _check_gff(self, gffs, type_):
        for gff in os.listdir(gffs):
            if type_ == "term":
                if gff.endswith("_term.gff"):
                    check_gff_attributes._check_uni_attributes(gffs + "/" + gff)
            else:
                if gff.endswith(".gff"):
                    check_gff_attributes._check_uni_attributes(gffs + "/" + gff)

    def run_Operon(self, tsss, gffs, trans, utr5s, utr3s, terms, tss_fuzzy, term_fuzzy, 
                   length, stat, output_folder, combine_gff, stat_folder, bin_path):
        ######################################################
        # parse the files and set the name of folders.       #
        ######################################################
#        self._check_gff(tsss, "tss")
#        self._check_gff(gffs, "gff")
#        self._check_gff(trans, "tran")
#        self._check_gff(utr5s, "utr")
#        self._check_gff(utr3s, "utr")
#        self._check_gff(terms, "term")
        tss_path = tsss + "/tmp/"
        tran_path = trans + "/tmp/"
        utr5_path = utr5s + "/tmp/"
        utr3_path = utr3s + "/tmp/"
        gff_path = gffs + "/"
        stat_path = stat_folder + "/"
        table_path = output_folder + "/tables/"
        if terms is not False:
            term_path = terms + "/tmp/"
        multiparser._parser_gff(gffs, None)
        multiparser._parser_gff(tsss, "TSS")
        multiparser._combine_gff(gffs, tss_path, None, "TSS")
        multiparser._parser_gff(trans, "transcript")
        multiparser._combine_gff(gffs, tran_path, None, "transcript")
        multiparser._parser_gff(utr5s, "5UTR")
        multiparser._combine_gff(gffs, utr5_path, None, "5UTR")
        multiparser._parser_gff(utr3s, "3UTR")
        multiparser._combine_gff(gffs, utr3_path, None, "3UTR")
        multiparser._parser_gff(terms, "term")
        multiparser._combine_gff(gffs, term_path, None, "term")
        prefixs = []
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                prefixs.append(gff.replace(".gff", ""))
        ########################################################
        # Detect the operons and generate the table of Operon. #
        ########################################################
        for prefix in prefixs:
            out_table = open(table_path + "operon_" + prefix + ".csv", "w")
            print("Detection operons of " + prefix)
            tss = self._get_correct_file(tss_path, "_TSS.gff", prefix)
            tran = self._get_correct_file(tran_path, "_transcript.gff", prefix)
            gff = self._get_correct_file(gff_path, ".gff", prefix)
            if terms is False:
                call(["python", bin_path + "/detect_operon.py",
                      "-tss", tss, "-ta", tran,
                      "-g", gff, "-sf", str(tss_fuzzy),
                      "-tf", str(term_fuzzy), 
                      "-l", str(length)], stdout=out_table)
            else:
                term = self._get_correct_file(term_path, "_term.gff", prefix)
                call(["python", bin_path + "/detect_operon.py",
                      "-tss", tss, "-ta", tran,
                      "-g", gff, "-term", term,
                      "-sf", str(tss_fuzzy),
                      "-tf", str(term_fuzzy), 
                      "-l", str(length)], stdout=out_table)
        ###############
        # statistics  #
        ###############
        if stat is not False:
            for table in os.listdir(table_path):
                if table.startswith("operon_") and \
                   table.endswith(".csv"):
                    filename = "stat_" + table
                    out_stat = open(stat_path + filename, "w")
                    call(["python", bin_path + "/stat_operon.py",
                          "-i", table_path + table], stdout=out_stat)
        ########################################################
        # Combine all gff files of features to one gff file.   #
        ########################################################
        if combine_gff is not False:
            for prefix in prefixs:
                out_file = open(output_folder + "/gffs/" + prefix + "_all_features.gff", "w")
                print("Combine all features of " + prefix)
                tss = self._get_correct_file(tss_path, "_TSS.gff", prefix)
                tran = self._get_correct_file(tran_path, "_transcript.gff", prefix)
                gff = self._get_correct_file(gff_path, ".gff", prefix)
                utr5 = self._get_correct_file(utr5_path, "_5UTR.gff", prefix)
                utr3 = self._get_correct_file(utr3_path, "_3UTR.gff", prefix)
                if terms is False:
                    call(["python", bin_path + "/combine_gff.py",
                          "-ts", tss, "-ta", tran,
                          "-g", gff, "-u5", utr5,
                          "-u3", utr3,
                          "-fts", str(tss_fuzzy)],
                          stdout=out_file)
                else:
                    term = self._get_correct_file(term_path, "_term.gff", prefix)
                    call(["python", bin_path + "/combine_gff.py",
                          "-ts", tss, "-ta", tran,
                          "-g", gff, "-te", term,
                          "-u5", utr5, "-u3", utr3,
                          "-fts", str(tss_fuzzy),
                          "-fte", str(term_fuzzy)],
                          stdout=out_file)
