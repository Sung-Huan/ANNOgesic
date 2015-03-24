#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
from subprocess import call
import transaplib.multiparser
import transaplib.helper
class UTR_detection(object):

    def _check_make_folder(self, path, folder):
        if folder in os.listdir(path):
            call(["rm", "-rf", path + folder])
        call(["mkdir", path + folder])

    def _check_folder(self, folder):
        if folder is None:
            print("Error: lack required files!!!")
            sys.exit()

    def _check_gff(self, gffs):
        for gff in os.listdir(gffs):
            if type_ == "term":
                if gff.endswith("_term.gff"):
                    check_gff_attributes._check_uni_attributes(gffs + "/" + gff)
            else:
                if gff.endswith(".gff"):
                    check_gff_attributes._check_uni_attributes(gffs + "/" + gff)

    def _get_correct_file(self, folder, prefix, feature):
        detect = False
        for file_ in os.listdir(folder):
            if prefix == file_.replace("_" + feature + ".gff", ""):
                detect = True
                break
        if detect:
            detect = False
            return folder + file_
        else:
            if feature == "term":
                print("Warning: No terminator file of " + \
                      prefix + "_" + feature + ".gff!!")
                return False
            else:
                print("Error: No " + feature + " file of " + \
                      prefix + "_" + feature + ".gff!!")
                sys.exit()

    def _remove_tmp(self, folder):
        if folder is not False:
            call(["rm", "-rf", folder + "/tmp"])
            os.system("rm -rf " + folder + "/*_folder")

    def run_UTR_detection(self, bin_path, tsss, gffs, trans, terms, 
                          fuzzy, out_folder, source):
        self._check_folder(tsss)
        self._check_folder(gffs)
        self._check_folder(trans)
        self._check_gff(tsss, "tss")
        self._check_gff(gffs, "gff")
        self._check_gff(trans, "tran")
        self._check_gff(terms, "term")
        tss_path = tsss + "/tmp/"
        tran_path = trans + "/tmp/"
        multiparser._parser_gff(gffs, None)
        multiparser._parser_gff(tsss, "TSS")
        multiparser._combine_gff(gffs, tss_path, None, "TSS")
        multiparser._parser_gff(trans, "transcript")
        multiparser._combine_gff(gffs, tran_path, None, "transcript")
        utr5_path = out_folder + "/5UTR/"
        utr3_path = out_folder + "/3UTR/"
        if terms is not False:
            term_path = terms + "/tmp/"
            multiparser._parser_gff(terms, "term")
            multiparser._combine_gff(gffs, term_path, None, "term")
        detect = False
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                prefix = gff[:-4]
                out_5 = open(utr5_path + "gffs/" + prefix + "_5UTR.gff", "w")
                out_3 = open(utr3_path + "gffs/" + prefix + "_3UTR.gff", "w")
                tss = self._get_correct_file(tss_path, prefix, "TSS")
                tran = self._get_correct_file(tran_path, prefix, "transcript")
                if terms is not False:
                    term = self._get_correct_file(term_path, prefix, "term")
                print("computing 5'UTR of " + prefix + " .....")
                if source is True:
                    call(["python", bin_path + "/detect_5utr.py",
                          "-t", tss, "-g", gffs + "/" + gff, "-ta", tran], stdout=out_5)
                else:
                    call(["python", bin_path + "/detect_5utr.py",
                          "-t", tss, "-g", gffs + "/" + gff, "-ta", tran, "-s"], stdout=out_5)
                print("computing 3'UTR of " + prefix + " .....")
                if term is not False:
                    call(["python", bin_path + "/detect_3utr.py",
                          "-a", tran, "-g", gffs + "/" + gff, "-t", term, "-f", str(fuzzy)], stdout=out_3)
                else:
                    call(["python", bin_path + "/detect_3utr.py",
                          "-a", tran, "-g", gffs + "/" + gff, "-f", str(fuzzy)], stdout=out_3)
                out_5.close()
                out_3.close()
                os.system("mv *_5utr_length.png " + utr5_path + "statistics/")
                os.system("mv *_3utr_length.png " + utr3_path + "statistics/")
        self._remove_tmp(gffs)
        self._remove_tmp(tsss)
        self._remove_tmp(trans)
        self._remove_tmp(terms)
