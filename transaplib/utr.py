#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
from subprocess import call
import transaplib.multiparser
from transaplib.helper import Helper
from transaplib.detect_utr import Detect_3UTR, Detect_5UTR
from transaplib.multiparser import Multiparser


class UTR_detection(object):

    def __init__(self):
        self.helper = Helper()
        self.multiparser = Multiparser()

    def _check_folder(self, folder):
        if folder is None:
            print("Error: lack required files!!!")
            sys.exit()

    def _check_gff(self, folder):
        for gff in os.listdir(folder):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(folder, gff))

    def _compute_utr(self, gffs, utr5_path, utr3_path, tss_path, tran_path, 
                     source, terms, fuzzy):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                prefix = gff[:-4]
                tss = self.helper.get_correct_file(tss_path, "_TSS.gff", prefix, None)
                tran = self.helper.get_correct_file(tran_path, "_transcript.gff", prefix, None)
                if terms is not False:
                    term = self.helper.get_correct_file(os.path.join(terms, "tmp"), "_term.gff", prefix, None)
                else:
                    term = None
                print("computing 5'UTR of {0} .....".format(prefix))
                Detect_5UTR(tss, os.path.join(gffs, gff), tran, source,
                            os.path.join(utr5_path, "gffs", "_".join([prefix, "5UTR.gff"])))                
                print("computing 3'UTR of {0} .....".format(prefix))
                Detect_3UTR(tran, os.path.join(gffs, gff), term, fuzzy,
                            os.path.join(utr3_path, "gffs", "_".join([prefix, "3UTR.gff"])))
                os.helper.move_all_content(os.getcwd(), os.path.join(utr5_path, "statistics"), "_5utr_length.png")
                os.helper.move_all_content(os.getcwd(), os.path.join(utr3_path, "statistics"), "_3utr_length.png")

    def run_UTR_detection(self, bin_path, tsss, gffs, trans, terms, 
                          fuzzy, out_folder, source):
        self._check_folder(tsss)
        self._check_folder(gffs)
        self._check_folder(trans)
        self._check_gff(tsss)
        self._check_gff(gffs)
        self._check_gff(trans)
        self._check_gff(terms)
        tss_path = os.path.join(tsss, "tmp")
        tran_path = os.path.join(trans, "tmp")
        self.multiparser._parser_gff(gffs, None)
        self.multiparser._parser_gff(tsss, "TSS")
        self.multiparser._combine_gff(gffs, tss_path, None, "TSS")
        self.multiparser._parser_gff(trans, "transcript")
        self.multiparser._combine_gff(gffs, tran_path, None, "transcript")
        utr5_path = os.path.join(out_folder, "5UTR")
        utr3_path = os.path.join(out_folder, "3UTR")
        if terms is not False:
            self.multiparser._parser_gff(terms, "term")
            self.multiparser._combine_gff(gffs, os.path.join(terms, "tmp"), None, "term")
        detect = False
        self._compute_utr(gffs, utr5_path, utr3_path, tss_path, 
                          tran_path, source, terms, fuzzy)
        self.helper._remove_tmp(gffs)
        self.helper._remove_tmp(tsss)
        self.helper._remove_tmp(trans)
        self.helper._remove_tmp(terms)
