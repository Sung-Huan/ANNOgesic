#!/usr/bin/python

import os
import sys
from subprocess import call
import annogesiclib.multiparser
from annogesiclib.helper import Helper
from annogesiclib.detect_utr import detect_3utr, detect_5utr
from annogesiclib.multiparser import Multiparser


class UTRDetection(object):

    def __init__(self, tsss, trans, out_folder):
        self.helper = Helper()
        self.multiparser = Multiparser()
        self.tss_path = os.path.join(tsss, "tmp")
        self.tran_path = os.path.join(trans, "tmp")
        self.utr5_path = os.path.join(out_folder, "5UTR")
        self.utr3_path = os.path.join(out_folder, "3UTR")
        self.utr5_stat_path = os.path.join(self.utr5_path, "statistics")
        self.utr3_stat_path = os.path.join(self.utr3_path, "statistics")

    def _check_folder(self, folder):
        if folder is None:
            print("Error: lack required files!!!")
            sys.exit()

    def _check_gff(self, folder):
        for gff in os.listdir(folder):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(folder, gff))

    def _compute_utr(self, gffs, utr5_path, utr3_path, tss_path, tran_path,
                     source, terms, fuzzy, base_5utr):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                prefix = gff[:-4]
                tss = self.helper.get_correct_file(tss_path, "_TSS.gff",
                                                   prefix, None)
                tran = self.helper.get_correct_file(tran_path, "_transcript.gff",
                                                    prefix, None)
                if terms:
                    term = self.helper.get_correct_file(
                                os.path.join(terms, "tmp"), "_term.gff",
                                prefix, None)
                else:
                    term = None
                print("computing 5'UTR of {0} .....".format(prefix))
                detect_5utr(tss, os.path.join(gffs, gff), tran, source,
                            base_5utr, os.path.join(utr5_path, "gffs",
                            "_".join([prefix, "5UTR.gff"])))
                print("computing 3'UTR of {0} .....".format(prefix))
                detect_3utr(tran, os.path.join(gffs, gff), term, fuzzy,
                            os.path.join(utr3_path, "gffs",
                            "_".join([prefix, "3UTR.gff"])))
                self.helper.move_all_content(os.getcwd(),
                     self.utr5_stat_path, "_5utr_length.png")
                self.helper.move_all_content(os.getcwd(),
                     self.utr3_stat_path, "_3utr_length.png")

    def run_utr_detection(self, tsss, gffs, trans, terms,
                          fuzzy, out_folder, source, base_5utr):
        self._check_folder(tsss)
        self._check_folder(gffs)
        self._check_folder(trans)
        self._check_gff(tsss)
        self._check_gff(gffs)
        self._check_gff(trans)
        self._check_gff(terms)
        self.multiparser.parser_gff(gffs, None)
        self.multiparser.parser_gff(tsss, "TSS")
        self.multiparser.combine_gff(gffs, self.tss_path, None, "TSS")
        self.multiparser.parser_gff(trans, "transcript")
        self.multiparser.combine_gff(gffs, self.tran_path, None, "transcript")
        if terms:
            self.multiparser.parser_gff(terms, "term")
            self.multiparser.combine_gff(gffs, os.path.join(terms, "tmp"),
                                          None, "term")
        self._compute_utr(gffs, self.utr5_path, self.utr3_path, self.tss_path,
                          self.tran_path, source, terms, fuzzy, base_5utr)
        self.helper.remove_tmp(gffs)
        self.helper.remove_tmp(tsss)
        self.helper.remove_tmp(trans)
        self.helper.remove_tmp(terms)
        self.helper.remove_tmp(self.utr5_path)
        self.helper.remove_tmp(self.utr3_path)
