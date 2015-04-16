#!/usr/bin/python
from Bio import SeqIO
import os
import sys
from subprocess import call
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.detect_operon import operon
from annogesiclib.stat_operon import stat
from annogesiclib.combine_gff import combine_gff


class Operon_detection(object):

    def __init__(self):
        self.multiparser = Multiparser()
        self.helper = Helper()

    def _check_gff(self, gffs, type_):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(gffs, gff))

    def _detect_operon(self, prefixs, paths, tss_fuzzy, term_fuzzy, length):
        for prefix in prefixs:
            out_table = os.path.join(paths["table"], "_".join(["operon", prefix + ".csv"]))
            print("Detection operons of {0}".format(prefix))
            tss = self.helper.get_correct_file(paths["tss"], "_TSS.gff", prefix, None)
            tran = self.helper.get_correct_file(paths["tran"], "_transcript.gff", prefix, None)
            gff = self.helper.get_correct_file(paths["gff"], ".gff", prefix, None)
            if paths["term"] is None:
                term = False
            else:
                term = self.helper.get_correct_file(paths["term"], "_term.gff", prefix, None)
            operon(tran, tss, gff, term, tss_fuzzy, term_fuzzy, length, out_table)

    def _check_and_parser_gff(self, tsss, gffs, trans, utr5s, utr3s, output_folder,
                              terms):
#        self._check_gff(tsss, "tss")
#        self._check_gff(gffs, "gff")
#        self._check_gff(trans, "tran")
#        self._check_gff(utr5s, "utr")
#        self._check_gff(utr3s, "utr")
        tss_path = os.path.join(tsss, "tmp")
        tran_path = os.path.join(trans, "tmp")
        utr5_path = os.path.join(utr5s, "tmp")
        utr3_path = os.path.join(utr3s, "tmp")
        table_path = os.path.join(output_folder, "tables")
        if terms is not False:
#            self._check_gff(terms, "term")
            term_path = os.path.join(terms, "tmp")
        else:
            term_path = None
        self.multiparser._parser_gff(gffs, None)
        self.multiparser._parser_gff(tsss, "TSS")
        self.multiparser._combine_gff(gffs, tss_path, None, "TSS")
        self.multiparser._parser_gff(trans, "transcript")
        self.multiparser._combine_gff(gffs, tran_path, None, "transcript")
        self.multiparser._parser_gff(utr5s, "5UTR")
        self.multiparser._combine_gff(gffs, utr5_path, None, "5UTR")
        self.multiparser._parser_gff(utr3s, "3UTR")
        self.multiparser._combine_gff(gffs, utr3_path, None, "3UTR")
        self.multiparser._parser_gff(terms, "term")
        self.multiparser._combine_gff(gffs, term_path, None, "term")
        paths = {"gff": gffs, "tss": tss_path, "tran": tran_path, "utr5": utr5_path,
                 "utr3": utr3_path, "term": term_path, "table": table_path}
        return paths

    def _stat(self, table_path, stat_folder):
        for table in os.listdir(table_path):
            if table.startswith("operon_") and \
               table.endswith(".csv"):
                filename = "_".join(["stat", table])
                out_stat = os.path.join(stat_folder, filename)
                stat(os.path.join(table_path, table), out_stat)

    def _combine_gff(self, prefixs, output_folder, paths, tss_fuzzy, term_fuzzy):
        for prefix in prefixs:
            out_file = os.path.join(output_folder, "gffs", 
                                    "_".join([prefix, "all_features.gff"]))
            print("Combine all features of {0}".format(prefix))
            tss = self.helper.get_correct_file(paths["tss"], "_TSS.gff", prefix, None)
            tran = self.helper.get_correct_file(paths["tran"], "_transcript.gff", prefix, None)
            gff = self.helper.get_correct_file(paths["gff"], ".gff", prefix, None)
            utr5 = self.helper.get_correct_file(paths["utr5"], "_5UTR.gff", prefix, None)
            utr3 = self.helper.get_correct_file(paths["utr3"], "_3UTR.gff", prefix, None)
            if paths["term"] is None:
                term = None
            else:
                term = self.helper.get_correct_file(paths["term"], "_term.gff", prefix, None)
            combine_Gff(gff, tran, tss, utr5, utr3, term, tss_fuzzy, term_fuzzy, out_file)

    def run_operon(self, tsss, gffs, trans, utr5s, utr3s, terms, tss_fuzzy, term_fuzzy, 
                   length, stat, output_folder, combine_gff, stat_folder, bin_path):
        paths = self._check_and_parser_gff(tsss, gffs, trans, utr5s, utr3s, 
                                           output_folder, terms)
        prefixs = []
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                prefixs.append(gff.replace(".gff", ""))
        self._detect_operon(prefixs, paths, tss_fuzzy, term_fuzzy, length)
        if stat is not False:
            self._stat(paths["table"], stat_folder)
        if combine_gff is not False:
            self._combine_gff(prefixs, output_folder, paths, tss_fuzzy, term_fuzzy)
