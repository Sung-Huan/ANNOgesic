import os
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.detect_operon import operon
from annogesiclib.stat_operon import stat
from annogesiclib.combine_gff import combine_gff


class OperonDetection(object):

    def __init__(self, tsss, trans, utr5s, utr3s, output_folder, terms):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.tss_path = os.path.join(tsss, "tmp")
        self.tran_path = os.path.join(trans, "tmp")
        self.utr5_path = os.path.join(utr5s, "tmp")
        self.utr3_path = os.path.join(utr3s, "tmp")
        self.table_path = os.path.join(output_folder, "tables")
        if terms is not None:
            self._check_gff(terms, "term")
            self.term_path = os.path.join(terms, "tmp")
        else:
            self.term_path = None

    def _check_gff(self, gffs, type_):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(gffs, gff))

    def _detect_operon(self, prefixs, gffs, tss_fuzzy, term_fuzzy, length):
        for prefix in prefixs:
            out_table = os.path.join(self.table_path,
                        "_".join(["operon", prefix + ".csv"]))
            print("Detection operons of {0}".format(prefix))
            tss = self.helper.get_correct_file(self.tss_path, "_TSS.gff",
                                               prefix, None)
            tran = self.helper.get_correct_file(self.tran_path,
                               "_transcript.gff", prefix, None)
            gff = self.helper.get_correct_file(gffs, ".gff", prefix, None)
            if self.term_path is None:
                term = False
            else:
                term = self.helper.get_correct_file(self.term_path, "_term.gff",
                                                    prefix, None)
            operon(tran, tss, gff, term, tss_fuzzy,
                   term_fuzzy, length, out_table)

    def _check_and_parser_gff(self, tsss, gffs, trans, utr5s, utr3s, terms):
        self._check_gff(tsss, "tss")
        self._check_gff(gffs, "gff")
        self._check_gff(trans, "tran")
        self._check_gff(utr5s, "utr")
        self._check_gff(utr3s, "utr")
        if terms is not False:
            self._check_gff(terms, "term")
        self.multiparser.parser_gff(gffs, None)
        self.multiparser.parser_gff(tsss, "TSS")
        self.multiparser.combine_gff(gffs, self.tss_path, None, "TSS")
        self.multiparser.parser_gff(trans, "transcript")
        self.multiparser.combine_gff(gffs, self.tran_path, None, "transcript")
        self.multiparser.parser_gff(utr5s, "5UTR")
        self.multiparser.combine_gff(gffs, self.utr5_path, None, "5UTR")
        self.multiparser.parser_gff(utr3s, "3UTR")
        self.multiparser.combine_gff(gffs, self.utr3_path, None, "3UTR")
        if terms is not None:
            self._check_gff(terms, "term")
            self.multiparser.parser_gff(terms, "term")
            self.multiparser.combine_gff(gffs, self.term_path, None, "term")

    def _stat(self, table_path, stat_folder):
        for table in os.listdir(table_path):
            if table.startswith("operon_") and table.endswith(".csv"):
                filename = "_".join(["stat", table])
                out_stat = os.path.join(stat_folder, filename)
                stat(os.path.join(table_path, table), out_stat)

    def _combine_gff(self, prefixs, output_folder, gffs, tss_fuzzy, term_fuzzy):
        for prefix in prefixs:
            out_file = os.path.join(output_folder, "gffs",
                                    "_".join([prefix, "all_features.gff"]))
            print("Combine all features of {0}".format(prefix))
            tss = self.helper.get_correct_file(self.tss_path, "_TSS.gff",
                                               prefix, None)
            tran = self.helper.get_correct_file(self.tran_path,
                               "_transcript.gff", prefix, None)
            gff = self.helper.get_correct_file(gffs, ".gff", prefix, None)
            utr5 = self.helper.get_correct_file(self.utr5_path,
                                                "_5UTR.gff", prefix, None)
            utr3 = self.helper.get_correct_file(self.utr3_path,
                                                "_3UTR.gff", prefix, None)
            if self.term_path is None:
                term = None
            else:
                term = self.helper.get_correct_file(self.term_path,
                                                    "_term.gff", prefix, None)
            combine_gff(gff, tran, tss, utr5, utr3, term,
                        tss_fuzzy, term_fuzzy, out_file)

    def run_operon(self, tsss, gffs, trans, utr5s, utr3s, terms,
                   tss_fuzzy, term_fuzzy, length, statistics, output_folder,
                   combine, stat_folder):
        self._check_and_parser_gff(tsss, gffs, trans, utr5s, utr3s, terms)
        prefixs = []
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                prefixs.append(gff.replace(".gff", ""))
        self._detect_operon(prefixs, gffs, tss_fuzzy, term_fuzzy, length)
        if statistics:
            self._stat(self.table_path, stat_folder)
        if combine:
            self._combine_gff(prefixs, output_folder,
                              gffs, tss_fuzzy, term_fuzzy)
