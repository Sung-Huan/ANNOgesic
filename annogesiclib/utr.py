import os
import sys
from annogesiclib.helper import Helper
from annogesiclib.detect_utr import detect_3utr, detect_5utr
from annogesiclib.multiparser import Multiparser


class UTRDetection(object):
    '''detection of UTR'''

    def __init__(self, args_utr):
        self.helper = Helper()
        self.multiparser = Multiparser()
        self.tss_path = os.path.join(args_utr.tsss, "tmp")
        self.tran_path = os.path.join(args_utr.trans, "tmp")
        self.utr5_path = os.path.join(args_utr.out_folder, "5UTRs")
        self.utr3_path = os.path.join(args_utr.out_folder, "3UTRs")
        self.utr5_stat_path = os.path.join(self.utr5_path, "statistics")
        self.utr3_stat_path = os.path.join(self.utr3_path, "statistics")

    def _check_folder(self, folder):
        if folder is None:
            print("Error: Lack required files!")
            sys.exit()

    def _check_gff(self, folder):
        for gff in os.listdir(folder):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(folder, gff))

    def _compute_utr(self, args_utr, log):
        log.write("Running detect_utr.py to detect UTRs.\n")
        for gff in os.listdir(args_utr.gffs):
            if gff.endswith(".gff"):
                prefix = gff[:-4]
                tss = self.helper.get_correct_file(
                        self.tss_path, "_TSS.gff", prefix, None, None)
                tran = self.helper.get_correct_file(
                        self.tran_path, "_transcript.gff", prefix, None, None)
                if args_utr.terms:
                    term = self.helper.get_correct_file(
                                os.path.join(args_utr.terms, "tmp"),
                                "_term.gff", prefix, None, None)
                else:
                    term = None
                print("Computing 5'UTRs of {0}".format(prefix))
                detect_5utr(tss, os.path.join(args_utr.gffs, gff),
                            tran, os.path.join(self.utr5_path, "gffs",
                            "_".join([prefix, "5UTR.gff"])), args_utr)
                print("Computing 3'UTRs of {0}".format(prefix))
                detect_3utr(tran, os.path.join(args_utr.gffs, gff),
                            term, os.path.join(self.utr3_path, "gffs",
                            "_".join([prefix, "3UTR.gff"])), args_utr)
                self.helper.move_all_content(
                    os.getcwd(), self.utr5_stat_path, ["_5utr_length.png"])
                self.helper.move_all_content(
                    os.getcwd(), self.utr3_stat_path, ["_3utr_length.png"])
        log.write("The following files are generated:\n")
        for folder in (os.path.join(self.utr5_path, "gffs"),
                       os.path.join(self.utr3_path, "gffs"),
                       self.utr5_stat_path, self.utr3_stat_path):
            for file_ in os.listdir(folder):
                log.write("\t" + os.path.join(folder, file_) + "\n")

    def run_utr_detection(self, args_utr, log):
        self._check_folder(args_utr.tsss)
        self._check_folder(args_utr.gffs)
        self._check_folder(args_utr.trans)
        self._check_gff(args_utr.tsss)
        self._check_gff(args_utr.gffs)
        self._check_gff(args_utr.trans)
        self._check_gff(args_utr.terms)
        self.multiparser.parser_gff(args_utr.gffs, None)
        self.multiparser.parser_gff(args_utr.tsss, "TSS")
        self.multiparser.combine_gff(args_utr.gffs, self.tss_path, None, "TSS")
        self.multiparser.parser_gff(args_utr.trans, "transcript")
        self.multiparser.combine_gff(args_utr.gffs, self.tran_path,
                                     None, "transcript")
        if args_utr.terms:
            self.multiparser.parser_gff(args_utr.terms, "term")
            self.multiparser.combine_gff(args_utr.gffs,
                                         os.path.join(args_utr.terms, "tmp"),
                                         None, "term")
        self._compute_utr(args_utr, log)
        self.helper.remove_tmp_dir(args_utr.gffs)
        self.helper.remove_tmp_dir(args_utr.tsss)
        self.helper.remove_tmp_dir(args_utr.trans)
        self.helper.remove_tmp_dir(args_utr.terms)
        self.helper.remove_tmp(self.utr5_path)
        self.helper.remove_tmp(self.utr3_path)
