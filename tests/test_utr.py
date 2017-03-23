import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file, import_data, extract_info
import annogesiclib.utr as ut
from annogesiclib.utr import UTRDetection
from mock_args_container import MockClass


class Mock_func(object):

    def __init__(self):
        self.example = Example()

    def mock_detect_5utr(self, tss, gff, tran, source, utr5):
        gen_file("test_5utr_length.png", "test")

    def mock_detect_3utr(self, tran, gff, term, fuzzy, utr3):
        gen_file("test_3utr_length.png", "test")

    def mock_check_gff(self, folder):
        pass


class Mock_Multiparser(object):

    def parser_wig(self, merge_wigs):
        pass

    def combine_wig(self, gffs, wig_path, test):
        pass

    def parser_gff(self, trans, type_):
        pass

    def combine_gff(self, gffs, tran_path, test, type_):
        pass

    def parser_fasta(self, fastas):
        pass

    def combine_fasta(self, gffs, fasta_path, test):
        pass


class TestsTSSpredator(unittest.TestCase):

    def setUp(self):
        self.mock_args = MockClass()
        self.mock = Mock_func()
        self.mock_parser = Mock_Multiparser()
        self.example = Example()
        self.test_folder = "test_folder"
        self.trans = "test_folder/trans"
        self.out = "test_folder/output"
        self.gffs = "test_folder/gffs"
        self.tsss = "test_folder/tsss"
        self.terms = "test_folder/terms"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
            os.mkdir(self.trans)
            os.mkdir(os.path.join(self.trans, "tmp"))
            os.mkdir(self.out)
            os.mkdir(self.gffs)
            os.mkdir(self.tsss)
            os.mkdir(os.path.join(self.tsss, "tmp"))
            os.mkdir(self.terms)
        args = self.mock_args.mock()
        args.tsss = self.tsss
        args.trans = self.trans
        args.out_folder = self.out
        self.utr = UTRDetection(args)

#    def tearDown(self):
#        if os.path.exists(self.test_folder):
#            shutil.rmtree(self.test_folder)

    def test_compute_utr(self):
        ut.detect_5utr = self.mock.mock_detect_5utr
        ut.detect_3utr = self.mock.mock_detect_3utr
        term_path = os.path.join(self.terms, "tmp")
        os.mkdir(term_path)
        utr5_path = os.path.join(self.out, "5UTRs")
        utr3_path = os.path.join(self.out, "3UTRs")
        os.mkdir(utr5_path)
        os.mkdir(utr3_path)
        utr5_stat_path = os.path.join(utr5_path, "statistics")
        utr3_stat_path = os.path.join(utr3_path, "statistics")
        os.mkdir(utr5_stat_path)
        os.mkdir(utr3_stat_path)
        gen_file(os.path.join(self.gffs, "test.gff"), self.example.gff_file)
        gen_file(os.path.join(self.trans, "test_transcript.gff"),
                 self.example.tran_file)
        gen_file(os.path.join(self.tsss, "test_TSS.gff"),
                 self.example.tss_file)
        gen_file(os.path.join(term_path, "test_term.gff"),
                 self.example.term_file)
        args = self.mock_args.mock()
        args.gffs = self.gffs
        args.tsss = self.tsss
        args.trans = self.trans
        args.terms = self.terms
        self.utr._compute_utr(args)
        self.assertTrue(os.path.exists(os.path.join(
            utr5_stat_path, "test_5utr_length.png")))
        self.assertTrue(os.path.exists(os.path.join(
            utr3_stat_path, "test_3utr_length.png")))
        shutil.rmtree(utr5_path)
        shutil.rmtree(utr3_path)

    def test_run_utr_detection(self):
        self.utr._check_gff = self.mock.mock_check_gff
        ut.detect_5utr = self.mock.mock_detect_5utr
        ut.detect_3utr = self.mock.mock_detect_3utr
        utr5_path = os.path.join(self.out, "5UTRs")
        utr3_path = os.path.join(self.out, "3UTRs")
        os.mkdir(utr5_path)
        os.mkdir(utr3_path)
        utr5_stat_path = os.path.join(utr5_path, "statistics")
        utr3_stat_path = os.path.join(utr3_path, "statistics")
        os.mkdir(utr5_stat_path)
        os.mkdir(utr3_stat_path)
        gen_file(os.path.join(self.gffs, "test.gff"), self.example.gff_file)
        gen_file(os.path.join(self.trans, "test_transcript.gff"),
                 self.example.tran_file)
        gen_file(os.path.join(self.tsss, "test_TSS.gff"),
                 self.example.tss_file)
        gen_file(os.path.join(self.terms, "test_term.gff"),
                 self.example.term_file)
        args = self.mock_args.mock()
        args.tsss = self.tsss
        args.gffs = self.gffs
        args.trans = self.trans
        args.terms = self.terms
        args.out_folder = self.out
        self.utr.run_utr_detection(args)
        self.assertTrue(os.path.exists(os.path.join(
            utr5_stat_path, "test_5utr_length.png")))
        self.assertTrue(os.path.exists(os.path.join(
            utr3_stat_path, "test_3utr_length.png")))


class Example(object):

    tran_file = """test	Transcript	Transcript	3	25	.	+	.	ID=tran0;Name=Transcript_0"""
    gff_file = """test	RefSeq	CDS	5	10	.	+	.	ID=cds0;Name=CDS_0"""
    tss_file = """test	RefSeq	TSS	3	3	.	+	.	ID=tss0;Name=TSS_0"""
    term_file = """test	RefSeq	terminator	12	16 	.	+	.	ID=term0;Name=Term_0"""

if __name__ == "__main__":
    unittest.main()
