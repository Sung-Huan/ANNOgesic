import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file
import annogesiclib.operon as op
from annogesiclib.operon import OperonDetection
from mock_args_container import MockClass


class Mock_func(object):

    def mock_operon(self, tran, tss, gff, term, tss_fuzzy,
                    term_fuzzy, length, out_table):
        gen_file(out_table, "test")

    def mock_stat(self, table, out_stat):
        gen_file(out_stat, "test")

    def mock_combine_gff(self, gff, tran, tss, utr5, utr3, term,
                        tss_fuzzy, term_fuzzy, out_file):
        gen_file(out_file, "test")

class TestOperonDetection(unittest.TestCase):

    def setUp(self):
        self.test_folder = "test_folder"
        self.mock_args = MockClass()
        self.mock = Mock_func()
        self.tsss = os.path.join(self.test_folder, "tsss")
        self.trans = os.path.join(self.test_folder, "trans")
        self.utr5s = os.path.join(self.test_folder, "utr5s")
        self.utr3s = os.path.join(self.test_folder, "utr3s")
        self.output = os.path.join(self.test_folder, "output")
        self.gffs = os.path.join(self.test_folder, "gffs")
        self.out_gff = os.path.join(self.output, "gffs")
        self.stat = os.path.join(self.test_folder, "stat")
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
            os.mkdir(self.gffs)
            os.mkdir(self.tsss)
            os.mkdir(self.stat)
            os.mkdir(os.path.join(self.tsss, "tmp"))
            os.mkdir(self.trans)
            os.mkdir(os.path.join(self.trans, "tmp"))
            os.mkdir(self.utr5s)
            os.mkdir(os.path.join(self.utr5s, "tmp"))
            os.mkdir(self.utr3s)
            os.mkdir(os.path.join(self.utr3s, "tmp"))
            os.mkdir(self.output)
            os.mkdir(self.out_gff)
            os.mkdir(os.path.join(self.output, "tables"))
        args = self.mock_args.mock()
        args.tsss = self.tsss
        args.trans = self.trans
        args.utr5s = self.utr5s
        args.utr3s = self.utr3s
        args.output_folder = self.output
        args.terms = None
        self.operon = OperonDetection(args)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_detect_operon(self):
        op.operon = self.mock.mock_operon
        gen_file(os.path.join(self.tsss, "tmp", "test_TSS.gff"), "test")
        gen_file(os.path.join(self.trans, "tmp", "test_transcript.gff"), "test")
        gen_file(os.path.join(self.gffs, "test.gff"), "test")
        args = self.mock_args.mock()
        args.gffs = self.out_gff
        args.term_fuzzy = 3
        args.tss_fuzzy = 3
        args.length = 100
        self.operon._detect_operon(["test"], args)
        self.assertTrue(os.path.exists(os.path.join(self.output, "tables",
                        "operon_test.csv")))
    
    def test_stat(self):
        op.stat = self.mock.mock_stat
        table_file = os.path.join(self.output, "tables", "operon_test.csv")
        if not os.path.exists(table_file):
            gen_file(table_file, "test")
        self.operon._stat(os.path.join(self.output, "tables"), self.stat)
        self.assertTrue(os.path.exists(os.path.join(self.stat, "stat_operon_test.csv")))

    def test_combine_gff(self):
        op.combine_gff = self.mock.mock_combine_gff
        gen_file(os.path.join(self.tsss, "tmp", "test_TSS.gff"), "test")
        gen_file(os.path.join(self.trans, "tmp", "test_transcript.gff"), "test")
        gen_file(os.path.join(self.gffs, "test.gff"), "test")
        gen_file(os.path.join(self.utr5s, "tmp", "test_5UTR.gff"), "test")
        gen_file(os.path.join(self.utr3s, "tmp", "test_3UTR.gff"), "test")
        args = self.mock_args.mock()
        args.gffs = self.out_gff
        args.term_fuzzy = 3
        args.tss_fuzzy = 3
        self.operon._combine_gff(["test"], args)
        self.assertTrue(os.path.exists(os.path.join(self.out_gff, "test_all_features.gff")))

if __name__ == "__main__":
    unittest.main()

