import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file, import_data, extract_info
import annogesiclib.sorf as so
from annogesiclib.sorf import sORFDetection
from mock_args_container import MockClass


class Mock_func(object):

    def mock_sorf_detection(self, fasta, srna_file, gff, tss_file, utr_length,
                            utr_detect, test3, test2):
        pass

    def mock_get_intergenic(self, gff, tran, inter, utr_detect, hypo, extend_5, extend_3):
        pass

    def mock_check_gff(self, gffs):
        pass

    def mock_remove_tmp(self, out_folder):
        pass

    def mock_check_necessary_files(self, gffs, log):
        pass
    

class Mock_Multiparser(object):

    def parser_wig(self, merge_wigs):
        pass

    def combine_wig(self, gffs, wig_path, test, libs):
        pass

    def parser_gff(self, trans, type_):
        pass

    def combine_gff(self, gffs, tran_path, test, type_):
        pass

    def parser_fasta(self, fastas):
        pass

    def combine_fasta(self, gffs, fasta_path, test):
        pass


class TestsORFDetection(unittest.TestCase):

    def setUp(self):
        self.mock = Mock_func()
        self.mock_args = MockClass()
        self.test_folder = "test_folder"
        self.tsss = "test_folder/tsss"
        self.srnas = "test_folder/sRNA"
        self.out = "test_folder/output"
        self.trans = "test_folder/trans"
        self.fastas = "test_folder/fastas"
        self.tex = "test_folder/tex"
        self.frag = "test_folder/frag"
        self.gffs = "test_folder/gffs"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
            os.mkdir(self.tsss)
            os.mkdir(self.out)
            os.mkdir(self.trans)
            os.mkdir(self.fastas)
            os.mkdir(self.tex)
            os.mkdir(self.frag)
            os.mkdir(self.srnas)
            os.mkdir(self.gffs)
        args = self.mock_args.mock()
        args.tsss = self.tsss
        args.srnas = self.srnas
        args.out_folder = self.out
        args.trans = self.trans
        args.fastas = self.fastas
        self.sorf = sORFDetection(args)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_start_stop_codon(self):
        gff_path = os.path.join(self.out, "gffs")
        table_path = os.path.join(self.out, "tables")
        os.mkdir(gff_path)
        os.mkdir(table_path)
        os.mkdir(os.path.join(gff_path, "all_candidates"))
        os.mkdir(os.path.join(table_path, "all_candidates"))
        os.mkdir(os.path.join(gff_path, "best_candidates"))
        os.mkdir(os.path.join(table_path, "best_candidates"))
        gen_file(os.path.join(gff_path, "all_candidates/test_sORF_all.gff"),
                 "test")
        gen_file(os.path.join(gff_path, "all_candidates/test_sORF_all.csv"),
                 "test")
        gen_file(os.path.join(gff_path, "all_candidates/test_sORF_best.gff"),
                 "test")
        gen_file(os.path.join(gff_path, "all_candidates/test_sORF_best.csv"),
                 "test")
        so.sorf_detection = self.mock.mock_sorf_detection
        args = self.mock_args.mock()
        args.libs = "libs"
        args.tex_notex = "tex_notex"
        args.replicates = "replicates"
        args.start_codon = ["ATG"]
        args.stop_codon = ["TTA"]
        args.background = "background"
        args.wig_path = "wig_path"
        args.merge_wigs = "merge_wigs"
        log = open(os.path.join(self.test_folder, "test.log"), "w")
        self.sorf._start_stop_codon(["test"], args, log)
        self.assertTrue(os.path.exists(os.path.join(
            gff_path, "best_candidates/test_sORF.gff")))
        self.assertTrue(os.path.exists(os.path.join(
            gff_path, "all_candidates/test_sORF.gff")))
        self.assertTrue(os.path.exists(os.path.join(
            table_path, "best_candidates/test_sORF.csv")))
        self.assertTrue(os.path.exists(os.path.join(
            table_path, "all_candidates/test_sORF.csv")))        
        log.close()

    def test_compare_tran_cds(self):
        so.get_intergenic = self.mock.mock_get_intergenic
        gen_file(os.path.join(self.test_folder, "test.gff"), "test")
        args = self.mock_args.mock()
        args.out_folder = self.out
        args.gffs = self.test_folder
        args.hypo = False
        args.utr_detect = True
        args.extend_5 = 5
        args.extend_3 = 75
        log = open(os.path.join(self.test_folder, "test.log"), "w")
        prefixs = self.sorf._compare_tran_cds(args, log)
        self.assertListEqual(prefixs, ["test"])
        log.close()

    def test_run_sorf_detection(self):
        gff_path = os.path.join(self.out, "gffs")
        table_path = os.path.join(self.out, "tables")
        os.mkdir(gff_path)
        os.mkdir(table_path)
        os.mkdir(os.path.join(gff_path, "all_candidates"))
        os.mkdir(os.path.join(table_path, "all_candidates"))
        os.mkdir(os.path.join(gff_path, "best"))
        os.mkdir(os.path.join(table_path, "best"))
        so.get_intergenic = self.mock.mock_get_intergenic
        so.sorf_detection = self.mock.mock_sorf_detection
        self.sorf._remove_tmp = self.mock.mock_remove_tmp
        self.sorf._check_gff = self.mock.mock_check_gff
        self.sorf._check_necessary_files = self.mock.mock_check_necessary_files
        self.sorf.multiparser = Mock_Multiparser()
        args = self.mock_args.mock()
        args.trans = self.trans
        args.gffs = self.gffs
        args.tsss = self.tsss
        args.out_folder = self.out
        args.libs = "libs"
        args.tex_notex = "tex_notex"
        args.replicates = "replicates"
        args.start_codon = ["ATG"]
        args.stop_codon = ["TTA"]
        args.background = "background"
        args.wig_path = "wig_path"
        args.merge_wigs = "merge_wigs"
        args.fuzzy_rbs = 2
        log = open(os.path.join(self.test_folder, "test.log"), "w")
        self.sorf.run_sorf_detection(args, log)
        log.close()

if __name__ == "__main__":
    unittest.main()
