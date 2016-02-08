import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file, import_data, extract_info
import annogesiclib.sorf as so
from annogesiclib.sorf import sORFDetection

class Mock_func(object):

    def mock_sorf_detection(self, fasta, srna_file, gff, tss_file, utr_length,
                            utr_detect, libs, tex_notex, replicate,
                            cutoff_inter, cutoff_3utr, cutoff_5utr,
                            cutoff_intercds, test3, test2, merge_wigs,
                            start_codon, stop_codon, table_best, max_len,
                            min_len, test1, background,
                            fuzzy_rbs, print_all, no_srna, noafter_tss, no_tss,
                            min_rbs, max_rbs):
        pass

    def mock_get_intergenic(self, gff, tran, inter, utr_detect, hypo):
        pass

    def mock_check_gff(self, gffs):
        pass

    def mock_remove_tmp(self, out_folder, fastas, gffs, tsss,
                        trans, srnas, tex_wigs, frag_wigs):
        pass

    def mock_check_necessary_files(self, gffs, trans, tex_wigs,
                                   frag_wigs, utr_detect, tsss, srnas):
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
        self.sorf = sORFDetection(self.tsss, self.srnas, self.out,
                                  self.trans, self.fastas)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_merge_wigs(self):
        frag_wig = os.path.join(self.frag, "frag_1.wig")
        tex_wig = os.path.join(self.tex, "tex_1.wig")
        gen_file(frag_wig, "test")
        gen_file(tex_wig, "test")
        merge = self.sorf._merge_wigs(self.frag, self.tex)
        self.assertEqual(merge, "test_folder/merge_wigs")

    def test_combine_libs(self):
        lib = self.sorf._combine_libs("tlibs", "flibs")
        self.assertEqual(lib, "tlibsflibs")

    def test_start_stop_codon(self):
        gff_path = os.path.join(self.out, "gffs")
        table_path = os.path.join(self.out, "tables")
        os.mkdir(gff_path)
        os.mkdir(table_path)
        os.mkdir(os.path.join(gff_path, "all_candidates"))
        os.mkdir(os.path.join(table_path, "all_candidates"))
        os.mkdir(os.path.join(gff_path, "best"))
        os.mkdir(os.path.join(table_path, "best"))
        gen_file(os.path.join(gff_path, "all_candidates/test_sORF_all.gff"), "test")
        gen_file(os.path.join(gff_path, "all_candidates/test_sORF_all.csv"), "test")
        gen_file(os.path.join(gff_path, "all_candidates/test_sORF_best.gff"), "test")
        gen_file(os.path.join(gff_path, "all_candidates/test_sORF_best.csv"), "test")
        so.sorf_detection = self.mock.mock_sorf_detection
        self.sorf._start_stop_codon(["test"], self.out, 300, "libs",
                                    "tex_notex", "replicate", 10, 10,
                                    10, 10, "wig_path", "merge_wigs",
                                    ["ATG"], ["TTA"], 300, 30, True,
                                    True, "background", 2,
                                    False, True, True, True, 0, 20)
        self.assertTrue(os.path.exists(os.path.join(gff_path, "best/test_sORF.gff")))
        self.assertTrue(os.path.exists(os.path.join(gff_path, "all_candidates/test_sORF.gff")))
        self.assertTrue(os.path.exists(os.path.join(table_path, "best/test_sORF.csv")))
        self.assertTrue(os.path.exists(os.path.join(table_path, "all_candidates/test_sORF.csv")))        

    def test_compare_tran_cds(self):
        so.get_intergenic = self.mock.mock_get_intergenic
        gen_file(os.path.join(self.test_folder, "test.gff"), "test")
        prefixs = self.sorf._compare_tran_cds(self.test_folder, self.out, True, False)
        self.assertListEqual(prefixs, ["test"])

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
        self.sorf.run_sorf_detection(self.out, True, self.trans, self.gffs, self.tsss,
                                     300, 30, 300, self.tex, self.frag, 10,
                                     10, 10, 10, self.fastas, "tlibs", "flibs",
                                     "tex_notex", 2, 3, True, self.srnas,
                                     ["ATG"], ["TTA"], "background",
                                     2, True, False, True, True, False, 0, 20)


if __name__ == "__main__":
    unittest.main()
