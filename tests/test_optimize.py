import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file, import_data
import annogesiclib.optimize as opt


class Mock_helper(object):

    def __init__(self):
        pass

    def check_uni_attributes(self, gff_file):
        pass

    def remove_all_content(self, filename, feature, type_):
        pass

    def remove_tmp(self, wigs):
        pass


class Mock_multiparser(object):

    def __init__(self):
        pass
    
    def parser_wig(self, wigs):
        pass

    def parser_gff(self, gffs, feature):
        pass

    def parser_fasta(self, fastas):
        pass

class Mock_func(object):

    def mock_optimization(self, tsspredator_path, max_height, max_height_reduction,
                          max_factor, max_factor_reduction, max_base_height,
                          max_enrichment, max_processing, output_folder, core,
                          wig_path, project_name, fasta_file, replicate_name, steps,
                          gff_file, program, manual, libs, length, cluster,
                          utr_length, replicate):
        gen_file(os.path.join(output_folder, "test.csv"), "test")

class TestOptimizeTSS(unittest.TestCase):

    def setUp(self):
        self.test_folder = "test_folder"
        self.fastas = os.path.join(self.test_folder, "fasta")
        self.wigs = os.path.join(self.test_folder, "wigs")
        self.gffs = os.path.join(self.test_folder, "gffs")
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
            os.mkdir(self.fastas)
            os.mkdir(os.path.join(self.fastas, "tmp"))
            os.mkdir(self.wigs)
            os.mkdir(os.path.join(self.wigs, "tmp"))
            os.mkdir(self.gffs)
            os.mkdir(os.path.join(self.gffs, "tmp"))

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_optimize_tss(self):
        opt.Helper = Mock_helper
        opt.Multiparser = Mock_multiparser
        opt.optimization = Mock_func().mock_optimization
        gen_file(os.path.join(self.gffs, "tmp", "test.gff"), "test")
        gen_file(os.path.join(self.fastas, "tmp", "test.fa"), "test")
        opt.optimize_tss("test", self.fastas, self.gffs, self.wigs, "test", self.test_folder,
                         "test", 9, 9, 9, 9, 9, 9, 9, 200, "test", "test", 3,
                         100, 4, "TSS", "test", 5000)
        self.assertTrue(os.path.exists(os.path.join(self.test_folder, "test.csv")))

if __name__ == "__main__":
    unittest.main()

