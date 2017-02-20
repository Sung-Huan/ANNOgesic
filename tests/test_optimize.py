import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file, import_data
from mock_args_container import MockClass
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

    def remove_tmp_dir(self, folder):
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

    def mock_optimization(self, wig_path, fasta_file, gff_file, args, strain, manual, length):
        gen_file(os.path.join(args.output_folder, "test.csv"), "test")

class TestOptimizeTSS(unittest.TestCase):

    def setUp(self):
        self.mock_args = MockClass()
        self.test_folder = "test_folder"
        self.fastas = os.path.join(self.test_folder, "fasta")
        self.wigs = os.path.join(self.test_folder, "wigs")
        self.gffs = os.path.join(self.test_folder, "gffs")
        self.manuals = os.path.join(self.test_folder, "manuals")
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
            os.mkdir(self.fastas)
            os.mkdir(os.path.join(self.fastas, "tmp"))
            os.mkdir(self.wigs)
            os.mkdir(os.path.join(self.wigs, "tmp"))
            os.mkdir(self.gffs)
            os.mkdir(os.path.join(self.gffs, "tmp"))
            os.mkdir(self.manuals)
            os.mkdir(os.path.join(self.manuals, "tmp"))

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_optimize_tss(self):
        opt.Helper = Mock_helper
        opt.Multiparser = Mock_multiparser
        opt.optimization = Mock_func().mock_optimization
        gen_file(os.path.join(self.gffs, "tmp", "test.gff"), "test")
        gen_file(os.path.join(self.fastas, "tmp", "test.fa"), "test")
        args = self.mock_args.mock()
        args.fastas = self.fastas
        args.gffs = self.gffs
        args.wigs = self.wigs
        args.tsspredator_path = "test"
        args.manuals = self.manuals
        gen_file(os.path.join(self.manuals, "tmp", "test.gff"), "test")
        args.output_folder = self.test_folder
        args.project_strain = "test"
        args.height = 9
        args.height_reduction = 9
        args.factor = 9
        args.factor_reduction = 9
        args.base_height = 9
        args.enrichment = 9
        args.processing = 9
        args.utr = 200
        args.libs = "test"
        args.replicate_name = "test"
        args.cluster = 2
        args.strain_lengths = {"test": 100}
        args.cores = 4
        args.program = "TSS"
        args.replicate = 2
        args.steps = 2000
        opt.optimize_tss(args)
        self.assertTrue(os.path.exists(os.path.join(self.test_folder, "test.csv")))

if __name__ == "__main__":
    unittest.main()

