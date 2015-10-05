import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file, import_data, extract_info
import annogesiclib.sublocal as su
from annogesiclib.sublocal import SubLocal

class Mock_func(object):

    def mock_psortb(self, psortb_path, strain_type, prot_seq_file,
                   out_raw, out_err):
        pass

    def mock_stat_sublocal(self, merge_table, stat, table):
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


class TestSubLocal(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.mock = Mock_func()
        self.test_folder = "test_folder"
        self.out = "test_folder/output"
        self.fastas = "test_folder/fastas"
        self.gffs = "test_folder/gffs"
        self.stat = "test_folder/stat"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
            os.mkdir(self.out)
            os.mkdir(self.fastas)
            os.mkdir(self.gffs)
            os.mkdir(self.stat)
        self.sub = SubLocal(self.gffs, self.fastas, self.out)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_get_protein_seq(self):
        gen_file(os.path.join(self.fastas, "aaa.fa"), self.example.fasta_file)
        gff = "aaa.gff"
        gen_file(os.path.join(self.gffs, gff), self.example.gff_file)
        prefix = self.sub._get_protein_seq(self.fastas, gff, self.test_folder, self.gffs)
        self.assertEqual(prefix, "aaa")

    def test_run_psortb(self):
        self.sub._psortb = self.mock.mock_psortb
        tmp_result = os.path.join(self.out, "tmp_results")
        os.mkdir(tmp_result)
        self.sub._run_psortb("aaa", self.out, "positive", "psortb_path", self.test_folder)
        self.assertTrue(os.path.exists(os.path.join(self.out, "tmp_log")))
        self.assertTrue(os.path.exists(os.path.join(tmp_result,
                       "_".join(["aaa", "raw.txt"]))))

    def test_merge_and_stat(self):
        su.stat_sublocal = self.mock.mock_stat_sublocal
        os.mkdir(os.path.join(self.gffs, "aaa.gff_folder"))
        gen_file(os.path.join(self.gffs, "aaa.gff_folder/aaa.gff"), "test")
        os.mkdir(os.path.join(self.out, "psortb_results"))
        gen_file(os.path.join(self.test_folder, "aaa_raw.txt"), "test")
        gen_file(os.path.join(self.test_folder, "aaa_table.csv"), "test")
        self.sub._merge_and_stat(self.gffs, self.out, self.test_folder, self.stat)
        self.assertTrue(os.path.exists(os.path.join(self.stat, "aaa")))
        self.assertTrue(os.path.exists(os.path.join(self.out, "psortb_results/aaa")))


class Example(object):

    fasta_file = """aaa
AATGTGTATACCCACAGTCTCAGCTACACATACAGTT"""
    gff_file = """aaa	RefSeq	CDS	3	17	.	+	.	ID=cds0;Name=CDS_0;locus_tag=AAA_00001"""

if __name__ == "__main__":
    unittest.main()
