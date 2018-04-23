import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file, import_data
from annogesiclib.seq_editer import SeqEditer


class TestSeqEditer(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        self.fasta = os.path.join(self.test_folder, "fasta")
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
            os.mkdir(self.fasta)
        self.seq = SeqEditer()

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_import_data(self):
        mod_table = os.path.join(self.test_folder, "mod")
        gen_file(mod_table, self.example.mutation)
        datas = self.seq._import_data(mod_table, "test")
        self.assertListEqual(datas, [{'target_id': 'test_NC_000915.1',
             'datas': [{'ref_nt': 'c', 'tar_nt': '', 'position': '3'},
             {'ref_nt': '-', 'tar_nt': 'deletion', 'position': '6'}], 'ref_id': 'NC_000915.1'}])

    def test_modify_seq(self):
        mod_table = os.path.join(self.test_folder, "mod")
        gen_file(mod_table, self.example.mutation)
        gen_file(os.path.join(self.fasta, "NC_000915.1.fa"),
                 self.example.fasta)
        self.seq.modify_seq(self.fasta, mod_table, self.test_folder, "test")
        datas = import_data(os.path.join(self.test_folder, "test_NC_000915.1.fa"))
        self.assertEqual("\n".join(datas), self.example.out_1)

    def test_modify_header(self):
        input_file = os.path.join(self.test_folder, "test.fa")
        gen_file(input_file, ">AAA|BBB|CCC|DDD|EEE\nACATACAAGTACAGTT")
        self.seq.modify_header(input_file)
        datas = import_data(input_file)
        self.assertEqual("\n".join(datas), ">DDD\nACATACAAGTACAGTT")


class Example(object):

    fasta = """>NC_000915.1
ATAGATAACCCAAGTACGACTCAGGTCCCTCACA"""
    out_1 = """>test_NC_000915.1
ATGATdeletionAACCCAAGTACGACTCAGGTCCCTCACA"""
    out_2 = """>test_case2
ATAGAgTAACCCAAGTACGACTCAGGTCCCTCACA"""
    mutation = """#refernce_id	target_id	reference_nt	position	target_nt	impact of correction	locus tag	gene	Description
NC_000915.1	3	a	c		SAOUHSC_00002	dnaA	XXXXXX
NC_000915.1	6	a	-	deletion			YYYYYY"""

if __name__ == "__main__":
    unittest.main()
