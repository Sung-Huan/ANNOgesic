import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
import annogesiclib.blast_class as blast_class

class Mock_func(object):

    def mock_read_file(self, blast_file, nums):
        nums['total']['dnaA'] = 2
        nums['aaa'] = {}
        nums['aaa']['dnaA'] = 2

class TestBlastClass(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
        self.blast_file = os.path.join(self.test_folder, "test.csv")
        with open(self.blast_file, "w") as rh:
            rh.write(self.example.blast)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_read_file(self):
        nums = {}
        nums["total"] = {}
        blast_class.read_file(self.blast_file, nums)
        self.assertDictEqual(nums, {'aaa': {'dnaA': 2}, 'total':{'dnaA': 2}})

    def test_blast_class(self):
        blast_class.read_file = Mock_func().mock_read_file
        lines = []
        out_file = os.path.join(self.test_folder, "test.out")
        blast_class.blast_class(self.blast_file, out_file)
        with open(out_file) as fh:
            for line in fh:
                line = line.strip()
                lines.append(line)
        self.assertEqual(set(lines), set(self.example.blast_table.split("\n")))

class Example(object):
    blast = """1\taaa\tdnaA\t2377296\t2377454\t-\tTSS:2377454_-\tNA\t2377296-2377454\tTEX+/-;Fragmented\t260123.91873361162\t446839.7634471806\t-0.0\tpMEM_t2_TEX_reverse(avg=155022.7050613754;high=266113.8349051722;low=0.6611942741842581)\t-0.2075\tIntergenic\tNA\tNA\t6\tNA\tsrn_4390|S._aureus_NCTC8325|dnaA|3e-55\tNA\tNA\tNA"""
    read_out = [{'ID': 'hit_1', 'srna_name': 'dnaA', 'blast_strain': 'strain_b', 'strain': 'aaa', 'start': '100', 'name': 'RNA_test1', 'strand': '+', 'end': '200', 'e': '0.0005'},
                {'ID': 'hit_2', 'srna_name': 'dnaa', 'blast_strain': 'strain_c', 'strain': 'aaa', 'start': '100', 'name': 'RNA_test1', 'strand': '+', 'end': '200', 'e': '0.0007'},
                {'ID': 'hit_3', 'srna_name': 'dnaC', 'blast_strain': 'strain_b', 'strain': 'aaa', 'start': '400', 'name': 'RNA_test2', 'strand': '-', 'end': '450', 'e': '0.000002'}]

    blast_table = """aaa:
sRNA_name	amount
dnaA	2"""

if __name__ == "__main__":
    unittest.main()

