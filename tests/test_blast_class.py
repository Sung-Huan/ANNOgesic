import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
import annogesiclib.blast_class as blast_class

class Mock_func(object):

    def mock_read_file(self, blast_file):
        repeats = {'RNA_test1': ['hit_1|strain_b|dnaA', 'hit_2|strain_c|dnaa']}
        return Example().read_out, repeats


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
        datas, repeats = blast_class.read_file(self.blast_file)
        for index in range(0, 3):
            self.assertDictEqual(datas[index], self.example.read_out[index])
        self.assertDictEqual(repeats, {'RNA_test1': ['hit_1|strain_b|dnaA', 'hit_2|strain_c|dnaa']})

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
    blast = """aaa	RNA_test1	+	100	200	hit_1|strain_b|dnaA	0.0005
aaa	RNA_test1	+	100	200	hit_2|strain_c|dnaa	0.0007
aaa	RNA_test2	-	400	450	hit_3|strain_b|dnaC	0.000002"""
    read_out = [{'ID': 'hit_1', 'srna_name': 'dnaA', 'blast_strain': 'strain_b', 'strain': 'aaa', 'start': '100', 'name': 'RNA_test1', 'strand': '+', 'end': '200', 'e': '0.0005'},
                {'ID': 'hit_2', 'srna_name': 'dnaa', 'blast_strain': 'strain_c', 'strain': 'aaa', 'start': '100', 'name': 'RNA_test1', 'strand': '+', 'end': '200', 'e': '0.0007'},
                {'ID': 'hit_3', 'srna_name': 'dnaC', 'blast_strain': 'strain_b', 'strain': 'aaa', 'start': '400', 'name': 'RNA_test2', 'strand': '-', 'end': '450', 'e': '0.000002'}]

    blast_table = """All strain:
sRNA_name	amount
dnaA	2
dnaC	1

repeat counting:
RNA_test1:hit_1|strain_b|dnaA;hit_2|strain_c|dnaa"""

if __name__ == "__main__":
    unittest.main()

