import sys
import os
import unittest
import shutil
import copy
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data, gen_file, extract_info
import annogesiclib.stat_operon as so


class TestStatOpern(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_boolean(self):
        data = "True"
        result = so._boolean(data)
        self.assertTrue(result)
        data = "False"
        result = so._boolean(data)
        self.assertFalse(result)

    def test_row_to_location(self):
        row = self.example.operon
        data = so.row_to_location(row)
        self.assertDictEqual(data, {'start with tss': True,
                                    'stop with terminator': True,
                                    'have sub-operons': True,
                                    'have no sub-operons': False})

    def test_plus_num(self):
        num_total = {"total": {"total": 2, "mono": 0, "multi": 2},
                     "aaa": {"mono": 0, "multi": 2, "total": 2}}
        so.plus_num(num_total, "aaa", "mono")
        self.assertDictEqual(num_total, {
            'aaa': {'mono': 1, 'multi': 2, 'total': 3},
            'total': {'mono': 1, 'multi': 2, 'total': 3}})

    def test_print_stat(self):
        out = StringIO()
        class_operon = {"na": 20, "mono": 70, "poly": 10, "total": 100}
        operons = [{'start with tss': True,
                    'stop with terminator': True,
                    'have sub-operons': True,
                    'have no sub-operons': False}]
        so.print_stat(operons, 100, class_operon, out)
        na = out.getvalue().split("monocistronic: ")[-1].split(" (0.7)")[0]
        self.assertEqual(na, "70")
        na = out.getvalue().split("polycistronic: ")[-1].split(" (0.1)")[0]
        self.assertEqual(na, "10")

    def test_stat(self):
        input_file = os.path.join(self.test_folder, "input_file")
        output_file = os.path.join(self.test_folder, "output_file")
        gen_file(input_file, self.example.operon)
        so.stat(input_file, output_file)
        datas = import_data(output_file)
        for data in datas:
            if "the number of operons which have sub-operons and " in data:
                self.assertEqual(data,
                                 ("\tthe number of operons which have "
                                  "sub-operons and start with tss = 1 (1.0)"))
            if "the number of operons which have sub-operons = " in data:
                self.assertEqual(data,
                                 ("\tthe number of operons which have "
                                  "sub-operons = 1 (1.0)"))
            if "the number of operons which start with tss = " in data:
                self.assertEqual(data,
                                 ("\tthe number of operons which start "
                                  "with tss = 1 (1.0)"))
            if "no associated with CDS: " in data:
                self.assertEqual(data, "\tno associated with CDS: 0 (0.0)")
            if "monocistronic: " in data:
                self.assertEqual(data, "\tmonocistronic: 0 (0.0)")
            if "polycistronic: " in data:
                self.assertEqual(data, "\tpolycistronic: 1 (1.0)")
            

class Example(object):
    operon = """Operon1	Staphylococcus_aureus_HG003	313-3344	+	3	313-2127	True	4	False	0	1	2	SAOUHSC_00001	SAOUHSC_00001, SAOUHSC_00002"""
    test_print = """Total number of Operons is 100
The sub operon and features:
	the number of operons which start with tss = 1 (0.01)
	the number of operons which have sub-operons and stop with terminator = 1 (0.01)
	the number of operons which stop with terminator = 1 (0.01)
	the number of operons which start with tss and have sub-operons and stop with terminator = 1 (0.01)
	the number of operons which have sub-operons = 1 (0.01)
	the number of operons which start with tss and stop with terminator = 1 (0.01)
	the number of operons which start with tss and have sub-operons = 1 (0.01)
mono/polycistronic:
	no associated with CDS: 20 (20)
	monocistronic: 70 (0.7)
	polycistronic: 10 (0.1)"""

if __name__ == "__main__":
    unittest.main()

