import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import gen_file, import_data
import annogesiclib.change_db_format as cdf


class TestChangeDBFormat(unittest.TestCase):

    def setUp(self):
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_change_format(self):
        input_file = os.path.join(self.test_folder, "input")
        output_file = os.path.join(self.test_folder, "output")
        gen_file(input_file, ">srna_1|Staphylococcus|Aar|12314|12444|forward\nATAGATTCCCGCGTATAGTCATCATTGTAC")
        cdf.change_format(input_file, output_file)
        data = import_data(output_file)
        self.assertListEqual(data, ['>srna_1|Staphylococcus|Aar', 'ATAGATTCCCGCGTATAGTCATCATTGTAC'])

if __name__ == "__main__":
    unittest.main()

