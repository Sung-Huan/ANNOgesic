import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file, import_data
import annogesiclib.modify_rbs_table as mrt


class TestGenSvg(unittest.TestCase):

    def setUp(self):
        self.test_folder = "test_folder"
        self.example = Example()
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_modify_table(self):
        result = """#ID|strain|strand|associated_CDS|start|end	Rfam	e_value	start_align	end_align
riboswitch_5|Staphylococcus_aureus_HG003|+|SAOUHSC_00013|15948|16046	RF00162	1.6e-18	1	99
riboswitch_11|Staphylococcus_aureus_HG003|-|SAOUHSC_00007|27955|28053	RF00162	1.6e-18	1	99
riboswitch_183|Staphylococcus_aureus_HG003|+|SAOUHSC_00372|377996|378098	RF00167	2.2e-18	1	103"""
        table = os.path.join(self.test_folder, "test")
        gen_file(table, self.example.ribos)
        mrt.modify_table(table, True)
        data = import_data(table)
        self.assertEqual("\n".join(data), result)
        gen_file(table, self.example.ribos)
        mrt.modify_table(table, False)
        data = import_data(table)
        self.assertEqual("\n".join(data), result)

class Example(object):

    ribos = """riboswitch_5|Staphylococcus_aureus_HG003|+|SAOUHSC_00013|15948|16046	RF00162	1.6e-18	1	99
riboswitch_11|Staphylococcus_aureus_HG003|-|SAOUHSC_00007|27955|28053	RF00162	1.6e-18	1	99
riboswitch_183|Staphylococcus_aureus_HG003|+|SAOUHSC_00372|377996|378098	RF00167	2.2e-18	1	103"""

if __name__ == "__main__":
    unittest.main()
