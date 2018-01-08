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
        result = """#ID\tGenome\tStrand\tAssociated_CDS\tStart_genome\tEnd_genome\tRfam\tE_value\tScore\tStart_align\tEnd_align
riboswitch_5\tStaphylococcus_aureus_HG003\t+\tSAOUHSC_00013\t15948\t16046\tRF00162\t1.6e-18\t74\t1\t99
riboswitch_11\tStaphylococcus_aureus_HG003\t-\tSAOUHSC_00007\t27955\t28053\tRF00162\t1.6e-18\t74\t1\t99
riboswitch_183\tStaphylococcus_aureus_HG003\t+\tSAOUHSC_00372\t377996\t378098\tRF00167\t2.2e-18\t45\t1\t103"""
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

    ribos = """riboswitch_5\tStaphylococcus_aureus_HG003\t+\tSAOUHSC_00013\t15948\t16046\tRF00162\t1.6e-18\t74\t1\t99
riboswitch_11\tStaphylococcus_aureus_HG003\t-\tSAOUHSC_00007\t27955\t28053\tRF00162\t1.6e-18\t74\t1\t99
riboswitch_183\tStaphylococcus_aureus_HG003\t+\tSAOUHSC_00372\t377996\t378098\tRF00167\t2.2e-18\t45\t1\t103"""

if __name__ == "__main__":
    unittest.main()
