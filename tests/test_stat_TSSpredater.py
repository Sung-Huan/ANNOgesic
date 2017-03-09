import sys
import os
import unittest
import shutil
import copy
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data, gen_file, extract_info
import annogesiclib.stat_TSSpredator as st


class TestStatTSSpredator(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_stat(self):
        detect = False
        out_stat = StringIO()
        out_lib = StringIO()
        st.stat(self.example.tsss, "aaa", "TSS", out_stat, "test", out_lib)
        for data in out_stat.getvalue().split("\n"):
            if "Primary" in data:
                self.assertEqual(data.split(" = ")[-1], "1 (1.0)")
        if ("TSB_OD_0.2" in out_lib.getvalue()) and (
            "pMEM_OD_0.5" in out_lib.getvalue()) and (
            "pMEM_t2" in out_lib.getvalue()):
            detect = True
        self.assertTrue(detect)
        os.remove("test_class_aaa.png")

    def test_plot(self):
        st.plot(20, 23, 10, 13, 5, 100, 200, "name",
                "TSS", os.path.join(self.test_folder, "test"))
        self.assertTrue(os.path.exists(os.path.join(
            self.test_folder, "test_class_name.png")))

    def test_stat_tsspredator(self):
        detect = False
        tss_file = os.path.join(self.test_folder, "aaa_TSS.gff")
        stat_file = os.path.join(self.test_folder, "stat")
        lib_file = os.path.join(self.test_folder, "lib")
        gen_file(tss_file, self.example.tss)
        st.stat_tsspredator(tss_file, "TSS", stat_file, lib_file)
        datas = import_data(stat_file)
        for data in datas:
            if "Primary" in data:
                self.assertEqual(data.split(" = ")[-1], "1 (1.0)")
        datas = import_data(lib_file)
        line = "\n".join(datas)
        if ("TSB_OD_0.2" in line) and (
            "pMEM_OD_0.5" in line) and (
            "pMEM_t2" in line):
            detect = True
        self.assertTrue(detect)
        self.assertTrue(os.path.exists("TSS_class_aaa.png"))
        os.remove("TSS_class_aaa.png")

class Example(object):

    tss = """aaa	TSSpredator	TSS	2131	2131	.	+	.	UTR_length=Primary_25;type=Primary;ID=tss3;libs=TSB_OD_0.2&pMEM_OD_0.5&pMEM_t2;associated_gene=SAOUHSC_00002;Name=TSS:2131_f"""
    tss_dict = [{"seq_id": "aaa", "source": "TSSpredator",
                 "feature": "TSS", "start": 2131,
                 "end": 2131, "phase": ".", "strand": "+", "score": "."}]
    attributes_tss = [{
        "ID": "tss3", "Name": "TSS:2131_f", "UTR_length": "Primary_25",
        "type": "Primary", "associated_gene": "SAOUHSC_00002",
        "libs": "TSB_OD_0.2&pMEM_OD_0.5&pMEM_t2"}]
    tsss = []
    tsss.append(Create_generator(tss_dict[0], attributes_tss[0], "gff"))

if __name__ == "__main__":
    unittest.main()

