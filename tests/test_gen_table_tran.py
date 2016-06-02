import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
import annogesiclib.gen_table_tran as gtt


class Mock_func(object):

    def __init__(self):
        self.example = Example()

class TestGenTableTran(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_detect_coverage(self):
        infos = {}
        gtt.detect_coverage(self.example.wigs, self.example.tas[0], infos)
        self.assertDictEqual(infos, {'track_1': {'high': 100, 'avg': 33.411764705882355, 'low': 2}})

    def test_print_coverage(self):
        out = StringIO()
        out_gff = StringIO()
        gtt.print_coverage(self.example.tas, out, out_gff, self.example.wigs, self.example.wigs, True)
        self.assertEqual(out.getvalue(), "aaa\tTranscript_0\t4\t20\t+\tfragmented&TEX+/-\tNA\tNA\tNA\ttrack_1(avg=33.411764705882355)\n")
        self.assertListEqual(out_gff.getvalue().split("\t")[:-1], ["aaa", "ANNOgesic", "Transcript", "4", "20", ".", "+", "."])
        self.assertEqual(set(out_gff.getvalue().split("\t")[-1].strip().split(";")), set(["Name=Transcript_0", "detect_lib=fragmented&TEX+/-",
                                                                                          "best_avg_coverage=33.411764705882355", "ID=tran0"]))
class Example(object):
    wigs = {"aaa": {"frag_1": {"track_1": [{"strand": "+", "pos": 1, "coverage": 100, "type": "frag"},
                                           {"strand": "+", "pos": 2, "coverage": 30, "type": "frag"},
                                           {"strand": "+", "pos": 3, "coverage": 23, "type": "frag"},
                                           {"strand": "+", "pos": 4, "coverage": 21, "type": "frag"},
                                           {"strand": "+", "pos": 5, "coverage": 21, "type": "frag"},
                                           {"strand": "+", "pos": 6, "coverage": 2, "type": "frag"},
                                           {"strand": "+", "pos": 7, "coverage": 100, "type": "frag"},
                                           {"strand": "+", "pos": 8, "coverage": 30, "type": "frag"},
                                           {"strand": "+", "pos": 9, "coverage": 23, "type": "frag"},
                                           {"strand": "+", "pos": 10, "coverage": 21, "type": "frag"},
                                           {"strand": "+", "pos": 11, "coverage": 21, "type": "frag"},
                                           {"strand": "+", "pos": 12, "coverage": 2, "type": "frag"},
                                           {"strand": "+", "pos": 13, "coverage": 100, "type": "frag"},
                                           {"strand": "+", "pos": 14, "coverage": 30, "type": "frag"},
                                           {"strand": "+", "pos": 15, "coverage": 23, "type": "frag"},
                                           {"strand": "+", "pos": 16, "coverage": 21, "type": "frag"},
                                           {"strand": "+", "pos": 17, "coverage": 21, "type": "frag"},
                                           {"strand": "+", "pos": 18, "coverage": 2, "type": "frag"},
                                           {"strand": "+", "pos": 19, "coverage": 100, "type": "frag"},
                                           {"strand": "+", "pos": 20, "coverage": 30, "type": "frag"},
                                           {"strand": "+", "pos": 21, "coverage": 23, "type": "frag"},
                                           {"strand": "+", "pos": 22, "coverage": 21, "type": "frag"},
                                           {"strand": "+", "pos": 23, "coverage": 21, "type": "frag"},
                                           {"strand": "+", "pos": 24, "coverage": 2, "type": "frag"}]}}}

    ta_dict = [{"seq_id": "aaa", "source": "ANNOgesic", "feature": "Transcript", "start": 4,
                "end": 20, "phase": ".", "strand": "+", "score": "."}]
    attributes_tas = [{"ID": "tran0", "Name": "Transcript_0", "detect_lib": "fragmented&tex_notex"}]
    tas = []
    for index in range(0, 1):
        tas.append(Create_generator(ta_dict[index], attributes_tas[index], "gff"))


if __name__ == "__main__":
    unittest.main()

