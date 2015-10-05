#!/usr/bin/python

import unittest
import os
import sys
import shutil
sys.path.append(".")
from io import StringIO
from mock_gff3 import Create_generator
import annogesiclib.filter_low_expression as fle


class TestLowExpress(unittest.TestCase):

    def setUp(self):
        self.test_folder = "test_project"
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)
        os.mkdir(self.test_folder)
        self.example = Example()

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_get_coverage(self):
        tar_dist = {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 3,
               "end": 3, "phase": ".", "strand": "+", "score": "."}
        attributes_tar = {"type": "Primary", "ID": "tss1", "Name": "TSS:3_+"}
        tar = Create_generator(tar_dist, attributes_tar, "gff")
        wigs = {"aaa": {"frag_1": {"track_1": [{"pos": 1, "coverage": 100, "type": "frag"},
                                               {"pos": 2, "coverage": 30, "type": "frag"},
                                               {"pos": 3, "coverage": 23, "type": "frag"},
                                               {"pos": 4, "coverage": 21, "type": "frag"},
                                               {"pos": 5, "coverage": 21, "type": "frag"}]},
                        "frag_2": {"track_2": [{"pos": 1, "coverage": 100, "type": "frag"},
                                               {"pos": 2, "coverage": 30, "type": "frag"},
                                               {"pos": 3, "coverage": 40, "type": "frag"},
                                               {"pos": 4, "coverage": 21, "type": "frag"},
                                               {"pos": 5, "coverage": 21, "type": "frag"}]}}}
        coverage = fle.get_coverage(tar, wigs)
        self.assertEqual(coverage, 40)

    def test_stat(self):
        stats, num = fle.stat(self.example.tars, self.example.refs, 5, 3000, 3)
        self.assertDictEqual(stats, {'fp': 1, 'fp_rate': 0.000333667000333667,
                                     'tp_rate': 0.3333333333333333, 'tp': 1,
                                     'miss_rate': 0.6666666666666666, 'miss': 2})
        self.assertEqual(num, 3)

    def test_change_best(self):
        best = {"tp": 100, "fp": 20, "miss": 10, "fp_rate": 0.0003,
                "tp_rate": 0.93, "miss_rate": 0.00001}
        stats1 = {"tp": 120, "fp": 10, "miss": 5, "fp_rate": 0.0002,
                 "tp_rate": 0.96, "miss_rate": 0.00001}
        stats2 = {"tp": 90, "fp": 30, "miss": 15, "fp_rate": 0.0004,
                  "tp_rate": 0.90, "miss_rate": 0.00011}
        best, change = fle.change_best(1000, best, stats1)
        self.assertDictEqual(best, stats1)
        self.assertTrue(change)
        best, change = fle.change_best(1000, best, stats2)
        self.assertDictEqual(best, stats1)
        self.assertFalse(change)

class Example(object):
    tar_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 3,
                 "end": 3, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 24,
                 "end": 24, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 1243,
                 "end": 1243, "phase": ".", "strand": "+", "score": "."}]
    attributes_tar = [{"coverage": "3", "ID": "tss1", "Name": "TSS:3_+"},
                      {"coverage": "340", "ID": "tss2", "Name": "TSS:24_+"},
                      {"coverage": "4440", "ID": "tss3", "Name": "TSS:1243_+"}]
    ref_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 3,
                 "end": 3, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 333,
                 "end": 333, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 1242,
                 "end": 1242, "phase": ".", "strand": "+", "score": "."}]
    attributes_ref = [{"coverage": "3", "ID": "tss1", "Name": "TSS:3_+"},
                      {"coverage": "330", "ID": "tss2", "Name": "TSS:333_+"},
                      {"coverage": "1230", "ID": "tss3", "Name": "TSS:1242_+"}]
    tars = []
    refs = []
    for index in range(0, 3):
        tars.append(Create_generator(tar_dict[index], attributes_tar[index], "gff"))
        tars[-1].attributes["print"] = False
        refs.append(Create_generator(ref_dict[index], attributes_ref[index], "gff"))
        refs[-1].attributes["print"] = False

if __name__ == "__main__":
    unittest.main()
