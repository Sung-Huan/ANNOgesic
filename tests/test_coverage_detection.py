import sys
import math
import csv
import os
import unittest
from io import StringIO
sys.path.append(".")
import annogesiclib.coverage_detection as cover_detect

class Mock_func(object):

    def mock_check_tex(self, template_texs, covers, cutoff,
                       target_datas, tex_notex, type_,
                       poss, median, coverages, utr_type):
        for srna in Example().cover_datas:
            target_datas.append(srna)
        poss["start"] = 100
        poss["end"] = 202
        return 2

class TestCoverageDetection(unittest.TestCase):

    def setUp(self):
        self.example = Example()

    def test_coverage_comparison_first(self):
        first = True
        cover_sets = {"high": -1, "low": -1}
        poss = {"high": -1, "low": -1}
        cover = {"coverage": 100, "pos": 50}
        cover_detect.coverage_comparison(cover, cover_sets, poss, first, "+")
        self.assertDictEqual(cover_sets, {"high": 100, "low": 100})
        self.assertDictEqual(poss, {"high": 50, "low": 50})
        
    def test_coverage_comparison_forward(self):
        first = False
        cover_sets = {"high": 50, "low": 20}
        poss = {"high": 10, "low": 30}
        cover = {"coverage": 100, "pos": 50}
        cover_detect.coverage_comparison(cover, cover_sets, poss, first, "+")
        self.assertDictEqual(cover_sets, {"high": 100, "low": 100})
        self.assertDictEqual(poss, {"high": 50, "low": 50})
        cover = {"coverage": 30, "pos": 51}
        cover_detect.coverage_comparison(cover, cover_sets, poss, first, "+")
        self.assertDictEqual(cover_sets, {"high": 100, "low": 30})
        self.assertDictEqual(poss, {"high": 50, "low": 51})

    def test_coverage_comparison_reverse(self):
        first = False
        cover_sets = {"high": 50, "low": 20}
        poss = {"high": 30, "low": 10}
        cover = {"coverage": 100, "pos": 50}
        cover_detect.coverage_comparison(cover, cover_sets, poss, first, "-")
        self.assertDictEqual(cover_sets, {"high": 100, "low": 100})
        self.assertDictEqual(poss, {"high": 50, "low": 50})
        cover = {"coverage": 30, "pos": 49}
        cover_detect.coverage_comparison(cover, cover_sets, poss, first, "-")
        self.assertDictEqual(cover_sets, {"high": 100, "low": 30})
        self.assertDictEqual(poss, {"high": 50, "low": 49})

    def test_define_cutoff_median(self):
        coverages = {"3utr": "mean", "5utr": "median"}
        median = {"track_a": {"median": 100, "mean": 200}, "track_b": {"median": 30, "mean": 80}}
        cutoff = cover_detect.define_cutoff(coverages, median, "5utr")
        self.assertDictEqual(cutoff, {'track_a': 100, 'track_b': 30})
        cutoff = cover_detect.define_cutoff(coverages, median, "3utr")
        self.assertDictEqual(cutoff, {'track_a': 200, 'track_b': 80})

    def test_check_tex(self):
        template_texs = self.example.texs
        covers = self.example.cover_datas
        coverages = {"3utr": 100, "5utr": 600}
        poss = {"high": 30, "low": 10}
        median = {"track1_tex": {"median": 100, "mean": 200}, "track1_notex": {"median": 30, "mean": 80},
                  "track2_tex": {"median": 150, "mean": 200}, "track2_notex": {"median": 10, "mean": 20},
                  "frag": {"median": 80, "mean": 100}}
        target_datas = []
        detect_num_lib = cover_detect.check_tex(template_texs, covers, 200, target_datas, 2,
                                            None, poss, median, coverages, "3utr")
        self.assertEqual(detect_num_lib, 2)
        num_frag = 0
        num_tex = 0
        for target in target_datas:
            if target["type"] == "frag":
                num_frag += 1
            else:
                num_tex += 1
        self.assertEqual(num_frag, 1)
        self.assertEqual(num_tex, 2)
        detect_num_lib = cover_detect.check_tex(template_texs, covers, 200, target_datas, 2,
                                                "sRNA_utr_derived", poss, median, coverages, "5utr")
        self.assertEqual(detect_num_lib, 0)
        self.assertDictEqual(poss, {"high": 30, "low": 10})

    def test_replicate_comparison(self):
        cover_detect.check_tex = Mock_func().mock_check_tex
        template_texs = self.example.texs
        srna_covers = {"texnotex": self.example.cover_datas}
        coverages = {"3utr": 100, "5utr": 600}
        median = {"track1_tex": {"median": 100, "mean": 200}, "track1_notex": {"median": 30, "mean": 80},
                  "track2_tex": {"median": 150, "mean": 200}, "track2_notex": {"median": 10, "mean": 20},
                  "frag": {"median": 80, "mean": 100}}
        srna_datas = cover_detect.replicate_comparison(srna_covers, template_texs, "+",
                                                      200, 2, {"tex": 2, "frag": 1},
                                                      "sRNA_utr_derived", median, coverages, "3utr", False)
        self.assertEqual(srna_datas["best"], 500)
        self.assertEqual(srna_datas["track"], "frag")
        self.assertEqual(srna_datas["high"], 700)
        self.assertEqual(srna_datas["low"], 400)
        self.assertEqual(srna_datas["start"], 100)
        self.assertEqual(srna_datas["end"], 202)

class Example(object):
    cover_datas = [{"avg": 300, "type": "tex", "track": "track1_tex",
                    "final_start": 100, "final_end": 200, "high": 500, "low": 100},
                   {"avg": 300, "type": "notex", "track": "track1_notex",
                    "final_start": 100, "final_end": 200, "high": 500, "low": 100},
                   {"avg": 10, "type": "tex", "track": "track2_tex",
                    "final_start": 101, "final_end": 200, "high": 50, "low": 0},
                   {"avg": 10, "type": "notex", "track": "track2_notex",
                    "final_start": 101, "final_end": 200, "high": 50, "low": 0},
                   {"avg": 500, "type": "frag", "track": "frag",
                    "final_start": 100, "final_end": 202, "high": 700, "low": 400}]
    texs = {"track1_tex_track1_notex": 0,
            "track2_tex_track2_notex": 0}

    cutoff = {"track1_tex": 200, "track1_notex": 200,
              "track2_tex": 200, "track2_notex": 200,
              "frag": 200}

if __name__ == "__main__":
    unittest.main()
