#!/usr/bin/python

import unittest
import os
import sys
import shutil
sys.path.append(".")
from io import StringIO
from mock_gff3 import Create_generator
from mock_helper import gen_file, import_data
import annogesiclib.gene_express_analysis as gea


class MockFunc(object):

    def __init__(self):
        self.example = Example()

    def mock_read_libs(self, input_libs, wigs):
        return None, {"tex1_tex2": 0}

    def mock_read_wig(self, wig_file, strand, libs):
        return self.example.wig_texs

    def mock_read_data(self, gff, features):
        stats = {"CDS": {}}
        outs = {"CDS": {"all": [], "least_one": [], "none": []}}
        gff_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 3,
                    "end": 5, "phase": ".", "strand": "+", "score": "."}
        attributes_gff = {"ID": "CDS0", "Name": "CDS_0", "locus_tag": "AAA_00001"}
        gff = Create_generator(gff_dict, attributes_gff, "gff")
        gff_list = {"CDS": [gff]}
        return gff_list, stats, outs

class TestGeneExpress(unittest.TestCase):

    def setUp(self):
        self.test_folder = "test_project"
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)
        os.mkdir(self.test_folder)
        self.example = Example()

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_set_cutoff(self):
        detects = {}
        detects["express"] = 100
        gff_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 3,
                    "end": 102, "phase": ".", "strand": "+", "score": "."}
        attributes_gff = {"ID": "CDS0", "Name": "CDS_0", "locus_tag": "AAA_00001"}
        gff = Create_generator(gff_dict, attributes_gff, "gff")
        diff, cutoff_percent = gea.set_cutoff("tex", "all", "all", detects, gff)
        self.assertEqual(diff, 100)
        self.assertEqual(cutoff_percent, 0)
        diff, cutoff_percent = gea.set_cutoff("frag", "all", "n_50", detects, gff)
        self.assertEqual(diff, 100)
        self.assertEqual(cutoff_percent, 50)
        diff, cutoff_percent = gea.set_cutoff("tex", "p_0.5", "n_50", detects, gff)
        self.assertEqual(diff, 1.0)
        self.assertEqual(cutoff_percent, 0.5)

    def test_detect_express(self):
        gff_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 3,
                    "end": 5, "phase": ".", "strand": "+", "score": "."}
        attributes_gff = {"ID": "CDS0", "Name": "CDS_0", "locus_tag": "AAA_00001"}
        gff = Create_generator(gff_dict, attributes_gff, "gff")
        texs = {"tex1_tex2": 0}
        detects = {"cond": 0, "track": 0, "import": False, "express": 0}
        gea.detect_express(self.example.wig_frags["aaa"]["frag"]["track_1"], gff, 5,
                           detects, "all", "all", texs, "frag", 2, "track_1")
        self.assertDictEqual({'track': 1, 'import': False, 'cond': 0, 'express': 2}, detects)

    def test_compare_wigs(self):
        gff_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 3,
                    "end": 5, "phase": ".", "strand": "+", "score": "."}
        attributes_gff = {"ID": "CDS0", "Name": "CDS_0", "locus_tag": "AAA_00001"}
        gff = Create_generator(gff_dict, attributes_gff, "gff")
        texs = {"tex1_tex2": 0}
        replicates = {"tex": 1, "frag": 1}
        stats = {"CDS": {"total": {"total": 0, "least_one": 0, "all": 0, "none": 0},
                         "aaa": {"total": 0, "least_one": 0, "all": 0, "none": 0}}}
        outs = {"CDS": {"least_one": [], "all": [], "none": []}}
        gea.compare_wigs(self.example.wig_texs, gff, 2, texs, replicates, stats["CDS"], outs["CDS"],
                         5, "all", "all")
        self.assertDictEqual(stats, {'CDS': {'total': {'tex': 1, 'none': 0, 'total': 0, 'all': 1, 'least_one': 1},
                                     'aaa': {'tex': 1, 'none': 0, 'total': 0, 'all': 1, 'least_one': 1}}})

    def test_gene_expression(self):
        gea.read_wig = MockFunc().mock_read_wig
        gea.read_libs = MockFunc().mock_read_libs
        gea.read_data = MockFunc().mock_read_data
        replicates = {"tex": 1, "frag": 1}
        stat_folder = os.path.join(self.test_folder, "stat")
        gff_folder = os.path.join(self.test_folder, "gff")
        if os.path.exists(gff_folder):
            shutil.rmtree(gff_folder)
        os.mkdir(gff_folder)
        gen_file(os.path.join(gff_folder, "aaa.gff"), "test")
        if not os.path.exists(stat_folder):
            os.mkdir(stat_folder)
        out_gff_folder = os.path.join(self.test_folder, "out_gff")
        if not os.path.exists(out_gff_folder):
            os.mkdir(out_gff_folder)
        gea.gene_expression(None, gff_folder, "all", "all", "test_wig", "test_wig", ["CDS"],
                           "test_wig_folder", 5, 2, replicates, stat_folder, out_gff_folder)
        datas = import_data(os.path.join(stat_folder, "aaa_CDS.csv"))
        dicts = {}
        for data in datas:
            dicts[data] = data
        refs = {}
        for data in self.example.out_stat.split("\n"):
            refs[data] = data
        self.assertDictEqual(dicts, refs)

class Example(object):
    wig_frags = {"aaa": {"frag": {"track_1": [{"coverage": 100, "pos": 1},
                                              {"coverage": 100, "pos": 2},
                                              {"coverage": 200, "pos": 3},
                                              {"coverage": 0, "pos": 4},
                                              {"coverage": 230, "pos": 5},
                                              {"coverage": 230, "pos": 6}]}}}
    wig_texs = {"aaa": {"tex": {"tex1": [{"coverage": 100, "pos": 1},
                                         {"coverage": 100, "pos": 2},
                                         {"coverage": 0, "pos": 3},
                                         {"coverage": 0, "pos": 4},
                                         {"coverage": 230, "pos": 5},
                                         {"coverage": 230, "pos": 6}],
                                "tex2": [{"coverage": 100, "pos": 1},
                                         {"coverage": 100, "pos": 2},
                                         {"coverage": 100, "pos": 3},
                                         {"coverage": 0, "pos": 4},
                                         {"coverage": 230, "pos": 5},
                                         {"coverage": 230, "pos": 6}]}}}

    out_stat = """aaa:
total input:	1
expression at all conditions:	1 (1.0)
expression at lease one condition:	1 (1.0)
no expression:	0 (0.0)
condition tex:	1 (1.0)"""

if __name__ == "__main__":
    unittest.main()
