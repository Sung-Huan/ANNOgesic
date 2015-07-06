import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
import annogesiclib.get_polyT as gpt


class TestGetPolyT(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_filter_term(self):
        cands = [{"r_stem": 6, "l_stem": 6, "miss": 2, "print": False, "strain": "aaa", "end": 500, "start": 30},
                 {"r_stem": 6, "l_stem": 6, "miss": 0, "print": False, "strain": "aaa", "end": 400, "start": 50},
                 {"r_stem": 3, "l_stem": 3, "miss": 2, "print": False, "strain": "aaa", "end": 50, "start": 10},
                 {"r_stem": 6, "l_stem": 6, "miss": 0, "print": False, "strain": "bbb", "end": 450, "start": 60},]
        terms = []
        gpt.filter_term(cands, terms)
        self.assertDictEqual(terms[0], {'end': 500, 'l_stem': 6, 'miss': 0, 'print': False,
                                        'strain': 'aaa', 'r_stem': 6, 'start': 30})
        self.assertDictEqual(terms[1], {'end': 450, 'l_stem': 6, 'miss': 0, 'print': False,
                                        'strain': 'bbb', 'r_stem': 6, 'start': 60})

    def test_check_sec(self):
        sec = "...((((((((..))))))..))...."
        features, detect = gpt.check_sec(sec, 26)
        self.assertDictEqual(features, {'rights': 8, 'loop': 2, 'real_miss': 2, 'lefts': 8,
                                        'l_stem': 0, 'st_pos': 23, 'r_stem': 10, 'tmp_miss': 2})
        self.assertDictEqual(detect, {'detect_l': True, 'conflict': False, 'detect_r': True})
        sec = "........))))))..))...."
        features, detect = gpt.check_sec(sec, 21)
        self.assertDictEqual(features, {'r_stem': 0, 'loop': 0, 'lefts': 0, 'l_stem': 0,
                                        'st_pos': 21, 'tmp_miss': 10, 'rights': 8, 'real_miss': 2})
        self.assertDictEqual(detect, {'detect_r': True, 'conflict': False, 'detect_l': False})

    def test_detect_candidates(self):
        seq = "GATCGGCAGTATTAAACGTACTTTTTTTTTT"
        sec = "...((((((((....))))))..))......"
        cands = gpt.detect_candidates(seq, sec, "test", "aaa", 30, 58,
                                      "AAA_00001", "AAA_00002", "+")
        refs = [{'l_stem': 4, 'name': 'test', 'strain': 'aaa', 'detect_m': False,
                 'parent_p': 'AAA_00001', 'miss': 0, 'ut': 4, 'parent_m': 'AAA_00002',
                 'detect_p': False, 'loop': 4, 'end': 58, 'strand': '+',
                 'print': False, 'length': 12, 'r_stem': 4, 'start': 27},
                {'l_stem': 5, 'name': 'test', 'strain': 'aaa', 'detect_m': False,
                 'parent_p': 'AAA_00001', 'miss': 0, 'ut': 5, 'parent_m': 'AAA_00002',
                 'detect_p': False, 'loop': 4, 'end': 59, 'strand': '+',
                 'print': False, 'length': 14, 'r_stem': 5, 'start': 26},
                {'l_stem': 6, 'name': 'test', 'strain': 'aaa', 'detect_m': False,
                 'parent_p': 'AAA_00001', 'miss': 0, 'ut': 6, 'parent_m': 'AAA_00002',
                 'detect_p': False, 'loop': 4, 'end': 60, 'strand': '+',
                 'print': False, 'length': 16, 'r_stem': 6, 'start': 25},
                {'l_stem': 7, 'name': 'test', 'strain': 'aaa', 'detect_m': False,
                 'parent_p': 'AAA_00001', 'miss': 2, 'ut': 6, 'parent_m': 'AAA_00002',
                 'detect_p': False, 'loop': 4, 'end': 63, 'strand': '+',
                 'print': False, 'length': 20, 'r_stem': 9, 'start': 24}]
        for index in range(len(cands)):
            self.assertDictEqual(cands[index], refs[index])

    def test_check_parent(self):
        term = {"strain": "aaa", "start": 11, "end": 14}
        detects = {"parent_p": False, "parent_m": False}
        parent = gpt.check_parent(self.example.cdss, term, detects, "+", 3, 3, "parent_p")
        self.assertEqual(parent, "AAA_00001")
        parent = gpt.check_parent(self.example.cdss, term, detects, "-", 3, 3, "parent_m")
        self.assertEqual(parent, "AAA_00002")

    def test_parents(self):
        terms = [{"strain": "aaa", "start": 11, "end": 14, "parent_p": "AAA_00001", "parent_m": "AAA_00002"},
                 {"strain": "aaa", "start": 12, "end": 15, "parent_p": "tran0:1-11_+", "parent_m": "tran1:16-30_-"}]
        gpt.parents(terms, self.example.cdss, 10, 10, 10, 10)
        self.assertDictEqual(terms[0], {'parent_p': 'AAA_00001', 'parent_m': 'AAA_00002',
                                        'start': 11, 'strain': 'aaa', 'end': 14})
        self.assertDictEqual(terms[1], {'parent_p': 'tran0:1-11_+,AAA_00001',
                                        'parent_m': 'tran1:16-30_-,AAA_00002', 'start': 12,
                                        'strain': 'aaa', 'end': 15})

    def test_compare_anno(self):
        terms = [{"strain": "aaa", "start": 11, "end": 14, "strand": "+"},
                 {"strain": "aaa", "start": 9, "end": 18, "strand": "-"},
                 {"strain": "aaa", "start": 209, "end": 218, "strand": "-"}]
        cands = gpt.compare_anno(self.example.cdss, terms, 3, 3)
        self.assertDictEqual(terms[0], {'strand': '+', 'start': 11, 'strain': 'aaa', 'end': 14})
        self.assertDictEqual(terms[1], {"strain": "aaa", "start": 9, "end": 18, "strand": "-"})

class Example(object):
    cds_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 3,
                 "end": 10, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 16,
                 "end": 30, "phase": ".", "strand": "-", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 38,
                 "end": 50, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 54,
                 "end": 60, "phase": ".", "strand": "+", "score": "."}]
    attributes_cds = [{"ID": "CDS0", "Name": "CDS_0", "locus_tag": "AAA_00001"},
                      {"ID": "CDS1", "Name": "CDS_1", "protein_id": "AAA_00002"},
                      {"ID": "CDS2", "Name": "CDS_2"},
                      {"ID": "CDS3", "Name": "CDS_3"}]
    cdss = []
    for index in range(0, 4):
        cdss.append(Create_generator(cds_dict[index], attributes_cds[index], "gff"))

if __name__ == "__main__":
    unittest.main()
