import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_args_container import MockClass
import annogesiclib.get_polyT as gpt


class TestGetPolyT(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.mock_args = MockClass()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_filter_term(self):
        cands = [{"r_stem": 6, "l_stem": 6, "miss": 2, "print": False,
                  "strain": "aaa", "end": 500, "start": 30},
                 {"r_stem": 6, "l_stem": 6, "miss": 0, "print": False,
                  "strain": "aaa", "end": 400, "start": 50},
                 {"r_stem": 3, "l_stem": 3, "miss": 2, "print": False,
                  "strain": "aaa", "end": 50, "start": 10},
                 {"r_stem": 6, "l_stem": 6, "miss": 0, "print": False,
                  "strain": "bbb", "end": 450, "start": 60},]
        terms = []
        gpt.filter_term(cands, terms, 0.25)
        self.assertDictEqual(terms[0], {
            'end': 500, 'l_stem': 6, 'miss': 0, 'print': False,
            'strain': 'aaa', 'r_stem': 6, 'start': 30})
        self.assertDictEqual(terms[1], {
            'end': 450, 'l_stem': 6, 'miss': 0, 'print': False,
            'strain': 'bbb', 'r_stem': 6, 'start': 60})

    def test_check_sec(self):
        sec = "...((((((((..))))))..))...."
        features, detect = gpt.check_sec(sec, 26)
        self.assertDictEqual(features, {
            'rights': 8, 'loop': 2, 'real_miss': 2, 'lefts': 8,
            'l_stem': 0, 'st_pos': 23, 'r_stem': 10, 'tmp_miss': 2})
        self.assertDictEqual(detect, {
            'detect_l': True, 'conflict': False, 'detect_r': True})
        sec = "........))))))..))...."
        features, detect = gpt.check_sec(sec, 21)
        self.assertDictEqual(features, {
            'r_stem': 0, 'loop': 0, 'lefts': 0, 'l_stem': 0,
            'st_pos': 21, 'tmp_miss': 10, 'rights': 8, 'real_miss': 2})
        self.assertDictEqual(detect, {
            'detect_r': True, 'conflict': False, 'detect_l': False})

    def test_detect_candidates(self):
        seq = "GATCGGCAGTATTAAACGTACTTTTTTTTTT"
        sec = "...((((((((....))))))..))......"
        args = self.mock_args.mock()
        args.max_loop = 10
        args.min_loop = 3
        args.max_stem = 20
        args.min_stem = 4
        args.miss_rate = 0.25
        args.at_tail = 3
        args.range_u = 6
        cands = gpt.detect_candidates(
            seq, sec, "test", "aaa", 30, 58,
            "gene_0", "gene_1", "+", args, "10-24", "70-100")
        refs = {'strand': '+', 'parent_m': 'gene_1', 'start': 32,
                'm_pos': '70-100', 'end': 60, 'detect_p': False,
                'name': 'test', 'detect_m': False, 'strain': 'aaa',
                'length': 23, 'print': False, 'ut': 6, 'loop': 4,
                'p_pos': '10-24', 'r_stem': 10, 'parent_p': 'gene_0',
                'miss': 2, 'l_stem': 9}
        self.assertDictEqual(cands[0], refs)

    def test_check_parent(self):
        term = {"strain": "aaa", "start": 11, "end": 14}
        detects = {"parent_p": False, "parent_m": False}
        parent = gpt.check_parent(
            self.example.cdss, term, detects, "+", 3, 3, "parent_p")
        self.assertEqual(parent, "gene_0")
        parent = gpt.check_parent(
            self.example.cdss, term, detects, "-", 3, 3, "parent_m")
        self.assertEqual(parent, "gene_1")

    def test_parents(self):
        terms = [{"strain": "aaa", "start": 11, "end": 14,
                  "parent_p": "gene_0", "parent_m": "gene_1",
                  "p_pos": "3-5", "m_pos": "20-50"},
                 {"strain": "aaa", "start": 12, "end": 15,
                  "parent_p": "tran0:1-11_+", "parent_m": "tran1:16-30_-",
                  "p_pos": "1-11", "m_pos": "16-30"}]
        args = self.mock_args.mock()
        args.fuzzy_up_gene = 10
        args.fuzzy_up_ta = 10
        args.fuzzy_down_gene = 10
        args.fuzzy_down_ta = 10
        gpt.parents(terms, self.example.cdss, args)
        self.assertDictEqual(terms[0], {
            'parent_p': 'gene_0', 'parent_m': 'gene_1',
            'start': 11, 'strain': 'aaa', 'end': 14,
            "p_pos": "3-5", "m_pos": "20-50"})
        self.assertDictEqual(terms[1], {
            'parent_p': 'tran0:1-11_+,gene_0',
            'parent_m': 'tran1:16-30_-,gene_1', 'start': 12,
            'strain': 'aaa', 'end': 15,
            "p_pos": "1-11", "m_pos": "16-30"})

    def test_compare_anno(self):
        terms = [{"strain": "aaa", "start": 11, "end": 14, "strand": "+"},
                 {"strain": "aaa", "start": 9, "end": 18, "strand": "-"},
                 {"strain": "aaa", "start": 209, "end": 218, "strand": "-"}]
        cands = gpt.compare_anno(self.example.cdss, terms, 3, 3)
        self.assertDictEqual(terms[0], {
            'strand': '+', 'start': 11, 'strain': 'aaa', 'end': 14})
        self.assertDictEqual(terms[1], {
            "strain": "aaa", "start": 9, "end": 18, "strand": "-"})

class Example(object):
    cds_dict = [
        {"seq_id": "aaa", "source": "Refseq", "feature": "gene", "start": 3,
         "end": 10, "phase": ".", "strand": "+", "score": "."},
        {"seq_id": "aaa", "source": "Refseq", "feature": "gene", "start": 16,
         "end": 30, "phase": ".", "strand": "-", "score": "."},
        {"seq_id": "aaa", "source": "Refseq", "feature": "gene", "start": 38,
         "end": 50, "phase": ".", "strand": "+", "score": "."},
        {"seq_id": "aaa", "source": "Refseq", "feature": "gene", "start": 54,
         "end": 60, "phase": ".", "strand": "+", "score": "."}]
    attributes_cds = [
        {"ID": "gene0", "Name": "gene_0", "locus_tag": "AAA_00001"},
        {"ID": "gene1", "Name": "gene_1", "protein_id": "AAA_00002"},
        {"ID": "gene2", "Name": "gene_2"},
        {"ID": "gene3", "Name": "gene_3"}]
    cdss = []
    for index in range(0, 4):
        cdss.append(Create_generator(cds_dict[index],
                                     attributes_cds[index], "gff"))
if __name__ == "__main__":
    unittest.main()
