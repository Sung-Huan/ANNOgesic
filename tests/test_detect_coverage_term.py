import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import read_dict
import annogesiclib.detect_coverage_term as dct

class Mock_coverage(object):

    def __init__(self):
        self.example = Example()

    def coverage_comparison(self, cover, cover_sets, poss, first, strand):
        cover_sets["low"] = 30
        cover_sets["high"] = 100
        poss["low"] = 21 
        poss["high"] = 20
        first = False
        return False
        
class Mock_Gff_parser(object):

    def __init__(self):
        self.example = Example()

    def entries(self, fh):
        for line in fh:
            if "gff" in line:
                lists = self.example.gff_dict
                attributes = self.example.attributes_gff
                num = 3
            elif "tran" in line:
                lists = self.example.tran_dict
                attributes = self.example.attributes_tran
                num = 3
            elif "term" in line:
                lists = self.example.term_dict
                attributes = self.example.attributes_term
                num = 3
            elif "tss" in line:
                lists = self.example.tss_dict
                attributes = self.example.attributes_tss
                num = 3
            elif "utr5" in line:
                lists = self.example.utr5_dict
                attributes = self.example.attributes_utr5
                num = 2
            elif "utr3" in line:
                lists = self.example.utr3_dict
                attributes = self.example.attributes_utr3
                num = 2
        for index in range(0, num):
            yield Create_generator(lists[index], attributes[index], "gff")

class TestCoverageTerminator(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_compare_ta(self):
        trans = read_dict(3, self.example.tran_dict, self.example.attributes_tran)
        dct.compare_ta(self.example.term_dict, trans, 5)
        express = []
        for term in self.example.term_dict:
            express.append(term["express"])
        self.assertListEqual(express, ["True", "True", "False"])

    def test_compare_transtermhp(self):
        hps = read_dict(3, self.example.hp_dict, self.example.attributes_term)
        terms = dct.compare_transtermhp(hps, self.example.term_dict)
        terms = sorted(terms, key=lambda x: (x["strain"], x["start"]))
        poss = []
        methods = []
        for term in terms:
            poss.append("_".join([str(term["start"]), str(term["end"])]))
            methods.append(term["method"])
        self.assertListEqual(poss, ['30_40', '350_367', '420_432', '1420_2429'])
        self.assertListEqual(methods, ['TransTermHP', 'forward_reverse&TransTermHP',
                                   'forward_reverse&TransTermHP', 'forward_reverse'])

    def test_compare_replicates(self):
        texs = {"track_tex_track_notex": 0}
        replicates = {"tex": 1, "frag": 1}
        cond = "texnotex"
        term_covers = [{"track": "track_tex", "high": 300,
                        "low": 50, "detect": "True",
                        "diff": 250, "type": "tex"},
                       {"track": "track_notex", "high": 200,
                        "low": 50, "detect": "True",
                        "diff": 150, "type": "notex"}]
        diff_cover, diff, term_datas, detect_num = \
            dct.compare_replicates(term_covers, texs, cond, 2, replicates)
        self.assertEqual(diff_cover, 250)
        self.assertDictEqual(diff, {'track': 'track_tex', 'detect': 'True',
                                    'high': 300, 'low': 50, 'type': 'tex', 'diff': 250})
        ref_datas = [{'track': 'track_notex', 'detect': 'True', 'high': 200, 'low': 50, 'type': 'notex', 'diff': 150},
                     {'track': 'track_tex', 'detect': 'True', 'high': 300, 'low': 50, 'type': 'tex', 'diff': 250}]
        for index in range(0, 2):
            self.assertDictEqual(ref_datas[index], term_datas[index])
        self.assertEqual(detect_num, 1)
        replicates = {"tex": 1, "frag": 1}
        cond = "frag"
        term_covers = [{"track": "frag", "high": 10,
                        "low": 0, "detect": "False",
                        "diff": 10, "type": "frag"}]
        diff_cover, diff, term_datas, detect_num = \
            dct.compare_replicates(term_covers, texs, cond, 2, replicates)
        self.assertEqual(diff_cover, 10)
        self.assertDictEqual(diff, {'detect': 'False', 'type': 'frag', 'low': 0, 'diff': 10, 'track': 'frag', 'high': 10})
        self.assertDictEqual(term_datas[0], {'detect': 'False', 'type': 'frag', 'low': 0, 'diff': 10, 'track': 'frag', 'high': 10})
        self.assertEqual(detect_num, 1)

    def test_coverage2term(self):
        dct.coverage_comparison = Mock_coverage().coverage_comparison
        hl_covers = {"low": 20, "high": 30}
        hl_poss = {"low": 1, "high": 2}
        term = {"start": 2, "end": 4}
        covers = [{"coverage": 100, "pos": 1, "type": "frag"},
                  {"coverage": 30, "pos": 2, "type": "frag"},
                  {"coverage": 23, "pos": 3, "type": "frag"},
                  {"coverage": 21, "pos": 4, "type": "frag"},
                  {"coverage": 21, "pos": 5, "type": "frag"},]
        term_covers = []
        dct.coverage2term(covers, term, 1, hl_covers, hl_poss, "+",
                          0.5, term_covers, "track_1")
        self.assertDictEqual(term_covers[0], {'diff': 70, 'track': 'track_1', 'type': 'frag', 'high': 100, 'low': 30, 'detect': 'True'})

    def test_get_coverage(self):
        term = {"start": 2, "end": 4, "strain": "aaa", "strand": "+"}
        texs = {"track_tex_track_notex": 0}
        replicates = {"tex": 1, "frag": 1}
        wigs = {"aaa": {"frag_1": {"track_1": [{"pos": 1, "coverage": 100, "type": "frag"},
                                               {"pos": 2, "coverage": 30, "type": "frag"},
                                               {"pos": 3, "coverage": 23, "type": "frag"},
                                               {"pos": 4, "coverage": 21, "type": "frag"},
                                               {"pos": 5, "coverage": 21, "type": "frag"}]}}}
        diff_cover, diff, term_datas, detect_nums = dct.get_coverage(
                                                    term, wigs, "+", texs, 1,
                                                    0.5, replicates, 2)
        self.assertEqual(diff_cover, 70)
        self.assertDictEqual(diff, {'track': 'track_1', 'high': 100, 'type': 'frag', 'detect': 'True', 'diff': 70, 'low': 30})
        self.assertDictEqual(term_datas["frag_1"][0],
                             {'track': 'track_1', 'high': 100, 'type': 'frag', 'detect': 'True', 'diff': 70, 'low': 30})
        self.assertDictEqual(detect_nums, {'frag_1': 1})

    def test_compare_term(self):
        terms = []
        term = {"miss": 5, "diff_cover": 30, "ut": 4}
        terms = dct.compare_term(term, terms)
        self.assertDictEqual(terms[0], term)
        term = {"miss": 4, "diff_cover": 30, "ut": 4}
        terms = dct.compare_term(term, terms)
        self.assertDictEqual(terms[0], term)
        term = {"miss": 6, "diff_cover": 80, "ut": 4}
        terms = dct.compare_term(term, terms)
        self.assertDictEqual(terms[0], {"miss": 4, "diff_cover": 30, "ut": 4})
        term = {"miss": 4, "diff_cover": 80, "ut": 4}
        terms = dct.compare_term(term, terms)
        self.assertDictEqual(terms[0], term)
        term = {"miss": 4, "diff_cover": 80, "ut": 6}
        terms = dct.compare_term(term, terms)
        self.assertDictEqual(terms[0], term)
        terms = dct.compare_term(term, terms)
        self.assertDictEqual(terms[0], term)
        self.assertDictEqual(terms[1], term)

    def test_first_term(self):
        detect_terms = {"detect": [], "undetect": []}
        detect = False
        term = {"detect_p": True, "detect_m": False}
        detect = dct.first_term("+", term, detect_terms, detect)
        self.assertTrue(detect)
        self.assertDictEqual(detect_terms["detect"][0], term)
        detect = False
        detect = dct.first_term("-", term, detect_terms, detect)
        self.assertFalse(detect)
        self.assertDictEqual(detect_terms["undetect"][0], term)

    def test_print_table(self):
        cutoff_coverage = 5
        out_t = StringIO()
        term = {"express": "True", "diff_cover": 70, "diff": {"high": 100, "low": 30, "track": "track_1"},
                "datas": {"data": [{"track": "track_1", "diff": 70, "high": 100, "low": 30},
                                   {"track": "track_2", "diff": 39, "high": 99, "low": 60}]}}
        dct.print_table(term, cutoff_coverage, out_t, True)
        self.assertEqual(set(out_t.getvalue().split("\n")), set(["	True	track_1(diff=70;high=100;low=30)"]))
        out_t.close()
        out_t = StringIO()
        dct.print_table(term, cutoff_coverage, out_t, False)
        self.assertEqual(set(out_t.getvalue().split("\n")), set(["	True	track_1(diff=70;high=100;low=30);track_2(diff=39;high=99;low=60)"]))
        term = {"express": "False", "diff_cover": 70, "diff": {"high": 100, "low": 30, "track": "track_1"},
                "datas": {"data": [{"track": "track_1", "diff": 70, "high": 100, "low": 30},
                                   {"track": "track_2", "diff": 39, "high": 99, "low": 60}]}}
        out_t.close()
        out_t = StringIO()
        dct.print_table(term, cutoff_coverage, out_t, False)
        self.assertEqual(set(out_t.getvalue().split("\n")), set(["	False	NA"]))
        term = {"express": "True", "diff_cover": -1, "diff": {"high": 100, "low": 30, "track": "track_1"},
                "datas": {"data": [{"track": "track_1", "diff": 70, "high": 100, "low": 30},
                                   {"track": "track_2", "diff": 39, "high": 99, "low": 60}]}}
        out_t.close()
        out_t = StringIO()
        dct.print_table(term, cutoff_coverage, out_t, False)
        self.assertEqual(set(out_t.getvalue().split("\n")), set(["	False	No_coverage_decreasing"]))
        out_t.close()

    def test_print2file(self):
        out = StringIO()
        out_t = StringIO()
        term = {"strain": "aaa", "express": "True", "diff_cover": 70,
                "strand": "+", "start": 2, "end": 4,
                "diff": {"high": 100, "low": 30, "track": "track_1"},
                "datas": {"data": [{"track": "track_1", "diff": 70, "high": 100, "low": 30},
                                   {"track": "track_2", "diff": 39, "high": 99, "low": 60}]}}
        dct.print2file(0, term, "70", "test", out, out_t,
                       "test_method", True, 5)
        self.assertEqual(set(out.getvalue().split("\n")[:-1]), set([self.example.gff_file]))
        self.assertEqual(set(out_t.getvalue().split("\n")[:-1]), set([self.example.table]))
        out.close()
        out_t.close()

class Example(object):

    term_dict = [{"strain": "aaa", "start": 350, "end": 367, "strand": "+",
                  "method": "forward_reverse", "express": "False"},
                 {"strain": "bbb", "start": 420, "end": 429, "strand": "-",
                  "method": "forward_reverse", "express": "False"},
                 {"strain": "bbb", "start": 1420, "end": 2429, "strand": "-",
                  "method": "forward_reverse", "express": "False"}]
    tran_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "Transcript", "start": 140,
                  "end": 367, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "Transcript", "start": 30,
                  "end": 40, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "bbb", "source": "Refseq", "feature": "Transcript", "start": 430,
                  "end": 567, "phase": ".", "strand": "-", "score": "."}]
    hp_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "Terminator", "start": 360,
                "end": 367, "phase": ".", "strand": "+", "score": "."},
               {"seq_id": "aaa", "source": "Refseq", "feature": "Terminator", "start": 30,
                "end": 40, "phase": ".", "strand": "+", "score": "."},
               {"seq_id": "bbb", "source": "Refseq", "feature": "Terminator", "start": 429,
                "end": 432, "phase": ".", "strand": "-", "score": "."}]
    attributes_term = [{"ID": "term0", "Name": "Term_0", "associated_gene": "AAA_00001"},
                       {"ID": "term1", "Name": "Term_1", "associated_gene": "AAA_00002"},
                       {"ID": "term2", "Name": "Term_2", "associated_gene": "AAA_00003"}]
    attributes_tran = [{"ID": "tran0", "Name": "Tran_0"},
                       {"ID": "tran1", "Name": "Tran_1"},
                       {"ID": "tran2", "Name": "Tran_2"}]
    gff_file = """aaa	test_method	terminator	2	4	.	+	.	ID=term_0;Name=Term_00000;associate=test;coverage_decrease=70;diff_coverage=track_1(high:100,low:30);express=True"""
    table = """aaa	Term_00000	2	4	+	True	track_1(diff=70;high=100;low=30)"""

if __name__ == "__main__":
    unittest.main()

