import sys
import os
import unittest
import shutil
import copy
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data, gen_file, extract_info
from mock_args_container import MockClass
import annogesiclib.sRNA_utr_derived as sud

get_coverage = sud.get_coverage

class Mock_func(object):

    def __init__(self):
        self.example = Example()

    def mock_read_libs(self, input_libs, wig_folder):
        return None, None

    def mock_read_wig(self, wig_file, strand, libs):
        return None

    def mock_coverage_comparison(self, cover, cover_sets, poss,
                                 first, strand):
        cover_sets['best'] = 20
        cover_sets['high'] = 50
        cover_sets['low'] = 10
        first = False
        return first

    def mock_replicate_comparison(self, srna_datas, template_texs, strand, test1,
                                  args, type_, median,
                                  coverages, srna, notex):
        datas = {"best": 500, "track": "frag", "high": 700, "low": 400,
                 "start": 100, "end": 202, "conds": {"frag_1": "track_1"}} 
        return datas

    def mock_get_coverage(self, wigs, inter, type_, pos, intercds, args):
        srna_covers = {'frag_1': [{'final_end': 20, 'final_start': 2,
                                   'track': 'track_1', 'high': 50, 'low': 10,
                                   'type': 'frag', 'avg': 41.36842105263158,
                                   'ori_avg': 27.52}]}
        utr_covers = srna_covers
        return srna_covers, utr_covers


class TestsRNAUTR(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.mock = Mock_func()
        self.mock_args = MockClass()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_import_data(self):
        pos = {"start": 4, "end": 40, "ori_start": 2, "ori_end": 3}
        datas = sud.import_data("+", "aaa", pos, "3UTR", "TSS", "cds", "srna_cover", "test")
        self.assertDictEqual(datas, {'start_cleavage': 'NA', 'strand': '+', 'end_cleavage': 'test',
                                     'start_tss': 'cds', 'end': 40, 'start': 4, 'utr': '3UTR',
                                     'strain': 'aaa', 'datas': 'srna_cover'})

    def test_read_data(self):
        args = self.mock_args.mock()
        args.gff_file = os.path.join(self.test_folder, "test.gff")
        args.ta_file = os.path.join(self.test_folder, "test.gff")
        args.tss_file = os.path.join(self.test_folder, "test.gff")
        args.pro_file = os.path.join(self.test_folder, "test.gff")
        args.seq_file = os.path.join(self.test_folder, "test.fa")
        gen_file(args.gff_file, self.example.gff_file)
        gen_file(args.seq_file, self.example.seq_file)
        args.hypo = False
        cdss, tas, tsss, pros, seq = sud.read_data(args)
        self.assertEqual(cdss[0].start, 4)
        self.assertEqual(tas[0].start, 4)
        self.assertEqual(tsss[0].start, 4)
        self.assertEqual(pros[0].start, 4)
        self.assertDictEqual(seq, {'aaa': 'ATATGACGATACGTAAACCGACCGAATATATCTTTTCACAACCAGATTACGATCGTCAT'})

    def test_get_terminal(self):
        inters = []
        seq = {"aaa": "ATATGACGATACGTAAACCGACCGAATATATCTTTTCACAACCAGATTACGATCGTCAT"}
        sud.get_terminal(self.example.gffs, inters, seq, "start")
        self.assertListEqual(inters, [{'end': 4, 'len_CDS': 0, 'strand': '+', 'strain': 'aaa', 'start': 1}])

    def test_get_inter(self):
        inters = []
        sud.get_inter(self.example.gffs, inters)
        self.assertListEqual(inters, [{'start': 14, 'strand': '+', 'end': 20, 'strain': 'aaa', 'len_CDS': 10}])

    def test_set_cover_and_point(self):
        covers = [2, 3, 4, 1, 6, 2, 8, 3, 5, 6, 7, 5, 2, 1]
        cover_results = {"covers": None, "check_point": None}
        pos = {"start": 2, "end": 6, "ori_start": 2, "ori_end": 3}
        sud.set_cover_and_point(cover_results, self.example.inters[0], covers, pos, 5)
        self.assertListEqual(cover_results["covers"], [2, 3, 4, 1, 6, 2, 8, 3, 5])
        self.assertDictEqual(cover_results["check_point"], {'srna_start': 0, 'utr_start': 2, 'utr_end': 3, 'srna_end': 12})

    def test_check_import_srna_covers(self):
        args = self.mock_args.mock()
        cover = {"type": "5utr"}
        datas = {"num": 0, "cover_tmp": {"total": 100, "ori_total": 200},
                  "checks": {"detect_decrease": True},
                 "final_poss": {"start": 3, "end": 23}}
        cover_results = {"cover_sets": {"high": 50, "low": 10},
                         "srna_covers": {"cond_1": []}, "utr_covers": {"cond_1": []},
                                         "type": "5utr", "intercds": "TSS"}
        args.min_len = 30
        args.max_len = 500
        pos = {"start": 1, "end": 25, "ori_start": 1, "ori_end": 25}
        sud.check_import_srna_covers(datas, cover_results, self.example.inters[0], "cond_1", "track",
                                     cover, pos, args)
        self.assertDictEqual(datas["final_poss"], {'end': 23, 'start': 3})
        self.assertDictEqual(cover_results["srna_covers"], {'cond_1': [{'final_start': 3, 'high': 50, 'ori_avg': 8.0,
                                                            'final_end': 23, 'low': 10, 'type': '5utr',
                                                            'avg': 4, 'track': 'track'}]})
        self.assertDictEqual(cover_results["utr_covers"], cover_results["srna_covers"])

        datas["checks"] = {"detect_decrease": False}
        cover_results["srna_covers"] = {"cond_1": []}
        cover_results["utr_covers"] = {"cond_1": []}
        sud.check_import_srna_covers(datas, cover_results, self.example.inters[0], "cond_1", "track",
                                     cover, pos, args)
        self.assertDictEqual(cover_results["srna_covers"], {'cond_1': []})

    def test_check_pos(self):
        cover = {"pos": 4}
        check_point = {"utr_start": 1, "utr_end": 29, "srna_start": 3, "srna_end": 11}
        checks = {"srna": False, "utr": False}
        sud.check_pos(cover, check_point, checks)
        self.assertDictEqual(checks, {'srna': True, 'utr': True})

    def test_get_cover_5utr(self):
        args = self.mock_args.mock()
        datas = {"num": 0, "cover_tmp": {"5utr": 0},
                 "checks": {"detect_decrease": True},
                 "final_poss": {"start": 1, "end": 26}}
        cover = {"coverage": 20, "pos": 10}
        cover_sets = {"high": 50, "low": 10}
        args.decrease_utr = 50
        args.fuzzy_utr = 2
        go_out = sud.get_cover_5utr(datas, cover_sets, cover, self.example.inters[0], args)
        self.assertDictEqual(datas["final_poss"], {'start': 1, 'end': 10})
        self.assertEqual(datas["num"], 0)
        self.assertTrue(go_out)
        self.assertDictEqual(datas["cover_tmp"], {'5utr': 0})
        self.assertDictEqual(cover_sets, {'high': 50, 'low': 10})
        cover = {"coverage": 20, "pos": 10}
        datas = {"num": 0, "cover_tmp": {"5utr": 30},
                 "checks": {"detect_decrease": True},
                 "final_poss": {"start": 1, "end": 26}}
        cover_sets = {"low": 10, "high": 50}
        args.decrease_utr = 0.5
        go_out = sud.get_cover_5utr(datas, cover_sets, cover, self.example.inters[0], args)
        self.assertEqual(datas["num"], 1)
        self.assertFalse(go_out)
        self.assertDictEqual(datas["final_poss"], {'start': 1, 'end': 26})
        self.assertDictEqual(datas["cover_tmp"], {'5utr': 20})
        self.assertDictEqual(cover_sets, {'low': 20, 'high': 50})

    def test_detect_cover_utr_srna(self):
        sud.coverage_comparison = self.mock.mock_coverage_comparison
        cover_results = {"cover_sets": {"low": 10, "high": 50}, "pos": {"low": 10, "high": 50},
                         "covers": [{"coverage": 20, "pos": 10, "type": "frag"}], "type": "5utr",
                         "srna_covers": {"frag_1": []}, "utr_covers": {"frag_1": []}, "intercds": "TSS",
                         "check_point": {"utr_start": 1, "utr_end": 29, "srna_start": 2, "srna_end": 25}}
        datas = {"num": 0, "cover_tmp": {"total": 100, "ori_total": 200},
                 "checks": {"detect_decrease": True},
                 "final_poss": {"start": 3, "end": 23}}
        pos = {"start": 2, "end": 20, "ori_start": 1, "ori_end": 23}
        args = self.mock_args.mock()
        args.min_len = 30
        args.max_len = 500
        args.decrease_utr = 0.5
        args.fuzzy_utr = 2
        sud.detect_cover_utr_srna(cover_results, pos, self.example.inters[0], "frag_1", "track_1", args)
        self.assertDictEqual(cover_results["srna_covers"], {'frag_1': [{'low': 20, 'high': 50, 'track': 'track_1',
                                           'final_start': 2, 'ori_avg': 0.8695652173913043,
                                           'type': 'frag', 'final_end': 20,
                                           'avg': 1.0526315789473684}]})
        self.assertDictEqual(cover_results["utr_covers"], cover_results["srna_covers"])
        self.assertDictEqual(cover_results["cover_sets"], {'best': 20, 'low': 20, 'high': 50})

    def test_get_coverage(self):
        sud.coverage_comparison = self.mock.mock_coverage_comparison
        pos = {"start": 2, "end": 20, "ori_start": 1, "ori_end": 25}
        args = self.mock_args.mock()
        args.min_len = 30
        args.max_len = 500
        args.decrease_utr = 0.5
        args.fuzzy_utr = 2
        srna_covers, utr_covers = sud.get_coverage(self.example.wigs, self.example.inters[0], pos, "3utr", "TSS", args)
        self.assertDictEqual(srna_covers, {'frag_1': [{'track': 'track_1', 'high': 50,
                                           'final_start': 2, 'type': 'frag',
                                           'avg': 8.052631578947368, 'low': 10,
                                           'final_end': 3, 'ori_avg': 2.12}]})
        self.assertDictEqual(utr_covers, srna_covers)

    def test_get_utr_cutoff(self):
        mediandict = {"aaa": {"5utr": {"bbb": {}}}}
        avgs = [30, 60, 550, 302, 44]
        sud.get_utr_cutoff("p_0.5", mediandict, avgs, "aaa", "5utr", "bbb")
        self.assertDictEqual(mediandict, {'aaa': {'5utr': {'bbb': {'mean': 197.2, 'median': 60}}}})

    def test_detect_normal(self):
        sud.get_coverage = self.mock.mock_get_coverage
        diff = 50
        pos = {"start": 2, "end": 20, "ori_start": 1, "ori_end": 25}
        args = self.mock_args.mock()
        args.min_len = 30
        args.max_len = 500
        args.decrease_utr = 0.5
        args.fuzzy_utr = 2
        args.utrs = []
        args.srnas = []
        sud.detect_normal(diff, self.example.wigs, self.example.inters[0], pos, "3utr", self.example.tsss[0], args)
        self.assertListEqual(args.srnas, [{'end': 20, 'strand': '+',
                                      'datas': {'frag_1': [{'track': 'track_1',
                                      'final_start': 2, 'avg': 41.36842105263158,
                                      'high': 50, 'type': 'frag', 'final_end': 20,
                                      'ori_avg': 27.52, 'low': 10}]}, 'end_cleavage': 'NA',
                                      'utr': '3utr', 'start_cleavage': 'NA', 'strain': 'aaa',
                                      'start': 2, 'start_tss': 'TSS:1_+'}])
        self.assertListEqual(args.utrs, [{'end': 20, 'strand': '+',
                                     'datas': {'frag_1': [{'track': 'track_1',
                                     'final_start': 2, 'avg': 41.36842105263158,
                                     'high': 50, 'type': 'frag', 'final_end': 20,
                                     'ori_avg': 27.52, 'low': 10}]}, 'end_cleavage': 'NA',
                                     'utr': '3utr', 'start_cleavage': 'NA', 'strain': 'aaa',
                                     'start': 2, 'start_tss': 'NA'}])
        args.utrs = []
        args.srnas = []
        args.pros = self.example.pros
        args.min_len = 3
        args.max_len = 20
        pos = {"start": 2, "end": 24, "ori_start": 1, "ori_end": 25}
        sud.detect_normal(diff, self.example.wigs, self.example.inters[0], pos, "3utr", self.example.tsss[0], args)
        self.assertListEqual(args.srnas, [{'start': 1, 'end': 18, 'start_tss': 'TSS:1_+',
                                      'datas': {'frag_1': [{'ori_avg': 27.52, 'track': 'track_1',
                                      'high': 50, 'low': 10, 'type': 'frag', 'final_end': 20,
                                      'avg': 41.36842105263158, 'final_start': 2}]},
                                      'start_cleavage': 'NA', 'end_cleavage': 'Cleavage:18_+',
                                      'utr': '3utr', 'strand': '+', 'strain': 'aaa'}])
        sud.get_coverage = get_coverage

    def test_detect_3utr_pro(self):
        sud.get_coverage = self.mock.mock_get_coverage
        args = self.mock_args.mock()
        args.min_len = 1
        args.max_len = 300
        args.decrease_utr = 0.5
        args.fuzzy_utr = 1
        args.fuzzy_tsss = {"3utr": 3}
        args.pros = self.example.pros
        args.utrs = []
        args.srnas = []
        pos = {"start": 2, "end": 20, "ori_start": 1, "ori_end": 25}
        sud.detect_3utr_pro(self.example.inters[0], pos, self.example.wigs, "3utr", args)
        self.assertListEqual(args.srnas, [{'end_cleavage': 'NA', 'end': 20,
                                      'start_cleavage': 'Cleavage:18_+', 'utr': '3utr',
                                      'datas': {'frag_1': [{'low': 10, 'final_start': 2,
                                      'track': 'track_1', 'type': 'frag', 'final_end': 20,
                                      'avg': 41.36842105263158, 'ori_avg': 27.52, 'high': 50}]},
                                      'strand': '+', 'start_tss': 'NA', 'start': 18, 'strain': 'aaa'}])
        self.assertListEqual(args.utrs, [{'end_cleavage': 'NA', 'end': 20, 'start_cleavage': 'NA',
                                     'utr': '3utr', 'datas': {'frag_1': [{'low': 10,
                                     'final_start': 2, 'track': 'track_1', 'type': 'frag',
                                     'final_end': 20, 'avg': 41.36842105263158, 'ori_avg': 27.52,
                                     'high': 50}]}, 'strand': '+', 'start_tss': 'NA',
                                     'start': 18, 'strain': 'aaa'}])
        sud.get_coverage = get_coverage

    def test_detect_twopro(self):
        sud.get_coverage = self.mock.mock_get_coverage
        pro_dict = [{"seq_id": "aaa", "source": "tsspredator", "feature": "processing", "start": 18,
                     "end": 18, "phase": ".", "strand": "+", "score": "."},
                    {"seq_id": "aaa", "source": "tsspredator", "feature": "processing", "start": 38,
                     "end": 38, "phase": ".", "strand": "+", "score": "."}]
        attributes_pro = [{"ID": "processing0", "Name": "Processing_0"},
                          {"ID": "processing1", "Name": "Processing_1"}]
        pros = []
        for index in range(0, 2):
            pros.append(Create_generator(pro_dict[index], attributes_pro[index], "gff"))
        args = self.mock_args.mock()
        args.min_len = 1
        args.max_len = 300
        args.decrease_utr = 0.5
        args.fuzzy_utr = 3
        args.fuzzy_tsss = {"3utr": 3}
        args.pros = pros
        args.utrs = []
        args.srnas = []
        pos = {"start": 2, "end": 50, "ori_start": 1, "ori_end": 25}
        sud.detect_twopro(self.example.inters[0], pos, self.example.wigs, "interCDS", "interCDS", args)
        self.assertListEqual(args.srnas, [{'start_cleavage': 'Cleavage:18_+', 'utr': 'interCDS',
                                      'datas': {'frag_1': [{'type': 'frag', 'low': 10,
                                      'final_start': 2, 'high': 50, 'avg': 41.36842105263158,
                                      'final_end': 20, 'track': 'track_1', 'ori_avg': 27.52}]},
                                      'start_tss': 'NA', 'end_cleavage': 'Cleavage:38_+',
                                      'strand': '+', 'end': 38, 'strain': 'aaa', 'start': 18}])
        self.assertListEqual(args.utrs, [{'start_cleavage': 'NA', 'utr': 'interCDS',
                                     'datas': {'frag_1': [{'type': 'frag', 'low': 10,
                                     'final_start': 2, 'high': 50, 'avg': 41.36842105263158,
                                     'final_end': 20, 'track': 'track_1', 'ori_avg': 27.52}]},
                                     'start_tss': 'NA', 'end_cleavage': 'Cleavage:38_+',
                                     'strand': '+', 'end': 38, 'strain': 'aaa', 'start': 18}])
        sud.get_coverage = get_coverage

    def test_run_utr_detection(self):
        args = self.mock_args.mock()
        args.min_len = 1
        args.max_len = 300
        args.decrease_utr = 0.5
        args.fuzzy_utr = 2
        args.fuzzy_tsss = {"5utr": "n_3"}
        args.utrs = []
        args.srnas = []
        args.tsss = self.example.tsss
        args.pros = self.example.pros
        sud.get_coverage = self.mock.mock_get_coverage
        sud.run_utr_detection(self.example.wigs, self.example.inters[0], 2, 50, "5utr", args)
        sud.get_coverage = get_coverage
        self.assertListEqual(args.srnas, [{'start': 1, 'end': 50, 'start_cleavage': 'NA',
                                      'datas': {'frag_1': [{'high': 50, 'final_end': 20,
                                      'avg': 41.36842105263158, 'low': 10, 'ori_avg': 27.52,
                                      'final_start': 2, 'type': 'frag', 'track': 'track_1'}]},
                                      'start_tss': 'TSS:1_+', 'strain': 'aaa', 'strand': '+',
                                      'utr': '5utr', 'end_cleavage': 'NA'}])
        self.assertListEqual(args.utrs, [{'start': 1, 'end': 50, 'start_cleavage': 'NA',
                                     'datas': {'frag_1': [{'high': 50, 'final_end': 20,
                                     'avg': 41.36842105263158, 'low': 10, 'ori_avg': 27.52,
                                     'final_start': 2, 'type': 'frag', 'track': 'track_1'}]},
                                     'start_tss': 'NA', 'strain': 'aaa', 'strand': '+',
                                     'utr': '5utr', 'end_cleavage': 'NA'}])

    def test_class_utr(self):
        args = self.mock_args.mock()
        args.min_len = 1
        args.max_len = 300
        args.decrease_utr = 0.5
        args.fuzzy_utr = 2
        args.fuzzy_tsss = {"3utr": "p_3"}
        args.utrs = []
        args.srnas = []
        args.tsss = self.example.tsss
        args.pros = self.example.pros
        args.wig_fs = self.example.wigs
        sud.get_coverage = self.mock.mock_get_coverage
        sud.class_utr(self.example.inters[0], self.example.tas[0], args)
        sud.get_coverage = get_coverage
        self.assertListEqual(args.srnas, [{'end_cleavage': 'NA', 'start_tss': 'TSS:1_+',
                                      'utr': '3utr', 'start_cleavage': 'NA', 'end': 20,
                                      'start': 1, 'datas': {'frag_1': [{'ori_avg': 27.52,
                                      'final_start': 2, 'avg': 41.36842105263158,
                                      'track': 'track_1', 'type': 'frag', 'final_end': 20,
                                      'low': 10, 'high': 50}]}, 'strain': 'aaa', 'strand': '+'},
                                     {'end_cleavage': 'NA', 'start_tss': 'NA', 'utr': '3utr',
                                      'start_cleavage': 'Cleavage:18_+', 'end': 20, 'start': 18,
                                      'datas': {'frag_1': [{'ori_avg': 27.52, 'final_start': 2,
                                      'avg': 41.36842105263158, 'track': 'track_1', 'type': 'frag',
                                      'final_end': 20, 'low': 10, 'high': 50}]}, 'strain': 'aaa', 'strand': '+'}])
        self.assertListEqual(args.utrs, [{'end_cleavage': 'NA', 'start_tss': 'NA', 'utr': '3utr',
                                     'start_cleavage': 'NA', 'end': 20, 'start': 1,
                                     'datas': {'frag_1': [{'ori_avg': 27.52, 'final_start': 2,
                                     'avg': 41.36842105263158, 'track': 'track_1', 'type': 'frag',
                                     'final_end': 20, 'low': 10, 'high': 50}]}, 'strain': 'aaa', 'strand': '+'},
                                    {'end_cleavage': 'NA', 'start_tss': 'NA', 'utr': '3utr',
                                     'start_cleavage': 'NA', 'end': 20, 'start': 18,
                                     'datas': {'frag_1': [{'ori_avg': 27.52, 'final_start': 2,
                                     'avg': 41.36842105263158, 'track': 'track_1', 'type': 'frag',
                                     'final_end': 20, 'low': 10, 'high': 50}]}, 'strain': 'aaa', 'strand': '+'}])

    def test_get_utr_coverage(self):
        utrs = [{'strand': '+', 'utr': '3utr', 'end': 20, 'start': 18, 'start_tss': 'NA',
                 'datas': {'frag_1': [{'final_end': 20, 'track': 'track_1', 'final_start': 2,
                 'ori_avg': 27.52, 'avg': 41.36842105263158, 'type': 'frag', 'low': 10,
                 'high': 50}]}, 'end_cleavage': 'NA', 'strain': 'aaa', 'start_cleavage': 'NA'}]
        covers = sud.get_utr_coverage(utrs)
        self.assertDictEqual(covers, {'aaa': {'interCDS': {}, '3utr': {'track_1': [27.52]}, '5utr': {}}})

    def test_set_cutoff(self):
        args = self.mock_args.mock()
        args.texs = {"track_4@AND@track_6": 0}
        covers = {'aaa': {'5utr': {'track_4': [52, 11, 23]}, 'inter': {'track_3': [111]},
                  'total': {'track_1': [27.52, 111]}, '3utr': {'track_1': [27.52, 111]},
                  'interCDS': {'track_2': [12, 0]}}}
        args.coverages = {"5utr": "p_0.3", "3utr": "n_10", "interCDS": "p_0.5"}
        args.cover_notex = {"5utr": "p_0.3", "3utr": "n_10", "interCDS": "p_0.5"}
        mediandict = sud.set_cutoff(covers, args)
        self.assertDictEqual(mediandict, {'aaa': {'5utr': {'track_4': {'median': 11, 'mean': 28.666666666666668}},
                                                  'interCDS': {'track_2': {}}, '3utr': {'track_1': {}}}})
        args.cover_notex = None
        mediandict = sud.set_cutoff(covers, args)
        self.assertDictEqual(mediandict, {'aaa': {'3utr': {'track_1': {'mean': 69.26, 'median': 10.0}},
                                                  '5utr': {'track_4': {'mean': 28.666666666666668, 'median': 11}},
                                                  'interCDS': {'track_2': {'mean': 6.0, 'median': 0}}}})

    def test_mean_score(self):
        lst = [1, 3, 5, 6, 7, 8]
        mean = sud.mean_score(lst)
        self.assertEqual(mean, 5.0)

    def test_median_score(self):
        lst = [1, 3, 5, 6, 7, 8]
        median = sud.median_score(lst, 0.5)
        self.assertEqual(median, 5)

    def test_detect_srna(self):
        sud.replicate_comparison = self.mock.mock_replicate_comparison
        args = self.mock_args.mock()
        args.min_len = 1
        args.max_len = 300
        args.decrease_utr = 0.5
        args.fuzzy_utr = 2
        args.coverages = "cover"
        args.texs = "template_texs"
        args.tex_notex = "tex_notex"
        args.replicates = "rep"
        args.table_best = True
        args.out = StringIO()
        args.out_t = StringIO()
        median = {"aaa": {"3utr": 555}}
        args.srnas = [{'strand': '+', 'utr': '3utr', 'end': 20,
                  'start': 18, 'start_tss': 'NA',
                  'datas': {'frag_1': [{'final_end': 20, 'track': 'track_1',
                  'final_start': 2, 'ori_avg': 27.52, 'avg': 41.36842105263158,
                  'type': 'frag', 'low': 10, 'high': 50, "conds": ["frag"]}]}, 'end_cleavage': 'NA',
                  'strain': 'aaa', 'start_cleavage': 'Cleavage:18_+'}]
        sud.detect_srna(median, args)
        self.assertEqual(args.out.getvalue(), "aaa\tANNOgesic\tncRNA\t18\t20\t.\t+\t.\tID=srna_utr0;Name=UTR_sRNA_00000;sRNA_type=3utr;best_avg_coverage=500;best_high_coverage=700;best_low_coverage=400;with_TSS=NA;start_cleavage=Cleavage:18_+;end_cleavage=NA\n")
        self.assertEqual(args.out_t.getvalue(), "aaa\t00000\t18\t20\t+\tfrag_1\ttrack_1\t500\t700\t400\tfrag(avg=500;high=700;low=400)\n")

    def test_print_file(self):
        args = self.mock_args.mock()
        args.min_len = 1
        args.max_len = 300
        args.decrease_utr = 0.5
        args.fuzzy_utr = 2
        args.coverages = "cover"
        args.texs = "template_texs"
        args.tex_notex = "tex_notex"
        args.replicates = "rep"
        args.table_best = True
        args.out = StringIO()
        args.out_t = StringIO()
        srna = {'strand': '+', 'utr': '3utr', 'end': 20,
                'start': 18, 'start_tss': 'NA',
                'datas': {'frag_1': [{'final_end': 20, 'track': 'track_1',
                'final_start': 2, 'ori_avg': 27.52, 'avg': 41.36842105263158,
                'type': 'frag', 'low': 10, 'high': 50, "conds": ["frag"]}]}, 'end_cleavage': 'NA',
                'strain': 'aaa', 'start_cleavage': 'Cleavage:18_+'}
        srna_datas = {"best": 500, "track": "frag", "high": 700, "low": 400,
                      "start": 100, "end": 202, "conds": {"frag_1": "track_1"}}
        sud.print_file(0, srna, 2, 50, srna_datas, args)
        self.assertEqual(args.out.getvalue(), "aaa\tANNOgesic\tncRNA\t2\t50\t.\t+\t.\tID=srna_utr0;Name=UTR_sRNA_00000;sRNA_type=3utr;best_avg_coverage=500;best_high_coverage=700;best_low_coverage=400;with_TSS=NA;start_cleavage=Cleavage:18_+;end_cleavage=NA\n")
        self.assertEqual(args.out_t.getvalue(), "aaa\t00000\t2\t50\t+\tfrag_1\ttrack_1\t500\t700\t400\tfrag(avg=500;high=700;low=400)\n")


class Example(object):

    gff_file = """aaa	Refseq	CDS	4	14	.	+	.	ID=cds0;Name=CDS_00000"""
    seq_file = """>aaa
ATATGACGATACGTAAACCGACCGAATATATCTTTTCACAACCAGATTACGATCGTCAT"""
    gff_dict = [{"seq_id": "aaa", "source": "RefSeq", "feature": "CDS", "start": 4,
                 "end": 14, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "RefSeq", "feature": "CDS", "start": 20,
                 "end": 30, "phase": ".", "strand": "+", "score": "."}]
    attributes_gff = [{"ID": "cds0", "Name": "CDS_0"},
                      {"ID": "cds1", "Name": "CDS_1"}]
    gffs = []
    for index in range(0, 2):
        gffs.append(Create_generator(gff_dict[index], attributes_gff[index], "gff"))
    inters = [{'strand': '+', 'start': 14, 'end': 20, 'cleavage': 'NA',
               'strain': 'aaa', 'utr': '', 'datas': None, 'tss': 'NA', "len_CDS": 1000}]
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
    tss_dict = [{"seq_id": "aaa", "source": "tsspredator", "feature": "TSS", "start": 1,
                "end": 1, "phase": ".", "strand": "+", "score": "."}]
    attributes_tss = [{"ID": "tss0", "Name": "TSS_0"}]
    tsss = []
    tsss.append(Create_generator(tss_dict[0], attributes_tss[0], "gff"))
    pro_dict = [{"seq_id": "aaa", "source": "tsspredator", "feature": "processing", "start": 18,
                "end": 18, "phase": ".", "strand": "+", "score": "."}]
    attributes_pro = [{"ID": "processing0", "Name": "Processing_0"}]
    pros = []
    pros.append(Create_generator(pro_dict[0], attributes_pro[0], "gff"))
    ta_dict = [{"seq_id": "aaa", "source": "RefSeq", "feature": "transcript", "start": 2,
                "end": 20, "phase": ".", "strand": "+", "score": "."}]
    attributes_ta = [{"ID": "tran0", "Name": "Tran_0"}]
    tas = []
    tas.append(Create_generator(ta_dict[0], attributes_ta[0], "gff"))


if __name__ == "__main__":
    unittest.main()

