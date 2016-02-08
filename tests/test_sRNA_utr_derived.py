import sys
import os
import unittest
import shutil
import copy
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data, gen_file, extract_info
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
                                  tex_notex, replicates, type_, median,
                                  coverages, srna, texs, notex):
        datas = {"best": 500, "track": "frag", "high": 700, "low": 400,
                 "start": 100, "end": 202, "conds": {"frag_1": "track_1"}} 
        return datas

    def mock_get_coverage(self, wigs, inter, start, end, type_,
                          intercds_type, fuzzy_end, decrease,
                          ori_start, ori_end, max_len, min_len):
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
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_import_data(self):
        datas = sud.import_data("+", "aaa", 4, 40, "3UTR", "TSS", "cds", "srna_cover", "test")
        self.assertDictEqual(datas, {'start_cleavage': 'NA', 'strand': '+', 'end_cleavage': 'test',
                                     'start_tss': 'cds', 'end': 40, 'start': 4, 'utr': '3UTR',
                                     'strain': 'aaa', 'datas': 'srna_cover'})

    def test_read_data(self):
        gff_file = os.path.join(self.test_folder, "test.gff")
        seq_file = os.path.join(self.test_folder, "test.fa")
        gen_file(gff_file, self.example.gff_file)
        gen_file(seq_file, self.example.seq_file)
        cdss, tas, tsss, pros, seq = sud.read_data(gff_file, gff_file,
                                         gff_file, gff_file, seq_file, False)
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
        covers, check_point = sud.set_cover_and_point(self.example.inters[0], covers,
                                                      2, 6, 2, 3, 5)
        self.assertListEqual(covers, [2, 3, 4, 1, 6, 2, 8, 3, 5])
        self.assertDictEqual(check_point, {'srna_start': 1, 'utr_start': 2, 'srna_end': 8, 'utr_end': 6})

    def test_check_import_srna_covers(self):
        cover_tmp = {"total": 100, "ori_total": 200}
        checks = {"detect_decrease": True}
        cover_sets = {"high": 50, "low": 10}
        srna_covers = {"cond_1": []}
        utr_covers = {"cond_1": []}
        cover = {"type": "5utr"}
        final_poss = {"start": 3, "end": 23}
        sud.check_import_srna_covers(checks, final_poss, 20, 2, "5utr", cover_sets,
                                     cover_tmp, "TSS", self.example.inters[0],
                                     srna_covers, "cond_1", "track", cover, 1, 25,
                                     utr_covers, 30, 500)
        self.assertDictEqual(final_poss, {'end': 20, 'start': 3})
        self.assertDictEqual(srna_covers, {'cond_1': [{'final_start': 3, 'high': 50, 'ori_avg': 8.0,
                                                       'final_end': 20, 'low': 10, 'type': '5utr',
                                                       'avg': 5.2631578947368425, 'track': 'track'}]})
        self.assertDictEqual(utr_covers, srna_covers)
        checks = {"detect_decrease": False}
        srna_covers = {"cond_1": []}
        utr_covers = {"cond_1": []}
        sud.check_import_srna_covers(checks, final_poss, 20, 2, "5utr", cover_sets,
                                     cover_tmp, "TSS", self.example.inters[0],
                                     srna_covers, "cond_1", "track", cover, 1, 25,
                                     utr_covers, 30, 500)
        self.assertDictEqual(srna_covers, {'cond_1': []})

    def test_check_pos(self):
        cover = {"pos": 4}
        check_point = {"utr_start": 1, "utr_end": 29, "srna_start": 3, "srna_end": 11}
        checks = {"srna": False, "utr": False}
        sud.check_pos(cover, check_point, checks)
        self.assertDictEqual(checks, {'srna': True, 'utr': True})

    def test_get_cover_5utr(self):
        num = 0
        cover_tmp = {"5utr": 0}
        cover = {"coverage": 20, "pos": 10}
        cover_sets = {"low": 10, "high": 50}
        final_poss = {"start": 1, "end": 26}
        num, go_out = sud.get_cover_5utr(num, 2, cover_sets, cover_tmp,
                           final_poss, 50, cover, self.example.inters[0])
        self.assertDictEqual(final_poss, {'start': 1, 'end': 10})
        self.assertEqual(num, 0)
        self.assertTrue(go_out)
        self.assertDictEqual(cover_tmp, {'5utr': 0})
        self.assertDictEqual(cover_sets, {'high': 50, 'low': 10})
        cover_tmp = {"5utr": 30}
        cover = {"coverage": 20, "pos": 10}
        cover_sets = {"low": 10, "high": 50}
        final_poss = {"start": 1, "end": 26}
        num, go_out = sud.get_cover_5utr(num, 2, cover_sets, cover_tmp,
                           final_poss, 0.5, cover, self.example.inters[0])
        self.assertEqual(num, 1)
        self.assertFalse(go_out)
        self.assertDictEqual(final_poss, {'start': 1, 'end': 26})
        self.assertDictEqual(cover_tmp, {'5utr': 20})
        self.assertDictEqual(cover_sets, {'low': 20, 'high': 50})

    def test_detect_cover_utr_srna(self):
        sud.coverage_comparison = self.mock.mock_coverage_comparison
        cover_sets = {"low": 10, "high": 50}
        srna_covers = {"frag_1": []}
        utr_covers = {"frag_1": []}
        covers = [{"coverage": 20, "pos": 10, "type": "frag"}]
        check_point = {"utr_start": 1, "utr_end": 29, "srna_start": 2, "srna_end": 25}
        sud.detect_cover_utr_srna(covers, check_point, cover_sets, 2, 20,
                          "poss", self.example.inters[0], "5utr", "TSS", 0.5,
                          2, srna_covers, utr_covers, "frag_1",
                          "track_1", 1, 23, 30, 500)
        self.assertDictEqual(srna_covers, {'frag_1': [{'low': 20, 'high': 50, 'track': 'track_1',
                                           'final_start': 2, 'ori_avg': 0.8695652173913043,
                                           'type': 'frag', 'final_end': 20,
                                           'avg': 1.0526315789473684}]})
        self.assertDictEqual(utr_covers, srna_covers)
        self.assertDictEqual(cover_sets, {'best': 20, 'low': 20, 'high': 50})

    def test_get_coverage(self):
        sud.coverage_comparison = self.mock.mock_coverage_comparison
        srna_covers, utr_covers = sud.get_coverage(self.example.wigs,
                                      self.example.inters[0], 2, 20, "3utr",
                                      "TSS", 2, 0.5, 1, 25, 30, 500)
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
        utrs = []
        srnas = []
        diff = 50
        sud.detect_normal(diff, self.example.wigs, self.example.inters[0], 2, 20,
                          "3utr", "TSS", "srna_1", self.example.tsss[0],
                          self.example.pros, utrs, srnas, 2, 0.5, 30, 300, 1, 25)
        self.assertListEqual(srnas, [{'end': 20, 'strand': '+',
                                      'datas': {'frag_1': [{'track': 'track_1',
                                      'final_start': 2, 'avg': 41.36842105263158,
                                      'high': 50, 'type': 'frag', 'final_end': 20,
                                      'ori_avg': 27.52, 'low': 10}]}, 'end_cleavage': 'NA',
                                      'utr': '3utr', 'start_cleavage': 'NA', 'strain': 'aaa',
                                      'start': 2, 'start_tss': 'TSS:1_+'}])
        self.assertListEqual(utrs, [{'end': 20, 'strand': '+',
                                     'datas': {'frag_1': [{'track': 'track_1',
                                     'final_start': 2, 'avg': 41.36842105263158,
                                     'high': 50, 'type': 'frag', 'final_end': 20,
                                     'ori_avg': 27.52, 'low': 10}]}, 'end_cleavage': 'NA',
                                     'utr': '3utr', 'start_cleavage': 'NA', 'strain': 'aaa',
                                     'start': 2, 'start_tss': 'srna_1'}])
        utrs = []
        srnas = []
        sud.detect_normal(diff, self.example.wigs, self.example.inters[0], 2, 24,
                          "3utr", "cleavage", "srna_1", self.example.tsss[0],
                          self.example.pros, utrs, srnas, 0, 0.5, 3, 20, 1, 25)
        self.assertListEqual(srnas, [{'start': 1, 'end': 18, 'start_tss': 'TSS:1_+',
                                      'datas': {'frag_1': [{'ori_avg': 27.52, 'track': 'track_1',
                                      'high': 50, 'low': 10, 'type': 'frag', 'final_end': 20,
                                      'avg': 41.36842105263158, 'final_start': 2}]},
                                      'start_cleavage': 'NA', 'end_cleavage': 'Cleavage:18_+',
                                      'utr': '3utr', 'strand': '+', 'strain': 'aaa'}])
        sud.get_coverage = get_coverage

    def test_detect_3utr_pro(self):
        utrs = []
        srnas = []
        fuzzys = {"3utr": 3}
        sud.get_coverage = self.mock.mock_get_coverage
        sud.detect_3utr_pro(self.example.pros, self.example.inters[0], utrs, srnas,
                            2, 20, self.example.wigs, "feature", "srna_0",
                            "3utr", fuzzys, 3, 0.5, 1, 300)
        self.assertListEqual(srnas, [{'end_cleavage': 'NA', 'end': 20,
                                      'start_cleavage': 'Cleavage:18_+', 'utr': '3utr',
                                      'datas': {'frag_1': [{'low': 10, 'final_start': 2,
                                      'track': 'track_1', 'type': 'frag', 'final_end': 20,
                                      'avg': 41.36842105263158, 'ori_avg': 27.52, 'high': 50}]},
                                      'strand': '+', 'start_tss': 'NA', 'start': 18, 'strain': 'aaa'}])
        self.assertListEqual(utrs, [{'end_cleavage': 'NA', 'end': 20, 'start_cleavage': 'NA',
                                     'utr': '3utr', 'datas': {'frag_1': [{'low': 10,
                                     'final_start': 2, 'track': 'track_1', 'type': 'frag',
                                     'final_end': 20, 'avg': 41.36842105263158, 'ori_avg': 27.52,
                                     'high': 50}]}, 'strand': '+', 'start_tss': 'NA',
                                     'start': 18, 'strain': 'aaa'}])
        sud.get_coverage = get_coverage

    def test_detect_twopro(self):
        utrs = []
        srnas = []
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
        sud.detect_twopro(pros, self.example.inters[0], utrs, srnas,
                          2, 50, self.example.wigs, "feature", "srna_0",
                          "interCDS", 3, 0.5, 1, 300, "interCDS")
        self.assertListEqual(srnas, [{'start_cleavage': 'Cleavage:18_+', 'utr': 'interCDS',
                                      'datas': {'frag_1': [{'type': 'frag', 'low': 10,
                                      'final_start': 2, 'high': 50, 'avg': 41.36842105263158,
                                      'final_end': 20, 'track': 'track_1', 'ori_avg': 27.52}]},
                                      'start_tss': 'NA', 'end_cleavage': 'Cleavage:38_+',
                                      'strand': '+', 'end': 38, 'strain': 'aaa', 'start': 18}])
        self.assertListEqual(utrs, [{'start_cleavage': 'NA', 'utr': 'interCDS',
                                     'datas': {'frag_1': [{'type': 'frag', 'low': 10,
                                     'final_start': 2, 'high': 50, 'avg': 41.36842105263158,
                                     'final_end': 20, 'track': 'track_1', 'ori_avg': 27.52}]},
                                     'start_tss': 'NA', 'end_cleavage': 'Cleavage:38_+',
                                     'strand': '+', 'end': 38, 'strain': 'aaa', 'start': 18}])
        sud.get_coverage = get_coverage


    def test_run_utr_detection(self):
        utrs = []
        srnas = []
        fuzzys = {"5utr": "n_3"}
        sud.get_coverage = self.mock.mock_get_coverage
        sud.run_utr_detection(self.example.wigs, self.example.inters[0],
                              2, 50, "5utr", utrs, self.example.tsss,
                              self.example.pros, srnas, "feature", "srna_0",
                              fuzzys, 2, 0.5, 1, 300)
        sud.get_coverage = get_coverage
        self.assertListEqual(srnas, [{'start': 1, 'end': 50, 'start_cleavage': 'NA',
                                      'datas': {'frag_1': [{'high': 50, 'final_end': 20,
                                      'avg': 41.36842105263158, 'low': 10, 'ori_avg': 27.52,
                                      'final_start': 2, 'type': 'frag', 'track': 'track_1'}]},
                                      'start_tss': 'TSS:1_+', 'strain': 'aaa', 'strand': '+',
                                      'utr': '5utr', 'end_cleavage': 'NA'}])
        self.assertListEqual(utrs, [{'start': 1, 'end': 50, 'start_cleavage': 'NA',
                                     'datas': {'frag_1': [{'high': 50, 'final_end': 20,
                                     'avg': 41.36842105263158, 'low': 10, 'ori_avg': 27.52,
                                     'final_start': 2, 'type': 'frag', 'track': 'track_1'}]},
                                     'start_tss': 'NA', 'strain': 'aaa', 'strand': '+',
                                     'utr': '5utr', 'end_cleavage': 'NA'}])

    def test_class_utr(self):
        utrs = []
        srnas = []
        fuzzys = {"3utr": "p_0.3"}
        sud.get_coverage = self.mock.mock_get_coverage
        sud.class_utr(self.example.inters[0], self.example.tas[0], utrs, srnas,
                      self.example.tsss, self.example.pros, self.example.wigs,
                      self.example.wigs, fuzzys, 30, 0.5, 1, 300)
        sud.get_coverage = get_coverage
        self.assertListEqual(srnas, [{'end_cleavage': 'NA', 'start_tss': 'TSS:1_+',
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
        self.assertListEqual(utrs, [{'end_cleavage': 'NA', 'start_tss': 'NA', 'utr': '3utr',
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
        texs = {"track_4@AND@track_6": 0}
        covers = {'aaa': {'5utr': {'track_4': [52, 11, 23]}, 'inter': {'track_3': [111]},
                  'total': {'track_1': [27.52, 111]}, '3utr': {'track_1': [27.52, 111]},
                  'interCDS': {'track_2': [12, 0]}}}
        coverages = {"5utr": "p_0.3", "3utr": "n_10", "interCDS": "p_0.5"}
        mediandict = sud.set_cutoff(covers, coverages, coverages, texs)
        self.assertDictEqual(mediandict, {'aaa': {'5utr': {'track_4': {'median': 11, 'mean': 28.666666666666668}},
                                                  'interCDS': {'track_2': {}}, '3utr': {'track_1': {}}}})
        mediandict = sud.set_cutoff(covers, coverages, None, texs)
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
        out = StringIO()
        out_t = StringIO()
        median = {"aaa": {"3utr": 555}}
        srnas = [{'strand': '+', 'utr': '3utr', 'end': 20,
                  'start': 18, 'start_tss': 'NA',
                  'datas': {'frag_1': [{'final_end': 20, 'track': 'track_1',
                  'final_start': 2, 'ori_avg': 27.52, 'avg': 41.36842105263158,
                  'type': 'frag', 'low': 10, 'high': 50, "conds": ["frag"]}]}, 'end_cleavage': 'NA',
                  'strain': 'aaa', 'start_cleavage': 'Cleavage:18_+'}]
        sud.detect_srna(srnas, median, out, out_t, "template_texs", "covers",
                        "tex_notex", "rep", True, 300, 1)
        self.assertEqual(out.getvalue(), "aaa\tANNOgesic\tsRNA\t18\t20\t.\t+\t.\tID=srna_utr0;Name=UTR_sRNA_00000;sRNA_type=3utr;best_avg_coverage=500;best_high_coverage=700;best_low_coverage=400;with_TSS=NA;start_cleavage=Cleavage:18_+;end_cleavage=NA\n")
        self.assertEqual(out_t.getvalue(), "aaa\t00000\t18\t20\t+\tfrag_1\ttrack_1\t500\t700\t400\tfrag(avg=500;high=700;low=400)\n")

    def test_print_file(self):
        out = StringIO()
        out_t = StringIO()
        srna = {'strand': '+', 'utr': '3utr', 'end': 20,
                'start': 18, 'start_tss': 'NA',
                'datas': {'frag_1': [{'final_end': 20, 'track': 'track_1',
                'final_start': 2, 'ori_avg': 27.52, 'avg': 41.36842105263158,
                'type': 'frag', 'low': 10, 'high': 50, "conds": ["frag"]}]}, 'end_cleavage': 'NA',
                'strain': 'aaa', 'start_cleavage': 'Cleavage:18_+'}
        srna_datas = {"best": 500, "track": "frag", "high": 700, "low": 400,
                      "start": 100, "end": 202, "conds": {"frag_1": "track_1"}}
        sud.print_file(0, out_t, out, srna, 2, 50, srna_datas, True)
        self.assertEqual(out.getvalue(), "aaa\tANNOgesic\tsRNA\t2\t50\t.\t+\t.\tID=srna_utr0;Name=UTR_sRNA_00000;sRNA_type=3utr;best_avg_coverage=500;best_high_coverage=700;best_low_coverage=400;with_TSS=NA;start_cleavage=Cleavage:18_+;end_cleavage=NA\n")
        self.assertEqual(out_t.getvalue(), "aaa\t00000\t2\t50\t+\tfrag_1\ttrack_1\t500\t700\t400\tfrag(avg=500;high=700;low=400)\n")


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

