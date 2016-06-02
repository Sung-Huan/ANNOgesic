import sys
import os
import unittest
import shutil
import copy
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data, gen_file
import annogesiclib.sORF_detection as sd
from mock_args_container import MockClass


get_coverage = copy.deepcopy(sd.get_coverage)

class Mock_func(object):

    def mock_get_coverage(self, inter_datas, wigs, strand, background, test1, test2):
        return "2"

    def mock_replicate_comparison(self, srna_covers, template_texs, strand, cutoff_coverage,
                                  tex_notex, type_, median, coverages,
                                  utr_type, notex):
        return {"best": 20, "high": 50, "low": 10, "start": 1, "end": 10,
                "track": "track_1", "detail": [], "conds": {"frag": "track_1"}}

    def mock_read_libs(self, input_libs, wig_folder):
        return None, None

    def mock_read_wig(self, wig_file, strand, libs):
        return None

    def mock_get_inter_coverage(self, inters, inter_covers):
        pass


class TestsORFDetection(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.mock_args = MockClass()
        self.mock = Mock_func()
        self.test_folder = "test_folder"
        self.fasta = "test_folder/fasta"
        self.wigs = "test_folder/wig"
        self.gff = "test_folder/gff"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
            os.mkdir(self.fasta)
            os.mkdir(self.wigs)
            os.mkdir(self.gff)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_get_coverage(self):
        coverages = {"3utr": "median", "5utr": "median",
                     "inter": 5, "interCDS": "median"}
        medianlist = {"aaa": {"3utr": {"track_1": {"median": 3}},
                              "5utr": {"track_1": {"median": 6}},
                              "interCDS": {"track_1": {"median": 2}},
                              "inter": {"track_1": {"median": 5}}}}
        cutoffs = {"track_1": 0}
        sorf = {"strain": "aaa", "strand": "+", "start": 2, "end": 6,
                "starts": [str(2)], "ends": [str(10)], "seq": "ATGTA",
                "type": "3utr", "print": False, "rbs": [1]}
        covers = sd.get_coverage(sorf, self.example.wigs, "+", coverages, medianlist, cutoffs)
        self.assertDictEqual(covers, {'frag_1': [{'avg': 19.4, 'pos': 2, 'type': 'frag',
                                                  'high': 30, 'track': 'track_1', 'low': 2}]})

    def test_detect_rbs_site(self):
        args = self.mock_args.mock()
        args.max_len = 20
        args.min_len = 3
        args.fuzzy_rbs = 2
        detect = sd.detect_rbs_site("AGGAGGCCGCTATGCCACACGT", 2, self.example.tas[0], args)
        self.assertListEqual(detect, [1])

    def test_detect_start_stop(self):
        seq = {"aaa": "TAGGAGGCCGCTATGCCATTA"}
        args = self.mock_args.mock()
        args.start_codon = ["ATG"]
        args.stop_codon = ["TTA"]
        args.max_len = 20
        args.min_len = 3
        args.fuzzy_rbs = 2
        sorf = sd.detect_start_stop(self.example.tas, seq, args)
        self.assertListEqual(sorf, [{'strand': '+', 'type': 'intergenic', 'starts': ['13'],
                                     'print': False, 'seq': 'ATGCCATTA', 'ends': ['21'],
                                     'end': 21, 'start': 13, 'rbs': [2], 'strain': 'aaa'}])
        seq = {"aaa": "TTAAAGGCATTATCCTCCTA"}
        self.example.tas[0].strand = "-"
        sorf = sd.detect_start_stop(self.example.tas, seq, args)
        self.assertListEqual(sorf, [{'end': 10, 'starts': ['2'], 'strain': 'aaa', 'ends': ['10'],
                                     'type': 'intergenic', 'print': False, 'seq': 'TAAAGGCAT',
                                     'rbs': [19], 'strand': '-', 'start': 2}])
        self.example.tas[0].strand = "+"

    def test_read_data(self):
        inter = os.path.join(self.test_folder, "inter")
        fasta = os.path.join(self.test_folder, "fa")
        gen_file(inter, self.example.inter)
        gen_file(fasta, ">aaa\nATATACCGATC")
        inters, tsss, srnas, seq = sd.read_data(inter, None, None, fasta, True)
        self.assertEqual(inters[0].start, 2)
        self.assertDictEqual(seq, {'aaa': 'ATATACCGATC'})

    def test_check_tss(self):
        sorf = {"strain": "aaa", "strand": "+", "start": 2, "end": 6,
                "starts": [str(2)], "ends": [str(10)], "seq": "ATGTA",
                "type": "3utr", "print": False, "rbs": [1], "with_TSS": []}
        checks = {"start": False, "rbs": False, "import": False}
        sd.check_tss(sorf, self.example.tsss[0], 300, checks)
        self.assertDictEqual(checks, {'start': True, 'rbs': [1], 'import': True})

    def test_compare_sorf_tss(self):
        sorfs = [{"strain": "aaa", "strand": "+", "start": 2, "end": 6,
                 "starts": [str(2)], "ends": [str(10)], "seq": "ATGTA",
                 "type": "3utr", "print": False, "rbs": [1]}]
        args = self.mock_args.mock()
        args.utr_length = 300
        args.noafter_tss = False
        args.no_tss = False
        sorfs_all, sorfs_best = sd.compare_sorf_tss(sorfs, self.example.tsss, "tss", args)
        self.assertListEqual(sorfs_all, [{'print': False, 'ends': ['10'], 'strand': '+',
                                          'end': 6, 'type': '3utr', 'starts': ['2'], 'seq': 'ATGTA',
                                          'strain': 'aaa', 'start': 2, 'rbs': [1],
                                          'start_TSS': '1+', 'with_TSS': ['TSS_1+']}])
        self.assertListEqual(sorfs_best, [{'print': False, 'ends': ['10'], 'strand': '+',
                                           'end': 6, 'type': '3utr', 'starts': ['2'], 'seq': 'ATGTA',
                                           'strain': 'aaa', 'start': 2, 'rbs': [1],
                                           'with_TSS': ['TSS_1+'], 'start_TSS': '1+'}])

    def test_compare_sorf_srna(self):
        sorfs = [{"strain": "aaa", "strand": "+", "start": 2, "end": 6,
                 "starts": [str(2)], "ends": [str(10)], "seq": "ATGTA",
                 "type": "3utr", "print": False, "rbs": [1]}]
        sd.compare_sorf_srna(sorfs, self.example.srnas, "test")
        self.assertListEqual(sorfs, [{'print': False, 'starts': ['2'], 'seq': 'ATGTA', 'strand': '+',
                                      'srna': ['srna0:5-8_f'], 'end': 6, 'rbs': [1], 'ends': ['10'],
                                      'start': 2, 'strain': 'aaa', 'type': '3utr'}])

    def test_import_overlap(self):
        sorf1 = {"strain": "aaa", "strand": "+", "start": 2, "end": 6,
                 "starts": [str(2)], "ends": [str(10)], "seq": "ATGTA",
                 "type": "3utr", "print": False, "rbs": [1], "start_TSS": "1"}
        sorf2 = {"strain": "aaa", "strand": "+", "start": 5, "end": 15,
                 "starts": [str(5)], "ends": [str(15)], "seq": "ATGTA",
                 "type": "3utr", "print": False, "rbs": [2], "start_TSS": "2"}
        final = {"strain": "aaa", "strand": "+", "start": 2, "end": 6,
                 "starts": [str(2)], "ends": [str(10)], "seq": "ATGTA",
                 "type": "3utr", "print": False, "rbs": [1], "start_TSS": "1"}
        sd.import_overlap(sorf2, final, sorf1, True)
        self.assertDictEqual(final, {'end': 15, 'candidate': ['2-6_TSS:1_RBS:1', '5-15_TSS:2_RBS:2'],
                                     'start': 2, 'rbs': [1, 2], 'strand': '+', 'strain': 'aaa',
                                     'print': False, 'seq': 'ATGTA', 'ends': ['10', '15'],
                                     'start_TSS': '1', 'type': '3utr', 'starts': ['2', '5']})

    def test_merge(self):
        seq = {"aaa": "TAGGAGGCCGCTATGCCATTA"}
        sorfs = [{"strain": "aaa", "strand": "+", "start": 2, "end": 6,
                  "starts": [str(2)], "ends": [str(10)], "seq": "ATGTA",
                  "type": "3utr", "print": False, "rbs": [1], "start_TSS": "1"},
                 {"strain": "aaa", "strand": "+", "start": 5, "end": 15,
                  "starts": [str(5)], "ends": [str(15)], "seq": "ATGTA",
                  "type": "3utr", "print": False, "rbs": [2], "start_TSS": "2"}]
        finals = sd.merge(sorfs, seq)
        self.assertDictEqual(finals[0], {'rbs': [1], 'end': 6, 'starts': ['2'], 'strand': '+',
                                         'start_TSS': '1', 'type': '3utr', 'start': 2, 'seq': 'AGGAG',
                                         'candidate': ['2-6_TSS:1_RBS:1'], 'ends': ['10', '6'],
                                         'strain': 'aaa'})

    def test_assign_utr_cutoff(self):
        coverages = {"3utr": "median", "5utr": 20, "interCDS": 11, "intergenic": 59}
        medians = {"median": 50, "mean": 20}
        cutoff =sd.assign_utr_cutoff(coverages, "3utr", medians)
        self.assertEqual(cutoff, 50)

    def test_get_cutoff(self):
        sorf = {"strain": "aaa", "strand": "+", "start": 2, "end": 6,
                "starts": [str(2)], "ends": [str(10)], "seq": "ATGTA",
                "type": "3utr", "print": False, "rbs": [1], "start_TSS": "1"}
        coverages = {"3utr": "median", "5utr": 20, "interCDS": 11, "intergenic": 59}
        medians = {"aaa": {"3utr": {"track_1": {"median": 50, "mean": 20}}}}
        cutoff = sd.get_cutoff(sorf, "track_1", coverages, medians)
        self.assertEqual(cutoff, 50)

    def test_get_attribute(self):
        sorf = {"strain": "aaa", "strand": "+", "start": 2, "end": 6,
                "starts": [str(2)], "ends": [str(10)], "seq": "ATGTA",
                "type": "3utr", "print": False, "rbs": ["1"], "start_TSS": "1",
                "with_TSS": "NA", "srna": "NA", "shift": 1}
        string = sd.get_attribute(1, "sORF_1", "4", sorf, "utr")
        self.assertEqual(string, "ID=sorf1;Name=sORF_sORF_1;start_TSS=4;with_TSS=N&A;sORF_type=3utr;sRNA=N&A;RBS=1;frame_shift=1")

    def test_print_file(self):
        out_g = StringIO()
        out_t = StringIO()
        sorf = {"strain": "aaa", "strand": "+", "start": 10, "end": 15,
                "starts": [str(10)], "ends": [str(15)], "seq": "ATGTA",
                "type": "3utr", "print": False, "rbs": ["3"], "start_TSS": "1",
                "with_TSS": ["NA"], "srna": ["NA"], "candidate": ["AAA"], "shift": 1}
        sorf_datas = {"best": 20, "high": 50, "low": 10, "start": 1,
                      "end": 10, "track": "track_1", "detail": [], "conds": {"frag": "track_1"}}
        args = self.mock_args.mock()
        args.table_best = True
        args.print_all = True
        sd.print_file(sorf, sorf_datas, 1, out_g, out_t, "best", args)
        self.assertEqual(out_g.getvalue(), "aaa\tANNOgesic\tsORF\t10\t15\t.\t+\t.\tID=sorf1;Name=sORF_00001;start_TSS=1;with_TSS=NA;sORF_type=3utr;sRNA=NA;RBS=RBS_3;frame_shift=1\n")
        self.assertEqual(out_t.getvalue(), "aaa\tsORF_00001\t10\t15\t+\t3'UTR_derived\tNA\tRBS_3\t10\t15\tNA\t1\tFragmented\t20\t50\t10\ttrack_1(avg=20;high=50;low=10)\tATGTA\tAAA\n")

    def test_print_table(self):
        out_t = StringIO()
        sorf = {"strain": "aaa", "strand": "+", "start": 2, "end": 6,
                "starts": [str(2)], "ends": [str(10)], "seq": "ATGTA",
                "type": "3utr", "print": False, "rbs": ["1"], "start_TSS": "1",
                "with_TSS": ["NA"], "srna": ["NA"], "candidate": ["AAA"], "shift": 1}
        sorf_datas = {"best": 20, "high": 50, "low": 10, "start": 1,
                      "end": 10, "track": "track_1", "detail": [], "conds": {"frag": "track_1"}}
        args = self.mock_args.mock()
        args.table_best = True
        args.print_all = True
        sd.print_table(out_t, sorf, "test", "3utr", "frag", sorf_datas, args)
        self.assertEqual(out_t.getvalue(), "aaa\tsORF_test\t2\t6\t+\t3utr\tNA\t1\t2\t10\tNA\t1\tfrag\t20\t50\t10\ttrack_1(avg=20;high=50;low=10)\tATGTA\tAAA\n")

    def test_get_inter_coverage(self):
        inter_covers = {}
        inters = [{"frag": [{"track": "track_1", "avg": 22}]}]
        sd.get_inter_coverage(inters, inter_covers)
        self.assertDictEqual(inter_covers, {'track_1': [22]})

    def test_detect_utr_type(self):
        ta_dict = [{"seq_id": "aaa", "source": "intergenic", "feature": "Transcript", "start": 1,
                    "end": 23, "phase": ".", "strand": "+", "score": "."}]
        attributes_tas = [{"ID": "tran0", "Name": "Transcript_0", "UTR_type": "intergenic"}]
        tas = []
        tas.append(Create_generator(ta_dict[0], attributes_tas[0], "gff"))
        sd.get_coverage = self.mock.mock_get_coverage
        med_inters = {"aaa": {"intergenic": []}}
        sd.detect_utr_type(tas[0], "intergenic", med_inters, "wigs", "+", "test")
        sd.get_coverage = get_coverage
        self.assertDictEqual(med_inters, {'aaa': {'intergenic': ["2"]}})

    def test_median_score(self):
        num = sd.median_score([1, 3, 11, 42, 2, 32, 111], "p_0.5")
        self.assertEqual(num, 11)

    def test_mean_score(self):
        num = sd.mean_score([1, 3, 11, 42, 2, 32, 111])
        self.assertEqual(num, 28.857142857142858)

    def test_validate_tss(self):
        sorf = {"strain": "aaa", "strand": "+", "start": 2, "end": 6,
                "starts": [str(2)], "ends": [str(10)], "seq": "ATGTA",
                "type": "3utr", "print": False, "rbs": ["1"], "start_TSS": "3",
                "with_TSS": ["TSS_3+"], "srna": ["NA"], "candidate": ["AAA"]}
        datas = sd.validate_tss([2], [6], sorf, 300)
        self.assertEqual(datas, (['TSS_3+'], 'NA'))

    def test_validate_srna(self):
        sorf = {"strain": "aaa", "strand": "+", "start": 2, "end": 6,
                "starts": [str(2)], "ends": [str(10)], "seq": "ATGTA",
                "type": "3utr", "print": False, "rbs": ["1"], "start_TSS": "1",
                "with_TSS": ["TSS_3+"], "srna": ["sRNA:2-5_+"], "candidate": ["AAA"]}
        srnas = sd.validate_srna([2], [6], sorf)
        self.assertListEqual(srnas, ['sRNA:2-5_+'])

    def test_get_best(self):
        sorfs = [{"strain": "aaa", "strand": "+", "start": 2, "end": 6,
                 "starts": [str(2)], "ends": [str(10)], "seq": "ATGTA",
                 "type": "3utr", "print": False, "rbs": ["1"], "start_TSS": "1",
                 "with_TSS": ["TSS_3+"], "srna": ["sRNA:2-5_+"], "candidate": ["2-6_TSS:3_RBS:1"]}]
        args = self.mock_args.mock()
        args.table_best = True
        args.no_srna = True
        args.utr_length = 300
        data = sd.get_best(sorfs, "tss", "srna", args)
        self.assertListEqual(data, [{'type': '3utr', 'strand': '+', 'print': False,
                                     'with_TSS': ['TSS_3+'], 'starts': ['2'], 'start': 2,
                                     'srna': ['sRNA:2-5_+'], 'rbs': ['1'], 'end': 6, 'seq': 'ATGTA',
                                     'start_TSS': '1', 'strain': 'aaa', 'ends': ['10'],
                                     'candidate': ['2-6_TSS:3_RBS:1']}])
    def test_coverage_and_output(self):
        out_t = StringIO()
        out_g = StringIO()
        sd.get_coverage = self.mock.mock_get_coverage
        sd.replicate_comparison = self.mock.mock_replicate_comparison
        sorfs = [{"strain": "aaa", "strand": "+", "start": 10, "end": 15,
                 "starts": [str(10)], "ends": [str(15)], "seq": "ATGTA",
                 "type": "3utr", "print": False, "rbs": [1], "start_TSS": "1",
                 "with_TSS": ["TSS_3+"], "srna": ["sRNA:2-5_+"], "candidate": ["2-6_TSS:3_RBS:1"]}]
        seq = {"aaa": "TAGGAGGCCGCTATGCCATTA"}
        wigs = {"forward": "wigs_f", "reverse": "wigs_r"}
        args = self.mock_args.mock()
        args.print_all = True
        args.min_rbs = 0
        args.max_rbs = 20
        args.min_len = 0
        args.max_len = 300
        args.table_best = True
        sd.coverage_and_output(sorfs, "median", wigs, out_g, out_t,
                               "best", seq, "cover", args, "texs")
        sd.get_coverage = copy.deepcopy(get_coverage)
        self.assertEqual(out_g.getvalue(), "##gff-version 3\naaa\tANNOgesic\tsORF\t10\t15\t.\t+\t.\tID=sorf0;Name=sORF_00000;start_TSS=1;with_TSS=TSS_3+;sORF_type=3utr;sRNA=NA;RBS=RBS_1;frame_shift=1\n")
        self.assertEqual(out_t.getvalue().split("\n")[1], "aaa\tsORF_00000\t10\t15\t+\t3'UTR_derived\tTSS_3+\tRBS_1\t10\t15\tNA\t1\tFragmented\t20\t50\t10\ttrack_1(avg=20;high=50;low=10)\tGCTATG\t10-15_TSS:3+_RBS:1")

    def test_detect_inter_type(self):
        inter_dict = [{"seq_id": "aaa", "source": "UTR_derived", "feature": "Transcript", "start": 1,
                       "end": 23, "phase": ".", "strand": "+", "score": "."}]
        attributes_inter = [{"ID": "tran0", "Name": "Transcript_0", "UTR_type": "3utr"}]
        inters = []
        inters.append(Create_generator(inter_dict[0], attributes_inter[0], "gff"))
        sd.get_coverage = self.mock.mock_get_coverage
        wigs = {"forward": "wigs_f", "reverse": "wigs_r"}
        data = sd.detect_inter_type(inters, wigs, "test")
        self.assertDictEqual(data, {'aaa': {'interCDS': [], '5utr': [], '3utr': ['2']}})
        sd.get_coverage = copy.deepcopy(get_coverage)

    def test_set_median(self):
        mediandict = {}
        covers = {"aaa": {"3utr": {"track_1": [1, 3, 4, 2, 55]}}}
        coverages = {"3utr": "p_0.5", "5utr": "p_0.5", "interCDS": "n_100"}
        sd.set_median(covers, mediandict, coverages)
        self.assertDictEqual(mediandict, {'aaa': {'5utr': {}, 'interCDS': {}, '3utr': {'track_1': {'median': 3}}}})

    def test_compute_candidate_best(self):
        sorfs = [{"strain": "aaa", "strand": "+", "start": 2, "end": 6,
                 "starts": [str(2)], "ends": [str(10)], "seq": "ATGTA",
                 "type": "3utr", "print": False, "rbs": ["1"], "start_TSS": "1",
                 "with_TSS": ["TSS_3+"], "srna": ["sRNA:2-5_+"]}]
        sd.compute_candidate_best(sorfs)
        self.assertListEqual(sorfs, [{'starts': ['2'], 'seq': 'ATGTA', 'strain': 'aaa',
                                      'ends': ['10'], 'print': False, 'rbs': ['1'], 'type': '3utr',
                                      'end': 6, 'start': 2, 'srna': ['sRNA:2-5_+'],
                                      'candidate': ['2-6_TSS:1_RBS:1'], 'start_TSS': '1',
                                      'strand': '+', 'with_TSS': ['TSS_3+']}])

    def test_sorf_detection(self):
        fasta = os.path.join(self.fasta, "fasta")
        gen_file(fasta, ">aaa\nTAGGAGGCCGCTATGCCATTA")
        srna_gff = os.path.join(self.gff, "srna.gff")
        inter_gff = os.path.join(self.gff, "inter.gff")
        tss_file = os.path.join(self.gff, "tss.gff")
        sd.get_coverage = self.mock.mock_get_coverage
        sd.read_libs = self.mock.mock_read_libs
        sd.read_wig = self.mock.mock_read_wig
        sd.get_inter_coverage = self.mock.mock_get_inter_coverage
        gen_file(srna_gff, self.example.srna)
        gen_file(inter_gff, self.example.inter)
        gen_file(tss_file, self.example.tss)
        args = self.mock_args.mock()
        args.start_codon = ["ATG"]
        args.stop_codon = ["TTA"]
        args.cutoff_5utr = "p_0.5"
        args.cutoff_intercds = "n_20"
        args.cutoff_3utr = "n_11"
        args.cutoff_inter = 50
        args.cutoff_anti = 50
        args.libs = ["frag:frag:1:a:+"]
        args.merge_wigs = "wig_folder"
        args.utr_detect = True
        args.background = 10
        args.print_all = True
        sd.sorf_detection(fasta, srna_gff, inter_gff, tss_file, "wig_f_file",
                          "wig_r_file", "test_folder/test", args)
        sd.get_coverage = copy.deepcopy(get_coverage)
        sd.replicate_comparison = self.mock.mock_replicate_comparison
        self.assertTrue(os.path.exists("test_folder/test_all.csv"))
        self.assertTrue(os.path.exists("test_folder/test_all.gff"))
        self.assertTrue(os.path.exists("test_folder/test_best.csv"))
        self.assertTrue(os.path.exists("test_folder/test_best.gff"))


class Example(object):
    inter = """aaa	UTR_derived	sORF	2	6	.	+	.	ID=inter0;Name=inter_00000;UTR_type=3utr"""
    srna = """aaa	UTR_derived	sRNA	5	8	.	+	.	ID=srna0;Name=srna_00000;UTR_type=3utr"""
    tss = """aaa	tsspredator	TSS	1	1	.	+	.	ID=tss0;Name=TSS_00000"""
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
    ta_dict = [{"seq_id": "aaa", "source": "intergenic", "feature": "Transcript", "start": 1,
                "end": 23, "phase": ".", "strand": "+", "score": "."}]
    attributes_tas = [{"ID": "tran0", "Name": "Transcript_0"}]
    tas = []
    tas.append(Create_generator(ta_dict[0], attributes_tas[0], "gff"))
    tss_dict = [{"seq_id": "aaa", "source": "tsspredator", "feature": "TSS", "start": 1,
                "end": 1, "phase": ".", "strand": "+", "score": "."}]
    attributes_tss = [{"ID": "tss0", "Name": "TSS_0"}]
    tsss = []
    tsss.append(Create_generator(tss_dict[0], attributes_tss[0], "gff"))
    srna_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 5,
                "end": 8, "phase": ".", "strand": "+", "score": "."}]
    attributes_srna = [{"ID": "srna0", "Name": "sRNA_0"}]
    srnas = []
    srnas.append(Create_generator(srna_dict[0], attributes_srna[0], "gff"))
if __name__ == "__main__":
    unittest.main()

