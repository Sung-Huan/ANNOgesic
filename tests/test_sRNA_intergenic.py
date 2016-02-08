import sys
import os
import unittest
import shutil
import copy
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data, gen_file, extract_info
import annogesiclib.sRNA_intergenic as si

get_coverage = si.get_coverage

class Mock_func(object):

    def __init__(self):
        self.example = Example()

    def mock_read_libs(self, input_libs, wig_folder):
        return None, None

    def mock_read_wig(self, wig_file, strand, libs):
        return self.example.wigs

    def mock_get_coverage(self, start, end, strain, wigs, strand, ta, nums, tss, output,
                 template_texs, out_table, cutoff_coverage, tex_notex,
                 replicates, decrease, fuzzy_end, table_best, tolerance):
        pass

    def mock_coverage_comparison(self, cover, cover_sets, poss, check, strand):
        pass

    def mock_replicate_comparison(self, srna_covers, template_texs,
                                  strand, cutoff_coverage, tex_notex,
                                  replicates, type_, test1, test2, test3, tex, notex):
        return {"best": 40, "high": 50, "low": 10, "pos": 5,
                "conds": {"cond1": "test1"}, "detail": None}

class TestsRNAIntergenic(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.mock = Mock_func()
        self.test_folder = "test_folder"
        self.wig_folder = "test_folder/wigs"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
            os.mkdir(self.wig_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_read_data(self):
        gff_file = os.path.join(self.test_folder, "anno.gff")
        tran_file = os.path.join(self.test_folder, "tran.gff")
        pro_file = os.path.join(self.test_folder, "pro.gff")
        gen_file(gff_file, self.example.gff_file)
        gen_file(tran_file, self.example.gff_file)
        gen_file(pro_file, self.example.gff_file)
        nums, cdss, tas, pros, genes = si.read_data(gff_file, tran_file, pro_file, True)
        self.assertDictEqual(nums, {'ta': 3, 'cds': 3, 'pro': 3, 'uni': 0} )
        self.assertEqual(cdss[0].start, 140)
        self.assertEqual(tas[0].start, 140)
        self.assertEqual(pros[0].start, 140)

    def test_read_tss(self):
        tss_file = os.path.join(self.test_folder, "tss.gff")
        gen_file(tss_file, self.example.gff_file)
        tsss, num_tss = si.read_tss(tss_file, "cutoff_coverage", True)
        self.assertEqual(tsss[0].start, 140)

    def test_compare_ta_cds(self):
        detects = {"overlap": False}
        si.compare_ta_cds(self.example.gffs, self.example.tas[0], detects)
        self.assertDictEqual(detects, {'overlap': True})

    def test_compare_ta_tss(self):
        out_table = StringIO()
        nums = {'pro': 3, 'tss': 3, 'uni': 0, 'cds': 3, 'ta': 3}
        output = StringIO()
        detects = {"overlap": False, "uni_with_tss": False}
        si.get_coverage = self.mock.mock_get_coverage
        si.compare_ta_tss(10, 2, 15, 30, 300, "", nums,
                   self.example.tas[0], self.example.tsss[0], output, out_table, "texs", detects, 50,
                   "cutoff", "texnotex", "rep", "decrease",
                   3, True, 5, "texs", 20)
        self.assertEqual(output.getvalue(), "aaa\tANNOgesic\tsRNA\t10\t15\t.\t+\t.\tID=srna0;Name=sRNA_00000;sRNA_type=in_CDS;with_TSS=TSS:170_+\n")
        self.assertEqual(out_table.getvalue(), "aaa\t00000\t10\t15\t+\tNA\tNA\tNA\tNA\tNA\tTSS:170_+\n")
        si.get_coverage = get_coverage

    def test_print_file(self):
        string = "aaa\tintergenic\tsRNA\t10\t15\t.\t+\t."
        nums = {'pro': 3, 'tss': 3, 'uni': 0, 'cds': 3, 'ta': 3}
        out_table = StringIO()
        output = StringIO()
        srna_datas = {"high": 20, "low": 5, "best": 13, "conds": {"cond1": "test1"},
                      "detail": [{"track": "test1", "high": 30, "low": 10, "avg": 15},
                                 {"track": "test2", "high": 25, "low": 13, "avg": 20}]}
        si.print_file(string, nums, "TSS_160+", output, out_table, srna_datas, False, "intergenic")
        self.assertEqual(out_table.getvalue(), "aaa\t00000\t10\t15\t+\tcond1\ttest1\t13\t20\t5\tTSS_160+\ttest1(avg=15;high=30;low=10);test2(avg=20;high=25;low=13)\n")
        self.assertEqual(output.getvalue(), "aaa\tintergenic\tsRNA\t10\t15\t.\t+\t.\tID=srna0;Name=sRNA_00000;sRNA_type=intergenic;with_TSS=TSS_160+;best_avg_coverage=13;best_high_coverage=20;best_low_coverage=5\n")

    def test_detect_include_tss(self):
        si.get_coverage = self.mock.mock_get_coverage
        nums = {'pro': 3, 'tss': 3, 'uni': 0, 'cds': 3, 'ta': 3}
        out_table = StringIO()
        output = StringIO()
        detects = {"overlap": False, "uni_with_tss": False}
        coverage = {"primary": 0, "secondary": 0, "internal": 0, "antisense": 50, "orphan": 10}
        si.detect_include_tss(nums, self.example.tsss, self.example.tas[0],
                              "", "", output, out_table,
                              "template_texs", 30, 100, detects, 5,
                              coverage, "tex_notex", "rep", "decrease",
                              5, True, 5, "texs", coverage)
        si.get_coverage = get_coverage
        self.assertEqual(output.getvalue(), "aaa\tintergenic\tsRNA\t180\t230\t.\t+\t.\tID=srna0;Name=sRNA_00000;sRNA_type=in_CDS\n")
        self.assertEqual(out_table.getvalue(), "aaa\t00000\t180\t230\t+\tNA\tNA\tNA\tNA\tNA\tFalse\n")

    def test_get_differential_cover(self):
        checks = {"detect_diff": True, "first": True}
        cover_sets = {"diff": 30, "low": 5, "high": 35}
        cover = {"coverage": 20, "pos": 80}
        poss = {"stop_point": 100}
        si.get_differential_cover(10, 0, checks, cover_sets,
                                  poss, cover, 200)
        self.assertDictEqual(cover_sets, {'diff': 20, 'low': 20, 'high': 35})
        cover = {"coverage": 50, "pos": 80}
        poss = {"stop_point": 100}
        num = 20
        fuzzy_end = 20
        si.get_differential_cover(fuzzy_end, num, checks, cover_sets,
                                  poss, cover, 200)
        self.assertDictEqual(poss, {"stop_point": 80})

    def test_check_coverage_pos(self):
        si.coverage_comparison = self.mock.mock_coverage_comparison
        cover_sets = {"low": 20, "high":30, "total": 90, "diff": 50}
        poss = {"high": 20, "low": 70, "stop_point": 70}
        tmps = {"total": 0, "toler": 10, "pos": 0}
        checks = {"detect_diff": True, "first": True}
        cover = {"coverage": 50, "pos": 80}
        detect = si.check_coverage_pos(30, 100, cover, 80, tmps, cover_sets,
                                       checks, poss, "+", 5)
        self.assertEqual(detect, (False, {'total': 0, 'pos': 0, 'toler': 0}))
        self.assertDictEqual(poss, {'high': 20, 'stop_point': 70, 'low': 70})
        detect = si.check_coverage_pos(30, 50, cover, 80, tmps, cover_sets,
                                       checks, poss, "+", 5)
        self.assertTrue(detect)

    def test_get_best(self):
        datas = si.get_best(self.example.wigs, "aaa", "+", 2, 23, "normal", 50,
                            10, 5, 5)
        self.assertDictEqual(datas, {'frag_1': [{'avg': 93.0, 'type': 'frag', 'track': 'track_1', 'low': -1, 'pos': 24, 'high': -1}]})

    def test_get_attribute_string(self):
        srna_datas = {'best': 23, 'low': 20, 'high': 35}
        data = si.get_attribute_string(srna_datas, "TSS_100+;Cleavage_150+", 1, "sRNA_00001", "3utr")
        self.assertEqual(data, "ID=srna1;Name=sRNA_sRNA_00001;sRNA_type=3utr;with_TSS=TSS_100+;end_cleavage=Cleavage_150+;best_avg_coverage=23;best_high_coverage=35;best_low_coverage=20")

    def test_check_pro(self):
        si.replicate_comparison = self.mock.mock_replicate_comparison
        srna_datas = {"pos": 50}
        texs = {"track_1@AND@track_2"}
        pro_pos, new_srna_datas, detect_pro = si.check_pro(
                                              self.example.pros, self.example.tas[0],
                                              20, 70, 30, 300, srna_datas, "within",
                                              5, self.example.wigs, 50, 5, "template_texs",
                                              "tex_notex", "replicates", 5, texs, 20)
        self.assertEqual(pro_pos, 190)
        self.assertDictEqual(new_srna_datas, {'best': 40, 'high': 50, 'low': 10, "pos": 5,
                                              "conds": {"cond1": "test1"}, "detail": None})
        self.assertEqual(detect_pro, "Cleavage:190_+")

    def test_exchange_to_pro(self):
        nums = {'pro': 3, 'tss': 3, 'uni': 0, 'cds': 3, 'ta': 3}
        out_table = StringIO()
        output = StringIO()
        srna_datas = {"pos": 50, "best": 10, "high": 12}
        si.replicate_comparison = self.mock.mock_replicate_comparison
        detect, srna_datas, pro = si.exchange_to_pro(srna_datas, 30, 300,
                                  self.example.pros, self.example.tas[0],
                                  20, 70, nums, 10, output, out_table, True, self.example.wigs,
                                  50, 5, "template_texs", "tex_notex", "replicates", 5, "texs", 20)
        self.assertTrue(detect)
        self.assertDictEqual(srna_datas, {'best': 40, 'high': 50, 'low': 10, 'pos': 190,
                                          "conds": {"cond1": "test1"}, "detail": None})
        self.assertEqual(pro, "Cleavage:190_+")

    def test_get_tss_type(self):
        coverage = {"primary": 0, "secondary": 0, "internal": 0, "antisense": 50, "orphan": 10}
        cover = si.get_tss_type(self.example.tsss[0], coverage)
        self.assertEqual(cover, 10)

    def test_detect_wig_pos(self):
        si.replicate_comparison = self.mock.mock_replicate_comparison
        nums = {'pro': 3, 'tss': 3, 'uni': 0, 'cds': 3, 'ta': 3}
        out_table = StringIO()
        output = StringIO()
        si.detect_wig_pos(self.example.wigs, self.example.tas[0], 20, 70, nums, output, "TSS_160+",
                          "template_texs", out_table, 10, 30, 300, 50,
                          5, True, "tex_notex", "replicates", self.example.pros, 5, "texs", 20)
        self.assertEqual(output.getvalue(), "aaa\tANNOgesic\tsRNA\t20\t190\t.\t+\t.\tID=srna0;Name=sRNA_00000;sRNA_type=in_CDS;with_TSS=TSS_160+;end_cleavage=Cleavage:190_+;best_avg_coverage=40;best_high_coverage=50;best_low_coverage=10\n")
        self.assertEqual(out_table.getvalue(), "aaa\t00000\t20\t190\t+\tcond1\ttest1\t40\t50\t10\t\n")

    def test_detect_longer(self):
        si.replicate_comparison = self.mock.mock_replicate_comparison
        si.coverage_comparison = self.mock.mock_coverage_comparison
        nums = {'pro': 3, 'tss': 3, 'uni': 0, 'cds': 3, 'ta': 3}
        out_table = StringIO()
        output = StringIO()
        detects = {"overlap": False, "uni_with_tss": False}
        coverage = {"primary": 0, "secondary": 0, "internal": 0, "antisense": 50, "orphan": 10}
        si.detect_longer(self.example.tsss, self.example.tas[0], nums, output, "", "",
                         "template_texs", out_table, 20, coverage, "tex_notex", "replicates",
                         50, 5, True, 30, 100, detects, self.example.pros, 5, "texs", coverage)
        self.assertEqual(output.getvalue(), "aaa\tANNOgesic\tsRNA\t170\t230\t.\t+\t.\tID=srna0;Name=sRNA_00000;sRNA_type=in_CDS;with_TSS=TSS:170_+\n")
        self.assertEqual(out_table.getvalue(), "aaa\t00000\t170\t230\t+\tNA\tNA\tNA\tNA\tNA\tTSS:170_+\n")

    def test_get_proper_tss(self):
        tss_file = os.path.join(self.test_folder, "tss.gff")
        gen_file(tss_file, self.example.gff_file)
        coverage = {"primary": 0, "secondary": 0, "internal": 0, "antisense": 50, "orphan": 10}
        tsss, num_tss = si.get_proper_tss(tss_file, coverage)
        self.assertEqual(tsss[0].start, 140)

    def test_check_srna_condition(self):
        si.replicate_comparison = self.mock.mock_replicate_comparison
        si.coverage_comparison = self.mock.mock_coverage_comparison
        nums = {'pro': 3, 'tss': 3, 'uni': 0, 'cds': 3, 'ta': 3}
        out_table = StringIO()
        output = StringIO()
        detects = {"overlap": False, "uni_with_tss": False}
        notex = {"primary": 0, "secondary": 0, "internal": 0, "antisense": 30, "orphan": 10}
        coverage = {"primary": 0, "secondary": 0, "internal": 0, "antisense": 50, "orphan": 10}
        si.check_srna_condition(self.example.tas[0], 30, 300, self.example.tsss, self.example.pros,
                                "", "", nums, output, out_table, "texs", 20, detects,
                                coverage, "tex_notex", "replicates", 50, 5, True, 5, notex)
        self.assertEqual(output.getvalue(), "aaa\tANNOgesic\tsRNA\t170\t230\t.\t+\t.\tID=srna0;Name=sRNA_00000;sRNA_type=intergenic;with_TSS=TSS:170_+\n")
        self.assertEqual(out_table.getvalue(), "aaa\t00000\t170\t230\t+\tNA\tNA\tNA\tNA\tNA\tTSS:170_+\n")

    def test_intergenic_srna(self):
        si.read_libs = self.mock.mock_read_libs
        si.read_wig = self.mock.mock_read_wig
        gff_file = os.path.join(self.test_folder, "aaa.gff")
        tss_file = os.path.join(self.test_folder, "aaa_TSS.gff")
        tran_file = os.path.join(self.test_folder, "aaa_tran.gff")
        pro_file = os.path.join(self.test_folder, "aaa_processing.gff")
        wig_f_file = os.path.join(self.wig_folder, "wig_f.wig")
        wig_r_file = os.path.join(self.wig_folder, "wig_r.wig")
        gen_file(gff_file, self.example.gff_file)
        gen_file(tss_file, self.example.gff_file)
        gen_file(tran_file, self.example.gff_file)
        gen_file(pro_file, self.example.gff_file)
        output_file = os.path.join(self.test_folder, "output")
        output_table = os.path.join(self.test_folder, "table")
        coverage = [0, 0, 0, 50, 10]
        si.replicate_comparison = self.mock.mock_replicate_comparison
        si.coverage_comparison = self.mock.mock_coverage_comparison
        si.intergenic_srna(gff_file, tran_file, tss_file, pro_file, 20, 300,
                           30, wig_f_file, wig_r_file, self.wig_folder, "input_libs",
                           "tex_notex", "replicates", output_file, output_table,
                           True, 50, 5, coverage, 5, False, True, False, "aaa", self.test_folder,
                           "frag", None)
        self.assertTrue(os.path.exists(output_file))
        self.assertTrue(os.path.exists(output_table))
        


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
    ta_dict = [{"seq_id": "aaa", "source": "intergenic", "feature": "Transcript", "start": 180,
                "end": 230, "phase": ".", "strand": "+", "score": "."}]
    attributes_tas = [{"ID": "tran0", "Name": "Transcript_0", "sRNA_type": "intergenic"}]
    tas = []
    tas.append(Create_generator(ta_dict[0], attributes_tas[0], "gff"))
    tss_dict = [{"seq_id": "aaa", "source": "intergenic", "feature": "TSS", "start": 170,
                "end": 170, "phase": ".", "strand": "+", "score": "."}]
    attributes_tsss = [{"ID": "tss0", "Name": "TSS_0", "type": "Orphan"}]
    tsss = []
    tsss.append(Create_generator(tss_dict[0], attributes_tsss[0], "gff"))
    pro_dict = [{"seq_id": "aaa", "source": "intergenic", "feature": "TSS", "start": 190,
                "end": 190, "phase": ".", "strand": "+", "score": "."}]
    attributes_pros = [{"ID": "tss0", "Name": "TSS_0"}]
    pros = []
    pros.append(Create_generator(pro_dict[0], attributes_pros[0], "gff"))
    gff_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 150,
                 "end": 200, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 1230,
                 "end": 1240, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 7100,
                 "end": 9167, "phase": ".", "strand": "-", "score": "."}]
    attributes_gff = [{"ID": "cds0", "Name": "CDS_0", "locus_tag": "AAA_00001"},
                      {"ID": "cds1", "Name": "CDS_1", "locus_tag": "AAA_00002"},
                      {"ID": "cds4", "Name": "CDS_4", "locus_tag": "BBB_00003", "product": "hypothetical protein"}]
    gffs = []
    for index in range(0, 3):
        gffs.append(Create_generator(gff_dict[index], attributes_gff[index], "gff"))
    gff_file = """aaa	RefSeq	CDS	140	160	.	+	.	ID=srna0;Name=sRNA_0
aaa	RefSeq	CDS	230	280	.	+	.	ID=srna1;Name=sRNA_1
aaa	RefSeq	CDS	5166	5266	.	-	.	ID=srna2;Name=sRNA_2"""

if __name__ == "__main__":
    unittest.main()

