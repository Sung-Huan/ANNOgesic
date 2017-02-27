import sys
import os
import unittest
import shutil
import copy
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data, gen_file
import annogesiclib.optimize_TSSpredator as ot
from mock_args_container import MockClass


class Mock_func(object):

    def __init__(self):
        self.example = Example()

    def mock_run_TSSpredator_paralle(self, config_files, tsspredator_path, processes):
        pass

    def mock_convert2gff(self, out_path, gff_files, args, test):
        if not os.path.exists("test_folder/gffs"):
            os.mkdir("test_folder/gffs")
        gen_file("test_folder/gffs/aaa.gff", self.example.gff_file)
        gff_files.append("test_folder/gffs/aaa.gff")

class TestOptimizeTSSpredator(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.mock_args = MockClass()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_initiate(self):
        args = self.mock_args.mock()
        args.height = 0.9
        args.height_reduction = 0.8
        args.factor = 0.9
        args.factor_reduction = 0.8
        args.base_height = 0.01
        args.enrichment = 0.5
        args.processing = 0.5
        max_num, best_para, current_para, indexs = ot.initiate(args)
        self.assertDictEqual(max_num, {'re_factor': 0.8, 'processing': 0.5, 'enrichment': 0.5,
                                       'height': 0.9, 'base_height': 0.01, 're_height': 0.8,
                                       'factor': 0.9})
        self.assertDictEqual(best_para, {'re_factor': 0, 'processing': 0, 'enrichment': 0,
                                         'height': 0, 'base_height': 0, 're_height': 0, 'factor': 0})
        self.assertDictEqual(current_para, {'re_factor': 0, 'processing': 0, 'enrichment': 0,
                                            'height': 0, 'base_height': 0, 're_height': 0,
                                            'factor': 0})
        self.assertDictEqual(indexs, {'step': 0, 'change': False, 'num': 0, 'first': True,
                                      'length': 0, 'exist': False, 'switch': 0, 'extend': False,
                                      'count': 0})

    def test_read_predict_manual_gff(self):
        gff = os.path.join(self.test_folder, "test.gff")
        gen_file(gff, self.example.gff_file)
        args = self.mock_args.mock()
        args.gene_length = 1000
        num, gffs = ot.read_predict_manual_gff(gff, 1000)
        self.assertEqual(num, 1)
        self.assertEqual(gffs[0].start, 633)

    def test_scoring_function(self):
        stat_value = {"tp_rate": 0.8, "fp_rate": 0.0003, "tp": 100, "fp": 3}
        best = {"tp_rate": 0.8, "fp_rate": 0.0005, "tp": 100, "fp": 31, "fn": 45, "missing_ratio": 0.004}
        ot.scoring_function(best, stat_value, self.example.indexs, 1000)
        self.assertTrue(self.example.indexs["change"])
        self.example.indexs["change"] = False
        stat_value = {"tp_rate": 0.8, "fp_rate": 0.0004, "tp": 100, "fp": 13}
        best = {"tp_rate": 0.8, "fp_rate": 0.0003, "tp": 100, "fp": 3}
        ot.scoring_function(best, stat_value, self.example.indexs, 1000)
        self.assertFalse(self.example.indexs["change"])

    def test_load_stat_csv(self):
        stat_file = os.path.join(self.test_folder, "stat.csv")
        gen_file(stat_file, self.example.stat)
        list_num = []
        best_para = {}
        datas = ot.load_stat_csv(self.test_folder, list_num, self.example.best,
                                 best_para, self.example.indexs, 1000, stat_file)
        self.assertEqual(datas[0], 2)
        self.assertDictEqual(datas[1], {'fp': 230.0, 'tp': 789.0, 'missing_ratio': 0.29991126885536823,
                                        'fp_rate': 8.15542105020548e-05, 'tp_rate': 0.7000887311446318,
                                        'fn': 338.0})
        self.assertDictEqual(datas[2], {'processing': 5.2, 'base_height': 0.086, 'factor': 7.6,
                                        're_height': 2.3, 're_factor': 5.5, 'enrichment': 3.1,
                                        'height': 2.4})

    def test_reload_data(self):
        stat_file = os.path.join(self.test_folder, "stat.csv")
        gen_file(stat_file, self.example.stat)
        list_num = []
        best_para = {}
        datas = ot.reload_data(self.test_folder, list_num, self.example.best, best_para, self.example.indexs, 1000, stat_file)
        self.assertDictEqual(datas[0], {'base_height': 0.086, 'processing': 5.2,
                                        'height': 2.4, 'enrichment': 3.1, 're_factor': 5.5,
                                        're_height': 2.3, 'factor': 7.6})
        self.assertDictEqual(datas[1], {'tp_rate': 0.7000887311446318, 'tp': 789.0,
                                        'fn': 338.0, 'fp': 230.0, 'fp_rate': 8.15542105020548e-05,
                                        'missing_ratio': 0.29991126885536823})

    def test_extend_data(self):
        best_para = copy.deepcopy(self.example.best_para)
        current_para = ot.extend_data(self.test_folder, self.example.best, best_para, 100, "aaa")
        self.assertDictEqual(current_para, best_para)

    def test_run_random_part(self):
        list_num = []
        current_para = copy.deepcopy(self.example.ref_para)
        para = ot.run_random_part(current_para, list_num, self.example.max_nums, 1000, self.example.indexs)
        self.assertTrue(para != self.example.ref_para)

    def test_run_large_change_part(self):
        list_num = []
        seeds = {"seed": 0, "pre_seed": []}
        features = {"feature": "r", "pre_feature": ""}
        current_para = copy.deepcopy(self.example.ref_para)
        best_para = copy.deepcopy(self.example.best_para)
        para = ot.run_large_change_part(seeds, features, self.example.indexs, current_para,
                                        self.example.max_nums, best_para, list_num)
        self.assertTrue(para != self.example.ref_para)
        self.assertTrue(para != best_para)

    def test_gen_large_random(self):
        list_num = []
        index_large = {0: "height", 1: "re_height", 2: "factor", 3: "re_factor",
                       4:"base_height", 5: "enrichment", 6: "processing"}
        best_para = copy.deepcopy(self.example.best_para)
        para = ot.gen_large_random(self.example.max_nums, "height", 0.2, list_num, 0.3,
                                   best_para, index_large, self.example.indexs)
        self.assertTrue(para != best_para)
        self.assertTrue(para["height"] > para["re_height"])

    def test_run_small_change_part(self):
        seeds = {"seed": 0, "pre_seed": []}
        features = {"feature": "l", "pre_feature": ""}
        current_para = copy.deepcopy(self.example.ref_para)
        list_num = []
        best_para = copy.deepcopy(self.example.best_para)
        para = ot.run_small_change_part(seeds, features, self.example.indexs, current_para,
                                        best_para, list_num, self.example.max_nums)
        self.assertTrue(para != best_para)

    def test_small_change(self):
        list_num = []
        best_para = copy.deepcopy(self.example.best_para)
        para = ot.small_change(0.9, "height", 0.2, list_num, 0.5, best_para)
        self.assertTrue(para != 0.5)
        self.assertTrue(para > 0.2)

    def test_plus_process(self):
        list_num = []
        actions = {"plus": False, "minus": False}
        best_para = copy.deepcopy(self.example.best_para)
        para = ot.plus_process("height", best_para, 0.9, 0.5, actions, list_num, 0.2)
        self.assertEqual(para, 0.4)

    def test_minus_process(self):
        list_num = []
        actions = {"plus": False, "minus": False}
        best_para = copy.deepcopy(self.example.best_para)
        para = ot.minus_process("height", best_para, 0.9, 0.5, actions, list_num, 0.1)
        self.assertEqual(para, 0.2)

    def test_compare_manual_predict(self):
        out = StringIO()
        manual = os.path.join(self.test_folder, "manual.gff")
        predict = os.path.join(self.test_folder, "predict.gff")
        gen_file(manual, self.example.manual_file)
        gen_file(predict, self.example.gff_file)
        para_list = [copy.deepcopy(self.example.best_para)]
        args = self.mock_args.mock()
        args.manual = manual
        args.cores = 1
        args.gene_length = 2000
        args.cluster = 3
        ot.compare_manual_predict(1000, para_list, [predict], self.test_folder,
                                  out, args, self.example.mans, 3, 2000)
        self.assertEqual(out.getvalue(), "1000\the_0.3_rh_0.2_fa_0.7_rf_0.3_bh_0.0_ef_2.5_pf_3.3\tTP=0\tTP_rate=0.0\tFP=2\tFP_rate=0.00100150225338007\tFN=2\tmissing_ratio=0.6666666666666666\n")

    def test_compute_stat(self):
        list_num = [self.example.best_para]
        best_para = {'re_factor': 0.3, 'processing': 3.3, 'enrichment': 2.5,
                     'height': 0.5, 'base_height': 0.0, 're_height': 0.2, 'factor': 0.7}
        self.example.indexs["change"] = True
        best = {"tp_rate": 0.6, "fp_rate": 0.0025, "tp": 40, "fp": 32, "fn": 45, "missing_ratio": 0.004}
        datas = ot.compute_stat(self.example.best, best, best_para, 1, list_num, self.test_folder, self.example.indexs, "aaa")
        self.assertDictEqual(datas[0], self.example.best_para)
        self.assertDictEqual(datas[1], self.example.best)

    def test_run_tss_and_stat(self):
        list_num = [self.example.best_para]
        seeds = {"seed": 0, "pre_seed": []}
        features = {"feature": "l", "pre_feature": ""}
        best_para = {'re_factor': 0.3, 'processing': 3.3, 'enrichment': 2.5,
                     'height': 0.5, 'base_height': 0.0, 're_height': 0.2, 'factor': 0.7}
        current_para = {'re_factor': 0.3, 'processing': 2.3, 'enrichment': 2.5,
                        'height': 0.5, 'base_height': 0.2, 're_height': 0.2, 'factor': 0.7}
        stat_out = StringIO()
        wig = os.path.join(self.test_folder, "wig")
        fasta = os.path.join(self.test_folder, "aaa.fa")
        gff = os.path.join(self.test_folder, "aaa.gff")
        if not os.path.exists(wig):
            os.mkdir(wig)
        gen_file(fasta, self.example.fasta)
        gen_file(gff, self.example.gff_file)
        output_prefix = ["test_aaa"]
        ot.run_TSSpredator_paralle = Mock_func().mock_run_TSSpredator_paralle
        ot.convert2gff = Mock_func().mock_convert2gff
        args = self.mock_args.mock()
        args.steps = 2000
        args.cores = 1
        args.tsspredator_path = "test"
        args.libs = self.example.libs
        args.program = "TSS"
        args.cluster = 3
        args.utr = 200
        args.replicate = "all_2"
        args.replicate_name = "rep"
        args.project_strain = "aaa"
        args.manual = os.path.join(self.test_folder, "manual.gff")
        args.gene_length = 2000
        gen_file(args.manual, self.example.manual_file)
        datas = ot.run_tss_and_stat(self.example.indexs, list_num, seeds, 0.4, 0.3,
                     self.test_folder, stat_out, best_para, current_para,
                     wig, fasta, gff, self.example.best, 3, args, "aaa", self.example.mans, 2000)
        self.assertFalse(datas[0])

    def test_gen_config(self):
        wig = os.path.join(self.test_folder, "wig")
        if not os.path.exists(wig):
            os.mkdir(wig)
        fasta = os.path.join(self.test_folder, "aaa.fa")
        gff = os.path.join(self.test_folder, "aaa.gff")
        gen_file(fasta, self.example.fasta)
        gen_file(gff, self.example.gff_file)
        args = self.mock_args.mock()
        args.libs = self.example.libs
        args.cores = 1
        args.cluster = 3
        args.program = "TSS"
        args.project_strain = "aaa"
        args.replicate = "all_1"
        args.utr = 200
        args.replicate_name = "test"
        filename = ot.gen_config(self.example.best_para, self.test_folder, 1, wig, fasta, gff, args, "aaa")
        self.assertEqual(filename, "test_folder/config_1.ini")
        data = import_data("test_folder/config_1.ini")
        self.assertEqual("\n".join(data), self.example.config)       

    def test_comparison(self):
        nums = {"overlap": 0, "predict": 0, "manual": 0}
        for index in range(0, 3):
            self.example.mans[index].attributes["print"] = False
            self.example.gffs[index].attributes["print"] = False
        args = self.mock_args.mock()
        args.cluster = 3
        args.gene_length = 2000
        ot.comparison(self.example.mans, self.example.gffs, nums, args, 2000)
        self.assertDictEqual(nums, {'manual': 1, 'predict': 2, 'overlap': 1})

    def test_check_overlap(self):
        nums = {"overlap": 0, "predict": 0, "manual": 0}
        datas = ot.check_overlap(True, None, nums, 2000, self.example.mans[0], self.example.gffs[0], 100)
        self.assertFalse(datas[0])
        self.assertEqual(datas[1], 140)

    def test_print_lib(self):
        libs = [{"condition": 1, "replicate": "a", "wig": "test_1.wig"},
                {"condition": 2, "replicate": "a", "wig": "test_2.wig"}]
        out = StringIO()
        ot.print_lib(2, libs, out, self.test_folder, "aaa", ["a"])
        self.assertEqual(out.getvalue(), "aaa_1a = test_folder/test_1.wig\naaa_2a = test_folder/test_2.wig\n")

    def test_import_lib(self):
        out = StringIO()
        if not os.path.exists(os.path.join(self.test_folder, "wigs")):
            os.mkdir(os.path.join(self.test_folder, "wigs"))
        wig_folder = os.path.join(self.test_folder, "wigs", "tmp")
        if not os.path.exists(wig_folder):
            os.mkdir(wig_folder)
        lib_dict = {"fp": [], "fm": [], "np": [], "nm": []}
        gen_file(os.path.join(wig_folder, "GSM1649587_Hp26695_ML_B1_HS1_-TEX_forward_STRAIN_aaa.wig"), "test")
        gen_file(os.path.join(wig_folder, "GSM1649587_Hp26695_ML_B1_HS1_-TEX_reverse_STRAIN_aaa.wig"), "test")
        gen_file(os.path.join(wig_folder, "GSM1649588_Hp26695_ML_B1_HS1_-TEX_forward_STRAIN_aaa.wig"), "test")
        gen_file(os.path.join(wig_folder, "GSM1649588_Hp26695_ML_B1_HS1_-TEX_reverse_STRAIN_aaa.wig"), "test")
        args = self.mock_args.mock()
        args.project_strain = "aaa"
        args.program = "TSS"
        args.libs = self.example.libs
        lib_num = ot.import_lib(wig_folder, set(), lib_dict, out, "aaa.gff",
                                [], "aaa.fa", args, "aaa")
        self.assertEqual(lib_num, 1)

    def test_optimization_process(self):
        current_para = copy.deepcopy(self.example.ref_para)
        best_ref_para = copy.deepcopy(self.example.best_para)
        list_num = [best_ref_para]
        indexs = copy.deepcopy(self.example.indexs)
        best_para = {'re_factor': 0.3, 'processing': 3.3, 'enrichment': 2.5,
                     'height': 0.6, 'base_height': 0.0, 're_height': 0.2, 'factor': 0.7}
        stat_out = StringIO()
        output_prefix = ["test_1"]
        gen_file(os.path.join(self.test_folder, "manual.gff"), self.example.manual_file)
        if not os.path.exists(os.path.join(self.test_folder, "wigs")):
            os.mkdir(os.path.join(self.test_folder, "wigs"))
        wig_folder = os.path.join(self.test_folder, "wigs", "tmp")
        if not os.path.exists(wig_folder):
            os.mkdir(wig_folder)
        gen_file(os.path.join(wig_folder, "GSM1649587_Hp26695_ML_B1_HS1_-TEX_forward_STRAIN_aaa.wig"), "test")
        gen_file(os.path.join(wig_folder, "GSM1649587_Hp26695_ML_B1_HS1_-TEX_reverse_STRAIN_aaa.wig"), "test")
        gen_file(os.path.join(wig_folder, "GSM1649588_Hp26695_ML_B1_HS1_-TEX_forward_STRAIN_aaa.wig"), "test")
        gen_file(os.path.join(wig_folder, "GSM1649588_Hp26695_ML_B1_HS1_-TEX_reverse_STRAIN_aaa.wig"), "test")
        ot.run_TSSpredator_paralle = Mock_func().mock_run_TSSpredator_paralle
        ot.convert2gff = Mock_func().mock_convert2gff
        args = self.mock_args.mock()
        args.libs = self.example.libs
        args.cores = 1
        args.cluster = 3
        args.program = "TSS"
        args.project_strain = "aaa"
        args.replicate = "all_1"
        args.utr = 200
        args.replicate_name = "test"
        args.steps = 2
        args.tsspredator_path = "test"
        args.gene_length = 2000
        args.manual = os.path.join(self.test_folder, "manual.gff")
        ot.optimization_process(indexs, current_para, list_num, self.example.max_nums, best_para,
                                self.test_folder, stat_out, self.example.best, wig_folder, "aaa.fa", "aaa.gff",
                                2, True, args, "aaa", self.example.mans, 2000)
        self.assertDictEqual(best_para, {'re_height': 0.2, 'factor': 0.7, 'processing': 3.3,
                                         'height': 0.6, 'base_height': 0.0, 're_factor': 0.3,
                                         'enrichment': 2.5})
#        self.assertDictEqual(self.example.best, {'fp': 3, 'missing_ratio': 0.004, 'tp_rate': 0.8,
#                                                 'tp': 100, 'fn': 45, 'fp_rate': 0.0005})

    def test_optimization(self):
        ot.run_TSSpredator_paralle = Mock_func().mock_run_TSSpredator_paralle
        ot.convert2gff = Mock_func().mock_convert2gff
        if not os.path.exists(os.path.join(self.test_folder, "wigs")):
            os.mkdir(os.path.join(self.test_folder, "wigs"))
        wig_folder = os.path.join(self.test_folder, "wigs", "tmp")
        if not os.path.exists(wig_folder):
            os.mkdir(wig_folder)
        fasta = os.path.join(self.test_folder, "aaa.fa")
        gff = os.path.join(self.test_folder, "aaa.gff")
        gen_file(fasta, self.example.fasta)
        gen_file(gff, self.example.gff_file)
        output_prefix = ["test_1"]
        args = self.mock_args.mock()
        args.libs = self.example.libs
        args.cores = 1
        args.cluster = 3
        args.program = "TSS"
        args.project_strain = "aaa"
        args.replicate = "all_1"
        args.utr = 200
        args.steps = 2
        args.gene_length = 2000
        args.height = 0.9
        args.height_reduction = 0.8
        args.factor = 0.9
        args.factor_reduction = 0.8
        args.base_height = 0.01
        args.enrichment = 0.5
        args.processing = 0.5
        args.length = None
        args.replicate_name = "test"
        args.tsspredator_path = "test"
        args.manual = os.path.join(self.test_folder, "manual.gff")
        gen_file(args.manual, self.example.manual_file)
        args.output_folder = self.test_folder
        ot.optimization(wig_folder, fasta, gff, args, args.manual, 2000, "aaa")
        self.assertTrue(os.path.exists(os.path.join(self.test_folder, "optimized_TSSpredator", "stat_aaa.csv")))

class Example(object):

    libs = ["GSM1649587_Hp26695_ML_B1_HS1_-TEX_forward.wig:notex:1:a:+",
            "GSM1649587_Hp26695_ML_B1_HS1_-TEX_reverse.wig:notex:1:a:-",
            "GSM1649588_Hp26695_ML_B1_HS1_+TEX_forward.wig:tex:1:a:+",
            "GSM1649588_Hp26695_ML_B1_HS1_+TEX_reverse.wig:tex:1:a:-"]
    best = {"tp_rate": 0.8, "fp_rate": 0.0005, "tp": 100, "fp": 3, "fn": 45, "missing_ratio": 0.004}
    best_para = {'re_factor': 0.3, 'processing': 3.3, 'enrichment': 2.5,
                 'height': 0.3, 'base_height': 0.0, 're_height': 0.2, 'factor': 0.7}
    max_nums = {'re_factor': 0.9, 'processing': 9.0, 'enrichment': 9.0,
                'height': 0.9, 'base_height': 0.1, 're_height': 0.9, 'factor': 0.9}
    ref_para = {'re_factor': 0.4, 'processing': 0.3, 'enrichment': 2.0,
                'height': 0.5, 'base_height': 0.0, 're_height': 0.1, 'factor': 0.5}
    indexs = {'step': 0, 'change': False, 'num': 0, 'first': True, 'length': 0,
              'exist': False, 'switch': 0, 'extend': False, 'count': 0}
    fasta = """>aaa
AGACTTCCTGATAGTTAAACACATGAGATGTTGGCGTACACACCCGGTGTTTACGTATACGTTACTATGATATTTAGAAAAAACCCGTGTATACGTTCGTGA"""
    gff_file = """NC_000915.1	RefSeq	TSS	633	633	.	-	.	Name=nusB;gene=nusB;locus_tag=HP0001;ID=gene0;Dbxref=GeneID:898756;gbkey=Gene
NC_000915.1	RefSeq	TSS	1105	1105	.	-	.	Name=ribH;gene=ribH;locus_tag=HP0002;ID=gene1;Dbxref=GeneID:898768;gbkey=Gene"""
    manual_file = """NC_000915.1	RefSeq	TSS	633	633	.	-	.	Name=nusB;gene=nusB;locus_tag=HP0001;ID=gene0;Dbxref=GeneID:898756;gbkey=Gene
NC_000915.1	RefSeq	TSS	1125	1125	.	-	.	Name=ribH;gene=ribH;locus_tag=HP0002;ID=gene1;Dbxref=GeneID:898768;gbkey=Gene"""
    stat = """0	he_2.4_rh_2.3_fa_7.6_rf_5.5_bh_0.086_ef_3.1_pf_5.2	TP=789	TP_rate=0.7000887311446318	FP=230	FP_rate=8.15542105020548e-05	FN=338	missing_ratio=0.29991126885536823
1	he_1.4_rh_1.2_fa_7.5_rf_2.5_bh_0.149_ef_5.2_pf_5.0	TP=595	TP_rate=0.5279503105590062	FP=195	FP_rate=6.91437871647856e-05	FN=532	missing_ratio=0.4720496894409938"""
    gffs_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 140,
                  "end": 140, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 30,
                  "end": 30, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 430,
                  "end": 430, "phase": ".", "strand": "-", "score": "."}]
    attributes_gffs = [{"ID": "TSS0", "Name": "TSS_0", "locus_tag": "AAA_00001"},
                       {"ID": "TSS1", "Name": "TSS_1", "locus_tag": "AAA_00002"},
                       {"ID": "TSS2", "Name": "TSS_2", "locus_tag": "AAA_00003"}]
    mans_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 140,
                  "end": 142, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 40,
                  "end": 40, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 5167,
                  "end": 5167, "phase": ".", "strand": "-", "score": "."}]
    attributes_mans = [{"ID": "TSS0", "Name": "TSS_0", "locus_tag": "AAA_00001"},
                       {"ID": "TSS1", "Name": "TSS_1", "locus_tag": "AAA_00002"},
                       {"ID": "TSS2", "Name": "TSS_2", "locus_tag": "AAA_00003"}]
    gffs = []
    mans = []
    for index in range(0, 3):
        gffs.append(Create_generator(gffs_dict[index], attributes_gffs[index], "gff"))
        mans.append(Create_generator(mans_dict[index], attributes_mans[index], "gff"))
    config = """TSSinClusterSelectionMethod = HIGHEST
allowedCompareShift = 1
allowedRepCompareShift = 1
annotation_1 = test_folder/aaa.gff
fivePrimeMinus_1a = test_folder/wig/GSM1649588_Hp26695_ML_B1_HS1_+TEX_reverse.wig
fivePrimePlus_1a = test_folder/wig/GSM1649588_Hp26695_ML_B1_HS1_+TEX_forward.wig
genome_1 = test_folder/aaa.fa
idList = 1
maxASutrLength = 100
maxGapLengthInGene = 500
maxNormalTo5primeFactor = 3.3
maxTSSinClusterDistance = 4
maxUTRlength = 200
min5primeToNormalFactor = 2.5
minCliffFactor = 0.7
minCliffFactorDiscount = 0.3
minCliffHeight = 0.3
minCliffHeightDiscount = 0.2
minNormalHeight = 0.0
minNumRepMatches = 1
minPlateauLength = 0
mode = cond
normPercentile = 0.9
normalMinus_1a = test_folder/wig/GSM1649587_Hp26695_ML_B1_HS1_-TEX_reverse.wig
normalPlus_1a = test_folder/wig/GSM1649587_Hp26695_ML_B1_HS1_-TEX_forward.wig
numReplicates = 1
numberOfDatasets = 1
outputDirectory = test_folder/MasterTable_1
outputPrefix_1 = t
outputPrefix_2 = e
outputPrefix_3 = s
outputPrefix_4 = t
projectName = aaa
superGraphCompatibility = igb
texNormPercentile = 0.5
writeGraphs = 0
writeNocornacFiles = 0"""
    table = """SuperPos	SuperStrand	mapCount	detCount	Genome	detected	enriched	stepHeight	stepFactor	enrichmentFactor	classCount	Pos	Strand	Locus_tag	sRNA/asRNA	Product	UTRlength	GeneLength	Primary	Secondary	Internal	Antisense	Automated	Manual	Putative sRNA	Putative asRNA	Comment	Sequence -50 nt upstream + TSS (51nt)
179	-	1	1	test	1	1	4.45	31.93	8.69	1	179	-	orphan	orphan	NA	NA	0	0	0	0	1	0	0	0		ACCCTTGAATTGAGGGTGTTTTATACCTAAATTTAAAAAATGATGCTATAA
681	-	1	1	test	1	1	4.2	3.0	3.54	2	681	-	HP0001	transcription antitermination protein NusB	48	417	1	0	0	0	1	0	0	0		GATTGAAAGAGCGGGCAGTAAAGCCGGCAATAAGGGCTTTGAAGCGATGAG
681	-	1	1	test	1	1	4.2	3.0	3.54	2	681	-	HP0002	6%2C7-dimethyl-8-ribityllumazine synthase	NA	471	0	0	1	0	1	0	0	0		GATTGAAAGAGCGGGCAGTAAAGCCGGCAATAAGGGCTTTGAAGCGATGAG
1361	-	1	1	test	1	1	24.16	6.31	3.98	2	1361	-	HP0002	6%2C7-dimethyl-8-ribityllumazine synthase	256	471	1	0	0	0	1	0	0	0		TATGGGAATTTAGTGGTGGATATGCGCTCTTTAAAAATCATGCGAGAATTT
1361	-	1	1	test	1	1	24.16	6.31	3.98	2	1361	-	HP0003	2-dehydro-3-deoxyphosphooctonate aldolase	NA	831	0	0	1	0	1	0	0	0		TATGGGAATTTAGTGGTGGATATGCGCTCTTTAAAAATCATGCGAGAATTT
2689	-	1	1	test	1	1	10.15	10.67	11.97	2	2689	-	HP0004	carbonic anhydrase IcfA 92	666	1	0	0	0	1	0	0	1		CAATGCGACACATAATTGCATGAAAGCCCTTTAAAGTGTAAAATAACGCCA
2689	-	1	1	test	1	1	10.15	10.67	11.97	2	2689	-	HP0005	orotidine 5%27-phosphate decarboxylase	NA	684	0	0	0	1	1	0	0	1		CAATGCGACACATAATTGCATGAAAGCCCTTTAAAGTGTAAAATAACGCCA"""







if __name__ == "__main__":
    unittest.main()

