import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file, import_data, extract_info
import annogesiclib.tsspredator as ts
from annogesiclib.tsspredator import TSSpredator


class Mock_func(object):

    def __init__(self):
        self.example = Example()

    def mock_start_to_run(self, tsspredator_path, config_file, out_path, prefix):
        gen_file(os.path.join(out_path, "TSSstatistics.tsv"), "test")        

    def mock_check_orphan(sel, pre_tss, gff, wig_f, wig_r, tmp_tss):
        gen_file(tmp_tss, "test")

    def mock_filter_low_expression(self, gff, manual, wig_f, wig_r, input_libs,
                                   wig_folder, cluster, nt_length, without):
        gen_file("tmp/without_low_expression.gff", self.example.tss_file)
        return 100

    def mock_merge_manual_predict_tss(self, predict, manual, stat_file,
                                      tmp_tss, gff, cluster, length, libs,
                                      wig_path, feature):
        gen_file('tmp_TSS/test_TSS.gff', self.example.tss_file)
        gen_file('stat_compare_TSSpredator_manual_test.csv', "test")

    def mock_filter_tss_pro(self, tss, ref, overlap, cluster):
        gen_file(tss, "test")

    def mock_plot_venn(self, compare_file, feature):
        pass

    def mock_stat_tsspredator(self, compare_file, feature, stat_class, stat_lib):
        gen_file(stat_class, "test")
        gen_file(stat_lib, "test")
        gen_file(os.path.join(os.getcwd(), "test_venn.png"), "test")
        gen_file(os.path.join(os.getcwd(), "test_class.png"), "test")

    def mock_validate_gff(self, compare_file, gff,
                          stat_file, out_cds_file, utr_length):
        gen_file(out_cds_file, self.example.tss_file)

    def mock_stat_ta_tss(self, ta, compare_file, stat_out,
                         ta_tss, tss_ta, fuzzy):
        gen_file("tmp_ta_tss", self.example.tran_file)
        gen_file("tmp_tss", self.example.tss_file)

    def mock_start_to_run(self, tsspredator_path, config_file, out_path, prefix):
        gen_file("test_folder/output/MasterTables/MasterTable_test/MasterTable.tsv", self.example.master)


class Mock_Multiparser(object):

    def parser_wig(self, merge_wigs):
        pass

    def combine_wig(self, gffs, wig_path, test):
        pass

    def parser_gff(self, trans, type_):
        pass

    def combine_gff(self, gffs, tran_path, test, type_):
        pass

    def parser_fasta(self, fastas):
        pass

    def combine_fasta(self, gffs, fasta_path, test):
        pass


class TestsTSSpredator(unittest.TestCase):

    def setUp(self):
        self.mock = Mock_func()
        self.mock_parser = Mock_Multiparser()
        self.example = Example()
        self.test_folder = "test_folder"
        self.trans = "test_folder/trans"
        self.out = "test_folder/output"
        self.wigs = "test_folder/wigs"
        self.gffs = "test_folder/gffs"
        self.tsss = "test_folder/tsss"
        self.fastas = "test_folder/fastas"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
            os.mkdir(self.trans)
            os.mkdir(self.out)
            os.mkdir(self.wigs)
            os.mkdir(self.gffs)
            os.mkdir(self.tsss)
            os.mkdir(self.fastas)
        self.tss = TSSpredator(self.out, self.trans, self.gffs, self.wigs, self.fastas)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)
        if os.path.exists("tmp"):
            shutil.rmtree("tmp")

    def test_print_lib(self):
        out = StringIO()
        lib_list = [{"condition": 1, "replicate": "a", "wig": "test_1.wig"},
                    {"condition": 2, "replicate": "a", "wig": "test_2.wig"}]
        self.tss._print_lib(2, lib_list, out, self.wigs, "test")
        self.assertEqual(out.getvalue(), "test_1a = test_folder/wigs/test_1.wig\ntest_2a = test_folder/wigs/test_2.wig\n")

    def test_import_lib(self):
        out = StringIO()
        libs = ["test1_forward.wig:notex:1:a:+",
                "test1_reverse.wig:notex:1:a:-",
                "test1_TEX_forward.wig:tex:1:a:+",
                "test1_TEX_reverse.wig:tex:1:a:-"]
        gen_file(os.path.join(self.wigs, "test1_forward.wig_STRAIN_test.wig"), "test")
        gen_file(os.path.join(self.wigs, "test1_reverse.wig_STRAIN_test.wig"), "test")
        gen_file(os.path.join(self.wigs, "test1_TEX_forward.wig_STRAIN_test.wig"), "test")
        gen_file(os.path.join(self.wigs, "test1_TEX_reverse.wig_STRAIN_test.wig"), "test")
        self.tss._import_lib(libs, self.wigs, "test", out, "test.gff", "TSS", "test.fa")
        self.assertListEqual(out.getvalue().split("\n"), ["annotation_1 = test.gff",
                                                          "fivePrimeMinus_1a = test_folder/wigs/test1_TEX_reverse.wig",
                                                          "fivePrimePlus_1a = test_folder/wigs/test1_TEX_forward.wig",
                                                          "genome_1 = test.fa", ""])

    def test_gen_config(self):
        os.mkdir(os.path.join(self.out, "MasterTables"))
        os.mkdir(os.path.join(self.wigs, "tmp"))
        config_file = os.path.join(self.test_folder, "config")
        libs = ["test1_forward.wig:notex:1:a:+",
                "test1_reverse.wig:notex:1:a:-",
                "test1_TEX_forward.wig:tex:1:a:+",
                "test1_TEX_reverse.wig:tex:1:a:-"]
        self.tss._gen_config("test", self.out, libs, self.gffs + "/tmp/test.gff",
                             self.wigs + "/tmp", self.fastas + "/tmp/test.fa", 0.3, 2.0, 0.2, 0.5, 2.0, 1.5,
                             0.00, ["test1"], config_file, "TSS",
                             2, 3, 300)
        datas = import_data(config_file)
        self.assertEqual("\n".join(datas), self.example.config)

    def test_set_gen_config(self):
        os.mkdir(os.path.join(self.fastas, "tmp"))
        os.mkdir(os.path.join(self.gffs, "tmp"))
        os.mkdir(os.path.join(self.wigs, "tmp"))
        os.mkdir(os.path.join(self.out, "MasterTables"))
        gen_file(os.path.join(self.fastas, "tmp/test.fa"), "test")
        gen_file(os.path.join(self.gffs, "tmp/test.gff"), "test")
        gen_file(os.path.join(self.wigs, "tmp/test1_forward.wig_STRAIN_test.wig"), "test")
        gen_file(os.path.join(self.wigs, "tmp/test1_reverse.wig_STRAIN_test.wig"), "test")
        gen_file(os.path.join(self.wigs, "tmp/test1_TEX_forward.wig_STRAIN_test.wig"), "test")
        gen_file(os.path.join(self.wigs, "tmp/test1_TEX_reverse.wig_STRAIN_test.wig"), "test")
        libs = ["test1_forward.wig:notex:1:a:+",
                "test1_reverse.wig:notex:1:a:-",
                "test1_TEX_forward.wig:tex:1:a:+",
                "test1_TEX_reverse.wig:tex:1:a:-"]
        self.tss._set_gen_config("TSS", self.test_folder, self.out, libs, 0.3,
                        2.0, 0.2, 0.5, 2.0, 1.5, 300, 0.00, ["test1"], 2, 3)
        datas = import_data(os.path.join(self.test_folder, "config_test.ini"))
        self.assertEqual("\n".join(datas), self.example.config)

    def test_convert_gff(self):
        os.mkdir(os.path.join(self.out, "gffs"))
        os.mkdir(os.path.join(self.out, "MasterTables"))
        os.mkdir(os.path.join(self.out, "MasterTables/MasterTable_test"))
        gen_file(os.path.join(self.out, "MasterTables/MasterTable_test/MasterTable.tsv"), self.example.master)
        self.tss._convert_gff(["test"], self.out, "TSS")
        datas = import_data(os.path.join(self.out, "gffs/test_TSS.gff"))
        self.assertEqual("\n".join(datas), self.example.master_gff)

    def test_merge_wigs(self):
        gen_file(os.path.join(self.wigs, "test1_forward.wig"), "test_f")
        gen_file(os.path.join(self.wigs, "test1_reverse.wig"), "test_r")
        gen_file(os.path.join(self.wigs, "test1_TEX_forward.wig"), "test_f")
        gen_file(os.path.join(self.wigs, "test1_TEX_reverse.wig"), "test_r")
        libs = ["test1_forward.wig:notex:1:a:+",
                "test1_reverse.wig:notex:1:a:-",
                "test1_TEX_forward.wig:tex:1:a:+",
                "test1_TEX_reverse.wig:tex:1:a:-"]
        self.tss._merge_wigs(self.wigs, "test", libs)
        datas = import_data(os.path.join("tmp", "merge_forward.wig"))
        self.assertEqual("\n".join(datas), "test_ftest_f")
        datas = import_data(os.path.join("tmp", "merge_reverse.wig"))
        self.assertEqual("\n".join(datas), "test_rtest_r")
        shutil.rmtree("tmp")

    def test_check_orphan(self):
        os.mkdir(os.path.join(self.out, "gffs"))
        gen_file(os.path.join(self.wigs, "test1_forward.wig"), "test_f")
        gen_file(os.path.join(self.wigs, "test1_reverse.wig"), "test_r")
        gen_file(os.path.join(self.wigs, "test1_TEX_forward.wig"), "test_f")
        gen_file(os.path.join(self.wigs, "test1_TEX_reverse.wig"), "test_r")
        ts.check_orphan = self.mock.mock_check_orphan
        libs = ["test1_TEX_forward.wig:tex:1:a:+", "test1_TEX_reverse.wig:tex:1:a:-",
                "test1_forward.wig:notex:1:a:+", "test1_reverse.wig:notex:1:a:-"]
        self.tss._check_orphan(["test"], self.wigs, "TSS", self.gffs, libs)
        self.assertTrue(os.path.exists(os.path.join(self.out, "gffs/test_TSS.gff")))

    def test_low_expression(self):
        ts.filter_low_expression = self.mock.mock_filter_low_expression
        gen_file(os.path.join(self.wigs, "test1_forward.wig"), "test_f")
        gen_file(os.path.join(self.wigs, "test1_reverse.wig"), "test_r")
        gen_file(os.path.join(self.wigs, "test1_TEX_forward.wig"), "test_f")
        gen_file(os.path.join(self.wigs, "test1_TEX_reverse.wig"), "test_r")
        gen_file(os.path.join(self.gffs, "test_TSS.gff"), self.example.tss_file)
        os.mkdir(os.path.join(self.out, "statistics"))
        os.mkdir(os.path.join(self.out, "statistics/test"))
        libs = ["test1_TEX_forward.wig:tex:1:a:+", "test1_TEX_reverse.wig:tex:1:a:-",
                "test1_forward.wig:notex:1:a:+", "test1_reverse.wig:notex:1:a:-"]
        self.tss._low_expression(100, 3, "manual", libs,
                                 self.gffs, "TSS", self.wigs)
        shutil.rmtree("tmp")
        datas = import_data(os.path.join(self.out, "statistics/test/stat_test_low_expression_cutoff.csv"))
        self.assertEqual("\n".join(datas), "strain\tcutoff_coverage\ntest\t100")

    def test_merge_manual(self):
        gen_file(os.path.join(self.gffs, "test.gff"), self.example.tss_file)
        os.mkdir(os.path.join(self.out, "statistics"))
        os.mkdir(os.path.join(self.out, "statistics/test"))
        os.mkdir(os.path.join(self.out, "gffs"))
        ts.merge_manual_predict_tss = self.mock.mock_merge_manual_predict_tss
        self.tss._merge_manual(["test"], self.gffs, "manual", self.wigs,
                               self.out, 3, "TSS", 300, "libs")
        self.assertTrue(os.path.exists(os.path.join(self.out,
                        "statistics/test/stat_compare_TSSpredator_manual_test.csv")))
        self.assertTrue(os.path.exists(os.path.join(self.out, "gffs/test_TSS.gff")))

    def test_deal_with_overlap(self):
        ts.filter_tss_pro = self.mock.mock_filter_tss_pro
        gen_file(os.path.join(self.out, "test_TSS.gff"), self.example.tss_file)
        gen_file(os.path.join(self.test_folder, "test_processing.gff"), self.example.tss_file)
        self.tss._deal_with_overlap(self.out, "overlap", self.test_folder, "TSS", 3)
        self.assertTrue(os.path.exists(os.path.join(self.out, "test_TSS.gff")))

    def test_stat_tss(self):
        ts.stat_tsspredator = self.mock.mock_stat_tsspredator
        ts.plot_venn = self.mock.mock_plot_venn
        os.mkdir(os.path.join(self.out, "statistics"))
        os.mkdir(os.path.join(self.out, "statistics/test"))
        self.tss._stat_tss(["test"], "TSS")
        self.assertTrue(os.path.exists(os.path.join(self.out, "statistics/test/test_venn.png")))
        self.assertTrue(os.path.exists(os.path.join(self.out, "statistics/test/test_class.png")))

    def test_validate(self):
        gen_file(os.path.join(self.gffs, "test.gff"), self.example.tss_file)
        os.mkdir(os.path.join(self.out, "gffs"))
        ts.validate_gff = self.mock.mock_validate_gff
        self.tss._validate(["test"], self.gffs, 300, self.out)

    def test_compare_ta(self):
        self.tss.multiparser = self.mock_parser
        ts.stat_ta_tss = self.mock.mock_stat_ta_tss
        ta_path = os.path.join(self.trans, "tmp")
        os.mkdir(ta_path)
        os.mkdir(os.path.join(self.out, "gffs"))
        gen_file(os.path.join(ta_path, "test_transcript.gff"), self.example.tran_file)
        self.tss._compare_ta(self.trans, self.gffs, ["test"], 3)
        self.assertTrue(os.path.exists(os.path.join(self.trans, "test_transcript.gff")))
        self.assertTrue(os.path.exists(os.path.join(self.out, "gffs/test_TSS.gff")))

    def test_run_tsspredator(self):
        os.mkdir(os.path.join(self.out, "gffs"))
        os.mkdir(os.path.join(self.out, "statistics"))
        os.mkdir(os.path.join(self.out, "statistics/test"))
        os.mkdir(os.path.join(self.out, "MasterTables"))
        os.mkdir(os.path.join(self.out, "MasterTables/MasterTable_test"))
        os.mkdir(os.path.join(self.out, "configs"))
        ts.stat_tsspredator = self.mock.mock_stat_tsspredator
        ts.plot_venn = self.mock.mock_plot_venn
        ts.validate_gff = self.mock.mock_validate_gff
        ts.stat_ta_tss = self.mock.mock_stat_ta_tss
        ts.filter_tss_pro = self.mock.mock_filter_tss_pro
        ts.merge_manual_predict_tss = self.mock.mock_merge_manual_predict_tss
        ts.filter_low_expression = self.mock.mock_filter_low_expression
        ts.check_orphan = self.mock.mock_check_orphan
        self.tss._start_to_run = self.mock.mock_start_to_run
        gen_file(os.path.join(self.wigs, "test1_forward.wig"), self.example.wig_f)
        gen_file(os.path.join(self.wigs, "test1_reverse.wig"), self.example.wig_r)
        gen_file(os.path.join(self.wigs, "test1_TEX_forward.wig"), self.example.wig_f)
        gen_file(os.path.join(self.wigs, "test1_TEX_reverse.wig"), self.example.wig_r)
        gen_file(os.path.join(self.trans, "test_transcript.gff"), self.example.tran_file)
        gen_file(os.path.join(self.gffs, "test.gff"), self.example.gff_file)
        gen_file(os.path.join(self.fastas, "test.fa"), ">test\nAAATATATATATATAAATTTATATATATATA")
        libs = ["test1_forward.wig:notex:1:a:+", "test1_reverse.wig:notex:1:a:-",
                "test1_TEX_forward.wig:tex:1:a:+", "test1_TEX_reverse.wig:tex:1:a:-"]
        self.tss.run_tsspredator("tsspredator_path", "TSS", self.fastas, 
            self.gffs, self.wigs, libs, "test", 0.3, 0.2,
            2.0, 0.5, 0.00, 2.0, 1.5, 2, self.out, True, True, "manual",
            self.trans, 2, 300, 2, 100, True, "TSS", self.gffs, True)
        self.assertTrue(os.path.exists(os.path.join(self.out, "gffs/test_TSS.gff")))
        self.assertTrue(os.path.exists(os.path.join(self.out,
                        "statistics/test/stat_compare_TSSpredator_manual_test.csv")))


class Example(object):

    tran_file = """test	Transcript	Transcript	3	25	.	+	.	ID=tran0;Name=Transcript_0"""
    gff_file = """test	RefSeq	CDS	5	10	.	+	.	ID=cds0;Name=CDS_0"""
    tss_file = """test	RefSeq	TSS	3	3	.	+	.	ID=tss0;Name=TSS_0"""
    wig_f = """track type=wiggle_0 name="test1_forward.wig"
variableStep chrom=test span=1
2 2.0
3 41.0
4 47.0
5 4.0
6 47.0
7 7.0
8 47.0
9 2.0
10 41.0
11 47.0
12 4.0
13 47.0
14 7.0
15 47.0
16 47.0
18 2.0
19 41.0
20 47.0
21 4.0
22 47.0
23 7.0
24 47.0
25 2.0
26 41.0
27 47.0
28 4.0
29 47.0
30 7.0"""
    wig_r = """track type=wiggle_0 name="test1_reverse.wig"
variableStep chrom=test span=1
2 -2.0
3 -41.0
4 -47.0
5 -4.0
6 -47.0
7 -7.0
8 -47.0
9 -2.0
10 -41.0"""
    config = """TSSinClusterSelectionMethod = HIGHEST
allowedCompareShift = 1
allowedRepCompareShift = 1
annotation_1 = test_folder/gffs/tmp/test.gff
fivePrimeMinus_1a = test_folder/wigs/tmp/test1_TEX_reverse.wig
fivePrimePlus_1a = test_folder/wigs/tmp/test1_TEX_forward.wig
genome_1 = test_folder/fastas/tmp/test.fa
idList = 1
maxASutrLength = 100
maxGapLengthInGene = 500
maxNormalTo5primeFactor = 1.5
maxTSSinClusterDistance = 4
maxUTRlength = 300
min5primeToNormalFactor = 2.0
minCliffFactor = 2.0
minCliffFactorDiscount = 0.5
minCliffHeight = 0.3
minCliffHeightDiscount = 0.2
minNormalHeight = 0.0
minNumRepMatches = 2
minPlateauLength = 0
mode = cond
normPercentile = 0.9
normalMinus_1a = test_folder/wigs/tmp/test1_reverse.wig
normalPlus_1a = test_folder/wigs/tmp/test1_forward.wig
numReplicates = 1
numberOfDatasets = 1
outputDirectory = test_folder/output/MasterTables/MasterTable_test
outputPrefix_1 = test1
projectName = test
superGraphCompatibility = igb
texNormPercentile = 0.5
writeGraphs = 0
writeNocornacFiles = 0"""
    master = """SuperPos	SuperStrand	mapCount	detCount	Genome	detected	enriched	stepHeight	stepFactor	enrichmentFactor	classCount	Pos	Strand	Locus_tag	sRNA/asRNA	Product	UTRlength	GeneLength	Primary	Secondary	Internal	Antisense	Automated	Manual	Putative sRNA	Putative asRNA	Comment	Sequence -50 nt upstream + TSS (51nt)
179	-	1	1	test	1	1	4.45	31.93	8.69	1	179	-	orphan		orphan	NA	NA	0	0	0	0	1	0	0	0		ACCCTTGAATTGAGGGTGTTTTATACCTAAATTTAAAAAATGATGCTATAA
681	-	1	1	test	1	1	4.2	3.0	3.54	2	681	-	HP0001		transcription antitermination protein NusB	48	417	1	0	0	0	1	0	0	0		GATTGAAAGAGCGGGCAGTAAAGCCGGCAATAAGGGCTTTGAAGCGATGAG
681	-	1	1	test	1	1	4.2	3.0	3.54	2	681	-	HP0002		6%2C7-dimethyl-8-ribityllumazine synthase	NA	471	0	0	1	0	1	0	0	0		GATTGAAAGAGCGGGCAGTAAAGCCGGCAATAAGGGCTTTGAAGCGATGAG"""
    master_gff = """##gff-version 3
test	TSSpredator	TSS	179	179	.	-	.	Name=TSS:179_r;ID=tss0;type=Orphan;UTR_length=Orphan_NA;associated_gene=orphan;libs=test
test	TSSpredator	TSS	681	681	.	-	.	Name=TSS:681_r;ID=tss1;type=Primary&Internal;UTR_length=Primary_48&Internal_NA;associated_gene=HP0001&HP0002;libs=test"""

if __name__ == "__main__":
    unittest.main()
