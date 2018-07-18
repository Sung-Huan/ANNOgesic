import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file, import_data, extract_info
import annogesiclib.srna as sr
from annogesiclib.srna import sRNADetection
from mock_args_container import MockClass


current_path = os.getcwd()

class Mock_func(object):

    def __init__(self):
        self.example = Example()

    def mock_run_filter_frag(self, merge_table, merge_gff):
        pass

    def mock_run_format(self, blast_path, database, type_, db_file, err, log):
        pass

    def mock_change_format(self, database, file_):
        gen_file(file_, "test")

    def mock_intergenic_srna(self, gff, tran, tss, pro, fuzzy_tss,
                        max_len, min_len, wig_f, wig_r,
                        merge_wigs, libs, tex_notex, replicates,
                        out_gff, out_table, table_best, decrease_inter,
                        fuzzy_inter, coverage, tolerance):
        pass

    def mock_check_gff(self, gffs):
        pass

    def mock_run_normal(self, import_info, tss_path,
                        pro_path, prefix, gff_path, log):
        return None, [None, None, 1, 2], [None, None, 1, 2]

    def mock_run_utrsrna(
            self, gff_path, gff, tran, fuzzy_tss, max_len, min_len,
            prefix, tex_wigs, frag_wigs, tss, pro, fasta_path,
            tex_notex, flibs, tlibs, replicates, table_best,
            decrease_utr, fuzzy_utr, utr_tex_cover, utr_frag_cover,
            out_folder, hypo, tex_path, frag_path, notex, min_utr):
        pass

    def mock_merge_tex_frag_datas(self, tex_datas, frag_datas):
        return [None, None, 1, 2]

    def mock_merge_srna_gff(self, utr_gff, normal_gff, in_cds, merge_gff, ex_srna):
        shutil.copy("test_folder/gffs/test.gff",
                    "test_folder/output/tmp_merge_test")

    def mock_merge_srna_table(self, merge_gff, normal_table, utr_table,
                              forward, reverse, merge_wig):
        pass

    def mock_get_seq(self, gff, fasta, seq_file):
        gen_file('test_folder/output/sRNA_index_test', ">test\nAAATTTGGGCCC")
        gen_file('test_folder/output/sRNA_2d_test', ">test\n...()()...")

    def mock_extract_energy(self, gff_file, sec_file, energy_file):
        pass

    def mock_run_RNAfold(self, seq_file, vienna_path, sec_file, log):
        pass

    def mock_run_replot(self, vienna_util, tmp_paths, file_,
                        dot_file, rel_file):
        pass

    def mock_convert_pdf(self, ps2pdf14_path, tmp_paths, file_, pdf_file):
        pass

    def mock_run_mountain(self, vienna_util, tmp_paths, dot_file, out, log):
        pass

    def mock_run_blast(self, program, database, e, seq_file, blast_file,
                       strand, para, test, log):
        gen_file('tmp_blast.txt', "test")

    def mock_extract_blast(self, blast_file, srna_file, out_file, csv_file,
                           database_type, score_a, score_b):
        pass

    def mock_classify_srna(self, gff_file, out_stat, in_cds, import_info):
        pass

    def mock_gen_srna_table(self, class_gff, merge, nr, srna,
                            max_len, min_len, out_table):
        pass
    
    def mock_blast_class(self, srna, out_srna_blast):
        pass

    def mock_srna_sorf_comparison(self, gff, sorf, tmp_srna, tmp_sorf):
        pass

    def mock_check_database(self, database, db_format):
        pass

    def mock_merge_blast_out(self, file1, file2):
        pass

class Mock_multiparser(object):

    def parser_gff(tsss, type_):
        pass

    def combine_gff(gffs, tss_path, detect, type_):
        pass

    def parser_fasta(fasta):
        pass

    def combine_fasta(gffs, fasta_path, detect):
        pass    


class TestsRNADetection(unittest.TestCase):

    def setUp(self):
        self.mock_args = MockClass()
        self.example = Example()
        self.mock = Mock_func()
        self.test_folder = "test_folder"
        self.gffs = "test_folder/gffs"
        self.tsss = "test_folder/tsss"
        self.sorf = "test_folder/sORF"
        self.out = "test_folder/output"
        self.trans = "test_folder/trans"
        self.fastas = "test_folder/fastas"
        self.tex = "test_folder/tex"
        self.frag = "test_folder/frag"
        self.pros = "test_folder/pros"
        self.terms = "test_folder/terms"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
            os.mkdir(self.gffs)
            os.mkdir(self.tsss)
            os.mkdir(os.path.join(self.tsss, "tmp"))
            os.mkdir(self.out)
            os.mkdir(self.trans)
            os.mkdir(os.path.join(self.trans, "tmp"))
            os.mkdir(self.fastas)
            os.mkdir(os.path.join(self.fastas, "tmp"))
            os.mkdir(self.tex)
            os.mkdir(self.frag)
            os.mkdir(self.pros)
            os.mkdir(os.path.join(self.pros, "tmp"))
            os.mkdir(self.sorf)
            os.mkdir(os.path.join(self.sorf, "tmp"))
            os.mkdir(self.terms)
        args = self.mock_args.mock()
        args.tss_folder = self.tsss
        args.pro_folder = self.pros
        args.out_folder = self.out
        args.sorf_file = self.sorf
        args.fastas = self.fastas
        args.trans = self.trans
        args.terms = self.terms
        self.srna = sRNADetection(args)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)
        os.chdir(current_path)
        if os.path.exists("tmp"):
            shutil.rmtree("tmp")
        if os.path.exists("tmp_srna.csv"):
            os.remove("tmp_srna.csv")
        if os.path.exists("tmp_srna.gff"):
            os.remove("tmp_srna.gff")
        if os.path.exists("tmp_blast.txt"):
            os.remove("tmp_blast.txt")

    def test_check_folder_exist(self):
        path_ = self.srna._check_folder_exist(self.sorf)
        self.assertEqual(path_, "test_folder/sORF/tmp")

    def test_formatdb(self):
        database = "test_folder/test.fa"
        gen_file(database, "test")
        log = open(os.path.join(self.test_folder, "test.log"), "w")
        sr.change_format = self.mock.mock_change_format
        self.srna._run_format = self.mock.mock_run_format
        self.srna._formatdb(database, "type_", self.out, "blast_path", "sRNA", log)
        self.assertTrue(os.path.exists(os.path.join(self.out, "log.txt")))

    def test_check_necessary_file(self):
        self.srna.multiparser = Mock_multiparser
        self.srna._check_gff = self.mock.mock_check_gff
        self.srna._check_database = self.mock.mock_check_database
        args = self.mock_args.mock()
        args.trans = self.trans
        args.tsss = self.tsss
        args.pros = self.pros
        args.import_info = ["tss", "blast_nr", "blast_srna", "sec_str", "sorf"]
        args.fastas = self.fastas
        args.terms = self.terms
        args.sorf_file = self.sorf
        args.gffs = self.gffs
        args.tex_wigs = self.tex
        args.frag_wigs = self.frag
        args.utr_srna = True
        args.nr_format = True
        args.srna_format = True
        args.nr_database = "test"
        args.srna_database = "test"
        log = open(os.path.join(self.test_folder, "test.log"), "w")
        self.srna._check_necessary_file(args, log)

    def test_run_program(self):
        self.srna.multiparser = Mock_multiparser
        self.srna._check_gff = self.mock.mock_check_gff
        self.srna._run_normal = self.mock.mock_run_normal
        self.srna._run_utrsrna = self.mock.mock_run_utrsrna
        self.srna._merge_tex_frag_datas = self.mock.mock_merge_tex_frag_datas
        sr.filter_frag = self.mock.mock_run_filter_frag
        sr.merge_srna_gff = self.mock.mock_merge_srna_gff
        sr.merge_srna_table = self.mock.mock_merge_srna_table
        gen_file(os.path.join(self.gffs, "test.gff"), self.example.sorf_file)
        gen_file(os.path.join(self.trans, "test_transcript.gff"),
                 self.example.sorf_file)
        gen_file(os.path.join(self.tsss, "test_TSS.gff"),
                 self.example.sorf_file)
        gen_file(os.path.join(self.tsss, "test_processing.gff"),
                 self.example.sorf_file)
        fuzzy_tsss = {"inter": 3}
        args = self.mock_args.mock()
        args.import_info = ["tss", "blast_nr", "blast_srna", "sec_str", "sorf"]
        args.trans = self.trans
        args.tsss = self.tsss
        args.pros = self.pros
        args.max_len = 300
        args.min_len = 30
        args.tex_notex = "tex_notex"
        args.fuzzy_tsss = fuzzy_tsss
        args.out_folder = self.out
        args.table_best = True
        args.wig_path = "wig_path"
        args.merge_wigs = "merge"
        args.libs = "libs"
        args.gffs = self.gffs
        args.in_cds = False
        args.utr_srna = True
        args.ex_srna = False
        args.cutoff_overlap = 0.5
        args.source = True
        log = open(os.path.join(self.test_folder, "test.log"), "w")
        prefixs = self.srna._run_program(args, log)
        self.assertListEqual(prefixs, ['test'])

    def test_get_seq_sec(self):
        sr.extract_energy = self.mock.mock_extract_energy
        self.srna.helper.get_seq = self.mock.mock_get_seq
        self.srna._run_RNAfold = self.mock.mock_run_RNAfold
        os.mkdir(os.path.join(self.out, "tmp_srna"))
        gen_file(os.path.join(self.fastas, "test.fa"), ">test\nAAATTTGGGCCC")
        log = open(os.path.join(self.test_folder, "test.log"), "w")
        datas = self.srna._get_seq_sec(
            self.fastas, self.out, "test", self.test_folder,
            self.test_folder, "vienna_path", log)
        self.assertEqual(datas["sec"].split("/")[-1], "test_folder")
        self.assertEqual(datas["dot"].split("/")[-1], "test_folder")
        self.assertEqual(datas["main"].split("/")[-1],
                         datas["tmp"].split("/")[-4])
        self.assertEqual(datas["tmp"].split("/")[-1], "tmp_srna")

    def test_replot_sec(self):
        self.srna._run_replot = self.mock.mock_run_replot
        self.srna._convert_pdf = self.mock.mock_convert_pdf
        gen_file(os.path.join(self.tsss, "test.rss.ps"), "test")
        gen_file(os.path.join(self.tsss, "test.dp.ps"), "test")
        log = open(os.path.join(self.test_folder, "test.log"), "w")
        tmp_paths = {"dot": self.out, "sec": self.fastas, "tmp": self.tsss}
        self.srna._replot_sec("vienna_util", tmp_paths, "test", log)
        self.assertTrue(os.path.exists(os.path.join(
            tmp_paths["dot"], "test/test.dp.ps")))
        self.assertTrue(os.path.exists(os.path.join(
            tmp_paths["sec"], "test/test.rss.ps")))

    def test_plot_mountain(self):
        self.srna._run_mountain = self.mock.mock_run_mountain
        tmp_paths = {"main": self.test_folder, "tmp": self.tsss,
                     "dot": self.sorf}
        moun_path = "fastas"
        log = open(os.path.join(self.test_folder, "test.log"), "w")
        os.mkdir(os.path.join(tmp_paths["dot"], "test"))
        gen_file(os.path.join(tmp_paths["dot"], "test/test.dp.ps"), "test")
        self.srna._plot_mountain(True, moun_path,
                                 tmp_paths, "test", "vienna_util", log)
        self.assertTrue("test_folder/fastas/test/test.mountain.pdf")

    def test_compute_2d_and_energy(self):
        sr.extract_energy = self.mock.mock_extract_energy
        sr.change_format = self.mock.mock_change_format
        self.srna._run_replot = self.mock.mock_run_replot
        self.srna._convert_pdf = self.mock.mock_convert_pdf
        self.srna._run_mountain = self.mock.mock_run_mountain
        sec_path = os.path.join(self.out, "figs")
        os.mkdir(sec_path)
        os.mkdir(os.path.join(sec_path, "sec_plots"))
        os.mkdir(os.path.join(sec_path, "dot_plots"))
        os.mkdir(os.path.join(sec_path, "mountain_plots"))
        tmp_paths = {"dot": self.out, "sec": self.fastas,
                     "tmp": self.tsss, "main": self.test_folder}
        gen_file(os.path.join(self.fastas, "tmp/test.fa"),
                 ">test\nAAATTTGGGCCC")
        gen_file(os.path.join(self.out, "tmp_basic_test"),
                 self.example.srna_file)
        gen_file(os.path.join(self.out, "tmp_energy_test"), "test")
        log = open(os.path.join(self.test_folder, "test.log"), "w")
        args = self.mock_args.mock()
        args.out_folder = self.out
        args.fastas = self.fastas
        args.rnafold = "test"
        args.relplot_pl = "test"
        args.mountain_pl = "test"
        args.mountain = True
        args.ps2pdf14_path = "test"
        self.srna._compute_2d_and_energy(args, ["test"], log)
        datas = import_data(os.path.join(self.out, "tmp_basic_test"))
        self.assertEqual("\n".join(datas), "test")

    def test_blast(self):
        self.srna.helper.merge_blast_out = self.mock.mock_merge_blast_out
        sr.extract_blast = self.mock.mock_extract_blast
        self.srna._run_blast = self.mock.mock_run_blast
        self.srna._run_format = self.mock.mock_run_format
        gen_file(os.path.join(self.out, "tmp_basic_test"),
                 self.example.srna_file)
        gen_file(os.path.join(self.out, "tmp_nr_test"), "test")
        gen_file(os.path.join(self.fastas, "tmp/test.fa"),
                 ">test\nAAATTTGGGCCC")
        args = self.mock_args.mock()
        args.blast_path = "test"
        args.para_blast = 1
        log = open(os.path.join(self.test_folder, "test.log"), "w")
        args.fastas = self.fastas
        args.out_folder = self.out
        args.blast_score_s = 0
        args.blast_score_n = 0
        self.srna._blast("database", False, "dna", args,
                         ["test"], "blast_all", "nr", 0.0001, "tss", log)
        datas = import_data(os.path.join(self.out, "tmp_basic_test"))
        self.assertEqual("\n".join(datas), "test")

    def test_class_srna(self):
        sr.classify_srna = self.mock.mock_classify_srna
        sr.gen_srna_table = self.mock.mock_gen_srna_table
        gff_out = os.path.join(self.out, "gffs")
        table_out = os.path.join(self.out, "tables")
        stat_out = os.path.join(self.out, "stat")
        os.mkdir(gff_out)
        os.mkdir(table_out)
        os.mkdir(stat_out)
        log = open(os.path.join(self.test_folder, "test.log"), "w")
        os.mkdir(os.path.join(table_out, "for_classes"))
        os.mkdir(os.path.join(gff_out, "for_classes"))
        args = self.mock_args.mock()
        args.max_len = 300
        args.min_len = 30
        args.import_info = ["tss", "blast_nr", "blast_srna", "sec_str", "sorf"]
        self.srna._class_srna(["test"], args, log)
        self.assertTrue(os.path.exists(os.path.join(
            gff_out, "for_classes/test")))
        self.assertTrue(os.path.exists(os.path.join(
            table_out, "for_classes/test")))

    def test_filter_srna(self):
        self.srna.helper.merge_blast_out = self.mock.mock_merge_blast_out
        sr.classify_srna = self.mock.mock_classify_srna
        sr.gen_srna_table = self.mock.mock_gen_srna_table
        sr.extract_blast = self.mock.mock_extract_blast
        self.srna._run_blast = self.mock.mock_run_blast
        self.srna._run_format = self.mock.mock_run_format
        sr.extract_energy = self.mock.mock_extract_energy
        sr.change_format = self.mock.mock_change_format
        self.srna._run_replot = self.mock.mock_run_replot
        self.srna._convert_pdf = self.mock.mock_convert_pdf
        self.srna._run_mountain = self.mock.mock_run_mountain
        self.srna.multiparser = Mock_multiparser
        self.srna._check_gff = self.mock.mock_check_gff
        self.srna._run_normal = self.mock.mock_run_normal
        self.srna._run_utrsrna = self.mock.mock_run_utrsrna
        sr.merge_srna_gff = self.mock.mock_merge_srna_gff
        sr.merge_srna_table = self.mock.mock_merge_srna_table
        sr.extract_energy = self.mock.mock_extract_energy
        self.srna.helper.get_seq = self.mock.mock_get_seq
        self.srna._run_RNAfold = self.mock.mock_run_RNAfold
        stat_out = os.path.join(self.out, "stat")
        if "mountain_plot" not in os.listdir(self.out):
            os.mkdir(os.path.join(self.out, "mountain_plot"))
        sec_path = os.path.join(self.out, "sec_structure")
        if "sec_structure" not in os.listdir(self.out):
            os.mkdir(sec_path)
            os.mkdir(os.path.join(sec_path, "sec_plot"))
            os.mkdir(os.path.join(sec_path, "dot_plot"))
        gen_file(os.path.join(self.fastas, "tmp/test.fa"),
                 ">test\nAAATTTGGGCCC")
        gen_file(os.path.join(self.out, "sRNA_seq_test"),
                 ">test\nAAATTTGGGCCC")
        gen_file(os.path.join(self.out, "sRNA_index_test"),
                 ">test\nAAATTTGGGCCC")
        gen_file(os.path.join(self.out, "tmp_basic_test"),
                 self.example.srna_file)
        gen_file(os.path.join(self.out, "tmp_energy_test"), "test")
        gen_file(os.path.join(self.out, "tmp_nr_test"), "test")
        gen_file(os.path.join(self.out, "tmp_sRNA_test"), "test")
        gen_file(os.path.join(self.out, "tmp_sRNA_test.csv"), "test")
        gen_file(os.path.join(self.test_folder, "srna"), "test")
        gen_file(os.path.join(self.test_folder, "nr"), "test")
        sr.blast_class = self.mock.mock_blast_class
        sr.srna_sorf_comparison = self.mock.mock_srna_sorf_comparison
        args = self.mock_args.mock()
        args.import_info = ["tss", "blast_nr", "blast_srna", "sec_str", "sorf"]
        args.out_folder = self.out
        args.fastas = self.fastas
        args.rnafold = "test"
        args.relplot_pl = "test"
        args.mountain_pl = "test"
        args.table_best = True
        args.in_cds = False
        args.ps2pdf14_path = "test"
        args.sorf_file = self.sorf
        args.mountain = True
        args.nr_database = os.path.join(self.test_folder, "nr")
        args.srna_database = os.path.join(self.test_folder, "srna")
        args.blastx = "blast_path"
        args.blastn = "blast_path"
        args.nr_format = False
        log = open(os.path.join(self.test_folder, "test.log"), "w")
        args.srna_format = False
        args.compute_sec_str = False
        args.e_nr = 0
        args.e_srna = 0
        args.para_blast = 1
        args.blast_score_s = 0
        args.blast_score_n = 0
        self.srna._filter_srna(args, ["test"], log)
        datas = import_data(os.path.join(self.out, "tmp_basic_test"))
        self.assertEqual("\n".join(datas), "test")


class Example(object):
    sorf_file = """aaa	RefSeq	sORF	140	160	.	+	.	ID=sorf0;Name=sORF_0
aaa	RefSeq	sORF	230	280	.	+	.	ID=sorf1;Name=sORF_1"""
    srna_file = """test	RefSeq	sRNA	1	5	.	+	.	ID=srna0;Name=sRNA_0"""

if __name__ == "__main__":
    unittest.main()
