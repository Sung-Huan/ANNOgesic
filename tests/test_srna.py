import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file, import_data, extract_info
import annogesiclib.srna as sr
from annogesiclib.srna import sRNADetection

current_path = os.getcwd()

class Mock_func(object):

    def __init__(self):
        self.example = Example()

    def mock_run_format(self, blast_path, database, type_, db_file, err):
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

    def mock_run_normal(self, import_info, tss_path, pro_path, prefix,
                        gff_path, gff, tran, fuzzy_tsss,
                        max_len, min_len, wig_path, coverage,
                        merge_wigs, libs, tex_notex, replicates,
                        table_best, decrease_inter, fuzzy_inter,
                        out_folder, tolerance):
        pass

    def mock_run_utrsrna(self, gff_path, gff, tran, fuzzy_tsss,
                          merge_wigs, max_len, min_len, wig_path, prefix,
                          tss, pro, fasta_path, libs, tex_notex,
                          replicates, table_best, decrease_utr,
                          fuzzy_utr, utr5_coverage, utr3_coverage,
                          intercds_coverage, out_folder):
        pass

    def mock_merge_srna_gff(self, utr_gff, normal_gff, merge_gff):
        shutil.copy("test_folder/test.gff", merge_gff)

    def mock_merge_srna_table(self, merge_gff, normal_table, utr_table,
                              wig_f, wig_r, merge_wigs, libs, tex_notex,
                              replicates, table_best, merge_table):
        pass

    def mock_get_seq(self, gff, fasta, seq_file):
        pass

    def mock_extract_energy(self, gff_file, sec_file, energy_file):
        pass

    def mock_run_RNAfold(self, seq_file, vienna_path, sec_file):
        pass

    def mock_run_replot(self, vienna_util, tmp_paths, file_, dot_file, rel_file):
        pass

    def mock_convert_pdf(self, ps2pdf14_path, tmp_paths, file_, pdf_file):
        pass

    def mock_run_mountain(self, vienna_util, tmp_paths, dot_file, out):
        pass

    def mock_run_blast(self, blast_path, program, database, e, seq_file, blast_file):
        pass

    def mock_extract_blast(self, blast_file, srna_file, out_file, csv_file, database_type):
        pass

    def mock_classify_srna(self, gff_file, class_gff, energy, nr_hit_num, out_stat):
        pass

    def mock_gen_srna_table(self, class_gff, merge, nr, srna, max_len, min_len, out_table):
        pass
    
    def mock_blast_class(self, srna, out_srna_blast):
        pass

    def mock_srna_sorf_comparison(self, gff, sorf, tmp_srna, tmp_sorf):
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
        self.example = Example()
        self.mock = Mock_func()
        self.test_folder = "test_folder"
        self.tsss = "test_folder/tsss"
        self.sorf = "test_folder/sORF"
        self.out = "test_folder/output"
        self.trans = "test_folder/trans"
        self.fastas = "test_folder/fastas"
        self.tex = "test_folder/tex"
        self.frag = "test_folder/frag"
        self.pros = "test_folder/pros"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
            os.mkdir(self.tsss)
            os.mkdir(self.out)
            os.mkdir(self.trans)
            os.mkdir(self.fastas)
            os.mkdir(self.tex)
            os.mkdir(self.frag)
            os.mkdir(self.pros)
            os.mkdir(self.sorf)
        self.srna = sRNADetection(self.out, self.tsss, self.pros,
                                  self.sorf, self.fastas, self.trans)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)
        os.chdir(current_path)
        if os.path.exists("tmp"):
            shutil.rmtree("tmp")

    def test_check_folder_exist(self):
        path_ = self.srna._check_folder_exist(self.sorf)
        self.assertEqual(path_, "test_folder/sORF")

    def test_formatdb(self):
        database = "test_folder/test.fa"
        gen_file(database, "test")
        sr.change_format = self.mock.mock_change_format
        self.srna._run_format = self.mock.mock_run_format
        self.srna._formatdb(database, "type_", self.out, "blast_path", "sRNA")
        self.assertTrue(os.path.exists(os.path.join(self.out, "log.txt")))

    def test_check_necessary_file(self):
        self.srna.multiparser = Mock_multiparser
        self.srna._check_gff = self.mock.mock_check_gff
        self.srna._check_necessary_file(self.test_folder, self.trans, self.tex,
                                        self.frag, True, self.tsss, self.pros,
                                        self.sorf, ["1", "2", "3", "4", "5"],
                                        self.fastas)

    def test_merge_wig(self):
        gen_file(os.path.join(self.tex, "tex.wig"), "test")
        gen_file(os.path.join(self.frag, "frag.wig"), "test")
        merge = self.srna._merge_wig(self.tex, self.frag, self.out)
        self.assertEqual(merge, "test_folder/output/merge_wigs")

    def test_merge_libs(self):
        libs = self.srna._merge_libs("tlibs", "flibs")
        self.assertEqual(libs, "tlibsflibs")

    def test_run_program(self):
        self.srna.multiparser = Mock_multiparser
        self.srna._check_gff = self.mock.mock_check_gff
        self.srna._run_normal = self.mock.mock_run_normal
        self.srna._run_utrsrna = self.mock.mock_run_utrsrna
        sr.merge_srna_gff = self.mock.mock_merge_srna_gff
        sr.merge_srna_table = self.mock.mock_merge_srna_table
        gen_file(os.path.join(self.test_folder, "test.gff"), self.example.sorf_file)
        gen_file(os.path.join(self.trans, "test_transcript.gff"), self.example.sorf_file)
        gen_file(os.path.join(self.tsss, "test_TSS.gff"), self.example.sorf_file)
        gen_file(os.path.join(self.tsss, "test_processing.gff"), self.example.sorf_file)
        fuzzy_tsss = {"inter": 3}
        prefixs = self.srna._run_program(self.test_folder, ["1", "2", "3", "4", "5"],
                                         self.trans, self.tsss, self.pros,
                                         "wig_path", 300, 30, "coverage", "merge_wigs",
                                         "libs", "tex_notex", "replicates", True,
                                         50, 50, 5, 5, fuzzy_tsss, 30, 30, 30, self.out, True,
                                         "fasta_path", 5)
        self.assertListEqual(prefixs, ['test'])

    def test_get_seq_sec(self):
        sr.extract_energy = self.mock.mock_extract_energy
        self.srna.helper.get_seq = self.mock.mock_get_seq
        self.srna._run_RNAfold = self.mock.mock_run_RNAfold
        os.mkdir(os.path.join(self.out, "tmp_srna"))
        gen_file(os.path.join(self.fastas, "test.fa"), ">test\nAAATTTGGGCCC")
        datas = self.srna._get_seq_sec(self.fastas, self.out, "test", self.test_folder,
                                       self.test_folder, "vienna_path")
        self.assertEqual(datas["sec"].split("/")[-1], "test_folder")
        self.assertEqual(datas["dot"].split("/")[-1], "test_folder")
        self.assertEqual(datas["main"].split("/")[-1], datas["tmp"].split("/")[-4])
        self.assertEqual(datas["tmp"].split("/")[-1], "tmp_srna")

    def test_replot_sec_to_pdf(self):
        self.srna._run_replot = self.mock.mock_run_replot
        self.srna._convert_pdf = self.mock.mock_convert_pdf
        gen_file(os.path.join(self.tsss, "test.rss.pdf"), "test")
        gen_file(os.path.join(self.tsss, "test.dp.pdf"), "test")
        tmp_paths = {"dot": self.out, "sec": self.fastas, "tmp": self.tsss}
        self.srna._replot_sec_to_pdf("vienna_util", tmp_paths, "ps2pdf14_path", "test")
        self.assertTrue(os.path.exists(os.path.join(tmp_paths["dot"], "test/test.dp.pdf")))
        self.assertTrue(os.path.exists(os.path.join(tmp_paths["sec"], "test/test.rss.pdf")))

    def test_plot_mountain(self):
        self.srna._run_mountain = self.mock.mock_run_mountain
        tmp_paths = {"main": self.test_folder, "tmp": self.tsss}
        moun_path = "fastas"
        gen_file(os.path.join(tmp_paths["tmp"], "test.dp.ps"), "test")
        self.srna._plot_mountain(True, moun_path,
                                 tmp_paths, "test", "vienna_util")
        self.assertTrue("test_folder/fastas/test/test.mountain.pdf")

    def test_compute_2d_and_energy(self):
        sr.extract_energy = self.mock.mock_extract_energy
        sr.change_format = self.mock.mock_change_format
        self.srna._run_replot = self.mock.mock_run_replot
        self.srna._convert_pdf = self.mock.mock_convert_pdf
        self.srna._run_mountain = self.mock.mock_run_mountain
        os.mkdir(os.path.join(self.out, "mountain_plot"))
        sec_path = os.path.join(self.out, "sec_structure")
        os.mkdir(sec_path)
        os.mkdir(os.path.join(sec_path, "sec_plot"))
        os.mkdir(os.path.join(sec_path, "dot_plot"))
        tmp_paths = {"dot": self.out, "sec": self.fastas, "tmp": self.tsss, "main": self.test_folder}
        gen_file(os.path.join(self.fastas, "test.fa"), ">test\nAAATTTGGGCCC")
        gen_file(os.path.join(self.out, "tmp_basic_test"), self.example.srna_file)
        gen_file(os.path.join(self.out, "tmp_energy_test"), "test")
        self.srna._compute_2d_and_energy(self.out, ["test"], self.fastas, "vienna_path",
                                         "vienna_util", True, "ps2pdf14_path")
        datas = import_data(os.path.join(self.out, "tmp_basic_test"))
        self.assertEqual("\n".join(datas), "test")

    def test_blast(self):
        sr.extract_blast = self.mock.mock_extract_blast
        self.srna._run_blast = self.mock.mock_run_blast
        self.srna._run_format = self.mock.mock_run_format
        gen_file(os.path.join(self.out, "tmp_basic_test"), self.example.srna_file)
        gen_file(os.path.join(self.out, "tmp_nr_test"), "test")
        gen_file(os.path.join(self.fastas, "test.fa"), ">test\nAAATTTGGGCCC")
        self.srna._blast("database", False, "dna", self.out,
                         "blast_path", ["test"], self.fastas, "blast_all", "nr", 0.0001)
        datas = import_data(os.path.join(self.out, "tmp_basic_test"))
        self.assertEqual("\n".join(datas), "test")

    def test_class_srna(self):
        sr.classify_srna = self.mock.mock_classify_srna
        sr.gen_srna_table = self.mock.mock_gen_srna_table
        gff_out = os.path.join(self.out, "gff")
        table_out = os.path.join(self.out, "table")
        stat_out = os.path.join(self.out, "stat")
        os.mkdir(gff_out)
        os.mkdir(table_out)
        os.mkdir(stat_out)
        os.mkdir(os.path.join(table_out, "for_class"))
        os.mkdir(os.path.join(gff_out, "for_class"))
        self.srna._class_srna(["1", "2", "3", "4", "5"], ["test"], gff_out, table_out,
                              stat_out, 0, 0, 300, 30)
        self.assertTrue(os.path.exists(os.path.join(gff_out, "for_class/test")))
        self.assertTrue(os.path.exists(os.path.join(table_out, "for_class/test")))

    def test_filter_srna(self):
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
        os.mkdir(os.path.join(self.out, "mountain_plot"))
        sec_path = os.path.join(self.out, "sec_structure")
        os.mkdir(sec_path)
        os.mkdir(os.path.join(sec_path, "sec_plot"))
        os.mkdir(os.path.join(sec_path, "dot_plot"))
        gen_file(os.path.join(self.fastas, "test.fa"), ">test\nAAATTTGGGCCC")
        gen_file(os.path.join(self.out, "tmp_basic_test"), self.example.srna_file)
        gen_file(os.path.join(self.out, "tmp_energy_test"), "test")
        gen_file(os.path.join(self.out, "tmp_nr_test"), "test")
        gen_file(os.path.join(self.out, "tmp_sRNA_test"), "test")
        gen_file(os.path.join(self.out, "tmp_sRNA_test.csv"), "test")
        gen_file(os.path.join(self.test_folder, "srna"), "test")
        gen_file(os.path.join(self.test_folder, "nr"), "test")
        sr.blast_class = self.mock.mock_blast_class
        sr.srna_sorf_comparison = self.mock.mock_srna_sorf_comparison
        self.srna._filter_srna(["1", "2", "3", "4", "5"], self.out, ["test"],
                     self.fastas, "vienna_path", "vienna_util", True, "ps2pdf14_path",
                     os.path.join(self.test_folder, "nr"), True, True, "blast_path",
                     os.path.join(self.test_folder, "srna"),
                     stat_out, self.sorf, 0, 0)
        datas = import_data(os.path.join(self.out, "tmp_basic_test"))
        self.assertEqual("\n".join(datas), "test\tRefSeq\tsRNA\t1\t5\t.\t+\t.\tID=srna0;Name=sRNA_0")

    def test_get_replicates(self):
        reps = self.srna._get_replicates(5, 10)
        self.assertDictEqual(reps, {'frag': 10, 'tex': 5})


class Example(object):
    sorf_file = """aaa	RefSeq	sORF	140	160	.	+	.	ID=sorf0;Name=sORF_0
aaa	RefSeq	sORF	230	280	.	+	.	ID=sorf1;Name=sORF_1"""
    srna_file = """test	RefSeq	sRNA	1	5	.	+	.	ID=srna0;Name=sRNA_0"""

if __name__ == "__main__":
    unittest.main()
