import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file, import_data, extract_info
import annogesiclib.terminator as te
from annogesiclib.terminator import Terminator

class Mock_func(object):

    def mock_check_gff_file(self, folder):
        pass

    def mock_remove_tmp_file(self, gffs, fastas, sRNAs, tex_wigs, frag_wigs,
                             trans, term_outfolder, out_folder, merge_wigs):
        pass

    def mock_TransTermHP(self, TransTermHP_path, expterm_path, fasta,
                         combine_path, file_, out_path, prefix, out):
        pass

    def mock_intergenic_seq(self, fasta, tran_file, gff_file, tmp_seq):
        pass

    def mock_poly_t(self, tmp_seq, tmp_sec, gff_file, tran_file, fuzzy_up_ta,
                    fuzzy_down_ta, fuzzy_up_cds, fuzzy_down_cds, tmp_cand):
        pass

    def mock_detect_coverage(self, tmp_cand, gff, tran, fasta, wig_f, wig_r,
                             fuzzy, cutoff_coverage, hp, merge_wigs, libs,
                             tex_notex, replicates, term_gff, term_table,
                             table_best, decrease):
        pass

    def mock_run_rnafold(self, RNAfold_path, tmp_seq, tmp_sec, prefix):
        gen_file(tmp_seq, "test")
        gen_file(tmp_sec, "test")

    def mock_stat_term(self, gff, table, stat, detect, express):
        pass


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


class TestTerminator(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.mock = Mock_func()
        self.mock_parser = Mock_Multiparser()
        self.test_folder = "test_folder"
        self.out = "test_folder/output"
        self.fastas = "test_folder/fastas"
        self.gffs = "test_folder/gffs"
        self.srnas = "test_folder/srnas"
        self.trans = "test_folder/trans"
        if (not os.path.exists(self.test_folder)):
            print("AAAA")
            os.mkdir(self.test_folder)
            os.mkdir(self.out)
            os.mkdir(self.fastas)
            os.mkdir(self.gffs)
            os.mkdir(self.srnas)
            os.mkdir(self.trans)
            os.mkdir(os.path.join(self.out, "tables"))
            os.mkdir(os.path.join(self.out, "gffs"))
        self.term = Terminator(self.gffs, self.fastas,
                               self.trans, self.out, self.srnas)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)
        if os.path.exists("tmp_transterm"):
            shutil.rmtree("tmp_transterm")

#    def test_convert_gff2rntptt(self):
#        os.mkdir(os.path.join(self.srnas, "tmp"))
#        gen_file(os.path.join(self.gffs, "aaa.gff"), self.example.gff_file)
#        gen_file(os.path.join(self.srnas, "aaa_sRNA.gff"), self.example.srna_file)
#        gen_file(os.path.join(self.fastas, "aaa.fa"), self.example.fasta_file)
#        file_types, prefixs = self.term._convert_gff2rntptt(self.gffs, self.fastas, self.srnas)
#        self.assertDictEqual(file_types, {'aaa': 'srna'})
#        self.assertListEqual(prefixs, ['aaa'])
#
#    def test_combine_annotation(self):
#        test1 = os.path.join(self.test_folder, "test1.ptt")
#        test2 = os.path.join(self.test_folder, "test2.ptt")
#        gen_file(test1, self.example.ptt)
#        gen_file(test2, self.example.ptt)
#        files = [test1, test2]
#        combine_file = os.path.join(self.test_folder, "combine")
#        self.term._combine_annotation(combine_file, files)
#        datas = import_data(combine_file)
#        result = self.example.ptt.split("\n")[3:]
#        self.assertEqual("\n".join(datas), "\n".join(result + result))
#
#    def test_run_TransTermHP(self):
#        self.term._TransTermHP = self.mock.mock_TransTermHP
#        gen_file(os.path.join(self.gffs, "aaa.ptt"), self.example.ptt)
#        gen_file(os.path.join(self.fastas, "aaa.fa"), self.example.fasta_file)
#        self.term._run_TransTermHP("TransTermHP_path", self.gffs,
#                                   self.fastas, self.out, "expterm_path")
#        self.assertTrue(os.path.exists(os.path.join(self.out, "aaa")))
#
    def test_convert_to_gff(self):
        self.term.multiparser = self.mock_parser
        hp_folder = os.path.join(self.out, "aaa")
        os.mkdir(hp_folder)
        gen_file(os.path.join(hp_folder, "aaa_best_terminator_after_gene.bag"), self.example.bag)
        os.mkdir("tmp_transterm")
        self.term._convert_to_gff(["aaa"], self.out, self.gffs)
        datas = import_data("/home/silas/ANNOgesic/tmp_transterm/aaa_transtermhp.gff")
        self.assertEqual("\n".join(datas), self.example.gff_bag)

#    def test_combine_libs_wigs(self):
#        tex_wigs = os.path.join(self.test_folder, "tex")
#        frag_wigs = os.path.join(self.test_folder, "frag")
#        os.mkdir(tex_wigs)
#        os.mkdir(frag_wigs)
#        gen_file(os.path.join(frag_wigs, "frag.wig"), "text")
#        gen_file(os.path.join(tex_wigs, "tex.wig"), "text")
#        merge, libs = self.term._combine_libs_wigs("tlibs", "flibs", tex_wigs, frag_wigs)
#        self.assertEqual(merge, "test_folder/merge_wigs")
#        self.assertEqual(libs, "tlibsflibs")
#
#    def test_merge_sRNA(self):
#        os.mkdir(os.path.join(self.srnas, "tmp"))
#        self.term.multiparser = self.mock_parser
#        gen_file(os.path.join(self.gffs, "aaa.gff"), self.example.gff_file)
#        gen_file(os.path.join(self.srnas, "tmp/aaa_sRNA.gff"), self.example.srna_file)
#        merge = self.term._merge_sRNA(self.srnas, ["aaa"], self.gffs)
#        self.assertEqual(merge.split("/")[-1], "tmp_merge_gff")
#        shutil.rmtree("tmp_merge_gff")
#
#    def test_move_file(self):
#        term_outfolder = self.gffs
#        csv_outfolder = self.out
#        gen_file(os.path.join(term_outfolder, "aaa_term.gff"), self.example.term_file)
#        os.mkdir("tmp_term_table")
#        gen_file("tmp_term_table/aaa_term_raw.csv", "test")
#        self.term._move_file(term_outfolder, csv_outfolder)
#        shutil.rmtree("tmp_term_table")
#        self.assertTrue("test_folder/output/gffs/all_candidates/aaa_term_all.gff")
#        self.assertTrue("test_folder/output/tables/all_candidates/aaa_term_all.csv")
#
#    def test_compute_intersection_forward_reverse(self):
#        self.term.multiparser = self.mock_parser
#        te.intergenic_seq = self.mock.mock_intergenic_seq
#        te.poly_t = self.mock.mock_poly_t
#        te.detect_coverage = self.mock.mock_detect_coverage
#        self.term._run_rnafold = self.mock.mock_run_rnafold
#        term_outfolder = os.path.join(self.out, "gffs")
#        csv_outfolder = os.path.join(self.out, "tables")
#        self.term._compute_intersection_forward_reverse("RNAfold_path", ["aaa"],
#                self.trans, self.test_folder, self.fastas, "cutoff_coverage", "fuzzy",
#                "wig_path", "merge_wigs", "libs", "tex_notex", "replicates", 0.5,
#                term_outfolder, csv_outfolder, True, self.gffs, self.out,
#                2, 2, 2, 2)
#        self.assertTrue(os.path.join(self.out, "inter_seq_aaa"))
#        self.assertTrue(os.path.join(self.out, "inter_sec_aaa"))
#
#    def test_compute_stat(self):
#        term_outfolder = os.path.join(self.out, "gffs")
#        csv_outfolder = os.path.join(self.out, "tables")
#        te.stat_term = self.mock.mock_stat_term
#        gen_file(os.path.join(term_outfolder, "all_candidates/aaa_term_all.gff"), self.example.term_file)
#        gen_file(os.path.join(term_outfolder, "detect/aaa_term.csv"), self.example.term_file)
#        gen_file(os.path.join(term_outfolder, "express/aaa_term.csv"), self.example.term_file)
#        self.term._compute_stat(term_outfolder, csv_outfolder, True, self.out)
#        self.assertTrue(os.path.exists(os.path.join(csv_outfolder, "express/aaa_term.csv")))
#        self.assertTrue(os.path.exists(os.path.join(csv_outfolder, "detect/aaa_term.csv")))
#
#    def test_run_terminator(self):
#        te.stat_term = self.mock.mock_stat_term
#        te.intergenic_seq = self.mock.mock_intergenic_seq
#        te.poly_t = self.mock.mock_poly_t
#        te.detect_coverage = self.mock.mock_detect_coverage
#        self.term.multiparser = self.mock_parser
#        self.term._run_rnafold = self.mock.mock_run_rnafold
#        self.term._TransTermHP = self.mock.mock_TransTermHP
#        os.mkdir(os.path.join(self.gffs, "tmp"))
#        os.mkdir(os.path.join(self.fastas, "tmp"))
#        os.mkdir(os.path.join(self.srnas, "tmp"))
#        os.mkdir(os.path.join(self.trans, "tmp"))
#        gen_file(os.path.join(self.gffs, "tmp/aaa.gff"), self.example.gff_file)
#        gen_file(os.path.join(self.fastas, "tmp/aaa.fa"), self.example.fasta_file)
#        gen_file(os.path.join(self.srnas, "tmp/aaa_sRNA.gff"), self.example.srna_file)
#        gen_file(os.path.join(self.trans, "tmp/aaa_transcript.gff"), self.example.tran_file)
#        tex_wigs = os.path.join(self.test_folder, "tex")
#        frag_wigs = os.path.join(self.test_folder, "frag")
#        os.mkdir(tex_wigs)
#        os.mkdir(frag_wigs)
#        gen_file(os.path.join(frag_wigs, "frag.wig"), "text")
#        gen_file(os.path.join(tex_wigs, "tex.wig"), "text")
#        self.term.run_terminator("TransTermHP_path", "expterm_path", "RNAfold_path",
#                       self.out, self.fastas, self.gffs, self.trans, self.srnas, True,
#                       tex_wigs, frag_wigs, 0.5, 2, 2, 2, 2, 2, 2,
#                       self.test_folder, "tlibs", "flibs", "tex_notex",
#                       2, 2, True)
#        self.assertTrue(os.path.exists(os.path.join(self.out, "tables/all_candidates")))
#        self.assertTrue(os.path.exists(os.path.join(self.out, "tables/express")))
#        self.assertTrue(os.path.exists(os.path.join(self.out, "tables/detect")))
#        self.assertTrue(os.path.exists(os.path.join(self.out, "gffs/all_candidates")))
#        self.assertTrue(os.path.exists(os.path.join(self.out, "gffs/express")))
#        self.assertTrue(os.path.exists(os.path.join(self.out, "gffs/detect")))
#

class Example(object):

    srna_file = """aaa	RefSeq	sRNA	1	5	.	+	.	ID=srna0;Name=sRNA_0"""
    gff_file = """aaa	RefSeq	CDS	10	15	.	+	.	ID=cds0;Name=CDS_0"""
    term_file = """aaa	RefSeq	CDS	17	22	.	+	.	ID=term0;Name=Term_0"""
    tran_file = """aaa	RefSeq	CDS	6	22	.	+	.	ID=tran0;Name=Transcirpt_0"""
    fasta_file = """>aaa
ATAGATAAGTTCACGATCACACATATCACTTACAATTAATTGACAGACA"""
    ptt = """Staphylococcus_aureus_HG003 - 1..2821337
2772 proteins
Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
517..1878	+	1362	YP_498609.1	dnaA	SAOUHSC_00001	-	-	chromosomal replication initiation protein
2156..3289	+	1134	YP_498610.1	-	SAOUHSC_00002	-	-	DNA polymerase III subunit beta
3670..3915	+	246	YP_498611.1	-	SAOUHSC_00003	-	-	hypothetical protein"""
    bag = """        recF NONE
SAOUHSC_00005 NONE
SAOUHSC_00006    9676 ..    9705 + -14.5  -4.91109   AAGAATAATAAAAAA       TAAGACTTCCCTA TATG TAGGGGAGTCTTA        TTTTTATGCTAGAAA 100 33
SAOUHSC_00008   12682 ..   12708 + -17.4  -6.14618   GGGTGGCAACGCGTA         GACCACGTCCCT TGT AGGGATGTGGTC         TTTTTTTATTTTCTA 100 300
SAOUHSC_00009   14156 ..   14177 +  -8.7  -5.69473   ATAATATTTTAAAAA           GTGGTGACG AAGC TGTCGCCAC            TTTTTTTGTGCTGTA 100 109"""
    gff_bag = """##gff-version 3
aaa	TransTermHP	terminator	9676	9705	.	+	.	associated_gene=SAOUHSC_00006;ID=term0;Name=terminator_00000
aaa	TransTermHP	terminator	12682	12708	.	+	.	associated_gene=SAOUHSC_00008;ID=term1;Name=terminator_00001
aaa	TransTermHP	terminator	14156	14177	.	+	.	associated_gene=SAOUHSC_00009;ID=term2;Name=terminator_00002"""

if __name__ == "__main__":
    unittest.main()
