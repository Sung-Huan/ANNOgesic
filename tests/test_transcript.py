import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file, import_data, extract_info
import annogesiclib.transcript as tr
from annogesiclib.transcript import TranscriptAssembly

class Mock_func(object):

    def mock_assembly(self, wig_f, wig_r, height, width, tolerance, low_cutoff,
                      wig_folder, tex, libs, replicates, out, wig_type):
        pass

    def mock_combine(self, frag, tex, tolerance, out):
        gen_file(out, "test")
        pass

    def mock_fill_gap(self, gff, tran, overlap, tmp_overlap):
        pass

    def mock_longer_ta(self, tmp_out, length, final_out):
        pass

    def mock_stat_ta_gff(self, ta_file, cds_file, stat_gff_out,
                         tmp_ta, tmp_gff, c_feature):
        pass

    def mock_stat_ta_tss(self, ta_file, tss_file, stat_tss_out, ta, tss, fuzzy):
        pass

    def mock_gen_table_tran(self, gff_outfolder, frag_wigs, tex_wigs,
                            tlibs, flibs, table_best, out_folder, gffs):
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


class TestsTranscriptAssembly(unittest.TestCase):

    def setUp(self):
        self.mock = Mock_func()
        self.mock_parser = Mock_Multiparser()
        self.example = Example()
        self.test_folder = "test_folder"
        self.trans = "test_folder/trans"
        self.out = "test_folder/output"
        self.tex = "test_folder/tex"
        self.frag = "test_folder/frag"
        self.gffs = "test_folder/gffs"
        self.tsss = "test_folder/tsss"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
            os.mkdir(self.trans)
            os.mkdir(self.out)
            os.mkdir(self.tex)
            os.mkdir(self.frag)
            os.mkdir(self.gffs)
            os.mkdir(self.tsss)
        self.tran = TranscriptAssembly(self.out)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_combine_wigs(self):
        gen_file(os.path.join(self.frag, "aaa_forward.wig_STRAIN_test.wig"), "test")
        gen_file(os.path.join(self.frag, "aaa_reverse.wig_STRAIN_test.wig"), "test")
        libs = ["aaa_forward.wig_STRAIN_test.wig:frag:1:a:+", "aaa_reverse.wig_STRAIN_test.wig:frag:1:a:-"]
        self.tran._combine_wigs(self.frag, self.test_folder, "bbb_forward.wig", libs, "test")
        self.assertTrue(os.path.exists(os.path.join(self.test_folder, "bbb_forward.wig")))

    def test_compute(self):
        tr.assembly = self.mock.mock_assembly
        os.mkdir(os.path.join(self.frag, "tmp"))
        gen_file(os.path.join(self.frag, "tmp/aaa_forward_STRAIN_test.wig"), "test")
        strains = self.tran._compute(self.test_folder, self.frag, "tex", 5, 20, 5,
                                     "replicates", self.out, "frag", "libs", 0)
        self.assertListEqual(strains, ['test'])

    def test_for_one_wig(self):
        tr.assembly = self.mock.mock_assembly
        self.tran.multiparser = self.mock_parser
        os.mkdir(os.path.join(self.frag, "tmp"))
        gen_file(os.path.join(self.frag, "tmp/aaa_forward_STRAIN_test.wig"), "test")
        gen_file(os.path.join(self.out, "test_frag"), self.example.tran_file)
        strains = self.tran._for_one_wig("frag", self.frag, "tex", 5, 20, 5,
                                         "replicates", self.out, "libs", self.gffs, 0)
        self.assertListEqual(strains, ['test'])
        datas = import_data(os.path.join(self.gffs, "test_transcript_assembly_frag.gff"))
        self.assertEqual("\n".join(datas), "##gff-version 3\n" + self.example.tran_file)

    def test_for_two_wigs(self):
        tr.combine = self.mock.mock_combine
        gen_file(os.path.join(self.gffs, "test_transcript_assembly_fragment.gff"), "test")
        gen_file(os.path.join(self.gffs, "test_transcript_assembly_tex_notex.gff"), "test")
        
        self.tran._for_two_wigs(self.frag, self.tex, ["test"],
                                self.gffs, 5)
        self.assertTrue(os.path.exists(os.path.join(self.gffs, "test_transcript.gff")))

    def test_post_modify(self):
        tr.longer_ta = self.mock.mock_longer_ta
        tr.fill_gap = self.mock.mock_fill_gap
        gen_file(os.path.join(self.gffs, "test.gff"), self.example.gff_file)
        gff_out = os.path.join(self.out, "gffs")
        os.mkdir(gff_out)
        os.mkdir(os.path.join(self.out, "tmp_tran"))
        gen_file(os.path.join(gff_out, "tmp_uni"), self.example.tran_file)
        gen_file(os.path.join(gff_out, "tmp_overlap"), self.example.tran_file)
        gen_file(os.path.join(gff_out, "final_test"), self.example.tran_file)
        self.tran._post_modify(["test"], self.gffs, self.trans, 20,
                               gff_out, self.out)
        self.assertTrue(os.path.exists(os.path.join(gff_out, "test_transcript.gff")))

    def test_compare_cds(self):
        tr.stat_ta_gff = self.mock.mock_stat_ta_gff
        self.tran.multiparser = self.mock_parser
        os.mkdir(os.path.join(self.gffs, "tmp"))
        gen_file(os.path.join(self.gffs, "test.gff"), self.example.gff_file)
        gen_file(os.path.join(self.gffs, "tmp/test.gff"), self.example.gff_file)
        gen_file(os.path.join(self.out, "test_transcript.gff"), self.example.tran_file)
        gff_out = os.path.join(self.out, "gffs")
        os.mkdir(gff_out)
        gen_file(os.path.join(gff_out, "tmp_ta_gff"), self.example.tran_file)
        gen_file(os.path.join(gff_out, "tmp_gff_ta"), self.example.gff_file)
        self.tran._compare_cds(["test"], self.out, self.gffs, self.trans, "CDS")
        datas = import_data(os.path.join(self.gffs, "test.gff"))
        self.assertEqual("\n".join(datas), "##gff-version 3\n" + self.example.gff_file)
        datas = import_data(os.path.join(self.out, "test_transcript.gff"))
        self.assertEqual("\n".join(datas), "##gff-version 3\n" + self.example.tran_file)

    def test_compare_tss(self):
        tr.stat_ta_tss = self.mock.mock_stat_ta_tss
        self.tran.multiparser = self.mock_parser
        os.mkdir(os.path.join(self.gffs, "tmp"))
        gen_file(os.path.join(self.gffs, "test_TSS.gff"), self.example.gff_file)
        gen_file(os.path.join(self.gffs, "tmp/test_TSS.gff"), self.example.gff_file)
        gen_file(os.path.join(self.out, "test_transcript.gff"), self.example.tran_file)
        gff_out = os.path.join(self.out, "gffs")
        os.mkdir(gff_out)
        gen_file(os.path.join(gff_out, "tmp_ta_tss"), self.example.tran_file)
        gen_file(os.path.join(gff_out, "tmp_tss_ta"), self.example.gff_file)
        self.tran._compare_tss(["test"], self.out, self.gffs, self.trans, 2)
        datas = import_data(os.path.join(self.gffs, "test_TSS.gff"))
        self.assertEqual("\n".join(datas), "##gff-version 3\n" + self.example.gff_file)
        datas = import_data(os.path.join(self.out, "test_transcript.gff"))
        self.assertEqual("\n".join(datas), "##gff-version 3\n" + self.example.tran_file)

    def test_run_transcript_assembly(self):
        tr.stat_ta_tss = self.mock.mock_stat_ta_tss
        tr.stat_ta_gff = self.mock.mock_stat_ta_gff
        tr.longer_ta = self.mock.mock_longer_ta
        tr.fill_gap = self.mock.mock_fill_gap
        tr.combine = self.mock.mock_combine
        tr.assembly = self.mock.mock_assembly
        self.tran._gen_table = self.mock.mock_gen_table_tran
        gen_file(os.path.join(self.frag, "test1_forward.wig"), self.example.wig_f)
        gen_file(os.path.join(self.frag, "test1_reverse.wig"), self.example.wig_r)
        gen_file(os.path.join(self.gffs, "test.gff"), self.example.gff_file)
        gen_file(os.path.join(self.tsss, "test_TSS.gff"), self.example.tss_file)
        gen_file("test_folder/output/test_fragment", self.example.tran_file)
        gff_out = os.path.join(self.out, "gffs")
        os.mkdir(gff_out)
        gen_file(os.path.join(gff_out, "test_transcript_assembly_fragment.gff"), self.example.tran_file)
        gen_file(os.path.join(gff_out, "tmp_uni"), self.example.tran_file)
        gen_file(os.path.join(gff_out, "tmp_overlap"), self.example.tran_file)
        gen_file(os.path.join(gff_out, "final_test"), self.example.tran_file)
        self.tran.run_transcript_assembly(self.frag, None, True, "tex", 1,
                                          self.gffs, 5, 1, 5, 0, 1, 1, self.out,
                                          None, None, 5, "tlibs", "flibs", "CDS", False)
        self.assertTrue(os.path.exists(os.path.join(gff_out, "test_transcript.gff")))


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


if __name__ == "__main__":
    unittest.main()
