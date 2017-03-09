import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file, import_data
from annogesiclib.screen import Screen
from mock_args_container import MockClass


class TestScreen(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.mock_args = MockClass()
        self.test_folder = "test_folder"
        self.output = os.path.join(self.test_folder, "output")
        self.tex_wig = os.path.join(self.test_folder, "tex")
        self.frag_wig = os.path.join(self.test_folder, "frag")
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
            os.mkdir(self.tex_wig)
            os.mkdir(self.frag_wig)
            os.mkdir(self.output)
        self.fasta = os.path.join(self.test_folder, "aaa.fa")
        gen_file(self.fasta, self.example.fasta)
        args = self.mock_args.mock()
        args.output_folder = self.output
        args.fasta = self.fasta
        self.screen = Screen(args)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_screenshot(self):
        gen_file(os.path.join(self.tex_wig, "tex_1_f.wig"), self.example.wig_f)
        gen_file(os.path.join(self.tex_wig, "notex_1_f.wig"),
                 self.example.wig_f)
        gen_file(os.path.join(self.frag_wig, "frag_f.wig"), self.example.wig_f)
        gen_file(os.path.join(self.tex_wig, "tex_1_r.wig"), self.example.wig_r)
        gen_file(os.path.join(self.tex_wig, "notex_1_r.wig"),
                 self.example.wig_r)
        gen_file(os.path.join(self.frag_wig, "frag_r.wig"), self.example.wig_r)
        args = self.mock_args.mock()
        args.fasta = self.fasta
        args.main_gff = os.path.join(self.test_folder, "main.gff")
        gen_file(args.main_gff, self.example.main_gff)
        side_gff = os.path.join(self.test_folder, "side.gff")
        args.side_gffs = [side_gff]
        gen_file(side_gff, self.example.side_gff)
        args.frag_wigs = self.frag_wig
        args.tex_wigs = self.tex_wig
        args.height = 1000
        args.tlibs = ["test_folder/tex/tex_1_f.wig:tex:1:a:+",
                      "test_folder/tex/tex_1_r.wig:tex:1:a:-",
                      "test_folder/tex/notex_1_f.wig:notex:1:a:+",
                      "test_folder/tex/notex_1_r.wig:notex:1:a:-"]
        args.flibs = ["test_folder/frag/frag_f.wig:frag:1:a:+",
                      "test_folder/frag/frag_r.wig:frag:1:a:-"]
        args.present = "expand"
        args.output_folder = self.output
        self.screen.screenshot(args)
        self.assertTrue(os.path.exists(os.path.join(
            self.output, "screenshots", "aaa", "forward")))
        self.assertTrue(os.path.exists(os.path.join(
            self.output, "screenshots", "aaa", "reverse")))
        datas = import_data(os.path.join(
            self.output, "screenshots", "aaa", "forward.txt"))
        datas = import_data(os.path.join(
            self.output, "screenshots", "aaa", "reverse.txt"))
        self.assertEqual("\n".join(datas), self.example.out_r)

    def test_import_libs(self):
        texs = [[os.path.join(self.tex_wig, "tex_1.wig"),
                 "tex", "1", "a", "+"],
                [os.path.join(self.tex_wig, "notex_1.wig"),
                 "notex", "1", "a", "+"]]
        lib_dict = {"ft": [], "fn": [], "rt": [],
                    "rn": [], "ff": [], "rf": []}
        self.screen._import_libs(texs, "+", lib_dict)
        self.assertDictEqual(lib_dict, {
            'fn': ['test_folder/tex/notex_1.wig'], 'rn': [],
            'rt': [], 'ft': ['test_folder/tex/tex_1.wig'],
            'rf': [], 'ff': []})


class Example(object):

    main_gff = """aaa\tRefseq\tCDS\t14\t18\t.\t+\t.\tID=cds0;Name=CDS_0"""
    side_gff = """aaa\tRefseq\tTSS\t13\t13\t.\t+\t.\tID=tss0;Name=TSS_0"""
    fasta = """>aaa
ATATAGATCCCGCTACAATGCACATCGAGTAACTGGA"""
    wig_f = """track type=wiggle_0 name="Hp26695_ML_B1_HS1_-TEX_forward"
variableStep chrom=aaa span=1
4 0.003744270207251674
5 0.003744270207251674
6 0.003744270207251674
7 0.008670941532582825
8 0.008670941532582825
9 0.008670941532582825
10 0.008670941532582825
11 0.0124152117398345
12 0.0124152117398345
13 0.0124152117398345
14 0.01734188306516565
15 0.01734188306516565
16 0.01734188306516565
17 0.01734188306516565
18 0.021086153272417325
19 0.021086153272417325
20 0.01734188306516565"""
    wig_r = """track type=wiggle_0 name="Hp26695_ML_B1_HS1_-TEX_forward"
variableStep chrom=aaa span=1
4 -0.003744270207251674
5 -0.003744270207251674
6 -0.003744270207251674
7 -0.008670941532582825
8 -0.008670941532582825
9 -0.008670941532582825
10 -0.008670941532582825
11 -0.0124152117398345
12 -0.0124152117398345
13 -0.0124152117398345
14 -0.01734188306516565
15 -0.01734188306516565
16 -0.01734188306516565
17 -0.01734188306516565
18 -0.021086153272417325
19 -0.021086153272417325
20 -0.01734188306516565"""
    out_f = """new
genome /home/silas/ANNOgesic/test_folder/aaa.fa
load /home/silas/ANNOgesic/test_folder/main.gff
expand main.gff
load /home/silas/ANNOgesic/test_folder/side.gff
expand side.gff
load /home/silas/ANNOgesic/test_folder/tex/tex_1_f.wig
load /home/silas/ANNOgesic/test_folder/tex/notex_1_f.wig
load /home/silas/ANNOgesic/test_folder/frag/frag_f.wig
maxPanelHeight 1000
snapshotDirectory /home/silas/ANNOgesic/test_folder/output/aaa/forward
goto aaa:-186-218
setDataRange 0,10
snapshot aaa:14-18.png"""
    out_f = out_f.replace("/home/silas/ANNOgesic", os.getcwd())
    out_r = """new
genome /home/silas/ANNOgesic/test_folder/aaa.fa
load /home/silas/ANNOgesic/test_folder/main.gff
expand main.gff
load /home/silas/ANNOgesic/test_folder/side.gff
expand side.gff
load /home/silas/ANNOgesic/test_folder/tex/tex_1_r.wig
load /home/silas/ANNOgesic/test_folder/tex/notex_1_r.wig
load /home/silas/ANNOgesic/test_folder/frag/frag_r.wig
maxPanelHeight 1000
snapshotDirectory /home/silas/ANNOgesic/test_folder/output/aaa/reverse"""
    out_r = out_r.replace("/home/silas/ANNOgesic", os.getcwd())

if __name__ == "__main__":
    unittest.main()
