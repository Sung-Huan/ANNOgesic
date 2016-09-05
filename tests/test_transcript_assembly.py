import sys
import math
import csv
import os
import shutil
import unittest
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data, gen_file
import annogesiclib.transcript_assembly as ta
from mock_args_container import MockClass


class TestTranscriptAssembly(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.mock_args = MockClass()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_read_wig(self):
        libs = [{"name": "test1", "type": "frag",
                 "cond": "1", "strand": "+", "rep": "a"}]
        filename = os.path.join(self.test_folder, "test_f.wig")
        gen_file(filename, self.example.wig_f)
        wigs = ta.read_wig(filename, libs, "+")
        self.assertDictEqual(wigs, self.example.wigs_f)

    def test_detect_hight_toler(self):
        cover = {"coverage": 100, "track": "test_1"}
        height = 5
        tmp_covers = {"best": 10, "toler": 2}
        tracks = []
        ta.detect_hight_toler(cover, height, tmp_covers, tracks)
        self.assertDictEqual(tmp_covers, {'best': 100, 'toler': 2})

    def test_check_tex_conds(self):
        check_tex = []
        tracks = ["test1", "test2"]
        libs = [{"name": "test1", "type": "frag",
                 "cond": "1", "strand": "+", "rep": "a"},
                {"name": "test2", "type": "tex",
                 "cond": "2", "strand": "+", "rep": "a"}]
        texs = {"test1": 2, "test2": 2}
        conds = {}
        ta.check_tex_conds(tracks, libs, texs, check_tex, conds, 1)
        self.assertDictEqual(conds, {'1_frag': 1, '2_tex': 1})

    def test_elongation(self):
        covers = [{"coverage": 10, "pos": 10, "track": "test1"},
                  {"coverage": 1, "pos": 10, "track": "test2"},
                  {"coverage": 100, "pos": 11, "track": "test1"},
                  {"coverage": 20, "pos": 11, "track": "test2"}]
        libs = [{"name": "test1", "type": "tex",
                 "cond": "1", "strand": "+", "rep": "a"},
                {"name": "test2", "type": "notex",
                 "cond": "1", "strand": "+", "rep": "a"}]
        reps = {"tex": 1, "frag": 1}
        tmp_texs = {"test1_test2": 2}
        tolers = []
        trans = {"aaa": []}
        args = self.mock_args.mock()
        args.replicates = reps
        args.height = 5
        args.tex = 2
        cover_best, conds, tracks, texs, pos = ta.elongation(covers, tmp_texs,
                                               libs, "+", trans, args, "aaa", tolers)
        self.assertEqual(cover_best, 100)
        self.assertListEqual(tracks, ['test1', 'test2'])
        self.assertDictEqual(texs, {'test1_test2': 2})
        self.assertEqual(pos, 11)
        self.assertDictEqual(trans, {'aaa': [{'coverage': 10, 'cond': 1, 'strand': '+', 'pos': 10}]})

    def test_transfer_to_tran(self):
        reps = {"tex": 1, "frag": 1}
        tmp_texs = {"test1": 2}
        libs = [{"name": "test1", "type": "frag",
                 "cond": "1", "strand": "+", "rep": "a"}]
        args = self.mock_args.mock()
        args.height = 10
        args.tex = 1
        args.replicates = reps
        tolers, trans = ta.transfer_to_tran(self.example.wigs_f, libs, tmp_texs, "+", args)
        self.assertDictEqual(tolers, {'aaa': [0.0, 2.0, 20, 20, 4.0, 20, 7.0]})
        self.assertDictEqual(trans, {'aaa': [{'pos': 3, 'cond': 1, 'strand': '+', 'coverage': 41.0},
                                             {'pos': 4, 'cond': 1, 'strand': '+', 'coverage': 47.0},
                                             {'pos': 6, 'cond': 1, 'strand': '+', 'coverage': 47.0},
                                             {'pos': 8, 'cond': 1, 'strand': '+', 'coverage': 47.0}]})

    def test_fill_gap_and_print(self):
        trans = {'aaa': [{'pos': 3, 'cond': 1, 'strand': '+', 'coverage': 41.0},
                         {'pos': 4, 'cond': 1, 'strand': '+', 'coverage': 47.0},
                         {'pos': 6, 'cond': 1, 'strand': '+', 'coverage': 47.0},
                         {'pos': 8, 'cond': 1, 'strand': '+', 'coverage': 47.0}]}
        out = StringIO()
        tolers = {'aaa': [0.0, 2.0, 20, 20, 4.0, 20, 20]}
        args = self.mock_args.mock()
        args.tolerance = 3
        args.low_cutoff = 5
        args.width = 1
        ta.fill_gap_and_print(trans, "+", out, tolers, "TEX", args)
        self.assertEqual(out.getvalue(), self.example.out_tran + "\n")

    def test_print_transctipt(self):
        out = StringIO()
        ta.print_transctipt(100, 200, 20, 1, 40, "TEX",
                            20, out, "aaa", "+")
        self.assertEqual(out.getvalue(), "aaa\tANNOgesic\ttranscript\t100\t200\t.\t+\t.\tID=tran_1;Name=transcript_00001;high_coverage=40;low_coverage=20;detect_lib=TEX\n")


    def test_assembly(self):
        wig_f_file = os.path.join(self.test_folder, "aaa_forward.wig")
        wig_r_file = os.path.join(self.test_folder, "aaa_reverse.wig")
        wig_f2_file = os.path.join(self.test_folder, "aaa2_forward.wig")
        wig_r2_file = os.path.join(self.test_folder, "aaa2_reverse.wig")
        gen_file(wig_f_file, self.example.wig_f)
        gen_file(wig_r_file, self.example.wig_r)
        gen_file(wig_f2_file, self.example.wig_f)
        gen_file(wig_r2_file, self.example.wig_r)
        reps = {"tex": 1, "frag": 1}
        out_file = os.path.join(self.test_folder, "out")
        input_lib = ["aaa_forward.wig:frag:1:a:+",
                     "aaa_reverse.wig:frag:1:a:-",
                     "aaa2_forward.wig:tex:1:a:+",
                     "aaa2_reverse.wig:tex:1:a:-"]
        args = self.mock_args.mock()
        args.replicates = reps
        args.height = 10
        args.width = 1
        args.tolerance = 3
        args.tex = 2
        args.low_cutoff = 5
        ta.assembly(wig_f_file, wig_r_file, self.test_folder, input_lib,
                    out_file, "TEX", args)
        datas = import_data(out_file)
        self.assertEqual("\n".join(datas), "##gff-version 3\n" + self.example.out_tran)


class Example(object):

    wig_f = """track type=wiggle_0 name="test1"
variableStep chrom=aaa span=1
2 2.0
3 41.0
4 47.0
5 4.0
6 47.0
7 7.0
8 47.0
track type=wiggle_0 name="test2"
variableStep chrom=aaa span=1
2 2.0
3 41.0
4 47.0
5 4.0
6 47.0
7 7.0
8 47.0"""
    wig_r = """track type=wiggle_0 name="test1"
variableStep chrom=aaa span=1
2 -2.0
3 -1.0
4 -7.0
5 -4.0
6 -7.0
7 -7.0
8 -7.0
track type=wiggle_0 name="test2"
variableStep chrom=aaa span=1
2 -2.0
3 -1.0
4 -7.0
5 -4.0
6 -7.0
7 -7.0
8 -7.0"""
    wigs_f = {'aaa': [{'cond': '1', 'track': 'test1', 'pos': 1, 'strand': '+', 'coverage': 0.0},
                      {'cond': '1', 'track': 'test1', 'pos': 2, 'strand': '+', 'coverage': 2.0},
                      {'cond': '1', 'track': 'test1', 'pos': 3, 'strand': '+', 'coverage': 41.0},
                      {'cond': '1', 'track': 'test1', 'pos': 4, 'strand': '+', 'coverage': 47.0},
                      {'cond': '1', 'track': 'test1', 'pos': 5, 'strand': '+', 'coverage': 4.0},
                      {'cond': '1', 'track': 'test1', 'pos': 6, 'strand': '+', 'coverage': 47.0},
                      {'cond': '1', 'track': 'test1', 'pos': 7, 'strand': '+', 'coverage': 7.0},
                      {'cond': '1', 'track': 'test1', 'pos': 8, 'strand': '+', 'coverage': 47.0}]}
    out_tran = """aaa	ANNOgesic	transcript	3	4	.	+	.	ID=tran_0;Name=transcript_00000;high_coverage=47.0;low_coverage=41.0;detect_lib=TEX
aaa	ANNOgesic	transcript	6	8	.	+	.	ID=tran_1;Name=transcript_00001;high_coverage=47.0;low_coverage=47.0;detect_lib=TEX"""


if __name__ == "__main__":
    unittest.main()
