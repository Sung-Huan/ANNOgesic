#!/usr/bin/python

import os
import sys
import csv
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from annogesiclib.parser_wig import WigParser


class TestParserWig(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.wig_parser = WigParser()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_parser(self):
        wigs = []
        wig_f_fh = StringIO(self.example.wig_forward_file)
        for entry in self.wig_parser.parser(wig_f_fh, "+"):
            self.assertEqual(entry.strain, "aaa")
            self.assertEqual(entry.track, "TSB_t0_TEX_forward")
            wigs.append(entry)
        self.assertEqual(wigs[2].pos, 3)
        self.assertEqual(wigs[2].coverage, 1.4041251228308191)
        wigs = []
        wig_r_fh = StringIO(self.example.wig_reverse_file)
        for entry in self.wig_parser.parser(wig_r_fh, "-"):
            self.assertEqual(entry.strain, "aaa")
            self.assertEqual(entry.track, "TSB_t0_TEX_reverse")
            wigs.append(entry)
        self.assertEqual(wigs[2].pos, 3)
        self.assertEqual(wigs[2].coverage, 1.4041251228308191)
 
class Example(object):
    wig_forward_file = """track type=wiggle_0 name="TSB_t0_TEX_forward"
variableStep chrom=aaa span=1
3 1.4041251228308191
4 56.867067474648174
5 56.867067474648174"""

    wig_reverse_file = """track type=wiggle_0 name="TSB_t0_TEX_reverse"
variableStep chrom=aaa span=1
3 -1.4041251228308191
4 -56.867067474648174
5 -56.867067474648174"""

if __name__ == "__main__":
    unittest.main()

