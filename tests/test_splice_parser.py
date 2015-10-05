#!/usr/bin/python

import os
import sys
import csv
import shutil
import unittest
from io import StringIO
sys.path.append(".")
from annogesiclib.splice_parser import SpliceParser


class TestGff3Parser(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.s_parser = SpliceParser()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_parser(self):
        splice_fh = StringIO(self.example.splice)
        starts = []
        splices = []
        for entry in self.s_parser.parser(splice_fh):
            starts.append(entry.start)
            splices.append(entry.splice)
        self.assertListEqual(starts, [17647, 20734, 43490, 49952])
        self.assertListEqual(splices, ['splits:1:1:1:N:F', 'splits:1:1:1:C:P',
                                       'splits:1:1:1:N:P', 'splits:2:2:2:N:P'])

class Example(object):

    splice = """Staphylococcus_aureus_HG003	17647	17667	splits:1:1:1:N:F	0	+
Staphylococcus_aureus_HG003	20734	21396	splits:1:1:1:C:P	0	+
Staphylococcus_aureus_HG003	43490	43644	splits:1:1:1:N:P	0	+
Staphylococcus_aureus_HG003	49952	50016	splits:2:2:2:N:P	0	+"""

if __name__ == "__main__":
    unittest.main()

