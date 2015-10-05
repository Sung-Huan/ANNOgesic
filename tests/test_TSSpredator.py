#!/usr/bin/python

import os
import sys
import csv
import shutil
import unittest
from io import StringIO
sys.path.append(".")
from annogesiclib.TSSpredator import TSSPredatorReader


class TestTSSPredatorReader(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
        self.tss = TSSPredatorReader()

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_entries(self):
        input_fh = StringIO(self.example.master)
        tsss = []
        for entry in self.tss.entries(input_fh):
            tsss.append(entry)
        self.assertEqual(tsss[0].pos, 179)
        self.assertTrue(tsss[1].is_primary)
        self.assertTrue(tsss[2].is_internal)


class Example(object):

    master = """SuperPos	SuperStrand	mapCount	detCount	Genome	detected	enriched	stepHeight	stepFactor	enrichmentFactor	classCount	Pos	Strand	Locus_tag	sRNA/asRNA	Product	UTRlength	GeneLength	Primary	Secondary	Internal	Antisense	Automated	Manual	Putative sRNA	Putative asRNA	Comment	Sequence -50 nt upstream + TSS (51nt)
179	-	1	1	test	1	1	4.45	31.93	8.69	1	179	-	orphan		orphan	NA	NA	0	0	0	0	1	0	0	0		ACCCTTGAATTGAGGGTGTTTTATACCTAAATTTAAAAAATGATGCTATAA
681	-	1	1	test	1	1	4.2	3.0	3.54	2	681	-	HP0001		transcription antitermination protein NusB	48	417	1	0	0	0	1	0	0	0		GATTGAAAGAGCGGGCAGTAAAGCCGGCAATAAGGGCTTTGAAGCGATGAG
681	-	1	1	test	1	1	4.2	3.0	3.54	2	681	-	HP0002		6%2C7-dimethyl-8-ribityllumazine synthase	NA	471	0	0	1	0	1	0	0	0		GATTGAAAGAGCGGGCAGTAAAGCCGGCAATAAGGGCTTTGAAGCGATGAG"""
if __name__ == "__main__":
    unittest.main()

