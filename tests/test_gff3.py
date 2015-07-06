#!/usr/bin/python

import os
import sys
import csv
import shutil
import unittest
from io import StringIO
sys.path.append(".")
from annogesiclib.gff3 import Gff3Parser


class TestGff3Parser(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.gff_parser = Gff3Parser()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_gff_parser(self):
        strains = []
        features = []
        starts = []
        ends = []
        strands = []
        IDs = []
        fh = StringIO(self.example.gff_file)
        for entry in self.gff_parser.entries(fh):
            strains.append(entry.seq_id)
            features.append(entry.feature)
            starts.append(entry.start)
            ends.append(entry.end)
            strands.append(entry.strand)
            IDs.append(entry.attributes["ID"])
        self.assertListEqual(strains, ["aaa", "aaa", "aaa", "aaa", "bbb", "bbb"])
        self.assertListEqual(features, ["gene", "CDS", "gene", "CDS", "gene", "tRNA"])
        self.assertListEqual(starts, [517, 517, 2156, 2156, 4444, 4444])
        self.assertListEqual(ends, [1878, 1878, 3289, 3289, 5444, 5444])
        self.assertListEqual(strands, ["+", "+", "-", "-", "+", "+"])
        self.assertListEqual(IDs, ["gene0", "cds0", "gene1", "cds1", "gene2", "rna0"])

class Example(object):

    gff_file = """#gff3
aaa	Refseq	gene	517	1878	.	+	.	Name=dnaA;locus_tag=AAA_00001;gene=dnaA;ID=gene0;db_xref=GeneID:3919798
aaa	Refseq	CDS	517	1878	.	+	.	protein_id=YP_498609.1;ID=cds0;Name=YP_498609.1;product=chromosomal replication initiation protein;Parent=gene0
aaa	Refseq	gene	2156	3289	.	-	.	Name=AAA_00002;locus_tag=AAA_00002;ID=gene1;db_xref=GeneID:3919799
aaa	Refseq	CDS	2156	3289	.	-	.	protein_id=YP_498610.1;ID=cds1;Name=YP_498610.1;locus_tag=AAA_00002
bbb	Refseq	gene	4444	5444	.	+	.	Name=AAA_T00004;locus_tag=AAA_T00004;ID=gene2
bbb	Refseq	tRNA	4444	5444	.	+	.	Name=AAA_T00018;locus_tag=AAA_T00004;ID=rna0"""

if __name__ == "__main__":
    unittest.main()

