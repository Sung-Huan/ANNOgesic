#!/usr/bin/python

import os
import sys
import copy
import shutil
import unittest
from io import StringIO
sys.path.append(".")
from mock_helper import import_data
from annogesiclib.helper import Helper


class TestHelper(unittest.TestCase):

    def setUp(self):
        self.example = ExampleData()
        self.helper = Helper()
        self.gff_out = self.example.gff_out
        self.rev_seq = self.example.rev_seq.replace("\n", "")
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)        
        self.gff_file = os.path.join(self.test_folder, "test.gff")
        with open(self.gff_file, "w") as rh:
            rh.write(self.example.gff_file)        
        self.seq_file = os.path.join(self.test_folder, "test.fa")
        with open(self.seq_file, "w") as rh:
            rh.write(self.example.seq)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_remove_all_content(self):
        tmp1 = os.path.join(self.test_folder, "tmp1.gff")
        tmp2 = os.path.join(self.test_folder, "tmp2")
        shutil.copyfile(self.gff_file, tmp1)
        os.mkdir(tmp2)
        self.helper.remove_all_content(self.test_folder, "tmp", "file")
        self.assertFalse(os.path.exists(tmp1))
        self.assertTrue(os.path.exists(tmp2))
        self.helper.remove_all_content(self.test_folder, "tmp", "dir")
        self.assertFalse(os.path.exists(tmp2))
        self.assertTrue(os.path.exists(self.gff_file))

    def test_remove_tmp(self):
        tmp1 = os.path.join(self.test_folder, "tmp")
        tmp2 = os.path.join(self.test_folder, "test.gff_folder")
        os.mkdir(tmp1)
        os.mkdir(tmp2)
        self.helper.remove_tmp(self.test_folder)
        self.assertFalse(os.path.exists(tmp1))
        self.assertFalse(os.path.exists(tmp2))

    def test_get_correct_file(self):
        gff_file = os.path.join(self.test_folder, "test.gff")
        wig_f_file = os.path.join(self.test_folder, "test_forward.wig_STRAIN_aaa.wig")
        wig_r_file = os.path.join(self.test_folder, "test_reverse.wig_STRAIN_aaa.wig")
        shutil.copyfile(gff_file, wig_f_file)
        shutil.copyfile(gff_file, wig_r_file)
        filename = self.helper.get_correct_file(self.test_folder, ".gff", "test", None)
        self.assertEqual(filename, gff_file)
        
    def test_sorf_gff(self):
        out_file = os.path.join(self.test_folder, "test.out")
        self.helper.sort_gff(self.gff_file, out_file)
        datas = import_data(out_file)
        self.assertEqual(set(datas), set(self.gff_out.split("\n")))

    def test_extract_gene(self):
        seq = self.example.seq.replace("\n", "")
        new_seq = self.helper.extract_gene(seq, 1, 70, "+")
        self.assertEqual(new_seq,
        "CGCAGGTTGAGTTCCTGTTCCCGATAGATCCGATAAACCCGCTTATGATTCCAGAGCTGTCCCTGCACAT")
        new_seq = self.helper.extract_gene(seq, 1, 140, "-")
        self.assertEqual(new_seq, self.rev_seq)

    def test_get_seq(self):
        gff_file = os.path.join(self.test_folder, "test.gff")
        out_file = os.path.join(self.test_folder, "test.cds")
        lines = self.example.gff_out.split("\n")
        with open(gff_file, "w") as gh:
            gh.write(lines[1])
        self.helper.get_seq(self.gff_file, self.seq_file, out_file)
        datas = import_data(out_file)
        self.assertEqual(set(datas), set([">cds0|aaa|1|10|+", "CGCAGGTTGA"]))

class ExampleData(object):

    gff_file = """##gff-version 3
#!gff-spec-version 1.20
#!processor NCBI annotwriter
##sequence-region ddd 1 2329769
##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=591001
ddd	Refseq	CDS	921	3443	.	+	.	locus_tag=SAOUHSC_00001;Name=dnaA;db_xref=GeneID:3919798;ID=cds0;gene=dnaA
aaa	Refseq	CDS	1	10	.	+	.	locus_tag=SAOUHSC_00001;Name=dnaR;db_xref=GeneID:3919434;ID=cds0;gene=dnaR
aaa	Refseq	CDS	1	5	.	-	.	locus_tag=SAOUHSC_00001;Name=dnaR;db_xref=GeneID:3919434;gene=dnaR
ddd	Refseq	CDS	517	1878	.	-	.	gene=dnaA;transl_table=11;note=binds to the dnaA-box as an ATP-bound complex at the origin of replication during the initiation of chromosomal replication, can also affect transcription of multiple genes including itself.;ID=cds1;locus_tag=SAOUHSC_00004;db_xref=GeneID:3919798"""

    gff_out = """##gff-version 3
aaa	Refseq	CDS	1	10	.	+	.	locus_tag=SAOUHSC_00001;Name=dnaR;db_xref=GeneID:3919434;ID=cds0;gene=dnaR
aaa	Refseq	CDS	1	5	.	-	.	locus_tag=SAOUHSC_00001;Name=dnaR;db_xref=GeneID:3919434;gene=dnaR
ddd	Refseq	CDS	517	1878	.	-	.	gene=dnaA;transl_table=11;note=binds to the dnaA-box as an ATP-bound complex at the origin of replication during the initiation of chromosomal replication, can also affect transcription of multiple genes including itself.;ID=cds1;locus_tag=SAOUHSC_00004;db_xref=GeneID:3919798
ddd	Refseq	CDS	921	3443	.	+	.	locus_tag=SAOUHSC_00001;Name=dnaA;db_xref=GeneID:3919798;ID=cds0;gene=dnaA"""

    seq = """
CGCAGGTTGAGTTCCTGTTCCCGATAGATCCGATAAACCCGCTTATGATTCCAGAGCTGTCCCTGCACAT
GCCGCTTCATTTTTTTCCAAGGGCTTCCTTCAGGATATCCGTCTGCATGCTCAAATCCGCATACATGCGC"""

    rev_seq = """
GCGCATGTATGCGGATTTGAGCATGCAGACGGATATCCTGAAGGAAGCCCTTGGAAAAAAATGAAGCGGC
ATGTGCAGGGACAGCTCTGGAATCATAAGCGGGTTTATCGGATCTATCGGGAACAGGAACTCAACCTGCG"""

if __name__ == "__main__":
    unittest.main()
