#!/usr/bin/python

import os
import csv
import shutil
import sys
import unittest
from io import StringIO
sys.path.append(".")
from annogesiclib.multiparser import Multiparser


class TestMultiparser(unittest.TestCase):

    def setUp(self):
        self.multiparser = Multiparser()
        self.example = Example()
        self.ref_folder = "ref_folder"
        if (not os.path.exists(self.ref_folder)):
            os.mkdir(self.ref_folder)
        self.tar_folder = "tar_folder"
        if (not os.path.exists(self.tar_folder)):
            os.mkdir(self.tar_folder)

    def tearDown(self):
        if os.path.exists(self.ref_folder):
            shutil.rmtree(self.ref_folder)
        if os.path.exists(self.tar_folder):
            shutil.rmtree(self.tar_folder)

    def test_combine_fasta(self):
        tmp_tar = os.path.join(self.tar_folder, "tmp")
        tmp_ref = os.path.join(self.ref_folder, "test.gff_folder")
        os.mkdir(tmp_ref)
        os.mkdir(tmp_tar)
        sub_fasta1 = os.path.join(tmp_tar, "aaa.fa")
        with open(sub_fasta1, "w") as rh:
            rh.write(self.example.sub_fasta1)
        sub_fasta2 = os.path.join(tmp_tar, "bbb.fa")
        with open(sub_fasta2, "w") as rh:
            rh.write(self.example.sub_fasta2)
        sub_gff1 = os.path.join(tmp_ref, "aaa.gff")
        with open(sub_gff1, "w") as rh:
            rh.write(self.example.sub_gff1)
        sub_gff2 = os.path.join(tmp_ref, "bbb.gff")
        with open(sub_gff2, "w") as rh:
            rh.write(self.example.sub_gff2)
        self.multiparser.combine_fasta(self.ref_folder, tmp_tar, None)
        self.assertTrue(os.path.exists(os.path.join(tmp_tar, "test.fa")))

    def test_combine_wig(self):
        tmp_tar = os.path.join(self.tar_folder, "tmp")
        tmp_ref = os.path.join(self.ref_folder, "test.fa_folder")
        os.mkdir(tmp_ref)
        os.mkdir(tmp_tar)
        sub_fasta1 = os.path.join(tmp_ref, "aaa.fa")
        with open(sub_fasta1, "w") as rh:
            rh.write(self.example.sub_fasta1)
        sub_fasta2 = os.path.join(tmp_ref, "bbb.fa")
        with open(sub_fasta2, "w") as rh:
            rh.write(self.example.sub_fasta2)
        sub_wig1 = os.path.join(tmp_tar, "test_forward.wig_STRAIN_aaa.wig")
        sub_wig2 = os.path.join(tmp_tar, "test_forward.wig_STRAIN_bbb.wig")
        sub_wig3 = os.path.join(tmp_tar, "test_reverse.wig_STRAIN_aaa.wig")
        sub_wig4 = os.path.join(tmp_tar, "test_reverse.wig_STRAIN_bbb.wig")
        wig_files = [sub_wig1, sub_wig2, sub_wig3, sub_wig4]
        example_wigs = [self.example.sub_f_wig1, self.example.sub_f_wig2,
                        self.example.sub_r_wig1, self.example.sub_r_wig2]
        for index in range(0, 4):
            with open(wig_files[index], "w") as fh:
                fh.write(example_wigs[index])
        self.multiparser.combine_wig(self.ref_folder, tmp_tar, "fasta")
        self.assertTrue(os.path.exists(os.path.join(tmp_tar, "test_forward.wig")))
        self.assertTrue(os.path.exists(os.path.join(tmp_tar, "test_reverse.wig")))

    def test_combine_gff(self):
        tmp_tar = os.path.join(self.tar_folder, "tmp")
        tmp_ref = os.path.join(self.ref_folder, "test.fa_folder")
        os.mkdir(tmp_ref)
        os.mkdir(tmp_tar)
        sub_fasta1 = os.path.join(tmp_ref, "aaa.fa")
        with open(sub_fasta1, "w") as rh:
            rh.write(self.example.sub_fasta1)
        sub_fasta2 = os.path.join(tmp_ref, "bbb.fa")
        with open(sub_fasta2, "w") as rh:
            rh.write(self.example.sub_fasta2)
        sub_gff1 = os.path.join(tmp_tar, "aaa.gff")
        with open(sub_gff1, "w") as rh:
            rh.write(self.example.sub_gff1)
        sub_gff2 = os.path.join(tmp_tar, "bbb.gff")
        with open(sub_gff2, "w") as rh:
            rh.write(self.example.sub_gff2)
        self.multiparser.combine_gff(self.ref_folder, tmp_tar, "fasta", None)
        self.assertTrue(os.path.exists(os.path.join(tmp_tar, "test.gff")))

    def test_parser_fasta(self):
        fasta_file = os.path.join(self.ref_folder, "test.fa")
        with open(fasta_file, "w") as rh:
            rh.write(self.example.fasta_file)
        self.multiparser.parser_fasta(self.ref_folder)
        self.assertTrue(os.path.exists(os.path.join(self.ref_folder, "tmp/aaa.fa")))
        self.assertTrue(os.path.exists(os.path.join(self.ref_folder, "tmp/bbb.fa")))
        self.assertTrue(os.path.exists(os.path.join(self.ref_folder, "test.fa_folder/aaa.fa")))
        self.assertTrue(os.path.exists(os.path.join(self.ref_folder, "test.fa_folder/bbb.fa")))

    def test_parser_gff(self):
        gff_file = os.path.join(self.ref_folder, "test.gff")
        with open(gff_file, "w") as rh:
            rh.write(self.example.gff_file)
        self.multiparser.parser_gff(self.ref_folder, None)
        self.assertTrue(os.path.exists(os.path.join(self.ref_folder, "tmp/aaa.gff")))
        self.assertTrue(os.path.exists(os.path.join(self.ref_folder, "tmp/bbb.gff")))
        self.assertTrue(os.path.exists(os.path.join(self.ref_folder, "test.gff_folder/aaa.gff")))
        self.assertTrue(os.path.exists(os.path.join(self.ref_folder, "test.gff_folder/bbb.gff")))
        tss_file = os.path.join(self.ref_folder, "test_TSS.gff")
        os.rename(gff_file, tss_file)
        tss_file = os.path.join(self.ref_folder, "test_TSS.gff")
        with open(tss_file, "w") as rh:
            rh.write(self.example.gff_file)
        self.multiparser.parser_gff(self.ref_folder, "TSS")
        self.assertTrue(os.path.exists(os.path.join(self.ref_folder, "tmp/aaa_TSS.gff")))
        self.assertTrue(os.path.exists(os.path.join(self.ref_folder, "tmp/bbb_TSS.gff")))
        self.assertTrue(os.path.exists(os.path.join(self.ref_folder, "test_TSS.gff_folder/aaa_TSS.gff")))
        self.assertTrue(os.path.exists(os.path.join(self.ref_folder, "test_TSS.gff_folder/bbb_TSS.gff")))

    def test_parser_wig(self):
        wig_f_file = os.path.join(self.ref_folder, "test_forward.wig")
        with open(wig_f_file, "w") as rh:
            rh.write(self.example.wig_f_file)
        wig_r_file = os.path.join(self.ref_folder, "test_reverse.wig")
        with open(wig_r_file, "w") as rh:
            rh.write(self.example.wig_r_file)
        self.multiparser.parser_wig(self.ref_folder)
        self.assertTrue(os.path.exists(os.path.join(self.ref_folder, "tmp/test_forward_STRAIN_aaa.wig")))
        self.assertTrue(os.path.exists(os.path.join(self.ref_folder, "tmp/test_forward_STRAIN_bbb.wig")))
        self.assertTrue(os.path.exists(os.path.join(self.ref_folder, "tmp/test_reverse_STRAIN_aaa.wig")))
        self.assertTrue(os.path.exists(os.path.join(self.ref_folder, "tmp/test_reverse_STRAIN_bbb.wig")))
        self.assertTrue(os.path.exists(
            os.path.join(self.ref_folder, "test_forward.wig_folder/test_forward_STRAIN_aaa.wig")))
        self.assertTrue(os.path.exists(
            os.path.join(self.ref_folder, "test_forward.wig_folder/test_forward_STRAIN_bbb.wig")))
        self.assertTrue(os.path.exists(
            os.path.join(self.ref_folder, "test_reverse.wig_folder/test_reverse_STRAIN_aaa.wig")))
        self.assertTrue(os.path.exists(
            os.path.join(self.ref_folder, "test_reverse.wig_folder/test_reverse_STRAIN_bbb.wig")))

class Example(object):

    fasta_file = """>aaa
CGCAGGTTGAGTTCCTGTTCCCGATAGATCCGATAAACCCGCTTATGATTCCAGAGCTGTCCCTGCACAT
>bbb
CGACAGGCCAGTGCTGACCCCATGATGCGCCACAGCTTGTGCGGCCAGTTCCCGGCGCTGGGCTG"""

    sub_fasta1 = """>aaa
CGCAGGTTGAGTTCCTGTTCCCGATAGATCCGATAAACCCGCTTATGATTCCAGAGCTGTCCCTGCACAT"""

    sub_fasta2 = """>bbb
CGACAGGCCAGTGCTGACCCCATGATGCGCCACAGCTTGTGCGGCCAGTTCCCGGCGCTGGGCTG"""

    gff_file = """##gff3
aaa	Refseq	gene	517	1878	.	+	.	db_xref=GeneID:3919798;locus_tag=SAOUHSC_00001;gene=dnaA
aaa	Refseq	gene	1234	2344	.	+	.	db_xref=GeneID:3919798;locus_tag=SAOUHSC_00001;gene=dnaA
bbb	Refseq	gene	5755	8456	.	-	.	db_xref=GeneID:3919180;locus_tag=SAOUHSC_00007
bbb	Refseq	gene	9755	10456	.	-	.	db_xref=GeneID:3919180;locus_tag=SAOUHSC_00007"""

    sub_gff1 = """aaa	Refseq	gene	517	1878	.	+	.	db_xref=GeneID:3919798;locus_tag=SAOUHSC_00001;gene=dnaA
aaa	Refseq	gene	1234	2344	.	+	.	db_xref=GeneID:3919798;locus_tag=SAOUHSC_00001;gene=dnaA"""

    sub_gff2 = """bbb	Refseq	gene	5755	8456	.	-	.	db_xref=GeneID:3919180;locus_tag=SAOUHSC_00007
bbb	Refseq	gene	9755	10456	.	-	.	db_xref=GeneID:3919180;locus_tag=SAOUHSC_00007"""

    wig_f_file = """track type=wiggle_0 name="test_forward"
variableStep chrom=aaa span=1
312 1.4041251228308191
313 56.867067474648174
314 56.867067474648174
315 56.867067474648174
variableStep chrom=bbb span=1
32 1.4041251228308191
33 56.867067474648174
34 56.867067474648174
35 56.867067474648174"""

    wig_r_file = """track type=wiggle_0 name="test_reverse"
variableStep chrom=aaa span=1
312 -1.4041251228308191
313 -56.867067474648174
314 -56.867067474648174
315 -56.867067474648174
variableStep chrom=bbb span=1
32 -1.4041251228308191
33 -56.867067474648174
34 -56.867067474648174
35 -56.867067474648174"""

    sub_f_wig1 = """track type=wiggle_0 name="test_forward"
variableStep chrom=aaa span=1
312 1.4041251228308191
313 56.867067474648174
314 56.867067474648174
315 56.867067474648174
"""

    sub_f_wig2 = """track type=wiggle_0 name="test_forward"
variableStep chrom=bbb span=1
32 1.4041251228308191
33 56.867067474648174
34 56.867067474648174
35 56.867067474648174
"""

    sub_r_wig1 = """track type=wiggle_0 name="test_reverse"
variableStep chrom=aaa span=1
312 -1.4041251228308191
313 -56.867067474648174
314 -56.867067474648174
315 -56.867067474648174
"""

    sub_r_wig2 = """track type=wiggle_0 name="test_reverse"
variableStep chrom=bbb span=1
32 -1.4041251228308191
33 -56.867067474648174
34 -56.867067474648174
35 -56.867067474648174
"""

if __name__ == "__main__":
    unittest.main()
