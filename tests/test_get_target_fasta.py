#!/usr/bin/python

import os
import shutil
import csv
import unittest
import sys
sys.path.append(".")
from io import StringIO
from annogesiclib.get_target_fasta import TargetFasta


class Mock_multiparser(object):

    def parser_fasta(folder):
        tmp_folder = os.path.join(folder, "tmp")
        os.mkdir(tmp_folder)
        fasta_file1 = """>aaa
CGCAGGTTGAGTTCCTGTTCCCGATAGATCCGATAAACCCGCTTATGATTCCAGAGCTGTCCCTGCACAT
TGCGCAGATACAGGAAACACAGACCAAATCCCCATCTCCTGTGAGCCTGGGTCAGTCCCACCAGAAGAGC
GGCAATCCTGTCGTTCTCCGCTGCCAGTCGCGGACGATAGCGAAAGCAGGTCTCGGATATCCCAAAAATC
CGACAGGCCAGCGCAATGCTGACCCCATGATGCGCCACAGCTTGTGCGGCCAGTTCCCGGCGCTGGGCTG
GCCGCTTCATTTTTTTCCAAGGGCTTCCTTCAGGATATCCGTCTGCATGCTCAAATCCGCATACATGCGC"""
        fasta_file2 = """>bbb
CGCAGGTTGAGTTCCTGTTCCCGATAGATCCGATAAACCCGCTTATGATTCCAGAGCTGTCCCTGCACAT
TGCGCAGATACAGGAAACACAGACCAAATCCCCATCTCCTGTGAGCCTGGGTCAGTCCCACCAGAAGAGC
GGCAATCCTGTCGTTCTCCGCTGCCAGTCGCGGACGATAGCGAAAGCAGGTCTCGGATATCCCAAAAATC
CGACAGGCCAGCGCAATGCTGACCCCATGATGCGCCACAGCTTGTGCGGCCAGTTCCCGGCGCTGGGCTG
GCCGCTTCATTTTTTTCCAAGGGCTTCCTTCAGGATATCCGTCTGCATGCTCAAATCCGCATACATGCGC"""
        f_h = open(os.path.join(tmp_folder, "aaa.fa"), "w")
        f_h.write(fasta_file1)
        f_h.close()
        f_h = open(os.path.join(tmp_folder, "bbb.fa"), "w")
        f_h.write(fasta_file2)
        f_h.close()

class Mock_seq_editer(object):

    def modify_seq(ref, mut_table, tar):
        shutil.copyfile("a_test_project/ref/tmp/aaa.fa", os.path.join(tar, "aaa.fa"))
        shutil.copyfile("a_test_project/ref/tmp/bbb.fa", os.path.join(tar, "bbb.fa"))

class Mock_helper(object):

    def remove_all_content(file_, type_, folder):
        pass
        

class TestTargetFasta(unittest.TestCase):

    def setUp(self):
        self.root_folder = "a_test_project"
        if os.path.exists(self.root_folder):
            shutil.rmtree(self.root_folder)
        self.ref_folder = os.path.join(self.root_folder, "ref")
        self.tar_folder = os.path.join(self.root_folder, "tar")
        os.mkdir(self.root_folder)
        os.mkdir(self.ref_folder)
        os.mkdir(self.tar_folder)
        self.target_fasta = TargetFasta(self.tar_folder, self.ref_folder)
        self.target_fasta.seq_editer = Mock_seq_editer
        self.target_fasta.multiparser = Mock_multiparser
        self.target_fasta.helper = Mock_helper
        self.example = ExampleData()

    def tearDown(self):
         if os.path.exists(self.root_folder):
            shutil.rmtree(self.root_folder)

    def test_get_target_fasta(self):
        self.target_fasta.get_target_fasta(None,
                                           self.tar_folder,
                                           self.ref_folder,
                                           ["ccc:aaa,bbb"])
        self.assertTrue(os.path.exists(os.path.join(self.tar_folder, "ccc.fa")))

class ExampleData(object):

    fasta_file1 = """>aaa
CGCAGGTTGAGTTCCTGTTCCCGATAGATCCGATAAACCCGCTTATGATTCCAGAGCTGTCCCTGCACAT
TGCGCAGATACAGGAAACACAGACCAAATCCCCATCTCCTGTGAGCCTGGGTCAGTCCCACCAGAAGAGC
GGCAATCCTGTCGTTCTCCGCTGCCAGTCGCGGACGATAGCGAAAGCAGGTCTCGGATATCCCAAAAATC
CGACAGGCCAGCGCAATGCTGACCCCATGATGCGCCACAGCTTGTGCGGCCAGTTCCCGGCGCTGGGCTG
GCCGCTTCATTTTTTTCCAAGGGCTTCCTTCAGGATATCCGTCTGCATGCTCAAATCCGCATACATGCGC"""

    fasta_file2 = """>bbb
CGCAGGTTGAGTTCCTGTTCCCGATAGATCCGATAAACCCGCTTATGATTCCAGAGCTGTCCCTGCACAT
TGCGCAGATACAGGAAACACAGACCAAATCCCCATCTCCTGTGAGCCTGGGTCAGTCCCACCAGAAGAGC
GGCAATCCTGTCGTTCTCCGCTGCCAGTCGCGGACGATAGCGAAAGCAGGTCTCGGATATCCCAAAAATC
CGACAGGCCAGCGCAATGCTGACCCCATGATGCGCCACAGCTTGTGCGGCCAGTTCCCGGCGCTGGGCTG
GCCGCTTCATTTTTTTCCAAGGGCTTCCTTCAGGATATCCGTCTGCATGCTCAAATCCGCATACATGCGC"""

    mutation_table = """#target_id	reference_id	ref_nt	pos	tar_nt	relationship	locus_tag	gene	description
bbb	aaa	g	2	-	Frame shift	AAAAA_00001	dnaA	Hypothetical protein
bbb	aaa	a	4	t		AAAAA_00004	dnaR	
bbb	aaa	c	14	a				
bbb	aaa	-	15	g		AAAAA_00023		Test protein
bbb	aaa	-	23	a						
bbb	aaa	-	25	t					
bbb	aaa	t	43	-					
"""

if __name__ == "__main__":
    unittest.main()
