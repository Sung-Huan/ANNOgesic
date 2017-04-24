#!/usr/bin/python

import os
import shutil
import csv
import unittest
import sys
sys.path.append(".")
from io import StringIO
from annogesiclib.get_target_fasta import TargetFasta
from mock_helper import gen_file


class Mock_multiparser(object):

    def parser_fasta(folder):
        if not isinstance(folder, str):
            folder = "a_test_project/ref"
        tmp_folder = os.path.join(folder ,"tmp")
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
        Mock_multiparser().parser_fasta()
        shutil.copyfile("a_test_project/ref/tmp/aaa.fa",
                        os.path.join(tar, "aaa.fa"))
        shutil.copyfile("a_test_project/ref/tmp/bbb.fa",
                        os.path.join(tar, "bbb.fa"))

class Mock_helper(object):

    def remove_all_content(file_, type_, folder):
        pass

    def check_make_folder(folder):
        path = "/".join(folder.split("/")[:-1])
        folder = folder.split("/")[-1]
        if folder in os.listdir(path):
            shutil.rmtree(os.path.join(path, folder))
        os.mkdir(os.path.join(path, folder))
        

class TestTargetFasta(unittest.TestCase):

    def setUp(self):
        self.root_folder = "a_test_project"
        if os.path.exists(self.root_folder):
            shutil.rmtree(self.root_folder)
        self.ref_folder = os.path.join(self.root_folder, "ref")
        self.tar_folder = os.path.join(self.root_folder, "tar")
        self.fasta_folder = os.path.join(self.root_folder, "fasta_files")
        os.mkdir(self.root_folder)
        os.mkdir(self.ref_folder)
        os.mkdir(self.tar_folder)
        os.mkdir(os.path.join(self.root_folder, "fasta_files"))
        self.target_fasta = TargetFasta(self.tar_folder, self.ref_folder)
        self.target_fasta.seq_editer = Mock_seq_editer
        self.target_fasta.multiparser = Mock_multiparser
        self.target_fasta.helper = Mock_helper
        self.example = ExampleData()

    def tearDown(self):
         if os.path.exists(self.root_folder):
            shutil.rmtree(self.root_folder)

    def test_get_target_fasta(self):
        gen_file("a_test_project/ref/bbb.fa", self.example.fasta_file2)
        gen_file("a_test_project/mut", self.example.mutation_table)
        self.target_fasta.get_target_fasta("a_test_project/mut",
                                           self.tar_folder,
                                           ["a_test_project/ref/bbb.fa"],
                                           False,
                                           self.root_folder)
        self.assertTrue(os.path.exists(
            os.path.join(self.fasta_folder, "aaa.fa")))

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
