#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
from subprocess import call, Popen
from transaplib.helper import Helper
from transaplib.multiparser import Multiparser
from transaplib.optimize_TSSpredator import Optimization


def Optimize_TSS(tsspredator_path, fasta, gff, wigs, manual, output_folder,
                 project_name, max_height, max_height_reduction, max_factor, 
                 max_factor_reduction, max_base_height, utr_length, libs, 
                 replicate_name, cluster, length, core, program, replicate, steps):
    Multiparser()._parser_wig(wigs)
    wig_path = wigs + "/tmp/"
    Helper().check_uni_attributes(gff)
    Optimization(tsspredator_path, max_height, max_height_reduction, max_factor,
                 max_factor_reduction, max_base_height, output_folder, core,
                 wig_path, project_name, fasta, replicate_name, steps, gff,
                 program, manual, libs, length)
