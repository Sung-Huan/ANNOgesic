#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
from subprocess import call, Popen
from annogesiclib.helper import Helper
from annogesiclib.multiparser import Multiparser
from annogesiclib.optimize_TSSpredator import optimization


def optimize_tss(tsspredator_path, fasta, gff, wigs, manual, output_folder,
                 project_name, max_height, max_height_reduction, max_factor, 
                 max_factor_reduction, max_base_height, utr_length, libs, 
                 replicate_name, cluster, length, core, program, replicate, steps):
    Multiparser()._parser_wig(wigs)
    wig_path = os.path.join(wigs, "tmp")
    Helper().check_uni_attributes(gff)
    optimization(tsspredator_path, max_height, max_height_reduction, max_factor,
                 max_factor_reduction, max_base_height, output_folder, core,
                 wig_path, project_name, fasta, replicate_name, steps, gff,
                 program, manual, libs, length)
