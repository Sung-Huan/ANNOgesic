import os
import sys
from annogesiclib.helper import Helper
from annogesiclib.multiparser import Multiparser
from annogesiclib.optimize_TSSpredator import optimization


def optimize_tss(args_ops):
    if len(os.listdir(args_ops.gffs)) == 0:
        print("Error: there is no gff files!!!")
        sys.exit()
    if len(os.listdir(args_ops.fastas)) == 0:
        print("Error: there is no fasta files!!!")
        sys.exit()
    if len(os.listdir(args_ops.wigs)) == 0:
        print("Error: there is no wiggle files!!!")
        sys.exit()
    Multiparser().parser_wig(args_ops.wigs)
    Multiparser().parser_gff(args_ops.gffs, None)
    Multiparser().parser_fasta(args_ops.fastas)
    gff_path = os.path.join(args_ops.gffs, "tmp")
    wig_path = os.path.join(args_ops.wigs, "tmp")
    fasta_path = os.path.join(args_ops.fastas, "tmp")
    for gff in os.listdir(gff_path):
        if args_ops.project_strain in gff:
            gff_file = os.path.join(gff_path, gff)
            break
    for fa in os.listdir(fasta_path):
        if args_ops.project_strain in fa:
            fasta_file = os.path.join(fasta_path, fa)
            break
    Helper().check_uni_attributes(gff_file)
    optimization(wig_path, fasta_file, gff_file, args_ops)
    Helper().remove_all_content(os.path.join(args_ops.output_folder,
                                "optimized_TSSpredator"), "config", "file")
    Helper().remove_all_content(os.path.join(args_ops.output_folder,
                                "optimized_TSSpredator"), "Master", "dir")
    Helper().remove_tmp(args_ops.wigs)
    Helper().remove_tmp(args_ops.gffs)
    Helper().remove_tmp(args_ops.fastas)
