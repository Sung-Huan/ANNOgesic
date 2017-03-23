import os
import sys
from annogesiclib.helper import Helper
from annogesiclib.multiparser import Multiparser
from annogesiclib.optimize_TSSpredator import optimization


def get_length(fasta_file):
    length = 0
    with open(fasta_file) as fh:
        for line in fh:
            line = line.strip()
            if not line.startswith(">"):
                length = length + len(line)
    return length 


def optimize_tss(args_ops):
    if len(os.listdir(args_ops.gffs)) == 0:
        print("Error: There is no gff file!")
        sys.exit()
    if len(os.listdir(args_ops.fastas)) == 0:
        print("Error: There is no fasta file!")
        sys.exit()
    if len(os.listdir(args_ops.wigs)) == 0:
        print("Error: There is no wiggle file!")
        sys.exit()
    Multiparser().parser_wig(args_ops.wigs)
    Multiparser().parser_gff(args_ops.gffs, None)
    Multiparser().parser_fasta(args_ops.fastas)
    Multiparser().parser_gff(args_ops.manuals, None)
    gff_path = os.path.join(args_ops.gffs, "tmp")
    wig_path = os.path.join(args_ops.wigs, "tmp")
    fasta_path = os.path.join(args_ops.fastas, "tmp")
    manual_path = os.path.join(args_ops.manuals, "tmp")
    if "all" not in args_ops.strain_lengths.keys():
        for strain in args_ops.strain_lengths.keys():
            detect = False
            for man in os.listdir(manual_path):
                if strain == man.replace(".gff", ""):
                    detect = True
            if not detect:
                print("Error: There are genomes in --genome_lengths "
                      "which is not contained in manual-detected "
                      "TSS gff files!")
                sys.exit()
    for man in os.listdir(manual_path):
        run = False
        prefix = man.replace(".gff", "")
        man_file = os.path.join(manual_path, man)
        if (prefix in args_ops.strain_lengths.keys()):
            length = args_ops.strain_lengths[prefix]
            run = True
        elif("all" in args_ops.strain_lengths.keys()):
            length = "all"
            run = True
        if run:
            for gff in os.listdir(gff_path):
                if (gff[:-4] == prefix) and (".gff" in gff):
                    gff_file = os.path.join(gff_path, gff)
                    break
            for fa in os.listdir(fasta_path):
                if (".".join(fa.split(".")[:-1]) == prefix) and (
                        ".fa" in fa):
                    fasta_file = os.path.join(fasta_path, fa)
                    break
            if length == "all":
                length = get_length(fasta_file)
            Helper().check_uni_attributes(gff_file)
            optimization(wig_path, fasta_file, gff_file, args_ops,
                         man_file, length, prefix)
            Helper().remove_all_content(os.path.join(
                args_ops.output_folder,
                "optimized_TSSpredator"), "config", "file")
            Helper().remove_all_content(os.path.join(
                args_ops.output_folder,
                "optimized_TSSpredator"), "Master", "dir")
    Helper().remove_tmp_dir(args_ops.wigs)
    Helper().remove_tmp_dir(args_ops.gffs)
    Helper().remove_tmp_dir(args_ops.fastas)
    Helper().remove_tmp_dir(args_ops.manuals)
