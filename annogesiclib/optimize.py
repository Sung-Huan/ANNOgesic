import os
from annogesiclib.helper import Helper
from annogesiclib.multiparser import Multiparser
from annogesiclib.optimize_TSSpredator import optimization


def optimize_tss(tsspredator_path, fastas, gffs, wigs, manual, output_folder,
                 project_name, max_height, max_height_reduction, max_factor,
                 max_factor_reduction, max_base_height, max_enrichment,
                 max_processing, utr_length, libs, replicate_name, cluster,
                 length, core, program, replicate, steps):
    Multiparser().parser_wig(wigs)
    Multiparser().parser_gff(gffs, None)
    Multiparser().parser_fasta(fastas)
    gff_path = os.path.join(gffs, "tmp")
    wig_path = os.path.join(wigs, "tmp")
    fasta_path = os.path.join(fastas, "tmp")
    for gff in os.listdir(gff_path):
        if project_name in gff:
            gff_file = os.path.join(gff_path, gff)
            break
    for fa in os.listdir(fasta_path):
        if project_name in fa:
            fasta_file = os.path.join(fasta_path, fa)
            break
    Helper().check_uni_attributes(gff_file)
    optimization(tsspredator_path, max_height, max_height_reduction,
                 max_factor, max_factor_reduction, max_base_height,
                 max_enrichment, max_processing, output_folder, core,
                 wig_path, project_name, fasta_file, replicate_name, steps,
                 gff_file, program, manual, libs, length, cluster,
                 utr_length, replicate)
    Helper().remove_all_content(os.path.join(output_folder, "optimized_TSSpredator"), "config", "file")
    Helper().remove_all_content(os.path.join(output_folder, "optimized_TSSpredator"), "Master", "dir")
