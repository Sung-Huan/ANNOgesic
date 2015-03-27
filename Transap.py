#!/usr/bin/python

"""transcript annotation and analysis pipeline"""
import argparse
import os
from transaplib.controller import Controller

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"
__version__ = "0.1.0"

def main():
    home_path = os.environ["HOME"]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--version", "-v", default=False, action="store_true",
        help="show version")
    parser.add_argument(
        "--bin_path", "-B", default="/home/silas/Transap/bin",
        help="path of scripts")
    subparsers = parser.add_subparsers(help="commands")
    # Arguments for project creation
    create_project_parser = subparsers.add_parser(
        "create", help="Create a project")
    create_project_parser.add_argument(
        "project_path", default=".", help="Name/path of the project.")
    create_project_parser.set_defaults(func=create_project)
    # Parameters for get input files
    get_input_parser = subparsers.add_parser(
        "get_input_files", help="Get required files (i.e. annotation files, fasta files)")
    get_input_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    get_input_parser.add_argument(
        "--FTP_path", "-F",
        help="Path of website where can download the required files")
    get_input_parser.add_argument(
        "--ref_fasta", "-f", default=False, action="store_true",
        help="Download fasta files of reference.")
    get_input_parser.add_argument(
        "--ref_gff", "-g", default=False, action="store_true",
        help="Download gff files of reference.")
    get_input_parser.add_argument(
        "--ref_gbk", "-k", default=False, action="store_true",
        help="Download genbank files of reference.")
    get_input_parser.add_argument(
        "--convert_embl", "-e", default=False, action="store_true",
        help="convert gbk to embl files of reference.")
    get_input_parser.set_defaults(func=get_input)
    # get target fasta
    get_target_fasta_parser = subparsers.add_parser(
        "get_target_fasta", help="Get target fasta.")
    get_target_fasta_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    get_target_fasta_parser.add_argument(
        "--ref_fasta_folder", "-r", default="Transap/input/reference/target/fasta",
        help="Folder of fasta files.")
    get_target_fasta_parser.add_argument(
        "--mutation_table", "-m", default="Transap/input/mutation_table",
        help="The path of mutation table.")
    get_target_fasta_parser.add_argument(
        "--combine_all_fasta", "-c", default=False,
        help="Conbime all fasta files to one fasta file.")
    get_target_fasta_parser.add_argument(
        "--output_format", "-o", default=False, nargs="+",
        help="Please assign the output filename and which strain should be included in it. "
        "For example: FILE1:strain1,strain2. FILE1 is a output fasta file which include the information of strain1 and strain2.")
    get_target_fasta_parser.set_defaults(func=get_target_fasta)    

    # run RATT
    RATT_parser = subparsers.add_parser(
        "annotation_transfer", help="Run RATT to transfer the annotation files from reference to target")
    RATT_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    RATT_parser.add_argument(
        "--PAGIT_folder", default=home_path + "/PAGIT",
        help="Path of the PAGIT folder. If your PAGIT folder is not located at " + home_path + "/PAGIT, please assign here.")
    RATT_parser.add_argument(
        "--RATT_path", default=home_path + "/PAGIT/RATT",
        help="Path of the RATT folder. If your RATT folder is not located at " + home_path + "/PAGIT/RATT, please assign here.")
    RATT_parser.add_argument(
        "--target_ref_comparison", "-tr", nargs="+",
        help="Please assign the fasta file of terget one and reference one. Use ':' to separate target and reference. "
        "For example, target1:ref1 target2:ref2 . "
        "ATTENTION:please make sure the target and ref names are the same as the header of target and reference fasta files. ")
    RATT_parser.add_argument(
        "--element", "-e",
        help="element for running RATT.")
    RATT_parser.add_argument(
        "--transfer_type", "-t",
        help="the transfer type for running RATT.(details can refer to the manual of RATT.)")
    RATT_parser.add_argument(
        "--ref_embl", "-re", default="Transap/input/reference/annotation",
        help="the folder which stores every reference embl folders.")
    RATT_parser.add_argument(
        "--ref_fasta", "-rf", default="Transap/input/reference/fasta",
        help="the folder which stores reference fasta files.")
    RATT_parser.add_argument(
        "--target_fasta", "-tf", default="Transap/output/target/fasta",
        help="the folder which stores target fasta files.")
    RATT_parser.add_argument(
        "--convert_to_gff_rnt_ptt", "-g", default=False, action="store_true",
        help="Do you want to convert to gff, rnt and ptt ?")
    RATT_parser.add_argument(
        "--combine_all_gff", "-c", default=False,
        help="Combine all gff files to be one gff file. please keyin the file name")
    RATT_parser.set_defaults(func=run_RATT)
    # Parameters of TSSpredator
    TSSpredator_parser = subparsers.add_parser(
        "tsspredator", help="Run TSSpredator to predict TSSs.")
    TSSpredator_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    TSSpredator_parser.add_argument(
        "--TSSpredator_path", default="TSSpredator.jar",
        help="If you want to assign the path of TSSpredator, please assign here.")
    TSSpredator_parser.add_argument(
        "--fasta_folder", "-f", default="Transap/output/target/fasta",
        help="Path of the target fasta folder.")
    TSSpredator_parser.add_argument(
        "--annotation_folder", "-g", default="Transap/output/target/annotation",
        help="Path of the target gff folder.")
    TSSpredator_parser.add_argument(
        "--wig_folder", "-w", default="Transap/input/wigs/tex_notex",
        help="The folder of the wig folder.")
    TSSpredator_parser.add_argument(
        "--height", "-he",
        help="This value relates to the minimal number of read starts at a "
	"certain genomic position to be considered as a TSS candidate.")
    TSSpredator_parser.add_argument(
        "--height_reduction", "-rh",
        help="When comparing different strains/conditions and "
	"the step height threshold is reached in at least one strain/condition, "
	"the threshold is reduced for the other strains/conditions "
	"by the value set here. "
        "This value must be smaller than the step height threshold.")
    TSSpredator_parser.add_argument(
        "--factor", "-fa",
        help="This is the minimal factor by which the TSS height has to "
	"exceed the local expression background.")
    TSSpredator_parser.add_argument(
        "--factor_reduction", "-rf",
        help="When comparing different strains/conditions and "
        "the step factor threshold is reached in at least one strain/condition, "
        "the threshold is reduced for the other strains/conditions "
        "by the value set here. "
        "This value must be smaller than the step factor threshold.")
    TSSpredator_parser.add_argument(
        "--base_height", "-bh",
        help="This is the minimal number of reads should be mapped on TSS. ")
    TSSpredator_parser.add_argument(
        "--replicate_match", "-rm",
        help="The TSS candidates should match to how many number of the replicates")
    TSSpredator_parser.add_argument(
        "--utr_length", "-u", default=300, type=int,
        help="The length of UTR. It is for Primary and Secondary definition.")
    TSSpredator_parser.add_argument(
        "--lib", "-l", nargs='+',
        help="The libraries of wig files for TSSpredator. The format is: "
	"wig_file_name:tex_treat_or_not(tex or notex):condition_id(integer):replicate_id(alphabet):strand(+ or -).")
    TSSpredator_parser.add_argument(
        "--output_prefix", "-p", nargs='+',
        help="The output prefix of all conditions.")
    TSSpredator_parser.add_argument(
        "--merge_manual", "-m", default=False,
        help="If you have gff file of manual checked TSS, you can use this function to merge manual ones and predicted ones.")
    TSSpredator_parser.add_argument(
        "--statistics", "-s", default=False, action="store_true",
        help="Doing statistics for TSS candidates. "
        "it will store in statistics folder.")
    TSSpredator_parser.add_argument(
        "--validate_gene", "-v", default=False, action="store_true",
        help="Using TSS candidates to validate genes in annotation file. "
        "it will store in statistics folder."
        "It also will output a gff file of annotation informatiom which only start from TSSs. "
        "The gff file will store in the folder - annotation_with_TSS")
    TSSpredator_parser.add_argument(
        "--compute_program", "-t", default="TSS",
        help="Which program do you want to predict. (TSS or processing_site)")
    TSSpredator_parser.add_argument(
        "--compare_transcript_assembly", "-ta", default=False,
        help="compare with transcriptome assembly")
    TSSpredator_parser.add_argument(
        "--fuzzy", "-fu", default=5, type=int,
        help="the fuzzy for comparing TSS and transcript assembly.")
    TSSpredator_parser.add_argument(
        "--cluster", "-c", default=2, type=int,
        help="This number if for compare manual detected TSS and prediced one. "
        "If the position between manual one and predicted one is smaller or equal than this value, "
        "it will only print one of them.")
    TSSpredator_parser.add_argument(
        "--length", "-le", default=False,
        help="The length that you want to compare with manual check for statistics. ")
    TSSpredator_parser.set_defaults(func=run_TSSpredator)
    # Parameter of opimization of TSSpredator
    op_TSSpredator_parser = subparsers.add_parser(
        "optimize_tsspredator", help="Optimize TSSpredator based on (partial)manual detect one.")
    op_TSSpredator_parser.add_argument(
        "--TSSpredator_path", default="TSSpredator.jar",
        help="If you want to assign the path of TSSpredator, please assign here.")
    op_TSSpredator_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    op_TSSpredator_parser.add_argument(
        "--fasta_file", "-fs",
        help="Path of the target fasta file.")
    op_TSSpredator_parser.add_argument(
        "--annotation_file", "-g",
        help="Path of the target gff file.")
    op_TSSpredator_parser.add_argument(
        "--wig_folder", "-w", default="Transap/input/wigs/tex_notex",
        help="The folder of the wig folder.")
    op_TSSpredator_parser.add_argument(
        "--manual", "-m",
        help="The file of manual check gff file.")
    op_TSSpredator_parser.add_argument(
        "--strain_name", "-n",
        help="The name of the strain you want to optimize.")
    op_TSSpredator_parser.add_argument(
        "--max_height", "-he", default=2.5, type=float,
        help="This value relates to the minimal number of read starts at a "
        "certain genomic position to be considered as a TSS candidate."
        "During optimization will never larger than this value.")
    op_TSSpredator_parser.add_argument(
        "--max_height_reduction", "-rh", default=2.4, type=float,
        help="When comparing different strains/conditions and "
        "the step height threshold is reached in at least one strain/condition, "
        "the threshold is reduced for the other strains/conditions "
        "by the value set here. "
        "This value must be smaller than the step height threshold."
        "During optimization will never larger than this value.")
    op_TSSpredator_parser.add_argument(
        "--max_factor", "-fa", default=10, type=float,
        help="This is the minimal factor by which the TSS height has to "
        "exceed the local expression background."
        "During optimization will never larger than this value.")
    op_TSSpredator_parser.add_argument(
        "--max_factor_reduction", "-rf", default=9.9, type=float,
        help="When comparing different strains/conditions and "
        "the step factor threshold is reached in at least one strain/condition, "
        "the threshold is reduced for the other strains/conditions "
        "by the value set here. "
        "This value must be smaller than the step factor threshold."
        "During optimization will never larger than this value.")
    op_TSSpredator_parser.add_argument(
        "--max_base_height", "-bh", default=0.2, type=float,
        help="This is the minimal number of reads should be mapped on TSS. "
        "During optimization will never larger than this value.")
    op_TSSpredator_parser.add_argument(
        "--utr_length", "-u", default=300, type=int,
        help="The length of UTR. It is for Primary and Secondary definition.")
    op_TSSpredator_parser.add_argument(
        "--lib", "-l", nargs='+',
        help="The libraries of wig files for TSSpredator. The format is: "
        "wig_file_name:tex_treat_or_not(tex or notex):condition_id(integer):replicate_id(alphabet):strand(+ or -).")
    op_TSSpredator_parser.add_argument(
        "--output_prefix", "-p", nargs='+',
        help="The output prefix of all conditions.")
    op_TSSpredator_parser.add_argument(
        "--cluster", "-cu", default=2, type=int,
        help="This number if for compare manual detected TSS and prediced one. "
        "If the position between manual one and predicted one is smaller or equal than this value, "
        "it will only print one of them.")
    op_TSSpredator_parser.add_argument(
        "--length", "-le", default=False,
        help="The length of nts for running optimization."
        "Default is compare whole genome")
    op_TSSpredator_parser.add_argument(
        "--core", "-c", type=int,
        help="How many paralle running do you want to use.")
    op_TSSpredator_parser.add_argument(
        "--program", "-r", default="TSS",
        help="The type which you want to run TSSpredator (TSS or Processing_site).")
    op_TSSpredator_parser.add_argument(
        "--replicate_match", "-rm",
        help="The TSS candidates should match to how many number of the replicates")
    op_TSSpredator_parser.add_argument(
        "--steps", "-s", default=6000, type=int,
        help="how many steps do you want to run.")
    op_TSSpredator_parser.set_defaults(func=optimize_TSSpredator)
    # Parameter of generating color png
    color_parser = subparsers.add_parser(
        "color_png", help="Generating color screenshots of TSS or processing site. "
        "It only works after running batch script. "
        "If the color bands are not located at the proper positions, "
        "please increase --figure_height. ")
    color_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    color_parser.add_argument(
        "--screenshot_type", "-s",
        help="please assign the screenshots is \"TSS\" or \"Processing_site\". ")
    color_parser.add_argument(
        "--track_number", "-t", type=int,
        help="how many number of tracks do you have. ")
    color_parser.set_defaults(func=color_png)
    # Parameter of TransTermHP
    TransTerm_parser = subparsers.add_parser(
        "terminator", help="run TransTermHP for detect Terminators.")
    TransTerm_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    TransTerm_parser.add_argument(
        "--TransTermHP_folder",
        help="Please assign the path of the folder of TransTermHP.")
    TransTerm_parser.add_argument(
        "--RNAfold_path", default="RNAfold",
        help="If you want to assign the path of RNAfold, please assign here.")
    TransTerm_parser.add_argument(
        "--fasta_folder", "-f", default="Transap/output/target/fasta",
        help="the path of fasta folder.")
    TransTerm_parser.add_argument(
        "--annotation_folder", "-g", default="Transap/output/target/annotation",
        help="the path of annotation gff folder.")
    TransTerm_parser.add_argument(
        "--transcript_folder", "-a", default="Transap/output/transcriptome_assembly/gffs",
        help="the path of transcript gff folder.")
    TransTerm_parser.add_argument(
        "--sRNA", "-sr", default=False,
        help="If you want to include sRNA information, please assign the gff file of sRNA.")
    TransTerm_parser.add_argument(
        "--statistics", "-s", default=False, action="store_true",
        help="Doing statistics for TransTermHP. "
        "the name of statistics file is - stat_terminator_$STRAIN_NAME.csv.")
    TransTerm_parser.add_argument(
        "--tex_wig_folder", "-tw", default=False,
        help="If you want to use tex +/- libraries, please assign fragmented wig folder.")
    TransTerm_parser.add_argument(
        "--frag_wig_folder", "-fw", default=False,
        help="If you want to use fragmented libraries, please assign fragmented wig folder.")
    TransTerm_parser.add_argument(
        "--decrease", "-d", default=0.5, type=float,
        help="If the (lowest coverage / highest coverage) in the terminator is smaller than this number, "
        "it will consider this terminator have dramatic coverage decrease in it. (this only for statistics)")
    TransTerm_parser.add_argument(
        "--fuzzy", "-u", default=10, type=int,
        help="It will elongate the number of nt(you assign here) from both terminal site. "
        "If it can found the coverage dramatic decrease, it will consider the terminator have dramatic coverage decrease in it.")
    TransTerm_parser.add_argument(
        "--highest_coverage", "-hc", default=0, type=float,
        help="If the highest coverage of the region of terminator is below to this number, "
        "the terminator will be classify to non-detect.")
    TransTerm_parser.add_argument(
        "-tl","--tex_notex_libs", default=False, nargs="+", help="library name of tex and notex library.")
    TransTerm_parser.add_argument(
        "-fl","--frag_libs", default=False, nargs="+", help="library name of fragment library.")
    TransTerm_parser.add_argument(
        "-te","--tex_notex", default=2, type=int,
        help="For tex +/- library, terminators should be detected by both or just one.(1/2)")
    TransTerm_parser.add_argument(
        "-r","--replicates", type=int, help="how many replicates should detect the same terminator.")
    TransTerm_parser.add_argument(
        "-tb","--table_best", default=False, action="store_true",
        help="output sRNA table only most decreasing track")
    TransTerm_parser.set_defaults(func=run_TransTermHP)

    # Parameter of Transcript assembly
    Transcript_parser = subparsers.add_parser(
        "transcript_assembly", help="run Transcript for doing transcriptome assembly.")
    Transcript_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    Transcript_parser.add_argument(
        "--annotation_folder", "-g", default=False,
        help="It is for comparing transcript assembly and annotation gff file. "
        "It can use annotation gff file as reference and modify transcript assembly file. "
        "If you don't want to do it, you don't need to turn it on.")
    Transcript_parser.add_argument(
        "--sort_annotation", "-s", default=False, action="store_true",
        help="The annotation gff files in annotation folder are sorted or not. "
        "If they didn't be sorted, please turn it on. ")
    Transcript_parser.add_argument(
        "--length", "-l", default=20, type=int,
        help="the minimum width of region to be a transcript. It is for refer to annotation file. "
        "If you want to compare with annotation files, it will be the final output. "
        "The default is 20.")
    Transcript_parser.add_argument(
        "--normal_wig_path", "-nw", default=False,
        help="The path of normal wig folder.")
    Transcript_parser.add_argument(
        "--frag_wig_path", "-fw", default=False,
        help="The path of fragment wig folder.")
    Transcript_parser.add_argument(
        "--replicates", "-p", nargs='+', default=False,
        help="If there are replicates, please assign to it. "
        "i.e. <replicate a1>:<replicate a2>:<replicate a3> <replicate b1>:<replicate b2>:<replicate b3>. "
        "It will use average coverage to compute it. "
        "default is no replicates.")
    Transcript_parser.add_argument(
        "--height", "-he", default=5, type=int,
        help="the minimum height of coverage to be a transcript. "
        "The default is 5.")
    Transcript_parser.add_argument(
        "--width", "-w", default=20, type=int,
        help="the minimum width of region to be a transcript. It is for without annotation to be reference. "
        "If you don't want to compare with annotation files (--length), it will be the final output. "
        "The default is 20.")
    Transcript_parser.add_argument(
        "--tolerance", "-t", default=5, type=int,
        help="This number indicates how willing the algorithm is to ignore a temporary drop below either cut-off. "
        "The default is 5.")
    Transcript_parser.add_argument(
        "--supported", "-r", type=int,
        help="the position is included in the transcript if there are more than the replicate which you assign here to supported it.")
    Transcript_parser.add_argument(
        "--tex_notex", "-te", default=2, type=int,
        help="If you use tex +/- libraries to run transcript assembly, please assign the tex +/- should both consider or just one. (1 or 2)")
    Transcript_parser.add_argument(
        "--compare_TSS", "-ct", default=False,
        help="Compare with TSS. Please assign TSS folder.")
    Transcript_parser.add_argument(
        "--compare_CDS", "-cg", default=False,
        help="Compare with annotation gff file. Please assign annotation folder.")
    Transcript_parser.add_argument(
        "--TSS_fuzzy", "-fu", default=5, type=int,
        help="the fuzzy for comparing TSS and transcript assembly.")
    Transcript_parser.add_argument(
        "--Tex_treated_libs", "-tl", nargs="+", default=False,
        help="input of tex +/- library.")
    Transcript_parser.add_argument(
        "--fragmented_libs", "-fl", nargs="+", default=False,
        help="input of fragmented library.")
    
    Transcript_parser.set_defaults(func=run_Transcript_Assembly)

    # Parameter of UTR detection
    UTR_parser = subparsers.add_parser(
        "utr", help="run UTR detection to detecting 5'UTR and 3'UTR.")
    UTR_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    UTR_parser.add_argument(
        "--annotation_folder", "-g", default="Transap/output/target/annotation",
        help="The path of annotation gff folder.")
    UTR_parser.add_argument(
        "--TSS_folder", "-t", default="Transap/output/TSS/gffs",
        help="The path of TSS folder.")
    UTR_parser.add_argument(
        "--transcript_assembly_folder", "-a", default="Transap/output/transcriptome_assembly/gffs",
        help="The path of transcriptome assembly folder.")
    UTR_parser.add_argument(
        "--terminator_folder", "-e", default=False,
        help="If you want to add the information of terminator, you can assign the path of terminator folder here.")
    UTR_parser.add_argument(
        "--terminator_fuzzy", "-f", default=10, type=int,
        help="If the distance(nt) of between terminator and the end of transcript assembly belows to this value, "
        "it will assign the terminator associated with the 3'UTR.")
    UTR_parser.add_argument(
        "--TSS_source", "-s", default=True, action="store_false",
        help="If you generate TSS from other method, please turn it on.")
    UTR_parser.set_defaults(func=run_UTR_detection)

    # Parameter of sRNA detectopn
    sRNA_parser = subparsers.add_parser(
        "srna", help="run sRNA detection to detecting sRNA candidates.")
    sRNA_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    sRNA_parser.add_argument(
        "--UTR_derived_sRNA", "-u", default=False, action="store_true",
        help="If you want to detect UTR derived sRNA, please turn it on.")
    sRNA_parser.add_argument(
        "--import_info", "-d", nargs="+",
        help="There are several types of information you can import to detect sRNA: "
        "TSS(1), energy of secondary structure(2), blast to nr(3), blast to sRNA (4), sORF(5). "
        "without any information, only detect sRNA (without any information) by transcriptome assembly(6)"
        "Please assign the number of type you want to import, i.e. 1 2 4 - means it used TSS, energy and blast result to detect sRNA. "
        "Besides these information, it will also consider the sequence length of sRNA.")
    sRNA_parser.add_argument(
        "--transcript_assembly_folder", "-a", default="Transap/output/transcriptome_assembly/gffs",
        help="The path of transcriptome assembly folder.")
    sRNA_parser.add_argument(
        "--annotation_folder", "-g", default="Transap/output/target/annotation",
        help="The path of annotation gff folder.")
    sRNA_parser.add_argument(
        "--TSS_folder", "-t", default=False,
        help="If you want to import TSS information, please assign the path of gff folder of TSS.")
    sRNA_parser.add_argument(
        "--processing_site_folder", "-p", default=False,
        help="If you want to import processing site information (only for UTR derived sRNA), please assign the path of gff folder of processing site.")
    sRNA_parser.add_argument(
        "--TSS_fuzzy", "-ft", default=2, type=int,
        help="If you want to import TSS information, the fuzzy for comparing TSS and transcript assembly. "
        "The default number is 2.")
    sRNA_parser.add_argument(
        "--min_length", "-lm",default=30,
        help="Please assign the minium length of sRNA. "
        "It will classify sRNA candidates based on the value.")
    sRNA_parser.add_argument(
        "--max_length", "-lM",default=500,
        help="Please assign the maxium length of sRNA. "
        "It will classify sRNA candidates based on the value.")
    sRNA_parser.add_argument(
        "--tex_wig_folder", "-tw", default=False,
        help="The path of tex+/- wig folder.")
    sRNA_parser.add_argument(
        "--frag_wig_folder", "-fw", default=False,
        help="The path of fragment wig folder.")
    sRNA_parser.add_argument(
        "--cutoff_intergenic_coverage", "-ci",default=5, type=float,
        help="The cutoff of minimal coverage of intergenic sRNA candidates.")
    sRNA_parser.add_argument(
        "--cutoff_5utr_coverage", "-cu5",default="median",
        help="The cutoff of minimal coverage of 5'UTR derived sRNA candidates.")
    sRNA_parser.add_argument(
        "--cutoff_3utr_coverage", "-cu3",default="median",
        help="The cutoff of minimal coverage of 3'UTR derived sRNA candidates.")
    sRNA_parser.add_argument(
        "--cutoff_interCDS_coverage", "-cuf",default="median",
        help="The cutoff of minimal coverage of inter CDS sRNA candidates.")
    sRNA_parser.add_argument(
        "--fasta_folder", "-f", default=False,
        help="If you want to import secondary structure information, please assign the path of fasta folder. ")
    sRNA_parser.add_argument(
        "--cutoff_energy", "-e", default=0, type=float,
        help="If you want to import secondary structure information, please assign the cutoff of energy. "
        "It will classify sRNA candidates based on the value.")    
    sRNA_parser.add_argument(
        "--mountain_plot", "-m", default=False, action="store_true",
        help="If you want to generate mountain plots of sRNA candidates, please turn it on. ")
    sRNA_parser.add_argument(
        "--database_format", "-fd", default=False, action="store_true",
        help="If you already format your database, you don't need to turn it on.")
    sRNA_parser.add_argument(
        "--sRNA_database_path", "-sd", default=False,
        help="If you want to import blast results of sRNA, please assign the path of sRNA database.")
    sRNA_parser.add_argument(
        "--nr_database_path", "-nd", default=False,
        help="If you want to import blast results of nr, please assign the path of nr database.")
    sRNA_parser.add_argument(
        "--tex_notex_libs", "-tl", default=False, nargs="+", help="library name of tex and notex library.")
    sRNA_parser.add_argument(
        "--frag_libs", "-fl", default=False, nargs="+", help="library name of fragment library.")
    sRNA_parser.add_argument(
        "--tex_notex", "-te", default=2, type=int,
        help="For tex +/- library, sRNA candidates should be detected by both or just one.(1/2)")
    sRNA_parser.add_argument(
        "--replicates", "-r", type=int, help="how many replicates should detect the same sRNA candidates.")
    sRNA_parser.add_argument(
        "--table_best", "-tb", default=False, action="store_true",
        help="the output table of sRNA candidates only print the best track.")
    sRNA_parser.add_argument(
        "--decrease_intergenic","-di", default=0.5, type=float,
        help="if the intergenic region is longer than the max_length, it will based on coverage to check the sRNA candidates. "
             "if the ratio of lowest coverage of intergenic region and the highest coverage of intergenic region is smaller than this number, "
             "it will consider the the point of lowest coverage to be end of sRNA. If the length of sRNA candidate is properly, "
             "it also assign to be one of sRNA candidates.")
    sRNA_parser.add_argument(
        "--decrease_utr","-du", default=0.5, type=float,
        help="if the kind of utr derived is 5'UTR, you have to consider the end of it's end."
             "if the ratio of lowest coverage of it and the highest coverage of it is smaller than this number, "
             "it will consider the the point of lowest coverage to be end of sRNA. If the length of sRNA candidate is properly, "
             "it also assign to be one of sRNA candidates.")
    sRNA_parser.add_argument(
        "--fuzzy_intergenic", "-fi", default=10, type=int,
        help="If the situation is like decrease_intergenic mentioned, the value would be fuzzy between the end of sRNA")
    sRNA_parser.add_argument(
        "--fuzzy_utr", "-fu", default=10, type=int,
        help="If the situation is like decrease_utr mentioned, the value would be fuzzy between the end of sRNA")
    sRNA_parser.add_argument(
        "--cutoff_nr_hit", "-cn", default=0, type=int,
        help="The cutoff of number of hits in nr database. "
        "If the number of nr hits more than this cutoff, program will exclude it during classification.")
    sRNA_parser.add_argument(
        "--sORF", "-O", default=False,
        help="If you want to compare sORF and sRNA, please assign the path of sORF gff folder.")
    sRNA_parser.add_argument(
        "--best_with_all_sRNAhit", "-ba", default=False, action="store_true",
        help="When you wnat to generate the files which store the best sRNA candidates, "
             "it should include all the sRNA candidates which can find the homology from blast sRNA database no matter other information(ex. TSS, blast in nr...)."
             "please turn it on. Or it will just select the the best candidates based on all filter conditions.")
    sRNA_parser.add_argument(
        "--best_without_sORF_candidate", "-bs", default=False, action="store_true",
        help="When you wnat to generate the files which store the best sRNA candidates, "
             "it shouldexclude all the sRNA candidates which also can be detected by sORF file."
             "please turn it on. Or it will just select the the best candidates based on all filter conditions.")
    sRNA_parser.set_defaults(func=run_sRNA_detection)
    # Parameters of small ORF
    sORF_parser = subparsers.add_parser(
        "sorf", help="run sORF detection to detecting sORF candidates which has expression.")
    sORF_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    sORF_parser.add_argument(
        "--UTR_derived_sORF", "-u", default=False, action="store_true",
        help="If you want to detect UTR derived sORF, please turn it on.")
    sORF_parser.add_argument(
        "--transcript_assembly_folder", "-a", default="Transap/output/transcriptome_assembly/gffs",
        help="The path of transcriptome assembly folder.")
    sORF_parser.add_argument(
        "--annotation_folder", "-g", default="Transap/output/target/annotation",
        help="The path of annotation gff folder.")
    sORF_parser.add_argument(
        "--TSS_folder", "-t", default=False,
        help="If you want to import TSS information, please assign the path of gff folder of TSS.")
    sORF_parser.add_argument(
        "--utr_length", "-ul", default=300, type=int,
        help="If you want to import TSS information, please assign the utr length for comparing TSS and sORF. "
        "The default number is 300.")
    sORF_parser.add_argument(
        "--min_length", "-lm",default=30,
        help="Please assign the minium length of sORF. "
        "It will classify sORF candidates based on the value.")
    sORF_parser.add_argument(
        "--max_length", "-lM",default=500,
        help="Please assign the maxium length of sORF. "
        "It will classify sORF candidates based on the value.")
    sORF_parser.add_argument(
        "--tex_wig_folder", "-tw", default=False,
        help="The path of tex+/- wig folder.")
    sORF_parser.add_argument(
        "--frag_wig_folder", "-fw", default=False,
        help="The path of fragment wig folder.")
    sORF_parser.add_argument(
        "--cutoff_intergenic_coverage", "-ci",default=5, type=float,
        help="The cutoff of minimal coverage of sORF candidates.")
    sORF_parser.add_argument(
        "--cutoff_5utr_coverage", "-cu5",default="median",
        help="The cutoff of minimal coverage of sORF candidates.")
    sORF_parser.add_argument(
        "--cutoff_3utr_coverage", "-cu3",default="median",
        help="The cutoff of minimal coverage of sORF candidates.")
    sORF_parser.add_argument(
        "--cutoff_interCDS_coverage", "-cuf",default="median",
        help="The cutoff of minimal coverage of inter CDS sORF candidates.")
    sORF_parser.add_argument(
        "--fasta_folder", "-f", default=False,
        help="If you want to import secondary structure information, please assign the path of fasta folder. ")
    sORF_parser.add_argument(
        "--tex_notex_libs", "-tl", default=False, nargs="+", help="library name of tex and notex library.")
    sORF_parser.add_argument(
        "--frag_libs", "-fl", default=False, nargs="+", help="library name of fragment library.")
    sORF_parser.add_argument(
        "--tex_notex", "-te", default=2, type=int,
        help="For tex +/- library, sORF candidates should be detected by both or just one.(1/2)")
    sORF_parser.add_argument(
        "--replicates", "-r", type=int, help="how many replicates should detect the same sORF candidates.")
    sORF_parser.add_argument(
        "--table_best", "-tb", default=False, action="store_true",
        help="the output table of sORF candidates only print the best track.")
    sORF_parser.add_argument(
        "--sRNA_folder", "-s", default=False,
        help="If you want to compare sORF and sRNA, please assign the path of sORF gff folder.")
    sORF_parser.add_argument(
        "--start_coden", "-ac", default=["ATG"], nargs="+", 
        help="what kinds of start coden ATG/GTG/TTG")
    sORF_parser.add_argument(
        "--stop_coden", "-oc", default=["TTA", "TAG", "TGA"], nargs="+", 
        help="what kinds of stop coden TTA/TAG/TGA")
    sORF_parser.add_argument(
        "--condition_best", "-c", default="TSS",
        help="For generate the result of best sORF, please assign which information you want to consider (TSS/sRNA/both). default is TSS")
    sORF_parser.set_defaults(func=run_sORF_detection)
    # Parameters of promoter detection
    promoter_parser = subparsers.add_parser(
        "promoter", help="run MEME to dicover promoter.")
    promoter_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    promoter_parser.add_argument(
        "--fasta_folder", "-f", default="Transap/output/target/fasta",
        help="please assign the path of gemonic fasta file. ")
    promoter_parser.add_argument(
        "--TSS_folder", "-t", default="Transap/output/TSS/gffs",
        help="The path of TSS gff file.")
    promoter_parser.add_argument(
        "--num_motif", "-n", default=10,
        help="How many of motifs you want to produce.")
    promoter_parser.add_argument(
        "--motif_width", "-w", nargs="+",
        help="motif length - it will refer the value to find the motif. "
        "if you want to detect a range of width, you can insert \"-\" between two values. "
        "for example, 2-10. ")
    promoter_parser.add_argument(
        "--parallel", "-p", default=4,
        help="How many process you want to use to run paralle.")
    promoter_parser.add_argument(
        "--TSS_source", "-s", default=True, action="store_false",
        help="If you generate TSS from other method, please turn it on.")
    promoter_parser.add_argument(
        "--tex_libs", "-tl", default=False, nargs="+", 
        help="library name of tex+/- library. If your TSS is not from Transap, please assign the libs of tex+/- too.")
    promoter_parser.add_argument(
        "--tex_wig_path", "-tw", default=False,
        help="The path of tex+/- wig folder. If your TSS is not from Transap, please assign the wig path too.")
    promoter_parser.add_argument(
        "--annotation_folder", "-g", default=False,
        help="The path of annotation gff folder. If your TSS is not from Transap, please assign the annotation gff path too.")
    promoter_parser.set_defaults(func=run_MEME)
    # Parameters of operon detection
    operon_parser = subparsers.add_parser(
        "Operon", help="detect operon and combine features together.")
    operon_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    operon_parser.add_argument(
        "--TSS_folder", "-t", default="Transap/output/TSS/gffs",
        help="The path of TSS gff folder.")
    operon_parser.add_argument(
        "--annotation_folder", "-g", default="Transap/output/target/annotation",
        help="The path of annotation gff folder.")
    operon_parser.add_argument(
        "--transcript_folder", "-a", default="Transap/output/transcriptome_assembly/gffs",
        help="The path of transcript assembly gff folder.")
    operon_parser.add_argument(
        "--UTR5_folder", "-u5", default="Transap/output/UTR/5UTR/gffs",
        help="The path of 5'UTR gff folder.")
    operon_parser.add_argument(
        "--UTR3_folder", "-u3", default="Transap/output/UTR/3UTR/gffs",
        help="The path of 3'UTR gff folder.")
    operon_parser.add_argument(
        "--term_folder", "-e", default=False, default="Transap/output/terminator/detect",
        help="The path of terminator gff folder.")
    operon_parser.add_argument(
        "--TSS_fuzzy", "-tf", default=3, type=int,
        help="the fuzzy for comparing TSS and transcript assembly. "
        "The default number is 3.")
    operon_parser.add_argument(
        "--term_fuzzy", "-ef", default=10, type=int,
        help="the fuzzy for comparing terminator and transcript assembly. "
        "The default number is 10.")
    operon_parser.add_argument(
        "--min_length", "-l", default=20, type=int,
        help="the minimum length of operon. "
        "The default number is 20.")
    operon_parser.add_argument(
        "--statistics", "-s", default=False, action="store_true",
        help="Doing statistics for Operon analysis. "
        "the name of statistics file is - stat_operon_$STRAIN_NAME.csv.")
    operon_parser.add_argument(
        "--combine_gff", "-c", default=False, action="store_true",
        help="convert the operon and all features you assigned to one gff file. ")
    operon_parser.set_defaults(func=run_operon)
    # Parameters of CircRNA detection
    circrna_parser = subparsers.add_parser(
        "circrna", help="detect circular RNA.")
    circrna_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    circrna_parser.add_argument(
        "--align", "-a", default=False, action="store_true",
        help="Using segemehl to align read (turn on splice detection). "
        "If you already usd segemehl with -S to align your reads, "
        "you can skip this step, don't need to turn it on. "
        "But remember put your results in Transap/output/CircRNA/segemehl_align. "
        "Please be attention, it only use default parameters of segemehl to align your reads. "
        "Moreover, it will align all read files in Transap/input/reads. "
        "If you want to run some specific functions of segemehl, "
        "please run it seperately and copy the files to Transap/output/CircRNA/segemehl_align.")
    circrna_parser.add_argument(
        "--process", "-p", default=10, type=int,
        help="How many processes for --align.  ")
    circrna_parser.add_argument(
        "--fasta_path", "-f", default="Transap/output/target/fasta",
        help="The genome fasta path. ")
    circrna_parser.add_argument(
        "--annotation_path", "-g", default="Transap/output/target/annotation",
        help="The path of the folder of annotation gff files.")
    circrna_parser.add_argument(
        "--convert_to_gff", "-cg", default=False, action="store_true",
        help="If you want to convert circRNA candidates to gff file, please turnn it on.")
    circrna_parser.add_argument(
        "--support_reads", "-s", default=5, type=int,
        help="If you want to convert circRNA candidates to gff file, please also assign the cut off of supported reads.")
    circrna_parser.set_defaults(func=run_circrna)
    # Parameters of Go term
    goterm_parser = subparsers.add_parser(
        "go_term", help="extract and find Go terms.")
    goterm_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    goterm_parser.add_argument(
        "--annotation_path", "-g", default="Transap/output/target/annotation",
        help="The path of annotation gff folder. ")
    goterm_parser.add_argument(
        "--UniProt_id", "-u", default=False,
        help="The path of UniProt ID mapping database. "
        "default is Transap/input/database/idmapping_selected.tab")
    goterm_parser.set_defaults(func=run_goterm)
    # Parameters of sRNA target prediction
    starget_parser = subparsers.add_parser(
        "srna_target", help="sRNA target prediction.")
    starget_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    starget_parser.add_argument(
        "--annotation_path", "-g", default="Transap/output/target/annotation",
        help="The path of annotation gff folder. ")
    starget_parser.add_argument(
        "--fasta_path", "-f", default="Transap/output/target/fasta",
        help="The path of genome fasta folder. ")
    starget_parser.add_argument(
        "--sRNA_path", "-r", default="Transap/output/sRNA/gffs",
        help="The path of sRNA gff folder. ")
    starget_parser.add_argument(
        "--program", "-p", default="both",
        help="Using RNAplex, RNAup or both. Default is both. ")
    starget_parser.add_argument(
        "--interaction_length", "-i", default=30, type=int,
        help="Maximal length of an interaction. Default is 30. ")
    starget_parser.add_argument(
        "--window_size_target", "-wt", default=240, type=int,
        help="only work when --program is RNAplex or both. "
        "Average the pair probabilities over windows of given size. "
        "Default is 240. ")
    starget_parser.add_argument(
        "--span_target", "-st", default=160, type=int,
        help="only work when --program is RNAplex or both. "
        "Set the maximum allowed separation of a base pair to span. "
        "Default is 160. ")
    starget_parser.add_argument(
        "--window_size_srna", "-ws", default=30, type=int,
        help="only work when --program is RNAplex or both. "
        "Average the pair probabilities over windows of given size. "
        "Default is 30. ")
    starget_parser.add_argument(
        "--span_srna", "-ss", default=30, type=int,
        help="only work when --program is RNAplex or both. "
        "Set the maximum allowed separation of a base pair to span. "
        "Defaults to winsize if parameter is omitted. ")
    starget_parser.add_argument(
        "--unstructured_region_RNAplex_target", "-ut", default=30, type=int,
        help="only work when --program is RNAplex or both. "
        "Compute the mean probability that regions of "
        "length 1 to a given length are unpaired. Default is 30. ")
    starget_parser.add_argument(
        "--unstructured_region_RNAplex_srna", "-us", default=30, type=int,
        help="only work when --program is RNAplex or both. "
        "Compute the mean probability that regions of "
        "length 1 to a given length are unpaired. Default is 30. ")
    starget_parser.add_argument(
        "--unstructured_region_RNAup", "-uu", default=30, type=int,
        help="only work when --program is RNAup or both. "
        "Compute the mean probability that regions of "
        "length 1 to a given length are unpaired. Default is 30. ")
    starget_parser.add_argument(
        "--energy_threshold", "-e", default=-8, type=float,
        help="only work when --program is RNAplex or both. "
        "Minimal energy for a duplex to be returned. "
        "Default is -8. ")
    starget_parser.add_argument(
        "--duplex_distance", "-d", default=20, type=int,
        help="only work when --program is RNAplex or both. "
        "Distance between target 3' ends of two consecutive duplexes. "
        "Default is 20. ")
    starget_parser.add_argument(
        "--process_rnaplex", "-pp", default=5, type=int,
        help="How many processes for running RNAplex prediction. ")
    starget_parser.add_argument(
        "--process_rnaup", "-pu", default=20, type=int,
        help="How many processes for running RNAup prediction. ")
    starget_parser.set_defaults(func=sRNA_target)
    # Parameters of SNP transcript
    snp_parser = subparsers.add_parser(
        "snp", help="Detection of SNP of transcripts.")
    snp_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    snp_parser.add_argument(
        "--samtools_path", default="samtools",
        help="If you want to assign the path of samtools, please assign here.")
    snp_parser.add_argument(
        "--bcftools_path", default="bcftools",
        help="If you want to assign the path of bcftools, please assign here.")
    snp_parser.add_argument(
        "--bam_type", "-t",
        help="Please assign the type of BAM. "
             "If your BAM file is mapping to reference genome and you want to know the difference between refenece genome and target genome, plase keyin 'reference'. "
             "If your BAM file already mapped to target genome and you want to check the genome sequence has SNP or not, please keyin 'target'")
    snp_parser.add_argument(
        "--program", "-p", nargs="+",
        help="Please assign the program for detect SNP of transcript: "
             "1: calculate with BAQ, 2: calculate without BAQ, 3: calculate with extend BAQ. "
             "You can assign more than 1 program. For example: 1 2 3")
    snp_parser.add_argument(
        "--fasta_path", "-f", default="Transap/output/target/fasta",
        help="The path of fasta folder. ")
    snp_parser.add_argument(
        "--normal_bam_path", "-nw", default=False,
        help="The path of normal wig folder. If you want to use tex treated and untreated bam files, "
             "please assign the path. Or it will not combine the bam files")
    snp_parser.add_argument(
        "--frag_bam_path", "-fw", default=False,
        help="The path of fragment wig folder. If you want to use tex treated and untreated bam files, "
             "please assign the path. Or it will not combine the bam files")
    snp_parser.add_argument(
        "--quality", "-q", default=20, type=int,
        help="the min quality which consider a real snp.")
    snp_parser.add_argument(
        "--read_depth", "-d", default=False, type=int,
        help="the minimum read depth, below to it will be excluded. default is 5 * number of BAM files,"
        "if the cutoff higher than 40, it will use 40.")
    snp_parser.add_argument(
        "--indel_fraction", "-imf", default=0.5, type=float,
        help="the fraction of maximum read depth, which support insertion of deletion.")
    snp_parser.set_defaults(func=SNP)
    # Parameters of protein-protein interaction network
    ppi_parser = subparsers.add_parser(
        "ppi_network", help="generate protein-protein interaction with literature supported.")
    ppi_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    ppi_parser.add_argument(
        "--ptt_path", "-p", default="Transap/output/target/annotation",
        help="the path of .ptt annotation folder. ")
    ppi_parser.add_argument(
        "--proteinID_strains", "-s", nargs="+",
        help="This is for assigning protein Id which you want to predict. "
             "In order to retrieve the data from STRING and Pubmed, you also have to assign the similar reference. "
             "For example, if you want to run all proteins in Staphylococcus aureus HG003, "
             "you can assign the strain_comparison like "
             "all:Staphylococcus_aureus_HG003.ptt:Staphylococcus_aureus_HG003:\"Staphylococcus aureus 8325\":\"Staphylococcus aureus\". "
             "or all:Staphylococcus_aureus_HG003.ptt:Staphylococcus_aureus_HG003:\"93061\":\"Staphylococcus aureus\". "
             "or all:Staphylococcus_aureus_HG003.ptt:Staphylococcus_aureus_HG003:\"Staphylococcus aureus NCTC 8325\":\"Staphylococcus aureus\". "
             "(ptt_name:STRING_name:Pubmed_name). "
             "First one is the locus tag of ptt. If you want to run all proteins, just assign 'all'. "
             "Therefore, if you want to run SAOUHSC_00001, just replace SAOUHSC_00001 to all like previous example. "
             "Second one is the ptt file name. Third one is the seq name in ptt files. "
             "Fourth one is for STRING database, and the fourth one is for Pubmed. "
             "Of course, you can run the script for several strains at the same time. "
             "Before running it, please check the species file which located in Transap/input/database ."
             "If you didn't download the file, you can download it by yourself or use get_package_database. "
             "You can use taxon_id, STRING_name_compact or official_name_NCBI to represent STRING_name."
             "BE CAREFUL, if the name which you assigned has spaces, please put \"\" at two ends. "
             "For the name of Pubmed, you can assign the name not so specific. "
             "If you assign a specific name, it may not be able to find the related literatures.")
    ppi_parser.add_argument(
        "--without_strain_pubmed", "-n", default=False, action="store_true",
        help="If you want to retrieve pubmed without assign any strains, please turn it on. ")
    ppi_parser.add_argument(
        "--species_STRING", "-d",
        help="Please assign the path of species file of STRING.")
    ppi_parser.add_argument(
        "--score", "-ps", default=0.0, type=float,
        help="Please assign the threadhold of score. The value is from -1 to 1")
    ppi_parser.add_argument(
        "--node_size", "-ns", default=4000, type=int,
        help="Please size of the nodes in figure, default is 4000.")
    ppi_parser.set_defaults(func=PPI)
    # Parameters of subcellular localization
    sub_local_parser = subparsers.add_parser(
        "subcellular_localization", help="prediction of subcellular localization of genomic CDS.")
    sub_local_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    sub_local_parser.add_argument(
        "--gff_path", "-g", default="Transap/output/target/annotation",
        help="The path of annotation gff folder. ")
    sub_local_parser.add_argument(
        "--fasta_path", "-f", default="Transap/output/target/fasta",
        help="The path of fasta folder. ")
    sub_local_parser.add_argument(
        "--bacteria_type", "-b",
        help="Is Gram-positive or Gram-negative. Please assign 'positive' or 'negative'. ")
    sub_local_parser.add_argument(
        "--merge_to_gff", "-m",  default=False, action="store_true",
        help="If you want to merge the information to annotation gff file, please turn it on. ")
    sub_local_parser.set_defaults(func=Sub_Local)
    # Parameters of riboswitch
    ribos_parser = subparsers.add_parser(
        "riboswitch", help="prediction of riboswitch.")
    ribos_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    ribos_parser.add_argument(
        "--riboswitch_ID", "-i",
        help="the path of the riboswitch ID of Rfam. ")
    ribos_parser.add_argument(
        "--gff_path", "-g", default="Transap/output/target/annotation",
        help="The path of annotation gff folder. ")
    ribos_parser.add_argument(
        "--fasta_path", "-f", default="Transap/output/target/fasta",
        help="The path of fasta folder. ")
    ribos_parser.add_argument(
        "--Rfam", "-R",
        help="The path of Rfam CM database. ")
    ribos_parser.add_argument(
        "--re_scan", "-r", default=False, action="store_true",
        help="Based on the results of first scaning, it will modify the input of sequence and re scan again. ")
    ribos_parser.add_argument(
        "--fuzzy", "-z", default=10, type=int,
        help="It will extend some nts of 3' ans 5' end. ")
    ribos_parser.set_defaults(func=Ribos)
    # Parameters of generate screenshots
    screen_parser = subparsers.add_parser(
        "screenshot", help="generate screenshot for selected feature.")
    screen_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    screen_parser.add_argument(
        "--main_gff", "-mg",
        help="It will based on the position of this gff file to generate screenshot. ")
    screen_parser.add_argument(
        "--side_gffs", "-sg", default=False, nargs="+",
        help="If you have more than one gff want to plot, please assign here. ")
    screen_parser.add_argument(
        "--fasta", "-f", default="Transap/output/target/fasta",
        help="The path of genome fasta folder. ")
    screen_parser.add_argument(
        "--frag_wig_folder", "-fw", default=False,
        help="If you want to include the information of fragmented wig file, please assign the folder. ")
    screen_parser.add_argument(
        "--tex_wig_folder", "-tw", default=False,
        help="If you want to include the information of tex treated wig file, please assign the folder. ")
    screen_parser.add_argument(
        "--height", "-he", default=1500, type=int,
        help="you can assign the height of screenshot. ")
    screen_parser.add_argument(
        "--tex_libs", "-tl", default=False, nargs="+",
        help="If you want to include the tex+/- library, please assign here.")
    screen_parser.add_argument(
        "--frag_libs", "-fl", default=False, nargs="+",
        help="If you want to include the fragmented library, please assign here.")
    screen_parser.add_argument(
        "--present", "-p", default="expand",
        help="Which type you want to present in the screen shot. expand/collapse/squish. ")
    screen_parser.add_argument(
        "--output_folder", "-o",
        help="Please assign the output folder. If the folder does not exist, it will generate automatically. ")
    screen_parser.set_defaults(func=Screen)

    args = parser.parse_args()
    
    if args.version is True:
        print("Transap version " + __version__)
    elif "func" in dir(args):
        controller = Controller(args)
        args.func(controller)
    else:
        parser.print_help()



def create_project(controller):
    controller.create_project(__version__)

def get_input(controller):
    controller.get_input()

def get_target_fasta(controller):
    controller.get_target_fasta()

def run_RATT(controller):
    controller.ratt()

def multiparser(controller):
    controller.multiparser()

def run_TSSpredator(controller):
    controller.tsspredator()

def optimize_TSSpredator(controller):
    controller.optimize()

def color_png(controller):
    controller.color()

def run_TransTermHP(controller):
    controller.transtermhp()

def run_Transcript_Assembly(controller):
    controller.transcript()

def run_UTR_detection(controller):
    controller.utr_detection()

def run_sRNA_detection(controller):
    controller.srna_detection()

def run_sORF_detection(controller):
    controller.sorf_detection()

def run_MEME(controller):
    controller.meme()

def run_operon(controller):
    controller.operon()

def run_circrna(controller):
    controller.circrna()

def run_goterm(controller):
    controller.goterm()

def sRNA_target(controller):
    controller.srna_target()

def SNP(controller):
    controller.snp()

def PPI(controller):
    controller.ppi()

def Sub_Local(controller):
    controller.sublocal()

def Ribos(controller):
    controller.ribos()

def Screen(controller):
    controller.screen()
if __name__ == "__main__":
    main()
