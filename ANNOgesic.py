#!/usr/bin/python

"""transcript annotation and analysis pipeline"""
import argparse
import os
from annogesiclib.controller import Controller

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"
__version__ = "0.1.0"

def main():
    home_path = os.environ["HOME"]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--version", "-v", default=False, action="store_true",
        help="show version")
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
        help="Download fasta files of reference. Default is False")
    get_input_parser.add_argument(
        "--ref_gff", "-g", default=False, action="store_true",
        help="Download gff files of reference. Default is False")
    get_input_parser.add_argument(
        "--ref_gbk", "-k", default=False, action="store_true",
        help="Download genbank files of reference. Default is False")
    get_input_parser.add_argument(
        "--convert_embl", "-e", default=False, action="store_true",
        help="Convert gbk to embl files of reference. Default is False")
    get_input_parser.set_defaults(func=get_input)
    # get target fasta
    get_target_fasta_parser = subparsers.add_parser(
        "get_target_fasta", help="Get target fasta.")
    get_target_fasta_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    get_target_fasta_parser.add_argument(
        "--ref_fasta_folder", "-r", default="ANNOgesic/input/reference/target/fasta",
        help="Folder of fasta files.")
    get_target_fasta_parser.add_argument(
        "--mutation_table", "-m", default="ANNOgesic/input/mutation_table",
        help="The path of mutation table.")
    get_target_fasta_parser.add_argument(
        "--output_format", "-o", default=None, nargs="+",
        help="Please assign the output filename and which strain should be included in it. "
        "For example: FILE1:strain1,strain2. FILE1 is a output fasta file which include the information of strain1 and strain2.")
    get_target_fasta_parser.set_defaults(func=get_target_fasta)    

    # run RATT
    RATT_parser = subparsers.add_parser(
        "annotation_transfer", help="Run RATT to transfer the annotation files from reference to target.")
    RATT_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    RATT_parser.add_argument(
        "--PAGIT_folder", default="/tools/PAGIT",
        help="Path of the PAGIT folder. If your PAGIT folder is not located at /tools/PAGIT, please assign here.")
    RATT_parser.add_argument(
        "--RATT_path", default="start.ratt.sh",
        help="Path of the start.ratt.sh file of RATT folder. Default is start.ratt.sh.")
    RATT_parser.add_argument(
        "--compare_pair", "-p", nargs="+",
        help="Please assign the name of pairs. ex. $REFERENCE.fa:$TARGET.fa. "
        "ATTENTION:please make sure the ref name is the same as embl file. ")
    RATT_parser.add_argument(
        "--element", "-e",
        help="It will become the prefix of all output file.")
    RATT_parser.add_argument(
        "--transfer_type", "-t", default="Strain",
        help="The transfer type for running RATT.(details can refer to the manual of RATT.) default is Strain.")
    RATT_parser.add_argument(
        "--ref_gbk", "-re", default="ANNOgesic/input/reference/annotation",
        help="The folder which stores every reference embl folders."
        "If you have no embl folder, you can assign the folder of genbank.")
    RATT_parser.add_argument(
        "--ref_fasta", "-rf",
        help="The folder of reference fasta files.")
    RATT_parser.add_argument(
        "--target_fasta", "-tf", default="ANNOgesic/output/target/fasta",
        help="The folder which stores target fasta files.")
    RATT_parser.add_argument(
        "--convert_to_gff_rnt_ptt", "-g", default=False, action="store_true",
        help="Do you want to convert to gff, rnt and ptt? Default is False")
    RATT_parser.set_defaults(func=run_RATT)
    # gene expression analysis
    expression_parser = subparsers.add_parser(
        "expression_analysis", help="Run gene expression analysis to compare which CDS is expressing in which libraries ")
    expression_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    expression_parser.add_argument(
        "-g","--annotation_folder", 
        help="The folder of annotation file which you want to analyze.")
    expression_parser.add_argument(
        "-tl","--tex_notex_libs", default=None, nargs="+",
        help="Library name of tex and notex library. The format is: "
        "wig_file_name:tex_treat_or_not(tex or notex):condition_id(integer):replicate_id(alphabet):strand(+ or -).")
    expression_parser.add_argument(
        "-fl","--frag_libs", default=None, nargs="+",
        help="Library name of fragment library. The format is: "
        "wig_file_name:fragmented(frag):condition_id(integer):replicate_id(alphabet):strand(+ or -).")
    expression_parser.add_argument(
        "-te","--tex_notex", default=2, type=int,
        help="For tex +/- library, expressing CDS should be detected by both or just one.(1 or 2)")
    expression_parser.add_argument(
        "-rt","--replicates_tex", default=None,
        help="The expressing CDS of tex +/- library should be detected more than this number of replicates.")
    expression_parser.add_argument(
        "-rf","--replicates_frag", default=None,
        help="The expressing CDS of fragmented library should be detected more than this number of replicates.")
    expression_parser.add_argument(
        "--tex_wig_folder", "-tw", default=None,
        help="The folder of TEX+/- wigge files.")
    expression_parser.add_argument(
        "--frag_wig_folder", "-fw", default=None,
        help="The folder of fragmented wigge files.")
    expression_parser.add_argument(
        "--cutoff_overlap_tex", "-ot", default="all",
        help="This value is for decision of CDS which is expressing or not in TEX+/- library. If the expressing nts more than this value, "
        "it will consider the CDS is expressing one. You can assign by percentage or nucleotide. ex: "
        "p_0.5 means the percentage of expressing nts should higher 0.5. n_100 means there should be 100 nts which are expressing. "
        "Default is \"all\" which means as long as there is a nt's coverage higher than cutoff_coverage, it would consider the CDS which is expressing.")
    expression_parser.add_argument(
        "--cutoff_overlap_frag", "-of", default="all",
        help="This value is for decision of CDS which is expressing or not in fragmented library. If the expressing nts more than this value, "
        "it will consider the CDS is expressing. You can assign by percentage or nucleotide. ex: "
        "p_0.5 means the percentage of expressing nts should higher 0.5. n_100 means there should be 100 nts which are expressing. "
        "Default is \"all\" which means as long as there is a nt's coverage higher than cutoff_coverage, it would consider the CDS which is expressing.")
    expression_parser.add_argument(
        "--cutoff_coverage", "-c", default=5, type=float,
        help="If the coverage is higher than this value, it will consider the nt is expressing")
    expression_parser.add_argument(
        "--features", "-f", nargs="+",
        help="The features which you want to compute, ex: CDS tRNA")
    expression_parser.set_defaults(func=run_expression)
    # Parameters of TSSpredator
    TSSpredator_parser = subparsers.add_parser(
        "tsspredator", help="Run TSSpredator to predict TSSs or processing sites.")
    TSSpredator_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    TSSpredator_parser.add_argument(
        "--TSSpredator_path", default="TSSpredator.jar",
        help="If you want to assign the path of TSSpredator, please assign here.")
    TSSpredator_parser.add_argument(
        "--fasta_folder", "-f", default="ANNOgesic/output/target/fasta",
        help="Path of the target fasta folder.")
    TSSpredator_parser.add_argument(
        "--annotation_folder", "-g", default="ANNOgesic/output/target/annotation",
        help="Path of the target gff folder.")
    TSSpredator_parser.add_argument(
        "--wig_folder", "-w", default="ANNOgesic/input/wigs/tex_notex",
        help="The folder of the wig folder.")
    TSSpredator_parser.add_argument(
        "--height", "-he", default=0.3,
        help="This value relates to the minimal number of read starts at a "
	"certain genomic position to be considered as a TSS candidate. Default is 0.3")
    TSSpredator_parser.add_argument(
        "--height_reduction", "-rh", default=0.2,
        help="When comparing different strains/conditions and "
	"the step height threshold is reached in at least one strain/condition, "
	"the threshold is reduced for the other strains/conditions "
	"by the value set here. "
        "This value must be smaller than the step height threshold. Default is 0.2")
    TSSpredator_parser.add_argument(
        "--factor", "-fa", default=2.0,
        help="This is the minimal factor by which the TSS height has to "
	"exceed the local expression background. Default is 2.0")
    TSSpredator_parser.add_argument(
        "--factor_reduction", "-rf", default=0.5,
        help="When comparing different strains/conditions and "
        "the step factor threshold is reached in at least one strain/condition, "
        "the threshold is reduced for the other strains/conditions "
        "by the value set here. "
        "This value must be smaller than the step factor threshold. Default is 0.5")
    TSSpredator_parser.add_argument(
        "--enrichment_factor", "-ef", default=2.0,
        help="This is the minimal enrichment factor. "
        "During optimization will never larger than this value. Default is 2.0")
    TSSpredator_parser.add_argument(
        "--processing_factor", "-pf", default=1.5,
        help="This is the minimal processing factor. If untreated library is higher than the treated library "
        "and above which the TSS candidate is considered as a processing site and not annotated as detected. "
        "During optimization will never larger than this value. Default is 1.5")
    TSSpredator_parser.add_argument(
        "--base_height", "-bh", default=0,
        help="This is the minimal number of reads should be mapped on TSS. Default is 0.0")
    TSSpredator_parser.add_argument(
        "--replicate_match", "-rm", default=1,
        help="The TSS candidates should match to how many number of the replicates. Default is 1.")
    TSSpredator_parser.add_argument(
        "--utr_length", "-u", default=300, type=int,
        help="The length of UTR. It is for Primary and Secondary definition. Default is 300")
    TSSpredator_parser.add_argument(
        "--lib", "-l", nargs='+',
        help="The libraries of wig files for TSSpredator. The format is: "
	"wig_file_name:tex_treat_or_not(tex or notex):condition_id(integer):replicate_id(alphabet):strand(+ or -).")
    TSSpredator_parser.add_argument(
        "--output_prefix", "-p", nargs='+',
        help="The output prefix of all conditions.")
    TSSpredator_parser.add_argument(
        "--merge_manual", "-m", default=None,
        help="If you have gff file of manual checked TSS, you can use this function to merge manual checked ones and predicted ones.")
    TSSpredator_parser.add_argument(
        "--statistics", "-s", default=False, action="store_true",
        help="Doing statistics for TSS candidates. "
        "it will store in statistics folder. Default is False")
    TSSpredator_parser.add_argument(
        "--validate_gene", "-v", default=False, action="store_true",
        help="Using TSS candidates to validate genes in annotation file. "
        "it will store in statistics folder. Default is False")
    TSSpredator_parser.add_argument(
        "--compute_program", "-t", default="TSS",
        help="Which program do you want to predict. (TSS or processing_site)")
    TSSpredator_parser.add_argument(
        "--compare_transcript_assembly", "-ta", default=None,
        help="If you want to compare with transcriptome assembly, please assign the folder of gff file of transcript assembly."
        "Default is False. ")
    TSSpredator_parser.add_argument(
        "--fuzzy", "-fu", default=5, type=int,
        help="The fuzzy for comparing TSS and transcript assembly. Default is 5")
    TSSpredator_parser.add_argument(
        "--cluster", "-c", default=2, type=int,
        help="This number is for compare manual detected TSS and prediced one. "
        "If the position between manual checked one and predicted one is smaller or equal than this value, "
        "It will only print one of them. Default is 2")
    TSSpredator_parser.add_argument(
        "--length", "-le", default=None,
        help="The length that you want to compare with manual check for statistics. If you want to compare whole genome, "
        "please don't turn it on. The default is comparing whole genome")
    TSSpredator_parser.add_argument(
        "--re_check_orphan", "-ro", default=False, action="store_true",
        help="If your annotation file lack information of gene or locus_tag, you can turn it on. "
        "It will try to compare with CDS. Default is False")
    TSSpredator_parser.add_argument(
        "--overlap_feature", "-of", default="both",
        help="If processing site and TSS are overlap, you can keep \"TSS\" or \"processing_site\" or \"both\". "
        "Default is both.")
    TSSpredator_parser.add_argument(
        "--reference_gff_folder", "-rg", default=None,
        help="For overlap_feature, if you want to only keep \"TSS\" or \"processing_site\", "
        "please assign the reference_gff_folder. If you are running TSS, please assign the folder of processing site. "
        "If you are running processing_site, please assign the folder of TSS. If you want to keep \"both\" at overlap position, "
        "please don't turn it on.")
    TSSpredator_parser.add_argument(
        "--remove_low_expression", "-rl", default=None,
        help="If you want to remove low expressed TSS/processing site, please assign the file of manual checked gff file here. "
        "Please Be ATTENTION: this parameter may remove some True positive, too. "
        "So, please make sure you want to do it.")
    TSSpredator_parser.set_defaults(func=run_TSSpredator)
    # Parameter of opimization of TSSpredator
    op_TSSpredator_parser = subparsers.add_parser(
        "optimize_tsspredator", help="Optimize TSSpredator based on (partial)manual detect one.")
    op_TSSpredator_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    op_TSSpredator_parser.add_argument(
        "--TSSpredator_path", default="TSSpredator.jar",
        help="If you want to assign the path of TSSpredator, please assign here.")
    op_TSSpredator_parser.add_argument(
        "--fasta_file", "-fs",
        help="Path of one target fasta file which you want to opimize it.")
    op_TSSpredator_parser.add_argument(
        "--annotation_file", "-g",
        help="Path of one target gff file which you want to opimize it.")
    op_TSSpredator_parser.add_argument(
        "--wig_folder", "-w", default="ANNOgesic/input/wigs/tex_notex",
        help="The folder of the wig folder.")
    op_TSSpredator_parser.add_argument(
        "--manual", "-m",
        help="The file of manual checked gff file.")
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
        "--max_enrichment_factor", "-ef", default=6.0, type=float,
        help="This is the minimal enrichment factor. "
        "During optimization will never larger than this value.")
    op_TSSpredator_parser.add_argument(
        "--max_processing_factor", "-pf", default=6.0, type=float,
        help="This is the minimal processing factor. If untreated library is higher than the treated library " 
        "and above which the TSS candidate is considered as a processing site and not annotated as detected. "
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
        "--length", "-le", default=None,
        help="The length of nts for running optimization."
        "Default is compare whole genome")
    op_TSSpredator_parser.add_argument(
        "--core", "-c", type=int, default=4,
        help="How many paralle running do you want to use. Default is 4")
    op_TSSpredator_parser.add_argument(
        "--program", "-r", default="TSS",
        help="The type which you want to run TSSpredator (TSS or Processing_site). Default is TSS")
    op_TSSpredator_parser.add_argument(
        "--replicate_match", "-rm", default=1,
        help="The TSS candidates should match to how many number of the replicates. Default is 1")
    op_TSSpredator_parser.add_argument(
        "--steps", "-s", default=4000, type=int,
        help="How many steps do you want to run. Default is 4000 runs.")
    op_TSSpredator_parser.set_defaults(func=optimize_TSSpredator)
    # Parameter of generating color png
    color_parser = subparsers.add_parser(
        "color_png", help="Generating color screenshots of TSS or processing site. "
        "It only works after running batch script. ")
    color_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    color_parser.add_argument(
        "--output_folder", "-f",
        help="The folder which stores the screenshots. ")
    color_parser.add_argument(
        "--track_number", "-t", type=int,
        help="How many number of tracks do you have. ")
    color_parser.add_argument(
        "--ImageMagick_covert_path", "-m", default="convert",
        help="Please assign the path of covert in ImageMagick package. ")
    color_parser.set_defaults(func=color_png)
    # Parameter of Terminator
    Terminator_parser = subparsers.add_parser(
        "terminator", help="Detect Terminators.")
    Terminator_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    Terminator_parser.add_argument(
        "--TransTermHP_path", default="transterm",
        help="Please assign the path of transterm in TransTermHP.")
    Terminator_parser.add_argument(
        "--expterm_path", default="expterm.dat",
        help="Please assign the path of your expterm.dat.")
    Terminator_parser.add_argument(
        "--RNAfold_path", default="RNAfold",
        help="If you want to assign the path of RNAfold of Vienna package, please assign here.")
    Terminator_parser.add_argument(
        "--fasta_folder", "-f", default="ANNOgesic/output/target/fasta",
        help="The path of fasta folder.")
    Terminator_parser.add_argument(
        "--annotation_folder", "-g", default="ANNOgesic/output/target/annotation",
        help="The path of annotation gff folder.")
    Terminator_parser.add_argument(
        "--transcript_folder", "-a", default="ANNOgesic/output/transcriptome_assembly/gffs",
        help="The path of the folder which store gff files of transcript assembly.")
    Terminator_parser.add_argument(
        "--sRNA", "-sr", default=None,
        help="If you want to include sRNA information, please assign the folder of gff files of sRNA.")
    Terminator_parser.add_argument(
        "--statistics", "-s", default=False, action="store_true",
        help="Doing statistics for TransTermHP. "
        "the name of statistics file is - stat_terminator_$STRAIN_NAME.csv.")
    Terminator_parser.add_argument(
        "--tex_wig_folder", "-tw", default=None,
        help="If you want to use tex +/- libraries, please assign tex +/- wig folder.")
    Terminator_parser.add_argument(
        "--frag_wig_folder", "-fw", default=None,
        help="If you want to use fragmented libraries, please assign fragmented wig folder.")
    Terminator_parser.add_argument(
        "--decrease", "-d", default=0.5, type=float,
        help="If the (lowest coverage / highest coverage) in the terminator is smaller than this number, "
        "it will consider this terminator have dramatic coverage decrease in it.")
    Terminator_parser.add_argument(
        "--fuzzy_detect_coverage", "-fc", default=30, type=int,
        help="It will elongate the number of nt(you assign here) from both terminal site. "
        "If it can found the coverage dramatic decrease within this range, it will consider the terminator have dramatic coverage decrease in it.")
    Terminator_parser.add_argument(
        "--fuzzy_upstream_transcript", "-fut", default=30, type=int,
        help="If the candidates are upstream of transcript and the distance between the end of gene and terminator candidate is within this number, "
             "it will be consider as terminator.")
    Terminator_parser.add_argument(
        "--fuzzy_downstream_transcript", "-fdt", default=30, type=int,
        help="If the candidates are downstream of transcript and the distance is within this number, it will be consider as terminator.")
    Terminator_parser.add_argument(
        "--fuzzy_upstream_cds", "-fuc", default=10, type=int,
        help="If the candidates are upstream of CDS/tRNA/rRNA/sRNA and the distance between the end of gene and terminator candidate is within this number, "
             "it will be consider as terminator.")
    Terminator_parser.add_argument(
        "--fuzzy_downstream_cds", "-fdg", default=310, type=int,
        help="If the candidates are downstream of CDS/tRNA/rRNA/sRNA and the distance is within this number, it will be consider as terminator.")
    Terminator_parser.add_argument(
        "--highest_coverage", "-hc", default=2.5, type=float,
        help="If the highest coverage of the region of terminator is below to this number, "
        "the terminator will be classify to non-detect.")
    Terminator_parser.add_argument(
        "-tl","--tex_notex_libs", default=None, nargs="+",
        help="Library name of tex and notex library. The format is: "
        "wig_file_name:tex_treat_or_not(tex or notex):condition_id(integer):replicate_id(alphabet):strand(+ or -).")
    Terminator_parser.add_argument(
        "-fl","--frag_libs", default=None, nargs="+", 
        help="Library name of fragment library. The format is: "
        "wig_file_name:fragmented(frag):condition_id(integer):replicate_id(alphabet):strand(+ or -).")
    Terminator_parser.add_argument(
        "-te","--tex_notex", default=2, type=int,
        help="For tex +/- library, terminators should be detected by both or just one.(1/2)")
    Terminator_parser.add_argument(
        "-rt","--replicates_tex", default=None,
        help="The terminator of tex +/- library should be detected more than this number of replicates.")
    Terminator_parser.add_argument(
        "-rf","--replicates_frag", default=None,
        help="The terminator of fragmented library should be detected more than this number of replicates.")
    Terminator_parser.add_argument(
        "-tb","--table_best", default=False, action="store_true",
        help="Output sRNA table only most decreasing track.")
    Terminator_parser.set_defaults(func=run_Terminator)

    # Parameter of Transcript assembly
    Transcript_parser = subparsers.add_parser(
        "transcript_assembly", help="Run Transcript for doing transcriptome assembly.")
    Transcript_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    Transcript_parser.add_argument(
        "--annotation_folder", "-g", default=None,
        help="It is for comparing transcript assembly and annotation gff file. "
        "It can use annotation gff file as reference and modify transcript assembly file. "
        "If you want to do it, please assign the annotation gff folder. Otherwise, don't turn it on.")
    Transcript_parser.add_argument(
        "--sort_annotation", "-s", default=False, action="store_true",
        help="The annotation gff files in annotation folder are sorted or not. "
        "If they didn't be sorted, please turn it on. Default is False")
    Transcript_parser.add_argument(
        "--length", "-l", default=20, type=int,
        help="The minimum width of region to be a transcript. It is for refer to annotation file. "
        "If you want to compare with annotation files, it will be the final output. "
        "If you don't want to compare with annotation files, --width would be the length for the final output. "
        "The default is 20.")
    Transcript_parser.add_argument(
        "--normal_wig_path", "-nw", default=None,
        help="The path of normal wig folder.")
    Transcript_parser.add_argument(
        "--frag_wig_path", "-fw", default=None,
        help="The path of fragment wig folder.")
    Transcript_parser.add_argument(
        "--height", "-he", default=5, type=int,
        help="The minimum height of coverage to be a transcript. "
        "The default is 5.")
    Transcript_parser.add_argument(
        "--width", "-w", default=20, type=int,
        help="The minimum width of region to be a transcript. It is for without annotation to be reference. "
        "If you don't want to compare with annotation files (--length), it will be the final output. "
        "Otherwise, --length would be the length of transcript for the final output."
        "The default is 20.")
    Transcript_parser.add_argument(
        "--tolerance", "-t", default=5, type=int,
        help="This number indicates how willing the algorithm is to ignore a temporary drop below this number. "
        "The default is 5.")
    Transcript_parser.add_argument(
        "--tolerance_coverage", "-tc", default=0, type=int,
        help="If the coverage is lower than tolerance_coverage, even the range is within tolerance, it will terminate the current transcript. "
        "The default is 0.")
    Transcript_parser.add_argument(
        "--replicates_tex", "-rt", default=None,
        help="The position is included in the transcript if there are more than the replicate which you assign here to supported it. (for tex +/- library)")
    Transcript_parser.add_argument(
        "--replicates_frag", "-rf", default=None,
        help="The position is included in the transcript if there are more than the replicate which you assign here to supported it. (for fragmented library)")
    Transcript_parser.add_argument(
        "--tex_notex", "-te", default=2, type=int,
        help="If you use tex +/- libraries to run transcript assembly, please assign the tex +/- should both consider or just one. (1 or 2). Default is 2")
    Transcript_parser.add_argument(
        "--compare_TSS", "-ct", default=None,
        help="If you want to compare with TSS, please assign TSS folder.")
    Transcript_parser.add_argument(
        "--compare_CDS", "-cg", default=None,
        help="If you want to compare with annotation file, please assign annotation folder.")
    Transcript_parser.add_argument(
        "--TSS_fuzzy", "-fu", default=5, type=int,
        help="The fuzzy for comparing TSS and transcript assembly. Default is 5")
    Transcript_parser.add_argument(
        "--Tex_treated_libs", "-tl", nargs="+", default=None,
        help="Input of tex +/- library. The format is: "
        "wig_file_name:tex_treat_or_not(tex or notex):condition_id(integer):replicate_id(alphabet):strand(+ or -).")
    Transcript_parser.add_argument(
        "--fragmented_libs", "-fl", nargs="+", default=None,
        help="Input of fragmented library. The format is: "
        "wig_file_name:fragmented(frag):condition_id(integer):replicate_id(alphabet):strand(+ or -).")
    
    Transcript_parser.set_defaults(func=run_Transcript_Assembly)

    # Parameter of UTR detection
    UTR_parser = subparsers.add_parser(
        "utr", help="Run UTR detection to detect 5'UTR and 3'UTR.")
    UTR_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    UTR_parser.add_argument(
        "--annotation_folder", "-g", default="ANNOgesic/output/target/annotation",
        help="The path of annotation gff folder.")
    UTR_parser.add_argument(
        "--TSS_folder", "-t", default="ANNOgesic/output/TSS/gffs",
        help="The path of TSS folder.")
    UTR_parser.add_argument(
        "--transcript_assembly_folder", "-a", default="ANNOgesic/output/transcriptome_assembly/gffs",
        help="The path of transcriptome assembly folder.")
    UTR_parser.add_argument(
        "--terminator_folder", "-e", default=None,
        help="If you want to add the information of terminator, you can assign the path of terminator folder here.")
    UTR_parser.add_argument(
        "--terminator_fuzzy", "-f", default=30, type=int,
        help="If the distance(nt) between terminator and the end of transcript assembly belows to this value, "
        "it will assign the terminator associated with the 3'UTR.")
    UTR_parser.add_argument(
        "--TSS_source", "-s", default=True, action="store_false",
        help="If you generate TSS from other method not from ANNOgesic, please turn it on.")
    UTR_parser.add_argument(
        "--base_5UTR", "-b", default="both",
        help="Which kind of information that you want to use for generating 5'UTR. TSS/transcript/both. Default is both.")
    UTR_parser.set_defaults(func=run_UTR_detection)

    # Parameter of sRNA detectopn
    sRNA_parser = subparsers.add_parser(
        "srna", help="Run sRNA detection to detect sRNA candidates.")
    sRNA_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    sRNA_parser.add_argument(
        "--Vienna_folder", default="",
        help="Please assign the folder of Vienna package. "
        "It should include RNAfold.")
    sRNA_parser.add_argument(
        "--Vienna_utils", default="",
        help="Please assign the folder of Utils of Vienna package. "
        "It should include relplot.pl and mountain.pl.")
    sRNA_parser.add_argument(
        "--blast_plus_folder", default="",
        help="Please assign the folder of blast+ which include blastn, blastx, makeblastdb.")
    sRNA_parser.add_argument(
        "--ps2pdf14_path", default="ps2pdf14",
        help="Please assign the path of ps2pdf14.")
    sRNA_parser.add_argument(
        "--UTR_derived_sRNA", "-u", default=False, action="store_true",
        help="If you want to detect UTR derived sRNA, please turn it on. Default is False.")
    sRNA_parser.add_argument(
        "--import_info", "-d", nargs="+",
        help="There are several types of information you can import to detect and filter sRNA: "
        "TSS(1), energy of secondary structure(2), blast to nr(3), blast to sRNA(4), sORF(5), "
        "without any information, only detect sRNA (without any information) by transcriptome assembly(6)."
        "Please assign the number of type you want to import, i.e. 1 2 4 - means it used TSS, energy and blast result to detect sRNA. "
        "Besides these information, it will also consider the sequence length of sRNA.")
    sRNA_parser.add_argument(
        "--transcript_assembly_folder", "-a", default="ANNOgesic/output/transcriptome_assembly/gffs",
        help="The path of transcriptome assembly folder.")
    sRNA_parser.add_argument(
        "--annotation_folder", "-g", default="ANNOgesic/output/target/annotation",
        help="The path of annotation gff folder.")
    sRNA_parser.add_argument(
        "--TSS_folder", "-t", default=None,
        help="If you want to import TSS information, please assign the path of gff folder of TSS. "
        "If you want to detect UTR derived sRNA, you must assign the folder of TSS.")
    sRNA_parser.add_argument(
        "--processing_site_folder", "-p", default=None,
        help="If you want to import processing site information, please assign the path of gff folder of processing site."
        "If you want to detect UTR derived sRNA, you must assign the folder of processing site.")
    sRNA_parser.add_argument(
        "--TSS_intergenic_fuzzy", "-ft", default=2, type=int,
        help="If you want to import TSS information, you need to assign the fuzzy for comparing TSS and transcript assembly. It is for intergenic."
        "The default number is 2.")
    sRNA_parser.add_argument(
        "--TSS_5UTR_fuzzy", "-f5", default=2, type=int,
        help="If you want to import TSS information, you need to assign the fuzzy for comparing TSS and transcript assembly. It is for 5'UTR of UTR derived sRNA."
        "The default number is 2.")
    sRNA_parser.add_argument(
        "--TSS_3UTR_fuzzy", "-f3", default=10, type=int,
        help="If you want to import TSS information, you need to assign the fuzzy for comparing TSS and transcript assembly. It is for 3'UTR of UTR derived sRNA."
        "The default number is 10.")
    sRNA_parser.add_argument(
        "--TSS_interCDS_fuzzy", "-fc", default=10, type=int,
        help="If you want to import TSS information, you need to assign the fuzzy for comparing TSS and transcript assembly. It is for interCDS derived sRNA."
        "The default number is 10.")
    sRNA_parser.add_argument(
        "--min_length", "-lm",default=30,
        help="Please assign the minium length of sRNA. "
        "It will classify sRNA candidates based on the value. Default is 30.")
    sRNA_parser.add_argument(
        "--max_length", "-lM",default=500,
        help="Please assign the maxium length of sRNA. "
        "It will classify sRNA candidates based on the value. Default is 500.")
    sRNA_parser.add_argument(
        "--tex_wig_folder", "-tw", default=None,
        help="The path of tex+/- wig folder.")
    sRNA_parser.add_argument(
        "--frag_wig_folder", "-fw", default=None,
        help="The path of fragment wig folder.")
    sRNA_parser.add_argument(
        "--cutoff_intergenic_coverage", "-ci",default=5, type=float,
        help="The cutoff of minimal coverage of intergenic sRNA candidates. Default is 5")
    sRNA_parser.add_argument(
        "--cutoff_5utr_coverage", "-cu5",default="median",
        help="The cutoff of minimal coverage of 5'UTR derived sRNA candidates. "
        "You can also assign median or mean. Default is median")
    sRNA_parser.add_argument(
        "--cutoff_3utr_coverage", "-cu3",default="median",
        help="The cutoff of minimal coverage of 3'UTR derived sRNA candidates. "
        "You can also assign median or mean. Default is median")
    sRNA_parser.add_argument(
        "--cutoff_interCDS_coverage", "-cuf",default="median",
        help="The cutoff of minimal coverage of inter CDS sRNA candidates. "
        "You can also assign median or mean. Default is median")
    sRNA_parser.add_argument(
        "--fasta_folder", "-f", default=None,
        help="If you want to import secondary structure information, please assign the path of fasta folder. ")
    sRNA_parser.add_argument(
        "--cutoff_energy", "-e", default=0, type=float,
        help="If you want to import secondary structure information, please assign the cutoff of folding energy. "
        "It will classify sRNA candidates based on the value. Default is 0.")
    sRNA_parser.add_argument(
        "--mountain_plot", "-m", default=False, action="store_true",
        help="If you want to generate mountain plots of sRNA candidates, please turn it on. Default is False.")
    sRNA_parser.add_argument(
        "--database_format", "-fd", default=False, action="store_true",
        help="If you already format your database, you don't need to turn it on. Default is False")
    sRNA_parser.add_argument(
        "--sRNA_database_path", "-sd", default=None,
        help="If you want to import blast results of sRNA, please assign the path of sRNA database.")
    sRNA_parser.add_argument(
        "--nr_database_path", "-nd", default=None,
        help="If you want to import blast results of nr, please assign the path of nr database.")
    sRNA_parser.add_argument(
        "--sRNA_blast_stat", "-sb", default=False, action="store_true",
        help="If the sRNA database which you used are the same format as our default sRNA database, "
        "you can run sRNA_blast_stat for do statistics of the result of sRNA blast."
        "If your format is not the same as our default database, please don't turn it on. "
        "Out default format of header is ID|strain|srna_name")
    sRNA_parser.add_argument(
        "--tex_notex_libs", "-tl", default=None, nargs="+",
        help="library name of tex and notex library. The format is: "
        "wig_file_name:tex_treat_or_not(tex or notex):condition_id(integer):replicate_id(alphabet):strand(+ or -).")
    sRNA_parser.add_argument(
        "--frag_libs", "-fl", default=None, nargs="+",
        help="library name of fragment library. The format is: "
        "wig_file_name:fragmented(frag):condition_id(integer):replicate_id(alphabet):strand(+ or -).")
    sRNA_parser.add_argument(
        "--tex_notex", "-te", default=2, type=int,
        help="For tex +/- library, sRNA candidates should be detected by both or just one.(1/2) Default is 2.")
    sRNA_parser.add_argument(
        "--replicates_tex", "-rt", default=None, 
        help="The sRNA of tex +/- library should be detected more than this number of replicates.")
    sRNA_parser.add_argument(
        "--replicates_frag", "-rf", default=None, 
        help="The sRNA of fragmented library should be detected more than this number of replicates.")
    sRNA_parser.add_argument(
        "--table_best", "-tb", default=False, action="store_true",
        help="The output table of sRNA candidates only print the best track. Default is False")
    sRNA_parser.add_argument(
        "--decrease_intergenic","-di", default=0.5, type=float,
        help="If the intergenic region is longer than the max_length, it will based on coverage to check the sRNA candidates. "
             "If the ratio of lowest coverage of intergenic region and the highest coverage of intergenic region is smaller than this number, "
             "it will consider the the point of lowest coverage to be end of sRNA. If the length of sRNA candidate is properly, "
             "it also assign the transcript to be one of sRNA candidates. Default is 0.5.")
    sRNA_parser.add_argument(
        "--decrease_utr","-du", default=0.5, type=float,
        help="If the kind of utr derived is 5'UTR, you have to consider the end of it's end."
             "If the ratio of lowest coverage of it and the highest coverage of it is smaller than this number, "
             "it will consider the the point of lowest coverage to be end of sRNA. If the length of sRNA candidate is properly, "
             "it also assign the transcript to be one of sRNA candidates. Default is 0.5.")
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
        "--blast_e_nr", "-en", default=0.0001, type=float,
        help="The cutoff of blast e value for nr alignment. ")
    sRNA_parser.add_argument(
        "--blast_e_srna", "-es", default=0.0001, type=float,
        help="The cutoff of blast e value for sRNA alignment. ")
    sRNA_parser.add_argument(
        "--sORF", "-O", default=None,
        help="If you want to compare sORF and sRNA, please assign the path of sORF gff folder.")
    sRNA_parser.add_argument(
        "--best_with_all_sRNAhit", "-ba", default=False, action="store_true",
        help="When you want to generate the files which store the best sRNA candidates, "
             "it should include all the sRNA candidates which can find the homology from blast sRNA database no matter other information(ex. TSS, blast in nr...)."
             "Please turn it on. Or it will just select the the best candidates based on all filter conditions. Default is False.")
    sRNA_parser.add_argument(
        "--best_without_sORF_candidate", "-bs", default=False, action="store_true",
        help="If you want to generate the files which store the best sRNA candidates excluded all the sRNA candidates which also can be detected by sORF file."
             "Please turn it on. Or it will just select the the best candidates based on all filter conditions. Default is False.")
    sRNA_parser.set_defaults(func=run_sRNA_detection)
    # Parameters of small ORF
    sORF_parser = subparsers.add_parser(
        "sorf", help="Run sORF detection to detect sORF candidates which has expression.")
    sORF_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    sORF_parser.add_argument(
        "--UTR_derived_sORF", "-u", default=False, action="store_true",
        help="If you want to detect UTR derived sORF, please turn it on. Default is False.")
    sORF_parser.add_argument(
        "--transcript_assembly_folder", "-a", default="ANNOgesic/output/transcriptome_assembly/gffs",
        help="The path of transcriptome assembly folder.")
    sORF_parser.add_argument(
        "--annotation_folder", "-g", default="ANNOgesic/output/target/annotation",
        help="The path of annotation gff folder.")
    sORF_parser.add_argument(
        "--TSS_folder", "-t", default=None,
        help="If you want to import TSS information, please assign the path of gff folder of TSS.")
    sORF_parser.add_argument(
        "--utr_length", "-ul", default=300, type=int,
        help="If you want to import TSS information, please assign the utr length for comparing TSS and sORF. "
        "The default number is 300.")
    sORF_parser.add_argument(
        "--min_length", "-lm",default=30,
        help="Please assign the minium length of sORF. "
        "It will classify sORF candidates based on the value. Default is 30.")
    sORF_parser.add_argument(
        "--max_length", "-lM",default=500,
        help="Please assign the maxium length of sORF. "
        "It will classify sORF candidates based on the value. Default is 500.")
    sORF_parser.add_argument(
        "--tex_wig_folder", "-tw", default=None,
        help="The path of tex+/- wig folder.")
    sORF_parser.add_argument(
        "--frag_wig_folder", "-fw", default=None,
        help="The path of fragment wig folder.")
    sORF_parser.add_argument(
        "--cutoff_intergenic_coverage", "-ci",default=5, type=float,
        help="The cutoff of minimal coverage of intergenic sORF candidates.")
    sORF_parser.add_argument(
        "--cutoff_5utr_coverage", "-cu5",default="median",
        help="The cutoff of minimal coverage of 5'UTR derived sORF candidates. You also can assign median or mean. Default is median.")
    sORF_parser.add_argument(
        "--cutoff_3utr_coverage", "-cu3",default="median",
        help="The cutoff of minimal coverage of 3'UTR derived sORF candidates. You also can assign median or mean. Default is median.")
    sORF_parser.add_argument(
        "--cutoff_interCDS_coverage", "-cuf",default="median",
        help="The cutoff of minimal coverage of interCDS derived sORF candidates. You also can assign median or mean. Default is median.")
    sORF_parser.add_argument(
        "--cutoff_background", "-cub",default=5, type=float,
        help="The cutoff of minimal coverage of all sORF candidates. Default is 5.")
    sORF_parser.add_argument(
        "--fasta_folder", "-f", default="ANNOgesic/output/target/fasta",
        help="The folder of fasta file. ")
    sORF_parser.add_argument(
        "--tex_notex_libs", "-tl", default=None, nargs="+", 
        help="Library name of tex and notex library. The format is: "
        "wig_file_name:tex_treat_or_not(tex or notex):condition_id(integer):replicate_id(alphabet):strand(+ or -).")
    sORF_parser.add_argument(
        "--frag_libs", "-fl", default=None, nargs="+", 
        help="Library name of fragment library The format is: "
        "wig_file_name:fragmented(frag):condition_id(integer):replicate_id(alphabet):strand(+ or -)..")
    sORF_parser.add_argument(
        "--tex_notex", "-te", default=2, type=int,
        help="For tex +/- library, sORF candidates should be detected by both or just one.(1/2) Default is 2.")
    sORF_parser.add_argument(
        "--replicates_tex", "-rt", default=None, 
        help="The sORF of tex +/- library should be detected more than this number of replicates.")
    sORF_parser.add_argument(
        "--replicates_frag", "-rf", default=None, 
        help="The sORF of fragmented library should be detected more than this number of replicates.")
    sORF_parser.add_argument(
        "--table_best", "-tb", default=False, action="store_true",
        help="The output table of sORF candidates only print the best track. Default is False.")
    sORF_parser.add_argument(
        "--sRNA_folder", "-s", default=None,
        help="If you want to compare sORF and sRNA, please assign the path of sORF gff folder.")
    sORF_parser.add_argument(
        "--start_coden", "-ac", default=["ATG"], nargs="+", 
        help="What kinds of start coden ATG/GTG/TTG you want to use.")
    sORF_parser.add_argument(
        "--stop_coden", "-oc", default=["TTA", "TAG", "TGA"], nargs="+", 
        help="What kinds of stop coden TTA/TAG/TGA you want to use.")
    sORF_parser.add_argument(
        "--condition_best", "-c", default="TSS",
        help="For generating the result of best sORF, please assign which information you want to consider (TSS/sRNA/both). default is TSS")
    sORF_parser.set_defaults(func=run_sORF_detection)
    # Parameters of promoter detection
    promoter_parser = subparsers.add_parser(
        "promoter", help="Run MEME to dicover promoter.")
    promoter_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    promoter_parser.add_argument(
        "--MEME_path", default="meme",
        help="path of MEME")
    promoter_parser.add_argument(
        "--fasta_folder", "-f", default="ANNOgesic/output/target/fasta",
        help="Please assign the folder of gemonic fasta file. ")
    promoter_parser.add_argument(
        "--TSS_folder", "-t", default="ANNOgesic/output/TSS/gffs",
        help="The folder of TSS gff file.")
    promoter_parser.add_argument(
        "--num_motif", "-n", default=10,
        help="How many of motifs you want to produce.")
    promoter_parser.add_argument(
        "--motif_width", "-w", nargs="+",
        help="Motif length - it will refer the value to find the motif. "
        "if you want to detect a range of width, you can insert \"-\" between two values. "
        "for example, 50 2-10. It means the range of width which you want to detect is 50 and within 2 to 10. ")
    promoter_parser.add_argument(
        "--parallel", "-p", default=4,
        help="How many process you want to use to run paralle.")
    promoter_parser.add_argument(
        "--TSS_source", "-s", default=True, action="store_false",
        help="If you generate TSS from other method, please turn it on.")
    promoter_parser.add_argument(
        "--tex_libs", "-tl", default=None, nargs="+", 
        help="Library name of tex+/- library. If your TSS is not from ANNOgesic, please assign the libs of tex+/- too."
        "The format is: "
        "wig_file_name:tex_treat_or_not(tex or notex):condition_id(integer):replicate_id(alphabet):strand(+ or -).")
    promoter_parser.add_argument(
        "--tex_wig_path", "-tw", default=None,
        help="The path of tex+/- wig folder. If your TSS is not from ANNOgesic, please assign the wig path too.")
    promoter_parser.add_argument(
        "--annotation_folder", "-g", default=None,
        help="The path of annotation gff folder. If your TSS is not from ANNOgesic, please assign the annotation gff path too.")
    promoter_parser.add_argument(
        "--combine_all", "-c", default=False, action="store_true",
        help="If you want to combine all TSS in TSS output folder to generate a overall promoter motif, please turn it on. Default is False.")
    promoter_parser.set_defaults(func=run_MEME)
    # Parameters of operon detection
    operon_parser = subparsers.add_parser(
        "operon", help="Detect operon and combine features together.")
    operon_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    operon_parser.add_argument(
        "--TSS_folder", "-t", default="ANNOgesic/output/TSS/gffs",
        help="The path of TSS gff folder.")
    operon_parser.add_argument(
        "--annotation_folder", "-g", default="ANNOgesic/output/target/annotation",
        help="The path of annotation gff folder.")
    operon_parser.add_argument(
        "--transcript_folder", "-a", default="ANNOgesic/output/transcriptome_assembly/gffs",
        help="The path of transcript assembly gff folder.")
    operon_parser.add_argument(
        "--UTR5_folder", "-u5", default="ANNOgesic/output/UTR/5UTR/gffs",
        help="The path of 5'UTR gff folder.")
    operon_parser.add_argument(
        "--UTR3_folder", "-u3", default="ANNOgesic/output/UTR/3UTR/gffs",
        help="The path of 3'UTR gff folder.")
    operon_parser.add_argument(
        "--term_folder", "-e", default=None,
        help="If you want to import the information of terminator, please assign the path of terminator gff folder.")
    operon_parser.add_argument(
        "--TSS_fuzzy", "-tf", default=5, type=int,
        help="The fuzzy for comparing TSS and transcript assembly. "
        "The default number is 5.")
    operon_parser.add_argument(
        "--term_fuzzy", "-ef", default=30, type=int,
        help="The fuzzy for comparing terminator and transcript assembly. "
        "The default number is 30.")
    operon_parser.add_argument(
        "--min_length", "-l", default=20, type=int,
        help="The minimum length of operon. "
        "The default number is 20.")
    operon_parser.add_argument(
        "--statistics", "-s", default=False, action="store_true",
        help="Doing statistics for Operon analysis. Default is False."
        "The name of statistics file is - stat_operon_$STRAIN_NAME.csv.")
    operon_parser.add_argument(
        "--combine_gff", "-c", default=False, action="store_true",
        help="Convert the operon and all features you assigned to one gff file. Default is False.")
    operon_parser.set_defaults(func=run_operon)
    # Parameters of CircRNA detection
    circrna_parser = subparsers.add_parser(
        "circrna", help="Detect circular RNA.")
    circrna_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    circrna_parser.add_argument(
        "--segemehl_folder", "-sg", default="",
        help="Please assign the folder of segemehl.")
    circrna_parser.add_argument(
        "--samtools_path", "-st", default="samtools",
        help="Please assign the path of samtools.")
    circrna_parser.add_argument(
        "--align", "-a", default=False, action="store_true",
        help="Using segemehl to align read (included splice detection). "
        "If you already usd segemehl with -S to align your reads, "
        "you can skip this step, don't need to turn it on. "
        "Please be attention, it only use default parameters of segemehl to align your reads. "
        "Moreover, it will align all read files in ANNOgesic/input/reads. "
        "If you want to run some specific functions of segemehl, "
        "please run it seperately.")
    circrna_parser.add_argument(
        "--normal_bam_path", "-nb", default=None,
        help="If you already has Bam files, Please assign the normal Bam path or fragmented_Bam_path.")
    circrna_parser.add_argument(
        "--fragmented_bam_path", "-fb", default=None,
        help="If you already has Bam files, Please assign the fragmented Bam path or normal Bam path.")
    circrna_parser.add_argument(
        "--process", "-p", default=10, type=int,
        help="How many parallels processes for --align.  ")
    circrna_parser.add_argument(
        "--fasta_path", "-f", default="ANNOgesic/output/target/fasta",
        help="The folder of genome fasta. ")
    circrna_parser.add_argument(
        "--annotation_path", "-g", default="ANNOgesic/output/target/annotation",
        help="The folder of annotation gff files.")
    circrna_parser.add_argument(
        "--convert_to_gff", "-cg", default=False, action="store_true",
        help="If you want to convert circRNA candidates to gff file, please turnn it on.")
    circrna_parser.add_argument(
        "--support_reads", "-s", default=5, type=int,
        help="If you want to convert circRNA candidates to gff file, please also assign the cut off of supported reads. Default is 5.")
    circrna_parser.add_argument(
        "--start_ratio", "-sr", default=0.25, type=float, 
        help="The ratio of (read support circ / all read) at starting point. The ratio of candidates should higher than this cutoff. Default is 0.25.")
    circrna_parser.add_argument(
        "--end_ratio", "-er", default=0.25, type=float, 
        help="The ratio of (read support circ / all read) at end point. The ratio of candidates should higher than this cutoff. Default is 0.25")
    circrna_parser.set_defaults(func=run_circrna)
    # Parameters of Go term
    goterm_parser = subparsers.add_parser(
        "go_term", help="Extract and find Go terms.")
    goterm_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    goterm_parser.add_argument(
        "--annotation_path", "-g", default="ANNOgesic/output/target/annotation",
        help="The path of annotation gff folder. ")
    goterm_parser.add_argument(
        "--UniProt_id", "-u", default="ANNOgesic/input/database/idmapping_selected.tab",
        help="The path of UniProt ID mapping database. "
        "default is ANNOgesic/input/database/idmapping_selected.tab")
    goterm_parser.add_argument(
        "--go_obo", "-go", default="ANNOgesic/input/database/go.obo",
        help="The path of go.obo. "
        "Default is ANNOgesic/input/database/go.obo")
    goterm_parser.add_argument(
        "--goslim_obo", "-gs", default="ANNOgesic/input/database/goslim.obo",
        help="The path of goslim.obo. "
        "Default is ANNOgesic/input/database/goslim.obo")
    goterm_parser.set_defaults(func=run_goterm)
    # Parameters of sRNA target prediction
    starget_parser = subparsers.add_parser(
        "srna_target", help="sRNA target prediction.")
    starget_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    starget_parser.add_argument(
        "--Vienna_folder", default="",
        help="Please assign the folder of Vienna package. It should include RNAplfold, RNAup and RNAplex.")
    starget_parser.add_argument(
        "--annotation_path", "-g", default="ANNOgesic/output/target/annotation",
        help="The path of annotation gff folder. ")
    starget_parser.add_argument(
        "--fasta_path", "-f", default="ANNOgesic/output/target/fasta",
        help="The path of genome fasta folder. ")
    starget_parser.add_argument(
        "--sRNA_path", "-r", default="ANNOgesic/output/sRNA/gffs",
        help="The path of sRNA gff folder. ")
    starget_parser.add_argument(
        "--query_sRNA", "-q", nargs="+", default="all",
        help="Please assign the query sRNA. If you want to compute all sRNA in gff file, please keyin 'all'."
        "the input format should be like, $STRAIN:$STRAND:$START:$END."
        "For example, NC_007795:+:200:534 NC_007795:-:6767:6900")
    starget_parser.add_argument(
        "--program", "-p", default="both",
        help="Using RNAplex, RNAup or both. Default is both. ")
    starget_parser.add_argument(
        "--interaction_length", "-i", default=30, type=int,
        help="Maximal length of an interaction. Default is 30. ")
    starget_parser.add_argument(
        "--window_size_target", "-wt", default=240, type=int,
        help="Only work when --program is RNAplex or both. "
        "Average the pair probabilities over windows of given size for RNAplex. "
        "Default is 240. ")
    starget_parser.add_argument(
        "--span_target", "-st", default=160, type=int,
        help="Only work when --program is RNAplex or both. "
        "Set the maximum allowed separation of a base pair to span for RNAplex. "
        "Default is 160. ")
    starget_parser.add_argument(
        "--window_size_srna", "-ws", default=30, type=int,
        help="Only work when --program is RNAplex or both. "
        "Average the pair probabilities over windows of given size for RNAplex. "
        "Default is 30. ")
    starget_parser.add_argument(
        "--span_srna", "-ss", default=30, type=int,
        help="Only work when --program is RNAplex or both. "
        "Set the maximum allowed separation of a base pair to span for RNAplex. "
        "Default is 30. ")
    starget_parser.add_argument(
        "--unstructured_region_RNAplex_target", "-ut", default=30, type=int,
        help="Only work when --program is RNAplex or both. "
        "Compute the mean probability that regions of "
        "length 1 to a given length are unpaired for RNAplex. Default is 30. ")
    starget_parser.add_argument(
        "--unstructured_region_RNAplex_srna", "-us", default=30, type=int,
        help="Only work when --program is RNAplex or both. "
        "Compute the mean probability that regions of "
        "length 1 to a given length are unpaired for RNAplex. Default is 30. ")
    starget_parser.add_argument(
        "--unstructured_region_RNAup", "-uu", default=30, type=int,
        help="Only work when --program is RNAup or both. "
        "Compute the mean probability that regions of "
        "length 1 to a given length are unpaired for RNAplex. Default is 30. ")
    starget_parser.add_argument(
        "--energy_threshold", "-e", default=-8, type=float,
        help="Only work when --program is RNAplex or both. "
        "Minimal energy for a duplex to be returned for RNAplex. "
        "Default is -8. ")
    starget_parser.add_argument(
        "--duplex_distance", "-d", default=20, type=int,
        help="Only work when --program is RNAplex or both. "
        "Distance between target 3' ends of two consecutive duplexes for RNAplex. "
        "Default is 20. ")
    starget_parser.add_argument(
        "--top", "-t", default=20, type=int,
        help="The output file only include top one(default is 20). ")
    starget_parser.add_argument(
        "--process_rnaplex", "-pp", default=5, type=int,
        help="How many parallel processes for running RNAplex prediction. ")
    starget_parser.add_argument(
        "--process_rnaup", "-pu", default=20, type=int,
        help="How many parallel processes for running RNAup prediction. ")
    starget_parser.add_argument(
        "--continue_rnaup", "-cr", default=False, action="store_true",
        help="RNAup will take a long time for running if you want to compute a lot of sRNA. "
        "If the process crush, you can turn it on. This flag will continue running RNAup based on your previous running.")
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
        help="Please assign the program for detecting SNP of transcript: "
        "1: calculate with BAQ, 2: calculate without BAQ, 3: calculate with extend BAQ. "
        "You can assign more than 1 program. For example: 1 2 3")
    snp_parser.add_argument(
        "--fasta_path", "-f", default="ANNOgesic/output/target/fasta",
        help="The path of fasta folder. ")
    snp_parser.add_argument(
        "--tex_bam_path", "-tw", default=None,
        help="The path of tex+/- wig folder. If you want to use tex treated and untreated bam files, "
             "please assign the path. Or it will not combine the bam files")
    snp_parser.add_argument(
        "--frag_bam_path", "-fw", default=None,
        help="The path of fragmented wig folder. If you want to use fragmented bam files, "
             "please assign the path. Or it will not combine the bam files")
    snp_parser.add_argument(
        "--quality", "-q", default=20, type=int,
        help="The min quality which consider a real snp. Default is 20")
    snp_parser.add_argument(
        "--read_depth", "-d", default=None, type=int,
        help="The minimum read depth, below to it will be excluded. default is 5 * number of BAM files,"
        "if the cutoff higher than 40, it will use 40.")
    snp_parser.add_argument(
        "--indel_fraction", "-imf", default=0.5, type=float,
        help="The fraction of maximum read depth, which support insertion of deletion. Default is 0.5")
    snp_parser.set_defaults(func=SNP)
    # Parameters of protein-protein interaction network
    ppi_parser = subparsers.add_parser(
        "ppi_network", help="Generate protein-protein interaction with literature supported.")
    ppi_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    ppi_parser.add_argument(
        "--ptt_path", "-p", default="ANNOgesic/output/target/annotation",
        help="The path of .ptt annotation folder. ")
    ppi_parser.add_argument(
        "--proteinID_strains", "-s", nargs="+",
        help="This is for assigning protein Id which you want to predict. "
             "In order to retrieve the data from STRING and Pubmed, you also have to assign the similar reference. "
             "For example, if you want to run all proteins in Staphylococcus aureus HG003, "
             "you can assign the it like "
             "Staphylococcus_aureus_HG003.ptt:Staphylococcus_aureus_HG003:\"Staphylococcus aureus 8325\":\"Staphylococcus aureus\". "
             "or Staphylococcus_aureus_HG003.ptt:Staphylococcus_aureus_HG003:\"93061\":\"Staphylococcus aureus\". "
             "or Staphylococcus_aureus_HG003.ptt:Staphylococcus_aureus_HG003:\"Staphylococcus aureus NCTC 8325\":\"Staphylococcus aureus\". "
             "(ptt_filename:strain_ptt:STRING_name:Pubmed_name). "
             "First one is the ptt file name. Second one is the seq name in ptt files. "
             "Third one is for STRING database, and the fourth one is for Pubmed. "
             "Of course, you can run the script for several strains at the same time. "
             "Before running it, please check the species file which located in ANNOgesic/input/database ."
             "If you didn't download the file, please download it. "
             "You can use taxon_id, STRING_name_compact or official_name_NCBI to represent STRING_name."
             "BE CAREFUL, if the name which you assigned has spaces, please put \"\" at two ends. "
             "For the name of Pubmed, you can assign the name not so specific. "
             "If you assign a specific name, it may not be able to find the related literatures.")
    ppi_parser.add_argument(
        "--without_strain_pubmed", "-n", default=False, action="store_true",
        help="If you want to retrieve pubmed without assign any strains, please turn it on. Default is False.")
    ppi_parser.add_argument(
        "--species_STRING", "-d",
        help="Please assign the path of species file of STRING.")
    ppi_parser.add_argument(
        "--score", "-ps", default=0.0, type=float,
        help="Please assign the cutoff of score. The value is from -1 to 1")
    ppi_parser.add_argument(
        "--node_size", "-ns", default=4000, type=int,
        help="Please size of the nodes in figure, default is 4000.")
    ppi_parser.add_argument(
        "--query", "-q", default="all", nargs="+",
        help="Please assign the query protein here. The format is $STRAIN:$START_POINT:$END_POINT:$STRAND."
        "For example, NC_007795:345:456:+ NC_007795:2000:3211:-. If you want to compute all protein, just type all."
        "The default is all.")
    ppi_parser.set_defaults(func=PPI)
    # Parameters of subcellular localization
    sub_local_parser = subparsers.add_parser(
        "subcellular_localization", help="prediction of subcellular localization of genomic CDS.")
    sub_local_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    sub_local_parser.add_argument(
        "--Psortb_path", default="psort",
        help="If you want to assign the path of Psortb, please assign here.")
    sub_local_parser.add_argument(
        "--gff_path", "-g", default="ANNOgesic/output/target/annotation",
        help="The path of annotation gff folder. ")
    sub_local_parser.add_argument(
        "--fasta_path", "-f", default="ANNOgesic/output/target/fasta",
        help="The path of fasta folder. ")
    sub_local_parser.add_argument(
        "--bacteria_type", "-b",
        help="Is Gram-positive or Gram-negative. Please assign 'positive' or 'negative'. ")
    sub_local_parser.add_argument(
        "--difference_multi", "-d", default=0.5,
        help="If the protein may have multiple location, "
        "it will calculte the difference of scores(psortb) between best one and others. "
        "If the difference is within this value, it will print it out, too. Default is 0.5. "
        "The maximum value is 10. ")
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
        "--infernal_path", "-if", default="",
        help="Please assign the folder of Infernal(where is cmscan and cmsearch located).")
    ribos_parser.add_argument(
        "--riboswitch_ID", "-i",
        help="The path of the riboswitch ID of Rfam. ")
    ribos_parser.add_argument(
        "--gff_path", "-g", default="ANNOgesic/output/target/annotation",
        help="The path of annotation gff folder. ")
    ribos_parser.add_argument(
        "--fasta_path", "-f", default="ANNOgesic/output/target/fasta",
        help="The path of fasta folder. ")
    ribos_parser.add_argument(
        "--Rfam", "-R",
        help="The path of Rfam CM database. ")
    ribos_parser.add_argument(
        "--re_scan", "-r", default=False, action="store_true",
        help="Based on the results of first scaning, it will modify the input of sequence and re scan again. "
        "Default is False.")
    ribos_parser.add_argument(
        "--e_value", "-e", default=0.001, type=float,
        help="The cutoff of e value. Default si 0.001.")
    ribos_parser.add_argument(
        "--output_all", "-o", default=False, action="store_true",
        help="One sequence may fit multiple riboswitches. If you want to output all of them, please turn it on. "
        "Or it will only print the best one.")
    ribos_parser.add_argument(
        "--fuzzy", "-z", default=10, type=int,
        help="It will extend some nts of 3' ans 5' end. Default is 10.")
    ribos_parser.set_defaults(func=Ribos)
    # Parameters of generate screenshots
    screen_parser = subparsers.add_parser(
        "screenshot", help="Generate screenshot for selected feature.")
    screen_parser.add_argument(
        "project_path", default=".", nargs="?",
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    screen_parser.add_argument(
        "--main_gff", "-mg",
        help="Screenshot will based on the position of main_gff file to generate screenshot. ")
    screen_parser.add_argument(
        "--side_gffs", "-sg", default=None, nargs="+",
        help="If you have more than one gff want to plot together, please assign here. ")
    screen_parser.add_argument(
        "--fasta", "-f", default="ANNOgesic/output/target/fasta",
        help="The path of genome fasta folder. ")
    screen_parser.add_argument(
        "--frag_wig_folder", "-fw", default=None,
        help="If you want to include the information of fragmented wig file, please assign the folder. ")
    screen_parser.add_argument(
        "--tex_wig_folder", "-tw", default=None,
        help="If you want to include the information of tex treated wig file, please assign the folder. ")
    screen_parser.add_argument(
        "--height", "-he", default=1500, type=int,
        help="You can assign the height of screenshot. Default is 1500.")
    screen_parser.add_argument(
        "--tex_libs", "-tl", default=None, nargs="+",
        help="If you want to include the tex treated wig file, please also assign proper format here. The format is: "
        "wig_file_name:tex_treat_or_not(tex or notex):condition_id(integer):replicate_id(alphabet):strand(+ or -).")
    screen_parser.add_argument(
        "--frag_libs", "-fl", default=None, nargs="+",
        help="If you want to include the fragmented wig file, please also assign proper format here. The format is: "
        "wig_file_name:fragmented(frag):condition_id(integer):replicate_id(alphabet):strand(+ or -)..")
    screen_parser.add_argument(
        "--present", "-p", default="expand",
        help="Which type you want to present in the screen shot. expand/collapse/squish. ")
    screen_parser.add_argument(
        "--output_folder", "-o",
        help="Please assign the output folder. If the folder does not exist, it will generate automatically. ")
    screen_parser.set_defaults(func=Screen)

    args = parser.parse_args()
    
    if args.version is True:
        print("ANNOgesic version " + __version__)
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

def run_expression(controller):
    controller.expression()

def multiparser(controller):
    controller.multiparser()

def run_TSSpredator(controller):
    controller.tsspredator()

def optimize_TSSpredator(controller):
    controller.optimize()

def color_png(controller):
    controller.color()

def run_Terminator(controller):
    controller.terminator()

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
