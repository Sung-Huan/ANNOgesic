import os
import sys

from transaplib.projectcreator import ProjectCreator
from transaplib.paths import Paths
from transaplib.get_input import Get_file
from transaplib.converter import Converter
from transaplib.get_target_fasta import Target_fasta
from transaplib.ratt import RATT
from transaplib.tsspredator import TSSpredator
from transaplib.optimize import Optimize_TSS
from transaplib.color_png import Color_PNG
from transaplib.terminator import Terminator
from transaplib.transcript import Transcript_Assembly
from transaplib.utr import UTR_detection
from transaplib.srna import sRNA_detection
from transaplib.sorf import sORF_detection
from transaplib.meme import MEME
from transaplib.operon import Operon_detection
from transaplib.circrna import CircRNA_detection
from transaplib.goterm import Go_term_finding
from transaplib.srna_target import sRNA_target_prediction
from transaplib.snp import SNP_calling
from transaplib.ppi import PPI_network
from transaplib.sublocal import Sub_Local
from transaplib.ribos import Ribos
from transaplib.screen import Screen

project_creator = ProjectCreator()
class Controller(object):

    """Manage the actions of the subcommands.

    The Controller take care of providing the argumentes like path
    names and the parallel processing of tasks.

    """
    def __init__(self, args):
        """Create an instance."""
        self._args = args
        self._paths = Paths(args.project_path)
    def create_project(self, version):
        """Create a new project."""
        project_creator.create_root_folder(self._args.project_path)
        project_creator.create_subfolders(self._paths.required_folders("root"))
        project_creator.create_version_file(self._paths.version_path, version)
        sys.stdout.write("Created folder \"%s\" and required subfolders.\n" % (
            self._args.project_path))
        sys.stdout.write("Please copy tex and no tex wig files into \"%s\.\n" % (
            self._paths.tex_folder))
        sys.stdout.write("Please copy mutation table files into \"%s\.\n" % (
            self._paths.mutation_table_folder))
        sys.stdout.write("Please copy fragment wig files into \"%s\.\n" % (
            self._paths.frag_folder))
        sys.stdout.write("Please copy read files into \"%s\.\n" % (
            self._paths.read_folder))
        sys.stdout.write("Please copy BAM files which map to reference genome(if you start from them) into \"%s\.\n" % (
            self._paths.bam_ref_folder))
        sys.stdout.write("Please copy BAM files which map to target genome(if you already had them) into \"%s\.\n" % (
            self._paths.bam_tar_folder))
        sys.stdout.write("Please download sRNAdatabse files(i.e. BSRD) into \"%s\.\n" % (
            self._paths.database_folder))
        sys.stdout.write("Please copy the riboswitch information of Rfam into \"%s\.\n" % (
            self._paths.ribos_folder))
    def get_input(self):
        """Download required files from website."""
        if self._args.ref_gff is True:
            Get_file(self._args.FTP_path, self._paths.ref_annotation_folder, "gff")
        if self._args.ref_fasta is True:
            Get_file(self._args.FTP_path, self._paths.ref_fasta_folder, "fna")
        if self._args.ref_gbk is True:
            Get_file(self._args.FTP_path, self._paths.ref_annotation_folder, "gbk")
        if self._args.convert_embl is True:
            annotation_files = os.listdir(self._paths.ref_annotation_folder)
            if len(annotation_files) == 0:
                sys.stdout.write("No gbk files!!\n")
            else:
                Converter().Convert_gbk2embl(self._paths.ref_annotation_folder, annotation_files)
    def get_target_fasta(self):
        """Get target fasta"""
        project_creator.create_subfolders(self._paths.required_folders("get_target_fasta"))
        target = Target_fasta()
        target.get_target_fasta(self._args.mutation_table,
                self._paths.tar_fasta_folder, self._args.ref_fasta_folder,
                self._args.output_format)
    def ratt(self):
        """Run RATT to transfer annotation file from reference to target."""
        project_creator.create_subfolders(self._paths.required_folders("annotation_transfer"))
        ratt = RATT()
        ratt.run_RATT(self._args.RATT_path, self._args.PAGIT_folder, 
                      self._args.target_ref_comparison, 
                      self._args.element, self._args.transfer_type, 
                      self._args.ref_embl, self._args.target_fasta, 
                      self._args.ref_fasta, self._paths.ratt_folder, 
                      self._args.convert_to_gff_rnt_ptt, self._args.combine_all_gff,
                      self._paths.tar_annotation_folder)

    def tsspredator(self):
        """Run TSSpredator for predicting TSS candidates."""
        if self._args.compute_program == "TSS":
            project_creator.create_subfolders(self._paths.required_folders("TSS"))
            out_folder = self._paths.tsspredator_folder
        elif self._args.compute_program == "processing_site":
            out_folder = self._paths.processing_site_folder
            project_creator.create_subfolders(self._paths.required_folders("processing"))
        else:
            print("Error:No such program!!!!")
            sys.exit()
        tsspredator = TSSpredator()
        tsspredator.run_TSSpredator(
                      self._args.TSSpredator_path, self._args.compute_program,
                      self._paths.tsspredator_input_folder, 
                      self._args.fasta_folder, self._args.annotation_folder,
                      self._args.wig_folder, self._args.lib, 
                      self._args.output_prefix,
                      self._args.height, self._args.height_reduction,
                      self._args.factor, self._args.factor_reduction,
                      self._args.base_height, self._args.replicate_match,
                      out_folder, self._args.project_path, self._args.statistics,
                      self._args.validate_gene, self._args.merge_manual,
                      self._args.compare_transcript_assembly, self._args.fuzzy,
                      self._args.utr_length, self._args.cluster, self._args.length)
    
    def optimize(self):
        """opimize TSSpredator"""
        if self._args.program == "TSS":
            output_folder = self._paths.tsspredator_base_folder
        elif self._args.program == "Processing_site":
            output_folder = self._paths.processing_base_folder
        else:
            print("Error:No such program!!!!")
            sys.exit()
        Optimize_TSS(self._args.TSSpredator_path, self._args.fasta_file,
                     self._args.annotation_file, self._args.wig_folder,
                     self._args.manual, output_folder, 
                     self._args.strain_name, self._args.max_height, 
                     self._args.max_height_reduction, self._args.max_factor,
                     self._args.max_factor_reduction, self._args.max_base_height, 
                     self._args.utr_length,self._args.lib,
                     self._args.output_prefix, self._args.cluster,
                     self._args.length, self._args.core, 
                     self._args.program, self._args.replicate_match,
                     self._args.steps)

    def color(self):
        """color the screenshots"""
        color = Color_PNG()
        color.Generate_color_png(
                self._args.bin_path, self._args.track_number, 
                self._args.screenshot_type, self._paths.output_folder)

    def transtermhp(self):
        """Run TransTermHP for detecting terminators."""
        terminator = Terminator()
        if self._args.TransTermHP_folder is None:
            print("Please assign the folder where you install TransTermHP.")
        terminator.run_Terminator(
                self._args.TransTermHP_folder, self._args.RNAfold_path,
                self._paths.transterm_folder, self._args.fasta_folder, 
                self._args.annotation_folder, self._args.transcript_folder, 
                self._args.sRNA, self._args.statistics, 
                self._args.tex_wig_folder, self._args.frag_wig_folder, 
                self._args.decrease, self._args.highest_coverage,
                self._args.fuzzy, self._paths.transtermhp_folder,
                self._args.tex_notex_libs, self._args.frag_libs,
                self._args.tex_notex, self._args.replicates,
                self._args.table_best)

    def transcript(self):
        """Run Transcriptome assembly."""
        transcript = Transcript_Assembly()
        transcript.run_Transcript_Assembly(
                self._args.project_path, self._args.bin_path,
                self._args.frag_wig_path, self._args.normal_wig_path, 
                self._args.sort_annotation, self._args.tex_notex,
                self._args.length, self._args.annotation_folder, 
                self._args.height, self._args.width, 
                self._args.tolerance, self._args.supported,
                self._paths.transcript_assembly_output_folder,
                self._args.compare_TSS, self._args.compare_CDS, 
                self._args.TSS_fuzzy, self._args.replicates,
                self._args.Tex_treated_libs, self._args.fragmented_libs)
    def utr_detection(self):
        """Run UTR detection."""
        utr = UTR_detection()
        utr.run_UTR_detection(
                self._args.bin_path,
                self._args.TSS_folder, self._args.annotation_folder,
                self._args.transcript_assembly_folder,
                self._args.terminator_folder,
                self._args.terminator_fuzzy, self._paths.utr_folder, 
                self._args.TSS_source)

    def srna_detection(self):
        """sRNA_detection."""
        srna = sRNA_detection()
        srna.run_sRNA_detection(
                self._args.bin_path, self._paths.srna_folder,
                self._args.UTR_derived_sRNA,
                self._args.annotation_folder, self._args.TSS_folder,
                self._args.transcript_assembly_folder,
                self._args.TSS_fuzzy, self._args.import_info,
                self._args.tex_wig_folder, self._args.frag_wig_folder,
                self._args.processing_site_folder,
                self._args.fasta_folder, self._args.mountain_plot,
                self._args.database_format,
                self._args.sRNA_database_path, self._args.nr_database_path,
                self._args.cutoff_energy, self._args.cutoff_intergenic_coverage,
                self._args.cutoff_5utr_coverage, self._args.cutoff_3utr_coverage,
                self._args.cutoff_interCDS_coverage,
                self._args.max_length, self._args.min_length,
                self._args.tex_notex_libs, self._args.frag_libs,
                self._args.replicates, self._args.tex_notex,
                self._args.table_best, self._args.decrease_intergenic,
                self._args.decrease_utr, self._args.fuzzy_intergenic,
                self._args.fuzzy_utr, self._args.cutoff_nr_hit,
                self._args.sORF, self._args.best_with_all_sRNAhit,
                self._args.best_without_sORF_candidate)

    def sorf_detection(self):
        """sORF_detection."""
        sorf = sORF_detection()
        sorf.run_sORF_detection(
                self._args.bin_path, self._paths.sorf_folder,
                self._args.UTR_derived_sORF, self._args.transcript_assembly_folder,
                self._args.annotation_folder, self._args.TSS_folder,
                self._args.utr_length,
                self._args.min_length, self._args.max_length,
                self._args.tex_wig_folder, self._args.frag_wig_folder,
                self._args.cutoff_intergenic_coverage, self._args.cutoff_5utr_coverage, 
                self._args.cutoff_3utr_coverage, self._args.cutoff_interCDS_coverage,
                self._args.fasta_folder, self._args.tex_notex_libs,
                self._args.frag_libs, self._args.tex_notex,
                self._args.replicates, self._args.table_best,
                self._args.sRNA_folder, self._args.start_coden,
                self._args.stop_coden, self._args.condition_best)

    def meme(self):
        """promoter detectopn"""
        meme = MEME()
        meme.run_MEME(
                self._args.bin_path, self._paths.promoter_input_folder,
                self._paths.promoter_output_folder, self._args.tex_libs,
                self._args.TSS_folder, self._args.fasta_folder, 
                self._args.num_motif, self._args.motif_width, 
                self._args.parallel, self._args.TSS_source, 
                self._args.tex_wig_path, self._args.annotation_folder)

    def operon(self):
        """operon detection"""
        operon = Operon_detection()
        operon.run_Operon(
                self._args.TSS_folder, self._args.annotation_folder,
                self._args.transcript_folder, self._args.UTR5_folder,
                self._args.UTR3_folder, self._args.term_folder,
                self._args.TSS_fuzzy, self._args.term_fuzzy,
                self._args.min_length, self._args.statistics,
                self._paths.operon_output_folder, self._args.combine_gff,
                self._paths.operon_statistics_folder, self._args.bin_path)
    def circrna(self):
        """circRNA detection"""
        circ = CircRNA_detection()
        circ.run_CircRNA(
                self._args.align, self._args.process, self._args.fasta_path,
                self._args.annotation_path, self._args.bin_path,
                self._paths.circrna_output_folder, self._paths.read_folder,
                self._paths.circrna_stat_folder, self._args.convert_to_gff,
                self._args.support_reads)
    def goterm(self):
        """Go term discovery"""
        goterm = Go_term_finding()
        goterm.run_Go_term(
                self._args.bin_path, self._args.annotation_path,
                self._paths.database_folder, self._paths.goterm_output_folder,
                self._args.UniProt_id)

    def srna_target(self):
        """sRNA target prediction"""
        srnatarget = sRNA_target_prediction()
        srnatarget.run_sRNA_target_prediction(
                self._args.bin_path, self._args.annotation_path,
                self._args.fasta_path, self._args.sRNA_path,
                self._args.program, self._args.interaction_length,
                self._args.window_size_target, self._args.span_target, 
                self._args.window_size_srna, self._args.span_srna,
                self._args.unstructured_region_RNAplex_target,
                self._args.unstructured_region_RNAplex_srna,
                self._args.unstructured_region_RNAup,
                self._args.energy_threshold, self._args.duplex_distance,
                self._paths.starget_output_folder, self._args.process_rnaplex,
                self._args.process_rnaup)
    def snp(self):
        """SNP transcript detection"""
        project_creator.create_subfolders(self._paths.required_folders("snp"))
        snp = SNP_calling()
        snp.run_SNP_calling(
                self._args.samtools_path, self._args.bcftools_path,
                self._args.bam_type, 
                self._args.program, self._args.fasta_path,
                self._args.normal_bam_path, self._args.frag_bam_path,
                self._args.quality, self._args.read_depth,
                self._paths.snp_output_folder, self._args.indel_fraction)

    def ppi(self):
        """PPI network retrieve"""
        ppi = PPI_network()
        ppi.retrieve_PPI_network(
            self._args.bin_path, self._args.ptt_path,
            self._args.proteinID_strains,
            self._args.without_strain_pubmed,
            self._args.species_STRING,
            self._args.score,
            self._paths.ppi_output_folder,
            self._args.node_size)

    def sublocal(self):
        """Subcellular Localization prediction"""
        sublocal = Sub_Local()
        sublocal.run_sub_local(
            self._args.bin_path, self._args.gff_path,
            self._args.fasta_path, self._args.bacteria_type,
            self._args.merge_to_gff,
            self._paths.sublocal_output_folder)

    def ribos(self):
        """riboswitch prediction"""
        ribos = Ribos()
        ribos.run_ribos(
            self._args.bin_path, self._args.riboswitch_ID,
            self._args.gff_path, self._args.fasta_path,
            self._args.Rfam, self._paths.ribos_output_folder,
            self._args.re_scan, self._paths.database_folder,
            self._args.fuzzy)
    def screen(self):
        """generate screenshot"""
        screen = Screen()
        screen.gen_screenshot(
            self._args.bin_path, self._args.main_gff,
            self._args.side_gffs, self._args.fasta,
            self._args.frag_wig_folder, self._args.tex_wig_folder,
            self._args.height, self._args.tex_libs,
            self._args.frag_libs, self._args.present,
            self._args.output_folder)
