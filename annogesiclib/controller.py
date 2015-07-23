import os
import sys
from annogesiclib.projectcreator import ProjectCreator
from annogesiclib.paths import Paths
from annogesiclib.get_input import get_file
from annogesiclib.converter import Converter
from annogesiclib.get_target_fasta import TargetFasta
from annogesiclib.ratt import RATT
from annogesiclib.tsspredator import TSSpredator
from annogesiclib.optimize import optimize_tss
from annogesiclib.color_png import ColorPNG
from annogesiclib.terminator import Terminator
from annogesiclib.transcript import TranscriptAssembly
from annogesiclib.utr import UTRDetection
from annogesiclib.srna import sRNADetection
from annogesiclib.sorf import sORFDetection
from annogesiclib.meme import MEME
from annogesiclib.operon import OperonDetection
from annogesiclib.circrna import CircRNADetection
from annogesiclib.goterm import GoTermFinding
from annogesiclib.srna_target import sRNATargetPrediction
from annogesiclib.snp import SNPCalling
from annogesiclib.ppi import PPINetwork
from annogesiclib.sublocal import SubLocal
from annogesiclib.ribos import Ribos
from annogesiclib.screen import Screen
from annogesiclib.expression import Expression

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
        project_creator.create_subfolders(self._paths.required_folders("get_target_fasta"))
        project_creator.create_version_file(self._paths.version_path, version)
#        sys.stdout.write("Created folder \"%s\" and required subfolders.\n" % (
#            self._args.project_path))
#        sys.stdout.write("Please copy tex and no tex wig files into \"%s\.\n" % (
#            self._paths.tex_folder))
#        sys.stdout.write("Please copy mutation table files into \"%s\.\n" % (
#            self._paths.mutation_table_folder))
#        sys.stdout.write("Please copy fragment wig files into \"%s\.\n" % (
#            self._paths.frag_folder))
#        sys.stdout.write("Please copy read files into \"%s\.\n" % (
#            self._paths.read_folder))
#        sys.stdout.write("Please copy BAM files which map to reference genome(if you start from them) into \"%s\.\n" % (
#            self._paths.bam_ref_folder))
#        sys.stdout.write("Please copy BAM files which map to target genome(if you already had them) into \"%s\.\n" % (
#            self._paths.bam_tar_folder))
#        sys.stdout.write("Please download sRNAdatabse files(i.e. BSRD) into \"%s\.\n" % (
#            self._paths.database_folder))
#        sys.stdout.write("Please extract the riboswitch information of Rfam into \"%s\.\n" % (
#            self._paths.database_folder))
    def get_input(self):
        """Download required files from website."""
        if self._args.ref_gff is True:
            get_file(self._args.FTP_path, self._paths.ref_annotation_folder, "gff")
        if self._args.ref_fasta is True:
            get_file(self._args.FTP_path, self._paths.ref_fasta_folder, "fna")
        if self._args.ref_gbk is True:
            get_file(self._args.FTP_path, self._paths.ref_annotation_folder, "gbk")
        if self._args.ref_ptt is True:
            get_file(self._args.FTP_path, self._paths.ref_annotation_folder, "ptt")
        if self._args.ref_rnt is True:
            get_file(self._args.FTP_path, self._paths.ref_annotation_folder, "rnt")
        if self._args.convert_embl is True:
            annotation_files = os.listdir(self._paths.ref_annotation_folder)
            if len(annotation_files) == 0:
                sys.stdout.write("No gbk files!!\n")
            else:
                Converter().convert_gbk2embl(self._paths.ref_annotation_folder)
    def get_target_fasta(self):
        """Get target fasta"""
        project_creator.create_subfolders(self._paths.required_folders("get_target_fasta"))
        target = TargetFasta(self._paths.tar_fasta_folder, self._args.ref_fasta_folder)
        target.get_target_fasta(self._args.mutation_table,
                self._paths.tar_fasta_folder, self._args.ref_fasta_folder,
                self._args.output_format)
    def ratt(self):
        """Run RATT to transfer annotation file from reference to target."""
        project_creator.create_subfolders(self._paths.required_folders("annotation_transfer"))
        ratt = RATT(self._args.ref_gbk, self._paths.ratt_folder,
                    self._args.target_fasta, self._args.ref_fasta,
                    self._paths.tar_annotation_folder)
        ratt.annotation_transfer(self._args.RATT_path, 
                                 self._args.element, self._args.transfer_type, 
                                 self._args.ref_gbk, self._args.target_fasta, 
                                 self._args.ref_fasta, self._paths.ratt_folder, 
                                 self._args.convert_to_gff_rnt_ptt,
                                 self._paths.tar_annotation_folder, self._args.compare_pair)

    def expression(self):
        """Run gene expression analysis."""
        expression = Expression(self._args.annotation_folder)
        expression.gene_expression(self._args.tex_notex_libs, self._args.frag_libs,
                                   self._args.tex_notex, self._args.replicates_tex,
                                   self._args.replicates_frag, self._args.tex_wig_folder,
                                   self._args.frag_wig_folder, self._args.cutoff_overlap_tex,
                                   self._args.cutoff_overlap_frag, self._args.cutoff_coverage,
                                   self._args.annotation_folder, self._args.features)
    def tsspredator(self):
        """Run TSSpredator for predicting TSS candidates."""
        if self._args.compute_program.lower() == "tss":
            project_creator.create_subfolders(self._paths.required_folders("TSS"))
            out_folder = self._paths.tsspredator_folder
        elif self._args.compute_program.lower() == "processing_site":
            out_folder = self._paths.processing_site_folder
            project_creator.create_subfolders(self._paths.required_folders("processing"))
        else:
            print("Error:No such program!!!!")
            sys.exit()
        tsspredator = TSSpredator(out_folder, self._args.compare_transcript_assembly,
                                  self._args.annotation_folder, self._args.wig_folder,
                                  self._args.fasta_folder)
        tsspredator.run_tsspredator(
                      self._args.TSSpredator_path, self._args.compute_program,
                      self._paths.tsspredator_input_folder, 
                      self._args.fasta_folder, self._args.annotation_folder,
                      self._args.wig_folder, self._args.lib, 
                      self._args.output_prefix,
                      self._args.height, self._args.height_reduction,
                      self._args.factor, self._args.factor_reduction,
                      self._args.base_height, self._args.enrichment_factor, 
                      self._args.processing_factor, self._args.replicate_match,
                      out_folder, self._args.statistics,
                      self._args.validate_gene, self._args.merge_manual,
                      self._args.compare_transcript_assembly, self._args.fuzzy,
                      self._args.utr_length, self._args.cluster, self._args.length,
                      self._args.re_check_orphan, self._args.overlap_feature,
                      self._args.reference_gff_folder, self._args.remove_low_expression)
    
    def optimize(self):
        """opimize TSSpredator"""
        if self._args.program.lower() == "tss":
            project_creator.create_subfolders(self._paths.required_folders("TSS"))
            out_folder = self._paths.tsspredator_folder
        elif self._args.program.lower() == "processing_site":
            out_folder = self._paths.processing_site_folder
            project_creator.create_subfolders(self._paths.required_folders("processing"))
        else:
            print("Error:No such program!!!!")
            sys.exit()
        optimize_tss(self._args.TSSpredator_path, self._args.fasta_file,
                     self._args.annotation_file, self._args.wig_folder,
                     self._args.manual, out_folder, 
                     self._args.strain_name, self._args.max_height, 
                     self._args.max_height_reduction, self._args.max_factor,
                     self._args.max_factor_reduction, self._args.max_base_height, 
                     self._args.max_enrichment_factor, self._args.max_processing_factor,
                     self._args.utr_length, self._args.lib, 
                     self._args.output_prefix, self._args.cluster,
                     self._args.length, self._args.core, 
                     self._args.program, self._args.replicate_match,
                     self._args.steps)

    def color(self):
        """color the screenshots"""
        color = ColorPNG()
        color.generate_color_png(
                self._args.track_number, self._args.screenshot_folder, 
                self._args.ImageMagick_covert_path)

    def terminator(self):
        """Run TransTermHP for detecting terminators."""
        project_creator.create_subfolders(self._paths.required_folders("terminator"))
        terminator = Terminator(self._args.annotation_folder, self._args.fasta_folder,
                                self._args.transcript_folder, self._paths.transterm_folder,
                                self._args.sRNA)
        if self._args.TransTermHP_path is None:
            print("Please assign the folder where you install TransTermHP.")
        terminator.run_terminator(
                self._args.TransTermHP_path, self._args.expterm_path,
                self._args.RNAfold_path,
                self._paths.transterm_folder, self._args.fasta_folder, 
                self._args.annotation_folder, self._args.transcript_folder, 
                self._args.sRNA, self._args.statistics, 
                self._args.tex_wig_folder, self._args.frag_wig_folder, 
                self._args.decrease, self._args.highest_coverage,
                self._args.fuzzy_detect_coverage, self._args.fuzzy_upstream_transcript,
                self._args.fuzzy_downstream_transcript, self._args.fuzzy_upstream_cds,
                self._args.fuzzy_downstream_cds, self._paths.transtermhp_folder,
                self._args.tex_notex_libs, self._args.frag_libs,
                self._args.tex_notex, self._args.replicates_tex,
                self._args.replicates_frag, self._args.table_best)

    def transcript(self):
        """Run Transcriptome assembly."""
        project_creator.create_subfolders(self._paths.required_folders("transcript_assembly"))
        transcript = TranscriptAssembly(self._paths.transcript_assembly_output_folder)
        transcript.run_transcript_assembly(
                self._args.frag_wig_path, self._args.normal_wig_path, 
                self._args.sort_annotation, self._args.tex_notex,
                self._args.length, self._args.annotation_folder, 
                self._args.height, self._args.width, 
                self._args.tolerance, self._args.tolerance_coverage,
                self._args.replicates_tex, self._args.replicates_frag,
                self._paths.transcript_assembly_output_folder,
                self._args.compare_TSS, self._args.compare_CDS, 
                self._args.TSS_fuzzy, self._args.Tex_treated_libs, 
                self._args.fragmented_libs)
    def utr_detection(self):
        """Run UTR detection."""
        project_creator.create_subfolders(self._paths.required_folders("utr"))
        utr = UTRDetection(self._args.TSS_folder, self._args.transcript_assembly_folder,
                            self._paths.utr_folder)
        utr.run_utr_detection(
                self._args.TSS_folder, self._args.annotation_folder,
                self._args.transcript_assembly_folder,
                self._args.terminator_folder,
                self._args.terminator_fuzzy, self._paths.utr_folder, 
                self._args.TSS_source, self._args.base_5UTR)

    def srna_detection(self):
        """sRNA_detection."""
        project_creator.create_subfolders(self._paths.required_folders("srna"))
        srna = sRNADetection(self._paths.srna_folder, self._args.TSS_folder,
                              self._args.processing_site_folder, self._args.sORF, 
                              self._args.fasta_folder, self._args.transcript_assembly_folder)
        srna.run_srna_detection(
                self._args.Vienna_folder, self._args.Vienna_utils,
                self._args.blast_plus_folder,
                self._args.ps2pdf14_path, self._paths.srna_folder,
                self._args.UTR_derived_sRNA, self._args.annotation_folder, 
                self._args.TSS_folder, self._args.transcript_assembly_folder,
                self._args.TSS_intergenic_fuzzy, self._args.TSS_5UTR_fuzzy, 
                self._args.TSS_3UTR_fuzzy, self._args.TSS_interCDS_fuzzy, 
                self._args.import_info, self._args.tex_wig_folder, 
                self._args.frag_wig_folder, self._args.processing_site_folder,
                self._args.fasta_folder, self._args.mountain_plot,
                self._args.database_format,
                self._args.sRNA_database_path, self._args.nr_database_path,
                self._args.cutoff_energy, self._args.cutoff_intergenic_coverage,
                self._args.cutoff_5utr_coverage, self._args.cutoff_3utr_coverage,
                self._args.cutoff_interCDS_coverage,
                self._args.max_length, self._args.min_length,
                self._args.tex_notex_libs, self._args.frag_libs,
                self._args.replicates_tex, self._args.replicates_tex, 
                self._args.tex_notex, self._args.blast_e_nr,
                self._args.blast_e_srna,
                self._args.table_best, self._args.decrease_intergenic,
                self._args.decrease_utr, self._args.fuzzy_intergenic,
                self._args.fuzzy_utr, self._args.cutoff_nr_hit,
                self._args.sORF, self._args.best_with_all_sRNAhit,
                self._args.best_without_sORF_candidate)

    def sorf_detection(self):
        """sORF_detection."""
        project_creator.create_subfolders(self._paths.required_folders("sorf"))
        sorf = sORFDetection(self._args.TSS_folder, self._args.sRNA_folder,
                              self._paths.sorf_folder, self._args.transcript_assembly_folder,
                              self._args.fasta_folder)
        sorf.run_sorf_detection(
                self._paths.sorf_folder, self._args.UTR_derived_sORF, 
                self._args.transcript_assembly_folder,
                self._args.annotation_folder, self._args.TSS_folder,
                self._args.utr_length,
                self._args.min_length, self._args.max_length,
                self._args.tex_wig_folder, self._args.frag_wig_folder,
                self._args.cutoff_intergenic_coverage, self._args.cutoff_5utr_coverage, 
                self._args.cutoff_3utr_coverage, self._args.cutoff_interCDS_coverage,
                self._args.fasta_folder, self._args.tex_notex_libs,
                self._args.frag_libs, self._args.tex_notex,
                self._args.replicates_tex, self._args.replicates_frag,
                self._args.table_best, self._args.sRNA_folder, 
                self._args.start_coden, self._args.stop_coden, 
                self._args.condition_best, self._args.cutoff_background)

    def meme(self):
        """promoter detectopn"""
        project_creator.create_subfolders(self._paths.required_folders("promoter"))
        meme = MEME(self._args.TSS_folder, self._args.annotation_folder,
                    self._args.fasta_folder)
        meme.run_meme(
                self._args.MEME_path, self._paths.promoter_input_folder,
                self._paths.promoter_output_folder, self._args.tex_libs,
                self._args.TSS_folder, self._args.fasta_folder, 
                self._args.num_motif, self._args.motif_width, 
                self._args.parallel, self._args.TSS_source, 
                self._args.tex_wig_path, self._args.annotation_folder,
                self._args.combine_all)

    def operon(self):
        """operon detection"""
        project_creator.create_subfolders(self._paths.required_folders("operon"))
        operon = OperonDetection(self._args.TSS_folder, self._args.transcript_folder,
                                  self._args.UTR5_folder, self._args.UTR3_folder,
                                  self._paths.operon_output_folder, self._args.term_folder)
        operon.run_operon(
                self._args.TSS_folder, self._args.annotation_folder,
                self._args.transcript_folder, self._args.UTR5_folder,
                self._args.UTR3_folder, self._args.term_folder,
                self._args.TSS_fuzzy, self._args.term_fuzzy,
                self._args.min_length, self._args.statistics,
                self._paths.operon_output_folder, self._args.combine_gff,
                self._paths.operon_statistics_folder)
    def circrna(self):
        """circRNA detection"""
        project_creator.create_subfolders(self._paths.required_folders("circrna"))
        circ = CircRNADetection(self._paths.circrna_output_folder, self._args.annotation_path,
                                 self._args.fasta_path, self._args.align)
        circ.run_circrna(
                self._args.align, self._args.process, self._args.fasta_path,
                self._args.annotation_path, self._args.normal_bam_path,
                self._args.fragmented_bam_path,
                self._paths.read_folder, self._paths.circrna_stat_folder, 
                self._args.convert_to_gff, self._args.support_reads,
                self._args.segemehl_folder, self._args.samtools_path,
                self._args.start_ratio, self._args.end_ratio)
    def goterm(self):
        """Go term discovery"""
        project_creator.create_subfolders(self._paths.required_folders("go_term"))
        goterm = GoTermFinding(self._paths.goterm_output_folder, self._args.annotation_path)
        goterm.run_go_term(
                self._args.annotation_path,
                self._paths.goterm_output_folder, self._args.UniProt_id,
                self._args.go_obo, self._args.goslim_obo)

    def srna_target(self):
        """sRNA target prediction"""
        srnatarget = sRNATargetPrediction(self._paths.starget_output_folder, self._args.sRNA_path,
                                            self._args.fasta_path, self._args.annotation_path)
        project_creator.create_subfolders(self._paths.required_folders("srna_target"))
        srnatarget.run_srna_target_prediction(
                self._args.Vienna_folder, self._args.annotation_path,
                self._args.fasta_path, self._args.sRNA_path,
                self._args.query_sRNA, self._args.program, 
                self._args.interaction_length, self._args.window_size_target, 
                self._args.span_target, self._args.window_size_srna, 
                self._args.span_srna, self._args.unstructured_region_RNAplex_target,
                self._args.unstructured_region_RNAplex_srna,
                self._args.unstructured_region_RNAup, self._args.energy_threshold,
                self._args.duplex_distance, self._args.top,
                self._paths.starget_output_folder, self._args.process_rnaplex,
                self._args.process_rnaup, self._args.continue_rnaup)
    def snp(self):
        """SNP transcript detection"""
        project_creator.create_subfolders(self._paths.required_folders("snp"))
        snp = SNPCalling(self._args.bam_type, self._paths.snp_output_folder,
                          self._args.fasta_path)
        snp.run_snp_calling(
                self._args.samtools_path, self._args.bcftools_path,
                self._args.bam_type, 
                self._args.program, self._args.fasta_path,
                self._args.tex_bam_path, self._args.frag_bam_path,
                self._args.quality, self._args.read_depth,
                self._paths.snp_output_folder, self._args.indel_fraction)

    def ppi(self):
        """PPI network retrieve"""
        project_creator.create_subfolders(self._paths.required_folders("ppi_network"))
        ppi = PPINetwork(self._paths.ppi_output_folder)
        ppi.retrieve_ppi_network(
            self._args.ptt_path, self._args.proteinID_strains,
            self._args.without_strain_pubmed, self._args.species_STRING,
            self._args.score, self._paths.ppi_output_folder,
            self._args.node_size, self._args.query)

    def sublocal(self):
        """Subcellular Localization prediction"""
        project_creator.create_subfolders(self._paths.required_folders("subcellular_localization"))
        sublocal = SubLocal(self._args.gff_path, self._args.fasta_path, self._paths.sublocal_output_folder)
        sublocal.run_sub_local(
            self._args.Psortb_path, self._args.gff_path,
            self._args.fasta_path, self._args.bacteria_type,
            self._args.difference_multi, self._args.merge_to_gff,
            self._paths.sublocal_output_folder)

    def ribos(self):
        """riboswitch prediction"""
        project_creator.create_subfolders(self._paths.required_folders("riboswitch"))
        ribos = Ribos(self._args.gff_path, self._args.fasta_path,
                      self._paths.ribos_output_folder, self._paths.database_folder)
        ribos.run_ribos(
            self._args.infernal_path, self._args.riboswitch_ID,
            self._args.gff_path, self._args.fasta_path,
            self._args.Rfam, self._paths.ribos_output_folder,
            self._args.re_scan, self._args.e_value,
            self._args.output_all,
            self._paths.database_folder, self._args.fuzzy)
    def screen(self):
        """generate screenshot"""
        screen = Screen(self._args.output_folder, self._args.fasta)
        screen.screenshot(
            self._args.main_gff,
            self._args.side_gffs, self._args.fasta,
            self._args.frag_wig_folder, self._args.tex_wig_folder,
            self._args.height, self._args.tex_libs,
            self._args.frag_libs, self._args.present,
            self._args.output_folder)
