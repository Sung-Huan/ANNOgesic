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
from annogesiclib.crispr import Crispr
from annogesiclib.merge_feature import run_merge
from annogesiclib.screen import Screen
from annogesiclib.args_container import ArgsContainer
from annogesiclib.helper import Helper
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
        self.args_container = ArgsContainer()
        self.helper = Helper()

    def check_folder(self, folders):
        '''Check the emtpy or wrong assigned folder'''
        for folder in folders:
            if folder is None:
                print("Error: There is wrong path of folder assigned, "
                      "please check it!!")
                sys.exit()
            else:
                if os.path.exists(folder):
                    if len(os.listdir(folder)) == 0:
                        print("Error: There is empty folder, "
                              "please check it!!")
                        sys.exit()
                else:
                    print("Error: There is wrong folder, please check it!!")
                    sys.exit()

    def check_parameter(self, paras, names):
        '''Check the parameter is assigned correct or not'''
        for i in range(len(paras)):
            if paras[i] is None:
                print("Error: {0} is wrong, "
                      "please check it!!".format(names[i]))
                sys.exit()

    def check_no_require_folder(self, folders):
        '''Check the folders which are not necessary.
        It should not be assigned a empty or wrong folder'''
        for folder in folders:
            if folder is not None:
                if os.path.exists(folder):
                    if len(os.listdir(folder)) == 0:
                        print("Error: There is empty folder, "
                              "please check it!!")
                        sys.exit()
                else:
                    print("Error: There is wrong folder, "
                          "please check it!!")
                    sys.exit()

    def check_file(self, files, names, require):
        '''Check the path of file'''
        for i in range(len(files)):
            if require:
                if files[i] is None:
                    print("Error: {0} is wrong, "
                          "please check it!!".format(names[i]))
                    sys.exit()
                else:
                    if not os.path.isfile(files[i]):
                        print("Error: There is wrong path of {0}, "
                              "please check it!!".format(names[i]))
                        sys.exit()
            else:
                if files[i] is not None:
                    if not os.path.isfile(files[i]):
                        print("Error: There is wrong path of {0}, "
                              "please check it!!".format(names[i]))
                        sys.exit()

    def create_project(self, version):
        """Create a new project."""
        project_creator.create_root_folder(self._args.project_path)
        project_creator.create_subfolders(self._paths.required_folders("root"))
        project_creator.create_subfolders(
            self._paths.required_folders("get_target_fasta"))
        project_creator.create_version_file(
            self._paths.version_path, version)
        sys.stdout.write("Created folder \"%s\" and required subfolders.\n" % (
            self._args.project_path))

    def get_input(self):
        """Download required files from website."""
        print("Running get input files...")
        if self._args.FTP_path is None:
            print("Error: Please assign the path for downloading the data!!")
            sys.exit()
        if self._args.for_target:
            annotation_folder = self._paths.tar_annotation_folder
            fasta_folder = self._paths.tar_fasta_folder
        else:
            annotation_folder = self._paths.ref_annotation_folder
            fasta_folder = self._paths.ref_fasta_folder
        self.helper.check_make_folder(annotation_folder)
        self.helper.check_make_folder(fasta_folder)
        if self._args.ref_gff is True:
            get_file(self._args.FTP_path, annotation_folder,
                     "gff", self._args.for_target)
            get_file(self._args.FTP_path, annotation_folder,
                     "_genomic.gff.gz", self._args.for_target)
        if self._args.ref_fasta is True:
            get_file(self._args.FTP_path, fasta_folder,
                     "fna", self._args.for_target)
            get_file(self._args.FTP_path, fasta_folder,
                     "_genomic.fna.gz", self._args.for_target)
        if self._args.ref_gbk is True:
            get_file(self._args.FTP_path, annotation_folder,
                     "gbk", self._args.for_target)
            get_file(self._args.FTP_path, annotation_folder,
                     "gbff", self._args.for_target)
            get_file(self._args.FTP_path, annotation_folder,
                     "_genomic.gbff.gz", self._args.for_target)
        if self._args.ref_ptt is True:
            get_file(self._args.FTP_path, annotation_folder,
                     "ptt", self._args.for_target)
        if self._args.ref_rnt is True:
            get_file(self._args.FTP_path, annotation_folder,
                     "rnt", self._args.for_target)
        if self._args.convert_embl is True:
            annotation_files = os.listdir(annotation_folder)
            if len(annotation_files) == 0:
                sys.stdout.write("No gff files!!\n")
            else:
                Converter().convert_gbk2embl(annotation_folder)

    def get_target_fasta(self):
        """Get target fasta"""
        print("Running get target fasta...")
        self.check_parameter([self._args.output_format], ["--output_format"])
        self.check_folder([self._args.ref_fasta_folder])
        self.check_file([self._args.mutation_table], "--mutation_table", True)
        project_creator.create_subfolders(
            self._paths.required_folders("get_target_fasta"))
        outputs = self._args.output_format.split(",")
        for output in outputs:
            output = output.strip()
        target = TargetFasta(self._paths.tar_fasta_folder,
                             self._args.ref_fasta_folder)
        target.get_target_fasta(
                self._args.mutation_table, self._paths.tar_fasta_folder,
                self._args.ref_fasta_folder, outputs)

    def ratt(self):
        """Run RATT to transfer annotation file from reference to target."""
        print("Running annotation transfer...")
        if (self._args.transfer_type != "Strain") and (
                self._args.transfer_type != "Assembly") and (
                self._args.transfer_type != "Species") and (
                self._args.transfer_type != "Assembly.Repetitive") and (
                self._args.transfer_type != "Strain.Repetitive") and (
                self._args.transfer_type != "Species.Repetitive") and (
                self._args.transfer_type != "Multiple") and (
                self._args.transfer_type != "Free"):
            print("Error: please assign correct --transfer_type!!")
            sys.exit()
        if (self._args.ref_embl is None) and (self._args.ref_gbk is None):
            print("Error: please assign proper embl or genbank folder")
            sys.exit()
        elif (self._args.ref_embl is not None) and (
                self._args.ref_gbk is not None):
            print("Error: please choose embl as input or genbank as input")
            sys.exit()
        self.check_folder([self._args.target_fasta,
                           self._args.ref_fasta])
        self.check_parameter([self._args.element, self._args.compare_pair],
                             ["--element", "--compare_pair"])
        project_creator.create_subfolders(
            self._paths.required_folders("annotation_transfer"))
        args_ratt = self.args_container.container_ratt(
            self._args.RATT_path, self._args.element, self._args.transfer_type,
            self._args.ref_embl, self._args.ref_gbk, self._args.target_fasta,
            self._args.ref_fasta, self._paths.ratt_folder,
            self._args.convert_to_gff_rnt_ptt,
            self._paths.tar_annotation_folder, self._args.compare_pair)
        ratt = RATT(args_ratt)
        ratt.annotation_transfer(args_ratt)

    def tsspredator(self):
        """Run TSSpredator for predicting TSS candidates."""
        self.check_folder([self._args.fasta_folder,
                           self._args.annotation_folder,
                           self._args.wig_folder])
        self.check_parameter([self._args.lib, self._args.output_prefix],
                             ["--lib", "--output_prefix"])
        self.check_no_require_folder([self._args.compare_transcript_assembly,
                                      self._args.reference_gff_folder])
        self.check_file([self._args.merge_manual], ["--merge_manual"], False)
        if self._args.compute_program.lower() == "tss":
            print("Running TSS prediction...")
            project_creator.create_subfolders(
                self._paths.required_folders("TSS"))
            out_folder = self._paths.tsspredator_folder
        elif self._args.compute_program.lower() == "processing_site":
            print("Running processing site prediction...")
            out_folder = self._paths.processing_site_folder
            project_creator.create_subfolders(
                self._paths.required_folders("processing"))
        else:
            print("Error:No such program!!!!")
            sys.exit()
        args_tss = self.args_container.container_tsspredator(
            self._args.TSSpredator_path, self._args.compute_program,
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
            self._args.utr_length, self._args.cluster,
            self._args.partial_length, self._args.re_check_orphan,
            self._args.overlap_feature, self._args.reference_gff_folder,
            self._args.remove_low_expression)
        tsspredator = TSSpredator(args_tss)
        tsspredator.run_tsspredator(args_tss)

    def optimize(self):
        """opimize TSSpredator"""
        self.check_folder([self._args.wig_folder, self._args.fasta_file,
                           self._args.annotation_file])
        self.check_file([self._args.manual],
                        ["--manual"], True)
        self.check_parameter([self._args.strain_name, self._args.lib,
                              self._args.output_prefix],
                             ["--strain_name", "--lib", "--output_prefix"])
        if self._args.program.lower() == "tss":
            print("Running optimization of TSS prediction...")
            project_creator.create_subfolders(
                self._paths.required_folders("TSS"))
            out_folder = self._paths.tsspredator_folder
        elif self._args.program.lower() == "processing_site":
            print("Running optimization of processing site prediction...")
            out_folder = self._paths.processing_site_folder
            project_creator.create_subfolders(
                self._paths.required_folders("processing"))
        else:
            print("Error:No such program!!!!")
            sys.exit()
        args_ops = self.args_container.container_optimize(
            self._args.TSSpredator_path, self._args.fasta_file,
            self._args.annotation_file, self._args.wig_folder,
            self._args.manual, out_folder,
            self._args.strain_name, self._args.max_height,
            self._args.max_height_reduction, self._args.max_factor,
            self._args.max_factor_reduction, self._args.max_base_height,
            self._args.max_enrichment_factor, self._args.max_processing_factor,
            self._args.utr_length, self._args.lib,
            self._args.output_prefix, self._args.cluster,
            self._args.partial_length, self._args.core,
            self._args.program, self._args.replicate_match,
            self._args.steps)
        optimize_tss(args_ops)

    def color(self):
        """color the screenshots"""
        print("Running png files coloring...")
        self.check_parameter([self._args.track_number], ["--track_numer"])
        self.check_folder([self._args.screenshot_folder])
        color = ColorPNG()
        color.generate_color_png(
                self._args.track_number, self._args.screenshot_folder,
                self._args.ImageMagick_covert_path)

    def terminator(self):
        """Run TransTermHP for detecting terminators."""
        print("Running terminator prediction...")
        if self._args.TransTermHP_path is None:
            print("Please assign the folder where you install TransTermHP.")
        self.check_folder([self._args.fasta_folder,
                           self._args.annotation_folder,
                           self._args.transcript_folder])
        self.check_no_require_folder([
            self._args.sRNA, self._args.tex_wig_folder,
            self._args.frag_wig_folder])
        project_creator.create_subfolders(
            self._paths.required_folders("terminator"))
        args_term = self.args_container.container_terminator(
            self._args.TransTermHP_path, self._args.expterm_path,
            self._args.RNAfold_path,
            self._paths.transterm_folder, self._args.fasta_folder,
            self._args.annotation_folder, self._args.transcript_folder,
            self._args.sRNA, self._args.statistics,
            self._args.tex_wig_folder, self._args.frag_wig_folder,
            self._args.decrease, self._args.highest_coverage,
            self._args.fuzzy_detect_coverage,
            self._args.fuzzy_within_transcript,
            self._args.fuzzy_downstream_transcript,
            self._args.fuzzy_within_gene,
            self._args.fuzzy_downstream_gene, self._paths.transtermhp_folder,
            self._args.tex_notex_libs, self._args.frag_libs,
            self._args.tex_notex, self._args.replicates_tex,
            self._args.replicates_frag, self._args.table_best,
            self._args.min_loop_length, self._args.max_loop_length,
            self._args.min_stem_length, self._args.max_stem_length,
            self._args.min_U_tail_length, self._args.miss_rate,
            self._args.range_U_tail, self._args.keep_multi_term)
        terminator = Terminator(args_term)
        terminator.run_terminator(args_term)

    def transcript(self):
        """Run Transcriptome assembly."""
        print("Running transcript assembly...")
        self.check_folder([self._args.annotation_folder])
        self.check_no_require_folder([
            self._args.compare_TSS, self._args.compare_genome_annotation,
            self._args.terminator_folder, self._args.frag_wig_path,
            self._args.tex_wig_path])
        project_creator.create_subfolders(
            self._paths.required_folders("transcript_assembly"))
        print(self._args.tex_wig_path)
        args_tran = self.args_container.container_transcript(
            self._args.frag_wig_path, self._args.tex_wig_path,
            self._args.tex_notex,
            self._args.length, self._args.annotation_folder,
            self._args.height, self._args.width,
            self._args.tolerance, self._args.tolerance_coverage,
            self._args.replicates_tex, self._args.replicates_frag,
            self._paths.transcript_assembly_output_folder,
            self._args.compare_TSS, self._args.compare_genome_annotation,
            self._args.TSS_fuzzy, self._args.Tex_treated_libs,
            self._args.fragmented_libs, self._args.compare_feature_genome,
            self._args.table_best, self._args.terminator_folder,
            self._args.fuzzy_term, self._args.max_length_distribution)
        transcript = TranscriptAssembly(args_tran)
        transcript.run_transcript_assembly(args_tran)

    def utr_detection(self):
        """Run UTR detection."""
        print("Running UTR detection...")
        self.check_folder([self._args.annotation_folder,
                           self._args.transcript_assembly_folder,
                           self._args.TSS_folder])
        self.check_no_require_folder([self._args.terminator_folder])
        project_creator.create_subfolders(self._paths.required_folders("utr"))
        args_utr = self.args_container.container_utr(
                self._args.TSS_folder, self._args.annotation_folder,
                self._args.transcript_assembly_folder,
                self._args.terminator_folder,
                self._args.terminator_fuzzy, self._paths.utr_folder,
                self._args.TSS_source, self._args.base_5UTR,
                self._args.UTR_length, self._args.base_3UTR,
                self._args.fuzzy_3utr, self._args.fuzzy_5utr)
        utr = UTRDetection(args_utr)
        utr.run_utr_detection(args_utr)

    def srna_detection(self):
        """sRNA_detection."""
        print("Running sRNA prediction...")
        self.check_folder([self._args.annotation_folder,
                           self._args.transcript_assembly_folder])
        self.check_no_require_folder([
            self._args.fasta_folder, self._args.sORF,
            self._args.terminator_folder, self._args.tex_wig_folder,
            self._args.frag_wig_folder, self._args.processing_site_folder])
        self.check_file([self._args.promoter_table],
                        ["--promoter_table"], False)
        if self._args.UTR_derived_sRNA:
            self.check_folder([self._args.TSS_folder])
        else:
            self.check_no_require_folder([self._args.TSS_folder])
        project_creator.create_subfolders(self._paths.required_folders("srna"))
        args_srna = self.args_container.container_srna(
                self._args.Vienna_folder, self._args.Vienna_utils,
                self._args.blast_plus_folder,
                self._args.ps2pdf14_path, self._paths.srna_folder,
                self._args.UTR_derived_sRNA, self._args.annotation_folder,
                self._args.TSS_folder, self._args.transcript_assembly_folder,
                self._args.TSS_intergenic_fuzzy, self._args.TSS_5UTR_fuzzy,
                self._args.TSS_3UTR_fuzzy, self._args.TSS_interCDS_fuzzy,
                self._args.filter_info, self._args.tex_wig_folder,
                self._args.frag_wig_folder, self._args.processing_site_folder,
                self._args.fasta_folder, self._args.mountain_plot,
                self._args.nr_format, self._args.srna_format,
                self._args.sRNA_database_path, self._args.nr_database_path,
                self._args.cutoff_energy,
                self._args.run_intergenic_TEX_coverage,
                self._args.run_intergenic_noTEX_coverage,
                self._args.run_intergenic_fragmented_coverage,
                self._args.run_break_transcript,
                self._args.run_antisense_TEX_coverage,
                self._args.run_antisense_noTEX_coverage,
                self._args.run_antisense_fragmented_coverage,
                self._args.run_utr_TEX_coverage,
                self._args.run_utr_noTEX_coverage,
                self._args.run_utr_fragmented_coverage,
                self._args.max_length, self._args.min_length,
                self._args.tex_notex_libs, self._args.frag_libs,
                self._args.replicates_tex, self._args.replicates_frag,
                self._args.tex_notex, self._args.blast_e_nr,
                self._args.blast_e_srna, self._args.detect_sRNA_in_CDS,
                self._args.table_best, self._args.decrease_intergenic_antisense,
                self._args.decrease_utr, self._args.fuzzy_intergenic_antisense,
                self._args.fuzzy_utr, self._args.cutoff_nr_hit,
                self._args.sORF,
                self._args.overlap_percent_CDS,
                self._args.terminator_folder,
                self._args.terminator_fuzzy_in_sRNA,
                self._args.terminator_fuzzy_out_sRNA,
                self._args.ignore_hypothetical_protein, self._args.TSS_source,
                self._args.min_utr_coverage, self._args.promoter_table,
                self._args.ranking_time_promoter, self._args.promoter_name)
        srna = sRNADetection(args_srna)
        srna.run_srna_detection(args_srna)

    def sorf_detection(self):
        """sORF_detection."""
        print("Running sORF prediction...")
        self.check_folder([self._args.transcript_assembly_folder,
                           self._args.annotation_folder,
                           self._args.fasta_folder])
        self.check_no_require_folder([
            self._args.sRNA_folder, self._args.TSS_folder,
            self._args.tex_wig_folder, self._args.frag_wig_folder])
        project_creator.create_subfolders(
            self._paths.required_folders("sorf"))
        args_sorf = self.args_container.container_sorf(
            self._paths.sorf_folder, self._args.UTR_derived_sORF,
            self._args.transcript_assembly_folder,
            self._args.annotation_folder,
            self._args.TSS_folder, self._args.utr_length,
            self._args.min_length, self._args.max_length,
            self._args.tex_wig_folder, self._args.frag_wig_folder,
            self._args.cutoff_intergenic_coverage,
            self._args.cutoff_antisense_coverage,
            self._args.cutoff_5utr_coverage,
            self._args.cutoff_3utr_coverage,
            self._args.cutoff_interCDS_coverage,
            self._args.fasta_folder, self._args.tex_notex_libs,
            self._args.frag_libs, self._args.tex_notex,
            self._args.replicates_tex, self._args.replicates_frag,
            self._args.table_best, self._args.sRNA_folder,
            self._args.start_codon, self._args.stop_codon,
            self._args.cutoff_background, self._args.fuzzy_rbs,
            self._args.rbs_not_after_TSS, self._args.print_all_combination,
            self._args.best_no_sRNA, self._args.best_no_TSS,
            self._args.ignore_hypothetical_protein,
            self._args.min_rbs_distance, self._args.max_rbs_distance)
        sorf = sORFDetection(args_sorf)
        sorf.run_sorf_detection(args_sorf)

    def meme(self):
        """promoter detectopn"""
        print("Running promoter detection...")
        self.check_folder([self._args.TSS_folder, self._args.fasta_folder])
        self.check_no_require_folder([self._args.tex_wig_path])
        if not self._args.TSS_source:
            self.check_folder([self._args.annotation_folder])
        project_creator.create_subfolders(
            self._paths.required_folders("promoter"))
        args_pro = self.args_container.container_promoter(
            self._args.MEME_path,
            self._paths.promoter_output_folder, self._args.tex_libs,
            self._args.TSS_folder, self._args.fasta_folder,
            self._args.num_motif, self._args.nt_before_TSS,
            self._args.motif_width, self._args.TSS_source,
            self._args.tex_wig_path, self._args.annotation_folder,
            self._args.combine_all, self._args.e_value,
            self._args.parallel)
        meme = MEME(args_pro)
        meme.run_meme(args_pro)

    def operon(self):
        """operon detection"""
        print("Running operon detection...")
        self.check_folder([self._args.TSS_folder, self._args.annotation_folder,
                           self._args.transcript_folder,
                           self._args.UTR5_folder, self._args.UTR3_folder])
        self.check_no_require_folder([self._args.term_folder])
        project_creator.create_subfolders(
            self._paths.required_folders("operon"))
        args_op = self.args_container.container_operon(
            self._args.TSS_folder, self._args.annotation_folder,
            self._args.transcript_folder, self._args.UTR5_folder,
            self._args.UTR3_folder, self._args.term_folder,
            self._args.TSS_fuzzy, self._args.term_fuzzy,
            self._args.min_length, self._args.statistics,
            self._paths.operon_output_folder, self._args.combine_gff,
            self._paths.operon_statistics_folder)
        operon = OperonDetection(args_op)
        operon.run_operon(args_op)

    def circrna(self):
        """circRNA detection"""
        print("Running circular RNA prediction...")
        self.check_folder([self._args.fasta_path, self._args.annotation_path])
        self.check_no_require_folder([self._args.tex_bam_path,
                                      self._args.fragmented_bam_path])
        project_creator.create_subfolders(
            self._paths.required_folders("circrna"))
        args_circ = self.args_container.container_circrna(
            self._args.align, self._args.process, self._args.fasta_path,
            self._args.annotation_path, self._args.tex_bam_path,
            self._args.fragmented_bam_path,
            self._args.read_path, self._paths.circrna_stat_folder,
            self._args.support_reads,
            self._args.segemehl_folder, self._args.samtools_path,
            self._args.start_ratio, self._args.end_ratio,
            self._args.ignore_hypothetical_protein,
            self._paths.circrna_output_folder)
        circ = CircRNADetection(args_circ)
        circ.run_circrna(args_circ)

    def goterm(self):
        """Go term discovery"""
        print("Running GO term mapping...")
        self.check_folder([self._args.annotation_path])
        self.check_no_require_folder([self._args.transcript_path])
        self.check_file([self._args.UniProt_id, self._args.go_obo,
                         self._args.goslim_obo],
                        ["--UniProt_id", "--go.obo", "--goslim_obo"], True)
        project_creator.create_subfolders(
            self._paths.required_folders("go_term"))
        args_go = self.args_container.container_goterm(
            self._args.annotation_path,
            self._paths.goterm_output_folder, self._args.UniProt_id,
            self._args.go_obo, self._args.goslim_obo,
            self._args.transcript_path)
        goterm = GoTermFinding(args_go)
        goterm.run_go_term(args_go)

    def srna_target(self):
        """sRNA target prediction"""
        print("Running sRNA target prediction...")
        self.check_folder([self._args.fasta_path, self._args.sRNA_path,
                           self._args.annotation_path])
        project_creator.create_subfolders(
            self._paths.required_folders("srna_target"))
        args_tar = self.args_container.container_srna_target(
            self._args.Vienna_folder, self._args.annotation_path,
            self._args.fasta_path, self._args.sRNA_path,
            self._args.query_sRNA, self._args.program,
            self._args.interaction_length, self._args.window_size_target,
            self._args.span_target, self._args.window_size_srna,
            self._args.span_srna,
            self._args.unstructured_region_RNAplex_target,
            self._args.unstructured_region_RNAplex_srna,
            self._args.unstructured_region_RNAup, self._args.energy_threshold,
            self._args.duplex_distance, self._args.top,
            self._paths.starget_output_folder, self._args.process_rnaplex,
            self._args.process_rnaup, self._args.continue_rnaup,
            self._args.potential_target_start, self._args.potential_target_end,
            self._args.target_feature)
        srnatarget = sRNATargetPrediction(args_tar)
        srnatarget.run_srna_target_prediction(args_tar)

    def snp(self):
        """SNP transcript detection"""
        print("Running SNP/mutations calling...")
        self.check_folder([self._args.fasta_path])
        if (self._args.bam_type != "target") and (
                self._args.bam_type != "reference"):
            print("Error: please assign \"target\" or"
                  " \"reference\" to --bam_type!!")
            sys.exit()
        if (self._args.ploidy != "haploid") and (
                self._args.ploidy != "diploid"):
            print("Error: please assign \"haploid\" or"
                  " \"diploid\" to --chromosome_type!!")
        if (self._args.caller != "c") and (
                self._args.caller != "m"):
            print("Error: please assign \"c\" or"
                  " \"m\" to --caller!!")
        project_creator.create_subfolders(self._paths.required_folders("snp"))
        args_snp = self.args_container.container_snp(
            self._args.samtools_path, self._args.bcftools_path,
            self._args.bam_type, self._args.sample_number,
            self._args.program, self._args.fasta_path,
            self._args.tex_bam_path, self._args.frag_bam_path,
            self._args.quality, self._args.read_depth_range,
            self._paths.snp_output_folder, self._args.indel_fraction,
            self._args.ploidy, self._args.RG_tag, self._args.caller,
            self._args.filter_tag_info, self._args.DP4_cutoff)
        snp = SNPCalling(args_snp)
        snp.run_snp_calling(args_snp)

    def ppi(self):
        """PPI network retrieve"""
        print("Running protein-protein interaction networks prediction...")
        self.check_folder([self._args.gff_path])
        self.check_parameter([self._args.proteinID_strains,
                              self._args.species_STRING],
                             ["--proteinID_strains", "--species_STRING"])
        project_creator.create_subfolders(
            self._paths.required_folders("ppi_network"))
        args_ppi = self.args_container.container_ppi(
            self._args.gff_path, self._args.proteinID_strains,
            self._args.without_strain_pubmed, self._args.species_STRING,
            self._args.score, self._paths.ppi_output_folder,
            self._args.node_size, self._args.query)
        ppi = PPINetwork(self._paths.ppi_output_folder)
        ppi.retrieve_ppi_network(args_ppi)

    def sublocal(self):
        """Subcellular Localization prediction"""
        print("Running subcellular localization prediction...")
        self.check_folder([self._args.gff_path, self._args.fasta_path])
        self.check_no_require_folder([self._args.transcript_path])
        if (self._args.bacteria_type != "positive") and (
                self._args.bacteria_type != "negative"):
            print("Error: please assign \"positive\" or"
                  " \"negative\" to --bacteria_type!!")
            sys.exit()
        project_creator.create_subfolders(
            self._paths.required_folders("subcellular_localization"))
        args_sub = self.args_container.container_sublocal(
            self._args.Psortb_path, self._args.gff_path,
            self._args.fasta_path, self._args.bacteria_type,
            self._args.difference_multi, self._args.merge_to_gff,
            self._paths.sublocal_output_folder, self._args.transcript_path)
        sublocal = SubLocal(args_sub)
        sublocal.run_sub_local(args_sub)

    def ribos(self):
        """riboswitch and RNA thermometer prediction"""
        print("Running riboswitch and RNA thermometer prediction...")
        self.check_folder([self._args.gff_path, self._args.fasta_path,
                           self._args.tss_path, self._args.transcript_path])
        if (self._args.program == "both"):
            self.check_file([self._args.riboswitch_ID, self._args.Rfam],
                            ["--riboswitch_ID", "--Rfam"], True)
            self.check_file([self._args.RNA_thermometer_ID, self._args.Rfam],
                            ["--RNA_thermometer_ID", "--Rfam"], True)
            project_creator.create_subfolders(
                    self._paths.required_folders("riboswitch"))
            project_creator.create_subfolders(
                    self._paths.required_folders("thermometer"))
            ribos_path = self._paths.ribos_output_folder
            thermo_path = self._paths.thermo_output_folder
        elif (self._args.program == "thermometer"):
            self.check_file([self._args.RNA_thermometer_ID, self._args.Rfam],
                            ["--thermometer_ID", "--Rfam"], True)
            project_creator.create_subfolders(
                    self._paths.required_folders("thermometer"))
            ribos_path = None
            thermo_path = self._paths.thermo_output_folder
        elif (self._args.program == "riboswitch"):
            self.check_file([self._args.riboswitch_ID, self._args.Rfam],
                            ["--riboswitch_ID", "--Rfam"], True)
            project_creator.create_subfolders(
                    self._paths.required_folders("riboswitch"))
            ribos_path = self._paths.ribos_output_folder
            thermo_path = None
        else:
            print("Error: Please assign \"thermometer\", \"riboswitch\" "
                  "or \"both\" in --program!!")
            sys.exit()
        args_ribo = self.args_container.container_ribos(
            self._args.program, self._args.RNA_thermometer_ID,
            self._args.infernal_path, self._args.riboswitch_ID,
            self._args.gff_path, self._args.fasta_path,
            self._args.tss_path, self._args.transcript_path,
            self._args.Rfam, ribos_path,
            thermo_path, self._args.e_value,
            self._args.output_all, self._paths.database_folder,
            self._args.fuzzy, self._args.start_codon,
            self._args.min_dist_rbs, self._args.max_dist_rbs,
            self._args.fuzzy_rbs, self._args.UTR_length)
        ribos = Ribos(args_ribo)
        ribos.run_ribos(args_ribo)

    def crispr(self):
        """CRISPR prediction"""
        print("Running CRISPR prediction...")
        self.check_folder([self._args.fasta_path])
        self.check_no_require_folder([self._args.gff_path])
        project_creator.create_subfolders(
            self._paths.required_folders("crispr"))
        args_cris = self.args_container.container_cris(
            self._args.fasta_path, self._args.gff_path,
            self._args.CRT_path, self._args.window_size,
            self._args.min_number_repeat, self._args.min_length_repeat,
            self._args.Max_length_repeat, self._args.min_length_spacer,
            self._args.Max_length_spacer, self._paths.crispr_output_folder,
            self._args.ignore_hypothetical_protein)
        cris = Crispr(args_cris)
        cris.run_crispr(args_cris)

    def merge(self):
        """Merge all features"""
        print("Merging all features to one gff file...")
        merge_folder = os.path.join(self._paths.output_folder,
                                    "merge_all_features")
        self.helper.check_make_folder(merge_folder)
        other_features = self._args.other_features_path.split(",")
        self.check_file([self._args.transcript_path] + other_features,
                        ["--transcript_path", "--other_features_path"],
                        False)
        self.check_parameter([self._args.strain_name], ["--strain_name"])
        run_merge(merge_folder, self._args.transcript_path,
                  self._args.other_features_path, self._args.fuzzy_term,
                  self._args.fuzzy_TSS,
                  os.path.join(merge_folder, self._args.strain_name))

    def screen(self):
        """generate screenshot"""
        print("Running screenshot generating...")
        self.check_file([self._args.main_gff, self._args.fasta],
                        ["--main_gff", "--fasta"], True)
        if self._args.side_gffs is not None:
            for gff in (self._args.side_gffs.split(",")):
                gff = gff.strip()
                if not os.path.isfile(gff):
                    print("Error: The --side_gffs no exist!!")
                    sys.exit()
        if self._args.output_folder is None:
            print("Error: please assign --output_folder!!")
            sys.exit()
        if (self._args.present != "expand") and (
                self._args.present != "collapse") and (
                self._args.present != "squish"):
            print("Error: please assign \"expand\" or "
                  "\"collapse\" or \"squish\" to --present!!")
            sys.exit()
        args_sc = self.args_container.container_screen(
            self._args.main_gff, self._args.side_gffs,
            self._args.fasta, self._args.frag_wig_folder,
            self._args.tex_wig_folder, self._args.height,
            self._args.tex_libs, self._args.frag_libs,
            self._args.present, self._args.output_folder)
        screen = Screen(args_sc)
        screen.screenshot(args_sc)
