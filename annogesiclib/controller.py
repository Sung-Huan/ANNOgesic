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
from annogesiclib.transcript import TranscriptDetection
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

    def check_folder(self, folders, flags):
        '''Check the emtpy or wrong assigned folder'''
        for folder, flag in zip(folders, flags):
            if folder is None:
                print("Error: {0} is wrong. Please check it!".format(flag))
                sys.exit()
            else:
                if os.path.exists(folder):
                    if len(os.listdir(folder)) == 0:
                        print("Error: {0} is empty folder!".format(flag))
                        sys.exit()
                else:
                    print("Error: {0} is wrong. Please check it!".format(
                          flag))
                    sys.exit()

    def check_multi_files(self, input_files, flags):
        if input_files is not None:
            for files, flag in zip(input_files, flags):
                if files is not None:
                    for file_ in files:
                        if not os.path.exists(file_):
                            print("Error: Some files in {0} are "
                                  "no existed!".format(flag))
                            sys.exit()

    def check_parameter(self, paras, names):
        '''Check the parameter is assigned correct or not'''
        for i in range(len(paras)):
            if paras[i] is None:
                print("Error: {0} is wrong. "
                      "Please check it!".format(names[i]))
                sys.exit()

    def check_no_require_folder(self, folders):
        '''Check the folders which are not necessary.
        It should not be assigned a empty or wrong folder'''
        for folder in folders:
            if folder is not None:
                if os.path.exists(folder):
                    if len(os.listdir(folder)) == 0:
                        print("Error: There is empty folder. "
                              "Please check it!")
                        sys.exit()
                else:
                    print("Error: There is wrong folder. "
                          "Please check it!")
                    sys.exit()

    def check_execute_file(self, exe):
        detect = False
        if os.path.exists(exe):
            detect = True
        for folder in os.environ["PATH"].split(":"):
            if exe in os.listdir(folder):
                detect = True
        if not detect:
            print("Error: {0} can't be found!".format(exe))
            sys.exit()

    def check_file(self, files, names, require):
        '''Check the path of file'''
        for i in range(len(files)):
            if require:
                if files[i] is None:
                    print("Error: {0} is wrong. "
                          "Please check it!".format(names[i]))
                    sys.exit()
                else:
                    if not os.path.isfile(files[i]):
                        print("Error: There is wrong path of {0}. "
                              "Please check it!".format(names[i]))
                        sys.exit()
            else:
                if files[i] is not None:
                    if not os.path.isfile(files[i]):
                        print("Error: There is wrong path of {0}. "
                              "Please check it!".format(names[i]))
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
        print("Running get input files")
        if self._args.ftp_path is None:
            print("Error: Please assign the path for downloading the data!")
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
            get_file(self._args.ftp_path, annotation_folder,
                     "gff", self._args.for_target)
            get_file(self._args.ftp_path, annotation_folder,
                     "_genomic.gff.gz", self._args.for_target)
        if self._args.ref_fasta is True:
            get_file(self._args.ftp_path, fasta_folder,
                     "fna", self._args.for_target)
            get_file(self._args.ftp_path, fasta_folder,
                     "_genomic.fna.gz", self._args.for_target)
        if self._args.ref_gbk is True:
            get_file(self._args.ftp_path, annotation_folder,
                     "gbk", self._args.for_target)
            get_file(self._args.ftp_path, annotation_folder,
                     "gbff", self._args.for_target)
            get_file(self._args.ftp_path, annotation_folder,
                     "_genomic.gbff.gz", self._args.for_target)
        if self._args.ref_ptt is True:
            get_file(self._args.ftp_path, annotation_folder,
                     "ptt", self._args.for_target)
        if self._args.ref_rnt is True:
            get_file(self._args.ftp_path, annotation_folder,
                     "rnt", self._args.for_target)
        if self._args.convert_embl is True:
            annotation_files = os.listdir(annotation_folder)
            if len(annotation_files) == 0:
                sys.stdout.write("No gff files!!\n")
            else:
                Converter().convert_gbk2embl(annotation_folder)

    def get_target_fasta(self):
        """Get target fasta"""
        print("Running get target fasta")
        self.check_multi_files([self._args.ref_fasta_files], ["--ref_fasta_files"])
        self.check_parameter([self._args.output_format], ["--output_format"])
        self.check_file([self._args.mutation_table], "--mutation_table", True)
        project_creator.create_subfolders(
            self._paths.required_folders("get_target_fasta"))
        for output in self._args.output_format:
            output = output.strip()
        target = TargetFasta(self._paths.tar_fasta_folder,
                             self._args.ref_fasta_files)
        target.get_target_fasta(
                self._args.mutation_table, self._paths.tar_fasta_folder,
                self._args.ref_fasta_files, self._args.output_format,
                self._paths.target_base_folder)

    def ratt(self):
        """Run RATT to transfer annotation file from reference to target."""
        print("Running annotation transfer")
        if (self._args.transfer_type != "Strain") and (
                self._args.transfer_type != "Assembly") and (
                self._args.transfer_type != "Species") and (
                self._args.transfer_type != "Assembly.Repetitive") and (
                self._args.transfer_type != "Strain.Repetitive") and (
                self._args.transfer_type != "Species.Repetitive") and (
                self._args.transfer_type != "Multiple") and (
                self._args.transfer_type != "Free"):
            print("Error: please assign correct --transfer_type!")
            sys.exit()
        if (self._args.ref_embl_files is None) and (
                self._args.ref_gbk_files is None):
            print("Error: please assign proper embl or genbank folder")
            sys.exit()
        elif (self._args.ref_embl_files is not None) and (
                self._args.ref_gbk_files is not None):
            print("Error: please choose embl as input or genbank as input")
            sys.exit()
        self.check_execute_file(self._args.ratt_path)
        self.check_multi_files(
                [self._args.target_fasta_files, self._args.ref_fasta_files],
                ["--target_fasta_files", "--ref_fasta_files"])
        self.check_parameter([self._args.element, self._args.compare_pair],
                             ["--element", "--compare_pair"])
        project_creator.create_subfolders(
            self._paths.required_folders("annotation_transfer"))
        args_ratt = self.args_container.container_ratt(
            self._args.ratt_path, self._args.element, self._args.transfer_type,
            self._args.ref_embl_files, self._args.ref_gbk_files,
            self._args.target_fasta_files, self._args.ref_fasta_files,
            self._paths.ratt_folder, self._args.convert_to_gff_rnt_ptt,
            self._paths.tar_annotation_folder, self._args.compare_pair)
        ratt = RATT(args_ratt)
        ratt.annotation_transfer(args_ratt)

    def tsspredator(self):
        """Run TSSpredator for predicting TSS candidates."""
        self.check_multi_files(
                [self._args.fasta_files, self._args.annotation_files,
                 self._args.reference_gff_files, self._args.manual_files,
                 self._args.compare_transcript_files],
                ["--fasta_files", "--annotation_files", "--reference_gff_files",
                 "--manual_files", "--compare_transcript_files"])
        self.check_parameter([self._args.tex_notex_libs, self._args.condition_names],
                             ["--tex_notex_libs", "--condition_names"])
        self.check_execute_file(self._args.tsspredator_path)
        if self._args.compute_program.lower() == "tss":
            print("Running TSS prediction")
            project_creator.create_subfolders(
                self._paths.required_folders("TSS"))
            out_folder = self._paths.tsspredator_folder
        elif self._args.compute_program.lower() == "processing_site":
            print("Running processing site prediction")
            out_folder = self._paths.processing_site_folder
            project_creator.create_subfolders(
                self._paths.required_folders("processing"))
        else:
            print("Error: No such program!")
            sys.exit()
        args_tss = self.args_container.container_tsspredator(
            self._args.tsspredator_path, self._args.compute_program,
            self._args.fasta_files, self._args.annotation_files,
            self._args.tex_notex_libs, self._args.condition_names,
            self._args.height, self._args.height_reduction,
            self._args.factor, self._args.factor_reduction,
            self._args.base_height, self._args.enrichment_factor,
            self._args.processing_factor, self._args.replicate_tex,
            out_folder, self._args.statistics,
            self._args.validate_gene, self._args.manual_files,
            self._args.compare_transcript_files, self._args.fuzzy,
            self._args.utr_length, self._args.cluster,
            self._args.partial_length, self._args.re_check_orphan,
            self._args.overlap_feature, self._args.reference_gff_files,
            self._args.remove_low_expression)
        tsspredator = TSSpredator(args_tss)
        tsspredator.run_tsspredator(args_tss)

    def optimize(self):
        """opimize TSSpredator"""
        self.check_file([self._args.manual],
                        ["--manual"], True)
        self.check_execute_file(self._args.tsspredator_path)
        self.check_parameter([self._args.strain_name, self._args.tex_notex_libs,
                              self._args.condition_names],
                             ["--strain_name", "--tex_notex_lib",
                              "--condition_names"])
        if self._args.program.lower() == "tss":
            print("Running optimization of TSS prediction")
            project_creator.create_subfolders(
                self._paths.required_folders("TSS"))
            out_folder = self._paths.tsspredator_folder
        elif self._args.program.lower() == "processing_site":
            print("Running optimization of processing site prediction")
            out_folder = self._paths.processing_site_folder
            project_creator.create_subfolders(
                self._paths.required_folders("processing"))
        else:
            print("Error: No such program!")
            sys.exit()
        args_ops = self.args_container.container_optimize(
            self._args.tsspredator_path, self._args.fasta_file,
            self._args.annotation_file,
            self._args.manual, out_folder,
            self._args.strain_name, self._args.max_height,
            self._args.max_height_reduction, self._args.max_factor,
            self._args.max_factor_reduction, self._args.max_base_height,
            self._args.max_enrichment_factor, self._args.max_processing_factor,
            self._args.utr_length, self._args.tex_notex_libs,
            self._args.condition_names, self._args.cluster,
            self._args.partial_length, self._args.parallels,
            self._args.program, self._args.replicate_tex,
            self._args.steps)
        optimize_tss(args_ops)

    def color(self):
        """color the screenshots"""
        print("Running png files coloring")
        self.check_parameter([self._args.track_number], ["--track_numer"])
        self.check_folder([self._args.screenshot_folder], ["--screenshot_folder"])
        self.check_execute_file(self._args.imagemagick_covert_path)
        color = ColorPNG()
        color.generate_color_png(
                self._args.track_number, self._args.screenshot_folder,
                self._args.imagemagick_covert_path)

    def terminator(self):
        """Run TransTermHP and Gene converaged for detecting terminators"""
        print("Running terminator prediction")
        if self._args.transtermhp_path is None:
            print("Please assign the folder where you install TransTermHP.")
        self.check_multi_files(
                [self._args.fasta_files, self._args.annotation_files,
                 self._args.transcript_files, self._args.srna_files],
                ["--fasta_files", "--annotation_files",
                 "--transcript_files", "--srna_files"])
        for exe in (self._args.transtermhp_path, self._args.expterm_path,
                    self._args.rnafold_path):
            self.check_execute_file(exe)
        project_creator.create_subfolders(
            self._paths.required_folders("terminator"))
        args_term = self.args_container.container_terminator(
            self._args.transtermhp_path, self._args.expterm_path,
            self._args.rnafold_path,
            self._paths.transterm_folder, self._args.fasta_files,
            self._args.annotation_files, self._args.transcript_files,
            self._args.srna_files, self._args.statistics,
            self._args.decrease, self._args.highest_coverage,
            self._args.fuzzy_detect_coverage,
            self._args.fuzzy_within_transcript,
            self._args.fuzzy_downstream_transcript,
            self._args.fuzzy_within_gene,
            self._args.fuzzy_downstream_gene, self._paths.transtermhp_folder,
            self._args.tex_notex_libs, self._args.frag_libs,
            self._args.tex_notex, self._args.replicate_tex,
            self._args.replicate_frag, self._args.table_best,
            self._args.min_loop_length, self._args.max_loop_length,
            self._args.min_stem_length, self._args.max_stem_length,
            self._args.min_u_tail_length, self._args.miss_rate,
            self._args.range_u_tail, self._args.keep_multi_term,
            self._args.window_size, self._args.window_shift)
        terminator = Terminator(args_term)
        terminator.run_terminator(args_term)

    def transcript(self):
        """Run Transcript detection"""
        print("Running transcript detection")
        self.check_multi_files(
                [self._args.annotation_files, self._args.tss_files,
                 self._args.terminator_files],
                ["--annotation_files", "--tss_files", "--terminator_files"])
        project_creator.create_subfolders(
            self._paths.required_folders("transcript"))
        args_tran = self.args_container.container_transcript(
            self._args.tex_notex,
            self._args.length, self._args.annotation_files,
            self._args.height, self._args.width,
            self._args.tolerance, self._args.tolerance_coverage,
            self._args.replicate_tex, self._args.replicate_frag,
            self._paths.transcript_output_folder,
            self._args.tss_files,
            self._args.tss_fuzzy, self._args.tex_notex_libs,
            self._args.frag_libs, self._args.compare_feature_genome,
            self._args.table_best, self._args.terminator_files,
            self._args.fuzzy_term, self._args.max_length_distribution)
        transcript = TranscriptDetection(args_tran)
        transcript.run_transcript(args_tran)

    def utr_detection(self):
        """Run UTR detection."""
        print("Running UTR detection")
        self.check_multi_files(
            [self._args.annotation_files, self._args.terminator_files,
             self._args.transcript_files, self._args.tss_files],
            ["--annotation_folder", "--terminator_files",
             "--transcript_folder", "--tss_folder"])
        project_creator.create_subfolders(self._paths.required_folders("utr"))
        args_utr = self.args_container.container_utr(
                self._args.tss_files, self._args.annotation_files,
                self._args.transcript_files,
                self._args.terminator_files,
                self._args.terminator_fuzzy, self._paths.utr_folder,
                self._args.tss_source, self._args.base_5utr,
                self._args.utr_length, self._args.base_3utr,
                self._args.fuzzy_3utr, self._args.fuzzy_5utr)
        utr = UTRDetection(args_utr)
        utr.run_utr_detection(args_utr)

    def _check_filter_input(self, files, info, filters):
        if files is None:
            print("Error: The {0} has to be provided "
                  "if \"{1}\" in --filter_info!".format(info, filters))
            sys.exit()

    def srna_detection(self):
        """sRNA_detection."""
        print("Running sRNA prediction")
        self.check_multi_files(
                [self._args.annotation_files, self._args.transcript_files,
                 self._args.fasta_files, self._args.sorf_files,
                 self._args.terminator_files, self._args.promoter_tables,
                 self._args.processing_site_files],
                ["--annotation_folder", "--transcript_folder",
                 "--fasta_files", "--sorf_files", "--terminator_files",
                 "--promoter_tables", "--processing_site_files"])
        for info in self._args.filter_info:
            if "sec_str" == info:
                self._check_filter_input(
                        self._args.fasta_files, "fasta", "fasta")
                for exe in (self._args.rnafold_path, self._args.relplot_path,
                            self._args.mountain_path, self._args.ps2pdf14_path):
                    self.check_execute_file(exe)
            elif ("blast_nr" == info) or (
                    "blast_srna"== info):
                for exe in (self._args.blastn_path, self._args.blastx_path,
                            self._args.makeblastdb_path):
                    self.check_execute_file(exe)
            elif "sorf" == info:
                self._check_filter_input(
                        self._args.sorf_files, "sORF", "sorf")
            elif "term" == info:
                self._check_filter_input(self._args.terminator_files,
                                         "terminator", "term")
            elif "promoter" == info:
                self._check_filter_input(self._args.promoter_tables,
                                         "Promoter", "promoter")
            elif "tss" == info:
                self._check_filter_input(self._args.tss_files,
                                         "TSS", "tss")
            else:
                print("Error: Please check the --filter_info, "
                      "invalid value was assigned!")
                sys.exit()
        if self._args.utr_derived_srna:
            if self._args.tss_files is None:
                print("Error: The TSS has to be provided "
                      "if you want to compute UTR-derived sRNA!")
                sys.exit()
        project_creator.create_subfolders(self._paths.required_folders("srna"))
        args_srna = self.args_container.container_srna(
                self._args.rnafold_path, self._args.relplot_path,
                self._args.mountain_path, self._args.blastn_path,
                self._args.blastx_path, self._args.makeblastdb_path,
                self._args.ps2pdf14_path, self._paths.srna_folder,
                self._args.utr_derived_srna, self._args.annotation_files,
                self._args.tss_files, self._args.transcript_files,
                self._args.tss_intergenic_fuzzy, self._args.tss_5utr_fuzzy,
                self._args.tss_3utr_fuzzy, self._args.tss_intercds_fuzzy,
                self._args.filter_info, self._args.processing_site_files,
                self._args.fasta_files, self._args.mountain_plot,
                self._args.nr_format, self._args.srna_format,
                self._args.srna_database_path, self._args.nr_database_path,
                self._args.cutoff_energy, self._args.parallel_blast,
                self._args.run_intergenic_tex_coverage,
                self._args.run_intergenic_notex_coverage,
                self._args.run_intergenic_fragmented_coverage,
                self._args.run_break_transcript,
                self._args.run_antisense_tex_coverage,
                self._args.run_antisense_notex_coverage,
                self._args.run_antisense_fragmented_coverage,
                self._args.run_utr_tex_coverage,
                self._args.run_utr_notex_coverage,
                self._args.run_utr_fragmented_coverage,
                self._args.max_length, self._args.min_length,
                self._args.tex_notex_libs, self._args.frag_libs,
                self._args.replicate_tex, self._args.replicate_frag,
                self._args.tex_notex, self._args.blast_e_nr,
                self._args.blast_e_srna, self._args.detect_srna_in_cds,
                self._args.table_best, self._args.decrease_intergenic_antisense,
                self._args.decrease_utr, self._args.fuzzy_intergenic_antisense,
                self._args.fuzzy_utr, self._args.cutoff_nr_hit,
                self._args.sorf_files,
                self._args.overlap_percent_cds,
                self._args.terminator_files,
                self._args.terminator_fuzzy_in_srna,
                self._args.terminator_fuzzy_out_srna,
                self._args.ignore_hypothetical_protein, self._args.tss_source,
                self._args.min_utr_coverage, self._args.promoter_tables,
                self._args.ranking_time_promoter, self._args.promoter_name)
        srna = sRNADetection(args_srna)
        srna.run_srna_detection(args_srna)

    def sorf_detection(self):
        """sORF_detection."""
        print("Running sORF prediction")
        self.check_multi_files(
                [self._args.transcript_files, self._args.annotation_files,
                 self._args.fasta_files, self._args.srna_files,
                 self._args.tss_files],
                ["--transcript_files", "--annotation_files",
                 "--fasta_files", "--srna_files", "--tss_files"])
        project_creator.create_subfolders(
            self._paths.required_folders("sorf"))
        args_sorf = self.args_container.container_sorf(
            self._paths.sorf_folder, self._args.utr_derived_sorf,
            self._args.transcript_files,
            self._args.annotation_files,
            self._args.tss_files, self._args.utr_length,
            self._args.min_length, self._args.max_length,
            self._args.cutoff_intergenic_coverage,
            self._args.cutoff_antisense_coverage,
            self._args.cutoff_5utr_coverage,
            self._args.cutoff_3utr_coverage,
            self._args.cutoff_intercds_coverage,
            self._args.fasta_files, self._args.tex_notex_libs,
            self._args.frag_libs, self._args.tex_notex,
            self._args.replicate_tex, self._args.replicate_frag,
            self._args.table_best, self._args.srna_files,
            self._args.start_codon, self._args.stop_codon,
            self._args.cutoff_background, self._args.fuzzy_rbs,
            self._args.rbs_not_after_tss, self._args.print_all_combination,
            self._args.best_no_srna, self._args.best_no_tss,
            self._args.ignore_hypothetical_protein,
            self._args.min_rbs_distance, self._args.max_rbs_distance)
        sorf = sORFDetection(args_sorf)
        sorf.run_sorf_detection(args_sorf)

    def meme(self):
        """promoter detectopn"""
        print("Running promoter detection")
        self.check_multi_files(
                [self._args.tss_files, self._args.fasta_files],
                ["--tss_files", "--fasta_files"])
        if not self._args.tss_source:
            self.check_multi_files([self._args.annotation_files],
                                   ["--annotation_files"])
        if (self._args.program == "both") or (
                self._args.program == "meme"):
            self.check_execute_file(self._args.meme_path)
        elif (self._args.program == "both") or (
                self._args.program == "glam2"):
            self.check_execute_file(self._args.glam2_path)
        project_creator.create_subfolders(
            self._paths.required_folders("promoter"))
        args_pro = self.args_container.container_promoter(
            self._args.meme_path, self._args.glam2_path,
            self._paths.promoter_output_folder, self._args.tex_libs,
            self._args.tss_files, self._args.fasta_files,
            self._args.num_motif, self._args.nt_before_tss,
            self._args.motif_width, self._args.tss_source,
            self._args.annotation_files, self._args.end_run,
            self._args.combine_all, self._args.e_value,
            self._args.parallel, self._args.program)
        meme = MEME(args_pro)
        meme.run_meme(args_pro)

    def operon(self):
        """operon detection"""
        print("Running operon detection")
        self.check_multi_files(
                [self._args.tss_files, self._args.annotation_files,
                 self._args.transcript_files, self._args.utr5_files,
                 self._args.utr3_files, self._args.term_files],
                ["--tss_files", "--annotation_files",
                 "--transcript_files", "--utr5_files",
                 "--utr3_files", "--term_files"])
        project_creator.create_subfolders(
            self._paths.required_folders("operon"))
        args_op = self.args_container.container_operon(
            self._args.tss_files, self._args.annotation_files,
            self._args.transcript_files, self._args.utr5_files,
            self._args.utr3_files, self._args.term_files,
            self._args.tss_fuzzy, self._args.term_fuzzy,
            self._args.min_length, self._args.statistics,
            self._paths.operon_output_folder, self._args.combine_gff,
            self._paths.operon_statistics_folder)
        operon = OperonDetection(args_op)
        operon.run_operon(args_op)

    def circrna(self):
        """circRNA detection"""
        print("Running circular RNA prediction")
        if self._args.align:
            self.check_execute_file(self._args.segemehl_path)
            if self._args.read_files is None:
                print("Error: The read files are needed for --align.")
                sys.exit()
            self.check_multi_files([self._args.read_files], ["--read_files"])
        for exe in (self._args.testrealign_path, self._args.samtools_path):
            self.check_execute_file(exe)
        self.check_multi_files(
                [self._args.fasta_files, self._args.annotation_files,
                 self._args.bam_files],
                ["--fasta_files", "--annotation_files",
                 "--bam_files"])
        project_creator.create_subfolders(
            self._paths.required_folders("circrna"))
        args_circ = self.args_container.container_circrna(
            self._args.align, self._args.process, self._args.fasta_files,
            self._args.annotation_files, self._args.bam_files,
            self._args.read_files, self._paths.circrna_stat_folder,
            self._args.support_reads, self._args.segemehl_path,
            self._args.testrealign_path, self._args.samtools_path,
            self._args.start_ratio, self._args.end_ratio,
            self._args.ignore_hypothetical_protein,
            self._paths.circrna_output_folder)
        circ = CircRNADetection(args_circ)
        circ.run_circrna(args_circ)

    def goterm(self):
        """Go term discovery"""
        print("Running GO term mapping")
        self.check_multi_files(
                [self._args.annotation_files, self._args.transcript_files],
                ["--annotation_files", "--transcript_files"])
        self.check_file([self._args.uniprot_id, self._args.go_obo,
                         self._args.goslim_obo],
                        ["--uniprot_id", "--go.obo", "--goslim_obo"], True)
        project_creator.create_subfolders(
            self._paths.required_folders("go_term"))
        args_go = self.args_container.container_goterm(
            self._args.annotation_files,
            self._paths.goterm_output_folder, self._args.uniprot_id,
            self._args.go_obo, self._args.goslim_obo,
            self._args.transcript_files)
        goterm = GoTermFinding(args_go)
        goterm.run_go_term(args_go)

    def srna_target(self):
        """sRNA target prediction"""
        print("Running sRNA target prediction")
        self.check_multi_files(
                [self._args.fasta_files, self._args.srna_files,
                 self._args.annotation_files],
                ["--fasta_files", "--srna_files",
                 "--annotation_files"])
        if (self._args.program.lower() == "both"):
            for exe in (self._args.rnaplfold_path, self._args.rnaplex_path,
                        self._args.rnaup_path):
                self.check_execute_file(exe)
        elif self._args.program.lower() == "RNAup":
            self.check_execute_file(self._args.rnaup_path)
        elif self._args.program.lower() == "RNAplex":
            for exe in (self._args.rnaplfold_path, self._args.rnaplex_path):
                self.check_execute_file(exe)
        project_creator.create_subfolders(
            self._paths.required_folders("srna_target"))
        args_tar = self.args_container.container_srna_target(
            self._args.rnaplfold_path, self._args.rnaplex_path,
            self._args.rnaup_path, self._args.annotation_files,
            self._args.fasta_files, self._args.srna_files,
            self._args.query_srna, self._args.program,
            self._args.interaction_length, self._args.window_size_target,
            self._args.span_target, self._args.window_size_srna,
            self._args.span_srna,
            self._args.unstructured_region_rnaplex_target,
            self._args.unstructured_region_rnaplex_srna,
            self._args.unstructured_region_rnaup, self._args.energy_threshold,
            self._args.duplex_distance, self._args.top,
            self._paths.starget_output_folder, self._args.process_rnaplex,
            self._args.process_rnaup, self._args.continue_rnaup,
            self._args.potential_target_start, self._args.potential_target_end,
            self._args.target_feature)
        srnatarget = sRNATargetPrediction(args_tar)
        srnatarget.run_srna_target_prediction(args_tar)

    def snp(self):
        """SNP transcript detection"""
        print("Running SNP/mutations calling")
        self.check_multi_files(
                [self._args.fasta_files, self._args.bam_files],
                ["--fasta_files", "--bam_files"])
        if (self._args.bam_type != "target") and (
                self._args.bam_type != "reference"):
            print("Error: Please assign \"target\" or"
                  " \"reference\" to --bam_type!")
            sys.exit()
        if (self._args.ploidy != "haploid") and (
                self._args.ploidy != "diploid"):
            print("Error: Please assign \"haploid\" or"
                  " \"diploid\" to --chromosome_type!")
        if (self._args.caller != "c") and (
                self._args.caller != "m"):
            print("Error: Please assign \"c\" or"
                  " \"m\" to --caller!")
        self.check_parameter([self._args.sample_number],
                             ["--sample_number"])
        for exe in (self._args.bcftools_path, self._args.samtools_path):
            self.check_execute_file(exe)
        project_creator.create_subfolders(self._paths.required_folders("snp"))
        args_snp = self.args_container.container_snp(
            self._args.samtools_path, self._args.bcftools_path,
            self._args.bam_type, self._args.sample_number,
            self._args.program, self._args.fasta_files,
            self._args.bam_files,
            self._args.quality, self._args.read_depth_range,
            self._paths.snp_output_folder, self._args.indel_fraction,
            self._args.ploidy, self._args.rg_tag, self._args.caller,
            self._args.filter_tag_info, self._args.dp4_cutoff)
        snp = SNPCalling(args_snp)
        snp.run_snp_calling(args_snp)

    def ppi(self):
        """PPI network retrieve"""
        print("Running protein-protein interaction networks prediction")
        self.check_multi_files([self._args.annotation_files],
                               ["--annotation_files"])
        self.check_parameter([self._args.proteinid_strains,
                              self._args.species_string],
                             ["--proteinid_strains", "--species_string"])
        project_creator.create_subfolders(
            self._paths.required_folders("ppi_network"))
        args_ppi = self.args_container.container_ppi(
            self._args.annotation_files, self._args.proteinid_strains,
            self._args.without_strain_pubmed, self._args.species_string,
            self._args.score, self._paths.ppi_output_folder,
            self._args.node_size, self._args.query)
        ppi = PPINetwork(self._paths.ppi_output_folder)
        ppi.retrieve_ppi_network(args_ppi)

    def sublocal(self):
        """Subcellular Localization prediction"""
        print("Running subcellular localization prediction")
        self.check_multi_files(
                [self._args.annotation_files, self._args.fasta_files,
                 self._args.transcript_files],
                ["--annotation_files", "--fasta_files",
                 "--transcript_files"])
        if (self._args.bacteria_type != "positive") and (
                self._args.bacteria_type != "negative"):
            print("Error: Please assign \"positive\" or"
                  " \"negative\" to --bacteria_type!")
            sys.exit()
        self.check_execute_file(self._args.psortb_path)
        project_creator.create_subfolders(
            self._paths.required_folders("subcellular_localization"))
        args_sub = self.args_container.container_sublocal(
            self._args.psortb_path, self._args.annotation_files,
            self._args.fasta_files, self._args.bacteria_type,
            self._args.difference_multi, self._args.merge_to_gff,
            self._paths.sublocal_output_folder, self._args.transcript_files)
        sublocal = SubLocal(args_sub)
        sublocal.run_sub_local(args_sub)

    def ribos(self):
        """riboswitch and RNA thermometer prediction"""
        print("Running riboswitch and RNA thermometer prediction")
        self.check_multi_files(
                [self._args.annotation_files, self._args.fasta_files,
                 self._args.tss_files, self._args.transcript_files],
                ["--annotation_files", "--fasta_files", "--tss_files",
                 "--transcript_files"])
        if (self._args.program == "both"):
            self.check_file([self._args.riboswitch_id, self._args.rfam_path],
                            ["--riboswitch_id", "--rfam_path"], True)
            self.check_file([self._args.rna_thermometer_id, self._args.rfam_path],
                            ["--rna_thermometer_id", "--rfam_path"], True)
            project_creator.create_subfolders(
                    self._paths.required_folders("riboswitch"))
            project_creator.create_subfolders(
                    self._paths.required_folders("thermometer"))
            ribos_path = self._paths.ribos_output_folder
            thermo_path = self._paths.thermo_output_folder
        elif (self._args.program == "thermometer"):
            self.check_file([self._args.rna_thermometer_id, self._args.rfam_path],
                            ["--thermometer_id", "--rfam_path"], True)
            project_creator.create_subfolders(
                    self._paths.required_folders("thermometer"))
            ribos_path = None
            thermo_path = self._paths.thermo_output_folder
        elif (self._args.program == "riboswitch"):
            self.check_file([self._args.riboswitch_id, self._args.rfam_path],
                            ["--riboswitch_id", "--rfam_path"], True)
            project_creator.create_subfolders(
                    self._paths.required_folders("riboswitch"))
            ribos_path = self._paths.ribos_output_folder
            thermo_path = None
        else:
            print("Error: Please assign \"thermometer\", \"riboswitch\" "
                  "or \"both\" in --program!")
            sys.exit()
        self.check_execute_file(self._args.cmscan_path)
        self.check_execute_file(self._args.cmpress_path)
        args_ribo = self.args_container.container_ribos(
            self._args.program, self._args.rna_thermometer_id,
            self._args.cmscan_path, self._args.cmpress_path,
            self._args.riboswitch_id,
            self._args.annotation_files, self._args.fasta_files,
            self._args.tss_files, self._args.transcript_files,
            self._args.rfam_path, ribos_path,
            thermo_path, self._args.e_value,
            self._args.output_all, self._paths.database_folder,
            self._args.fuzzy, self._args.start_codon,
            self._args.min_dist_rbs, self._args.max_dist_rbs,
            self._args.fuzzy_rbs, self._args.utr_length)
        ribos = Ribos(args_ribo)
        ribos.run_ribos(args_ribo)

    def crispr(self):
        """CRISPR prediction"""
        print("Running CRISPR prediction")
        self.check_multi_files(
                [self._args.fasta_files, self._args.annotation_files],
                ["--fasta_files", "--annotation_files"])
        self.check_execute_file(self._args.crt_path)
        project_creator.create_subfolders(
            self._paths.required_folders("crispr"))
        args_cris = self.args_container.container_cris(
            self._args.fasta_files, self._args.annotation_files,
            self._args.crt_path, self._args.window_size,
            self._args.min_number_repeat, self._args.min_length_repeat,
            self._args.Max_length_repeat, self._args.min_length_spacer,
            self._args.Max_length_spacer, self._paths.crispr_output_folder,
            self._args.ignore_hypothetical_protein)
        cris = Crispr(args_cris)
        cris.run_crispr(args_cris)

    def merge(self):
        """Merge all features"""
        print("Merging all features to one gff file")
        merge_folder = os.path.join(self._paths.output_folder,
                                    "merge_all_features")
        self.helper.check_make_folder(merge_folder)
        other_features = self._args.other_features_files_path
        self.check_file([self._args.transcript_file_path] + other_features,
                        ["--transcript_path", "--other_features_path"],
                        False)
        self.check_parameter([self._args.strain_name], ["--strain_name"])
        run_merge(merge_folder, self._args.transcript_file_path,
                  self._args.other_features_files_path, self._args.fuzzy_term,
                  self._args.fuzzy_tss,
                  os.path.join(merge_folder, self._args.strain_name))

    def screen(self):
        """generate screenshot"""
        print("Running screenshot generation")
        self.check_file([self._args.main_gff, self._args.fasta],
                        ["--main_gff", "--fasta"], True)
        if self._args.side_gffs is not None:
            for gff in (self._args.side_gffs):
                gff = gff.strip()
                if not os.path.isfile(gff):
                    print("Error: The --side_gffs no exist!")
                    sys.exit()
        if self._args.output_folder is None:
            print("Error: Please assign --output_folder!")
            sys.exit()
        if (self._args.present != "expand") and (
                self._args.present != "collapse") and (
                self._args.present != "squish"):
            print("Error: Please assign \"expand\" or "
                  "\"collapse\" or \"squish\" to --present!")
            sys.exit()
        args_sc = self.args_container.container_screen(
            self._args.main_gff, self._args.side_gffs,
            self._args.fasta, self._args.height,
            self._args.tex_notex_libs, self._args.frag_libs,
            self._args.present, self._args.output_folder)
        screen = Screen(args_sc)
        screen.screenshot(args_sc)
