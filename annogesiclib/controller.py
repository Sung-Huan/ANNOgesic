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
from annogesiclib.overlap import deal_overlap
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
        if (len(args.__dict__) > 3):
            if not os.path.exists(args.project_path):
                print("Error: --project_path does not exists!")
                sys.exit()
        self._paths = Paths(args.project_path)
        self.args_container = ArgsContainer()
        self.helper = Helper()

    def check_folder(self, folders, flags, log):
        '''Check the emtpy or wrong assigned folder'''
        for folder, flag in zip(folders, flags):
            if folder is None:
                log.write("{0} of {1} is not found. Please check it!".format(
                          folder, flag))
                print("Error: {0} of {1} is not found. Please check it!".format(
                      folder, flag))
                sys.exit()
            else:
                if os.path.exists(folder):
                    if len(os.listdir(folder)) == 0:
                        log.write("{0} is a empty folder!".format(flag))
                        print("Error: {0} is a empty folder!".format(flag))
                        sys.exit()
                else:
                    log.write("{0} of {1} is not found. Please check it!".format(
                               folder, flag))
                    print("Error: {0} of {1} is not found. Please check it!".format(
                          folder, flag))
                    sys.exit()

    def check_multi_files(self, input_files, flags, log):
        if input_files is not None:
            for files, flag in zip(input_files, flags):
                if files is not None:
                    for file_ in files:
                        if not os.path.exists(file_):
                            print("Error: {0} in {1} does "
                                  "not exist!".format(file_, flag))
                            log.write(file_ + " does not exists\n")
                            sys.exit()
                        else:
                            log.write(file_ + " exists\n")

    def check_parameter(self, paras, names, log):
        '''Check the parameter is assigned correct or not'''
        for i in range(len(paras)):
            if paras[i] is None:
                print("Error: {0} can not be None. "
                      "Please check it!".format(names[i]))
                log.write(file_ + " need to be assigned.\n")
                sys.exit()

    def check_execute_file(self, exe, log):
        detect = False
        exe_folder = ""
        log.write("Checking " + exe + "\n")
        if os.path.exists(exe):
            detect = True
            full_exe = os.path.realpath(exe)
            log.write(full_exe + " is found.\n")
        else:
            exes = []
            for folder in os.environ["PATH"].split(":"):
                if os.path.isfile(os.path.join(folder, exe)):
                    exe_folder = folder
                    detect = True
                    full_exe = exe
                    if os.path.join(folder, exe) not in exes:
                        exes.append(os.path.join(folder, exe))
                        log.write(os.path.join(folder, exe) + " is found.\n")
        if not detect:
            if os.path.exists(os.path.realpath(exe)):
                full_exe = os.path.realpath(exe)
                log.write(full_exe + " is found.\n")
            else:
                print("Error: {0} can't be found!".format(exe))
                print("Please assign the correct path!")
                sys.exit()
        if (os.path.isfile(full_exe)) or (
                os.path.isfile(os.path.join(exe_folder, exe))):
            log.write("The execute path is " + os.popen(
                "which " + exe).read())
            return full_exe
        else:
            log.write(full_exe + " is not found.\n")
            print("Error: {0} is not a file!".format(exe))
            sys.exit()

    def check_file(self, files, names, require, log):
        '''Check the path of file'''
        for i in range(len(files)):
            if require:
                if files[i] is None:
                    print("Error: {0} can not be None. "
                          "Please check it!".format(names[i]))
                    log.write(names[i] + " is None.\n")
                    sys.exit()
                else:
                    if not os.path.isfile(files[i]):
                        print("Error: {0} is not found. "
                              "Please check it!".format(files[i]))
                        log.write(files[i] + " does not exist.\n")
                        sys.exit()
            else:
                if files[i] is not None:
                    if not os.path.isfile(files[i]):
                        print("Error: {0} is not found. "
                              "Please check it!".format(files[i]))
                        log.write(files[i] + " does not exist.\n")
                        sys.exit()

    def create_project(self, version):
        """Create a new project."""
        project_creator.create_root_folder(self._args.project_path)
        project_creator.create_subfolders(self._paths.required_folders("root"))
        project_creator.create_version_file(
            self._paths.version_path, version)
        sys.stdout.write("Created folder \"%s\" and required subfolders.\n" % (
            self._args.project_path))

    def get_input(self):
        """Download required files from website."""
        print("Running get input files")
        log = open(os.path.join(self._paths.reference_input_folder, "log.txt"), "w")
        if self._args.ftp_path is None:
            print("Error: Please assign the path for downloading the data!")
            sys.exit()
            annotation_folder = self._paths.ref_annotation_folder
            fasta_folder = self._paths.ref_fasta_folder
        self.helper.check_make_folder(self._paths.ref_annotation_folder)
        self.helper.check_make_folder(self._paths.ref_fasta_folder)
        if self._args.ref_gff is True:
            log.write("Get gff files\n")
            get_file(self._args.ftp_path, self._paths.ref_annotation_folder,
                     "gff", log)
            get_file(self._args.ftp_path, self._paths.ref_annotation_folder,
                     "_genomic.gff.gz", log)
        if self._args.ref_fasta is True:
            log.write("Get fasta files\n")
            get_file(self._args.ftp_path, self._paths.ref_fasta_folder,
                     "fna", log)
            get_file(self._args.ftp_path, self._paths.ref_fasta_folder,
                     "_genomic.fna.gz", log)
        if self._args.ref_gbk is True:
            log.write("Get gbk files\n")
            get_file(self._args.ftp_path, self._paths.ref_annotation_folder,
                     "gbk", log)
            get_file(self._args.ftp_path, self._paths.ref_annotation_folder,
                     "gbff", log)
            get_file(self._args.ftp_path, self._paths.ref_annotation_folder,
                     "_genomic.gbff.gz", log)
        if self._args.ref_ptt is True:
            log.write("Get ptt files\n")
            get_file(self._args.ftp_path, self._paths.ref_annotation_folder,
                     "ptt", log)
        if self._args.ref_rnt is True:
            log.write("Get rnt files\n")
            get_file(self._args.ftp_path, self._paths.ref_annotation_folder,
                     "rnt", log)
        if self._args.convert_embl is True:
            annotation_files = os.listdir(self._paths.ref_annotation_folder)
            if len(annotation_files) == 0:
                sys.stdout.write("No gff files!!\n")
            else:
                log.write("Running converter.py for converting gbk file "
                          "to embl formet\n")
                Converter().convert_gbk2embl(self._paths.ref_annotation_folder)
        log.close()

    def get_target_fasta(self):
        """Get target fasta"""
        print("Running update genome fasta")
        project_creator.create_subfolders(
            self._paths.required_folders("get_target_fasta"))
        log = open(os.path.join(self._paths.target_folder, "log.txt"), "w")
        self.check_multi_files([self._args.related_fasta_files],
                               ["--related_fasta_files"], log)
        self.check_file([self._args.mutation_table], ["--mutation_table"], True, log)
        target = TargetFasta(self._paths.tar_fasta_folder,
                             self._args.related_fasta_files)
        target.get_target_fasta(
                self._args.mutation_table, self._paths.tar_fasta_folder,
                self._args.related_fasta_files, self._args.updated_seq_name,
                self._paths.target_base_folder, log)
        log.close()

    def ratt(self):
        """Run RATT to transfer annotation file from reference to target."""
        print("Running annotation transfer")
        project_creator.create_subfolders(
            self._paths.required_folders("get_target_fasta"))
        project_creator.create_subfolders(
            self._paths.required_folders("annotation_transfer"))
        log = open(os.path.join(self._paths.ratt_folder, "log.txt"), "w")
        if (self._args.transfer_type != "Strain") and (
                self._args.transfer_type != "Assembly") and (
                self._args.transfer_type != "Species") and (
                self._args.transfer_type != "Assembly.Repetitive") and (
                self._args.transfer_type != "Strain.Repetitive") and (
                self._args.transfer_type != "Species.Repetitive") and (
                self._args.transfer_type != "Multiple") and (
                self._args.transfer_type != "Free"):
            log.write("Incorrect --transfer_type. Please assign 'Assembly', 'Species', "
                      "'Assembly.Repetitive', 'Strain.Repetitive', 'Species.Repetitive', "
                      "'Multiple' or 'Free'\n")
            print("Error: please assign correct --transfer_type!")
            sys.exit()
        if (self._args.related_embl_files is None) and (
                self._args.related_gbk_files is None):
            print("Error: please assign proper embl or genbank files")
            log.write("--related_gbk_files and --related_embl_files can not be both None.\n")
            sys.exit()
        elif (self._args.related_embl_files is not None) and (
                self._args.related_gbk_files is not None):
            log.write("Please choose --related_gbk_files as input or "
                      "--related_embl_files as input. Do not assign both.\n")
            print("Error: please choose embl as input or genbank as input")
            sys.exit()
        self._args.ratt_path = self.check_execute_file(self._args.ratt_path, log)
        self.check_multi_files(
                [self._args.target_fasta_files, self._args.related_fasta_files],
                ["--target_fasta_files", "--closed_fasta_files"], log)
        self.check_parameter([self._args.element, self._args.compare_pair],
                             ["--element", "--compare_pair"], log)
        args_ratt = self.args_container.container_ratt(
            self._args.ratt_path, self._args.element, self._args.transfer_type,
            self._args.related_embl_files, self._args.related_gbk_files,
            self._args.target_fasta_files, self._args.related_fasta_files,
            self._paths.ratt_folder,
            self._paths.tar_annotation_folder, self._args.compare_pair)
        ratt = RATT(args_ratt)
        ratt.annotation_transfer(args_ratt, log)
        log.close()

    def tsspredator(self):
        """Run TSSpredator for predicting TSS candidates."""
        if self._args.program.lower() == "tss":
            print("Running TSS prediction")
            project_creator.create_subfolders(
                self._paths.required_folders("TSS"))
            out_folder = self._paths.tsspredator_folder
            log = open(os.path.join(self._paths.tsspredator_folder,
                                    "log.txt"), "w")
            log.write("Running TSS prediction.\n")
        elif self._args.program.lower() == "ps":
            print("Running processing site prediction")
            out_folder = self._paths.processing_site_folder
            project_creator.create_subfolders(
                self._paths.required_folders("processing"))
            log = open(os.path.join(self._paths.processing_site_folder,
                                    "log.txt"), "w")
            log.write("Running PS prediction.\n")
        else:
            print("Error: No such program!")
            sys.exit()
        self.check_multi_files(
                [self._args.fasta_files, self._args.annotation_files,
                 self._args.compare_overlap_gff, self._args.manual_files,
                 self._args.compare_transcript_files],
                ["--fasta_files", "--annotation_files", "--compare_overlap_gff",
                 "--manual_files","--compare_transcript_files"], log)
        self.check_parameter([self._args.tex_notex_libs, self._args.condition_names],
                             ["--tex_notex_libs", "--condition_names"], log)
        self._args.tsspredator_path = self.check_execute_file(
                self._args.tsspredator_path, log)
        args_tss = self.args_container.container_tsspredator(
            self._args.tsspredator_path, self._args.program,
            self._args.fasta_files, self._args.annotation_files,
            self._args.tex_notex_libs, self._args.condition_names,
            self._args.height, self._args.height_reduction,
            self._args.factor, self._args.factor_reduction,
            self._args.base_height, self._args.enrichment_factor,
            self._args.processing_factor, self._args.replicate_tex,
            out_folder, self._args.validate_gene,
            self._args.manual_files, self._args.curated_sequence_length,
            self._args.compare_transcript_files, self._args.tolerance,
            self._args.utr_length, self._args.cluster,
            self._args.re_check_orphan,
            self._args.remove_overlap_feature, self._args.compare_overlap_gff,
            self._args.remove_low_expression)
        tsspredator = TSSpredator(args_tss)
        tsspredator.run_tsspredator(args_tss, log)
        log.close()

    def optimize(self):
        """opimize TSSpredator"""
        if self._args.program.lower() == "tss":
            print("Running optimization of TSS prediction")
            project_creator.create_subfolders(
                self._paths.required_folders("TSS"))
            out_folder = self._paths.tsspredator_folder
            if "optimized_TSSpredator" not in os.listdir(out_folder):
                os.mkdir(os.path.join(out_folder, "optimized_TSSpredator"))
            log = open(os.path.join(out_folder, "optimized_TSSpredator",
                                    "log.txt"), "w")
            log.write("Running optimization of TSS prediction\n")
        elif self._args.program.lower() == "ps":
            out_folder = self._paths.processing_site_folder
            project_creator.create_subfolders(
                self._paths.required_folders("processing"))
            if "optimized_TSSpredator" not in os.listdir(out_folder):
                os.mkdir(os.path.join(out_folder, "optimized_TSSpredator"))
            log = open(os.path.join(out_folder, "optimized_TSSpredator",
                                    "log.txt"), "w")
            log.write("Running optimization of PS prediction\n")
            print("Running optimization of processing site prediction")
        else:
            print("Error: No such program!")
            sys.exit()
        self.check_multi_files(
                [self._args.fasta_files, self._args.annotation_files,
                 self._args.manual_files],
                ["--fasta_files", "--annotation_files", "--manual_files"], log)
        self._args.tsspredator_path = self.check_execute_file(
                self._args.tsspredator_path, log)
        self.check_parameter([self._args.tex_notex_libs,
                              self._args.condition_names],
                             ["--tex_notex_lib",
                              "--condition_names"], log)
        args_ops = self.args_container.container_optimize(
            self._args.tsspredator_path, self._args.fasta_files,
            self._args.annotation_files,
            self._args.manual_files, out_folder, self._args.max_height,
            self._args.max_height_reduction, self._args.max_factor,
            self._args.max_factor_reduction, self._args.max_base_height,
            self._args.max_enrichment_factor, self._args.max_processing_factor,
            self._args.utr_length, self._args.tex_notex_libs,
            self._args.condition_names, self._args.cluster,
            self._args.curated_sequence_length, self._args.parallels,
            self._args.program, self._args.replicate_tex,
            self._args.steps)
        optimize_tss(args_ops, log)
        log.close()

    def color(self):
        """color the screenshots"""
        print("Running png files coloring")
        if not os.path.exists(os.path.join(self._args.screenshot_folder,
                                           "screenshots")):
            print("The folder -- screenshots needs to be found in "
                  "{0}.".format(self._args.screenshot_folder))
            sys.exit()
        log = open(os.path.join(self._args.screenshot_folder, "screenshots",
                                "color_log.txt"), "w")
        self.check_parameter([self._args.track_number], ["--track_numer"], log)
        self.check_folder([self._args.screenshot_folder], ["--screenshot_folder"], log)
        self._args.imagemagick_covert_path = self.check_execute_file(
                self._args.imagemagick_covert_path, log)
        color = ColorPNG()
        color.generate_color_png(
                self._args.track_number, self._args.screenshot_folder,
                self._args.imagemagick_covert_path, log)
        log.close()

    def terminator(self):
        """Run TransTermHP and Gene converaged for detecting terminators"""
        print("Running terminator prediction")
        project_creator.create_subfolders(
            self._paths.required_folders("terminator"))
        log = open(os.path.join(self._paths.transterm_folder, "log.txt"), "w")
        if self._args.transterm_path is None:
            print("Please assign the path of transterm in TransTermHP.")
        self.check_multi_files(
                [self._args.fasta_files, self._args.annotation_files,
                 self._args.transcript_files, self._args.srna_files],
                ["--fasta_files", "--annotation_files",
                 "--transcript_files", "--srna_files"], log)
        for prop in ("transterm_path", "expterm_path", "rnafold_path"):
            setattr(self._args, prop,
                    self.check_execute_file(getattr(self._args, prop), log))
        args_term = self.args_container.container_terminator(
            self._args.transterm_path, self._args.expterm_path,
            self._args.rnafold_path,
            self._paths.transterm_folder, self._args.fasta_files,
            self._args.annotation_files, self._args.transcript_files,
            self._args.srna_files, self._args.decrease,
            self._args.highest_coverage, self._args.tolerance_detect_coverage,
            self._args.tolerance_within_transcript,
            self._args.tolerance_downstream_transcript,
            self._args.tolerance_within_gene,
            self._args.tolerance_downstream_gene, self._paths.transtermhp_folder,
            self._args.tex_notex_libs, self._args.frag_libs,
            self._args.tex_notex, self._args.replicate_tex,
            self._args.replicate_frag,
            self._args.min_loop_length, self._args.max_loop_length,
            self._args.min_stem_length, self._args.max_stem_length,
            self._args.min_u_tail, self._args.miss_rate,
            self._args.mutation_u_tail, self._args.keep_multi_term,
            self._args.window_size, self._args.window_shift)
        terminator = Terminator(args_term)
        terminator.run_terminator(args_term, log)
        log.close()

    def transcript(self):
        """Run Transcript detection"""
        project_creator.create_subfolders(
            self._paths.required_folders("transcript"))
        log = open(os.path.join(self._paths.transcript_output_folder, "log.txt"), "w")
        print("Running transcript detection")
        self.check_multi_files(
                [self._args.annotation_files, self._args.tss_files,
                 self._args.terminator_files],
                ["--annotation_files", "--tss_files", "--terminator_files"], log)
        args_tran = self.args_container.container_transcript(
            self._args.tex_notex, self._args.modify_transcript,
            self._args.length, self._args.annotation_files,
            self._args.height, self._args.width,
            self._args.tolerance, self._args.tolerance_coverage,
            self._args.replicate_tex, self._args.replicate_frag,
            self._paths.transcript_output_folder,
            self._args.tss_files, self._args.tss_tolerance,
            self._args.tex_notex_libs, self._args.frag_libs,
            self._args.compare_feature_genome,
            self._args.terminator_files, self._args.terminator_tolerance,
            self._args.max_length_distribution)
        transcript = TranscriptDetection(args_tran)
        transcript.run_transcript(args_tran, log)

    def utr_detection(self):
        """Run UTR detection."""
        print("Running UTR detection")
        project_creator.create_subfolders(self._paths.required_folders("utr"))
        log = open(os.path.join(self._paths.utr_folder, "log.txt"), "w")
        self.check_multi_files(
            [self._args.annotation_files, self._args.terminator_files,
             self._args.transcript_files, self._args.tss_files],
            ["--annotation_files", "--terminator_files",
             "--transcript_files", "--tss_files"], log)
        args_utr = self.args_container.container_utr(
                self._args.tss_files, self._args.annotation_files,
                self._args.transcript_files, self._args.terminator_files,
                self._args.terminator_tolerance, self._paths.utr_folder,
                self._args.tss_source, self._args.base_5utr,
                self._args.utr_length, self._args.base_3utr)
        utr = UTRDetection(args_utr)
        utr.run_utr_detection(args_utr, log)

    def _check_filter_input(self, files, info, filters, log):
        if files is None:
            print("Error: The {0} has to be provided "
                  "if \"{1}\" in --filter_info!".format(info, filters))
            log.write("The {0} has to be provided "
                      "if \"{1}\" in --filter_info!\n".format(info, filters))
            sys.exit()

    def _check_database(self, database, flag, info, log):
        wrong = False
        if database is None:
            wrong = True
        elif not os.path.isfile(database):
            if (os.path.isfile(database + ".fa")) or (
                    os.path.isfile(database + ".fna")) or (
                    os.path.isfile(database + ".fasta")):
                wrong = False
            else:
                wrong = True
        if wrong:
            print("Error: {0} is required if {1} is in --filter_info. "
                  "But the assignment of {0} is empty or wrong. "
                  "Please check the {0} or remove {1} from "
                  "--filter_info!".format(flag, info))
            log.write("{0} is required if {1} is in --filter_info. "
                      "But the assignment of {0} is empty or wrong. "
                      "Please check the {0} or remove {1} from "
                      "--filter_info!\n".format(flag, info))
            sys.exit()

    def srna_detection(self):
        """sRNA_detection."""
        print("Running sRNA prediction")
        project_creator.create_subfolders(self._paths.required_folders("srna"))
        log = open(os.path.join(self._paths.srna_folder, "log.txt"), "w")
        self.check_multi_files(
                [self._args.annotation_files, self._args.transcript_files,
                 self._args.fasta_files, self._args.sorf_files,
                 self._args.terminator_files, self._args.promoter_tables,
                 self._args.processing_site_files],
                ["--annotation_files", "--transcript_files",
                 "--fasta_files", "--sorf_files", "--terminator_files",
                 "--promoter_tables", "--processing_site_files"], log)
        for info in self._args.filter_info:
            if "sec_str" == info:
                if not self._args.compute_sec_structures:
                    log.write("If you want to use secondary structure to "
                              "filter the false positive, "
                              "--compute_sec_structure need to be switch on.\n")
                    print("Error: --compute_sec_structures is not switch on, "
                          "but sec_str is still in --filter_info.")
                    sys.exit()
                self._check_filter_input(
                        self._args.fasta_files, "fasta file", "sec_str", log)
                for prop in ("rnafold_path", "relplot_path",
                             "mountain_path"):
                    setattr(self._args, prop,
                            self.check_execute_file(getattr(self._args, prop), log))
            elif ("blast_nr" == info) or (
                    "blast_srna"== info):
                for prop in ("blastn_path", "blastx_path", "makeblastdb_path"):
                    setattr(self._args, prop,
                            self.check_execute_file(getattr(self._args, prop), log))
                if ("blast_nr" == info):
                    self._check_database(self._args.nr_database_path,
                                         "--nr_database_path", "blast_nr", log)
                if ("blast_srna" == info):
                    self._check_database(self._args.srna_database_path,
                                         "--srna_database_path", "blast_srna", log)
            elif "sorf" == info:
                self._check_filter_input(
                        self._args.sorf_files, "sORF", "sorf", log)
            elif "term" == info:
                self._check_filter_input(self._args.terminator_files,
                                         "terminator", "term", log)
            elif "promoter" == info:
                self._check_filter_input(self._args.promoter_tables,
                                         "Promoter", "promoter", log)
            elif "tss" == info:
                self._check_filter_input(self._args.tss_files,
                                         "TSS", "tss", log)
            else:
                if "none" != info.lower():
                    print("Error: Please check the --filter_info, "
                          "invalid value was assigned!")
                    log.write("invalid value was assigned to --filter_info.\n")
                    sys.exit()
        log.write("--filter_info and databases are assigned correctly.\n")
        if self._args.utr_derived_srna:
            if self._args.tss_files is None:
                print("Error: The TSS has to be provided "
                      "if you want to compute UTR-derived sRNA!")
                log.write("The TSS has to be provided "
                          "if you want to compute UTR-derived sRNA!\n")
                sys.exit()
        if self._args.search_poly_u != 0:
            if self._args.fasta_files is None:
                print("Error: The fasta files have to be provided "
                      "if you want to extend 3'end of sRNA by "
                      "searching poly U tail!")
                log.write("The fasta files have to be provided "
                          "if you want to extend 3'end of sRNA by "
                          "searching poly U tail!\n")
                sys.exit()
        args_srna = self.args_container.container_srna(
                self._args.rnafold_path, self._args.relplot_path,
                self._args.mountain_path, self._args.blastn_path,
                self._args.blastx_path, self._args.makeblastdb_path,
                self._paths.srna_folder, self._args.utr_derived_srna,
                self._args.annotation_files, self._args.tss_files,
                self._args.transcript_files,
                self._args.tss_intergenic_antisense_tolerance,
                self._args.tss_5utr_tolerance, self._args.tss_3utr_tolerance,
                self._args.tss_intercds_tolerance, self._args.filter_info,
                self._args.processing_site_files, self._args.fasta_files,
                self._args.mountain_plot, self._args.nr_format,
                self._args.srna_format, self._args.srna_database_path,
                self._args.nr_database_path, self._args.cutoff_energy,
                self._args.parallel_blast, self._args.blast_score_srna,
                self._args.blast_score_nr,
                self._args.min_intergenic_tex_coverage,
                self._args.min_intergenic_notex_coverage,
                self._args.min_intergenic_fragmented_coverage,
                self._args.min_complete_5utr_transcript_coverage,
                self._args.min_antisense_tex_coverage,
                self._args.min_antisense_notex_coverage,
                self._args.min_antisense_fragmented_coverage,
                self._args.min_utr_tex_coverage,
                self._args.min_utr_notex_coverage,
                self._args.min_utr_fragmented_coverage,
                self._args.max_length, self._args.min_length,
                self._args.tex_notex_libs, self._args.frag_libs,
                self._args.replicate_tex, self._args.replicate_frag,
                self._args.tex_notex, self._args.blast_e_nr,
                self._args.blast_e_srna, self._args.detect_srna_in_cds,
                self._args.decrease_intergenic_antisense,
                self._args.decrease_utr, self._args.tolerance_intergenic_antisense,
                self._args.tolerance_utr, self._args.cutoff_nr_hit,
                self._args.sorf_files, self._args.overlap_percent_cds,
                self._args.terminator_files,
                self._args.terminator_tolerance_in_srna,
                self._args.terminator_tolerance_out_srna,
                self._args.ignore_hypothetical_protein, self._args.tss_source,
                self._args.min_all_utr_coverage, self._args.promoter_tables,
                self._args.ranking_time_promoter, self._args.promoter_names,
                self._args.compute_sec_structures, self._args.search_poly_u,
                self._args.min_u_poly_u, self._args.mutation_poly_u,
                self._args.exclude_srna_in_annotation_file)
        srna = sRNADetection(args_srna)
        srna.run_srna_detection(args_srna, log)

    def sorf_detection(self):
        """sORF_detection."""
        print("Running sORF prediction")
        project_creator.create_subfolders(
            self._paths.required_folders("sorf"))
        log = open(os.path.join(self._paths.sorf_folder, "log.txt"), "w")
        self.check_multi_files(
                [self._args.transcript_files, self._args.annotation_files,
                 self._args.fasta_files, self._args.srna_files,
                 self._args.tss_files],
                ["--transcript_files", "--annotation_files",
                 "--fasta_files", "--srna_files", "--tss_files"], log)
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
            self._args.srna_files,
            self._args.start_codon, self._args.stop_codon,
            self._args.cutoff_base_coverage, self._args.rbs_seq,
            self._args.tolerance_rbs,
            self._args.rbs_not_after_tss, self._args.print_all_combination,
            self._args.best_no_srna, self._args.best_no_tss,
            self._args.ignore_hypothetical_protein,
            self._args.min_rbs_distance, self._args.max_rbs_distance,
            self._args.tolerance_3end, self._args.tolerance_5end,
            self._args.contain_multi_stop)
        sorf = sORFDetection(args_sorf)
        sorf.run_sorf_detection(args_sorf, log)

    def meme(self):
        """promoter detectopn"""
        print("Running promoter detection")
        project_creator.create_subfolders(
            self._paths.required_folders("promoter"))
        log = open(os.path.join(self._paths.promoter_output_folder, "log.txt"), "w")
        self.check_multi_files(
                [self._args.tss_files, self._args.fasta_files],
                ["--tss_files", "--fasta_files"], log)
        if not self._args.tss_source:
            self.check_multi_files([self._args.annotation_files],
                                   ["--annotation_files"], log)
        if (self._args.program == "both") or (
                self._args.program == "meme"):
            self._args.meme_path = self.check_execute_file(self._args.meme_path, log)
        elif (self._args.program == "both") or (
                self._args.program == "glam2"):
            self._args.glam2_path = self.check_execute_file(self._args.glam2_path, log)
        args_pro = self.args_container.container_promoter(
            self._args.meme_path, self._args.glam2_path,
            self._paths.promoter_output_folder, self._args.tex_libs,
            self._args.tss_files, self._args.fasta_files,
            self._args.num_motifs, self._args.nt_before_tss,
            self._args.motif_width, self._args.tss_source,
            self._args.annotation_files, self._args.end_run,
            self._args.combine_all, self._args.e_value,
            self._args.parallels, self._args.program, self._args.use_tss_type)
        meme = MEME(args_pro)
        meme.run_meme(args_pro, log)

    def operon(self):
        """operon detection"""
        print("Running operon detection")
        project_creator.create_subfolders(
            self._paths.required_folders("operon"))
        log = open(os.path.join(self._paths.operon_output_folder, "log.txt"), "w")
        self.check_multi_files(
                [self._args.tss_files, self._args.annotation_files,
                 self._args.transcript_files, self._args.terminator_files],
                ["--tss_files", "--annotation_files",
                 "--transcript_files", "--terminator_files"], log)
        args_op = self.args_container.container_operon(
            self._args.tss_files, self._args.annotation_files,
            self._args.transcript_files, self._args.terminator_files,
            self._args.tss_tolerance, self._args.terminator_tolerance,
            self._args.min_length, self._paths.operon_output_folder,
            self._paths.operon_statistics_folder)
        operon = OperonDetection(args_op)
        operon.run_operon(args_op, log)

    def circrna(self):
        """circRNA detection"""
        print("Running circular RNA prediction")
        project_creator.create_subfolders(
            self._paths.required_folders("circrna"))
        log = open(os.path.join(self._paths.circrna_output_folder,
                                "log.txt"), "w")
        if self._args.read_files:
            self._args.segemehl_path = self.check_execute_file(
                    self._args.segemehl_path, log)
        for prop in ("testrealign_path", "samtools_path"):
            setattr(self._args, prop,
                    self.check_execute_file(getattr(self._args, prop), log))
        self.check_multi_files(
                [self._args.fasta_files, self._args.annotation_files],
                ["--fasta_files", "--annotation_files"], log)
        args_circ = self.args_container.container_circrna(
            self._args.parallels, self._args.fasta_files,
            self._args.annotation_files, self._args.bam_files,
            self._args.read_files, self._paths.circrna_stat_folder,
            self._args.support_reads, self._args.segemehl_path,
            self._args.testrealign_path, self._args.samtools_path,
            self._args.start_ratio, self._args.end_ratio,
            self._args.ignore_hypothetical_protein,
            self._paths.circrna_output_folder)
        circ = CircRNADetection(args_circ)
        circ.run_circrna(args_circ, log)

    def goterm(self):
        """Go term discovery"""
        print("Running GO term mapping")
        project_creator.create_subfolders(
            self._paths.required_folders("go_term"))
        log = open(os.path.join(self._paths.goterm_output_folder, "log.txt"), "w")
        self.check_multi_files(
                [self._args.annotation_files, self._args.transcript_files],
                ["--annotation_files", "--transcript_files"], log)
        self.check_file([self._args.uniprot_id, self._args.go_obo,
                         self._args.goslim_obo],
                        ["--uniprot_id", "--go.obo", "--goslim_obo"], True, log)
        args_go = self.args_container.container_goterm(
            self._args.annotation_files,
            self._paths.goterm_output_folder, self._args.uniprot_id,
            self._args.go_obo, self._args.goslim_obo,
            self._args.transcript_files)
        goterm = GoTermFinding(args_go)
        goterm.run_go_term(args_go, log)

    def srna_target(self):
        """sRNA target prediction"""
        print("Running sRNA target prediction")
        project_creator.create_subfolders(
            self._paths.required_folders("srna_target"))
        log = open(os.path.join(self._paths.starget_output_folder,
                                "log.txt"), "w")
        self.check_multi_files(
                [self._args.fasta_files, self._args.srna_files,
                 self._args.annotation_files],
                ["--fasta_files", "--srna_files",
                 "--annotation_files"], log)
        if "RNAup" in self._args.program:
            self._args.rnaup_path = self.check_execute_file(
                    self._args.rnaup_path, log)
        if "RNAplex" in self._args.program:
            for prop in ("rnaplfold_path", "rnaplex_path"):
                setattr(self._args, prop,
                        self.check_execute_file(getattr(self._args, prop), log))
        if "IntaRNA" in self._args.program:
            self._args.intarna_path = self.check_execute_file(
                    self._args.intarna_path, log)
            if self._args.mode_intarna is None:
                print("Error: --mode_IntaRNA need to be assigned!")
                sys.exit()
        args_tar = self.args_container.container_srna_target(
            self._args.rnaplfold_path, self._args.rnaplex_path,
            self._args.rnaup_path, self._args.intarna_path,
            self._args.annotation_files,
            self._args.fasta_files, self._args.srna_files,
            self._args.query_srnas, self._args.program,
            self._args.interaction_length,
            self._args.window_size_target_rnaplex,
            self._args.span_target_rnaplex,
            self._args.window_size_srna_rnaplfold,
            self._args.span_srna_rnaplfold,
            self._args.unstructured_region_rnaplex_target,
            self._args.unstructured_region_rnaplex_srna,
            self._args.unstructured_region_rnaup,
            self._args.energy_threshold_rnaplex,
            self._args.duplex_distance_rnaplex, self._args.top,
            self._paths.starget_output_folder, self._args.parallels_rnaplex,
            self._args.parallels_rnaup, self._args.parallels_intarna,
            self._args.continue_rnaup,
            self._args.slide_window_size_srna_intarna,
            self._args.max_loop_length_srna_intarna,
            self._args.slide_window_size_target_intarna,
            self._args.max_loop_length_target_intarna,
            self._args.mode_intarna, self._args.potential_target_start,
            self._args.potential_target_end, self._args.target_feature)
        srnatarget = sRNATargetPrediction(args_tar)
        srnatarget.run_srna_target_prediction(args_tar, log)

    def snp(self):
        """SNP transcript detection"""
        print("Running SNP/mutations calling")
        project_creator.create_subfolders(self._paths.required_folders("snp"))
        log = open(os.path.join(self._paths.snp_output_folder,
                                "log.txt"), "w")
        self.check_multi_files(
                [self._args.fasta_files],
                ["--fasta_files"], log)
        if (self._args.bam_type != "related_genome") and (
                self._args.bam_type != "reference_genome"):
            print("Error: Please assign \"related_genome\" or"
                  " \"reference_genome\" to --bam_type!")
            sys.exit()
        if (self._args.ploidy != "haploid") and (
                self._args.ploidy != "diploid"):
            print("Error: Please assign \"haploid\" or"
                  " \"diploid\" to --chromosome_type!")
        if (self._args.caller != "c") and (
                self._args.caller != "m"):
            print("Error: Please assign \"c\" or"
                  " \"m\" to --caller!")
        for prop in ("bcftools_path", "samtools_path"):
            setattr(self._args, prop,
                    self.check_execute_file(getattr(self._args, prop), log))
        args_snp = self.args_container.container_snp(
            self._args.samtools_path, self._args.bcftools_path,
            self._args.bam_type,
            self._args.program, self._args.fasta_files,
            self._args.bam_files,
            self._args.quality, self._args.read_depth_range,
            self._paths.snp_output_folder, self._args.indel_fraction,
            self._args.ploidy, self._args.rg_tag, self._args.caller,
            self._args.filter_tag_info, self._args.dp4_cutoff)
        snp = SNPCalling(args_snp)
        snp.run_snp_calling(args_snp, log)

    def ppi(self):
        """PPI network retrieve"""
        project_creator.create_subfolders(
            self._paths.required_folders("ppi_network"))
        log = open(os.path.join(self._paths.ppi_output_folder,
                                "log.txt"), "w")
        print("Running protein-protein interaction networks prediction")
        self.check_multi_files([self._args.annotation_files],
                               ["--annotation_files"], log)
        self.check_parameter([self._args.query_strains,
                              self._args.species_string],
                             ["--query_strains", "--species_string"], log)
        args_ppi = self.args_container.container_ppi(
            self._args.annotation_files, self._args.query_strains,
            self._args.without_strain_pubmed, self._args.species_string,
            self._args.score, self._paths.ppi_output_folder,
            self._args.node_size, self._args.query)
        ppi = PPINetwork(self._paths.ppi_output_folder)
        ppi.retrieve_ppi_network(args_ppi, log)

    def sublocal(self):
        """Subcellular Localization prediction"""
        print("Running subcellular localization prediction")
        project_creator.create_subfolders(
            self._paths.required_folders("subcellular_localization"))
        log = open(os.path.join(self._paths.sublocal_output_folder,
                                "log.txt"), "w")
        self.check_multi_files(
                [self._args.annotation_files, self._args.fasta_files,
                 self._args.transcript_files],
                ["--annotation_files", "--fasta_files",
                 "--transcript_files"], log)
        if (self._args.bacteria_type != "positive") and (
                self._args.bacteria_type != "negative"):
            print("Error: Please assign \"positive\" or"
                  " \"negative\" to --bacteria_type!")
            sys.exit()
        self._args.psortb_path = self.check_execute_file(self._args.psortb_path, log)
        args_sub = self.args_container.container_sublocal(
            self._args.psortb_path, self._args.annotation_files,
            self._args.fasta_files, self._args.bacteria_type,
            self._args.difference_multi,
            self._paths.sublocal_output_folder, self._args.transcript_files)
        sublocal = SubLocal(args_sub)
        sublocal.run_sub_local(args_sub, log)

    def ribos(self):
        """riboswitch and RNA thermometer prediction"""
        print("Running riboswitch and RNA thermometer prediction")
        log_t = None
        log_r = None
        if (self._args.program == "both"):
            project_creator.create_subfolders(
                    self._paths.required_folders("riboswitch"))
            project_creator.create_subfolders(
                    self._paths.required_folders("thermometer"))
            log_r = open(os.path.join(self._paths.ribos_output_folder,
                                "log.txt"), "w")
            log_t = open(os.path.join(self._paths.thermo_output_folder,
                                "log.txt"), "w")
            self.check_file([self._args.riboswitch_id_file, self._args.rfam_path],
                            ["--riboswitch_id_file", "--rfam_path"], True, log_r)
            self.check_file([self._args.rna_thermometer_id_file,
                             self._args.rfam_path],
                            ["--rna_thermometer_id_file", "--rfam_path"], True, log_t)
            ribos_path = self._paths.ribos_output_folder
            thermo_path = self._paths.thermo_output_folder
        elif (self._args.program == "thermometer"):
            project_creator.create_subfolders(
                    self._paths.required_folders("thermometer"))
            log_t = open(os.path.join(self._paths.thermo_output_folder,
                                "log.txt"), "w")
            self.check_file([self._args.rna_thermometer_id_file,
                             self._args.rfam_path],
                            ["--thermometer_id_file", "--rfam_path"], True, log_t)
            ribos_path = None
            thermo_path = self._paths.thermo_output_folder
        elif (self._args.program == "riboswitch"):
            project_creator.create_subfolders(
                    self._paths.required_folders("riboswitch"))
            log_r = open(os.path.join(self._paths.ribos_output_folder,
                                "log.txt"), "w")
            self.check_file([self._args.riboswitch_id_file, self._args.rfam_path],
                            ["--riboswitch_id_file", "--rfam_path"], True, log_r)
            ribos_path = self._paths.ribos_output_folder
            thermo_path = None
        else:
            log.write("Please assign \"thermometer\", \"riboswitch\" "
                      "or \"both\" in --program.\n")
            print("Error: Please assign \"thermometer\", \"riboswitch\" "
                  "or \"both\" in --program!")
            sys.exit()
        for log in (log_t, log_r):
            if log is not None:
                self.check_multi_files(
                        [self._args.annotation_files, self._args.fasta_files,
                         self._args.tss_files, self._args.transcript_files],
                        ["--annotation_files", "--fasta_files", "--tss_files",
                         "--transcript_files"], log)
                self._args.cmscan_path = self.check_execute_file(self._args.cmscan_path, log)
                self._args.cmpress_path = self.check_execute_file(self._args.cmpress_path, log)
        args_ribo = self.args_container.container_ribos(
            self._args.program, self._args.rna_thermometer_id_file,
            self._args.cmscan_path, self._args.cmpress_path,
            self._args.riboswitch_id_file,
            self._args.annotation_files, self._args.fasta_files,
            self._args.tss_files, self._args.transcript_files,
            self._args.rfam_path, ribos_path,
            thermo_path, self._args.cutoff,
            self._args.output_all, self._paths.database_folder,
            self._args.tolerance, self._args.without_rbs,
            self._args.rbs_seq,
            self._args.tolerance_rbs, self._args.utr_length)
        ribos = Ribos(args_ribo)
        ribos.run_ribos(args_ribo, log_t, log_r)

    def crispr(self):
        """CRISPR prediction"""
        print("Running CRISPR prediction")
        project_creator.create_subfolders(
            self._paths.required_folders("crispr"))
        log = open(os.path.join(self._paths.crispr_output_folder,
                                "log.txt"), "w")
        self.check_multi_files(
                [self._args.fasta_files, self._args.annotation_files],
                ["--fasta_files", "--annotation_files"], log)
        self._args.crt_path = self.check_execute_file(self._args.crt_path, log)
        args_cris = self.args_container.container_cris(
            self._args.fasta_files, self._args.annotation_files,
            self._args.crt_path, self._args.window_size,
            self._args.min_number_repeats, self._args.min_length_repeat,
            self._args.Max_length_repeat, self._args.min_length_spacer,
            self._args.Max_length_spacer, self._paths.crispr_output_folder,
            self._args.ignore_hypothetical_protein)
        cris = Crispr(args_cris)
        cris.run_crispr(args_cris, log)

    def merge(self):
        """Merge all features"""
        print("Merging all features to one gff file")
        merge_folder = os.path.join(self._paths.output_folder,
                                    "merge_all_features")
        self.helper.check_make_folder(merge_folder)
        log = open(os.path.join(merge_folder, "log.txt"), "w")
        other_features = self._args.other_features_files
        self.check_multi_files([[self._args.transcript_file], other_features],
                                ["--transcript_file", "--other_features_files"], log)
        self.check_parameter([self._args.output_prefix], ["--output_prefix"], log)
        run_merge(merge_folder, self._args.transcript_file,
                  self._args.other_features_files,
                  self._args.terminator_tolerance, self._args.tss_tolerance,
                  os.path.join(merge_folder, self._args.output_prefix), log)
        if self._args.source_for_overlapping is not None:
            deal_overlap(merge_folder, self._args.source_for_overlapping)

    def screen(self):
        """generate screenshot"""
        print("Running screenshot generation")
        out_folder = os.path.join(self._args.output_folder, "screenshots")
        if os.path.exists(out_folder):
            print("Error: The {0} already exists!".format(
                  out_folder))
            sys.exit()
        else:
            os.mkdir(out_folder)
        log = open(os.path.join(out_folder, "log.txt"), "w")
        self.check_file([self._args.main_gff, self._args.fasta_file],
                        ["--main_gff", "--fasta_file"], True, log)
        if self._args.side_gffs is not None:
            for gff in (self._args.side_gffs):
                gff = gff.strip()
                if not os.path.isfile(gff):
                    print("Error: The --side_gffs do not exist!")
                    sys.exit()
        if self._args.output_folder is None:
            log.write("No --output_folder can be found.\n")
            print("Error: Please assign --output_folder!")
            sys.exit()
        if (self._args.present != "expand") and (
                self._args.present != "collapse") and (
                self._args.present != "squish"):
            log.write("Please assign \"expand\" or "
                  "\"collapse\" or \"squish\" to --present.\n")
            print("Error: Please assign \"expand\" or "
                  "\"collapse\" or \"squish\" to --present!")
            sys.exit()
        args_sc = self.args_container.container_screen(
            self._args.main_gff, self._args.side_gffs,
            self._args.fasta_file, self._args.height,
            self._args.tex_notex_libs, self._args.frag_libs,
            self._args.present, self._args.output_folder)
        screen = Screen(args_sc, out_folder)
        screen.screenshot(args_sc, log)
