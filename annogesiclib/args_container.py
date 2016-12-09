import os
import sys
import shutil
from glob import glob
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from contextlib import redirect_stdout

class ArgsContainer(object):

    def __init__(self):
        self.multiparser = Multiparser()
        self.helper = Helper()

    def _create_working_wigs(self, out_folder, libs, wig_folder):
        new_libs = []
        if libs is not None:
            self.helper.check_make_folder(wig_folder)
            for lib in libs:
                shutil.copy(lib.split(":")[0], wig_folder)
                wig = lib.split(":")[0].split("/")[-1]
                new_libs.append(":".join([wig, ":".join(lib.split(":")[1:])]))
        else:
            new_libs = None
        return new_libs

    def _check_replicates(self, replicates_tex, replicates_frag,
                          tex_lib, frag_lib):
        '''Check the replicate of frag and tex libs'''
        if (tex_lib is not None) and (replicates_tex is None):
            print("Error: No replicates numbers for "
                  "TEX treated libraries are assigned!")
            sys.exit()
        if (frag_lib is not None) and (replicates_frag is None):
            print("Error: No replicates numbers for "
                  "fragmented libraries are assigned!")
            sys.exit()
        if (replicates_tex is not None) and (replicates_frag is not None):
            replicates = {"tex": replicates_tex,
                          "frag": replicates_frag}
        elif replicates_tex is not None:
            replicates = {"tex": replicates_tex, "frag": -1}
        elif replicates_frag is not None:
            replicates = {"tex": -1, "frag": replicates_frag}
        else:
            print("Error: No replicates number was assigned!")
            sys.exit()
        if replicates["tex"] != -1:
            for rep in replicates["tex"]:
                if ("_" not in rep):
                    print("Error: Please check the input format of replicate_tex! "
                          "It should also contain condition name.")
                    sys.exit()
        if replicates["frag"] != -1:
            for rep in replicates["frag"]:
                if ("_" not in rep):
                    print("Error: Please check the input format of replicate_frag! "
                          "It should also contain condition name.")
                    sys.exit()
        return replicates

    def _check_assign_info(self, infos, file_type, wig_type):
        if file_type == "cond":
            index = 1
        elif file_type == "rep":
            index = 97
        for info, num in sorted(infos.items()):
            if (file_type == "cond") or (file_type == "rep"):
                if file_type == "cond":
                    if info != index:
                        print("Error: The condition number and order "
                              "of --tex_notex_libs should follow 1, 2, 3.")
                        sys.exit()
                elif file_type == "rep":
                    if ord(info) != index:
                        print("Error: The replicate index and order "
                              "of --tex_notex_libs should follow a, b, c.")
                        sys.exit()
                if (wig_type == "tex") and (num % 4 != 0):
                    print("Error: The --tex_notex_libs was assinged incorrectly. "
                          "Please check it again.")
                    sys.exit()
                elif (wig_type == "frag") and (num % 2 != 0):
                    print("Error: The --frag_libs was assinged incorrectly. "
                          "Please check it again.")
                    sys.exit()
                index += 1
            elif file_type == "strand":
                if wig_type == "frag":
                    if (infos["+"] != infos["-"]):
                        print("Error: The --frag_libs was assinged incorrectly. "
                              "Please check it again.")
                        sys.exit()
                if wig_type == "tex":
                    if (num % 2 != 0) and (infos["+"] != infos["-"]):
                        print("Error: The --tex_notex_libs was assinged incorrectly.  "
                              "Please check it again.")
                        sys.exit()

    def _check_tex_frag(self, libs, wig_type):
        conds = {}
        reps = {}
        strands = {}
        for lib in libs:
            datas = lib.split(":")
            if not datas[0].endswith(".wig"):
                print("Error: The input wiggle file should end with .wig!")
                sys.exit()
            if (datas[1] != "notex") and (
                    datas[1] != "tex") and (
                    datas[1] != "frag"):
                print("Error: Please assign \"tex\", \"notex\" or "
                      "\"frag\" to your input libraries.")
                sys.exit()
            if int(datas[2]) not in conds.keys():
                conds[int(datas[2])] = 0
            conds[int(datas[2])] += 1
            if datas[3] not in reps.keys():
                reps[datas[3]] = 0
            reps[datas[3]] += 1
            if (datas[4] != "+") and (datas[4] != "-"):
                print("Error: Strand of libs should be assigned as + or -")
                sys.exit()
            if datas[4] not in strands.keys():
                strands[datas[4]] = 0
            strands[datas[4]] += 1
        self._check_assign_info(conds, "cond", wig_type)
        self._check_assign_info(reps, "rep", wig_type)
        self._check_assign_info(strands, "strand", wig_type)

    def _check_libs(self, tex_notex_libs, frag_libs):
        '''Check the libs of frag and tex'''
        if (tex_notex_libs is None) and (frag_libs is None):
            print("Error: please input proper libraries!!")
        elif (tex_notex_libs is not None) and (frag_libs is not None):
            libs = tex_notex_libs + frag_libs
            self._check_tex_frag(tex_notex_libs, "tex")
            self._check_tex_frag(frag_libs, "frag")
        elif (tex_notex_libs is not None):
            libs = tex_notex_libs
            self._check_tex_frag(tex_notex_libs, "tex")
        elif (frag_libs is not None):
            libs = frag_libs
            self._check_tex_frag(frag_libs, "frag")
        return libs

    def _check_condition_num(self, out_prefix, libs):
        high = 0
        for lib in libs:
            datas = lib.split(":")
            if int(datas[2]) > high:
                high = int(datas[2])
        if len(out_prefix) != high:
            print("Error: The number of --condition_names should be "
                  "the same to the condition of input libraries!")
            sys.exit()

    def _combine_files(self, ref_files, out_folder, filename):
        if ref_files is not None:
            tar_file = os.path.join(out_folder, filename)
            if os.path.exists(tar_file):
                os.remove(tar_file)
            for files in ref_files:
                for file_ in glob(files):
                    self.helper.merge_file(file_, tar_file)
            return tar_file
        else:
            return None

    def _merge_by_strain(self, wig_path, libs):
        strains = []
        merge_folder = os.path.join(wig_path, "merge_tmp")
        self.helper.check_make_folder(merge_folder)
        for wig in os.listdir(wig_path):
            if "_STRAIN_" in wig:
                strain = wig.split("_STRAIN_")[-1].replace(".wig", "")
                if strain not in strains:
                    strains.append(strain)
        for strain in strains:
            change_f = False
            change_r = False
            for wig in os.listdir(wig_path):
                filename = wig.split("_STRAIN_")
                if ("_STRAIN_" in wig) and (
                        filename[-1].replace(
                            ".wig", "") == strain):
                    for lib in libs:
                        if (filename[0] in lib) and (lib[-1] == "+"):
                            self.helper.merge_file(
                                os.path.join(wig_path, wig),
                                os.path.join(merge_folder,
                                             "tmp_forward.wig"))
                            change_f = True
                        elif (filename[0] in lib) and (lib[-1] == "-"):
                            self.helper.merge_file(
                                os.path.join(wig_path, wig),
                                os.path.join(merge_folder,
                                             "tmp_reverse.wig"))
                            change_r = True
            if change_f and change_r:
                change_f = False
                change_r = False
                shutil.move(os.path.join(merge_folder, "tmp_forward.wig"),
                            os.path.join(merge_folder,
                                         strain + "_forward.wig"))
                shutil.move(os.path.join(merge_folder, "tmp_reverse.wig"),
                            os.path.join(merge_folder,
                                         strain + "_reverse.wig"))
            else:
                print("Error: .wig files should be compose of "
                      "forward or reverse files.")
                sys.exit()
        self.helper.remove_all_content(wig_path, ".wig", "file")
        self.helper.move_all_content(merge_folder, wig_path, None)
        shutil.rmtree(merge_folder)

    def _parser_combine_wigs(self, subcommand):
        '''Check the wig folders of frag and tex, then merge them'''
        self.tex_path = None
        self.frag_path = None
        if subcommand == "transcript":
            if self.gffs is not None:
                self.multiparser.parser_gff(self.gffs, None)
                gff_path = self.gffs
        elif subcommand == "terminator":
            self.multiparser.parser_gff(self.gffs, None)
            gff_path = os.path.join(self.gffs, "tmp")
            tmp_file = os.path.join(self.out_folder, "tmp.txt")
            with open(tmp_file, 'w') as fh:
                with redirect_stdout(fh):
                    self.multiparser.parser_gff(gff_path, None)
            os.remove(tmp_file)
        else:
            self.multiparser.parser_gff(self.gffs, None)
            gff_path = self.gffs
        if self.tex_wigs is not None:
            self.tex_path = os.path.join(self.tex_wigs, "tmp")
            self.multiparser.parser_wig(self.tex_wigs)
            if self.gffs is not None:
                self.multiparser.combine_wig(gff_path, self.tex_path,
                                             None, self.libs)
            else:
                self._merge_by_strain(self.tex_path, self.libs)
            self.merge_wigs = self.tex_wigs
            self.wig_path = self.tex_path
        if self.frag_wigs is not None:
            self.frag_path = os.path.join(self.frag_wigs, "tmp")
            self.multiparser.parser_wig(self.frag_wigs)
            if self.gffs is not None:
                self.multiparser.combine_wig(gff_path, self.frag_path,
                                             None, self.libs)
            else:
                self._merge_by_strain(self.frag_path, self.libs)
            self.merge_wigs = self.frag_wigs
            self.wig_path = self.frag_path
        if (self.tex_path is not None) and (
                self.frag_path is not None):
            self = self._merge_wig()
        if (self.tex_path is None) and (
                self.frag_path is None):
            print("Error: There is no proper wig files assigned!!")
            sys.exit()
        return self

    def _merge_wig(self):
        '''Copy the wig files to one folder'''
        self.merge_wigs = os.path.join(self.out_folder, "merge_wigs")
        if (self.tex_wigs is not None) and (
                self.frag_wigs is not None):
            self.helper.check_make_folder(self.merge_wigs)
            self.wig_path = os.path.join(self.merge_wigs, "tmp")
            self.helper.check_make_folder(self.wig_path)
            for wig in os.listdir(self.tex_wigs):
                if os.path.isfile(os.path.join(self.tex_wigs, wig)):
                    shutil.copy(os.path.join(self.tex_wigs, wig),
                                self.merge_wigs)
            for wig in os.listdir(self.frag_wigs):
                if os.path.isfile(os.path.join(self.frag_wigs, wig)):
                    shutil.copy(os.path.join(self.frag_wigs, wig),
                                self.merge_wigs)
            for wig in os.listdir(self.tex_path):
                if os.path.isfile(os.path.join(self.tex_path, wig)):
                    shutil.copy(os.path.join(self.tex_path, wig),
                                self.wig_path)
            for wig in os.listdir(self.frag_path):
                if os.path.isfile(os.path.join(self.frag_path, wig)):
                    self.helper.merge_file(os.path.join(self.frag_path, wig),
                                           os.path.join(self.wig_path, wig))
        elif (self.tex_wigs is not None):
            self.merge_wigs = self.tex_wigs
        elif (self.frag_wigs is not None):
            self.merge_wigs = self.frag_wigs
        return self

    def _deal_multi_inputs(self, inputs, file_type, num, command):
        '''It is for split the input if it is assigned to multiple factors'''
        if inputs is not None:
            datas = inputs.split(",")
            if num is not None:
                if (len(datas) != num):
                    print("Error: the amount of {0} is not correct!!".format(
                        command))
            new_inputs = []
            for data in datas:
                if file_type == "float":
                    new_inputs.append(float(data.strip()))
                elif file_type == "int":
                    new_inputs.append(int(data.strip()))
                else:
                    new_inputs.append(data)
            return new_inputs
        else:
            return inputs

    def _gen_copy_new_folder(self, file_types, out_folder,
                             folder_name, ref_files, flag):
        if ref_files is not None:
            new_ref_folder = os.path.join(out_folder, folder_name)
            self.helper.check_make_folder(new_ref_folder)
            for files in ref_files:
                detect = False
                for file_ in glob(files):
                    for type_ in file_types:
                        if file_.endswith(type_):
                            detect = True
                if not detect:
                    print("Error: The {0} is not end with {1}!".format(
                          flag, " ".join(file_types)))
                    sys.exit()
                shutil.copy(file_, new_ref_folder)
            return new_ref_folder
        else:
            return None

    def container_ratt(self, ratt_path, element, transfer_type,
                       ref_embl, ref_gbk, target_fasta, ref_fasta, ratt_folder,
                       convert_to_gff_rnt_ptt, tar_annotation_folder,
                       compare_pair):
        self.ratt_path = ratt_path
        self.element = element
        self.transfer_type = transfer_type
        self.ref_embls = self._gen_copy_new_folder(
                [".embl"], ratt_folder, "temp_embl",
                ref_embl, ["--ref_embl_files"])
        self.ref_gbk = self._gen_copy_new_folder(
                [".gbk", ".gbff", ".gb"], ratt_folder, "temp_gbk",
                ref_gbk, ["--ref_gbk_files"])
        file_types = [".fa", ".fna", ".fasta"]
        self.tar_fastas = self._gen_copy_new_folder(
                file_types, ratt_folder, "temp_tar", target_fasta,
                ["--ref_fasta_files"])
        self.ref_fastas = self._gen_copy_new_folder(
                file_types, ratt_folder, "temp_ref", ref_fasta,
                ["--target_fasta_files"])
        self.output_path = ratt_folder
        self.convert = convert_to_gff_rnt_ptt
        self.gff_outfolder = tar_annotation_folder
        self.pairs = compare_pair
        return self

    def container_tsspredator(self, TSSpredator_path, compute_program,
                              fasta_files, annotation_files, lib,
                              output_prefix, height, height_reduction, factor,
                              factor_reduction, base_height, enrichment_factor,
                              processing_factor, replicate_match, out_folder,
                              statistics, validate_gene, merge_manual,
                              compare_transcript_assembly, fuzzy, utr_length,
                              cluster, length, re_check_orphan,
                              overlap_feature, reference_gff_files,
                              remove_low_expression):
        self.tsspredator_path = TSSpredator_path
        self.program = compute_program
        self.fastas = self._gen_copy_new_folder(
                [".fa", ".fna", ".fasta"], out_folder, "tmp_fasta", fasta_files,
                ["--fasta_files"])
        self.gffs = self._gen_copy_new_folder(
                [".gff"], out_folder, "tmp_anno", annotation_files,
                ["--annotation_files"])
        self.wig_folder = os.path.join(out_folder, "tmp_wig")
        self.helper.check_make_folder(self.wig_folder)
        self.libs = self._create_working_wigs(out_folder, lib, self.wig_folder)
        self.libs = self._check_libs(self.libs, None)
        self._check_condition_num(output_prefix, self.libs)
        self.output_prefixs = output_prefix
        self.height = height
        self.height_reduction = height_reduction
        self.factor = factor
        self.factor_reduction = factor_reduction
        self.base_height = base_height
        self.enrichment_factor = enrichment_factor
        self.processing_factor = processing_factor
        self.repmatch = replicate_match
        self.out_folder = out_folder
        self.stat = statistics
        self.validate = validate_gene
        self.manual = self._combine_files(merge_manual, out_folder,
                                          "tmp_manual_file")
        self.ta_files = self._gen_copy_new_folder(
                [".gff"], out_folder, "tmp_ta", compare_transcript_assembly,
                ["--compare_transcript_assembly"])
        self.fuzzy = fuzzy
        self.utr_length = utr_length
        self.cluster = cluster
        self.nt_length = length
        self.check_orphan = re_check_orphan
        self.overlap_feature = overlap_feature
        self.references = self._gen_copy_new_folder(
                [".gff"], out_folder, "tmp_reference", reference_gff_files,
                ["--eference_gff_files"])
        self.remove_low_expression = remove_low_expression
        return self

    def container_optimize(self, TSSpredator_path, fasta_file, annotation_file,
                           manual, out_folder, strain_name,
                           max_height, max_height_reduction, max_factor,
                           max_factor_reduction, max_base_height,
                           max_enrichment_factor, max_processing_factor,
                           utr_length, lib, output_prefix, cluster, length,
                           core, program, replicate_match, steps):
        self.tsspredator_path = TSSpredator_path
        self.fastas = self._gen_copy_new_folder(
                [".fa", ".fna", ".fasta"], out_folder, "tmp_fasta",
                [fasta_file], ["--fasta_file"])
        self.gffs = self._gen_copy_new_folder(
                [".gff"], out_folder, "tmp_anno",
                [annotation_file], ["--annotation_file"])
        self.wigs = os.path.join(out_folder, "tmp_wig")
        self.helper.check_make_folder(self.wigs)
        self.libs = self._create_working_wigs(out_folder, lib, self.wigs)
        self.libs = self._check_libs(self.libs, None)
        self.manual = manual
        self.output_folder = out_folder
        self.project_strain = strain_name
        self.height = max_height
        self.height_reduction = max_height_reduction
        self.factor = max_factor
        self.factor_reduction = max_factor_reduction
        self.base_height = max_base_height
        self.enrichment = max_enrichment_factor
        self.processing = max_processing_factor
        self.utr = utr_length
        self._check_condition_num(output_prefix, self.libs)
        self.replicate_name = output_prefix
        self.cluster = cluster
        self.length = length
        self.cores = core
        self.program = program
        self.replicate = replicate_match
        self.steps = steps
        return self

    def _create_wig_folder(self, folder, libs):
        if libs is not None:
            self.helper.check_make_folder(folder)
            return folder
        else:
            return None

    def container_terminator(
            self, TransTermHP_path, expterm_path, RNAfold_path, out_folder,
            fasta_files, annotation_files, transcript_files, srna, statistics,
            decrease, highest_coverage, fuzzy_detect_coverage,
            fuzzy_within_transcript, fuzzy_downstream_transcript,
            fuzzy_within_gene, fuzzy_downstream_gene, transtermhp_folder,
            tex_notex_libs, frag_libs, tex_notex, replicates_tex,
            replicates_frag, table_best, min_loop_length, max_loop_length,
            min_stem_length, max_stem_length, min_AT_tail_length, miss_rate,
            range_u, keep_multi, window, shift):
        self.TransTermHP_path = TransTermHP_path
        self.expterm_path = expterm_path
        self.RNAfold_path = RNAfold_path
        self.out_folder = out_folder
        self.fastas = self._gen_copy_new_folder(
                [".fa", ".fna", ".fasta"], out_folder, "tmp_fasta", fasta_files,
                ["--fasta_files"])
        self.gffs = self._gen_copy_new_folder(
                [".gff"], out_folder, "tmp_anno", annotation_files,
                ["--annotation_files"])
        self.trans = self._gen_copy_new_folder(
                [".gff"], out_folder, "tmp_ta", transcript_files,
                ["--transcript_files"])
        self.srnas = self._gen_copy_new_folder(
                [".gff"], out_folder, "tmp_srna", srna, ["srna_files"])
        self.stat = statistics
        self.helper.check_make_folder(os.path.join(out_folder, "tmp_wig"))
        self.tex_wigs = self._create_wig_folder(
                os.path.join(out_folder, "tmp_wig", "tex_notex"),
                tex_notex_libs)
        self.frag_wigs = self._create_wig_folder(
                os.path.join(out_folder, "tmp_wig", "frag"), frag_libs)
        self.decrease = decrease
        self.cutoff_coverage = highest_coverage
        self.fuzzy = fuzzy_detect_coverage
        self.fuzzy_up_ta = fuzzy_within_transcript
        self.fuzzy_down_ta = fuzzy_downstream_transcript
        self.fuzzy_up_gene = fuzzy_within_gene
        self.fuzzy_down_gene = fuzzy_downstream_gene
        self.hp_folder = transtermhp_folder
        self.tlibs = tex_notex_libs
        self.tlibs = self._create_working_wigs(
                out_folder, tex_notex_libs, self.tex_wigs)
        self.flibs = self._create_working_wigs(
                out_folder, frag_libs, self.frag_wigs)
        self.libs = self._check_libs(self.tlibs, self.flibs)
        self.tex_notex = tex_notex
        self.replicates_tex = replicates_tex
        self.replicates_frag = replicates_frag
        self.replicates = self._check_replicates(
                replicates_tex, replicates_frag, tex_notex_libs, frag_libs)
        self.table_best = table_best
        self.min_loop = min_loop_length
        self.max_loop = max_loop_length
        self.min_stem = min_stem_length
        self.max_stem = max_stem_length
        self.at_tail = min_AT_tail_length
        self.miss_rate = miss_rate
        self.range_u = range_u
        self.keep_multi = keep_multi
        self.window = window
        self.shift = shift
        self = self._parser_combine_wigs("terminator")
        return self

    def container_transcript(self, tex_notex, length, annotation_files, height,
                             width, tolerance, tolerance_coverage,
                             replicates_tex, replicates_frag, out_folder,
                             tss_files, TSS_fuzzy, tex_treated_libs,
                             fragmented_libs, compare_feature_genome,
                             table_best, terminator_files, fuzzy_term, max_dist):
        if (compare_feature_genome is not None) and (annotation_files is None):
            print("Error: --annotation_files needs to be assigned if "
                  "--compare_feature_genome is assigned.")
            sys.exit()
        self.helper.check_make_folder(os.path.join(out_folder, "tmp_wig"))
        self.tex_wigs = self._create_wig_folder(
                os.path.join(out_folder, "tmp_wig", "tex_notex"),
                tex_treated_libs)
        self.frag_wigs = self._create_wig_folder(
                os.path.join(out_folder, "tmp_wig", "frag"), fragmented_libs)
        self.tex = tex_notex
        self.length = length
        self.gffs = self._gen_copy_new_folder(
                [".gff"], out_folder, "tmp_anno", annotation_files,
                ["--annotation_files"])
        self.height = height
        self.width = width
        self.tolerance = tolerance
        self.low_cutoff = tolerance_coverage
        self.replicates_tex = replicates_tex
        self.replicates_frag = replicates_frag
        self.replicates = self._check_replicates(
                replicates_tex, replicates_frag,
                tex_treated_libs, fragmented_libs)
        self.out_folder = out_folder
        self.compare_tss = self._gen_copy_new_folder(
                [".gff"], out_folder, "tmp_tss", tss_files, ["--tss_files"])
        self.fuzzy = TSS_fuzzy
        self.tlibs = self._create_working_wigs(
                out_folder, tex_treated_libs, self.tex_wigs)
        self.flibs = self._create_working_wigs(
                out_folder, fragmented_libs, self.frag_wigs)
        self.libs = self._check_libs(self.tlibs, self.flibs)
        self.c_feature = compare_feature_genome
        self.table_best = table_best
        self.terms = self._gen_copy_new_folder(
                [".gff"], out_folder, "tmp_term", terminator_files,
                ["--terminator_files"])
        self.fuzzy_term = fuzzy_term
        self.max_dist = max_dist
        self = self._parser_combine_wigs("transcript")
        return self

    def container_utr(self, tss_files, annotation_files,
                      transcript_assembly_files, terminator_files,
                      terminator_fuzzy, utr_folder, tss_source, base_5utr,
                      length, base_3utr, fuzzy_3utr, fuzzy_5utr):
        self.tsss = self._gen_copy_new_folder(
                [".gff"], utr_folder, "tmp_tss", tss_files, ["--tss_files"])
        self.gffs = self._gen_copy_new_folder(
                [".gff"], utr_folder, "tmp_anno", annotation_files,
                ["--annotation_files"])
        self.trans = self._gen_copy_new_folder(
                [".gff"], utr_folder, "tmp_ta", transcript_assembly_files,
                ["--transcript_files"])
        self.terms = self._gen_copy_new_folder(
                [".gff"], utr_folder, "tmp_term", terminator_files,
                ["--terminator_files"])
        self.fuzzy = terminator_fuzzy
        self.out_folder = utr_folder
        self.source = tss_source
        self.base_5utr = base_5utr
        self.base_3utr = base_3utr
        self.length = length
        self.fuzzy_3utr = fuzzy_3utr
        self.fuzzy_5utr = fuzzy_5utr
        return self

    def container_srna(self, rnafold, relplot_pl, mountain_pl, blastn, blastx,
                       blastdb, ps2pdf14_path, srna_folder, UTR_derived_sRNA,
                       annotation_files, TSS_files, transcript_files,
                       TSS_intergenic_fuzzy, TSS_5UTR_fuzzy, TSS_3UTR_fuzzy,
                       TSS_interCDS_fuzzy, import_info, processing_site_files,
                       fasta_files, mountain_plot, nr_format, srna_format,
                       sRNA_database_path, nr_database_path, cutoff_energy,
                       para_blast, run_intergenic_TEX_coverage,
                       run_intergenic_noTEX_coverage,
                       run_intergenic_fragmented_coverage, break_tran,
                       run_antisense_TEX_coverage,
                       run_antisense_noTEX_coverage,
                       run_antisense_fragmented_coverage, run_utr_TEX_coverage,
                       run_utr_noTEX_coverage, run_utr_fragmented_coverage,
                       max_length, min_length, tex_notex_libs, frag_libs,
                       replicates_tex, replicates_frag, tex_notex, blast_e_nr,
                       blast_e_srna, detect_sRNA_in_CDS, table_best,
                       decrease_intergenic, decrease_utr, fuzzy_intergenic,
                       fuzzy_utr, cutoff_nr_hit, sORF, overlap_percent_CDS,
                       terminator_files, terminator_fuzzy_in_sRNA,
                       terminator_fuzzy_out_sRNA, ignore_hypothetical_protein,
                       TSS_source, min_utr_coverage, promoter_tables,
                       ranking_promoter, promoter_name):
        self.rnafold = rnafold
        self.para_blast = para_blast
        self.relplot_pl = relplot_pl
        self.mountain_pl = mountain_pl
        self.blastx = blastx
        self.blastn = blastn
        self.blastdb = blastdb
        self.ps2pdf14_path = ps2pdf14_path
        self.out_folder = srna_folder
        self.utr_srna = UTR_derived_sRNA
        self.gffs = self._gen_copy_new_folder(
                [".gff"], srna_folder, "temp_anno", annotation_files,
                ["--annotation_files"])
        self.tss_folder = self._gen_copy_new_folder(
                [".gff"], srna_folder, "temp_tss", TSS_files, ["--tss_files"])
        self.trans = self._gen_copy_new_folder(
                [".gff"], srna_folder, "temp_ta", transcript_files,
                ["--transcript_files"])
        self.fuzzy_inter_tss = TSS_intergenic_fuzzy
        self.fuzzy_5utr_tss = TSS_5UTR_fuzzy
        self.fuzzy_3utr_tss = TSS_3UTR_fuzzy
        self.fuzzy_intercds_tss = TSS_interCDS_fuzzy
        self.fuzzy_tsss = {"5utr": self.fuzzy_5utr_tss,
                           "3utr": self.fuzzy_3utr_tss,
                           "interCDS": self.fuzzy_intercds_tss,
                           "inter": self.fuzzy_inter_tss}
        self.import_info = import_info
        self.helper.check_make_folder(os.path.join(srna_folder, "temp_wig"))
        self.tex_wigs = self._create_wig_folder(
                os.path.join(srna_folder, "temp_wig", "tex_notex"),
                tex_notex_libs)
        self.frag_wigs = self._create_wig_folder(
                os.path.join(srna_folder, "temp_wig", "frag"), frag_libs)
        self.pro_folder = self._gen_copy_new_folder(
                [".gff"], srna_folder, "temp_pro", processing_site_files,
                ["--processing_site_files"])
        self.fastas = self._gen_copy_new_folder(
                [".fa", ".fna", ".fasta"], srna_folder,
                "temp_fasta", fasta_files, ["--fasta_files"])
        self.mountain = mountain_plot
        self.nr_format = nr_format
        self.srna_format = srna_format
        self.srna_database = sRNA_database_path
        self.nr_database = nr_database_path
        self.energy = cutoff_energy
        self.coverage_tex = self._deal_multi_inputs(
                run_intergenic_TEX_coverage, "float", 5,
                "--run_intergenic_TEX_coverage")
        self.coverage_notex = self._deal_multi_inputs(
                run_intergenic_noTEX_coverage, "float", 5,
                "--run_intergenic_noTEX_coverage")
        self.coverage_frag = self._deal_multi_inputs(
                run_intergenic_fragmented_coverage, "float", 5,
                "--run_intergenic_fragmented_coverage")
        self.anti_cover_tex = self._deal_multi_inputs(
                run_antisense_TEX_coverage, "float", 5,
                "--run_antisense_TEX_coverage")
        self.anti_cover_notex = self._deal_multi_inputs(
                run_antisense_noTEX_coverage, "float", 5,
                "--run_antisense_noTEX_coverage")
        self.anti_cover_frag = self._deal_multi_inputs(
                run_antisense_fragmented_coverage, "float", 5,
                "--run_antisense_fragmented_coverage")
        self.break_tran = self._deal_multi_inputs(
                break_tran, "float", 3,
                "--run_break_transcript")
        self.utr_tex_cover = self._deal_multi_inputs(
                run_utr_TEX_coverage, "str", 3, "--run_utr_TEX_coverage")
        self.utr_notex_cover = self._deal_multi_inputs(
                run_utr_noTEX_coverage, "str", 3, "--run_utr_TEX_coverage")
        self.utr_frag_cover = self._deal_multi_inputs(
                run_utr_fragmented_coverage, "str", 3,
                "--run_utr_fragmented_coverage")
        self.max_len = max_length
        self.min_len = min_length
        self.tlibs = self._create_working_wigs(
                srna_folder, tex_notex_libs, self.tex_wigs)
        self.flibs = self._create_working_wigs(
                srna_folder, frag_libs, self.frag_wigs)
        self.libs = self._check_libs(self.tlibs, self.flibs)
        self.replicates_tex = replicates_tex
        self.replicates_frag = replicates_frag
        self.replicates = self._check_replicates(
                replicates_tex, replicates_frag, tex_notex_libs, frag_libs)
        self.tex_notex = tex_notex
        self.e_nr = blast_e_nr
        self.e_srna = blast_e_srna
        self.in_cds = detect_sRNA_in_CDS
        self.table_best = table_best
        self.decrease_inter = decrease_intergenic
        self.decrease_utr = decrease_utr
        self.fuzzy_inter = fuzzy_intergenic
        self.fuzzy_utr = fuzzy_utr
        self.nr_hits_num = cutoff_nr_hit
        self.sorf_file = self._gen_copy_new_folder(
                [".gff"], srna_folder, "temp_sorf", sORF, ["--sorf_files"])
        self.cutoff_overlap = overlap_percent_CDS
        self.terms = self._gen_copy_new_folder(
                [".gff"], srna_folder, "temp_term", terminator_files,
                ["--terminator_files"])
        self.fuzzy_b = terminator_fuzzy_in_sRNA
        self.fuzzy_a = terminator_fuzzy_out_sRNA
        self.hypo = ignore_hypothetical_protein
        self.tss_source = TSS_source
        self.min_utr = min_utr_coverage
        self.promoter_table = self._combine_files(
                promoter_tables, srna_folder, "tmp_promoter_table")
        if ranking_promoter < 1:
            print("Error: --ranking_time_promoter must larger than 1...")
            sys.exit()
        self.rank_promoter = ranking_promoter
        self.promoter_name = promoter_name
        self = self._parser_combine_wigs("srna")
        return self

    def container_intersrna(self, file_type, files, args_srna, prefix,
                            gff_file, tran_file, tss_file, pro_file, fuzzy):
        '''Especially for intergenic and antisense sRNA'''
        args_srna.file_type = file_type
        args_srna.gff_file = gff_file
        args_srna.tran_file = tran_file
        args_srna.tss_file = tss_file
        args_srna.pro_file = pro_file
        args_srna.fuzzy = fuzzy
        args_srna.prefix = prefix
        if file_type == "frag":
            args_srna.wig_f_file = os.path.join(
                    args_srna.frag_path, "_".join([prefix, "forward.wig"]))
            args_srna.wig_r_file = os.path.join(
                    args_srna.frag_path, "_".join([prefix, "reverse.wig"]))
            args_srna.wig_folder = args_srna.frag_wigs
            args_srna.input_libs = args_srna.flibs
            args_srna.output_file = files["frag_gff"]
            args_srna.output_table = files["frag_csv"]
            args_srna.cutoffs = args_srna.coverage_frag
            args_srna.tss_source = True
            args_srna.cut_notex = None
            args_srna.anti_notex_cutoff = None
        else:
            args_srna.wig_f_file = os.path.join(
                    args_srna.tex_path, "_".join([prefix, "forward.wig"]))
            args_srna.wig_r_file = os.path.join(
                    args_srna.tex_path, "_".join([prefix, "reverse.wig"]))
            args_srna.wig_folder = args_srna.tex_wigs
            args_srna.input_libs = args_srna.tlibs
            args_srna.output_file = files["tex_gff"]
            args_srna.output_table = files["tex_csv"]
            args_srna.cutoffs = args_srna.coverage_tex
            args_srna.tss_source = args_srna.tss_source
            args_srna.cut_notex = args_srna.coverage_notex
            args_srna.anti_notex_cutoff = args_srna.anti_cover_notex
        return args_srna

    def container_utrsrna(self, gff, tran, tss, files, pro, fasta, file_type,
                          prefix, args_srna):
        '''Especially for UTR-derived sRNA'''
        args_srna.file_type = file_type
        args_srna.gff_file = gff
        args_srna.ta_file = tran
        args_srna.tss_file = tss
        args_srna.pro_file = pro
        args_srna.prefix = prefix
        args_srna.seq_file = fasta
        if file_type == "frag":
            args_srna.wig_f_file = os.path.join(
                    args_srna.frag_path, "_".join([prefix, "forward.wig"]))
            args_srna.wig_r_file = os.path.join(
                    args_srna.frag_path, "_".join([prefix, "reverse.wig"]))
            args_srna.wig_folder = args_srna.frag_wigs
            args_srna.input_libs = args_srna.flibs
            args_srna.output_file = files["frag_gff"]
            args_srna.output_table = files["frag_csv"]
            args_srna.utr_coverages = args_srna.utr_frag_cover
            args_srna.notex = None
        else:
            args_srna.wig_f_file = os.path.join(
                    args_srna.tex_path, "_".join([prefix, "forward.wig"]))
            args_srna.wig_r_file = os.path.join(
                    args_srna.tex_path, "_".join([prefix, "reverse.wig"]))
            args_srna.wig_folder = args_srna.tex_wigs
            args_srna.input_libs = args_srna.tlibs
            args_srna.output_file = files["tex_gff"]
            args_srna.output_table = files["tex_csv"]
            args_srna.utr_coverages = args_srna.utr_tex_cover
            args_srna.notex = args_srna.utr_notex_cover
        args_srna.coverages = {"5utr": args_srna.utr_coverages[0],
                               "3utr": args_srna.utr_coverages[1],
                               "interCDS": args_srna.utr_coverages[2]}
        if args_srna.notex is not None:
            args_srna.cover_notex = {"5utr": args_srna.notex[0],
                                     "3utr": args_srna.notex[1],
                                     "interCDS": args_srna.notex[2]}
        else:
            args_srna.cover_notex = None
        return args_srna

    def extend_inter_container(self, args_srna, tsss, pros,
                               nums, output, out_table, texs, detects,
                               cutoff_coverage, notex):
        '''Especially for intergenic and antisense sRNA'''
        args_srna.tsss = tsss
        args_srna.pros = pros
        args_srna.nums = nums
        args_srna.output = output
        args_srna.out_table = out_table
        args_srna.texs = texs
        args_srna.detects = detects
        args_srna.cutoff_coverage = cutoff_coverage
        args_srna.notex = notex
        return args_srna

    def extend_utr_container(self, args_srna, cdss, tsss, pros,
                             out, out_t, texs):
        '''Especially for UTR-derived sRNA'''
        args_srna.cdss = cdss
        args_srna.tsss = tsss
        args_srna.pros = pros
        args_srna.out = out
        args_srna.out_t = out_t
        args_srna.texs = texs
        args_srna.utrs = []
        args_srna.srnas = []
        return args_srna

    def container_sorf(self, sorf_folder, UTR_derived_sORF, transcript_files,
                       annotation_files, TSS_files, utr_length, min_length,
                       max_length, cutoff_intergenic_coverage,
                       cutoff_antisense_coverage, cutoff_5utr_coverage,
                       cutoff_3utr_coverage, cutoff_interCDS_coverage,
                       fasta_files, tex_notex_libs, frag_libs, tex_notex,
                       replicates_tex, replicates_frag, table_best,
                       sRNA_files, start_codon, stop_codon, cutoff_background,
                       fuzzy_rbs, rbs_not_after_TSS, print_all_combination,
                       best_no_sRNA, best_no_TSS, ignore_hypothetical_protein,
                       min_rbs_distance, max_rbs_distance):
        self.out_folder = sorf_folder
        self.utr_detect = UTR_derived_sORF
        self.trans = self._gen_copy_new_folder(
                [".gff"], sorf_folder, "temp_ta", transcript_files,
                ["--transcript_files"])
        self.gffs = self._gen_copy_new_folder(
                [".gff"], sorf_folder, "temp_anno", annotation_files,
                ["--annotation_files"])
        self.tsss = self._gen_copy_new_folder(
                [".gff"], sorf_folder, "temp_tss", TSS_files, ["--tss_files"])
        self.utr_length = utr_length
        self.min_len = min_length
        self.max_len = max_length
        self.helper.check_make_folder(os.path.join(sorf_folder, "temp_wig"))
        self.tex_wigs = self._create_wig_folder(
                os.path.join(sorf_folder, "temp_wig", "tex_notex"),
                tex_notex_libs)
        self.frag_wigs = self._create_wig_folder(
                os.path.join(sorf_folder, "temp_wig", "frag"), frag_libs)
        self.cutoff_inter = cutoff_intergenic_coverage
        self.cutoff_anti = cutoff_antisense_coverage
        self.cutoff_5utr = cutoff_5utr_coverage
        self.cutoff_3utr = cutoff_3utr_coverage
        self.cutoff_intercds = cutoff_interCDS_coverage
        self.fastas = self._gen_copy_new_folder(
                [".fa", ".fna", ".fasta"], sorf_folder,
                "temp_fasta", fasta_files, ["--fasta_files"])
        self.tlibs = tex_notex_libs
        self.flibs = frag_libs
        self.tlibs = self._create_working_wigs(
                sorf_folder, tex_notex_libs, self.tex_wigs)
        self.flibs = self._create_working_wigs(
                sorf_folder, frag_libs, self.frag_wigs)
        self.libs = self._check_libs(self.tlibs, self.flibs)
        self.tex_notex = tex_notex
        self.replicates_tex = replicates_tex
        self.replicates_frag = replicates_frag
        self.replicates = self._check_replicates(
                replicates_tex, replicates_frag, tex_notex_libs, frag_libs)
        self.table_best = table_best
        self.srnas = self._gen_copy_new_folder(
                [".gff"], sorf_folder, "temp_srna", sRNA_files,
                ["--srna_files"])
        self.start_codon = start_codon
        self.stop_codon = stop_codon
        self.background = cutoff_background
        self.fuzzy_rbs = fuzzy_rbs
        self.noafter_tss = rbs_not_after_TSS
        self.print_all = print_all_combination
        self.no_srna = best_no_sRNA
        self.no_tss = best_no_TSS
        self.hypo = ignore_hypothetical_protein
        self.min_rbs = min_rbs_distance
        self.max_rbs = max_rbs_distance
        self = self._parser_combine_wigs("sorf")
        return self

    def container_srna_target(
            self, rnaplfold_path, rnaplex_path, rnaup_path, annotation_files,
            fasta_files, sRNA_files, query_sRNA, program, interaction_length,
            window_size_target, span_target, window_size_srna, span_srna,
            unstructured_region_RNAplex_target,
            unstructured_region_RNAplex_srna, unstructured_region_RNAup,
            energy_threshold, duplex_distance, top, starget_output_folder,
            process_rnaplex, process_rnaup, continue_rnaup,
            potential_target_start, potential_target_end, target_feature):
        self.rnaplfold_path = rnaplfold_path
        self.rnaplex_path = rnaplex_path
        self.rnaup_path = rnaup_path
        self.gffs = self._gen_copy_new_folder(
                [".gff"], starget_output_folder, "tmp_anno",
                annotation_files, ["--annotation_files"])
        self.fastas = self._gen_copy_new_folder(
                [".fa", ".fna", ".fasta"], starget_output_folder,
                "tmp_fasta", fasta_files, ["--fasta_files"])
        self.srnas = self._gen_copy_new_folder(
                [".gff"], starget_output_folder, "tmp_srna", sRNA_files,
                ["--srna_files"])
        self.query = query_sRNA
        self.program = program
        self.inter_length = interaction_length
        self.win_size_t = window_size_target
        self.span_t = span_target
        self.win_size_s = window_size_srna
        self.span_s = span_srna
        self.unstr_region_rnaplex_t = unstructured_region_RNAplex_target
        self.unstr_region_rnaplex_s = unstructured_region_RNAplex_srna
        self.unstr_region_rnaup = unstructured_region_RNAup
        self.energy = energy_threshold
        self.duplex_dist = duplex_distance
        self.top = top
        self.out_folder = starget_output_folder
        self.core_plex = process_rnaplex
        self.core_up = process_rnaup
        self.continue_rnaup = continue_rnaup
        self.tar_start = potential_target_start
        self.tar_end = potential_target_end
        self.features = target_feature
        return self

    def container_goterm(self, annotation_files, goterm_output_folder,
                         UniProt_id, go_obo, goslim_obo, transcript_files):
        self.gffs = self._gen_copy_new_folder(
                [".gff"], goterm_output_folder, "tmp_anno", annotation_files,
                ["--annotation_files"])
        self.out_folder = goterm_output_folder
        self.uniprot = UniProt_id
        self.go = go_obo
        self.goslim = goslim_obo
        self.trans = self._gen_copy_new_folder(
                [".gff"], goterm_output_folder, "tmp_ta", transcript_files,
                ["--transcript_files"])
        return self

    def container_sublocal(self, Psortb_path, annotation_files, fasta_files,
                           bacteria_type, difference_multi, merge_to_gff,
                           sublocal_output_folder, transcript_files):
        self.psortb_path = Psortb_path
        self.gffs = self._gen_copy_new_folder(
                [".gff"], sublocal_output_folder, "tmp_anno",
                annotation_files, ["--annotation_files"])
        self.fastas = self._gen_copy_new_folder(
                [".fa", ".fna", ".fasta"], sublocal_output_folder,
                "tmp_fa", fasta_files, ["--fasta_files"])
        self.gram = bacteria_type
        self.fuzzy = difference_multi
        self.merge = merge_to_gff
        self.out_folder = sublocal_output_folder
        self.trans = self._gen_copy_new_folder(
                [".gff"], sublocal_output_folder, "tmp_ta",
                transcript_files, ["--transcript_files"])
        return self

    def container_ppi(self, annotation_files, proteinID_strains,
                      without_strain_pubmed, species_STRING, score,
                      ppi_output_folder, node_size, query):
        self.ptts = self._gen_copy_new_folder(
                [".gff"], ppi_output_folder, "temp_anno",
                annotation_files, ["--annotation_files"])
        self.strains = proteinID_strains
        self.no_specific = without_strain_pubmed
        self.species = species_STRING
        self.score = score
        self.out_folder = ppi_output_folder
        self.size = node_size
        self.querys = query
        return self

    def container_promoter(self, MEME_path, GLAM2_path, out_folder, tex_libs,
                           TSS_files, fasta_files, num_motif, nt_before_TSS,
                           motif_width, TSS_source, annotation_files, end_run,
                           combine_all, e_value, para, program):
        self.meme_path = MEME_path
        self.glam2_path = GLAM2_path
        self.program = program
        self.end_run = end_run
        if (program.lower() != "both") and (
                program.lower() != "meme") and (
                program.lower() != "glam2"):
            print("Error: Please assign meme or glam2 or both to --program.")
            sys.exit()
        self.output_folder = out_folder
        self.tsss = self._gen_copy_new_folder(
                [".gff"], out_folder, "tmp_tss", TSS_files, ["--tss_files"])
        self.fastas = self._gen_copy_new_folder(
                [".fa", ".fna", ".fasta"], out_folder, "tmp_fasta",
                fasta_files, ["--fasta_files"])
        self.num_motif = num_motif
        self.nt_before = nt_before_TSS
        self.widths = motif_width
        self.source = TSS_source
        self.tex_wigs = None
        self.frag_wigs = None
        self.gffs = self._gen_copy_new_folder(
                [".gff"], out_folder, "tmp_anno", annotation_files,
                ["--annotation_files"])
        self.combine = combine_all
        self.e_value = e_value
        self.para = para
        if tex_libs is not None:
            self.helper.check_make_folder(os.path.join(out_folder, "tmp_wig"))
            self.tex_wigs = self._create_wig_folder(
                    os.path.join(out_folder, "tmp_wig", "tex_notex"), tex_libs)
            self.input_libs = self._create_working_wigs(
                     out_folder, tex_libs, self.tex_wigs)
            self.libs = self.input_libs
            self = self._parser_combine_wigs("promoter")
        return self

    def container_operon(self, TSS_files, annotation_files,
                         transcript_files, UTR5_files, UTR3_files,
                         term_files, TSS_fuzzy, term_fuzzy, min_length,
                         statistics, operon_output_folder, combine_gff,
                         operon_statistics_folder):
        self.tsss = self._gen_copy_new_folder(
                [".gff"], operon_output_folder, "tmp_tss",
                TSS_files, ["--tss_files"])
        self.gffs = self._gen_copy_new_folder(
                [".gff"], operon_output_folder, "tmp_anno",
                annotation_files, ["--annotation_files"])
        self.trans = self._gen_copy_new_folder(
                [".gff"], operon_output_folder, "tmp_ta",
                transcript_files, ["--transcript_files"])
        self.utr5s = self._gen_copy_new_folder(
                [".gff"], operon_output_folder, "tmp_utr5",
                UTR5_files, ["--utr5_files"])
        self.utr3s = self._gen_copy_new_folder(
                [".gff"], operon_output_folder, "tmp_utr3",
                UTR3_files, ["--utr3_files"])
        self.terms = self._gen_copy_new_folder(
                [".gff"], operon_output_folder, "tmp_term",
                term_files, ["--term_files"])
        self.tss_fuzzy = TSS_fuzzy
        self.term_fuzzy = term_fuzzy
        self.length = min_length
        self.statistics = statistics
        self.output_folder = operon_output_folder
        self.combine = combine_gff
        self.stat_folder = operon_statistics_folder
        return self

    def container_snp(self, samtools_path, bcftools_path, bam_type, min_sample,
                      program, fasta_files, bam_files,
                      quality, read_depth_range, snp_output_folder,
                      indel_fraction, chrom, rg, caller, filters, DP4_cutoff):
        self.samtools_path = samtools_path
        self.bcftools_path = bcftools_path
        self.types = bam_type
        self.program = program 
        self.fastas = self._gen_copy_new_folder(
                [".fa", ".fna", ".fasta"], snp_output_folder,
                "tmp_fa", fasta_files, ["--fasta_files"])
        self.bams = bam_files
        self.quality = quality
        self.depth_s = read_depth_range.split(",")[0]
        self.depth_b = read_depth_range.split(",")[-1]
        self.out_folder = snp_output_folder
        self.idv = indel_fraction.split(",")[0]
        self.imf = indel_fraction.split(",")[-1]
        if chrom == "haploid":
            chrom = "1"
        elif chrom == "diploid":
            chrom = "2"
        self.chrom = chrom
        self.rg = rg
        self.caller = caller
        self.filters = filters
        self.dp4_sum = DP4_cutoff.split(",")[0]
        self.dp4_frac = DP4_cutoff.split(",")[-1]
        self.min_sample = min_sample
        return self

    def container_circrna(self, align, process, fasta_files, annotation_files,
                          bam_files, read_files,
                          circrna_stat_folder, support_reads, segemehl_path,
                          testrealign, samtools_path, start_ratio,
                          end_ratio, ignore_hypothetical_protein, out_folder):
        self.align = align
        self.cores = process
        self.fastas = self._gen_copy_new_folder(
                [".fa", ".fna", ".fasta"], out_folder, "tmp_fa", fasta_files,
                ["--fasta_files"])
        self.gffs = self._gen_copy_new_folder(
                [".gff"], out_folder, "tmp_anno", annotation_files,
                ["--annotation_files"])
        self.bams = bam_files
        self.read_files = read_files
        self.stat_folder = circrna_stat_folder
        self.support = support_reads
        self.segemehl_path = segemehl_path
        self.testrealign_path = testrealign
        self.samtools_path = samtools_path
        self.start_ratio = start_ratio
        self.end_ratio = end_ratio
        self.hypo = ignore_hypothetical_protein
        self.output_folder = out_folder
        return self

    def container_ribos(self, program, thermo_ID, cmscan_path, cmpress_path,
                        riboswitch_ID, annotation_files, fasta_files,
                        tss_files, transcript_files, Rfam, ribos_output_folder,
                        thermo_output_folder, e_value, output_all,
                        database_folder, fuzzy, start_codon, min_dist_rbs,
                        max_dist_rbs, fuzzy_rbs, UTR_length):
        self.program = program
        if (program.lower() == "riboswitch") or (
                program.lower() == "both"):
            output = ribos_output_folder
        elif (program.lower() == "thermometer"):
            output = thermo_output_folder
        self.thermo_id = thermo_ID
        self.cmscan_path = cmscan_path
        self.cmpress_path = cmpress_path
        self.ribos_id = riboswitch_ID
        self.gffs = self._gen_copy_new_folder(
                [".gff"], output, "temp_anno", annotation_files,
                ["--annotation_files"])
        self.fastas = self._gen_copy_new_folder(
                [".fa", ".fna", ".fasta"], output, "temp_fa",
                fasta_files, ["--fasta_files"])
        self.tsss = self._gen_copy_new_folder(
                [".gff"], output, "temp_tss", tss_files, ["--tss_files"])
        self.trans = self._gen_copy_new_folder(
                [".gff"], output, "temp_ta", transcript_files,
                ["--transcript_files"])
        self.rfam = Rfam
        self.ribos_out_folder = ribos_output_folder
        self.thermo_out_folder = thermo_output_folder
        self.e_value = e_value
        self.output_all = output_all
        self.database = database_folder
        self.fuzzy = fuzzy
        self.start_codons = start_codon
        self.start_rbs = min_dist_rbs
        self.end_rbs = max_dist_rbs
        self.fuzzy_rbs = fuzzy_rbs
        self.utr = UTR_length
        return self

    def container_cris(self, fasta_files, annotation_files, CRT_path,
                       window_size, min_number_repeat, min_length_repeat,
                       Max_length_repeat, min_length_spacer, Max_length_spacer,
                       cris_out_folder, ignore_hypo):
        self.gffs = self._gen_copy_new_folder(
                [".gff"], cris_out_folder, "tmp_anno", annotation_files,
                ["--annotation_files"])
        self.fastas = self._gen_copy_new_folder(
                [".fa", ".fna", ".fasta"], cris_out_folder, "tmp_fa",
                fasta_files, ["--fasta_files"])
        self.crt_path = CRT_path
        self.win_size = window_size
        self.out_folder = cris_out_folder
        self.min_num_r = min_number_repeat
        self.min_len_r = min_length_repeat
        self.max_len_r = Max_length_repeat
        self.min_len_s = min_length_spacer
        self.max_len_s = Max_length_spacer
        self.ignore_hypo = ignore_hypo
        return self

    def container_screen(self, main_gff, side_gffs, fasta, height, tex_libs,
                         frag_libs, present, output_folder):
        self.main_gff = main_gff
        self.side_gffs = side_gffs
        self.fasta = fasta
        self.height = height
        self.tlibs = tex_libs
        self.flibs = frag_libs
        self.present = present
        self.output_folder = output_folder
        return self
