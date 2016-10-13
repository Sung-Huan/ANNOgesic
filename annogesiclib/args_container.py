import os
import sys
import shutil
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper


class ArgsContainer(object):

    def __init__(self):
        self.multiparser = Multiparser()
        self.helper = Helper()

    def _check_replicates(self, replicates_tex, replicates_frag):
        '''Check the replicate of frag and tex libs'''
        if (replicates_tex is not None) and (replicates_frag is not None):
            replicates = {"tex": replicates_tex,
                          "frag": replicates_frag}
        elif replicates_tex is not None:
            replicates = {"tex": replicates_tex, "frag": -1}
        elif replicates_frag is not None:
            replicates = {"tex": -1, "frag": replicates_frag}
        else:
            print("Error:No replicates number assign!!!")
            sys.exit()
        return replicates

    def _check_libs(self, tex_notex_libs, frag_libs):
        '''Check the libs of frag and tex'''
        if (tex_notex_libs is None) and (frag_libs is None):
            print("Error: please input proper libraries!!")
        elif (tex_notex_libs is not None) and (frag_libs is not None):
            libs = tex_notex_libs + frag_libs
        elif (tex_notex_libs is not None):
            libs = tex_notex_libs
        elif (frag_libs is not None):
            libs = frag_libs
        return libs

    def _parser_combine_wigs(self, subcommand):
        '''Check the wig folders of frag and tex, then merge them'''
        self.tex_path = None
        self.frag_path = None
        self.multiparser.parser_gff(self.gffs, None)
        if subcommand == "terminator":
            gff_path = os.path.join(self.gffs, "tmp")
            self.multiparser.parser_gff(gff_path, None)
        else:
            gff_path = self.gffs
        if self.tex_wigs is not None:
            self.tex_path = os.path.join(self.tex_wigs, "tmp")
            self.multiparser.parser_wig(self.tex_wigs)
            self.multiparser.combine_wig(gff_path, self.tex_path,
                                         None, self.libs)
            self.merge_wigs = self.tex_wigs
            self.wig_path = self.tex_path
        if self.frag_wigs is not None:
            self.frag_path = os.path.join(self.frag_wigs, "tmp")
            self.multiparser.parser_wig(self.frag_wigs)
            self.multiparser.combine_wig(gff_path, self.frag_path,
                                         None, self.libs)
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

    def container_ratt(self, ratt_path, element, transfer_type,
                       ref_embl, ref_gbk, target_fasta, ref_fasta, ratt_folder,
                       convert_to_gff_rnt_ptt, tar_annotation_folder,
                       compare_pair):
        self.ratt_path = ratt_path
        self.element = element
        self.transfer_type = transfer_type
        self.ref_embls = ref_embl
        self.ref_gbk = ref_gbk
        self.tar_fastas = target_fasta
        self.ref_fastas = ref_fasta
        self.output_path = ratt_folder
        self.convert = convert_to_gff_rnt_ptt
        self.gff_outfolder = tar_annotation_folder
        self.pairs = self._deal_multi_inputs(compare_pair, "str", None, None)
        return self

    def container_tsspredator(self, TSSpredator_path, compute_program,
                              fasta_folder, annotation_folder, wig_folder, lib,
                              output_prefix, height, height_reduction, factor,
                              factor_reduction, base_height, enrichment_factor,
                              processing_factor, replicate_match, out_folder,
                              statistics, validate_gene, merge_manual,
                              compare_transcript_assembly, fuzzy, utr_length,
                              cluster, length, re_check_orphan,
                              overlap_feature, reference_gff_folder,
                              remove_low_expression):
        self.tsspredator_path = TSSpredator_path
        self.program = compute_program
        self.fastas = fasta_folder
        self.gffs = annotation_folder
        self.wig_folder = wig_folder
        self.libs = self._deal_multi_inputs(lib, "str", None, None)
        self.output_prefixs = self._deal_multi_inputs(output_prefix, "str",
                                                      None, None)
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
        self.manual = merge_manual
        self.ta_files = compare_transcript_assembly
        self.fuzzy = fuzzy
        self.utr_length = utr_length
        self.cluster = cluster
        self.nt_length = length
        self.check_orphan = re_check_orphan
        self.overlap_feature = overlap_feature
        self.references = reference_gff_folder
        self.remove_low_expression = remove_low_expression
        return self

    def container_optimize(self, TSSpredator_path, fasta_file, annotation_file,
                           wig_folder, manual, out_folder, strain_name,
                           max_height, max_height_reduction, max_factor,
                           max_factor_reduction, max_base_height,
                           max_enrichment_factor, max_processing_factor,
                           utr_length, lib, output_prefix, cluster, length,
                           core, program, replicate_match, steps):
        self.tsspredator_path = TSSpredator_path
        self.fastas = fasta_file
        self.gffs = annotation_file
        self.wigs = wig_folder
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
        self.libs = self._deal_multi_inputs(lib, "str", None, None)
        self.replicate_name = self._deal_multi_inputs(output_prefix, "str",
                                                      None, None)
        self.cluster = cluster
        self.length = length
        self.cores = core
        self.program = program
        self.replicate = replicate_match
        self.steps = steps
        return self

    def container_terminator(
            self, TransTermHP_path, expterm_path, RNAfold_path, out_folder,
            fasta_folder, annotation_folder, transcript_folder, srna,
            statistics, tex_wig_folder, frag_wig_folder, decrease,
            highest_coverage, fuzzy_detect_coverage, fuzzy_within_transcript,
            fuzzy_downstream_transcript, fuzzy_within_gene,
            fuzzy_downstream_gene, transtermhp_folder, tex_notex_libs,
            frag_libs, tex_notex, replicates_tex, replicates_frag, table_best,
            min_loop_length, max_loop_length, min_stem_length, max_stem_length,
            min_AT_tail_length, miss_rate, range_u, keep_multi):
        self.TransTermHP_path = TransTermHP_path
        self.expterm_path = expterm_path
        self.RNAfold_path = RNAfold_path
        self.out_folder = out_folder
        self.fastas = fasta_folder
        self.gffs = annotation_folder
        self.trans = transcript_folder
        self.srnas = srna
        self.stat = statistics
        self.tex_wigs = tex_wig_folder
        self.frag_wigs = frag_wig_folder
        self.decrease = decrease
        self.cutoff_coverage = highest_coverage
        self.fuzzy = fuzzy_detect_coverage
        self.fuzzy_up_ta = fuzzy_within_transcript
        self.fuzzy_down_ta = fuzzy_downstream_transcript
        self.fuzzy_up_gene = fuzzy_within_gene
        self.fuzzy_down_gene = fuzzy_downstream_gene
        self.hp_folder = transtermhp_folder
        self.tlibs = self._deal_multi_inputs(tex_notex_libs, "str", None, None)
        self.flibs = self._deal_multi_inputs(frag_libs, "str", None, None)
        self.libs = self._check_libs(self.tlibs, self.flibs)
        self.tex_notex = tex_notex
        self.replicates_tex = replicates_tex
        self.replicates_frag = replicates_frag
        self.replicates = self._check_replicates(
                replicates_tex, replicates_frag)
        self.table_best = table_best
        self.min_loop = min_loop_length
        self.max_loop = max_loop_length
        self.min_stem = min_stem_length
        self.max_stem = max_stem_length
        self.at_tail = min_AT_tail_length
        self.miss_rate = miss_rate
        self.range_u = range_u
        self.keep_multi = keep_multi
        self = self._parser_combine_wigs("terminator")
        return self

    def container_transcript(
            self, frag_wig_path, tex_wig_path, tex_notex, length,
            annotation_folder, height, width, tolerance, tolerance_coverage,
            replicates_tex, replicates_frag, transcript_assembly_output_folder,
            compare_TSS, compare_genome_annotation, TSS_fuzzy,
            tex_treated_libs, fragmented_libs, compare_feature_genome,
            table_best, terminator_folder, fuzzy_term, max_dist):
        self.frag_wigs = frag_wig_path
        self.tex_wigs = tex_wig_path
        self.tex = tex_notex
        self.length = length
        self.gffs = annotation_folder
        self.height = height
        self.width = width
        self.tolerance = tolerance
        self.low_cutoff = tolerance_coverage
        self.replicates_tex = replicates_tex
        self.replicates_frag = replicates_frag
        self.replicates = self._check_replicates(
                replicates_tex, replicates_frag)
        self.out_folder = transcript_assembly_output_folder
        self.compare_tss = compare_TSS
        self.compare_cds = compare_genome_annotation
        self.fuzzy = TSS_fuzzy
        self.tlibs = self._deal_multi_inputs(tex_treated_libs, "str", None,
                                             None)
        self.flibs = self._deal_multi_inputs(fragmented_libs, "str", None,
                                             None)
        self.libs = self._check_libs(self.tlibs, self.flibs)
        self.c_feature = self._deal_multi_inputs(compare_feature_genome, "str",
                                                 None, None)
        self.table_best = table_best
        self.terms = terminator_folder
        self.fuzzy_term = fuzzy_term
        self.max_dist = max_dist
        self = self._parser_combine_wigs("transcript")
        return self

    def container_utr(self, tss_folder, annotation_folder,
                      transcript_assembly_folder, terminator_folder,
                      terminator_fuzzy, utr_folder, tss_source, base_5utr,
                      length, base_3utr, fuzzy_3utr, fuzzy_5utr):
        self.tsss = tss_folder
        self.gffs = annotation_folder
        self.trans = transcript_assembly_folder
        self.terms = terminator_folder
        self.fuzzy = terminator_fuzzy
        self.out_folder = utr_folder
        self.source = tss_source
        self.base_5utr = base_5utr
        self.base_3utr = base_3utr
        self.length = length
        self.fuzzy_3utr = fuzzy_3utr
        self.fuzzy_5utr = fuzzy_5utr
        return self

    def container_srna(
            self, Vienna_folder, Vienna_utils, blast_plus_folder,
            ps2pdf14_path, srna_folder, UTR_derived_sRNA, annotation_folder,
            TSS_folder, transcript_assembly_folder, TSS_intergenic_fuzzy,
            TSS_5UTR_fuzzy, TSS_3UTR_fuzzy, TSS_interCDS_fuzzy, import_info,
            tex_wig_folder, frag_wig_folder, processing_site_folder,
            fasta_folder, mountain_plot, nr_format, srna_format,
            sRNA_database_path, nr_database_path, cutoff_energy,
            run_intergenic_TEX_coverage, run_intergenic_noTEX_coverage,
            run_intergenic_fragmented_coverage, break_tran,
            run_antisense_TEX_coverage, run_antisense_noTEX_coverage,
            run_antisense_fragmented_coverage, run_utr_TEX_coverage,
            run_utr_noTEX_coverage, run_utr_fragmented_coverage, max_length,
            min_length, tex_notex_libs, frag_libs, replicates_tex,
            replicates_frag, tex_notex, blast_e_nr, blast_e_srna,
            detect_sRNA_in_CDS, table_best, decrease_intergenic, decrease_utr,
            fuzzy_intergenic, fuzzy_utr, cutoff_nr_hit, sORF,
            overlap_percent_CDS, terminator_folder, terminator_fuzzy_in_sRNA,
            terminator_fuzzy_out_sRNA, ignore_hypothetical_protein, TSS_source,
            min_utr_coverage, promoter_table, ranking_promoter, promoter_name):
        self.vienna_path = Vienna_folder
        self.vienna_util = Vienna_utils
        self.blast_path = blast_plus_folder
        self.ps2pdf14_path = ps2pdf14_path
        self.out_folder = srna_folder
        self.utr_srna = UTR_derived_sRNA
        self.gffs = annotation_folder
        self.tss_folder = TSS_folder
        self.trans = transcript_assembly_folder
        self.fuzzy_inter_tss = TSS_intergenic_fuzzy
        self.fuzzy_5utr_tss = TSS_5UTR_fuzzy
        self.fuzzy_3utr_tss = TSS_3UTR_fuzzy
        self.fuzzy_intercds_tss = TSS_interCDS_fuzzy
        self.fuzzy_tsss = {"5utr": self.fuzzy_5utr_tss,
                           "3utr": self.fuzzy_3utr_tss,
                           "interCDS": self.fuzzy_intercds_tss,
                           "inter": self.fuzzy_inter_tss}
        self.import_info = self._deal_multi_inputs(import_info, "str",
                                                   None, None)
        self.tex_wigs = tex_wig_folder
        self.frag_wigs = frag_wig_folder
        self.pro_folder = processing_site_folder
        self.fastas = fasta_folder
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
        self.tlibs = self._deal_multi_inputs(tex_notex_libs, "str", None, None)
        self.flibs = self._deal_multi_inputs(frag_libs, "str", None, None)
        self.libs = self._check_libs(self.tlibs, self.flibs)
        self.replicates_tex = replicates_tex
        self.replicates_frag = replicates_frag
        self.replicates = self._check_replicates(
                replicates_tex, replicates_frag)
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
        self.sorf_file = sORF
        self.cutoff_overlap = overlap_percent_CDS
        self.terms = terminator_folder
        self.fuzzy_b = terminator_fuzzy_in_sRNA
        self.fuzzy_a = terminator_fuzzy_out_sRNA
        self.hypo = ignore_hypothetical_protein
        self.tss_source = TSS_source
        self.min_utr = min_utr_coverage
        self.promoter_table = promoter_table
        if ranking_promoter < 1:
            print("Error: --ranking_time_promoter must larger than 1...")
            sys.exit()
        self.rank_promoter = ranking_promoter
        self.promoter_name = self._deal_multi_inputs(promoter_name, "str",
                                                     None, None)
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
#        args_srna.wigs_f = wigs_f
#        args_srna.wigs_r = wigs_r
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
#        args_srna.wig_fs = wig_fs
#        args_srna.wig_rs = wig_rs
        args_srna.out = out
        args_srna.out_t = out_t
        args_srna.texs = texs
        args_srna.utrs = []
        args_srna.srnas = []
        return args_srna

    def container_sorf(self, sorf_folder, UTR_derived_sORF, transcript_folder,
                       annotation_folder, TSS_folder, utr_length, min_length,
                       max_length, tex_wig_folder, frag_wig_folder,
                       cutoff_intergenic_coverage, cutoff_antisense_coverage,
                       cutoff_5utr_coverage, cutoff_3utr_coverage,
                       cutoff_interCDS_coverage, fasta_folder, tex_notex_libs,
                       frag_libs, tex_notex, replicates_tex, replicates_frag,
                       table_best, sRNA_folder, start_codon, stop_codon,
                       cutoff_background, fuzzy_rbs, rbs_not_after_TSS,
                       print_all_combination, best_no_sRNA, best_no_TSS,
                       ignore_hypothetical_protein, min_rbs_distance,
                       max_rbs_distance):
        self.out_folder = sorf_folder
        self.utr_detect = UTR_derived_sORF
        self.trans = transcript_folder
        self.gffs = annotation_folder
        self.tsss = TSS_folder
        self.utr_length = utr_length
        self.min_len = min_length
        self.max_len = max_length
        self.tex_wigs = tex_wig_folder
        self.frag_wigs = frag_wig_folder
        self.cutoff_inter = cutoff_intergenic_coverage
        self.cutoff_anti = cutoff_antisense_coverage
        self.cutoff_5utr = cutoff_5utr_coverage
        self.cutoff_3utr = cutoff_3utr_coverage
        self.cutoff_intercds = cutoff_interCDS_coverage
        self.fastas = fasta_folder
        self.tlibs = self._deal_multi_inputs(tex_notex_libs, "str", None, None)
        self.flibs = self._deal_multi_inputs(frag_libs, "str", None, None)
        self.libs = self._check_libs(self.tlibs, self.flibs)
        self.tex_notex = tex_notex
        self.replicates_tex = replicates_tex
        self.replicates_frag = replicates_frag
        self.replicates = self._check_replicates(
                replicates_tex, replicates_frag)
        self.table_best = table_best
        self.srnas = sRNA_folder
        self.start_codon = self._deal_multi_inputs(start_codon, "str",
                                                   None, None)
        self.stop_codon = self._deal_multi_inputs(stop_codon, "str",
                                                  None, None)
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

    def container_srna_target(self, Vienna_folder, annotation_path, fasta_path,
                              sRNA_path, query_sRNA, program,
                              interaction_length, window_size_target,
                              span_target, window_size_srna, span_srna,
                              unstructured_region_RNAplex_target,
                              unstructured_region_RNAplex_srna,
                              unstructured_region_RNAup, energy_threshold,
                              duplex_distance, top, starget_output_folder,
                              process_rnaplex, process_rnaup, continue_rnaup,
                              potential_target_start, potential_target_end,
                              target_feature):
        self.vienna_path = Vienna_folder
        self.gffs = annotation_path
        self.fastas = fasta_path
        self.srnas = sRNA_path
        self.query = self._deal_multi_inputs(query_sRNA, "str", None, None)
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
        self.features = self._deal_multi_inputs(target_feature, "str",
                                                None, None)
        return self

    def container_goterm(self, annotation_path, goterm_output_folder,
                         UniProt_id, go_obo, goslim_obo, transcript_path):
        self.gffs = annotation_path
        self.out_folder = goterm_output_folder
        self.uniprot = UniProt_id
        self.go = go_obo
        self.goslim = goslim_obo
        self.trans = transcript_path
        return self

    def container_sublocal(self, Psortb_path, gff_path, fasta_path,
                           bacteria_type, difference_multi, merge_to_gff,
                           sublocal_output_folder, transcript_path):
        self.psortb_path = Psortb_path
        self.gffs = gff_path
        self.fastas = fasta_path
        self.gram = bacteria_type
        self.fuzzy = difference_multi
        self.merge = merge_to_gff
        self.out_folder = sublocal_output_folder
        self.trans = transcript_path
        return self

    def container_ppi(self, gff_path, proteinID_strains, without_strain_pubmed,
                      species_STRING, score, ppi_output_folder, node_size,
                      query):
        self.ptts = gff_path
        self.strains = self._deal_multi_inputs(proteinID_strains, "str",
                                               None, None)
        self.no_specific = without_strain_pubmed
        self.species = species_STRING
        self.score = score
        self.out_folder = ppi_output_folder
        self.size = node_size
        self.querys = self._deal_multi_inputs(query, "str", None, None)
        return self

    def container_promoter(self, MEME_path, promoter_output_folder, tex_libs,
                           TSS_folder, fasta_folder, num_motif, nt_before_TSS,
                           motif_width, TSS_source, tex_wig_path,
                           annotation_folder, combine_all, e_value, para):
        self.meme_path = MEME_path
        self.output_folder = promoter_output_folder
        self.input_libs = self._deal_multi_inputs(tex_libs, "str", None, None)
        self.libs = self.input_libs
        self.tsss = TSS_folder
        self.fastas = fasta_folder
        self.num_motif = num_motif
        self.nt_before = nt_before_TSS
        self.widths = self._deal_multi_inputs(motif_width, "str", None, None)
        self.source = TSS_source
        self.tex_wigs = tex_wig_path
        self.frag_wigs = None
        self.gffs = annotation_folder
        self.combine = combine_all
        self.e_value = e_value
        self.para = para
        self = self._parser_combine_wigs("promoter")
        return self

    def container_operon(self, TSS_folder, annotation_folder,
                         transcript_folder, UTR5_folder, UTR3_folder,
                         term_folder, TSS_fuzzy, term_fuzzy, min_length,
                         statistics, operon_output_folder, combine_gff,
                         operon_statistics_folder):
        self.tsss = TSS_folder
        self.gffs = annotation_folder
        self.trans = transcript_folder
        self.utr5s = UTR5_folder
        self.utr3s = UTR3_folder
        self.terms = term_folder
        self.tss_fuzzy = TSS_fuzzy
        self.term_fuzzy = term_fuzzy
        self.length = min_length
        self.statistics = statistics
        self.output_folder = operon_output_folder
        self.combine = combine_gff
        self.stat_folder = operon_statistics_folder
        return self

    def container_snp(self, samtools_path, bcftools_path, bam_type, min_sample,
                      program, fasta_path, tex_bam_path, frag_bam_path,
                      quality, read_depth_range, snp_output_folder,
                      indel_fraction, chrom, rg, caller, filters, DP4_cutoff):
        self.samtools_path = samtools_path
        self.bcftools_path = bcftools_path
        self.types = bam_type
        self.program = self._deal_multi_inputs(program, "str", None, None)
        self.fastas = fasta_path
        self.normal_bams = tex_bam_path
        self.frag_bams = frag_bam_path
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
        self.filters = filters.split(",")
        self.dp4_sum = DP4_cutoff.split(",")[0]
        self.dp4_frac = DP4_cutoff.split(",")[-1]
        self.min_sample = min_sample
        return self

    def container_circrna(self, align, process, fasta_path, annotation_path,
                          tex_bam_path, fragmented_bam_path, read_folder,
                          circrna_stat_folder, support_reads,
                          segemehl_folder, samtools_path, start_ratio,
                          end_ratio, ignore_hypothetical_protein, out_folder):
        self.align = align
        self.cores = process
        self.fastas = fasta_path
        self.gffs = annotation_path
        self.normal_bams = tex_bam_path
        self.frag_bams = fragmented_bam_path
        self.read_folder = read_folder
        self.stat_folder = circrna_stat_folder
        self.support = support_reads
        self.segemehl_path = segemehl_folder
        self.samtools_path = samtools_path
        self.start_ratio = start_ratio
        self.end_ratio = end_ratio
        self.hypo = ignore_hypothetical_protein
        self.output_folder = out_folder
        return self

    def container_ribos(self, program, thermo_ID, infernal_path, riboswitch_ID,
                        gff_path, fasta_path, tss_path, transcript_path, Rfam,
                        ribos_output_folder, thermo_output_folder, e_value,
                        output_all, database_folder, fuzzy, start_codon,
                        min_dist_rbs, max_dist_rbs, fuzzy_rbs, UTR_length):
        self.program = program
        self.thermo_id = thermo_ID
        self.infernal_path = infernal_path
        self.ribos_id = riboswitch_ID
        self.gffs = gff_path
        self.fastas = fasta_path
        self.tsss = tss_path
        self.trans = transcript_path
        self.rfam = Rfam
        self.ribos_out_folder = ribos_output_folder
        self.thermo_out_folder = thermo_output_folder
        self.e_value = e_value
        self.output_all = output_all
        self.database = database_folder
        self.fuzzy = fuzzy
        self.start_codons = self._deal_multi_inputs(start_codon, "str",
                                                    None, None)
        self.start_rbs = min_dist_rbs
        self.end_rbs = max_dist_rbs
        self.fuzzy_rbs = fuzzy_rbs
        self.utr = UTR_length
        return self

    def container_cris(self, fasta_path, gff_path, CRT_path, window_size,
                       min_number_repeat, min_length_repeat, Max_length_repeat,
                       min_length_spacer, Max_length_spacer, cris_out_folder,
                       ignore_hypo):
        self.fastas = fasta_path
        self.gffs = gff_path
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

    def container_screen(self, main_gff, side_gffs, fasta, frag_wig_folder,
                         tex_wig_folder, height, tex_libs, frag_libs, present,
                         output_folder):
        self.main_gff = main_gff
        self.side_gffs = self._deal_multi_inputs(side_gffs, "str", None, None)
        self.fasta = fasta
        self.frag_wigs = frag_wig_folder
        self.tex_wigs = tex_wig_folder
        self.height = height
        self.tlibs = self._deal_multi_inputs(tex_libs, "str", None, None)
        self.flibs = self._deal_multi_inputs(frag_libs, "str", None, None)
        self.present = present
        self.output_folder = output_folder
        return self
