class Paths(object):
    '''Setup the folders of ANNOgesic'''

    def __init__(self, base_path="."):
        self.base_path = base_path
        self._set_folder_names()

    def _set_folder_names(self):
        """Set the name of folders used in a project."""
        self.input_folder = "%s/input" % (self.base_path)
        self.output_folder = "%s/output" % (self.base_path)
        self._set_input_folder_names()
        self._set_reference_input_folder_names()
        self._set_wig_folder_names()
        self._set_bam_folder_names()
        self._set_output_folder_names()
        self._set_target_folder_names()
        self._set_tsspredator_folder_names()
        self._set_transterm_folder_names()
        self._set_processing_folder_names()
        self._set_transcript_folder_names()
        self._set_utr_folder_names()
        self._set_srna_folder_names()
        self._set_sorf_folder_names()
        self._set_operon_folder_names()
        self._set_circrna_folder_names()
        self._set_goterm_folder_names()
        self._set_starget_folder_names()
        self._set_snp_folder_names()
        self._set_ppi_folder_names()
        self._set_sublocal_folder_names()
        self._set_ribos_folder_names()
        self._set_thermo_folder_names()
        self._set_crispr_folder_names()
        self.version_path = "%s/used_annogesic_version.txt" % (self.base_path)

    def _set_input_folder_names(self):
        self.reference_input_folder = "%s/references" % self.input_folder
        self.wig_folder = "%s/wigs" % self.input_folder
        self.mutation_table_folder = "%s/mutation_tables" % self.input_folder
        self.database_folder = "%s/databases" % self.input_folder
        self.manual_TSS_folder = "%s/manual_TSSs" % self.input_folder
        self.manual_pro_folder = "%s/manual_processing_sites" % (
             self.input_folder)
        self.read_folder = "%s/reads" % self.input_folder
        self.bam_folder = "%s/BAMs" % self.input_folder
        self.riborfam_folder = "%s/riboswitch_ID_file" % self.input_folder
        self.thermorfam_folder = "%s/RNA_thermometer_ID_file" % self.input_folder

    def _set_output_folder_names(self):
        self.target_folder = "%s/updated_references" % self.output_folder
        self.ratt_folder = "%s/annotation_transfer" % self.output_folder
        self.tsspredator_folder = "%s/TSSs" % self.output_folder
        self.utr_folder = "%s/UTRs" % self.output_folder
        self.transterm_folder = "%s/terminators" % self.output_folder
        self.transcript_output_folder = (
                "%s/transcripts" % self.output_folder)
        self.processing_site_folder = "%s/processing_sites" % self.output_folder
        self.srna_folder = "%s/sRNAs" % self.output_folder
        self.sorf_folder = "%s/sORFs" % self.output_folder
        self.promoter_output_folder = "%s/promoters" % (
             self.output_folder)
        self.operon_output_folder = "%s/operons" % self.output_folder
        self.circrna_output_folder = "%s/circRNAs" % self.output_folder
        self.goterm_output_folder = "%s/GO_terms" % self.output_folder
        self.starget_output_folder = "%s/sRNA_targets" % self.output_folder
        self.snp_output_folder = "%s/SNP_calling" % self.output_folder
        self.ppi_output_folder = "%s/PPI_networks" % self.output_folder
        self.sublocal_output_folder = "%s/subcellular_localization" % (
             self.output_folder)
        self.ribos_output_folder = "%s/riboswitches" % self.output_folder
        self.thermo_output_folder = "%s/RNA_thermometers" % self.output_folder
        self.crispr_output_folder = "%s/crisprs" % self.output_folder

    def _set_transcript_folder_names(self):
        self.transcript_base_folder = "%s/transcripts" % (
             self.output_folder)
        self.transcript_gff_folder = "%s/gffs" % self.transcript_base_folder
        self.transcript_stat_folder = "%s/statistics" % (
             self.transcript_base_folder)
        self.transcript_table_folder = "%s/tables" % (
             self.transcript_base_folder)

    def _set_reference_input_folder_names(self):
        self.reference_base_folder = "%s/references" % self.input_folder
        self.ref_annotation_folder = "%s/annotations" % (
            self.reference_base_folder)
        self.ref_fasta_folder = "%s/fasta_files" % (
            self.reference_base_folder)

    def _set_wig_folder_names(self):
        self.wig_base_folder = "%s/wigs" % self.input_folder
        self.frag_folder = "%s/fragment" % (
            self.wig_base_folder)
        self.tex_folder = "%s/tex_notex" % (
            self.wig_base_folder)

    def _set_bam_folder_names(self):
        self.bam_base_folder = "%s/BAMs" % self.input_folder
        self.bam_ref_folder = "%s/BAMs_map_related_genomes" % self.bam_base_folder
        self.bam_tar_folder = "%s/BAMs_map_reference_genomes" % self.bam_base_folder
        self.bam_ref_frag_folder = "%s/fragment" % (
            self.bam_ref_folder)
        self.bam_tar_frag_folder = "%s/fragment" % (
            self.bam_tar_folder)
        self.bam_ref_tex_folder = "%s/tex_notex" % (
            self.bam_ref_folder)
        self.bam_tar_tex_folder = "%s/tex_notex" % (
            self.bam_tar_folder)

    def _set_target_folder_names(self):
        self.target_base_folder = "%s/updated_references" % self.output_folder
        self.tar_fasta_folder = "%s/fasta_files" % (
            self.target_base_folder)
        self.tar_annotation_folder = "%s/annotations" % (
            self.target_base_folder)

    def _set_tsspredator_folder_names(self):
        self.tsspredator_base_folder = "%s/TSSs" % self.output_folder
        self.tss_to_gff_folder = "%s/gffs" % (
            self.tsspredator_base_folder)
        self.tss_statistics_folder = "%s/statistics" % (
            self.tsspredator_base_folder)
        self.tss_Master_folder = "%s/MasterTables" % (
            self.tsspredator_base_folder)
        self.tss_config_folder = "%s/configs" % (
            self.tsspredator_base_folder)

    def _set_processing_folder_names(self):
        self.processing_base_folder = "%s/processing_sites" % self.output_folder
        self.processing_to_gff_folder = "%s/gffs" % (
            self.processing_base_folder)
        self.processing_statistics_folder = "%s/statistics" % (
            self.processing_base_folder)
        self.processing_screenshot_folder = "%s/screenshots" % (
            self.processing_base_folder)
        self.processing_Master_folder = "%s/MasterTables" % (
            self.processing_base_folder)
        self.processing_config_folder = "%s/configs" % (
            self.processing_base_folder)

    def _set_transterm_folder_names(self):
        self.transterm_base_folder = "%s/terminators" % self.output_folder
        self.term_to_gff_folder = "%s/gffs" % (
            self.transterm_base_folder)
        self.term_to_table_folder = "%s/tables" % (
            self.transterm_base_folder)
        self.transtermhp_folder = "%s/transtermhp_results" % (
            self.transterm_base_folder)
        self.term_statistics_folder = "%s/statistics" % (
            self.transterm_base_folder)

    def _set_utr_folder_names(self):
        self.utr_base_folder = "%s/UTRs" % self.output_folder
        self.utr5_folder = "%s/5UTRs" % (
            self.utr_base_folder)
        self.utr3_folder = "%s/3UTRs" % (
            self.utr_base_folder)
        self.utr3_stat_folder = "%s/statistics" % (
             self.utr3_folder)
        self.utr3_gff_folder = "%s/gffs" % (
             self.utr3_folder)
        self.utr5_stat_folder = "%s/statistics" % (
             self.utr5_folder)
        self.utr5_gff_folder = "%s/gffs" % (
             self.utr5_folder)

    def _set_srna_folder_names(self):
        self.srna_base_folder = "%s/sRNAs" % self.output_folder
        self.srna_gff_folder = "%s/gffs" % (
            self.srna_base_folder)
        self.srna_table_folder = "%s/tables" % (
            self.srna_base_folder)
        self.srna_plot_folder = "%s/figs" % (
            self.srna_base_folder)
        self.srna_sec_plot_folder = "%s/figs/sec_plots" % (
            self.srna_base_folder)
        self.srna_dot_plot_folder = "%s/figs/dot_plots" % (
            self.srna_base_folder)
        self.srna_mountain_folder = "%s/figs/mountain_plots" % (
            self.srna_base_folder)
        self.srna_blast_folder = "%s/blast_results_and_misc" % (
            self.srna_base_folder)
        self.srna_stat_folder = "%s/statistics" % (
            self.srna_base_folder)
        self.srna_gff_class_folder = "%s/for_classes" % (
            self.srna_gff_folder)
        self.srna_gff_best_folder = "%s/best_candidates" % (
            self.srna_gff_folder)
        self.srna_gff_all_folder = "%s/all_candidates" % (
            self.srna_gff_folder)
        self.srna_table_class_folder = "%s/for_classes" % (
            self.srna_table_folder)
        self.srna_table_best_folder = "%s/best_candidates" % (
            self.srna_table_folder)
        self.srna_table_all_folder = "%s/all_candidates" % (
            self.srna_table_folder)

    def _set_sorf_folder_names(self):
        self.sorf_base_folder = "%s/sORFs" % self.output_folder
        self.sorf_gff_folder = "%s/gffs" % (
            self.sorf_base_folder)
        self.sorf_table_folder = "%s/tables" % (
            self.sorf_base_folder)
        self.sorf_stat_folder = "%s/statistics" % (
            self.sorf_base_folder)
        self.sorf_gff_best_folder = "%s/best_candidates" % (
            self.sorf_gff_folder)
        self.sorf_gff_all_folder = "%s/all_candidates" % (
            self.sorf_gff_folder)
        self.sorf_table_best_folder = "%s/best_candidates" % (
            self.sorf_table_folder)
        self.sorf_table_all_folder = "%s/all_candidates" % (
            self.sorf_table_folder)

    def _set_operon_folder_names(self):
        self.operon_base_folder = "%s/operons" % self.output_folder
        self.operon_gff_folder = "%s/gffs" % (
             self.operon_base_folder)
        self.operon_table_folder = "%s/tables" % (
             self.operon_base_folder)
        self.operon_statistics_folder = "%s/statistics" % (
            self.operon_base_folder)

    def _set_circrna_folder_names(self):
        self.circrna_base_folder = "%s/circRNAs" % self.output_folder
        self.circrna_align_folder = "%s/segemehl_alignment_files" % (
             self.circrna_base_folder)
        self.circrna_splice_folder = "%s/segemehl_splice_results" % (
             self.circrna_base_folder)
        self.circrna_circ_folder = "%s/circRNA_tables" % (
             self.circrna_base_folder)
        self.circrna_stat_folder = "%s/statistics" % (
             self.circrna_base_folder)
        self.circrna_gff_folder = "%s/gffs" % (
             self.circrna_base_folder)

    def _set_goterm_folder_names(self):
        self.goterm_base_folder = "%s/GO_terms" % self.output_folder
        self.goterm_all_folder = "%s/all_CDSs" % self.goterm_base_folder
        self.goterm_express_folder = "%s/expressed_CDSs" % (
             self.goterm_base_folder)
        self.goterm_express_result_folder = "%s/GO_term_results" % (
             self.goterm_express_folder)
        self.goterm_express_stat_folder = "%s/statistics" % (
             self.goterm_express_folder)
        self.goterm_all_result_folder = "%s/GO_term_results" % (
             self.goterm_all_folder)
        self.goterm_all_stat_folder = "%s/statistics" % (
             self.goterm_all_folder)

    def _set_starget_folder_names(self):
        self.starget_base_folder = "%s/sRNA_targets" % self.output_folder
        self.starget_RNAplex_folder = "%s/RNAplex_results" % (
             self.starget_base_folder)
        self.starget_RNAup_folder = "%s/RNAup_results" % (
             self.starget_base_folder)
        self.starget_IntaRNA_folder = "%s/IntaRNA_results" % (
             self.starget_base_folder)
        self.starget_merge_folder = "%s/merged_results" % (
             self.starget_base_folder)
        self.starget_srna_seq_folder = "%s/sRNA_seqs" % (
             self.starget_base_folder)
        self.starget_target_seq_folder = "%s/target_seqs" % (
             self.starget_base_folder)

    def _set_snp_folder_names(self):
        self.snp_base_folder = "%s/SNP_calling" % self.output_folder
        self.ref_snp_folder = "%s/compare_related_and_reference_genomes" % self.snp_base_folder
        self.tar_snp_folder = "%s/mutations_of_reference_genomes" % self.snp_base_folder
        self.snp_ref_stat_folder = "%s/statistics" % (
             self.ref_snp_folder)
        self.snp_tar_stat_folder = "%s/statistics" % (
             self.tar_snp_folder)
        self.snp_ref_table_folder = "%s/SNP_tables" % (
             self.ref_snp_folder)
        self.snp_tar_table_folder = "%s/SNP_tables" % (
             self.tar_snp_folder)
        self.snp_ref_raw_folder = "%s/SNP_raw_outputs" % (
             self.ref_snp_folder)
        self.snp_tar_raw_folder = "%s/SNP_raw_outputs" % (
             self.tar_snp_folder)
        self.snp_ref_seq_folder = "%s/seqs" % (
             self.ref_snp_folder)
        self.snp_tar_seq_folder = "%s/seqs" % (
             self.tar_snp_folder)
        self.snp_ref_seq_extend_BAQ_folder = "%s/extend_BAQ" % (
             self.snp_ref_seq_folder)
        self.snp_tar_seq_extend_BAQ_folder = "%s/extend_BAQ" % (
             self.snp_tar_seq_folder)
        self.snp_ref_seq_with_BAQ_folder = "%s/with_BAQ" % (
             self.snp_ref_seq_folder)
        self.snp_tar_seq_with_BAQ_folder = "%s/with_BAQ" % (
             self.snp_tar_seq_folder)
        self.snp_ref_seq_without_BAQ_folder = "%s/without_BAQ" % (
             self.snp_ref_seq_folder)
        self.snp_tar_seq_without_BAQ_folder = "%s/without_BAQ" % (
             self.snp_tar_seq_folder)

    def _set_ppi_folder_names(self):
        self.ppi_base_folder = "%s/PPI_networks" % self.output_folder
        self.ppi_all_folder = "%s/all_results" % (
             self.ppi_base_folder)
        self.ppi_best_folder = "%s/best_results" % (
             self.ppi_base_folder)
        self.ppi_fig_folder = "%s/figures" % (
             self.ppi_base_folder)

    def _set_sublocal_folder_names(self):
        self.sublocal_base_folder = "%s/subcellular_localization" % (
             self.output_folder)
        self.sublocal_all_folder = "%s/all_CDSs" % self.sublocal_base_folder
        self.sublocal_express_folder = "%s/expressed_CDSs" % (
             self.sublocal_base_folder)
        self.sublocal_all_results_folder = "%s/psortb_results" % (
             self.sublocal_all_folder)
        self.sublocal_all_stat_folder = "%s/statistics" % (
             self.sublocal_all_folder)
        self.sublocal_express_results_folder = "%s/psortb_results" % (
             self.sublocal_express_folder)
        self.sublocal_express_stat_folder = "%s/statistics" % (
             self.sublocal_express_folder)

    def _set_ribos_folder_names(self):
        self.ribos_base_folder = "%s/riboswitches" % self.output_folder
        self.ribos_gff_folder = "%s/gffs" % (
             self.ribos_base_folder)
        self.ribos_stat_folder = "%s/statistics" % (
             self.ribos_base_folder)
        self.ribos_table_folder = "%s/tables" % (
             self.ribos_base_folder)
        self.ribos_rfam_folder = "%s/scan_Rfam_results" % (
             self.ribos_base_folder)

    def _set_thermo_folder_names(self):
        self.thermo_base_folder = "%s/RNA_thermometers" % self.output_folder
        self.thermo_gff_folder = "%s/gffs" % (
             self.thermo_base_folder)
        self.thermo_stat_folder = "%s/statistics" % (
             self.thermo_base_folder)
        self.thermo_table_folder = "%s/tables" % (
             self.thermo_base_folder)
        self.thermo_rfam_folder = "%s/scan_Rfam_results" % (
             self.thermo_base_folder)

    def _set_crispr_folder_names(self):
        self.crispr_base_folder = "%s/crisprs" % self.output_folder
        self.crispr_gff_folder = "%s/gffs" % (
             self.crispr_base_folder)
        self.crispr_stat_folder = "%s/statistics" % (
             self.crispr_base_folder)
        self.crispr_data_folder = "%s/CRT_results" % (
             self.crispr_base_folder)

    def required_folders(self, folder_type):
        if (folder_type == "root"):
            return (self.required_base_folders() +
                    self.required_input_folders() +
                    self.required_reference_input_folders() +
                    self.required_wig_folders() +
                    self.required_bam_folders())
        else:
            return (self.required_base_folders() +
                    self.required_input_folders() +
                    self.required_reference_input_folders() +
                    self.required_wig_folders() +
                    self.required_bam_folders() +
                    self.required_output_folders(folder_type))

    def required_base_folders(self):
        return [self.input_folder, self.output_folder]

    def required_input_folders(self):
        return [self.reference_input_folder, self.wig_folder,
                self.mutation_table_folder, self.read_folder,
                self.bam_folder,
                self.database_folder, self.manual_TSS_folder,
                self.manual_pro_folder, self.riborfam_folder,
                self.thermorfam_folder]

    def required_output_folders(self, folder_type):
        folder_dict = {"get_target_fasta": (
                           [self.target_folder] +
                           self.required_target_folders()),
                       "annotation_transfer": (
                           [self.ratt_folder] +
                           self.required_target_folders()),
                       "TSS": (
                           [self.tsspredator_folder] +
                           self.required_tsspredator_folders()),
                       "processing": (
                           [self.processing_site_folder] +
                           self.required_processing_folders()),
                       "terminator": (
                           [self.transterm_folder] +
                           self.required_transterm_folders()),
                       "transcript": (
                           [self.transcript_output_folder] +
                           self.required_transcript_folders()),
                       "utr": [self.utr_folder] + self.required_utr_folders(),
                       "srna": (
                           [self.srna_folder] + self.required_srna_folders()),
                       "sorf": (
                           [self.sorf_folder] + self.required_sorf_folders()),
                       "promoter": [self.promoter_output_folder],
                       "circrna": (
                           [self.circrna_output_folder] +
                           self.required_circrna_folders()),
                       "go_term": (
                           [self.goterm_output_folder] +
                           self.required_goterm_folders()),
                       "srna_target": (
                           [self.starget_output_folder] +
                           self.required_starget_folders()),
                       "snp": (
                           [self.snp_output_folder] +
                           self.required_snp_folders()),
                       "ppi_network": (
                           [self.ppi_output_folder] +
                           self.required_ppi_folders()),
                       "subcellular_localization": (
                           [self.sublocal_output_folder] +
                           self.required_sublocal_folders()),
                       "riboswitch": (
                           [self.ribos_output_folder] +
                           self.required_ribos_folders()),
                       "thermometer": (
                           [self.thermo_output_folder] +
                           self.required_thermo_folders()),
                       "crispr": (
                           [self.crispr_output_folder] +
                           self.required_crispr_folders()),
                       "operon": (
                           [self.operon_output_folder] +
                           self.required_operon_folders())}
        return folder_dict[folder_type]

    def required_reference_input_folders(self):
        return [self.ref_annotation_folder, self.ref_fasta_folder]

    def required_wig_folders(self):
        return [self.tex_folder, self.frag_folder]

    def required_bam_folders(self):
        return [self.bam_ref_folder, self.bam_tar_folder,
                self.bam_ref_tex_folder, self.bam_tar_tex_folder,
                self.bam_ref_frag_folder, self.bam_tar_frag_folder]

    def required_target_folders(self):
        return [self.tar_annotation_folder, self.tar_fasta_folder]

    def required_tsspredator_folders(self):
        return [self.tss_to_gff_folder, self.tss_statistics_folder,
                self.tss_Master_folder, self.tss_config_folder]

    def required_transterm_folders(self):
        return [self.term_to_gff_folder, self.term_statistics_folder,
                self.transtermhp_folder, self.term_to_table_folder]

    def required_processing_folders(self):
        return [self.processing_to_gff_folder,
                self.processing_statistics_folder,
                self.processing_Master_folder,
                self.processing_config_folder]

    def required_transcript_folders(self):
        return [self.transcript_gff_folder, self.transcript_stat_folder,
                self.transcript_table_folder]

    def required_utr_folders(self):
        return [self.utr5_folder, self.utr3_folder,
                self.utr5_stat_folder, self.utr5_gff_folder,
                self.utr3_stat_folder, self.utr3_gff_folder]

    def required_srna_folders(self):
        return [self.srna_gff_folder, self.srna_plot_folder,
                self.srna_sec_plot_folder, self.srna_dot_plot_folder,
                self.srna_mountain_folder, self.srna_table_folder,
                self.srna_blast_folder, self.srna_stat_folder,
                self.srna_gff_class_folder, self.srna_gff_best_folder,
                self.srna_gff_all_folder, self.srna_table_class_folder,
                self.srna_table_best_folder, self.srna_table_all_folder]

    def required_sorf_folders(self):
        return [self.sorf_gff_folder, self.sorf_table_folder,
                self.sorf_stat_folder, self.sorf_gff_best_folder,
                self.sorf_gff_all_folder, self.sorf_table_best_folder,
                self.sorf_table_all_folder]

    def required_operon_folders(self):
        return [self.operon_gff_folder, self.operon_table_folder,
                self.operon_statistics_folder]

    def required_circrna_folders(self):
        return [self.circrna_align_folder, self.circrna_splice_folder,
                self.circrna_circ_folder, self.circrna_stat_folder,
                self.circrna_gff_folder]

    def required_goterm_folders(self):
        return [self.goterm_all_folder, self.goterm_express_folder,
                self.goterm_express_result_folder,
                self.goterm_express_stat_folder,
                self.goterm_all_result_folder, self.goterm_all_stat_folder]

    def required_starget_folders(self):
        return [self.starget_RNAplex_folder, self.starget_RNAup_folder,
                self.starget_IntaRNA_folder, self.starget_merge_folder,
                self.starget_srna_seq_folder, self.starget_target_seq_folder]

    def required_snp_folders(self):
        return [self.ref_snp_folder, self.tar_snp_folder,
                self.snp_ref_stat_folder, self.snp_tar_stat_folder,
                self.snp_ref_table_folder, self.snp_tar_table_folder,
                self.snp_ref_raw_folder, self.snp_tar_raw_folder,
                self.snp_ref_seq_folder, self.snp_tar_seq_folder,
                self.snp_ref_seq_extend_BAQ_folder,
                self.snp_tar_seq_extend_BAQ_folder,
                self.snp_ref_seq_with_BAQ_folder,
                self.snp_tar_seq_with_BAQ_folder,
                self.snp_ref_seq_without_BAQ_folder,
                self.snp_tar_seq_without_BAQ_folder]

    def required_ppi_folders(self):
        return [self.ppi_all_folder, self.ppi_best_folder,
                self.ppi_fig_folder]

    def required_sublocal_folders(self):
        return [self.sublocal_all_folder, self.sublocal_express_folder,
                self.sublocal_all_results_folder,
                self.sublocal_all_stat_folder,
                self.sublocal_express_results_folder,
                self.sublocal_express_stat_folder]

    def required_ribos_folders(self):
        return [self.ribos_gff_folder, self.ribos_table_folder,
                self.ribos_stat_folder, self.ribos_rfam_folder]

    def required_thermo_folders(self):
        return [self.thermo_gff_folder, self.thermo_table_folder,
                self.thermo_stat_folder, self.thermo_rfam_folder]

    def required_crispr_folders(self):
        return [self.crispr_gff_folder, self.crispr_stat_folder,
                self.crispr_data_folder]
