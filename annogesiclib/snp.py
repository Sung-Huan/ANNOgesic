#!/usr/bin/python
import os	
import sys
import csv
from subprocess import call
from annogesiclib.multiparser import Multiparser
from annogesiclib.seq_editer import Seq_Editer
from annogesiclib.transcript_SNP import snp_detect
from annogesiclib.helper import Helper

class SNP_calling(object):

    def __init__(self, types, out_folder, fastas):
        self.multiparser = Multiparser()
        self.seq_editer = Seq_Editer()
        self.helper = Helper()
        if types == "reference":
            file_type = "compare_reference"
        else:
            file_type = "validate_target"
        self.seq_path = os.path.join(out_folder, file_type, "seqs")
        self.stat_path = os.path.join(out_folder, file_type, "statistics")
        self.fasta_path = os.path.join(fastas, "tmp")
        self.table_folder = os.path.join(out_folder, file_type, "SNP_table")
        self.raw_folder = os.path.join(out_folder, file_type, "SNP_raw_outputs")
        self.tmp_bcf = os.path.join(out_folder, "tmp_bcf")
#        if "whole_reads.bam" in os.listdir(out_folder):
#           self.helper.remove_all_content(out_folder, "whole_read", "file")
        self.whole_bam = os.path.join(out_folder, "whole_reads.bam")
        self.sort_bam = os.path.join(out_folder, "whole_reads_sorted.bam")
        self.header = os.path.join(out_folder, "header")
        self.with_baq = "with_BAQ"
        self.without_baq = "without_BAQ"
        self.extend_baq = "extend_BAQ"

    def _import_bam(self, bam_folder, bams):
        num_bam = 0
        for bam in os.listdir(bam_folder):
            if bam.endswith(".bam"):
                num_bam += 1
                bams.append(os.path.join(bam_folder, bam))
        return num_bam

    def _transcript_snp(self, fasta, snp, out_table_prefix, qual, 
                        seq_path, fasta_prefix, read_depth, stat_path, 
                        types, bam_number, fraction, table_path):
        stat_file = os.path.join(stat_path, 
                   "_".join(["stat", "_".join([fasta_prefix, types]), "SNP.csv"]))
        out_stat = open(stat_file, "w")
        snp_detect(fasta, snp, out_table_prefix, qual, seq_path + fasta_prefix, 
                   read_depth, fraction, bam_number, stat_file)
        self.helper.move_all_content(table_path, stat_path, [".png"]) 

    def _run_program(self, program, samtools_path, bcftools_path, fasta_file, out_bcf, 
                     out_raw_prefix, out_table_prefix, quality, seq_path, prefix, 
                     depth, stat_path, bam_number, fraction, out_folder, table_path):
        if "1" in program:
            print("Running SNP calling with BAQ...")
            call([samtools_path, "mpileup", 
                  "-t", "DP", "-ugf", fasta_file, self.sort_bam,
                  "--ignore-RG"], stdout=out_bcf)
            out_vcf = "_".join([out_raw_prefix, self.with_baq + ".vcf"])
            call([bcftools_path , "call", self.tmp_bcf,
                  "-vmO", "v", "-o", out_vcf])
            self.helper.check_make_folder(os.path.join(seq_path, self.with_baq, prefix))
            self._transcript_snp(fasta_file, out_vcf,
                                 "_".join([out_table_prefix, self.with_baq]), quality,
                                 os.path.join(seq_path, self.with_baq, prefix),
                                 prefix, depth, stat_path, self.with_baq,
                                 bam_number, fraction, table_path)
        if "2" in program:
            print("Running SNP calling without BAQ...")
            call([samtools_path, "mpileup", 
                  "-t", "DP", "-B", "-ugf", fasta_file,
                  self.sort_bam,
                  "--ignore-RG"], stdout=out_bcf)
            out_vcf = "_".join([out_raw_prefix, self.without_baq + ".vcf"])
            call([bcftools_path, "call", self.tmp_bcf,
                  "-vmO", "v", "-o", out_vcf])
            self.helper.check_make_folder(os.path.join(seq_path, self.without_baq, prefix))
            self._transcript_snp(fasta_file, out_vcf,
                                 "_".join([out_table_prefix, self.without_baq]), quality,
                                 os.path.join(seq_path, self.without_baq, prefix),
                                 prefix, depth, stat_path, self.without_baq,
                                 bam_number, fraction, table_path)
        if "3" in program:
            print("Running SNP calling extend BAQ...")
            call([samtools_path, "mpileup",
                  "-t", "DP", "-E", "-ugf", fasta_file,
                  self.sort_bam,
                  "--ignore-RG"], stdout=out_bcf)
            out_vcf = "_".join([out_raw_prefix, self.extend_baq + ".vcf"])
            call([bcftools_path, "call", self.tmp_bcf,
                  "-vmO", "v", "-o", out_vcf])
            self.helper.check_make_folder(os.path.join(seq_path, self.extend_baq, prefix))
            self._transcript_snp(fasta_file, out_vcf,
                                 "_".join([out_table_prefix, self.extend_baq]), quality,
                                 os.path.join(seq_path, self.extend_baq, prefix),
                                 prefix, depth, stat_path, self.extend_baq,
                                 bam_number, fraction, table_path)

    def _detect_fasta(self, fasta):
        detect = False
        if fasta.endswith(".fa"):
            prefix = fasta[:-3]
            detect = True
        elif fasta.endswith(".fna"):
            prefix = fasta[:-4]
            detect = True
        elif fasta.endswith(".fasta"):
            prefix = fasta[:-6]
            detect = True
        return (detect, prefix)

    def _merge_bams(self, normal_bams, frag_bams, samtools_path, out_folder):
        bams = []
        num_normal = 0
        num_frag = 0
        if (frag_bams is False) and \
           (normal_bams is False):
            print("Error: There is no BAMs folders!!")
            sys.exit()
        else:
            if normal_bams is not False:
                num_normal = self._import_bam(normal_bams, bams)
            if frag_bams is not False:
                num_frag = self._import_bam(frag_bams, bams)
        num_bam = num_normal + num_frag
        print("Merge BAM files now ...")
#        command = (" ".join([samtools_path, "merge", 
#                   self.whole_bam, " ".join(bams)]))
#        os.system(command)
        print("Sort BAM file now ...")
#        command = (" ".join([samtools_path, "sort", 
#                   self.whole_bam, self.sort_bam.replace(".bam", "")]))
#        os.system(command)
        return num_bam

    def _modify_header(self, fastas):
        for fasta in os.listdir(fastas):
            if fasta.endswith("fasta") or \
               fasta.endswith("fa") or \
               fasta.endswith("fna"):
                self.seq_editer.modify_header(os.path.join(fastas, fasta))

    def _get_genome_name(self, samtools_path, out_folder):
        command = " ".join([samtools_path, "view", "-H", self.sort_bam])
        os.system(">".join([command, self.header]))
        fh = open(self.header, "r");
        seq_names = []
        for row in csv.reader(fh, delimiter="\t"):
            if row[0] == "@SQ":
                seq_names.append(row[1].split(":")[1])
        return seq_names

    def run_snp_calling(self, samtools_path, bcftools_path, types, program, fastas, 
                        normal_bams, frag_bams, quality, depth, out_folder, fraction):
        if types == "reference":
            file_type = "compare_reference"
        else:
            file_type = "validate_target"
        self.multiparser._parser_fasta(fastas)
        self._modify_header(fastas)
        bam_number = self._merge_bams(normal_bams, frag_bams, samtools_path, out_folder)
        seq_names = self._get_genome_name(samtools_path, out_folder)
        #### running SNP calling
        if ("1" not in program ) and \
           ("2" not in program) and \
           ("3" not in program):
            print("Error:Please assign a correct BAQ type: '1' means 'with_BAQ', '2' means 'with_BAQ' or '3' means 'extend_BAQ'.")
            sys.exit()
        else:
            for fasta in os.listdir(self.fasta_path):
                if (fasta.split(".f")[0] in seq_names):
                    fasta_datas = self._detect_fasta(fasta)
                    detect = fasta_datas[0]
                    prefix = fasta_datas[1]
                    if detect:
                        detect = False
                        print("Computing {0} now ...".format(fasta))
                        out_bcf = open(self.tmp_bcf, "w")
                        self.helper.check_make_folder(os.path.join(self.table_folder, prefix))
                        self.helper.check_make_folder(os.path.join(self.raw_folder, prefix))
                        out_raw_prefix = os.path.join(self.raw_folder, prefix, prefix)
                        out_table_prefix = os.path.join(self.table_folder, prefix, prefix)
                        fasta_file = os.path.join(self.fasta_path, fasta)
                        table_path = os.path.join(self.table_folder, prefix)
                        self._run_program(program, samtools_path, bcftools_path, fasta_file, out_bcf, 
                                          out_raw_prefix, out_table_prefix, quality, self.seq_path, prefix, 
                                          depth, self.stat_path, bam_number, fraction, out_folder, table_path)
                        os.remove(self.tmp_bcf)
        self.helper.remove_tmp(fastas)
        os.remove(self.whole_bam)
        os.remove(self.sort_bam)
        os.remove(self.header)
