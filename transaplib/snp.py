#!/usr/bin/python
import os	
import sys
import csv
from subprocess import call
from transaplib.multiparser import Multiparser
from transaplib.seq_editer import Seq_Editer
from transaplib.transcript_SNP import SNP_detect
from transaplib.helper import Helper

class SNP_calling(object):

    def __init__(self):
        self.multiparser = Multiparser()
        self.seq_editer = Seq_Editer()
        self.helper = Helper()

    def _import_bam(self, bam_folder, bams):
        num_bam = 0
        for bam in os.listdir(bam_folder):
            if bam.endswith(".bam"):
                num_bam += 1
                bams.append(bam_folder + "/" + bam)
        return num_bam

    def _transcript_SNP(self, fasta, snp, out_table_prefix, 
                        qual, seq_path, fasta_prefix, read_depth, 
                        stat_path, types, bam_number, fraction):
        stat_file = os.path(stat_path, 
                   "".join(["stat_", fasta_prefix, types, "_SNP.csv"]))
        out_stat = open(stat_file, "w")
        SNP_detect(fasta, snp, out_table_prefix, qual, seq_path + fasta_prefix, 
                   read_depth, fraction, bam_number, stat_file)
        self.helper.move_all_content(out_table_prefix, stat_path, ".png") 

    def _run_program(self, program, samtools_path, bcftools_path, fasta_file, out_bcf, 
                     out_raw_prefix, out_table_prefix, quality, seq_path, prefix, 
                     depth, stat_path, bam_number, fraction, out_folder):
        if "1" in program:
            print("Running SNP calling with BAQ...")
            call([samtools_path, "mpileup", 
                  "-t", "DP", "-ugf", fasta_file,
                  os.path.join(out_folder, "whole_reads_sorted.bam"),
                  "--ignore-RG"], stdout=out_bcf)
            call([bcftools_path , "call", "tmp_bcf",
                  "-vmO", "v", "-o", out_raw_prefix + "_with_BAQ.vcf"])
            self.helper.check_make_folder(os.path.join(seq_path, "with_BAQ"), prefix)
            self._transcript_SNP(fasta_file, out_raw_prefix + "_with_BAQ.vcf",
                                 out_table_prefix + "_with_BAQ", quality,
                                 os.path.join(seq_path, "with_BAQ", prefix),
                                 prefix, depth, stat_path, "_with_BAQ",
                                 bam_number, fraction)
        if "2" in program:
            print("Running SNP calling without BAQ...")
            call([samtools_path, "mpileup", 
                  "-t", "DP", "-B", "-ugf", fasta_file,
                  os.path.join(out_folder, "whole_reads_sorted.bam"),
                  "--ignore-RG"], stdout=out_bcf)
            call([bcftools_path, "call", "tmp_bcf",
                  "-vmO", "v", "-o", out_raw_prefix + "_without_BAQ.vcf"])
            self.helper.check_make_folder(os.path.join(seq_path, "without_BAQ"), prefix)
            self._transcript_SNP(fasta_file, out_raw_prefix + "_without_BAQ.vcf",
                                 out_table_prefix + "_without_BAQ", quality,
                                 os.path.join(seq_path, "without_BAQ", prefix),
                                 prefix, depth, stat_path, "_without_BAQ",
                                 bam_number, fraction)
        if "3" in program:
            print("Running SNP calling extend BAQ...")
            call([samtools_path, "mpileup",
                  "-t", "DP", "-E", "-ugf", fasta_file,
                  os.path.join(out_folder, "whole_reads_sorted.bam"),
                  "--ignore-RG"], stdout=out_bcf)
            call([bcftools_path, "call", "tmp_bcf",
                  "-vmO", "v", "-o", out_raw_prefix + "_extend_BAQ.vcf"])
            self.helper.check_make_folder(os.path.join(seq_path, "extend_BAQ", prefix))
            self._transcript_SNP(fasta_file, out_raw_prefix + "_extend_BAQ.vcf",
                                 out_table_prefix + "_extend_BAQ", quality,
                                 os.path.join(seq_path, "extend_BAQ", prefix),
                                 prefix, depth, stat_path, "_extend_BAQ",
                                 bam_number, fraction)

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
        if "whole_reads.bam" in os.listdir(out_folder):
            self.helper.remove_all_content(out_folder, "whole_read", "file")
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
        command = (" ".join([samtools_path, "merge", 
                   os.path.join(out_folder, "whole_reads.bam"), " ".join(bams)]))
        os.system(command)
        print("Sort BAM file now ...")
        command = (" ".join([samtools_path, "sort", 
                   os.path.join(out_folder, "whole_reads.bam"),
                   os.path.join(out_folder, "whole_reads_sorted")]))
        os.system(command)
        return num_bam

    def _modify_header(self, fastas):
        for fasta in os.listdir(fastas):
            if fasta.endswith("fasta") or \
               fasta.endswith("fa") or \
               fasta.endswith("fna"):
                self.seq_editer.modify_header(os.path.join(fastas, fasta))

    def _get_genome_name(self, samtools_path, out_folder):
        command = " ".join([samtools_path, "view", "-H", 
                  os.path.join(out_folder, "whole_reads_sorted.bam")])
        os.system(command + "> header")
        fh = open("header", "r");
        seq_names = []
        for row in csv.reader(fh, delimiter="\t"):
            if row[0] == "@SQ":
                seq_names.append(row[1].split(":")[1])
        return seq_names

    def run_SNP_calling(self, samtools_path, bcftools_path, types, program, fastas, 
                        normal_bams, frag_bams, quality, depth, out_folder, fraction):
        if types == "reference":
            file_type = "compare_reference"
        else:
            file_type = "validate_target"
        seq_path = os.path.join(out_folder, file_type, "seqs")
        stat_path = os.path.join(out_folder, file_type, "statistics")
        self.multiparser._parser_fasta(fastas)
        fasta_path = os.path.join(fastas, "tmp")
        self._modify_header(fastas)
        bam_number = self._merge_bams(normal_bams, frag_bams, samtools_path, out_folder)
        seq_names = self._get_genome_name(samtools_path, out_folder)
        #### running SNP calling
#        if depth is False:
#            depth = -5
#        else:
#            depth = int(depth)
#        if ("1" not in program ) and \
#           ("2" not in program) and \
#           ("3" not in program):
#            print("Error:Please assign a correct BAQ type: '1' means 'with_BAQ', '2' means 'with_BAQ' or '3' means 'extend_BAQ'.")
#            sys.exit()
#        else:
#            for fasta in os.listdir(fasta_path):
#                if (fasta.split(".f")[0] in seq_names):
#                    fasta_datas = self._detect_fasta(fasta)
#                    detect = fasta_datas[0]
#                    prefix = fasta_datas[1]
#                    if detect:
#                        detect = False
#                        print("Computing " + fasta + " now ...")
#                        out_bcf = open("tmp_bcf", "w")
#                        self.helper.check_make_folder(os.path.join(out_folder, file_type, "SNP_table"), prefix)
#                        out_raw_prefix = os.path.join(out_folder, file_type, "SNP_raw_outputs", prefix, prefix)
#                        out_table_prefix = os.path.join(out_folder, file_type, "SNP_table", prefix, prefix)
#                        fasta_file = os.path.join(fasta_path, fasta)
#                        self._run_program(program, samtools_path, bcftools_path, fasta_file, out_bcf, 
#                                          out_raw_prefix, out_table_prefix, quality, seq_path, prefix, 
#                                          depth, stat_path, bam_number, fraction, out_folder)
#                        os.remove("tmp_bcf")
#        self.helper.remove_tmp(fastas)
#        os.remove(os.path.join(out_folder, "whole_reads.bam"))
#        os.remove(os.path.join(out_folder, "whole_reads_sorted.bam"))
#        os.remove("header")
