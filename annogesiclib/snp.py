import os
import sys
import csv
import shutil
from subprocess import call
from annogesiclib.multiparser import Multiparser
from annogesiclib.seq_editer import SeqEditer
from annogesiclib.transcript_SNP import snp_detect
from annogesiclib.helper import Helper

class SNPCalling(object):

    def __init__(self, types, out_folder, fastas):
        self.multiparser = Multiparser()
        self.seq_editer = SeqEditer()
        self.helper = Helper()
        if types == "reference":
            file_type = "compare_reference"
        else:
            file_type = "validate_target"
        self.seq_path = os.path.join(out_folder, file_type, "seqs")
        self.stat_path = os.path.join(out_folder, file_type, "statistics")
        self.fasta_path = os.path.join(fastas, "tmp")
        self.outputs = {"table": os.path.join(
                                 out_folder, file_type, "SNP_table"),
                        "raw": os.path.join(
                               out_folder, file_type, "SNP_raw_outputs"),
                        "tmp": os.path.join(out_folder, "tmp_bcf")}
        if "whole_reads.bam" in os.listdir(out_folder):
            self.helper.remove_all_content(out_folder, "whole_read", "file")
        self.bams = {"whole": os.path.join(out_folder, "whole_reads.bam"),
                     "sort": os.path.join(out_folder, "whole_reads_sorted.bam")}
        self.header = os.path.join(out_folder, "header")
        self.baqs = {"with": "with_BAQ", "without": "without_BAQ",
                     "extend": "extend_BAQ"}

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
        stat_file = os.path.join(stat_path, "_".join(["stat",
                    "_".join([fasta_prefix, types]), "SNP.csv"]))
        snp_detect(fasta, snp, out_table_prefix, qual,
                   os.path.join(seq_path, fasta_prefix),
                   read_depth, fraction, bam_number, stat_file)
        self.helper.move_all_content(table_path, stat_path, [".png"])

    def _run_program(self, program, samtools_path, bcftools_path, fasta_file,
                     out_raw_prefix, out_table_prefix, quality, seq_path,
                     prefix, depth, stat_path, bam_number, fraction,
                     out_folder, table_path):
        if "1" in program:
            print("Running SNP calling with BAQ...")
            out_bcf = open(self.outputs["tmp"], "w")
            call([samtools_path, "mpileup",
                  "-t", "DP", "-ugf", fasta_file, self.bams["sort"],
                  "--ignore-RG"], stdout=out_bcf)
            out_vcf = "_".join([out_raw_prefix, self.baqs["with"] + ".vcf"])
            call([bcftools_path, "call", self.outputs["tmp"],
                  "-vmO", "v", "-o", out_vcf])
            self.helper.check_make_folder(
                 os.path.join(seq_path, self.baqs["with"], prefix))
            self._transcript_snp(fasta_file, out_vcf,
                 "_".join([out_table_prefix, self.baqs["with"]]), quality,
                 os.path.join(seq_path, self.baqs["with"], prefix),
                 prefix, depth, stat_path, self.baqs["with"],
                 bam_number, fraction, table_path)
            out_bcf.close()
        if "2" in program:
            print("Running SNP calling without BAQ...")
            out_bcf = open(self.outputs["tmp"], "w")
            call([samtools_path, "mpileup",
                  "-t", "DP", "-B", "-ugf", fasta_file,
                  self.bams["sort"], "--ignore-RG"],
                  stdout=out_bcf)
            out_vcf = "_".join([out_raw_prefix, self.baqs["without"] + ".vcf"])
            call([bcftools_path, "call", self.outputs["tmp"],
                  "-vmO", "v", "-o", out_vcf])
            self.helper.check_make_folder(
                 os.path.join(seq_path, self.baqs["without"], prefix))
            self._transcript_snp(fasta_file, out_vcf,
                 "_".join([out_table_prefix, self.baqs["without"]]), quality,
                 os.path.join(seq_path, self.baqs["without"], prefix),
                 prefix, depth, stat_path, self.baqs["without"],
                 bam_number, fraction, table_path)
            out_bcf.close()
        if "3" in program:
            print("Running SNP calling extend BAQ...")
            out_bcf = open(self.outputs["tmp"], "w")
            call([samtools_path, "mpileup",
                  "-t", "DP", "-E", "-ugf", fasta_file,
                  self.bams["sort"], "--ignore-RG"], stdout=out_bcf)
            out_vcf = "_".join([out_raw_prefix, self.baqs["extend"] + ".vcf"])
            call([bcftools_path, "call", self.outputs["tmp"],
                  "-vmO", "v", "-o", out_vcf])
            self.helper.check_make_folder(
                 os.path.join(seq_path, self.baqs["extend"], prefix))
            self._transcript_snp(fasta_file, out_vcf,
                 "_".join([out_table_prefix, self.baqs["extend"]]), quality,
                 os.path.join(seq_path, self.baqs["extend"], prefix),
                 prefix, depth, stat_path, self.baqs["extend"],
                 bam_number, fraction, table_path)
            out_bcf.close()

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
        if (frag_bams is None) and (normal_bams is None):
            print("Error: There is no BAMs folders!!")
            sys.exit()
        else:
            if normal_bams is not None:
                num_normal = self._import_bam(normal_bams, bams)
            if frag_bams is not None:
                num_frag = self._import_bam(frag_bams, bams)
        num_bam = num_normal + num_frag
        if num_bam <= 1:
            shutil.copyfile(bams[0], self.bams["whole"])
            print("Sort BAM file now ...")
            command = (" ".join([samtools_path, "sort",
                       self.bams["whole"], self.bams["sort"].replace(".bam", "")]))
            os.system(command)
        else:
            print("Merge BAM files now ...")
            command = (" ".join([samtools_path, "merge",
                       self.bams["whole"], " ".join(bams)]))
            os.system(command)
            print("Sort BAM file now ...")
            command = (" ".join([samtools_path, "sort",
                       self.bams["whole"], self.bams["sort"].replace(".bam", "")]))
            os.system(command)
        return num_bam

    def _modify_header(self, fastas):
        for fasta in os.listdir(fastas):
            if fasta.endswith("fasta") or \
               fasta.endswith("fa") or \
               fasta.endswith("fna"):
                self.seq_editer.modify_header(os.path.join(fastas, fasta))

    def _get_genome_name(self, samtools_path, out_folder):
        command = " ".join([samtools_path, "view", "-H", self.bams["sort"]])
        os.system(">".join([command, self.header]))
        fh = open(self.header, "r");
        seq_names = []
        for row in csv.reader(fh, delimiter="\t"):
            if row[0] == "@SQ":
                seq_names.append(row[1].split(":")[1])
        return seq_names

    def run_snp_calling(self, samtools_path, bcftools_path, types, program,
                        fastas, normal_bams, frag_bams, quality, depth,
                        out_folder, fraction):
        if types == "reference":
            file_type = "compare_reference"
        else:
            file_type = "validate_target"
        self.multiparser.parser_fasta(fastas)
        self._modify_header(fastas)
        bam_number = self._merge_bams(normal_bams, frag_bams,
                                      samtools_path, out_folder)
        seq_names = self._get_genome_name(samtools_path, out_folder)
        #### running SNP calling
        if ("1" not in program) and (
            "2" not in program) and (
            "3" not in program):
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
                        self.helper.check_make_folder(
                             os.path.join(self.outputs["table"], prefix))
                        self.helper.check_make_folder(
                             os.path.join(self.outputs["raw"], prefix))
                        out_raw_prefix = os.path.join(
                                         self.outputs["raw"], prefix, prefix)
                        out_table_prefix = os.path.join(
                                         self.outputs["table"], prefix, prefix)
                        fasta_file = os.path.join(self.fasta_path, fasta)
                        table_path = os.path.join(self.outputs["table"], prefix)
                        self._run_program(program, samtools_path, bcftools_path,
                                          fasta_file, out_raw_prefix,
                                          out_table_prefix, quality,
                                          self.seq_path, prefix, depth,
                                          self.stat_path, bam_number, fraction,
                                          out_folder, table_path)
                        os.remove(self.outputs["tmp"])
        self.helper.remove_tmp(fastas)
        os.remove(self.bams["whole"])
        os.remove(self.bams["sort"])
        os.remove(self.header)
