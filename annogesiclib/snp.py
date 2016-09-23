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
    '''detection of SNP'''

    def __init__(self, args_snp):
        self.multiparser = Multiparser()
        self.seq_editer = SeqEditer()
        self.helper = Helper()
        if args_snp.types == "reference":
            file_type = "compare_reference"
        else:
            file_type = "validate_target"
        self.seq_path = os.path.join(args_snp.out_folder, file_type, "seqs")
        self.stat_path = os.path.join(args_snp.out_folder, file_type,
                                      "statistics")
        self.fasta_path = os.path.join(args_snp.fastas, "tmp")
        self.outputs = {"table": os.path.join(
                        args_snp.out_folder, file_type, "SNP_table"),
                        "raw": os.path.join(
                        args_snp.out_folder, file_type, "SNP_raw_outputs"),
                        "tmp": os.path.join(args_snp.out_folder, "tmp_bcf"),
                        "depth": os.path.join(args_snp.out_folder, "tmp_depth")}
        if "whole_reads.bam" in os.listdir(args_snp.out_folder):
            self.helper.remove_all_content(args_snp.out_folder,
                                           "whole_read", "file")
        self.bams = {"whole": os.path.join(args_snp.out_folder,
                                           "whole_reads.bam"),
                     "sort": os.path.join(args_snp.out_folder,
                                          "whole_reads_sorted.bam"),
                     "bams": []}
        self.header = os.path.join(args_snp.out_folder, "header")
        self.baqs = {"with": "with_BAQ", "without": "without_BAQ",
                     "extend": "extend_BAQ"}

    def _import_bam(self, bam_folder, bams):
        num_bam = 0
        for bam in os.listdir(bam_folder):
            if bam.endswith(".bam"):
                num_bam += 1
                bams.append(os.path.join(bam_folder, bam))
        return num_bam

    def _transcript_snp(self, fasta, snp, out_table_prefix, type_,
                        prefix, bam_number, table_path, args_snp):
        seq_path = os.path.join(self.seq_path, self.baqs[type_], prefix)
        stat_prefix = os.path.join(self.stat_path, "_".join([
            "stat", "_".join([prefix, self.baqs[type_]]), "SNP"]))
        snp_detect(fasta, snp, self.outputs["depth"], out_table_prefix,
                   os.path.join(seq_path, prefix), bam_number,
                   stat_prefix, args_snp)
        self.helper.move_all_content(table_path, self.stat_path, [".png"])

    def _get_para(self, args_snp):
        bams = self.bams["sort"]
        if args_snp.caller == "c":
            bcf_para = "-vcO"
        else:
            bcf_para = "-vmO"
        return bams, bcf_para

    def _run_tools(self, fasta_file, out_raw_prefix, type_, args_snp):
        bams, bcf_para = self._get_para(args_snp)
        if type_ == "with":
            command = [args_snp.samtools_path, "mpileup", "-t", "DP"]
        elif type_ == "without":
            command = [args_snp.samtools_path, "mpileup", "-t", "DP", "-B"]
        elif type_ == "extend":
            command = [args_snp.samtools_path, "mpileup", "-t", "DP", "-E"]
        if args_snp.rg:
            command = command + ["-ugf", fasta_file, bams]
        else:
            command = command + ["--ignore-RG", "-ugf", fasta_file, bams]
        os.system(" ".join(command) + ">" + self.outputs["tmp"])
        out_vcf = "_".join([out_raw_prefix, self.baqs[type_] + ".vcf"])
        if args_snp.chrom == "1":
            call([args_snp.bcftools_path, "call", "--ploidy", args_snp.chrom,
                  self.outputs["tmp"], bcf_para, "v", "-o", out_vcf])
        elif args_snp.chrom == "2":
            call([args_snp.bcftools_path, "call",
                  self.outputs["tmp"], bcf_para, "v", "-o", out_vcf])
        return out_vcf

    def _run_sub(self, args_snp, fasta_file, type_, file_prefixs, prefix,
                 table_path, bam_number):
        out_vcf = self._run_tools(fasta_file, file_prefixs["raw_prefix"],
                                  type_, args_snp)
        self.helper.check_make_folder(
             os.path.join(self.seq_path, self.baqs[type_], prefix))
        self._transcript_snp(
            fasta_file, out_vcf,
            "_".join([file_prefixs["table_prefix"], self.baqs[type_]]),
            type_, prefix, bam_number, table_path, args_snp)

    def _run_program(self, fasta_file, file_prefixs, prefix, bam_number,
                     table_path, args_snp):
        for index in args_snp.program:
            if index == "1":
                type_ = "with"
                print("Running SNP calling with BAQ...")
            elif index == "2":
                type_ = "without"
                print("Running SNP calling without BAQ...")
            elif index == "3":
                print("Running SNP calling extend BAQ...")
                type_ = "extend"
            else:
                print("Error: No correct program, please assign 1, 2, 3")
                sys.exit()
            self._run_sub(args_snp, fasta_file, type_, file_prefixs, prefix,
                          table_path, bam_number)

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

    def _run_bam(self, samtools_path, sub_command, bam_file):
        if sub_command == "merge":
            command = (" ".join([samtools_path, sub_command,
                       self.bams["whole"], bam_file]))
        elif sub_command == "sort":
            command = (" ".join([samtools_path, sub_command,
                                 "-o", bam_file, self.bams["whole"]]))
        os.system(command)
        self.bams["bams"].append(bam_file.replace(".bam", "_sort.bam"))

    def _merge_bams(self, args_snp):
        bams = []
        num_normal = 0
        num_frag = 0
        if (args_snp.frag_bams is None) and (args_snp.normal_bams is None):
            print("Error: There is no BAMs folders!!")
            sys.exit()
        else:
            if args_snp.normal_bams is not None:
                num_normal = self._import_bam(args_snp.normal_bams, bams)
            if args_snp.frag_bams is not None:
                num_frag = self._import_bam(args_snp.frag_bams, bams)
        num_bam = num_normal + num_frag
        if num_bam <= 1:
            shutil.copyfile(bams[0], self.bams["whole"])
            print("Sort BAM file now ...")
            self._run_bam(args_snp.samtools_path, "sort",
                          self.bams["sort"])
        else:
            print("Merge BAM files now ...")
            self._run_bam(args_snp.samtools_path, "merge",
                          " ".join(bams))
            print("Sort BAM file now ...")
            self._run_bam(args_snp.samtools_path, "sort",
                          self.bams["sort"])
        out_depth = open(self.outputs["depth"], "w")
        call([args_snp.samtools_path, "index",  self.bams["sort"]])
        call([args_snp.samtools_path, "depth",  self.bams["sort"]],
             stdout=out_depth)
        return num_bam

    def _modify_header(self, fastas):
        for fasta in os.listdir(fastas):
            if fasta.endswith("fasta") or \
               fasta.endswith("fa") or \
               fasta.endswith("fna"):
                self.seq_editer.modify_header(os.path.join(fastas, fasta))

    def _get_header(self, samtools_path, bam, seq_names):
        command = " ".join([samtools_path, "view", "-H", bam])
        os.system(">".join([command, self.header]))
        fh = open(self.header, "r")
        for row in csv.reader(fh, delimiter="\t"):
            if row[0] == "@SQ":
                seq_names.append(row[1].split(":")[1])
        fh.close()

    def _get_genome_name(self, args_snp):
        seq_names = []
        self._get_header(args_snp.samtools_path, self.bams["sort"],
                         seq_names)
        return seq_names

    def _remove_bams(self):
        if os.path.exists(self.bams["whole"]):
            os.remove(self.bams["whole"])
        if os.path.exists(self.bams["sort"]):
            os.remove(self.bams["sort"])
        if os.path.exists(self.header):
            os.remove(self.header)
        os.remove(self.outputs["depth"])

    def run_snp_calling(self, args_snp):
        self.multiparser.parser_fasta(args_snp.fastas)
        self._modify_header(args_snp.fastas)
        bam_number = self._merge_bams(args_snp)
        seq_names = self._get_genome_name(args_snp)
        if ("1" not in args_snp.program) and (
                "2" not in args_snp.program) and (
                "3" not in args_snp.program):
            print("Error:Please assign a correct BAQ type: "
                  "'1' means 'with_BAQ', '2' means 'with_BAQ' or "
                  "'3' means 'extend_BAQ'.")
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
                        file_prefixs = {"raw_prefix": os.path.join(
                                        self.outputs["raw"], prefix, prefix),
                                        "table_prefix": os.path.join(
                                        self.outputs["table"], prefix, prefix)}
                        fasta_file = os.path.join(self.fasta_path, fasta)
                        table_path = os.path.join(self.outputs["table"],
                                                  prefix)
                        self._run_program(fasta_file, file_prefixs, prefix,
                                          bam_number, table_path, args_snp)
                        os.remove(self.outputs["tmp"])
        self.helper.remove_tmp(args_snp.fastas)
        self._remove_bams()
