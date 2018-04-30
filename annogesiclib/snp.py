import os
import sys
import csv
import shutil
from glob import glob
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
        if args_snp.types == "related_genome":
            file_type = "compare_related_and_reference_genomes"
        else:
            file_type = "mutations_of_reference_genomes"
        self.seq_path = os.path.join(args_snp.out_folder, file_type, "seqs")
        self.stat_path = os.path.join(args_snp.out_folder, file_type,
                                      "statistics")
        self.fig_path = os.path.join(self.stat_path, "figs")
        self.helper.check_make_folder(self.fig_path)
        self.outputs = {"table": os.path.join(
                        args_snp.out_folder, file_type, "SNP_tables"),
                        "raw": os.path.join(
                        args_snp.out_folder, file_type, "SNP_raw_outputs"),
                        "tmp": os.path.join(args_snp.out_folder, "tmp_bcf"),
                        "depth": os.path.join(args_snp.out_folder, "tmp_depth")}
        self.bams = {"whole": os.path.join(args_snp.out_folder,
                                           "whole_reads.bam"),
                     "sort": os.path.join(args_snp.out_folder,
                                          "whole_reads_sorted.bam"),
                     "bams": []}
        self.header = os.path.join(args_snp.out_folder, "header")
        self.baqs = {"with": "with_BAQ", "without": "without_BAQ",
                     "extend": "extend_BAQ"}

    def _transcript_snp(self, fasta, out_table_prefix, type_,
                        prefix, bam_datas, table_path, args_snp):
        seq_path = os.path.join(self.seq_path, self.baqs[type_], prefix)
        for bam in bam_datas:
            stat_prefix = os.path.join(self.stat_path, "_".join([
                "stat", "_".join([prefix, self.baqs[type_], bam["sample"]]),
                "SNP"]))
            snp_file = os.path.join(self.outputs["raw"], prefix, "_".join(
                [prefix, self.baqs[type_], bam["sample"] + ".vcf"]))
            snp_detect(
                fasta, snp_file, self.outputs["depth"] + bam["sample"],
                "_".join([out_table_prefix, bam["sample"]]),
                os.path.join(seq_path, "_".join([prefix, bam["sample"]])),
                bam["bam_number"], stat_prefix, args_snp, bam["rep"])
            self.helper.move_all_content(table_path, self.fig_path, [".png"])

    def _get_para(self, args_snp):
        if args_snp.caller == "c":
            bcf_para = "-vcO"
        else:
            bcf_para = "-vmO"
        return bcf_para

    def _run_tools(self, fasta_file, type_, args_snp, bam_datas, log):
        bcf_para = self._get_para(args_snp)
        for bam in bam_datas:
            bam_file = os.path.join(args_snp.out_folder,
                                    bam["sample"] + ".bam")
            if type_ == "with":
                command = [args_snp.samtools_path, "mpileup", "-t", "DP"]
            elif type_ == "without":
                command = [args_snp.samtools_path, "mpileup", "-t", "DP", "-B"]
            elif type_ == "extend":
                command = [args_snp.samtools_path, "mpileup", "-t", "DP", "-E"]
            if args_snp.rg:
                command = command + ["-ugf", fasta_file, bam_file]
            else:
                command = command + ["--ignore-RG", "-ugf", fasta_file, bam_file]
            log.write(" ".join(command) + ">" + self.outputs["tmp"] + "\n")
            os.system(" ".join(command) + ">" + self.outputs["tmp"])
            bam["vcf"] = os.path.join(self.outputs["raw"], "_".join(
                [self.baqs[type_], bam["sample"] + ".vcf"]))
            if args_snp.chrom == "1":
                log.write(" ".join([
                      args_snp.bcftools_path, "call", "--ploidy", args_snp.chrom,
                      self.outputs["tmp"], bcf_para, "v", "-o", bam["vcf"]]) + "\n")
                call([args_snp.bcftools_path, "call", "--ploidy", args_snp.chrom,
                      self.outputs["tmp"], bcf_para, "v", "-o", bam["vcf"]])
            elif args_snp.chrom == "2":
                log.write(" ".join([args_snp.bcftools_path, "call",
                      self.outputs["tmp"], bcf_para, "v", "-o", bam["vcf"]]) + "\n")
                call([args_snp.bcftools_path, "call",
                      self.outputs["tmp"], bcf_para, "v", "-o", bam["vcf"]])
        log.write("Done!\n")
        log.write("The following files are generated:\n")
        for file_ in os.listdir(self.outputs["raw"]):
            log.write("\t" + os.path.join(self.outputs["raw"], file_) + "\n")

    def _parse_vcf_by_fa(self, args_snp, type_, num_prog, log):
        seq_names = []
        fa_prefixs = []
        log.write("Parsing Vcf files by comparing fasta information.\n")
        for fa in os.listdir(args_snp.fastas):
            if (fa != "all.fa") and (not fa.endswith(".fai")):
                with open(os.path.join(args_snp.fastas, fa)) as fh:
                    for line in fh:
                        line = line.strip()
                        if line.startswith(">"):
                            seq_names.append(line[1:])
                fa_prefix = ".".join(fa.split(".")[:-1])
                fa_prefixs.append(fa_prefix)
                vcf_folder = os.path.join(
                    self.outputs["raw"], fa_prefix)
                if num_prog == 0:
                    self.helper.check_make_folder(vcf_folder)
                    self.helper.check_make_folder(os.path.join(
                        self.outputs["table"], fa_prefix))
                self.helper.check_make_folder(
                    os.path.join(self.seq_path, self.baqs[type_], fa_prefix))
                for vcf in os.listdir(self.outputs["raw"]):
                    if vcf.endswith(".vcf"):
                        out = open(os.path.join(vcf_folder, "_".join(
                            [fa_prefix, vcf])), "w")
                        with open(os.path.join(self.outputs["raw"],
                                  vcf)) as vh:
                            for line in vh:
                                line = line.strip()
                                if line.startswith("#"):
                                    out.write(line + "\n")
                                else:
                                    if line.split("\t")[0] in seq_names:
                                        out.write(line + "\n")
                        out.close()
                        log.write("\t" + os.path.join(vcf_folder, "_".join(
                            [fa_prefix, vcf])) + " is generated.\n")
        for vcf in os.listdir(self.outputs["raw"]):
            if vcf.endswith(".vcf"):
                os.remove(os.path.join(self.outputs["raw"], vcf))
        return fa_prefixs

    def _run_sub(self, args_snp, all_fasta, type_, bam_datas, num_prog, log):
        self._run_tools(all_fasta, type_, args_snp, bam_datas, log)
        fa_prefixs = self._parse_vcf_by_fa(args_snp, type_, num_prog, log)
        log.write("Running transcript_SNP.py to do statistics, filter SNPs, "
                  "and generate potential sequences.\n")
        log.write("The following files are generated:\n")
        for fa_prefix in fa_prefixs:
            for fasta in os.listdir(args_snp.fastas):
                if fa_prefix in fasta:
                    fasta_file = os.path.join(args_snp.fastas, fasta)
            table_path = os.path.join(self.outputs["table"], fa_prefix)
            table_prefix = os.path.join(table_path, "_".join(
                [fa_prefix, self.baqs[type_]]))
            self._transcript_snp(
                fasta_file, table_prefix,
                type_, fa_prefix, bam_datas, table_path, args_snp)
            seq_path = os.path.join(self.seq_path, self.baqs[type_], fa_prefix)
            for folder in (table_path, self.stat_path, seq_path, self.fig_path):
                for file_ in os.listdir(folder):
                    if os.path.isfile(os.path.join(folder, file_)):
                        log.write("\t" + os.path.join(folder, file_) + "\n")

    def _run_program(self, all_fasta, bam_datas, args_snp, log):
        num_prog = 0
        log.write("Running Samtools to mpileup, and using Bcftools to "
                  "call snp.\n")
        log.write("Please make sure the version of Samtools and Bcftools "
                  "are both at least 1.3.1.\n")
        for index in args_snp.program:
            if index == "with_BAQ":
                type_ = "with"
                print("Running SNP calling with BAQ")
                log.write("Running SNP calling with BAQ.\n")
            elif index == "without_BAQ":
                type_ = "without"
                print("Running SNP calling without BAQ")
                log.write("Running SNP calling without BAQ.\n")
            elif index == "extend_BAQ":
                print("Running SNP calling extend BAQ")
                log.write("Running SNP calling extend BAQ.\n")
                type_ = "extend"
            else:
                print("Error: No correct program, please assign "
                      "\"with_BAQ\", \"without_BAQ\", \"extend_BAQ\"!")
                log.write("No valid program can be found, please assign"
                          "\"with_BAQ\", \"without_BAQ\", \"extend_BAQ\".\n")
                sys.exit()
            self._run_sub(args_snp, all_fasta, type_, bam_datas, num_prog, log)
            num_prog += 1

    def _run_bam(self, samtools_path, sub_command, bam_file, type_file, log):
        if sub_command == "merge":
            command = (" ".join([samtools_path, sub_command,
                       self.bams["whole"], bam_file]))
        elif sub_command == "sort":
            if type_file == "all":
                command = (" ".join([samtools_path, sub_command,
                                     "-o", bam_file, self.bams["whole"]]))
            else:
                command = (" ".join([samtools_path, sub_command,
                                     "-o",
                                     bam_file, type_file]))
        log.write(command + "\n")
        os.system(command)

    def _merge_bams(self, args_snp, bam_datas, log):
        bams = []
        num_normal = 0
        num_frag = 0
        log.write("Using Samtools to merge and sort BAM files.\n")
        log.write("Please make sure the version of Samtools is at least 1.3.1.\n")
        for bam in bam_datas:
            bam["bam_number"] = 0
            out_bam = os.path.join(args_snp.out_folder, bam["sample"] + ".bam")
            if len(bam["bams"]) == 1:
                print("Sorting BAM files of " + bam["sample"])
                self._run_bam(
                    args_snp.samtools_path, "sort",
                    out_bam, bam["bams"][0], log)
                bam["bam_number"] = 1
            else:
                print("Merging BAM files of " + bam["sample"])
                self._run_bam(args_snp.samtools_path, "merge",
                              " ".join(bam["bams"]), "all", log)
                print("Sorting BAM files of " + bam["sample"])
                self._run_bam(
                    args_snp.samtools_path, "sort",
                    out_bam, "all", log)
                bam["bam_number"] += 1
            if os.path.exists(self.bams["whole"]):
                os.remove(self.bams["whole"])
            out_depth = open(self.outputs["depth"] + bam["sample"], "w")
            log.write(" ".join([args_snp.samtools_path, "index",  out_bam]) + "\n")
            call([args_snp.samtools_path, "index",  out_bam])
            log.write(" ".join([args_snp.samtools_path, "depth",  out_bam]) + "\n")
            call([args_snp.samtools_path, "depth",  out_bam],
                 stdout=out_depth)
            out_depth.close()
        log.write("Done!\n")
        log.write("The following files are generated:\n")
        log.write("\t" + self.bams["whole"] + " is temporary generated "
                  "(be deleted afterward).\n")
        for file_ in os.listdir(args_snp.out_folder):
            if os.path.isfile(os.path.join(args_snp.out_folder, file_)):
                log.write("\t" + os.path.join(args_snp.out_folder, file_) + "\n")
        

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
                if row[1].split(":")[1] not in seq_names:
                    seq_names.append(row[1].split(":")[1])
        fh.close()

    def _get_genome_name(self, args_snp, bam_datas):
        seq_names = []
        for bam in bam_datas:
            bam_file = os.path.join(args_snp.out_folder,
                                    bam["sample"] + ".bam")
            self._get_header(args_snp.samtools_path,
                             bam_file, seq_names)
        return seq_names

    def _remove_bams(self, bam_datas, args_snp):
        for bam in bam_datas:
            bam_file = os.path.join(args_snp.out_folder,
                                    bam["sample"] + ".bam")
            if os.path.exists(bam_file):
                os.remove(bam_file)
            if os.path.exists(bam_file + ".bai"):
                os.remove(bam_file + ".bai")
            if os.path.exists(self.header):
                os.remove(self.header)
            os.remove(self.outputs["depth"] + bam["sample"])

    def _extract_bams(self, bams, log):
        bam_datas = []
        for bam in bams:
            datas = bam.split(":")
            if len(datas) != 2:
                log.write("the format of --bam_files is wrong!\n")
                print("Error: the format of --bam_files is wrong!")
                sys.exit()
            for file_ in datas[-1].split(","):
                if not os.path.exists(file_):
                    print("Error: there are some Bam files "
                          "which do not exist!")
                    log.write(file_ + " is not found.\n")
                    sys.exit()
            bam_datas.append({"sample": datas[0],
                              "rep": len(datas[-1].split(",")),
                              "bams": datas[-1].split(",")})
        return bam_datas

    def _merge_fasta(self, fastas, log):
        all_fasta = os.path.join(fastas, "all.fa")
        names = []
        out = open(all_fasta, "w")
        print_ = False
        for fasta in os.listdir(fastas):
            if (fasta.endswith(".fa")) or (
                    fasta.endswith(".fasta")) or (
                    fasta.endswith(".fna")):
                with open(os.path.join(fastas, fasta)) as fh:
                    for line in fh:
                        line = line.strip()
                        if line.startswith(">"):
                            if line not in names:
                                print_ = True
                                names.append(line)
                            else:
                                print_ = False
                        if print_:
                            out.write(line + "\n")
                log.write(os.path.join(fastas, fasta) + " is loaded.\n")
        out.close()
        return all_fasta

    def run_snp_calling(self, args_snp, log):
        self._modify_header(args_snp.fastas)
        all_fasta = self._merge_fasta(args_snp.fastas, log)
        bam_datas = self._extract_bams(args_snp.bams, log)
        self._merge_bams(args_snp, bam_datas, log)
        if ("with_BAQ" not in args_snp.program) and (
                "without_BAQ" not in args_snp.program) and (
                "extend_BAQ" not in args_snp.program):
            print("Error: Please assign a correct programs: "
                  "\"with_BAQ\", \"without_BAQ\", \"extend_BAQ\".")
            sys.exit()
        else:
            print("Detecting mutations now")
            self._run_program(all_fasta, bam_datas, args_snp, log)
            os.remove(self.outputs["tmp"])
            os.remove(all_fasta)
            os.remove(all_fasta + ".fai")
        self.helper.remove_tmp_dir(args_snp.fastas)
        self._remove_bams(bam_datas, args_snp)
        log.write("Remove all the temporary files.\n")
