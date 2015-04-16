#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
import time
import shutil
from subprocess import call, Popen
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.converter import Converter
from annogesiclib.circRNA import detect_circrna


class CircRNA_detection(object):

    def __init__(self):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.converter = Converter()

    def _remove_tmp(self, folder):
        if folder is not False:
            call(["rm", "-rf", folder + "/tmp"])
            os.system("rm -rf " + folder + "/*_folder")

    def _wait_process(self, processes):
        for p in processes:
            p.wait()
        if p.stdout:
            p.stdout.close()
        if p.stdin:
            p.stdin.close()
        if p.stderr:
            p.stderr.close()
        try:
            p.kill()
        except OSError:
            pass
        time.sleep(5)

    def _deal_zip_file(self, read_folder):
        for read in os.listdir(read_folder):
            if read.endswith(".bz2"):
                mod_read = read.replace(".bz2", "")
                if (".fa" not in mod_read) and \
                   (".fasta" not in mod_read) and \
                   (".fna" not in mod_read):
                    mod_read = mod_read + ".fa"
                read_out = open(os.path.join(read_folder, mod_read), "w")
                tmp_reads.append(os.path.join(read_folder, mod_read))
                print(" ".join(["unzip", read]))
                call(["bzcat", os.path.join(read_folder, read)], stdout=read_out)

    def _align(self, fasta_path, segemehl_path, prefixs, 
               alignment_path, read_folder, cores):
        for fasta in os.listdir(fasta_path):
            index = fasta.replace(".fa", ".idx")
            call([os.path.join(segemehl_path, "segemehl.x"),
                  "-x", os.path.join(fasta_path, index),
                  "-d", os.path.join(fasta_path, fasta)])
            processes = []
            num_process = 0
            fasta_prefix = fasta.replace(".fa", "")
            prefixs.append(fasta_prefix)
            self.helper.check_make_folder(alignment_path, fasta_prefix)
            for read in os.listdir(read_folder):
                num_process += 1
                if read.endswith(".fa") or \
                   read.endswith(".fna") or \
                   read.endswith("fasta"):
                    filename = read.split(".")
                    read_prefix = ".".join(filename[:-1])
                    sam_file = "_".join([read_prefix, fasta_prefix + ".sam"])
                    log_file = "_".join([read_prefix, fasta_prefix + ".log"])
                    print("mapping {0}".format(sam_file))
                    out = open(os.path.join(alignment_path, fasta_prefix, sam_file), "w")
                    log = open(os.path.join(alignment_path, fasta_prefix, log_file), "w")
                    p = Popen([os.path.join(segemehl_path, "segemehl.x"),
                               "-i", os.path.join(fasta_path, index),
                               "-d", os.path.join(fasta_path, fasta),
                               "-q", os.path.join(read_folder, read), "-S"], 
                               stdout=out, stderr=log)
                    processes.append(p)
                    if num_process == cores:
                        self._wait_process(processes)
                        num_process = 0
            self._wait_process(processes)

    def _convert_sam2bam(self, sub_alignment_path, samtools_path, convert_ones, bam_files):
        for sam in os.listdir(sub_alignment_path):
            if sam.endswith(".sam"):
                bam_file = sam.replace(".sam", ".bam")
                print("Convert {0} to {1}".format(sam, bam_file))
                out_bam = open(os.path.join(sub_alignment_path, bam_file), "w")
                call([samtools_path, "view", "-bS", 
                      os.path.join(sub_alignment_path, sam)], stdout=out_bam)
                bam_files.append(os.path.join(sub_alignment_path, bam_file))
                convert_ones.append(os.path.join(sub_alignment_path, bam_file))
            elif sam.endswith(".bam"):
                bam_files.append(os.path.join(sub_alignment_path, sam))

    def _merge_sort_aligment_file(self, bam_files, samtools_path, sub_alignment_path, 
                                  convert_ones, tmp_reads):
        print("Merge all bam files....")
        file_line = " ".join(bam_files)
        os.system(" ".join([samtools_path, "merge",
                            os.path.join(sub_alignment_path, "whole_reads.bam"),
                            file_line]))
        print("Sort bam files....")
        call([samtools_path, "sort", 
              os.path.join(sub_alignment_path, "whole_reads.bam"),
              os.path.join(sub_alignment_path, "whole_reads_sorted")])
        os.remove(os.path.join(sub_alignment_path, "whole_reads.bam"))
        print("Convert whole reads bam file to sam file....")
        call([samtools_path, "view", "-h", "-o",
              os.path.join(sub_alignment_path, "whole_reads_sorted.sam"),
              os.path.join(sub_alignment_path, "whole_reads_sorted.bam")])
        for bam in convert_ones:
            os.remove(bam)
        if len(tmp_reads) != 0:
            for read in tmp_reads:
                os.remove(read)

    def _run_testrealign(self, splice_path, prefix, segemehl_path, fasta_path, sub_alignment_path):
        self.helper.check_make_folder(splice_path, prefix)
        sub_splice_path = os.path.join(splice_path, prefix)
        err_log = os.path.join(sub_splice_path, prefix + ".log")
        print("Running testrealign.x for {0}".format(prefix))
        command = " ".join([os.path.join(segemehl_path, "testrealign.x"),
                            "-d", os.path.join(fasta_path, prefix + ".fa"),
                            "-q", os.path.join(sub_alignment_path, "whole_reads_sorted.sam"),
                            "-n"])
        os.system(command + " 2>" + err_log)
        self.helper.move_all_content(os.getcwd(), sub_splice_path, ".bed")
        self.helper.remove_all_content(sub_alignment_path, "whole_reads_sorted", "file")

    def _merge_bed(self, fastas, splice_path):
        tmp_prefixs = []
        for fasta in os.listdir(fastas):
            headers = []
            if fasta.endswith(".fa") or \
               fasta.endswith(".fna") or \
               fasta.endswith(".fasta"):
                with open(os.path.join(fastas, fasta), "r") as f_h:
                    for line in f_h:
                        line = line.strip()
                        if line.startswith(">"):
                            headers.append(line[1:])
                filename = fasta.split(".")
                fasta_prefix = filename[0]
                tmp_prefixs.append(fasta_prefix)
                self.helper.check_make_folder(os.getcwd(), fasta_prefix)
                for header in headers:
                    shutil.copyfile(os.path.join(splice_path, header, "splicesites.bed"),
                                    os.path.join(fasta_prefix, 
                                    "_".join(["splicesites", header + ".bed"])))
                    shutil.copyfile(os.path.join(splice_path, header, "transrealigned.bed"),
                                    os.path.join(fasta_prefix, 
                                    "_".join(["transrealigned", header + ".bed"])))
                out_splice = os.path.join(fasta_prefix, "splicesites_all.bed")
                out_trans = os.path.join(fasta_prefix, "transrealigned_all.bed")
                if len(headers) > 1:
                    for file_ in os.listdir(fasta_prefix):
                        if ("splicesites" in file_) and \
                           ("splicesites_all" not in file_):
                            self.helper.merge_file(fasta_prefix, file_, 
                                                   fasta_prefix, "splicesites_all.bed")
                        elif ("transrealigned" in file_) and \
                             ("transrealigned_all" not in file_):
                            self.helper.merge_file(fasta_prefix, file_, 
                                                   fasta_prefix, "transrealigned_all.bed")
                else:
                    os.rename(os.path.join(fasta_prefix, 
                              "_".join(["splicesites", headers[0] + ".bed"])), out_splice)
                    os.rename(os.path.join(fasta_prefix,
                              "_".join(["transrealigned", headers[0] + ".bed"])), out_trans)
        self.helper.remove_all_content(splice_path, None, "dir")
        return tmp_prefixs

    def _stat_and_gen_gff(self, tmp_prefixs, gff_folder, splice_path, candidate_path, gff_path,
                          stat_folder, support, convert, start_ratio, end_ratio):
        for prefix in tmp_prefixs:
            self.helper.check_make_folder(gff_folder, prefix)
            shutil.copytree(prefix, os.path.join(splice_path, prefix))
            self.helper.check_make_folder(candidate_path, prefix)
            print("comparing with annotation of {0}".format(prefix))
            if "splicesites_all.bed" in os.listdir(os.path.join(splice_path, prefix)):
                detect_circrna(os.path.join(splice_path, prefix, "splicesites_all.bed"), 
                               os.path.join(gff_path, prefix + ".gff"), 
                               os.path.join(candidate_path, prefix, 
                                            "_".join(["circRNA", prefix + ".txt"])), 
                               start_ratio, end_ratio, 
                               os.path.join(stat_folder, 
                                            "_".join(["stat_circRNA", prefix + ".csv"])))
                if convert is True:
                    self.converter.convert_circ2gff(
                         os.path.join(candidate_path, prefix, 
                                      "_".join(["circRNA", prefix + ".txt"])), 
                         support, start_ratio, end_ratio, 
                         os.path.join(gff_folder, prefix, "_".join([prefix, "circRNA_all.gff"])), 
                         os.path.join(gff_folder, prefix, "_".join([prefix, "circRNA_best.gff"])))

    def run_circrna(self, align, cores, fastas, gffs, output_folder, read_folder, 
                    stat_folder, convert, support, segemehl_path, samtools_path,
                    start_ratio, end_ratio):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(gffs, gff))
        if segemehl_path is None:
            print("Error: please assign segemehl folder!!")
            sys.exit()
        alignment_path = os.path.join(output_folder, "segemehl_align")
        splice_path = os.path.join(output_folder, "segemehl_splice")
        candidate_path = os.path.join(output_folder, "circRNA_tables")
        gff_folder = os.path.join(output_folder, "gffs")
        self.multiparser._parser_gff(gffs, None)
        gff_path = os.path.join(gffs, "tmp")
        self.multiparser._combine_gff(fastas, gff_path, "fasta", None)
        tmp_reads = []
        if align is True:
            if fastas is False:
                print("Error: There is no genome fasta file!!!")
                sys.exit()
            else:
                fasta_path = os.path.join(fastas, "tmp")
                self.multiparser._parser_fasta(fastas)
                prefixs = []
                self._deal_zip_file(read_folder)
                self._align(fasta_path, segemehl_path, prefixs, 
                            alignment_path, read_folder, cores)
        else:
            fasta_path = os.path.join(fastas, "tmp")
            self.multiparser._parser_fasta(fastas)
            prefixs = []
            for fasta in os.listdir(fasta_path):
                fasta_prefix = fasta.replace(".fa", "")
                prefixs.append(fasta_prefix)
        for prefix in prefixs:
            sub_alignment_path = os.path.join(alignment_path, prefix)
            bam_files = []
            convert_ones = []
            self._convert_sam2bam(sub_alignment_path, samtools_path, convert_ones, bam_files)
            self._merge_sort_aligment_file(bam_files, samtools_path, sub_alignment_path, 
                                           convert_ones, tmp_reads)
            self._run_testrealign(splice_path, prefix, segemehl_path, fasta_path, sub_alignment_path)
        tmp_prefixs = self._merge_bed(fastas, splice_path) ## based on fasta file
        self._stat_and_gen_gff(tmp_prefixs, gff_folder, splice_path, candidate_path, gff_path,
                               stat_folder, support, convert, start_ratio, end_ratio)
        self._remove_tmp(fastas)
        self._remove_tmp(gffs)
