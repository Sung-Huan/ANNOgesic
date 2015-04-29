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

    def __init__(self, output_folder, gffs, fastas, align):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.converter = Converter()
        self.alignment_path = os.path.join(output_folder, "segemehl_align")
        self.splice_path = os.path.join(output_folder, "segemehl_splice")
        self.candidate_path = os.path.join(output_folder, "circRNA_tables")
        self.gff_folder = os.path.join(output_folder, "gffs")
        self.gff_path = os.path.join(gffs, "tmp")
        self.splice_all_file = "splicesites_all.bed"
        self.splice_file = "splicesites.bed"
        self.splice_all = "splicesites_all"
        self.splice = "splicesites"
        self.trans_all_file = "transrealigned_all.bed"
        self.trans_file = "transrealigned.bed"
        self.trans_all = "transrealigned_all"
        self.trans = "transrealigned"
        self.whole_bam = "whole_reads.bam"
        self.sort_bam = "whole_reads_sort"
        if align is True:
            if fastas is False:
                print("Error: There is no genome fasta file!!!")
                sys.exit()
            else:
                self.fasta_path = os.path.join(fastas, "tmp")
        else:
            self.fasta_path = os.path.join(fastas, "tmp")

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
            self.helper.check_make_folder(os.path.join(alignment_path, fasta_prefix))
            align_files = []
            for read in os.listdir(read_folder):
                num_process += 1
                if read.endswith(".fa") or \
                   read.endswith(".fna") or \
                   read.endswith("fasta"):
                    filename = read.split(".")
                    read_prefix = ".".join(filename[:-1])
                    sam_file = "_".join([read_prefix, fasta_prefix + ".sam"])
                    log_file = "_".join([read_prefix, fasta_prefix + ".log"])
                    align_files.append("_".join([read_prefix, fasta_prefix]))
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
        return align_files

    def _convert_sam2bam(self, sub_alignment_path, samtools_path, 
                         convert_ones, bam_files, align_files, remove_ones):
        for sam in os.listdir(sub_alignment_path):
            if sam.endswith(".sam"):
                bam_file = sam.replace(".sam", ".bam")
                print("Convert {0} to {1}".format(sam, bam_file))
                out_bam = os.path.join(sub_alignment_path, bam_file)
                call([samtools_path, "view", "-bS", 
                      os.path.join(sub_alignment_path, sam), "-o", out_bam])
                bam_files.append(os.path.join(sub_alignment_path, bam_file))
                if align_files:
                    if bam_file.replace(".bam", "") not in align_files:
                        convert_ones.append(os.path.join(sub_alignment_path, bam_file))
                    else:
                        remove_ones.append(os.path.join(sub_alignment_path, sam))
            elif sam.endswith(".bam"):
                bam_files.append(os.path.join(sub_alignment_path, sam))
            elif sam.endswith(".log"):
                os.remove(os.path.join(sub_alignment_path, sam))

    def _merge_sort_aligment_file(self, bam_files, samtools_path, sub_alignment_path, 
                                  convert_ones, tmp_reads, remove_ones):
        print("Merge all bam files....")
        whole_bam = os.path.join(sub_alignment_path, self.whole_bam)
        file_line = " ".join(bam_files)
        os.system(" ".join([samtools_path, "merge",
                            whole_bam, file_line]))
        print("Sort bam files....")
        call([samtools_path, "sort", whole_bam,
              os.path.join(sub_alignment_path, self.sort_bam)])
        os.remove(os.path.join(sub_alignment_path, self.whole_bam))
        print("Convert whole reads bam file to sam file....")
        call([samtools_path, "view", "-h", "-o",
              os.path.join(sub_alignment_path, self.sort_bam + ".sam"),
              os.path.join(sub_alignment_path, self.sort_bam + ".bam")])
        for bam in convert_ones:
            os.remove(bam)
        for sam in remove_ones:
            os.remove(sam)
        if len(tmp_reads) != 0:
            for read in tmp_reads:
                os.remove(read)

    def _run_testrealign(self, splice_path, prefix, segemehl_path, fasta_path, sub_alignment_path):
        self.helper.check_make_folder(os.path.join(splice_path, prefix))
        sub_splice_path = os.path.join(splice_path, prefix)
        err_log = os.path.join(sub_splice_path, prefix + ".log")
        print("Running testrealign.x for {0}".format(prefix))
        command = " ".join([os.path.join(segemehl_path, "testrealign.x"),
                            "-d", os.path.join(fasta_path, prefix + ".fa"),
                            "-q", os.path.join(sub_alignment_path, self.sort_bam + ".sam"),
                            "-n"])
        os.system(command + " 2>" + err_log)
        self.helper.move_all_content(os.getcwd(), sub_splice_path, ".bed")
        self.helper.remove_all_content(sub_alignment_path, self.sort_bam, "file")

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
                self.helper.check_make_folder(os.path.join(os.getcwd(), fasta_prefix))
                for header in headers:
                    shutil.copyfile(os.path.join(splice_path, header, self.splice_file),
                                    os.path.join(fasta_prefix, 
                                    "_".join([self.splice, header + ".bed"])))
                    shutil.copyfile(os.path.join(splice_path, header, self.trans_file),
                                    os.path.join(fasta_prefix, 
                                    "_".join([self.trans, header + ".bed"])))
                out_splice = os.path.join(fasta_prefix, self.splice_all_file)
                out_trans = os.path.join(fasta_prefix, self.trans_all_file)
                if len(headers) > 1:
                    for file_ in os.listdir(fasta_prefix):
                        if (self.splice in file_) and \
                           (self.splice_all not in file_):
                            self.helper.merge_file(os.path.join(fasta_prefix, file_),
                                                   out_splice)
                        elif (self.trans in file_) and \
                             (self.trans_all not in file_):
                            self.helper.merge_file(os.path.join(fasta_prefix, file_), 
                                                   out.trans)
                else:
                    os.rename(os.path.join(fasta_prefix, 
                              "_".join([self.splice, headers[0] + ".bed"])), out_splice)
                    os.rename(os.path.join(fasta_prefix,
                              "_".join([self.trans, headers[0] + ".bed"])), out_trans)
        self.helper.remove_all_content(splice_path, None, "dir")
        return tmp_prefixs

    def _stat_and_gen_gff(self, tmp_prefixs, gff_folder, splice_path, candidate_path, gff_path,
                          stat_folder, support, convert, start_ratio, end_ratio):
        for prefix in tmp_prefixs:
            self.helper.check_make_folder(os.path.join(gff_folder, prefix))
            shutil.copytree(prefix, os.path.join(splice_path, prefix))
            self.helper.check_make_folder(os.path.join(candidate_path, prefix))
            print("comparing with annotation of {0}".format(prefix))
            if self.splice_all_file in os.listdir(os.path.join(splice_path, prefix)):
                detect_circrna(os.path.join(splice_path, prefix, self.splice_all_file), 
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
        self.multiparser._parser_gff(gffs, None)
        self.multiparser._combine_gff(fastas, self.gff_path, "fasta", None)
        tmp_reads = []
        if align is True:
            self.multiparser._parser_fasta(fastas)
            prefixs = []
            self._deal_zip_file(read_folder)
            align_files = self._align(self.fasta_path, segemehl_path, prefixs, 
                                      self.alignment_path, read_folder, cores)
        else:
            self.multiparser._parser_fasta(fastas)
            prefixs = []
            for fasta in os.listdir(self.fasta_path):
                fasta_prefix = fasta.replace(".fa", "")
                prefixs.append(fasta_prefix)
            align_files = None
        for prefix in prefixs:
            sub_alignment_path = os.path.join(self.alignment_path, prefix)
            bam_files = []
            convert_ones = []
            remove_ones = []
            self._convert_sam2bam(sub_alignment_path, samtools_path, convert_ones, 
                                  bam_files, align_files, remove_ones)
            self._merge_sort_aligment_file(bam_files, samtools_path, sub_alignment_path, 
                                           convert_ones, tmp_reads, remove_ones)
            self._run_testrealign(self.splice_path, prefix, segemehl_path, 
                                  self.fasta_path, sub_alignment_path)
        tmp_prefixs = self._merge_bed(fastas, self.splice_path) ## based on fasta file
        self.multiparser._parser_gff(gffs, None)
        self.multiparser._combine_gff(fastas, self.gff_path, "fasta", None)
        self._stat_and_gen_gff(tmp_prefixs, self.gff_folder, self.splice_path, 
                               self.candidate_path, self.gff_path,
                               stat_folder, support, convert, start_ratio, end_ratio)
        self.helper.remove_tmp(fastas)
        self.helper.remove_tmp(gffs)
        for tmp_prefix in tmp_prefixs:
            shutil.rmtree(tmp_prefix)
