#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
from subprocess import call, Popen
from sort_gff import Sort_GFF
import multiparser
import time
import check_gff_attributes
class CircRNA_detection(object):

    def _check_make_folder(self, path, folder):
        if folder in os.listdir(path):
            call(["rm", "-rf", path + folder])
        call(["mkdir", path + folder])

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

    def run_CircRNA(self, align, cores, fastas, gffs, bin_path, output_folder, 
                    read_folder, stat_folder, convert, support):
        segemehl_path = os.environ["SEGEMEHL_HOME"]
        samtools_path = os.environ["SAMTOOLS_HOME"]
#        for gff in os.listdir(gffs):
#            if gff.endswith(".gff"):
#                check_gff_attributes._check_uni_attributes(gffs + "/" + gff)
        alignment_path = output_folder + "/segemehl_align/"
        splice_path = output_folder + "/segemehl_splice/"
        candidate_path = output_folder + "/circRNA_tables/"
        gff_folder = output_folder + "/gffs/"
        multiparser._parser_gff(gffs, None)
        gff_path = gffs + "/tmp/"
        multiparser._combine_gff(fastas, gff_path, "fasta", None)
        stat_path = stat_folder + "/"
        tmp_reads = []
        if align is True:
            if fastas is False:
                print("Error: There is no genome fasta file!!!")
                sys.exit()
            ###################################
            # generate index of fasta file    #
            ###################################
            else:
                fasta_path = fastas + "/tmp/"
                read_path = read_folder + "/"
                multiparser._parser_fasta(fastas)
                prefixs = []
                for read in os.listdir(read_path):
                    if read.endswith(".bz2"):
                        mod_read = read.replace(".bz2", "")
                        if (".fa" not in mod_read) and \
                           (".fasta" not in mod_read) and \
                           (".fna" not in mod_read):
                            mod_read = mod_read + ".fa"
                        read_out = open(read_path + mod_read, "w")
                        tmp_reads.append(read_path + mod_read)
                        print("unzip " + read)
                        call(["bzcat", read_path + read], stdout=read_out)
                for fasta in os.listdir(fasta_path):
                    index = fasta.replace(".fa", ".idx")
                    call([segemehl_path + "/segemehl.x",
                          "-x", fasta_path + index,
                          "-d", fasta_path + fasta])
                    ################
                    # align reads  #
                    ################
                    processes = []
                    num_process = 0
                    fasta_prefix = fasta.replace(".fa", "")
                    prefixs.append(fasta_prefix)
                    self._check_make_folder(alignment_path, fasta_prefix)
                    for read in os.listdir(read_path):
                        num_process += 1
                        if read.endswith(".fa") or \
                           read.endswith(".fna") or \
                           read.endswith("fasta"):
                            filename = read.split(".")
                            read_prefix = ".".join(filename[:-1])
                            sam_file = read_prefix + "_" + fasta_prefix + ".sam"
                            log_file = read_prefix + "_" + fasta_prefix + ".log"
                            print("mapping " + sam_file)
                            out = open(alignment_path + fasta_prefix + "/" + sam_file, "w")
                            log = open(alignment_path + fasta_prefix + "/" + log_file, "w")
                            p = Popen([segemehl_path + "/segemehl.x",
                                       "-i", fasta_path + index,
                                       "-d", fasta_path + fasta,
                                       "-q", read_path + read, "-S"], stdout=out, stderr=log)
                            processes.append(p)
                            if num_process == cores:
                                self._wait_process(processes)
                                num_process = 0
                    self._wait_process(processes)
        else:
            fasta_path = fastas + "/tmp/"
            read_path = read_folder + "/"
            multiparser._parser_fasta(fastas)
            prefixs = []
            for fasta in os.listdir(fasta_path):
                fasta_prefix = fasta.replace(".fa", "")
                prefixs.append(fasta_prefix)
        #########################################
        # convert sam files to bam files,       #
        # in order to merge to one bam file.    #
        #########################################
#        for prefix in prefixs:
#            sub_alignment_path = alignment_path + prefix + "/"
#            bam_files = []
#            convert_ones = []
#            for sam in os.listdir(sub_alignment_path):
#                if sam.endswith(".sam"):
#                    bam_file = sam.replace(".sam", ".bam")
#                    print("Convert " + sam + " to " + bam_file)
#                    out_bam = open(sub_alignment_path + bam_file, "w")
#                    call([samtools_path + "/samtools",
#                          "view", "-bS", sub_alignment_path + sam], stdout=out_bam)
#                    bam_files.append(sub_alignment_path + bam_file)
#                    convert_ones.append(sub_alignment_path + bam_file)
#                elif sam.endswith(".bam"):
#                    bam_files.append(sub_alignment_path + sam)
#            ##########################################
#            # merge to one bam file and sort it.     #
#            # After that, convert to sam again,      #
#            # in order to run testrealign.x          #
#            ##########################################
#            file_line = " ".join(bam_files)
#            command = " ".join([samtools_path + "/samtools", "merge",
#                                sub_alignment_path + "whole_reads.bam",
#                                file_line])
#            print("Merge all bam files....")
#            os.system(command)
#            print("Sort bam files....")
#            call([samtools_path + "/samtools",
#                  "sort", sub_alignment_path + "whole_reads.bam",
#                  sub_alignment_path + "whole_reads_sorted"])
#            call(["rm", sub_alignment_path + "whole_reads.bam"])
#            print("Convert whole reads bam file to sam file....")
#            call([samtools_path + "/samtools",
#                  "view", "-h", "-o", 
#                  sub_alignment_path + "whole_reads_sorted.sam",
#                  sub_alignment_path + "whole_reads_sorted.bam"])
#            for bam in convert_ones:
#                call(["rm", bam])
#            if len(tmp_reads) != 0:
#                for read in tmp_reads:
#                    call(["rm", read])
#            self._check_make_folder(splice_path, prefix)
#            sub_splice_path = splice_path + prefix + "/"
#            err_log = sub_splice_path + prefix + ".log"
#            print("Running testrealign.x for " + prefix)
#            command = " ".join([segemehl_path + "/testrealign.x",
#                                "-d", fasta_path + prefix + ".fa",
#                                "-q", sub_alignment_path + "whole_reads_sorted.sam",
#                                "-n"])
#            print(command)
#            os.system(command + " 2>" + err_log)
#            os.system("mv *.bed " + sub_splice_path)
#            os.system("rm " + sub_alignment_path + "whole_reads_sorted*")
#        ###############################################
#        # merge splice bedfile based on fasta files   #
#        ###############################################
        tmp_prefixs = []
        for fasta in os.listdir(fastas):
            headers = []
            if fasta.endswith(".fa") or \
               fasta.endswith(".fna") or \
               fasta.endswith(".fasta"):
                with open(fastas + "/" + fasta, "r") as f_h:
                    for line in f_h:
                        line = line.strip()
                        if line.startswith(">"):
                            headers.append(line[1:])
                filename = fasta.split(".")
                fasta_prefix = filename[0]
                tmp_prefixs.append(fasta_prefix)
#                self._check_make_folder(os.getcwd() + "/", fasta_prefix)
#                tmp_fasta_path = fasta_prefix + "/"
#                for header in headers:
#                    call(["cp", splice_path + header + "/splicesites.bed",
#                          tmp_fasta_path + "splicesites_" + header + ".bed"])
#                    call(["cp", splice_path + header + "/transrealigned.bed",
#                          tmp_fasta_path + "transrealigned_" + header + ".bed"])
#                if len(headers) > 1:
#                    out_splice = open(tmp_fasta_path + "splicesites_all.bed", "w")
#                    out_trans = open(tmp_fasta_path + "/transrealigned_all.bed", "w")
#                    for file_ in os.listdir(tmp_fasta_path):
#                        if ("splicesites" in file_) and \
#                           ("splicesites_all" not in file_):
#                            call(["cat", tmp_fasta_path + file_], stdout=out_splice)
#                        elif ("transrealigned" in file_) and \
#                             ("transrealigned_all" not in file_):
#                            call(["cat", tmp_fasta_path + file_], stdout=out_trans)
#                else:
#                    call(["mv", tmp_fasta_path + "splicesites_" + headers[0] + ".bed",
#                          tmp_fasta_path + "splicesites_all.bed"])
#                    call(["mv", tmp_fasta_path + "transrealigned_" + headers[0] + ".bed",
#                          tmp_fasta_path + "transrealigned_all.bed"])
#        os.system("rm -rf " + splice_path + "*")
        for prefix in tmp_prefixs:
            self._check_make_folder(gff_folder, prefix)
#            call(["mv", prefix, splice_path])
            self._check_make_folder(candidate_path, prefix)
            print("comparing with annotation of " + prefix)
            if "splicesites_all.bed" in os.listdir(splice_path + prefix):
                call(["python", bin_path + "/circRNA.py",
                      "-i", splice_path + prefix + "/splicesites_all.bed",
                      "-g", gff_path + prefix + ".gff",
                      "-o", candidate_path + prefix + "/circRNA_" + prefix + ".txt",
                      "-t", stat_path + "stat_circRNA_" + prefix + ".csv"])
                if convert is True:
                    call(["python", bin_path + "/circ2gff.py",
                          "-c", candidate_path + prefix + "/circRNA_" + prefix + ".txt",
                          "-d", str(support),
                          "-oa", gff_folder + prefix + "/" + prefix + "_circRNA_all.gff",
                          "-of", gff_folder + prefix + "/" + prefix + "_circRNA_best.gff"])
        ######################################
        # remove temperary folders and files.#
        ######################################
#        self._remove_tmp(fastas)
#        self._remove_tmp(gffs)
