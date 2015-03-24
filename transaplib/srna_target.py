#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
from subprocess import call, Popen
from sort_gff import Sort_GFF
import multiparser
import time
import check_gff_attributes
class sRNA_target_prediction(object):

    def _remove_tmp(self, folder):
        if folder is not False:
            call(["rm", "-rf", folder + "/tmp"])
            os.system("rm -rf " + folder + "/*_folder")

    def _check_make_folder(self, path, folder):
        if folder in os.listdir(path):
            call(["rm", "-rf", path + folder])
        call(["mkdir", path + folder])

    def _check_gff(self, gffs):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                check_gff_attributes._check_uni_attributes(gffs + "/" + gff)

    def _run_rnaplfold(self, vienna_path, file_type, win_size, span,
                       unstr_region, seq_path, prefix, out_path):
        current = os.getcwd()
        os.chdir(out_path)
        command = " ".join([vienna_path + "/Progs/RNAplfold",
                            "-W", str(win_size),
                            "-L", str(span),
                            "-u", str(unstr_region),
                            "-O"])
        if file_type == "sRNA":
            os.system(command + " < " + current + "/" + seq_path + "tmp_" + prefix + "_" + file_type + ".fa")
        else:
            os.system(command + " < " + current + "/" + seq_path + prefix + "_" + file_type + ".fa")
#        os.system("mv *_openen " + out_path)
        os.chdir(current)
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

    def _organize_files(self, srna_prefix, target_path, file_type, folder):
        self._check_make_folder(os.getcwd() + "/", folder.replace(".gff_folder", "") + file_type)
        tmp_path = os.getcwd() + "/" + folder.replace(".gff_folder", "") + file_type
        for srna in os.listdir(target_path):
            if srna == srna_prefix:
                if file_type != "_merge":
                    os.system("cp " + target_path + srna + "/*.txt " + tmp_path)
                os.system("cp " + target_path + srna + "/*.csv " + tmp_path)

    def _sort_sRNA_fasta(self, fasta, prefix, path):
        out = open(path + "tmp_" + prefix + "_sRNA.fa", "w")
        srnas = []
        with open(fasta) as f_h:
            for line in f_h:
                line = line.strip()
                if line.startswith(">"):
                    name = line[1:]
                else:
                    srnas.append({"name": name, "seq": line, "len": len(line)})
        srnas = sorted(srnas, key = lambda x: (x["len"]))
        for srna in srnas:
            out.write(">" + srna["name"].split("|")[0] + "\n")
            out.write(srna["seq"] + "\n")

    def run_sRNA_target_prediction(self, bin_path, gffs, fastas, srnas,
                                   program, inter_length, win_size_t, span_t, 
                                   win_size_s, span_s, unstr_region_rnaplex_t,
                                   unstr_region_rnaplex_s, unstr_region_rnaup, 
                                   energy, duplex_dist, out_folder, core_plex, core_up):
        vienna_path = os.environ["Vienna_HOME"]
#        self._check_gff(gffs)
#        self._check_gff(srnas)
        target_seq_path = out_folder + "/target_seqs/"
        srna_seq_path = out_folder + "/sRNA_seqs/"
        rnaplex_path = out_folder + "/RNAplex/"
        rnaup_path = out_folder + "/RNAup/"
        merge_path = out_folder + "/merge/"
        multiparser._parser_gff(gffs, None)
        multiparser._parser_fasta(fastas)
        multiparser._parser_gff(srnas, "sRNA")
        srna_path = srnas + "/tmp/"
        fasta_path = fastas + "/tmp/"
        gff_path = gffs + "/tmp/"
        prefixs = []
        #############################################
        # extract sRNA seq                          #
        #############################################
        print("Generating sRNA fasta files...")
        for srna in os.listdir(srna_path):
            if srna.endswith("_sRNA.gff"):
                prefix = srna.replace("_sRNA.gff", "")
                prefixs.append(prefix)
                srna_out = open(srna_seq_path + prefix + "_sRNA.fa", "w")
                call(["python", bin_path + "/get_seq.py",
                      "-s", srna_path + srna, 
                      "-f", fasta_path + prefix + ".fa"],
                      stdout=srna_out)
                self._sort_sRNA_fasta(srna_seq_path + prefix + "_sRNA.fa", prefix, srna_seq_path)
       ###########################################
       # extract target seq                      #
       ###########################################
        print("Generating target fasta files...")
        for gff in os.listdir(gff_path):
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                call(["python", bin_path + "/potential_target.py",
                      "-g", gff_path + gff,
                      "-s", fasta_path + prefix + ".fa",
                      "-o", target_seq_path])
                file_num = 1
                num = 0
                sub_out = open(target_seq_path + prefix +
                               "_target_" + str(file_num) + ".fa", "w")
                with open(target_seq_path + prefix + "_target.fa", "r") as t_f:
                    for line in t_f:
                        line = line.strip()
                        if line.startswith(">"):
                            num += 1
                        if (num == 100):
                            num = 0
                            file_num += 1
                            sub_out.close()
                            sub_out = open(target_seq_path + prefix + 
                                           "_target_" + str(file_num) + ".fa", "w")
                        sub_out.write(line + "\n")
        if (program == "both") or \
           (program == "RNAplex"):
            for prefix in prefixs:
                ########################################
                # Run RNAplfold first for RNAplex      #
                ########################################
                print("Running RNAplfold of " + prefix)
                self._check_make_folder(rnaplex_path, prefix)
                rnaplfold_path = rnaplex_path + prefix + "/RNAplfold"
                call(["mkdir", rnaplfold_path])
                self._run_rnaplfold(vienna_path, "sRNA", win_size_s, span_s,
                                    unstr_region_rnaplex_s, srna_seq_path,
                                    prefix, rnaplfold_path)
                self._run_rnaplfold(vienna_path, "target", win_size_t, span_t,
                                    unstr_region_rnaplex_t, target_seq_path, 
                                    prefix, rnaplfold_path)
                ########################################
                # Run RNAplex                          #
                ########################################
                print("Running RNAplex of " + prefix)
                num_process = 0
                processes = []
                for seq in os.listdir(target_seq_path):
                    if (prefix in seq) and ("_target_" in seq):
                        print(seq)
                        out_rnaplex = open(rnaplex_path + prefix + "/" + prefix + "_RNAplex_" + str(num_process) + ".txt", "w")
                        num_process += 1
                        p = Popen([vienna_path + "/Progs/RNAplex",
                                   "-q", srna_seq_path + "tmp_" + prefix + "_sRNA.fa",
                                   "-t", target_seq_path + seq,
                                   "-l", str(inter_length),
                                   "-e", str(energy), "-z", str(duplex_dist),
                                   "-a", rnaplfold_path], stdout=out_rnaplex)
                        processes.append(p)
                        if num_process % core_plex == 0:
                            self._wait_process(processes)
                self._wait_process(processes)
                if prefix + "_RNAplex.txt" in os.listdir(rnaplex_path + prefix):
                    call(["rm", rnaplex_path + prefix + "/" + prefix + "_RNAplex.txt"])
                for index in range(0, num_process):
                    os.system("cat " + rnaplex_path + prefix + "/" + prefix + "_RNAplex_" + str(index) + ".txt >> " +
                              rnaplex_path + prefix + "/" + prefix + "_RNAplex.txt")
                os.system("rm " + rnaplex_path + prefix + "/*_RNAplex_*")
                out_tmp = open("tmp", "w")
                call(["python", bin_path + "/fix_rnaplex.py",
                      "-i", rnaplex_path + prefix + "/" + prefix + "_RNAplex.txt"],
                      stdout=out_tmp)
                call(["mv", "tmp", rnaplex_path + prefix + "/" + prefix + "_RNAplex.txt"])
        os.system("rm " + target_seq_path + "*_target_*")
        if (program == "both") or \
           (program == "RNAup"):
            for prefix in prefixs:
                #############################################
                # Run RNAup                                 #
                #############################################
                print("Running RNAup of " + prefix)
                self._check_make_folder(rnaup_path, prefix)
                num_up = 0
                processes = []
                if prefix + "_RNAup.txt" in os.listdir(rnaup_path + prefix):
                    call(["rm", rnaup_path + prefix + "/" + prefix + "_RNAup.txt"])
                out_rnaup = rnaup_path + prefix + "/" + prefix + "_RNAup.txt"
                out_log = rnaup_path + prefix + "/" + prefix + "_RNAup.log"
                with open(srna_seq_path + "tmp_" + prefix + "_sRNA.fa", "r") as s_f:
                    for line in s_f:
                        line = line.strip()
                        if line.startswith(">"):
                            print(line)
                            num_up += 1
                            out_up = open("tmp" + str(num_up) + ".fa", "w")
                            out_up.write(line + "\n")
                        else:
                            out_up.write(line + "\n")
                            out_up.close()
                            os.system("cat " + target_seq_path + prefix + "_target.fa >> tmp" + str(num_up) + ".fa")
                            if num_up == core_up:
                                for index in range(1, num_up + 1):
                                    out_tmp_up = open("tmp_rnaup" + str(index) + ".txt", "w")
                                    out_err = open("tmp_log" + str(index) + ".txt", "w")
                                    in_up = open("tmp" + str(index) + ".fa", "r")
                                    p = Popen([vienna_path + "/Progs/RNAup",
                                               "-u", str(unstr_region_rnaup),
                                               "-o", "--interaction_first"], 
                                               stdin = in_up, stdout=out_tmp_up, stderr=out_err)
                                    processes.append(p)
                                self._wait_process(processes)
                                os.system("rm tmp*.fa")
                                for index in range(1, num_up + 1):
                                    os.system("cat tmp_rnaup" + str(index) + ".txt >> " + out_rnaup)
                                    os.system("cat tmp_log" + str(index) + ".txt >> " + out_log)
                                os.system("rm tmp*.txt")
                                processes = []
                                num_up = 0
                for index in range(1, num_up + 1):
                    out_tmp_up = open("tmp_rnaup" + str(index) + ".txt", "w")
                    out_err = open("tmp_log" + str(index) + ".txt", "w")
                    in_up = open("tmp" + str(index) + ".fa", "r")
                    p = Popen([vienna_path + "/Progs/RNAup",
                               "-u", str(unstr_region_rnaup),
                               "-o", "--interaction_first"], 
                               stdin = in_up, stdout=out_tmp_up, stderr=out_err)
                    processes.append(p)
                if len(processes) != 0:
                    self._wait_process(processes)
                    os.system("rm tmp*.fa")
                    for index in range(1, num_up + 1):
                        os.system("cat tmp_rnaup" + str(index) + ".txt >> " + out_rnaup)
                        os.system("cat tmp_log" + str(index) + ".txt >> " + out_log)
                    os.system("rm tmp*.txt")
        ###############################################
        # merge and rank the sRNA target output       #
        ###############################################
        for prefix in prefixs:
            self._check_make_folder(merge_path, prefix)
            print("Ranking " + prefix + " now...")
            command = ["python", bin_path + "/merge_rnaplex_rnaup.py",
                       "-sg", srna_path + prefix + "_sRNA.gff",
                       "-g", gff_path + prefix + ".gff"]
            if program == "both":
                call(command + ["-p", rnaplex_path + "/" + prefix + "/" + prefix + "_RNAplex.txt",
                                "-u", rnaup_path + "/" + prefix + "/" + prefix + "_RNAup.txt",
                                "-op", rnaplex_path + prefix + "/" + prefix + "_RNAplex_rank.csv",
                                "-ou", rnaup_path + prefix + "/" + prefix + "_RNAup_rank.csv",
                                "-o", merge_path + prefix + "/" + prefix + "_merge.csv",
                                "-ol", merge_path + prefix + "/" + prefix + "_overlap.csv"])
            elif program == "RNAplex":
                call(command + ["-p", rnaplex_path + "/" + prefix + "/" + prefix + "_RNAplex.txt",
                                "-op", rnaplex_path + prefix + "/" + prefix + "_RNAplex_rank.csv"])
            elif program == "RNAup":
                call(command + ["-u", rnaup_path + "/" + prefix + "/" + prefix + "_RNAup.txt",
                                "-ou", rnaup_path + prefix + "/" + prefix + "_RNAup_rank.csv"])
        ########################################################
        # re-organize the files based on annotation gff files. #
        ########################################################
        srna_prefixs = []
        folders = []
        for folder in os.listdir(gffs):
            if folder.endswith(".gff_folder"):
                folders.append(folder.replace(".gff_folder", ""))
                for gff in os.listdir(gffs + "/" + folder):
                    if gff.endswith(".gff"):
                       srna_prefixs.append(gff.replace(".gff", ""))
                for srna_prefix in srna_prefixs:
                    self._organize_files(srna_prefix, rnaplex_path, "_RNAplex", folder)
                    self._organize_files(srna_prefix, rnaup_path, "_RNAup", folder)
                    self._organize_files(srna_prefix, merge_path, "_merge", folder)
        ############################
        # remove temperary files   #
        ############################
        if (program == "RNAplex") or \
           (program == "both"):
            for strain in os.listdir(out_folder + "/RNAplex"):
                os.system("rm -rf " + out_folder + "/RNAplex/" + strain + "/RNAplfold")
        os.system("rm -rf " + os.getcwd() + "/tmp*")
        self._remove_tmp(gffs)
        self._remove_tmp(srnas)
        self._remove_tmp(fastas)
