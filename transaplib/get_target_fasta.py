#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
import shutil
from transaplib.multiparser import Multiparser
from transaplib.seq_editer import Seq_Editer
from transaplib.helper import Helper


class Target_fasta(object):

    def __init__(self):
        self.multiparser = Multiparser()
        self.seq_editer = Seq_Editer()
        self.helper = Helper()

    def get_target_fasta(self, mut_table, tar_folder, ref_folder, output):
        self.multiparser._parser_fasta(ref_folder)
        tar_fastas = os.path.join(tar_folder, "tmp")
        ref_fastas = os.path.join(ref_folder, "tmp")
        if "tmp" in os.listdir(tar_folder):
            shutil.rmtree(tar_fastas)
        os.mkdir(tar_fastas)
        self.seq_editer.modify_seq(ref_fastas, mut_table, tar_fastas)
        print("transfer to target fasta...")
        if output is not False:
            for file_ in output:
                datas = file_.split(":")
                filename = datas[0]
                strains = datas[1].split(",")
                out = open(os.path.join(tar_folder, filename + ".fa"), "w")
                for strain in strains:
                    if strain + ".fa" in os.listdir(tar_fastas):
                        with open(os.path.join(tar_fastas, strain + ".fa")) as f_h:
                            for line in f_h:
                                out.write(line)
                    else:
                        print("Error:There is no fasta information of " + strain + ".fa")
                out.close()
        else:
            os.rename(tar_fastas + "*.fa ", tar_folder)
        shutil.rmtree(os.path.join(tar_folder, "tmp"))
        shutil.rmtree(os.path.join(ref_folder,"tmp"))
        self.helper.remove_all_content(ref_folder, "_folder", "dir")
        print("please use the new fasta file to remapping again.")
        print("Then copy BAMs and wigs back to input/align_results/BAMs and input/align_results/wigs")
