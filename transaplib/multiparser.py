#!/usr/bin/python
from Bio import SeqIO
import os        
import sys
from subprocess import call
import csv
import shutil
from transaplib.seq_editer import Seq_Editer
from transaplib.helper import Helper


class Multiparser(object):

    def __init__(self):
        self.seq_editer = Seq_Editer()
        self.helper = Helper()

    def _fix_folder_name(self, folder):
        if folder.endswith("/"):
            return folder[:-1]
        else:
            return folder
    
    def _combine_fasta(self, ref_folder, tar_folder, ref_feature):
        change = False
#        ref_folder = self._fix_folder_name(ref_folder)
#        tar_folder = self._fix_folder_name(tar_folder)
        if ref_feature is None:
            ref_feature = ""
        else:
            ref_feature = "_" + ref_feature
        self.helper.check_make_folder(tar_folder, "merge_tmp")
        for folder in os.listdir(ref_folder):
            files = []
            if "_folder" in folder:
                datas = folder.split("_folder")
                if ref_feature == "":
                    prefix = datas[0][:-4]
                elif ref_feature == "_fasta":
                    if datas[0].endswith(".fa"):
                        prefix = datas[0][:-3]
                    elif datas[0].endswith(".fna"):
                        prefix = datas[0][:-4]
                    elif datas[0].endswith(".fasta"):
                        prefix = datas[0][:-6]
                else:
                    datas = datas[0][:-4]
                    datas = datas.split(ref_feature)
                    prefix = datas[0]
                print("Merging fasta file of " + prefix)
                for file_ in os.listdir("/".join([ref_folder, folder])):
                    if ref_feature == "":
                        files.append(file_[:-4])
                    elif ref_feature == "_fasta":
                        files.append(file_[:-3])
                    else:
                        filename = file_.split(ref_feature)
                        files.append(filename[0])
                for tar in os.listdir(tar_folder):
                    if tar.endswith(".fa") or \
                       tar.endswith(".fna") or \
                       tar.endswith(".fasta"):
                        filename = ".".join((tar.split("."))[:-1])
                        for file_ in files:
                            if filename == file_:
                                self.helper.merge_file(tar_folder, tar, tar_folder, "tmp.fa")
                                change = True
                if change:
                    change = False
                    os.rename(os.path.join(tar_folder, "tmp.fa"), 
                              os.path.join(tar_folder, "merge_tmp", prefix + ".fa"))
        self.helper.remove_all_content(tar_folder, ".fa", "file")
        self.helper.move_all_content(os.path.join(tar_folder, "merge_tmp"), tar_folder, None)
        shutil.rmtree(os.path.join(tar_folder, "merge_tmp"))
    
    def _combine_wig(self, ref_folder, tar_folder, ref_feature):
        change_f = False
        change_r = False
#        ref_folder = self._fix_folder_name(ref_folder)
#        tar_folder = self._fix_folder_name(tar_folder)
        if ref_feature is None:
            ref_feature = ""
        else:
            ref_feature = "_" + ref_feature
        self.helper.check_make_folder(tar_folder, "merge_tmp")
        for folder in os.listdir(ref_folder):
            files = []
            if "_folder" in folder:
                datas = folder.split("_folder")
                if ref_feature == "":
                    prefix = datas[0][:-4]
                elif ref_feature == "_fasta":
                    if datas[0].endswith(".fa"):
                        prefix = datas[0][:-3]
                    elif datas[0].endswith(".fna"):
                        prefix = datas[0][:-4]
                    elif datas[0].endswith(".fasta"):
                        prefix = datas[0][:-6]
                else:
                    datas = datas[0][:-4]
                    datas = datas.split(ref_feature)
                    prefix = datas[0]
                print("Merging wig file of " + prefix)
                for file_ in os.listdir(os.path.join(ref_folder, folder)):
                    if ref_feature == "":
                        files.append(file_[:-4])
                    elif ref_feature == "_fasta":
                        files.append(file_[:-3])
                    else:
                        filename = file_.split(ref_feature)
                        files.append(filename[0])
                for tar in os.listdir(tar_folder):
                    filename = tar.split("_STRAIN_")
                    for file_ in files:
                        if (tar.endswith(".wig")) and (file_ == filename[-1][:-4]):
                            if ("forward" in tar) and ("reverse" in tar):
                                print("Error: Unclear wig file. It is reverse or forward!!!")
                            elif ("forward" in tar):
                                self.helper.merge_file(tar_folder, tar, tar_folder, "tmp_forward.wig")
                                change_f = True
                            elif ("reverse" in tar):
                                self.helper.merge_file(tar_folder, tar, tar_folder, "tmp_reverse.wig")
                                change_r = True
                if change_f and change_r:
                    change_f = False
                    change_r = False
                    os.rename(os.path.join(tar_folder, "tmp_forward.wig"), 
                              os.path.join(tar_folder, "merge_tmp", prefix + "_forward.wig"))
                    os.rename(os.path.join(tar_folder, "tmp_reverse.wig"), 
                              os.path.join(tar_folder, "merge_tmp", prefix + "_reverse.wig"))
        self.helper.remove_all_content(tar_folder, ".wig", "file")
        self.helper.move_all_content(os.path.join(tar_folder, "merge_tmp"),
                                     tar_folder, None)
        shutil.rmtree(os.path.join(tar_folder, "merge_tmp"))
    
    
    def _combine_gff(self, ref_folder, tar_folder, ref_feature, tar_feature):
#        ref_folder = self._fix_folder_name(ref_folder)
#        tar_folder = self._fix_folder_name(tar_folder)
        change = False
        if tar_feature is None:
            tar_feature = ""
        else:
            tar_feature = "_" + tar_feature
        if ref_feature is None:
            ref_feature = ""
        else:
            ref_feature = "_" + ref_feature
        self.helper.check_make_folder(tar_folder, "merge_tmp")
        for folder in os.listdir(ref_folder):
            files = []
            if "_folder" in folder:
                datas = folder.split("_folder")
                if ref_feature == "":
                    prefix = datas[0][:-4]
                elif ref_feature == "_fasta":
                    if datas[0].endswith(".fa"):
                        prefix = datas[0][:-3]
                    elif datas[0].endswith(".fna"):
                        prefix = datas[0][:-4]
                    elif datas[0].endswith(".fasta"):
                        prefix = datas[0][:-6]
                else:
                    datas = datas[0][:-4]
                    datas = datas.split(ref_feature)
                    prefix = datas[0]
                print("Merging gff file of " + prefix + tar_feature)
                for file_ in os.listdir(os.path.join(ref_folder, folder)):
                    if ref_feature == "":
                        files.append(file_[:-4])
                    elif ref_feature == "_fasta":
                        files.append(file_[:-3])
                    else:
                        filename = file_.split(ref_feature)
                        files.append(filename[0])
                for tar in os.listdir(tar_folder):
                    for file_ in files:
                        if (".gff" in tar) and (file_ + tar_feature == tar[:-4]):
                            self.helper.merge_file(tar_folder, tar, tar_folder, "tmp.gff")
                            change = True
                if change:
                    change = False
#                    print("/".join([tar_folder, "tmp.gff"]))
#                    print("/".join([tar_folder, "merge_tmp", prefix + tar_feature + ".gff"]))
                    os.rename(os.path.join(tar_folder, "tmp.gff"),
                              os.path.join(tar_folder, "merge_tmp", prefix + tar_feature + ".gff"))
        self.helper.remove_all_content(tar_folder, ".gff", "file")
        self.helper.move_all_content(os.path.join(tar_folder, "merge_tmp"), tar_folder, None)
        shutil.rmtree(os.path.join(tar_folder, "merge_tmp"))
    def _parser_fasta(self, fastas):
        first = True
#        fastas = self._fix_folder_name(fastas)
        ### fix header ###
        for fasta in os.listdir(fastas):
            if fasta.endswith("fasta") or \
               fasta.endswith("fa") or \
               fasta.endswith("fna"):
                self.seq_editer.modify_header(os.path.join(fastas, fasta))
        self.helper.check_make_folder(fastas, "tmp")
        for fasta in os.listdir(fastas):
            if ("_folder" not in fasta) and ("tmp" != fasta):
                if (fasta.endswith(".fa")) or \
                   (fasta.endswith(".fna")) or \
                   (fasta.endswith(".fasta")):
                    out_path = os.path.join(fastas, fasta + "_folder")
                    print("Parser " + fasta + "...")
                    self.helper.check_make_folder(fastas, fasta + "_folder")
                    with open(os.path.join(fastas, fasta), "r") as f_f:
                        for line in f_f:
                            if line[0] == ">":
                                line = line.strip()
                                if "|" in line:
                                    strain = line.split("|")
                                    name = strain[3]
                                else:
                                    name = line[1:]
                                if first:
                                    first = False
                                else:
                                    out.close()
                                    out_t.close()
                                out = open(os.path.join(out_path, name + ".fa"), "w")
                                out_t = open(os.path.join(fastas, "tmp", name + ".fa"), "w")
                                out.write(">" + name + "\n")
                                out_t.write(">" + name + "\n")
                            else:
                                out.write(line)
                                out_t.write(line)
    def _parser_gff(self, gff_folder, feature):
#        gff_folder = self._fix_folder_name(gff_folder)
        first = True
        if feature is None:
            feature = ""
        else:
            feature = "_" + feature
        self.helper.check_make_folder(gff_folder, "tmp")
        for filename in os.listdir(gff_folder):
            pre_seq_id = ""
            if ("_folder" not in filename) and ("tmp" != filename):
                out_path = os.path.join(gff_folder, filename + "_folder")
                if ".gff" in filename:
                    print("Parser " + filename + "...")
                    self.helper.check_make_folder(gff_folder, filename + "_folder")
                    fh = open(os.path.join(gff_folder, filename), "r")
                    for row in csv.reader(fh, delimiter="\t"):
                        if row[0].startswith("#"):
                            continue
                        else:
                            if pre_seq_id == row[0]:
                                out.write("\t".join(row) + "\n")
                                out_t.write("\t".join(row) + "\n")
                            else:
                                if first:
                                    first = False
                                else:
                                    out.close()
                                    out_t.close()
                                out = open(os.path.join(out_path, row[0] + feature + ".gff"), "w")
                                out_t = open(os.path.join(gff_folder, "tmp", row[0] + feature + ".gff"), "w")
                                pre_seq_id = row[0]
                                out.write("\t".join(row) + "\n")
                                out_t.write("\t".join(row) + "\n")
    
    def _parser_wig(self, wig_folder):
        first = True
#        wig_folder = self._fix_folder_name(wig_folder)
        self.helper.check_make_folder(wig_folder, "tmp")
        for filename in os.listdir(wig_folder):
            track_info = ""
            if ("_folder" not in filename) and ("tmp" != filename):
                out_path = os.path.join(wig_folder, filename + "_folder")
                if ".wig" in filename:
                    print("Parser " + filename + "...")
                    self.helper.check_make_folder(wig_folder, filename + "_folder")
                    with open(os.path.join(wig_folder, filename), "r") as w_f:
                        for line in w_f:
                            line = line.split(" ")
                            if (line[0] == "track"):
                                track_info = " ".join(line)
                            if (line[0] == "variableStep"):
                                strain = line[1].split("=")
                                if first:
                                    first = False
                                else:
                                    out.close()
                                    out_t.close()
                                out = open("".join([out_path, filename[:-4], 
                                           "_STRAIN_", strain[1], ".wig"]), "w")
                                out_t = open("".join([os.path.join(wig_folder, "tmp", filename[:-4]), 
                                             "_STRAIN_", strain[1], ".wig"]), "w")
                                if track_info != "":
                                    out.write(track_info)
                                    out_t.write(track_info)
                                out.write(" ".join(line))
                                out_t.write(" ".join(line))
                            if (line[0] != "track") and \
                               (line[0] != "variableStep"):
                                out.write(" ".join(line))
                                out_t.write(" ".join(line))
