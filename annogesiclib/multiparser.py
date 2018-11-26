import os
import sys
import csv
import shutil
from annogesiclib.seq_editer import SeqEditer
from annogesiclib.helper import Helper


class Multiparser(object):

    def __init__(self):
        self.seq_editer = SeqEditer()
        self.helper = Helper()
        self.tmp_fa = "tmp.fa"
        self.tmp_gff = "tmp.gff"
        self.tmp_wig_forward = "tmp_forward.wig"
        self.tmp_wig_reverse = "tmp_reverse.wig"

    def combine_fasta(self, ref_folder, tar_folder, ref_feature):
        '''combine multiple fasta files'''
        tar_merge = os.path.join(tar_folder, "merge_tmp")
        change = False
        if ref_feature is None:
            ref_feature = ""
        else:
            ref_feature = "_" + ref_feature
        self.helper.check_make_folder(tar_merge)
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
                print("Merging fasta files of " + prefix)
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
                                self.helper.merge_file(
                                     os.path.join(tar_folder, tar),
                                     os.path.join(tar_folder, self.tmp_fa))
                                change = True
                if change:
                    change = False
                    shutil.move(os.path.join(tar_folder, self.tmp_fa),
                                os.path.join(tar_merge, prefix + ".fa"))
        self.helper.remove_all_content(tar_folder, ".fa", "file")
        self.helper.move_all_content(tar_merge, tar_folder, None)
        shutil.rmtree(tar_merge)

    def get_prefix(self, folder, ref_feature):
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
        return prefix

    def combine_wig(self, ref_folder, tar_folder, ref_feature, libs):
        '''combine multiple wig files'''
        tar_merge = os.path.join(tar_folder, "merge_tmp")
        change_f = False
        change_r = False
        if ref_feature is None:
            ref_feature = ""
        else:
            ref_feature = "_" + ref_feature
        self.helper.check_make_folder(tar_merge)
        for folder in os.listdir(ref_folder):
            files = []
            if "_folder" in folder:
                prefix = self.get_prefix(folder, ref_feature)
                print("Merging wig files of " + prefix)
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
                        if (tar.endswith(".wig")) and (
                                file_ == filename[-1][:-4]):
                            for lib in libs:
                                if (filename[0] in lib) and (lib[-1] == "+"):
                                    self.helper.merge_file(
                                        os.path.join(tar_folder, tar),
                                        os.path.join(tar_folder,
                                                     self.tmp_wig_forward))
                                    change_f = True
                                elif (filename[0] in lib) and (lib[-1] == "-"):
                                    self.helper.merge_file(
                                        os.path.join(tar_folder, tar),
                                        os.path.join(tar_folder,
                                                     self.tmp_wig_reverse))
                                    change_r = True
                if change_f and change_r:
                    change_f = False
                    change_r = False
                    shutil.move(os.path.join(tar_folder, self.tmp_wig_forward),
                                os.path.join(tar_merge,
                                             prefix + "_forward.wig"))
                    shutil.move(os.path.join(tar_folder, self.tmp_wig_reverse),
                                os.path.join(tar_merge,
                                             prefix + "_reverse.wig"))
                else:
                    print("Error: comparing input files of {0} failed. "
                          "Please check the seq IDs of all gff and fasta "
                          "files, they should be the same.\nPlease "
                          "also check the wiggle files which should contain "
                          "forward and reverse files.".format(prefix))
                    sys.exit()
        self.helper.remove_all_content(tar_folder, ".wig", "file")
        self.helper.move_all_content(tar_merge, tar_folder, None)
        shutil.rmtree(tar_merge)

    def combine_gff(self, ref_folder, tar_folder, ref_feature, tar_feature):
        '''combine multiple gff files'''
        tar_merge = os.path.join(tar_folder, "merge_tmp")
        change = False
        if tar_feature is None:
            tar_feature = ""
        else:
            tar_feature = "_" + tar_feature
        if ref_feature is None:
            ref_feature = ""
        else:
            ref_feature = "_" + ref_feature
        self.helper.check_make_folder(tar_merge)
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
                print("Merging gff files of " + prefix + tar_feature)
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
                        if (".gff" in tar) and (
                                file_ + tar_feature == tar[:-4]):
                            self.helper.merge_file(
                                 os.path.join(tar_folder, tar),
                                 os.path.join(tar_folder, self.tmp_gff))
                            change = True
                if change:
                    change = False
                    shutil.move(os.path.join(tar_folder, self.tmp_gff),
                                os.path.join(tar_folder, "merge_tmp",
                                prefix + tar_feature + ".gff"))
        self.helper.remove_all_content(tar_folder, ".gff", "file")
        self.helper.move_all_content(tar_merge, tar_folder, None)
        shutil.rmtree(tar_merge)

    def parser_fasta(self, fastas):
        '''parser the fasta file based on strain'''
        par_tmp = os.path.join(fastas, "tmp")
        first = True
        out = None
        out_t = None
        detect = False
        for fasta in os.listdir(fastas):
            if (fasta.endswith(".fasta") or
                    fasta.endswith(".fa") or
                    fasta.endswith(".fna")):
                detect = True
                self.seq_editer.modify_header(os.path.join(fastas, fasta))
        self.helper.check_make_folder(par_tmp)
        if not detect:
            print("Error: there are folders which conatin no fasta files! "
                  "The files should end with .fa or .fna or .fasta!")
            sys.exit()
        for fasta in os.listdir(fastas):
            if ("_folder" not in fasta) and ("tmp" != fasta):
                if (fasta.endswith(".fa")) or \
                   (fasta.endswith(".fna")) or \
                   (fasta.endswith(".fasta")):
                    out_path = os.path.join(fastas, fasta + "_folder")
                    print("Parsing " + fasta)
                    self.helper.check_make_folder(out_path)
                    with open(os.path.join(fastas, fasta), "r") as f_f:
                        for line in f_f:
                            if line[0] == ">":
                                line = line.strip()
                                if ("|" in line) and (
                                        len(line.split("|")) > 4):
                                    strain = line.split("|")
                                    name = strain[3]
                                else:
                                    name = line[1:]
                                if first:
                                    first = False
                                else:
                                    out.close()
                                    out_t.close()
                                out = open(os.path.join(
                                           out_path, name + ".fa"), "w")
                                out_t = open(os.path.join(
                                             par_tmp, name + ".fa"), "w")
                                out.write(">" + name + "\n")
                                out_t.write(">" + name + "\n")
                            else:
                                out.write(line)
                                out_t.write(line)
        if out is not None:
            out.close()
        if out_t is not None:
            out_t.close()

    def parser_gff(self, gff_folder, feature):
        '''parser gff file based on strain'''
        par_tmp = os.path.join(gff_folder, "tmp")
        out = None
        out_t = None
        first = True
        detect = False
        if feature is None:
            feature = ""
        else:
            feature = "_" + feature
        self.helper.check_make_folder(par_tmp)
        for filename in os.listdir(gff_folder):
            pre_seq_id = ""
            if ("_folder" not in filename) and ("tmp" != filename):
                out_path = os.path.join(gff_folder, filename + "_folder")
                if ".gff" in filename:
                    detect = True
                    print("Parsing " + filename)
                    self.helper.check_make_folder(out_path)
                    self.helper.sort_gff(os.path.join(gff_folder, filename),
                                         os.path.join(gff_folder, "tmp.gff"))
                    f_h = open(os.path.join(gff_folder, "tmp.gff"), "r")
                    for row in csv.reader(f_h, delimiter="\t"):
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
                                out = open(os.path.join(out_path,
                                           row[0] + feature + ".gff"), "w")
                                out_t = open(os.path.join(par_tmp,
                                             row[0] + feature + ".gff"), "w")
                                pre_seq_id = row[0]
                                out.write("\t".join(row) + "\n")
                                out_t.write("\t".join(row) + "\n")
                    f_h.close()
        if not detect:
            print("Error: There are folders which contain no gff3 files! "
                  "The files should end with .gff!")
            sys.exit()
        if os.path.exists(os.path.join(gff_folder, "tmp.gff")):
            os.remove(os.path.join(gff_folder, "tmp.gff"))
        if out is not None:
            out.close()
        if out_t is not None:
            out_t.close()

    def parser_wig(self, wig_folder):
        '''parser the wig file based on strain'''
        par_tmp = os.path.join(wig_folder, "tmp")
        first = True
        out = None
        out_t = None
        detect = False
        self.helper.check_make_folder(par_tmp)
        for filename in os.listdir(wig_folder):
            track_info = ""
            if ("_folder" not in filename) and ("tmp" != filename):
                out_path = os.path.join(wig_folder, filename + "_folder")
                if ".wig" in filename:
                    detect = True
                    print("Parsing {0}".format(filename))
                    self.helper.check_make_folder(out_path)
                    with open(os.path.join(wig_folder, filename), "r") as w_f:
                        for line in w_f:
                            if (not line.startswith("#")) and (len(line) != 0):
                                line = line.split(" ")
                                if (line[0] == "track"):
                                    track_info = " ".join(line)
                                if (line[0] == "variableStep") or (line[0] == "fixedStep"):
                                    strain = line[1].split("=")
                                    if first:
                                        first = False
                                    else:
                                        out.close()
                                        out_t.close()
                                    out = open("".join([
                                        os.path.join(out_path, filename[:-4]),
                                        "_STRAIN_", strain[1], ".wig"]), "w")
                                    out_t = open("".join([
                                        os.path.join(wig_folder, "tmp",
                                                     filename[:-4]),
                                        "_STRAIN_", strain[1], ".wig"]), "w")
                                    if track_info != "":
                                        out.write(track_info)
                                        out_t.write(track_info)
                                    out.write(" ".join(line))
                                    out_t.write(" ".join(line))
                                if (line[0] != "track") and (
                                        line[0] != "variableStep"):
                                    out.write(" ".join(line))
                                    out_t.write(" ".join(line))
        if not detect:
            print("Error: There are folders which contain no wig files! "
                  "The files should end with .wig!")
            sys.exit()
        if out is not None:
            out.close()
        if out_t is not None:
            out_t.close()
