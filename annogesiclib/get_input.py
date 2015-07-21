import os
import sys
from annogesiclib.seq_editer import SeqEditer


def wget(input_folder, ftp, files_type):
    os.system(" ".join(["wget", "-cP", input_folder, ftp + "/*" + files_type]))

def get_file(ftp, input_folder, files_type):
    """Download required files from FTP."""
    detect = False
    wget(input_folder, ftp, files_type)
    for file_ in os.listdir(input_folder):
        input_file = os.path.join(input_folder, file_)
        if (file_[-3:] == "fna"):
            filename = file_[0:-3] + "fa"
            detect = True
            change = True
        elif (file_[-5:] == "fasta"):
            filename = file_[0:-5] + "fa"
            detect = True
            change = True
        elif (file_[-2:] == "fa"):
            filename = file_[0:-2] + "fa"
            detect = True
            change = False
        elif (file_[-3:] == "gff"):
            with open(input_file, "r") as g_f:
                for line in g_f:
                    if line[0] != "#":
                        line = line.strip()
                        line = line.split("\t")
                        break
            if line[0] != file_[:-4]:
                name = line[0]
                os.rename(input_file, os.path.join(input_folder, name + ".gff"))
        elif (file_[-3:] == "gbk"):
            with open("/".join([input_folder, file_]), "r") as g_f:
                for line in g_f:
                    if line[0:7] == "VERSION":
                        data = line[12:].split(" ")
                        break
            if data[0] != file_[:-4]:
                name = data[0]
                os.rename(input_file, os.path.join(input_folder, name + ".gbk"))
        if detect:
            detect = False
            if change:
                os.rename(input_file, os.path.join(input_folder, filename))
                change = False
            SeqEditer().modify_header(os.path.join(input_folder, filename))
