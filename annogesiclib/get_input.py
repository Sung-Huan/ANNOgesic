import os
import csv
import shutil
from subprocess import call
from annogesiclib.seq_editer import SeqEditer


def wget(input_folder, ftp, files_type, log):
    log.write("\t" + " ".join(["wget", "-cP", input_folder, ftp + "/*" + files_type]) + "\n")
    os.system(" ".join(["wget", "-cP", input_folder, ftp + "/*" + files_type]))
    log.write("Done!\n")

def deal_detect(input_file, file_path, change, input_folder):
    '''deal with the header of fasta file and 
    put the files to corresponding folders'''
    if change:
        shutil.move(input_file, file_path)
        change = False
    SeqEditer().modify_header(file_path)
    with open(os.path.join(file_path)) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                seq_name = line[1:]
    shutil.move(file_path,
                os.path.join(input_folder, seq_name + ".fa"))
    return change, seq_name


def get_file(ftp, input_folder, files_type, log):
    checks = {"detect": False, "change": None}
    filename = None
    files = []
    wget(input_folder, ftp, files_type, log)
    for file_ in os.listdir(input_folder):
        input_file = os.path.join(input_folder, file_)
        if (file_[-3:] == "fna"):
            filename = file_[0:-3] + "fa"
            checks = {"detect": True, "change": True}
        elif (file_[-5:] == "fasta"):
            filename = file_[0:-5] + "fa"
            checks = {"detect": True, "change": True}
        elif (file_[-2:] == "fa"):
            filename = file_[0:-2] + "fa"
            checks = {"detect": True, "change": True}
        elif (file_[-6:] == "fna.gz") and ("_genomic" in file_):
            if ("_cds_from_genomic" in file_) or (
                    "_rna_from_genomic" in file_):
                os.remove(input_file)
            else:
                filename = file_[0:-6] + "fa"
                checks = {"detect": True, "change": True}
                log.write("\tgunzip " + input_file + "\n")
                call(["gunzip", input_file])
                input_file = input_file[:-3]
        elif (file_[-6:] == "gff.gz") or (file_[-3:] == "gff"):
            if ("_genomic" in file_) and (file_[-6:] == "gff.gz"):
                log.write("\tgunzip " + input_file + "\n")
                call(["gunzip", input_file])
                input_file = input_file[:-3]
            fh = open(input_file, "r")
            for row in csv.reader(fh, delimiter='\t'):
                if not row[0].startswith("#"):
                    gff_name = row[0]
                    break
            shutil.move(input_file, os.path.join(input_folder,
                                               gff_name + ".gff"))
            fh.close()
        elif (file_[-3:] == "gbk") or (file_[-7:] == "gbff.gz") or (
                file_[-4:] == "gbff"):
            if (file_[-7:] == "gbff.gz") and ("_genomic" in file_):
                log.write("\tgunzip " + input_file + "\n")
                call(["gunzip", input_file])
                input_file = input_file[:-3]
            with open(input_file, "r") as g_f:
                for line in g_f:
                    line = line.strip()
                    if line.startswith("VERSION"):
                        for data in line.split(" "):
                            if (len(data) != 0) and (data != "VERSION"):
                                break
                        break
            print(os.path.join(input_folder, data + ".gbk"))
            shutil.move(input_file, os.path.join(input_folder, data + ".gbk"))
        if checks["detect"]:
            checks["detect"] = False
            checks["change"], seq_name = deal_detect(
                    input_file, filename, checks["change"], input_folder)
