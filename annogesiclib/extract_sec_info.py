import os
import shutil

def mod_file(input_file, out, indexs):
    with open(input_file) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                out.write(indexs[line] + "\n")
            else:
                out.write(line + "\n")
    out.close()

def extract_info_sec(sec_file, seq_file, index_file):
    out_sec = open(sec_file + "tmp", "w")
    out_seq = open(seq_file + "tmp", "w")
    indexs = {}
    with open(index_file) as hi:
        for line in hi:
            line = line.strip()
            if line.startswith(">"):
                tag = line.split("|")[0]
                indexs[tag] = line
    mod_file(sec_file, out_sec, indexs)
    mod_file(seq_file, out_seq, indexs)
    os.remove(sec_file)
    shutil.move(sec_file + "tmp", sec_file)
    os.remove(seq_file)
    shutil.move(seq_file + "tmp", seq_file)

def modify_header(seq_file, index_file):
    out = open(seq_file, "w")
    with open(index_file) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                tag = line.split("|")[0]
                out.write(tag + "\n")
            else:
                out.write(line + "\n")
