#!/usr/bin/python

import os	
import sys
import csv

def read_file(seq_file, seqs):
    with open(seq_file, "r") as s_h:
        for line in s_h:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]
            else:
                seqs.append({"name": header, "seq": line})

def import_ribo(line, ribos, seq_name):
    if line.startswith("("):
        num = 0
        datas = line.split(" ")
        for data in datas:
            if len(data) != 0:
                num += 1
                if (num == 2):
                    detect = data
                elif (num == 3):
                    e = float(data)
                elif (num == 6):
                    name = data
                elif (num == 7):
                    start = int(data)
                elif (num == 8):
                    end = int(data)
                elif (num == 9):
                    if data == "+":
                        ribos.append({"name": name, "detect": detect, 
                                      "e": e, "seq_name": seq_name, 
                                      "start": start, "end": end})
                    else:
                        ribos.append({"name": name, "detect": detect, 
                                      "e": e, "seq_name": seq_name, 
                                      "start": end, "end": start})

def print_file(ribos, out_t, out_s, seq_name):
    if len(ribos) != 0:
        for rbs in ribos:
            if rbs["detect"] == "!":
                out_t.write("\t".join([seq_name, rbs["name"], str(rbs["e"]), 
                            str(rbs["start"]), str(rbs["end"])]) + "\n")
            for seq in seqs:
                if rbs["seq_name"] == seq["name"]:
                    tags = seq["name"].split("|")
                    out_s.write(">" + "|".join(["|".join(tags[0:-2]),
                                      str(int(tags[-2]) + rbs["start"] - 1),
                                      str(int(tags[-2]) + rbs["end"] - 1)]) + "\n")
                    out_s.write(seq["seq"][(rbs["start"] - 1): (rbs["end"])] + "\n")

def regenerate_seq(align_file, seq_file, out_table, out_seq):
    hit = False
    seqs = []
    out_t = open(out_talbe, "w")
    out_s = open(out_seq, "w")
    read_file(seq_file, seqs)
    with open(align_file, "r") as a_h:
        for line in a_h:
            line = line.strip()
            if line.startswith("#"):
                continue
            else:
                if line.startswith("Query:"):
                    datas = line.split("[")[0]
                    ribos = []
                    seq_name = datas.strip().split(" ")[-1]
                    hit = False
                elif line.startswith("Hit scores:"):
                    hit = True
                elif hit:
                    import_ribo(line, ribos, seq_name)
                if line.startswith("Hit alignments:"):
                    hit = False
                    print_file(ribos, out_t, out_s, seq_name)

def compare_first_result(ribos, firsts, seq_name, out, extras):
    if len(ribos) != 0:
        for rbs in ribos:
            check = False
            same = False
            info = ""
            best = ""
            if rbs["detect"] == "!":
                for first in firsts:
                    if first["seq_name"] == "|".join(rbs["seq_name"].split("|")[0:4]):
                        same = True
                        if (first["acc"] == rbs["name"]) and \
                           (first["e"] > rbs["e"]):
                            first["print"] = True
                            first["e"] = rbs["e"]
                            out.write("\t".join([seq_name, rbs["name"], str(rbs["e"]), 
                                      str(rbs["start"]), str(rbs["end"])]) + "\n")
                            if len(info) != 0:
                                info = ""
                        elif (first["acc"] != rbs["name"]):
                            info = "\t".join([seq_name, rbs["name"], str(rbs["e"]), 
                                   str(rbs["start"]), str(rbs["end"])])
                if len(info) != 0:
                    out.write(info + "\n")
                if same is False:
                    if (len(extras) == 0):
                        extras.append(rbs)
                    else:
                        for extra in extras:
                            if ("|".join(extra["seq_name"].split("|")[0:4]) == \
                                "|".join(rbs["seq_name"].split("|")[0:4])):
                                check = True
                                if (extra["name"] == rbs["name"]):
                                    if extra["e"] > rbs["e"]:
                                        extra["seq_name"] = rbs["seq_name"]
                                        extra["e"] = rbs["e"]
                                        extra["start"] = rbs["start"]
                                        extra["end"] = rbs["end"]
                                else:
                                    extras.append(rbs)
                        if check is not True:
                            extras.append(rbs)

def reextract_rbs(align_file, first_file, output_file):
    hit = False
    seqs = []
    extras = []
    out = open(output_file, "w")
    fh = open(first_file, "r");
    firsts = []
    pre_seq_name = ""
    for row in csv.reader(fh, delimiter="\t"):
        datas = row[0].split("|")
        firsts.append({"seq_name": "|".join(datas[0:4]), "acc": row[1], "e": float(row[2]), 
                       "start": int(row[3]), "end": int(row[4]), "print": False, 
                       "pre_start": int(datas[4]), "pre_end": int(datas[5])})
    with open(align_file, "r") as a_h:
        for line in a_h:
            line = line.strip()
            if line.startswith("#"):
                continue
            else:
                if line.startswith("Query:"):
                    datas = line.split("[")[0]
                    seq_name = datas.strip().split(" ")[-1]
                    ribos = []
                    hit = False
                elif line.startswith("Hit scores:"):
                    hit = True
                elif hit:
                    import_ribo(line, ribos, seq_name)
                if line.startswith("Hit alignments:"):
                    hit = False
                    compare_first_result(ribos, firsts, seq_name, out, extras)
        if len(extras) != 0:
            for extra in extras:
                out.write("\t".join([extra["seq_name"], extra["name"], str(extra["e"]), str(extra["start"]), str(extra["end"])]) + "\n")
    for first in firsts:
        if first["print"] is False:
            out.write("\t".join(["|".join([first["seq_name"], str(first["pre_start"]), str(first["pre_end"])]), 
                                 first["acc"], str(first["e"]), str(first["start"]), str(first["end"])]) + "\n")
