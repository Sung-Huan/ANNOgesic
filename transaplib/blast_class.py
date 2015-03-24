#!/usr/bin/python

import os	
import sys
import csv
from transaplib.gff3 import Gff3Parser


def plus_data(id_name, name, e_value, srna, strain):
    if name not in srna[strain].keys():
        srna[strain][name] = {"id_name": [], "num": 1, "e": []}
        srna[strain][name]["id_name"].append(id_name)
        srna[strain][name]["e"].append(e_value)
    else:
        srna[strain][name]["id_name"].append(id_name)
        srna[strain][name]["e"].append(e_value)
        srna[strain][name]["num"] += 1

def import_data(srnas, row):
    srnas.append({"strain": row[0], "name": row[1], "strand": row[2],
                  "start": row[3], "end": row[4], "blast": row[5],
                  "ori_strain": row[6], "ID": row[7], "e": row[8]})

def read_file(sRNA_file):
    fh = open(sRNA_file, "r")
    first = True
    srna_names = []
    srnas = []
    for row in csv.reader(fh, delimiter="\t"):
        if row[5] != "NA":
            if first:
                first = False
            else:
                if (row[0] == pre[0]) and \
                   (row[1] == pre[1]) and \
                   (row[2] == pre[2]) and \
                   (row[3] == pre[3]) and \
                   (row[4] == pre[4]):
                    srna_names.append(pre)
                    if row[5] != pre[5]:
                        if row[1] in repeats.keys():
                            if row[5] not in repeats[row[1]]:
                                repeats[row[1]].append(row[5])
                        else:
                            repeats[row[1]] = []
                            repeats[row[1]].append(pre[5])
                            repeats[row[1]].append(row[5])
                else:
                    srna_names.append(pre)
                    srna_names = sorted(srna_names, key = lambda x: (x[5]))
                    pre_name = ""
                    for srna_name in srna_names:
                        if pre_name != srna_name[5]:
                            import_data(srnas, srna_name)
                            pre_name = srna_name[5]
                    srna_names = []
            pre = row
    if row[5] != "NA":
        srna_names.append(pre)
        srna_names = sorted(srna_names, key = lambda x: (x[5]))
        pre_name = ""
        for srna_name in srna_names:
            if pre_name != srna_name[5]:
                import_data(srnas, srna_name)
    srnas = sorted(srnas, key = lambda x: (x["strain"]))
    return srnas

def Blast_class(sRNA_file, out_file):
    nums = {}
    nums["total"] = {}
    pre_strain = ""
    repeats = {}
    srnas = read_file(sRNA_file, repeats)
    out = open(out_file, "w")
    for srna in srnas:
        if srna["strain"] != pre_strain:
            nums[srna["strain"]] = {}
            pre_strain = srna["strain"]
        if srna["blast"] not in nums[srna["strain"]].keys(): 
            nums[srna["strain"]][srna["blast"]] = {"num": 1, "ID": [], "e": []}
            nums["total"][srna["blast"]] = {"num": 1, "ID": [], "e": []}
        else:
            nums[srna["strain"]][srna["blast"]]["num"] += 1
            nums["total"][srna["blast"]]["num"] += 1
        nums[srna["strain"]][srna["blast"]]["ID"].append(srna["ID"])
        nums[srna["strain"]][srna["blast"]]["e"].append(srna["e"])
        nums["total"][srna["blast"]]["ID"].append(srna["ID"])
        nums["total"][srna["blast"]]["e"].append(srna["e"])
    if len(nums) > 1:
        out.write("All strain:\n")
        out.write("sRNA_name\tamount\n")
        for blast, details in nums["total"].items():
            out.write("%s\t%s\n" % (blast, details["num"]))
    else:
        for strain, srna_name in srna.items():
            out.write(strain + ":\n")
            out.write("sRNA_name\tamount\n")
            for blast, details in srna_name.items():
                srna_name = blast.split("_")
                out.write("%s\t%s\n" %(blast, details["num"]))
    out.write("\nrepeat counting:\n")
    for name, srna_name in repeats.items():
        out.write(name + ":" + ";".join(srna_name) + "\n")
