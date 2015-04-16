#!/usr/bin/python

import os        
import sys
import csv
from annogesiclib.gff3 import Gff3Parser


def line_to_dict(hit, strain, e_value):
    return {"hit": hit, "strain": strain,"e_value": e_value}

def get_proteins(datas, nums, checks, proteins, blast_f):
    for data in datas:
        if (nums["index"] % 4 == 0) and (nums["index"] != 0):
            if "[" in data:
                data1 = data.split("[")
                data2 = data1[1].split("]")
                strain = data2[0].strip()
            else:
                data1 = data.split("Length")
                strain = "NA"
            name = data1[0].strip()
            tag = datas[nums["index"] - 1]
            if ("hypothetical" in name) or \
               ("unknown" in name):
                nums["hypo"] += 1
            if checks["detect"] is False:
                for line in blast_f:
                    line = line.strip()
                    if "Expect =" in line:
                        e_value = line.split(",")[1].split(" ")[-1]
                        checks["detect"] = True
                        checks["print"] = True
                        proteins.append({"name": name, "strain": strain, "e": e_value, "tag": tag})
                        break
            else:
                proteins.append({"name": name, "strain": strain, "e": e_value, "tag": tag})
        nums["index"] += 1

def detect_nr(line, blast_f, out_t, blasts, prefix):
    checks = {"print": False, "detect": False}
    if line.startswith(">"):
        info = line.replace(">", "")
        for line in blast_f:
            line = line.strip()
            if len(line) != 0:
                info = " ".join([info, line])
            else:
                break
        datas = info.split("|")
        nums = {"index": 0, "hypo": 0}
        proteins = []
        get_proteins(datas, nums, checks, proteins, blast_f)
        if checks["print"] is True:
            if float(nums["hypo"]) / float(len(proteins)) <= 0.5:
                blasts["hit_num"] += 1
                if blasts["hit_num"] <= 3:
                    protein_names = {}
                    for protein in proteins:
                        if ("hypothetical" not in protein["name"]) and \
                           ("unknown" not in protein["name"]):
                            if (protein["name"] not in protein_names.keys()):
                                protein_names[protein["name"]] = []
                            protein_names[protein["name"]].append(protein["tag"])
                            blasts["blast"] = True
                    if len(protein_names) != 0:
                        for key, value in protein_names.items():
                            out_t.write("{0}\t{1}\t{2}\t{3}\n".format(prefix, key, ",".join(value), protein["e"]))

def detect_sRNA(line, blast_f, out_t, blasts, prefix):
    print_ = False
    blasts["name"] = ""
    if line[0] == ">":
        blasts["name"] = line[1:]
        blasts["hit_num"] += 1
        for line in blast_f:
            if "Expect =" in line:
                e_value = line.split(" ")[-1].strip()
                print_ = True
                break
    if print_:
        out_t.write("{0}\t{1}\t{2}\n".format(
                    prefix, blasts["name"], e_value))
        blasts["blast"] = True
    return blasts

def read_gff(sRNA_file, data_type):
    srnas = []
    srna_f = open(sRNA_file, "r")
    print(data_type)
    for entry in Gff3Parser().entries(srna_f):
        attributes = {}
        for key, value in entry.attributes.items():
            if data_type == "sRNA" and \
               "sRNA_hit" not in key:
                attributes[key] = value
            elif data_type == "nr" and \
                 "nr_hit" not in key:
                attributes[key] = value
        entry.attributes = attributes
        attribute_string = ";".join(
            ["=".join(items) for items in entry.attributes.items()])
        entry.info = "\t".join([entry.info_without_attributes, attribute_string])
        srnas.append(entry)
    srnas = sorted(srnas, key=lambda k: (k.seq_id, k.start))
    return srnas

def print_file(database, out_f, info, srna_hit, nr_hit):
    if database == "sRNA":
        out_f.write("{0};sRNA_hit={1}\n".format(info, srna_hit))
    elif database == "nr":
        out_f.write("{0};nr_hit={1}\n".format(info, nr_hit))

def extract_blast(blast_result, sRNA_file, output_file, output_table, database):
    out_f = open(output_file, "w")
    out_t = open(output_table, "w")
    out_f.write("##gff-version 3\n")
    srnas = read_gff(sRNA_file, database)
    for srna in srnas:
        blasts = {"hit_num": 0, "blast": False, "name": ""}
        names = []
        detect_name = False
        prefix = "\t".join([srna.seq_id, srna.attributes["ID"], 
                           srna.strand, str(srna.start), str(srna.end)])
        with open(blast_result, "r") as blast_f:
            for line in blast_f:
                line = line.strip()
                if line.startswith("Query= "):
                    query = line.split("=")[1].strip()
                    if (query == ("|".join([srna.attributes["ID"], srna.seq_id, \
                              str(srna.start), str(srna.end), srna.strand]))):
                        for line in blast_f:
                            line = line.strip()
                            if line.find("No hits found") != -1:
                                print_file(database, out_f, srna.info, "NA", "NA")
                                out_t.write("{0}\tNA\n".format(prefix))
                                break
                            elif line.find("Sequences producing significant alignments:") != -1:
                                for line in blast_f:
                                    line = line.strip()
                                    if line:
                                        if line.startswith("Effective search space"):
                                            break
                                        if database == "sRNA":
                                            detect_sRNA(line, blast_f, out_t, blasts, prefix)
                                            if (len(blasts["name"]) > 0):
                                                if blasts["name"] not in names:
                                                    names.append(blasts["name"])
                                        elif database == "nr":
                                            detect_nr(line, blast_f, out_t, blasts, prefix)
                                if blasts["blast"] is not True:
                                    out_t.write("{0}\tNA\n".format(prefix))
                                    print_file(database, out_f, srna.info, "NA", "NA")
                                else:
                                    print_file(database, out_f, srna.info,
                                               "&".join(names), blasts["hit_num"])
                                blasts["hit_num"] = 0
                                break
    out_f.close()

def extract_energy(sRNA_file, sec_file, out_file):
    s_f = open(sRNA_file, "r")
    srnas = []
    check = False
    get_length = False
    out = open(out_file, "w")
    out.write("##gff-version 3\n")
    for srna in Gff3Parser().entries(s_f):
        with open (sec_file, "r") as d_f:
            for structure in d_f:
                structure = structure.rstrip('\n')
                if get_length:
                    length = len(structure)
                    get_length = False
                if (structure.startswith(">")):
                    if (("|".join([srna.attributes["ID"], srna.seq_id, str(srna.start),
                         str(srna.end), srna.strand])) == structure[1:]) or \
                       (("|".join([srna.feature, srna.seq_id, str(srna.start),
                         str(srna.end), srna.strand])) == structure[1:]):
                        check = True
                        get_length = True
                if (check is True) and \
                   ((structure[0] == "(") or (structure[0] == ")") or (structure[0] == ".")) and \
                   (structure[-1] == ")"):
                    check = False
                    data=structure.split(" ")
                    if (data[-1].find("(") == -1):
                        energy = float(data[-1][0:-1])
                    else:
                        energy = float(data[-1][1:-1])
                    out.write("{0};2d_energy={1:.4f}\n".format(srna.info, (energy / float(length))))
                    break
    s_f.close()
    d_f.close()
