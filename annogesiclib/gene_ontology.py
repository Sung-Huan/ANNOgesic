#!/usr/bin/python

import os    
import sys
import random
import csv
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from annogesiclib.gff3 import Gff3Parser


def _import_uniprot_data(entry, name_list, feature):
    ref_name = entry.attributes[feature]
    if ref_name not in name_list:
        name_list.add(ref_name)

def retrieve_uniprot(database_file, gff_file, out_file):
    name_list = set()
    gffs = []
    out = open(out_file, "w")
    out.write("\t".join(["strain", "strand", "start", "end",
                         "protein_id", "Go_term"]) + "\n")
    for entry in Gff3Parser().entries(open(gff_file)):
        if entry.feature == "CDS":
            if ("Name" in entry.attributes.keys()) and \
               ("protein_id" in entry.attributes.keys()):
                if entry.attributes["Name"] == entry.attributes["protein_id"]:
                    _import_uniprot_data(entry, name_list, "Name")
                else:
                    _import_uniprot_data(entry, name_list, "Name")
                    _import_uniprot_data(entry, name_list, "protein_id")
            elif ("Name" in entry.attributes.keys()):
                _import_uniprot_data(entry, name_list, "Name")
            elif ("protein_id" in entry.attributes.keys()):
                _import_uniprot_data(entry, name_list, "protein_id")
            else:
                ref_name = ''
            gffs.append(entry)
    idmapping = open(database_file,"r")
    for uni_id in idmapping:
        uni_line = uni_id.rstrip("\n")
        uni_lines = uni_line.split("\t")
        if uni_lines[3] in name_list:
            for gff in gffs:
                if (uni_lines[3] == gff.attributes["Name"]) or \
                   (uni_lines[3] == gff.attributes["protein_id"]):
                    out.write("\t".join([gff.seq_id, gff.strand, str(gff.start),
                                   str(gff.end), uni_lines[3], uni_lines[6]]) + "\n")

def _get_label(class_, classes, percent, num):
    if percent >= 3:
        name = ""
        if len(class_) > 30:
            num_word = 0
            first = True
            words = class_.split(" ")
            for word in words:
                num_word = num_word + len(word)
                if num_word > 30:
                    name = name + "\n" + word
                    num_word = len(word)
                else:
                    if first:
                        first = False
                        name = word
                    else:
                        name = " ".join([name, word])
        else:
            name = class_
        classes.append(("{0} {1} ({2}%)").format(name, num, round(percent, 2)))

def _random_color(count, colors, color, color_elements):
    while True:
        color_index = random.randint(0, 5)
        if count % 3 == 0:
            color = color[:4] + color_elements[color_index]
        elif count % 3 == 1:
            color = color_elements[color_index] + color[2:]
        elif count % 3 == 2:
            color = color[:2] + color_elements[color_index] + color[4:]
        if (color != "000000") and (color != "FFFFFF") and \
           ("#" + color not in colors):
            colors.append("#" + color)
            break
    return color

def _plot(total_nums, strain, filename, total, out_folder):
    color_elements = ["FF", "CC", "99", "66", "33", "00"]
    classes = []
    nums = []
    explode = []
    colors = []
    color = "FFFFFF"
    ini_color = ["FFFF99", "99FFCC", "FFCCCC", "CCFF66", "FFCC33", "CC99FF"]
    fig = plt.figure(figsize=(12, 6))
    plt.subplot(121)
    count = 0
    sort_total_nums = sorted(total_nums.items(), key = lambda x: (x[1]), reverse=True)
    for total_num in sort_total_nums:
        class_ = total_num[0]
        num = total_num[1]
        if class_ != "total":
            percent = (float(num) / float(total)) * 100
            _get_label(class_, classes, percent, num)
            nums.append(num)
            explode.append(0)
            if count <= 5:
                colors.append("#" + ini_color[count])
            else:
                color = _random_color(count, colors, color, color_elements)
            count += 1
    plt.pie(nums, explode=explode, colors=colors, startangle=90)
    plt.axis('equal')
    plt.legend(classes, bbox_to_anchor=(1.7, 0.5), loc=10, shadow=True, prop={'size':15})
    plt.savefig(os.path.join(out_folder, "_".join([strain, filename + ".png"])))
    plt.clf()

def _import_obo(filename, obos):
    start = False
    with open(filename, "r") as o_h:
        for line in o_h:
            line = line.strip()
            if line == "[Term]":
                obo = {}
                start = True
            elif start:
                if len(line) == 0:
                    obos.append(obo.copy())
                    start = False
                else:
                    datas = line.split(": ")
                    if datas[0] == "is_a":
                        if "is_a" not in obo.keys():
                            obo["is_a"] = []
                        obo["is_a"].append(datas[1].strip())
                    else:
                        obo[datas[0]] = datas[1].strip()

def _import_class_name(slim_obo, classes, strain):
    if slim_obo["name"] not in classes[strain][slim_obo["namespace"]]:
        classes[strain][slim_obo["namespace"]][slim_obo["name"]] = 0
    classes[strain][slim_obo["namespace"]][slim_obo["name"]] += 1

def _import_total_num(slim_obo, total_nums, strain):
    total_nums[strain][slim_obo["namespace"]] += 1
    total_nums[strain]["total"] += 1

def _print_file(classes, total_nums, out_folder, stat):
    out_stat = open(stat, "w")
    printed = True
    for strain, datas in classes.items():
        if (strain == "All_strain") and len(classes) <= 2:
            printed = False
        if printed:
            _plot(total_nums[strain], strain, "three_roots", 
                  total_nums[strain]["total"], out_folder)
            out_stat.write("{0}:\n".format(strain))
            for origin, types in datas.items():
                _plot(types, strain, origin, total_nums[strain][origin], out_folder)
                out_stat.write("\t{0}: {1}(percentage in total: {2})\n".format(
                               origin, total_nums[strain][origin],
                               float(total_nums[strain][origin]) / float(total_nums[strain]["total"])))
                for type_, num in types.items():
                    out_stat.write("\t\t{0}: {1}(percentage in {2}: {3})\n".format(
                                   type_, num, origin, float(num) / float(total_nums[strain][origin])))
        else:
            printed = True

def _initiate_dict(classes, total_nums, index):
    classes[index] = {"biological_process": {},
                      "cellular_component": {},
                      "molecular_function": {}}
    total_nums[index] = {"biological_process": 0,
                         "cellular_component": 0,
                         "molecular_function": 0,
                         "total": 0}

def _compare_go_slim(gos, term_obos, slim_obos, classes, total_nums):
    detect = False
    for strain, go_ids in gos.items():
        for go_id in go_ids:
            target_terms = [go_id]
            for target_term in target_terms:
                for term_obo in term_obos:
                    if target_term == term_obo["id"]:
                        if "is_a" in term_obo.keys():
                            for is_a in term_obo["is_a"]:
                                go_a = is_a.split(" ! ")
                                if (go_a[1] != "biological_process") and \
                                   (go_a[1] != "cellular_component") and \
                                   (go_a[1] != "molecular_function"):
                                    target_terms.append(go_a[0])
                        elif ("is_obsolete" in term_obo.keys()):
                            if term_obo["is_obsolete"] == "true":
                                break
                        for slim_obo in slim_obos:
                            for target_term in target_terms:
                                if target_term == slim_obo["id"]:
                                    detect = True
                                    _import_class_name(slim_obo, classes, strain)
                                    _import_class_name(slim_obo, classes, "All_strain")
                                    _import_total_num(slim_obo, total_nums, strain)
                                    _import_total_num(slim_obo, total_nums, "All_strain")
                        break
                if detect:
                    detect = False
                    break

def map2goslim(slim_file, term_file, go_table, stat, out_folder):
    gos = {}
    term_obos = []
    slim_obos = []
    classes = {}
    total_nums = {}
    _initiate_dict(classes, total_nums, "All_strain")
    detect = False
    pre_strain = ""
    g_h = open(go_table, "r")
    print("load go tabel")
    for row in csv.reader(g_h, delimiter="\t"):
        if row[0] != "strain":
            if row[0] != pre_strain:
                gos[row[0]] = []
                _initiate_dict(classes, total_nums, row[0])
            go_terms = row[-1].split("; ")
            gos[row[0]] = gos[row[0]] + go_terms
            pre_strain = row[0]
    print("load obo file")
    _import_obo(term_file, term_obos)
    _import_obo(slim_file, slim_obos)
    print("start mapping")
    _compare_go_slim(gos, term_obos, slim_obos, classes, total_nums)
    print("statistics and ploting...")
    _print_file(classes, total_nums, out_folder, stat)
