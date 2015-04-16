#!/usr/bin/python

import os        
import sys
import math
import csv
import itertools
from annogesiclib.gff3 import Gff3Parser


def get_gene_name(data):
    if ("blast_hit" in data.attributes.keys()):
        return "NA"
    else:
        name = []
        for key, value in data.attributes.items():
            if "blast_hit" in key:
                if value.split(",")[2].split(":")[1] not in name:
                    name.append(value.split(",")[2].split(":")[1])
        return name

def print_intersection(datas, keys, num_srna, gff_name, type_):
    num = 0
    datas_merge = []
    datas_cover = []
    if type_ == "total":
        out = open(gff_name, "w")
        out.write("##gff-version 3\n")
    for data in datas[keys[0]]:
        check_same = []
        for key in keys[1:]:
            if data in datas[key]:
                check_same.append("True")
        if len(check_same) == (len(keys) - 1):
            if len(keys) <= 5:
                if type_ == "total":
                    out.write(data.info + "\t" + "\n")
            if "best_avg_coverage" in data.attributes.keys():
                datas_merge.append({"data": data, "wig": data.attributes["best_avg_coverage"]})
            num += 1
    if len(keys) == 5:
        datas_sort = sorted(datas_merge, key=lambda k: float(k['wig']), reverse=True)
        for data in datas_sort:
            if type_ == "total":
                out.write(data["data"].info + "\n")
    print("\t{0} = {1}({2})".format(" and ".join(keys), str(num), str(float(num)/float(num_srna))))

def initiate(key, key_list, class_name, class_num, index, out, content):
    if key in key_list:
        class_num += 1
        index[class_name] = class_num
        out.write(str(class_num) + content + "\n")
    return class_num

def import_class(class_num, datas_srna, datas, index, num_srna, strain, type_, utr):
    for num in range(1, class_num + 1):
        datas_srna["class_" + str(num)] = []
    for data in datas[strain]:
        detect = False
        if (data.source == type_) or (type_ == "total"):
            if type_ == "UTR_derived":
                if utr in data.attributes["UTR_type"]:
                    detect = True
            else:
                detect = True
            if detect:
                num_srna += 1
                if "2d_energy" in data.attributes.keys():
                    if float(data.attributes["2d_energy"]) < args.energy:
                        datas_srna["class_" + str(index["2d_energy"])].append(data)
                if "with_TSS" in data.attributes.keys():
                    if data.attributes["with_TSS"] != "NA":
                        datas_srna["class_" + str(index["with_TSS"])].append(data)
                    elif ((type_ == "UTR_derived") or (type_ == "total")) and \
                         (data.source == "UTR_derived"):
                        if (data.attributes["with_cleavage"] != "NA") and \
                           (("3utr" in data.attributes["UTR_type"]) or \
                            ("interCDS" in data.attributes["UTR_type"])):
                            datas_srna["class_" + str(index["with_TSS"])].append(data)
                if "nr_hit" in data.attributes.keys():
                    if ((data.attributes["nr_hit"] != "NA") and \
                       (int(data.attributes["nr_hit"]) <= args.hit_nr_num)) or \
                       (data.attributes["nr_hit"] == "NA"):
                        datas_srna["class_" + str(index["nr_no_hit"])].append(data)
                if "sORF" in data.attributes.keys():
                    if (data.attributes["sORF"] == "NA"):
                        datas_srna["class_" + str(index["sORF"])].append(data)
                if "sRNA_hit" in data.attributes.keys():
                    if data.attributes["sRNA_hit"] != "NA":
                        datas_srna["class_" + str(index["sRNA_hit"])].append(data)
                    else:
                        datas_srna["class_" + str(index["sRNA_no_hit"])].append(data)
    return num_srna

def import_data(class_num, datas_srna, datas, index, num_srna, strain, utr, inter):
    if utr:
        datas_srna["5'UTR_derived"] = {}
        num_srna["5'UTR_derived"] = import_class(class_num, datas_srna["5'UTR_derived"], datas, 
                                  index, num_srna["5'UTR_derived"], strain, "UTR_derived", "5utr")
        datas_srna["3'UTR_derived"] = {}
        num_srna["3'UTR_derived"] = import_class(class_num, datas_srna["3'UTR_derived"], datas,
                                  index, num_srna["3'UTR_derived"], strain, "UTR_derived", "3utr")
        datas_srna["interCDS"] = {}
        num_srna["interCDS"] = import_class(class_num, datas_srna["interCDS"], datas,
                               index, num_srna["interCDS"], strain, "UTR_derived", "interCDS")
    if inter:
        datas_srna["intergenic"] = {}
        num_srna["intergenic"] = import_class(class_num, datas_srna["intergenic"], datas, 
                                         index, num_srna["intergenic"], strain, "intergenic", None)
    if utr:
        datas_srna["total"] = {}
        num_srna["total"] = import_class(class_num, datas_srna["total"], datas,
                                         index, num_srna["total"], strain, "total", None)

def sort_keys(keys):
    nums = []
    final_keys = []
    for key in keys:
        nums.append(int(key.split("_")[1]))
    nums = sorted(nums)
    for num in nums:
        final_keys.append("_".join(["class", str(num)]))
    return final_keys

def print_stat_title(checks, out_stat, strain, srna_datas, index, energy, hit_nr_num):
    class_num = 0
    if checks["first"]:
        checks["first"] = False
    else:
        out_stat.write("\n")
    if checks["limit"] is True:
        return class_num
    elif strain == "all":
        out_stat.write("All strains:\n")
        if len(strains) <= 2:
            checks["limit"] = True
    else:
        out_stat.write(strain + ":\n")
    class_num = initiate("2d_energy", srna_datas[strain][0].attributes.keys(),
                         "2d_energy", class_num, index, out_stat,
                         " - free energy change of secondary structure below to " + str(energy))
    name = " ".join([" - sRNA candidates start with TSS", \
                     "(3'UTR derived sRNA also includes the sRNA candidates which start with processing site.", \
                     "interCDS includes start and end with processing site.)"])
    class_num = initiate("with_TSS", srna_datas[strain][0].attributes.keys(), "with_TSS",
                         class_num, index, out_stat, name)
    class_num = initiate("nr_hit", srna_datas[strain][0].attributes.keys(),
                         "nr_no_hit", class_num, index, out_stat,
                         "".join([" - blast can not find the homology from nr database (the cutoff is ",
                                  str(hit_nr_num), ")."]))
    class_num = initiate("sORF", srna_datas[strain][0].attributes.keys(),
                         "sORF", class_num, index, out_stat, " - have no confliction of sORF candidates.")
    class_num = initiate("sRNA_hit", srna_datas[strain][0].attributes.keys(),
                         "sRNA_no_hit", class_num, index, out_stat,
                         " - blast can not find the homology from sRNA database.")
    class_num = initiate("sRNA_hit", srna_datas[strain][0].attributes.keys(),
                         "sRNA_hit", class_num, index, out_stat,
                         " - blast can find the homology from sRNA database.")
    return class_num

def read_file(sRNA_file, checks, strains):
    srna_datas = {}
    srna_datas["all"] = []
    strains.append("all")
    pre_seq_id = ""
    for entry in Gff3Parser().entries(open(sRNA_file)):
        if entry.source == "UTR_derived":
            checks["utr"] = True
        elif entry.source == "intergenic":
            checks["inter"] = True
        if entry.seq_id != pre_seq_id:
            srna_datas[entry.seq_id] = []
            strains.append(entry.seq_id)
            pre_seq_id = entry.seq_id
        srna_datas[entry.seq_id].append(entry)
        srna_datas["all"].append(entry)
    for strain in datas.keys():
        srna_datas[strain] = sorted(srna_datas[strain], key=lambda k: (k.seq_id, k.start))
    return srna_datas

def classify_sRNA(sRNA_file, out_folder, energy, hit_nr_num, out_stat_file):
    strains = []
    checks = {"limit": False, "first": True, "utr": False, "inter": False}
    srna_datas = read_file(sRNA_file, checks, strains)
    out_stat = open(out_stat_file, "w")
    for strain in strains:
        class_num = 0
        if checks["utr"]:
            num_srna = {"total": 0, "intergenic": 0, "5'UTR_derived": 0, 
                        "3'UTR_derived": 0, "interCDS": 0}
        else:
            num_srna = {"intergenic": 0}
        index = {}
        srna_class = {}
        class_num = print_stat_title(checks, out_stat, strain, srna_datas, 
                                index, energy, hit_nr_num)
        if checks["limit"] is True:
            break
        import_data(class_num, srna_class, srna_datas, index, 
                    num_srna, strain, checks["utr"], checks["inter"])
        for type_, srna in num_srna.items():
            print("sRNA type - {0}:".format(type_))
            print("\ttotal sRNA candidates = {0}".format(srna))
            for num in range(1, class_num + 1):
                print("\tclass {0} = {1}({2})".format(
                         str(num), str(len(srna_class[type_]["class_" + str(num)])), 
                         str(float(len(srna_class[type_]["class_" + str(num)])) / float(srna))))
                if type_ == "total":
                    out = open(os.path.join(out_folder, "_".join(["class", str(num), strain + ".gff"])), "w")
                    out.write("##gff-version 3\n")
                    for data in srna_class[type_]["_".join(["class", str(num)])]:
                        out.write(data.info + "\n")
            if class_num >= 2:
                for comb in range(2 , class_num):
                    for keys in itertools.combinations(srna_class[type_].keys(), comb):
                            if (("class_" + str(index["sRNA_hit"])) in keys) and \
                               (("class_" + str(index["sRNA_no_hit"])) in keys):
                                continue
                            else:
                                keys = sort_keys(list(keys))
                                gff_name = os.path.join(out_folder, "_".join(sorted(list(keys)) + [strain]) + ".gff")
                                print_intersection(srna_class[type_], keys, srna, gff_name, type_)
