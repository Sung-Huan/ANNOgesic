#!/usr/bin/python

import os
import sys
import csv
from annogesiclib.gff3 import Gff3Parser


def create_dict(nums, strain, utr_detect):
    nums[strain] = {}
    if utr_detect:
        types = ["all", "5'UTR_derived", "3'UTR_derived", "interCDS", "intergenic"]
    else:
        types = ["all"]
    for type_ in types:
        nums[strain][type_] = {}
        for feature in ["TSS", "sRNA", "all", "both"]:
            nums[strain][type_][feature] = 0

def plus_data(nums, strain, sorf_types, features, utr_detect):
    for sorf_type in sorf_types:
        if ((not utr_detect) and \
           (sorf_type != "intergenic")) or \
           (utr_detect):
            for feature in features:
                nums[strain][sorf_type][feature] += 1

def print_stat(nums, strain, out, utr_detect):
    out.write(strain + ":\n")
    if utr_detect:
        out.write("\ttotal sORF in this strain are {}0\n".format(nums[strain]["all"]["all"]))
    for type_, features in nums[strain].items():
        out.write("\ttotal sORF of {0} sORF candidates are{1}".format(type_, nums[strain][type_]["all"]))
        out.write("(for this strain - {0})\n".format(float(nums[strain][type_]["all"]) / \
                                                     float(nums[strain]["all"]["all"])))
        for feature, num in features.items():
            if feature == "TSS":
                out.write("\t\ttotal sORF which start from TSS are {0}".format(num))
                out.write("(for strain {0};".format(float(num) / float(nums[strain]["all"]["all"])))
                out.write("for {0} - {1})\n".format(type_, (float(num) / float(nums[strain][type_]["all"]))))
            elif feature == "sRNA":
                out.write("\t\ttotal sORF without overlap with sRNA candidates are {0}".format(num))
                out.write("(for strain {0};".format(float(num) / float(nums[strain]["all"]["all"])))
                out.write("for {0} - {1})\n".format(type_, float(num) / float(nums[strain][type_]["all"])))
            elif feature == "both":
                out.write("\t\ttotal sORF which start from TSS and without overlap with sRNA candidates are {0}".format(num))
                out.write("(for strain {0}; ".format(float(num) / float(nums[strain]["all"]["all"])))
                out.write("for {0} - {1})\n".format(type_, float(num) / float(nums[strain][type_]["all"])))

def stat(sORF_gff, stat_file, utr_detect):
    nums = {}
    create_dict(nums, "total", utr_detect)
    sorfs = []
    for entry in Gff3Parser().entries(open(sORF_gff)):
        sorfs.append(entry)
    sorfs = sorted(sorfs, key=lambda k: (k.seq_id, k.start))
    strain = ""
    for sorf in sorfs:
        if strain != sorf.seq_id:
            create_dict(nums, sorf.seq_id, utr_detect)
            strain = sorf.seq_id
        if sorf.source == "intergenic":
            sorf_type = "intergenic"
        else:
            if "5utr" in sorf.attributes["UTR_type"]:
                sorf_type = "5'UTR_derived"
            elif "3utr" in sorf.attributes["UTR_type"]:
                sorf_type = "3'UTR_derived"
            elif "interCDS" in sorf.attributes["UTR_type"]:
                sorf_type = "interCDS"
        if (sorf.attributes["with_TSS"] != "NA") and \
           (sorf.attributes["sRNA"] == "NA"):
            plus_data(nums, "total", [sorf_type, "all"], ["all", "TSS", "sRNA", "both"], utr_detect)
            plus_data(nums, strain, [sorf_type, "all"], ["all", "TSS", "sRNA", "both"], utr_detect)
        elif sorf.attributes["with_TSS"] != "NA":
            plus_data(nums, "total", [sorf_type, "all"], ["all", "TSS"], utr_detect)
            plus_data(nums, strain, [sorf_type, "all"], ["all", "TSS"], utr_detect)
        elif sorf.attributes["sRNA"] != "NA":
            plus_data(nums, "total", [sorf_type, "all"], ["all", "sRNA"], utr_detect)
            plus_data(nums, strain, [sorf_type, "all"], ["all", "sRNA"], utr_detect)
        else:
            plus_data(nums, "total", [sorf_type, "all"], ["all"], utr_detect)
            plus_data(nums, strain, [sorf_type, "all"], ["all"], utr_detect)
    out = open(stat_file, "w")
    if len(nums) <= 2:
        for strain, types in nums.items():
            if strain != "total":
                print_stat(nums, strain, out, utr_detect)
    else:
        for strain, types in nums.items():
            print_stat(nums, strain, out, utr_detect)
