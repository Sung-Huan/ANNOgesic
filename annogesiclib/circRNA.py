#!/usr/bin/python

import os
import sys
import csv
from annogesiclib.splice_parser import SpliceParser
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.helper import Helper

def get_feature(cds):
    if "locus_tag" in cds.attributes.keys():
        feature = cds.attributes["locus_tag"]
    elif "protein_id" in cds.attributes.keys():
        feature = cds.attributes["protein_id"]
    elif "ID" in cds.attributes.keys():
        strand = Helper().get_strand_name(cds.strand)
        feature = "".join([cds.attributes["ID"], ":",
                  str(cds.start), "-", str(cds.end),
                  "_", strand])
    else:
        strand = Helper().get_strand_name(cds.strand)
        feature = "".join([cds.feature, ":",
                  str(cds.start), "-", str(cds.end),
                  "_", strand])
    return feature

def detect_conflict(gffs, circ, num, out):
    detect = False
    gff = None
    for gff in gffs:
        if (gff.seq_id == circ.strain) and (
            gff.strand == circ.strand):
            if ((gff.start < circ.start) and (
                 gff.end > circ.start) and (
                 gff.end < circ.end)) or (
                (gff.start > circ.start) and (
                 gff.end < circ.end)) or (
                (gff.start > circ.start) and (
                 gff.start < circ.end) and (
                 gff.end > circ.end)):
                detect = True
                break
    if detect:
        feature = get_feature(gff)
        out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(
                  "_".join(["circRNA", str(num)]), circ.strain, circ.strand,
                  circ.start, circ.end, feature, circ.supported_reads,
                  float(circ.supported_reads) / float(circ.start_site_reads),
                  float(circ.supported_reads) / float(circ.end_site_reads)))
    else:
        out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(
                  "_".join(["circRNA", str(num)]), circ.strain, circ.strand,
                  circ.start, circ.end, "NA", circ.supported_reads,
                  float(circ.supported_reads) / float(circ.start_site_reads),
                  float(circ.supported_reads) / float(circ.end_site_reads)))
    return detect

def import_num(support, nums, strain):
    if support not in nums[strain].keys():
        nums[strain][support] = 0
    if support not in nums["all"].keys():
        nums["all"][support] = 0
    nums[strain][support] += 1
    nums["all"][support] += 1

def print_file(nums, stat, strain):
    for key in sorted(nums[strain].keys()):
        stat.write("\tthe number of potential circular RNAs, ")
        stat.write("more than {0} supported it = {1}\n".format(
                   key, nums[strain][key]))

def read_file(input_file, gff_file):
    circs = []
    gffs = []
    ps = SpliceParser()
    high = 0
    splice_fh = open(input_file)
    for entry in ps.parser(splice_fh):
        if entry.supported_reads > high:
            high = entry.supported_reads
        circs.append(entry)
    gff_parser = Gff3Parser()
    for entry in gff_parser.entries(open(gff_file)):
        gffs.append(entry)
    gffs = sorted(gffs, key=lambda k: (k.seq_id, k.start))
    circs = sorted(circs, key=lambda x: (x.strain, x.supported_reads),
                   reverse=True)
    splice_fh.close()
    return circs, gffs, high

def get_circrna(circs, gffs, high, start_ratio, end_ratio, out):
    num_circular = {}
    num_circular["all"] = 0
    num_support = {}
    num_support["all"] = {}
    num_conflict = {}
    num_conflict["all"] = {}
    pre_seq_id = ""
    num = 0
    for circ in circs:
        if pre_seq_id != circ.strain:
            num_support[circ.strain] = {}
            num_conflict[circ.strain] = {}
            num_circular[circ.strain] = 0
        if (circ.situation != "F") and \
           (circ.splice_type == "C"):
            num_circular[circ.strain] += 1
            num_circular["all"] += 1
            detect = detect_conflict(gffs, circ, num, out)
            for support in range(0, high + 5, 5):
                if circ.supported_reads >= int(support):
                    import_num(support, num_support, circ.strain)
                    if detect is False:
                        if (float(circ.supported_reads) / float(
                            circ.start_site_reads) >= start_ratio) and (
                            float(circ.supported_reads) / float(
                            circ.end_site_reads) >= end_ratio):
                            import_num(support, num_conflict, circ.strain)
            num += 1
        pre_seq_id = circ.strain
    return {"circular": num_circular, "support": num_support,
            "conflict": num_conflict}

def detect_circrna(input_file, gff_file, output_file,
                   start_ratio, end_ratio, statistics):
    circs, gffs, high = read_file(input_file, gff_file)
    out = open(output_file, "w")
    out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(
               "ID", "strain", "strand", "start", "end", "annotation_overlap",
               "supported_reads", "supported_reads/reads_at_start",
               "supported_reads/reads_at_end"))
    nums = get_circrna(circs, gffs, high, start_ratio, end_ratio, out)
    stat = open(statistics, "w")
    stat.write("All strains:\n")
    stat.write("\tthe number of all circular RNAs = {0}\n".format(
               nums["circular"]["all"]))
    print_file(nums["support"], stat, "all")
    stat.write("\n\tthe circular RNAs:\n")
    stat.write("\t\twithout conflict with annotation\n")
    stat.write("\t\tsupport reat ratio of starting point is larger than {0}\n".format(
               start_ratio))
    stat.write("\t\tsupport reat ratio of end point is larger than {0}\n".format(
               end_ratio))
    print_file(nums["conflict"], stat, "all")
    if len(nums["circular"]) > 2:
        for strain in nums["circular"].keys():
            if strain != "all":
                stat.write("\n{0}:\n".format(strain))
                stat.write("\tthe number of all circular RNAs = {0}\n".format(
                           nums["circular"][strain]))
                print_file(nums["support"], stat, strain)
                stat.write("\n\tthe circular RNAs:\n")
                stat.write("\t\twithout conflict with annotation\n")
                stat.write("\t\tsupport reat ratio of starting point is larger than {0}\n".format(
                           start_ratio))
                stat.write("\t\tsupport reat ratio of end point is larger than {0}\n".format(
                           end_ratio))
                print_file(nums["conflict"], stat, strain)
    out.close()
    stat.close()
