#!/usr/bin/python

import os        
import sys
import csv
import transaplib.splice_parser
from transaplib.gff3 import Gff3Parser


def get_feature(cds):
    if "locus_tag" in cds.attributes.keys():
        feature = cds.attributes["locus_tag"]
    elif "protein_id" in cds.attributes.keys():
        feature = cds.attributes["protein_id"]
    elif "ID" in cds.attributes.keys():
        feature = cds.attributes["ID"] + "_" + \
                  str(cds.start) + "-" + str(cds.end) + \
                  "_" + cds.strand
    else:
        feature = cds.feature + "_" + \
                  str(cds.start) + "-" + str(cds.end) + \
                  "_" + cds.strand
    return feature

def detect_conflict(gffs, circ, num, out):
    detect = False
    for gff in gffs:
        if (gff.seq_id == circ.strain) and \
           (gff.strand == circ.strand):
            if ((gff.start < circ.start) and (gff.end > circ.start) and (gff.end < circ.end)) or \
               ((gff.start > circ.start) and (gff.end < circ.end)) or \
               ((gff.start > circ.start) and (gff.start < circ.end) and (gff.end > circ.end)):
                detect = True
                break
    if detect:
        feature = get_feature(gff)
        out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % \
                 ("circRNA_" + str(num), circ.strain, circ.strand, circ.start,
                  circ.end, feature, circ.supported_reads,
                  str(float(circ.supported_reads) / float(circ.start_site_reads)),
                  str(float(circ.supported_reads) / float(circ.end_site_reads))))
    else:
        out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % \
                 ("circRNA_" + str(num), circ.strain, circ.strand, circ.start,
                  circ.end, "NA", circ.supported_reads,
                  str(float(circ.supported_reads) / float(circ.start_site_reads)),
                  str(float(circ.supported_reads) / float(circ.end_site_reads))))
    return detect

def import_num(support, nums, strain):
    if support not in nums[strain].keys():
        nums[strain][support] = 0
    if support not in nums["all"].keys():
        nums["all"][support] = 0
    nums[strain][support] += 1
    nums["all"][support] += 1

def print_num(nums, stat, strain):
    for key in sorted(nums[strain].keys()):
        stat.write("\tthe number of potential circular RNAs, more than %s supported it = %s\n" %
                   (str(key), str(nums[strain][key])))

def print_stat(statistics, num_circular, num_support, num_conflict, start_ratio, end_ratio):
    stat = open(args.statistics, "w")
    stat.write("All strains:\n")
    stat.write("\tthe number of all circular RNAs = %s\n" % (num_circular["all"]))
    print_num(num_support, stat, "all")
    stat.write("\n\tthe circular RNAs:\n")
    stat.write("\t\twithout conflict with annotation and support read\n")
    stat.write("\t\tsupport reat ratio of starting point is larger than %s\n" % (start_ratio))
    stat.write("\t\tsupport reat ratio of end point is larger than %s\n" % (end_ratio))
    print_num(num_conflict, stat, "all")
    if len(num_circular) > 2:
        for strain in num_circular.keys():
            if strain != "all":
                stat.write("\n" + strain + ":\n")
                stat.write("\tthe number of all circular RNAs = %s\n" % (num_circular[strain]))
                print_num(num_support, stat, strain)
                stat.write("\n\tthe circular RNAs:\n")
                stat.write("\t\twithout conflict with annotation and support read\n")
                stat.write("\t\tsupport reat ratio of starting point is larger than %s\n" % (start_ratio))
                stat.write("\t\tsupport reat ratio of end point is larger than %s\n" % (end_ratio))
                print_num(num_conflict, stat, strain)

def read_file(input_file, gff_file, circs):
    ps = splice_parser.parser_splice()
    for entry in ps.parser(input_file):
        if entry.supported_reads > high:
            high = entry.supported_reads
        circs.append(entry)
    gffs = []
    for entry in Gff3Parser().entries(open(gff_file)):
        gffs.append(entry)
    gffs = sorted(gffs, key=lambda k: (k.seq_id, k.start))
    return gffs

def circRNA_detection(input_file, gff_file, output_file, start_ratio, end_ratio, statistics):
    num_circular = {}
    num_circular["all"] = 0
    num_support = {}
    num_support["all"] = {}
    num_conflict = {}
    num_conflict["all"] = {}
    pre_seq_id = ""
    ps = splice_parser.parser_splice()
    circs = []
    high = 0
    gffs = read_file(input_file, gff_file, circs)
    out = open(args.output_file, "w")
    circs = sorted(circs, key = lambda x: (x.strain, x.supported_reads), reverse=True)
    out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % \
              ("ID","strain","strand","start","end","annotation_overlap",
               "supported_reads","supported_reads/reads_at_start",
               "supported_reads/reads_at_end"))
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
                        if (float(circ.supported_reads) / float(circ.start_site_reads) >= start_ratio) and \
                           (float(circ.supported_reads) / float(circ.end_site_reads) >= end_ratio):
                            import_num(support, num_conflict, circ.strain)
            num += 1
        pre_seq_id = circ.strain
    print_stat(statistics, num_circular, num_support, num_conflict, start_ratio, end_ratio)
if __name__ == "__main__":
    main()
