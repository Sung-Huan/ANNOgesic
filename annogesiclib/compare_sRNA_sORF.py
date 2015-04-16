#!/usr/bin/python

import os
import sys
import csv
from annogesiclib.gff3 import Gff3Parser

def print_file(datas, out, type_, feature):
    for data in datas:
        if feature not in data.attributes.keys():
            data.attributes[feature] = "NA"
        else:
            data.attributes[feature] = "&".join(data.attributes[feature])
        data.attribute_string = ";".join(
            ["=".join(items) for items in data.attributes.items()])
        out.write("\t".join([data.info_without_attributes, data.attribute_string]) + "\n")

def srna_sorf_comparison(sRNA_file, sORF_file, sRNA_out, sORF_out):
    sorfs = []
    srnas = []
    out_r = open(args.sRNA_out, "w")
    out_o = open(args.sORF_out, "w")
    gff_parser = gff3.Gff3Parser()
    out_r.write("##gff-version 3\n")
    out_o.write("##gff-version 3\n")
    for entry in gff_parser.entries(open(args.sRNA_file)):
        if "sORF" in entry.attributes.keys():
            del entry.attributes["sORF"]
        srnas.append(entry)
    srnas = sorted(srnas, key=lambda k: (k.seq_id, k.start))    
    for entry in gff_parser.entries(open(args.sORF_file)):
        if "sRNA" in entry.attributes.keys():
            del entry.attributes["sRNA"]
        sorfs.append(entry)
    sorfs = sorted(sorfs, key=lambda k: (k.seq_id, k.start))
    for srna in srnas:
        for sorf in sorfs:
            if (srna.seq_id == sorf.seq_id) and \
               (srna.strand == sorf.strand):
                if (srna.start <= sorf.start) and \
                   (srna.end >= sorf.end):
                    if "sORF" not in srna.attributes.keys():
                        srna.attributes["sORF"] = []
                    srna.attributes["sORF"].append("".join([sorf.attributes["ID"], ":",
                                              str(sorf.start), "-", str(sorf.end),
                                              "_", sorf.strand]))
                    if "sRNA" not in sorf.attributes.keys():
                        sorf.attributes["sRNA"] = []
                    sorf.attributes["sRNA"].append("".join([srna.attributes["ID"], ":",
                                              str(srna.start), "-", str(srna.end),
                                              "_", srna.strand]))
    print_file(sorfs, out_o, "sORF", "sRNA")
    print_file(srnas, out_r, "sRNA", "sORF")
