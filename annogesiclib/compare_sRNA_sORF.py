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

def del_attributes(feature, entry):
    attributes = {}
    for key, value in entry.attributes.items():
        if feature not in key:
            attributes[key] = value
    return attributes

def srna_sorf_comparison(sRNA_file, sORF_file, sRNA_out, sORF_out):
    sorfs = []
    srnas = []
    out_r = open(sRNA_out, "w")
    out_o = open(sORF_out, "w")
    out_r.write("##gff-version 3\n")
    out_o.write("##gff-version 3\n")
    for entry in Gff3Parser().entries(open(sRNA_file)):
        entry.attributes = del_attributes("sORF", entry)
        srnas.append(entry)
    srnas = sorted(srnas, key=lambda k: (k.seq_id, k.start))    
    for entry in Gff3Parser().entries(open(sORF_file)):
        entry.attributes = del_attributes("sRNA", entry)
        sorfs.append(entry)
    sorfs = sorted(sorfs, key=lambda k: (k.seq_id, k.start))
    for srna in srnas:
        for sorf in sorfs:
            if (srna.seq_id == sorf.seq_id) and \
               (srna.strand == sorf.strand):
                if (srna.start <= sorf.start) and \
                   (srna.end >= sorf.end):
                    print(sorf.attributes)
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
