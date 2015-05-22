#!/usr/bin/python

import os	
import sys
import csv
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.helper import Helper

def get_feature(cds):
    if "locus_tag" in cds.attributes.keys():
        feature = cds.attributes["locus_tag"]
    elif "protein_id" in cds.attributes.keys():
        feature = cds.attributes["protein_id"]
    else:
        feature = cds.attributes["ID"]
    return feature

def import_data(seq, cds, start, end):
    feature = get_feature(cds)
    return {"seq": seq, "strain": cds.seq_id, "strand": cds.strand,
            "protein": feature, "start": start, "end": end}

def detect_site(inters, rbss):
    for inter in inters:
        for nts in range(0, len(inter["seq"]) - 6):
            num = 0
            miss = 0
            for nt in inter["seq"][nts:nts + 6]:
                if miss == 2:
                    break
                else:
                    if (num == 0) and (nt != "A"):
                        miss += 1
                    elif (num == 1) and (nt != "G"):
                        miss += 1
                    elif (num == 2) and (nt != "G"):
                        miss += 1
                    elif (num == 3) and (nt != "A"):
                        miss += 1
                    elif (num == 4) and (nt != "G"):
                        miss += 1
                    elif (num == 5) and (nt != "G"):
                        miss += 1
                    num += 1
            if miss < 2:
                if ("ATG" in inter["seq"][nts + 10:nts + 21]) or \
                   ("GTG" in inter["seq"][nts + 10:nts + 21]) or \
                   ("TTG" in inter["seq"][nts + 10:nts + 21]):
                    rbss.append(inter)
                    break

def read_file(seq_file, gff_file, seq):
    cdss = []
    with open(seq_file, "r") as f_h:
        for line in f_h:
            line = line.strip()
            if line.startswith(">"):
                strain = line[1:]
                seq[strain] = ""
            else:
                seq[strain] = seq[strain] + line
    for entry in Gff3Parser().entries(open(gff_file)):
        if (entry.feature == "CDS"):
            cdss.append(entry)
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start))
    return cdss

def extract_seq(cdss, inters, seq):
    first = True
    helper = Helper()
    for cds in cdss:
        if not first:
            if cds.seq_id != pre_cds.seq_id:
                first = True
                inter = helper.extract_gene(seq[cds.seq_id], 1, cds.start + 10, "+")
                inters.append(import_data(inter, cds, 1, cds.start + 10))
                inter = helper.extract_gene(seq[pre_minus.seq_id], pre_minus.end - 10, 
                                            len(seq[pre_minus.seq_id]), "-")
                inters.append(import_data(inter, pre_minus, pre_minus.end - 10, 
                                          len(seq[pre_minus.seq_id])))
                pre_cds = cds
                continue
        if cds.strand == "+":
            if first:
                inter = helper.extract_gene(seq[cds.seq_id], 1, cds.start + 10, "+")
                inters.append(import_data(inter, cds, 1, cds.start + 10))
                first = False
            else:
                inter = helper.extract_gene(seq[cds.seq_id], pre_cds.end, cds.start + 10, "+")
                inters.append(import_data(inter, cds, pre_cds.end, cds.start + 10))
        else:
            inter = helper.extract_gene(seq[pre_cds.seq_id], pre_cds.end - 10, cds.start, "-")
            inters.append(import_data(inter, pre_cds, pre_cds.end - 10, cds.start))
            pre_minus = cds
        pre_cds = cds
    inter = helper.extract_gene(seq[pre_minus.seq_id], pre_minus.end - 10, 
                                len(seq[pre_minus.seq_id]), "-")
    inters.append(import_data(inter, pre_minus, pre_minus.end - 10, 
                              len(seq[pre_minus.seq_id])))

def extract_potential_rbs(seq_file, gff_file, out_file):
    seq = {}
    rbss = []
    out = open(out_file, "w")
    cdss = read_file(seq_file, gff_file, seq)
    inters = []
    extract_seq(cdss, inters, seq)
    detect_site(inters, rbss)
    num = 0
    for rbs in rbss:
        out.write(">riboswitch_{0}\n".format(
                  "|".join([str(num), rbs["strain"], rbs["strand"],
                  rbs["protein"], str(rbs["start"]), str(rbs["end"])])))
        out.write(rbs["seq"] + "\n")
        num += 1
