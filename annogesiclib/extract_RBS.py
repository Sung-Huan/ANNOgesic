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

def detect_site(inters, start_codons, start_rbs, end_rbs, fuzzy_rbs):
    rbss = []
    for inter in inters:
        for nts in range(0, len(inter["seq"]) - 6):
            num = 0
            miss = 0
            detect = False
            for nt in inter["seq"][nts:nts + 6]:
                if miss > fuzzy_rbs:
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
            if miss <= fuzzy_rbs:
                for start_codon in start_codons:
                    if (start_codon in inter["seq"][nts + (6 + start_rbs - 1):\
                                       nts + (6 + end_rbs)]):
                        rbss.append(inter)
                        detect = True
                        break
            if detect:
                break
    return rbss

def read_file(seq_file, gff_file):
    cdss = []
    seq = {}
    with open(seq_file, "r") as f_h:
        for line in f_h:
            line = line.strip()
            if line.startswith(">"):
                strain = line[1:]
                seq[strain] = ""
            else:
                seq[strain] = seq[strain] + line
    g_h = open(gff_file)
    for entry in Gff3Parser().entries(g_h):
        if (entry.feature == "CDS"):
            cdss.append(entry)
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start))
    g_h.close()
    return cdss, seq

def extract_seq(cdss, seq):
    first = True
    helper = Helper()
    inters = []
    pre_positive = None
    pre_minus = None
    for cds in cdss:
        if not first:
            if (cds.seq_id != pre_cds.seq_id):
                first = True
                if cds.strand == "+":
                    inter = helper.extract_gene(seq[cds.seq_id], 1, cds.start + 10, "+")
                    inters.append(import_data(inter, cds, 1, cds.start + 10))
                if pre_minus is not None:
                    inter = helper.extract_gene(seq[pre_minus.seq_id], pre_minus.end - 10, 
                                                len(seq[pre_minus.seq_id]), "-")
                    inters.append(import_data(inter, pre_minus, pre_minus.end - 10, 
                                              len(seq[pre_minus.seq_id])))
                if cds.strand == "+":
                    pre_positive = cds
                else:
                    pre_minus = cds
                pre_cds = cds
                continue
        if cds.strand == "+":
            if first:
                inter = helper.extract_gene(seq[cds.seq_id], 1, cds.start + 10, "+")
                inters.append(import_data(inter, cds, 1, cds.start + 10))
                first = False
            else:
                inter = helper.extract_gene(seq[cds.seq_id], pre_positive.end, cds.start + 10, "+")
                inters.append(import_data(inter, cds, pre_positive.end, cds.start + 10))
            pre_positive = cds
        else:
            if pre_minus is not None:
                inter = helper.extract_gene(seq[pre_minus.seq_id], pre_minus.end - 10, cds.start, "-")
                inters.append(import_data(inter, pre_minus, pre_minus.end - 10, cds.start))
            pre_minus = cds
        pre_cds = cds
    if pre_minus is not None:
        inter = helper.extract_gene(seq[pre_minus.seq_id], pre_minus.end - 10, 
                                    len(seq[pre_minus.seq_id]), "-")
        inters.append(import_data(inter, pre_minus, pre_minus.end - 10, 
                                  len(seq[pre_minus.seq_id])))
    return inters

def extract_potential_rbs(seq_file, gff_file, out_file, start_codons,
                          start_rbs, end_rbs, fuzzy_rbs):
    out = open(out_file, "w")
    cdss, seq = read_file(seq_file, gff_file)
    inters = extract_seq(cdss, seq)
    rbss = detect_site(inters, start_codons, start_rbs, end_rbs, fuzzy_rbs)
    num = 0
    for rbs in rbss:
        out.write(">riboswitch_{0}\n".format(
                  "|".join([str(num), rbs["strain"], rbs["strand"],
                  rbs["protein"], str(rbs["start"]), str(rbs["end"])])))
        out.write(rbs["seq"] + "\n")
        num += 1
    out.close()
