#!/usr/bin/python

import os	
import sys
import csv
import math
from annogesiclib.gff3 import Gff3Parser

def import_candidate(cands, term_features, strain, start, end, ut, name, total_length, strand,
                     parent_p, parent_m):
    cands.append({"strain": strain, "start": start, "end": end, "print": False, "ut": ut, "name": name,
                  "miss": term_features["real_miss"], "loop": term_features["loop"],
                  "length": total_length, "r_stem": term_features["r_stem"], "strand": strand,
                  "l_stem": term_features["l_stem"], "parent_p": parent_p, "parent_m": parent_m,
                  "detect_p": False, "detect_m": False})

def get_feature(cds):
    if "locus_tag" in cds.attributes.keys():
        feature = cds.attributes["locus_tag"]
    elif "protein_id" in cds.attributes.keys():
        feature = cds.attributes["protein_id"]
    else:
        feature = "".join([cds.feature, ":", str(cds.start),
                           "-", str(cds.end), "_", cds.strand])
    return feature

def filter_term(cands, terms):
    cutoff_miss = 0.25
    for cand1 in cands:
        detect = False
        stem_len = (cand1["r_stem"] + cand1["l_stem"] - cand1["miss"])
        if (not cand1["print"]) and \
           ((float(cand1["miss"]) / float(stem_len)) <= cutoff_miss):
            tmp_term = cand1.copy()
            for cand2 in cands:
                if (tmp_term["strain"] == cand2["strain"]) and \
                   (tmp_term["miss"] >= cand2["miss"]):
                    detect = True
                    if (tmp_term["start"] >= cand2["start"]) and \
                       (tmp_term["start"] < cand2["end"]) and \
                       (tmp_term["end"] > cand2["end"]):
                        tmp_term["start"] = cand2["start"]
                        cand2["print"] = True
                    elif (cand2["start"] > tmp_term["start"]) and \
                         (cand2["start"] < tmp_term["end"]) and \
                         (cand2["end"] >= tmp_term["end"]):
                        tmp_term["end"] = cand2["end"]
                        cand2["print"] = True
                    elif (tmp_term["start"] >= cand2["start"]) and \
                         (tmp_term["end"] <= cand2["end"]):
                        tmp_term["start"] = cand2["start"]
                        tmp_term["end"] = cand2["end"]
                        cand2["print"] = True
                    elif (cand2["start"] >= tmp_term["start"]) and \
                         (cand2["end"] <= tmp_term["end"]):
                        cand2["print"] = True
            terms.append(tmp_term)

def check_sec(sec, term_features, detects, nts):
    for st in reversed(sec[0:nts]):
        term_features["st_pos"] += 1
        if st == ")":
            if not detects["detect_l"]:
                term_features["rights"] += 1
                term_features["real_miss"] = term_features["tmp_miss"]
                detects["detect_r"] = True
            else:
                detects["conflict"] = True
                break
        elif st == ".":
            if detects["detect_r"]:
                term_features["tmp_miss"] += 1
        elif st == "(":
            term_features["lefts"] += 1
            if not detects["detect_l"]:
                term_features["loop"] = term_features["tmp_miss"] - term_features["real_miss"]
                term_features["tmp_miss"] = term_features["real_miss"]
                term_features["r_stem"] = term_features["rights"] + term_features["real_miss"]
            else:
                term_features["real_miss"] = term_features["tmp_miss"]
            detects["detect_l"] = True
            if term_features["lefts"] == term_features["rights"]:
                break

def detect_candidates(seq, sec, cands, name, strain, start, end, parent_p, parent_m, strand):
    for nts in range(0, len(seq) - 6):
        ut = 0
        for nt in seq[nts:nts + 6]:
            if (nt == "U") or (nt == "T"):
                ut += 1
        if (ut >= 3) and (nts > 10):
            if sec[nts - 1] == ")":
                term_features = {"st_pos": 0, "rights": 0, "lefts": 0, "tmp_miss": 0, 
                                 "real_miss": 0, "loop": 0, "r_stem": 0, "l_stem": 0}
                detects = {"detect_r": False, "detect_l": False, "conflict": False}
                check_sec(sec, term_features, detects, nts)
                if detects["conflict"] is False:
                    total_length = (nts) - (nts - term_features["st_pos"] + 1) + 1
                    term_features["l_stem"] = total_length - term_features["r_stem"] - \
                                              term_features["loop"]
                    if (total_length <= 60) and (term_features["loop"] <= 10) and \
                       (term_features["loop"] >= 3) and \
                       (((term_features["r_stem"] + term_features["l_stem"] - \
                          term_features["real_miss"]) / 2) >= 4):
                        if strand == "+":
                            import_candidate(cands, term_features, strain, 
                                             start + (nts - term_features["st_pos"]) - 10, 
                                             start + nts - 1 + 10, ut, name, total_length, 
                                             strand, parent_p, parent_m)
                        else:
                            import_candidate(cands, term_features, strain, end - (nts - 1) - 10, 
                                             end - (nts - term_features["st_pos"]) + 10, 
                                             ut, name, total_length, strand, parent_p, parent_m)

def check_parent(cdss, term, detects, strand, fuzzy_up, fuzzy_down, type_):
    tmp = None
    for cds in cdss:
        if (term["strain"] == cds.seq_id) and \
           (cds.strand == strand):
            if type_ == "parent_p":
#                print(term["start"] - fuzzy_down)
#                print(cds.end)
#                print(term["start"])
                if ((term["start"] - fuzzy_down) <= cds.end) and \
                   (term["start"] >= cds.end):
                    detects[type_] = True
                    tmp = get_feature(cds)
                elif ((cds.end - term["end"]) <= fuzzy_up) and \
                     ((cds.end - term["end"]) >= 0):
                    detects[type_] = True
                    tmp = get_feature(cds)
                elif ((cds.end - term["start"]) > fuzzy_up) and \
                     (cds.end - term["start"] >= 0):
                    break
            elif type_ == "parent_m":
#                print(term["end"] + fuzzy_down)
#                print(cds.start)
#                print(term["end"])
                if ((term["end"] + fuzzy_down) >= cds.start) and \
                   (term["end"] <= cds.start):
                    detects[type_] = True
                    tmp = get_feature(cds)
                elif ((term["start"] - cds.start) <= fuzzy_up) and \
                     ((term["start"] - cds.start) >= 0):
                    detects[type_] = True
                    tmp = get_feature(cds)
                elif (cds.start - term["end"] > fuzzy_down):
                    break
    return tmp

def parents(terms, cdss, fuzzy_up_cds, fuzzy_down_cds, fuzzy_up_ta, fuzzy_down_ta, tas):
    for term in terms:
#        print(term)
        check = False
        detects = {"parent_p": False, "parent_m": False}
        if "tran" in term["parent_p"]:
            tmp_p = check_parent(cdss, term, detects, "+", fuzzy_up_cds, fuzzy_down_cds, "parent_p")
            pos = term["parent_p"].split(":")[-1].split("_")[0].split("-")[-1]
            if ((term["start"] - int(pos) <= fuzzy_down_ta) and \
                (term["start"] - int(pos) >= 0)) or \
               ((int(pos) - term["end"] <= fuzzy_up_ta) and \
                (int(pos) - term["end"] >= 0)):
                pass
            else:
                term["parent_p"] = ""
        if "tran" in term["parent_m"]:
            tmp_m = check_parent(cdss, term, detects, "-", fuzzy_up_cds, fuzzy_down_cds, "parent_m")
            pos = term["parent_m"].split(":")[-1].split("_")[0].split("-")[0]
            if ((int(pos) - term["end"] <= fuzzy_down_ta) and \
                (int(pos) - term["end"] >= 0)) or \
               ((term["start"] - int(pos) <= fuzzy_up_ta) and \
                (term["start"] - int(pos) >= 0)):
                pass
            else:
                term["parent_m"] = ""
#        print(pos)
#        print(detects["parent_p"])
        if detects["parent_p"]:
            if len(term["parent_p"]) == 0:
                term["parent_p"] = tmp_p
            else:
                term["parent_p"] = ",".join([term["parent_p"], tmp_p])
#        print(detects["parent_m"])
        if detects["parent_m"]:
            if len(term["parent_m"]) == 0:
                term["parent_m"] = tmp_m
            else:
                term["parent_m"] = ",".join([term["parent_m"], tmp_m])

def read_gff(cdss, genome, trans, seq_file, gff_file, tran_file):
    for entry in Gff3Parser().entries(open(gff_file)):
        if (entry.feature == "CDS") or \
           (entry.feature == "tRNA") or \
           (entry.feature == "rRNA") or \
           (entry.feature == "sRNA"):
            cdss.append(entry)
    for entry in Gff3Parser().entries(open(tran_file)):
        trans.append(entry)
    with open(seq_file, "r") as q_h:
        for line in q_h:
            line = line.strip()
            if line.startswith(">"):
                strain = line[1:]
                genome[strain] = ""
            else:
                genome[strain] = genome[strain] + line

def compare_anno(gffs, cands, fuzzy_up, fuzzy_down):
    detect = False
    new_cands = []
    for cand in cands:
        for gff in gffs:
            if (gff.seq_id == cand["strain"]) and \
               (gff.strand == cand["strand"]):
                if cand["strand"] == "+":
                    if (gff.start <= cand["start"]) and \
                       (gff.end >= cand["start"]) and \
                       (gff.end <= cand["end"]):
                        detect = True
                        break
                    elif (math.fabs(gff.end - cand["end"]) <= fuzzy_up) and \
                         (gff.end >= cand["end"]):
                        detect = True
                        break
                    elif (math.fabs(gff.end - cand["end"]) <= fuzzy_down) and \
                         (gff.end <= cand["start"]):
                        detect = True
                        break
                else:
                    if (gff.start >= cand["start"]) and \
                       (gff.start <= cand["end"]) and \
                       (gff.end >= cand["end"]):
                        detect = True
                        break
                    elif (math.fabs(gff.start - cand["end"]) <= fuzzy_down) and \
                         (gff.start >= cand["end"]):
                        detect = True
                        break
                    elif (math.fabs(gff.start - cand["start"]) <= fuzzy_up) and \
                         (cand["start"] >= gff.start):
                        detect = True
                        break
        if detect:
            detect = False
            new_cands.append(cand)
    return new_cands
        
def merge_cands(new_cands_cds, new_cands_ta):
    new_cands = []
    for cand_cds in new_cands_cds:
        new_cands.append(cand_cds)
    for cand_ta in new_cands_ta:
        if cand_ta not in new_cands:
            new_cands.append(cand_ta)
    new_cands = sorted(new_cands, key=lambda k: (k["strain"], k["start"]))
    return new_cands

def poly_t(seq_file, sec_file, gff_file, tran_file, fuzzy_up_ta, fuzzy_down_ta, 
           fuzzy_up_cds, fuzzy_down_cds, out_file):
    genome = {}
    terms = []
    cdss = []
    trans = []
    read_gff(cdss, genome, trans, seq_file, gff_file, tran_file)
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start))
    trans = sorted(trans, key=lambda k: (k.seq_id, k.start))
    out = open(out_file, "w")
    with open(sec_file, "r") as s_h:
        for line in s_h:
            line = line.strip()
            if line.startswith(">"):
                line = line[1:]
                name = line.split("|")[0]
                start = int(line.split("|")[1])
                end = int(line.split("|")[2])
                strain = line.split("|")[3]
                parent_p = line.split("|")[4]
                parent_m = line.split("|")[5]
                strand = line.split("|")[6]
                sec = ""
                seq = ""
                for line in s_h:
                    if ("(" in line) or ("." in line) or (")" in line):
                        line = line.split(" ")
                        sec = line[0]
                        break
                    else:
                        seq = line
                if len(seq) <= 6:
                    continue
                else:
                    cands = []
                    detect_candidates(seq, sec, cands, name, strain, start, end, parent_p, parent_m, strand)
                cands = sorted(cands, key = lambda x: (x["miss"], x["start"]))
                new_cands_cds = compare_anno(cdss, cands, fuzzy_up_cds, fuzzy_down_cds)
                new_cands_ta = compare_anno(trans, cands, fuzzy_up_ta, fuzzy_down_ta)
                new_cands = merge_cands(new_cands_cds, new_cands_ta)
                filter_term(new_cands, terms)
    parents(terms, cdss,  fuzzy_up_cds, fuzzy_down_cds, fuzzy_up_ta, fuzzy_down_ta, trans)
    for term in terms:
        print_ = False
        if (term["strand"] == "+") and (len(term["parent_p"]) != 0):
            print_ = True
        elif (term["strand"] == "-") and (len(term["parent_m"]) != 0):
            print_ = True
        if print_:
            out.write("\t".join([term["strain"], str(term["start"]), str(term["end"]),
                                 term["name"], str(term["miss"]), str(term["loop"]),
                                 str(term["length"]), str(term["r_stem"]), term["strand"],
                                 str(term["l_stem"]), term["parent_p"], term["parent_m"],
                                 str(term["ut"])]) + "\n")
#    print_file(terms, genes, fuzzy_up, fuzzy_down, out)
