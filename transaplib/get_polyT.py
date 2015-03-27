#!/usr/bin/python

import os	
import sys
import csv
from transaplib.gff3 import Gff3Parser

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
        if (cand1["print"] is False) and \
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
            if detects["detect_l"] is not True:
                term_features["rights"] += 1
                term_features["real_miss"] = term_features["tmp_miss"]
                detects["detect_r"] = True
            else:
                detects["conflict"] = True
                break
        elif st == ".":
            if detects["detect_r"] is True:
                term_features["tmp_miss"] += 1
        elif st == "(":
            term_features["lefts"] += 1
            if detects["detect_l"] is not True:
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

def check_parent(cdss, term, detects, strand, utr_length, type_):
    tmp = None
    for cds in cdss:
        if (term["strain"] == cds.seq_id) and \
           (cds.strand == strand):
            if type_ == "parent_p":
                check_point = cds.end
            elif type_ == "parent_m":
                check_point = cds.start
            # when getting intergenic region already get upstream/downstream 50 nt
            if ((term["start"] - utr_length) <= check_point) and \
               (term["start"] >= check_point):
                detects[type_] = True
                tmp = get_feature(cds)
            elif (term["start"] < check_point):
                break
    return tmp

def parents(terms, cdss):
    for term in terms:
        detects = {"parent_p": False, "parent_m": False}
        if "tran" in term["parent_p"]:
            tmp_p = check_parent(cdss, term, detects, "+", 250, "parent_p")
        if "tran" in term["parent_m"]:
            tmp_m = check_parent(cdss, term, detects, "-", -1 * 250, "parent_m")
        if detects["parent_p"] is True:
            term["parent_p"] = ",".join(term["parent_p"], tmp_p)
        if detects["parent_m"] is True:
            term["parent_m"] = ",".join(term["parent_m"], tmp_m)

def read_gff(cdss, genome, seq_file, gff_file):
    for entry in Gff3Parser().entries(open(gff_file)):
        if (entry.feature == "CDS") or \
           (entry.feature == "tRNA") or \
           (entry.feature == "rRNA") or \
           (entry.feature == "sRNA"):
            cdss.append(entry)
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start))
    with open(seq_file, "r") as q_h:
        for line in q_h:
            line = line.strip()
            if line.startswith(">"):
                strain = line[1:]
                genome[strain] = ""
            else:
                genome[strain] = genome[strain] + line
    return cdss

def Poly_T(seq_file, sec_file, gff_file, fuzzy, out_file):
    genome = {}
    terms = []
    cdss = []
    cdss = read_gff(cdss, genome, seq_file, gff_file)
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
                filter_term(cands, terms)
    parents(terms, cdss)
    for term in terms:
        out.write("\t".join([term["strain"], str(term["start"]), str(term["end"]), 
                             term["name"], str(term["miss"]), str(term["loop"]),
                             str(term["length"]), str(term["r_stem"]), term["strand"], 
                             str(term["l_stem"]), term["parent_p"], term["parent_m"], 
                             str(term["ut"])]) + "\n")
