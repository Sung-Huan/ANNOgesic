#!/usr/bin/python

import os	
import sys
import csv
from transaplib.gff3 import Gff3Parser
from transaplib.coverage_detection import Coverage_comparison, Check_Tex
from transaplib.lib_reader import Read_libs, Read_wig


def import_data(row):
    return{"method": "forward_reverse", "strain": row[0], "start": int(row[1]), "end": int(row[2]),
           "name": row[3], "miss": int(row[4]), "loop": int(row[5]), "diff": [],
           "length": int(row[6]), "r_stem": int(row[7]), "strand": row[8],
           "l_stem": int(row[9]), "parent_p": row[10], "parent_m": row[11], "ut": int(row[12]), 
           "print": False, "detect_p": False, "detect_m": False, "express": "False"}

def compare_gff(term, gffs, seq):
    first = True
    detect = False
    seq_start = -1
    seq_end = -1
    for gff in gffs:
        if (term["strain"] == gff.seq_id) and \
           (term["strand"] == gff.strand):
            if first:
                first = False
                if (term["strand"] == "-") and \
                   (term["end"] <= gff.end):
                    seq_start = 1
                    seq_end = gff.start
                    detect = True
            else:
                if term["strand"] == "+":
                    if (term["start"] >= pre_gff.start) and \
                       (term["end"] <= gff.start):
                        seq_start = pre_gff.end
                        seq_end = gff.start
                        detect = True
                else:
                    if (term["end"] <= gff.end) and \
                       (term["start"] >= pre_gff.end):
                        seq_start = pre_gff.end
                        seq_end = gff.start
                        detect = True
            pre_gff = gff
    if detect is not True:
        if (term["strand"] == "+"):
            for gff in reversed(gffs):
               if gff.strand == "+":
                   if (term["start"] >= gff.start):
                       seq_start = gff.end
                       seq_end = len(seq[term["strain"]])
                   break
    return (seq_start, seq_end)

def check_overlap(strain, strand, start, end, tas):
    for ta in tas:
        if (ta.seq_id == strain) and \
           (ta.strand == strand):
            if ((ta.start < start) and (ta.end > start) and (ta.end < end)) or \
               ((ta.start > start) and (ta.end < end)) or \
               ((ta.start > start) and (ta.start < end) and (ta.end > end)) or \
               ((ta.start < start) and (ta.end > end)):
                return "True"
    return "False"

def compare_ta(terms, tas, gffs, seq):
    for term in terms:
        pos = compare_gff(term, gffs, seq)
        start = pos[0]
        end = pos[1]
        if (start == -1):
            term["conflict"] = True
        else:
            term["express"] = check_overlap(term["strain"], term["strand"], start, end, tas)

def compare_transtermhp(hps, fr_terms, terms):
    for term in fr_terms:
        detect = False
        for hp in hps:
            if (hp.seq_id == term["strain"]) and \
               (hp.strand == term["strand"]):
                if ((hp.start < term["start"]) and (hp.end > term["start"]) and (hp.end < term["end"])) or \
                   ((hp.start > term["start"]) and (hp.end < term["end"])) or \
                   ((hp.start > term["start"]) and (hp.start < term["end"]) and (hp.end > term["end"])) or \
                   ((hp.start < term["start"]) and (hp.end > term["end"])):
                    if hp.start < term["start"]:
                        term["start"] = hp.start
                    if hp.end > term["end"]:
                        term["end"] = hp.end
                    hp.attributes["print"] = True
                    detect = True
        if detect:
            term["method"] = term["method"] + "&" + "TransTermHP"
        terms.append(term)
    for hp in hps:
        if "print" not in hp.attributes.keys():
            if hp.strand == "+":
                terms.append({"method": "TransTermHP", "strain": hp.seq_id, "start": hp.start, "end": hp.end,
                              "strand": hp.strand, "name": hp.attributes["ID"], "parent_p": hp.attributes["associated_gene"],
                              "print": False, "detect_p": False, "detect_m": False, "express": "False", "diff": []})
            else:
                terms.append({"method": "TransTermHP", "strain": hp.seq_id, "start": hp.start, "end": hp.end,
                              "strand": hp.strand, "name": hp.attributes["ID"], "parent_m": hp.attributes["associated_gene"],
                              "print": False, "detect_p": False, "detect_m": False, "express": "False", "diff": []})

def compare_replicates(term_covers, template_texs, cond, tex_notex, replicates):
    detect_num = 0
    term_datas = []
    diff_cover = -1
    diff = ""
    detect = False
    detect_num = Check_Tex(template_texs, term_covers, 0, term_datas, tex_notex,
                           detect_num, "terminator", None, None, None, None)
    if detect_num >= replicates:
        detect = True
        for term in term_datas:
            if (len(diff) == 0) or (diff_cover < term["diff"]):
                diff_cover = term["diff"]
                diff = term
    return (detect, diff_cover, diff, term_datas, detect_num)

def coverage2term(covers, term, fuzzy, hl_covers, hl_poss, strand, 
                  decrease, term_covers, track):
    first = True
    for cover in covers:
        if (term["start"] <= cover["pos"] + fuzzy) and \
           (term["end"] >= cover["pos"] - fuzzy):
            first = Coverage_comparison(cover, hl_covers, hl_poss, first, strand)
        else:
            if (strand == "+") and (cover["pos"] > term["end"] + fuzzy):
                break
            elif (strand == "-") and (cover["pos"] < term["start"] - fuzzy):
                break
        if (first is not True) and (hl_covers["high"] > 0):
            if ((hl_covers["low"] / hl_covers["high"]) < decrease) and \
               (hl_covers["low"] > -1):
                term_covers.append({"track": track, "high": hl_covers["high"],
                                    "low": hl_covers["low"], "detect": "True",
                                    "diff": (hl_covers["high"] - hl_covers["low"]),
                                    "type": cover["type"]})
                break

def detect_coverage(term, wigs, strand, template_texs, fuzzy, decrease, replicates, tex_notex):
    hl_poss = {"high": 0, "low": 0}
    hl_covers = {"high": 0, "low": 0}
    term_datas = {}
    detect_nums = {}
    diff_cover = -1
    diff = []
    for wig_strain, conds in wigs.items():
        if wig_strain == term["strain"]:
            for cond, tracks in conds.items():
                term_covers = []
                term_datas[cond] = []
                detect_nums[cond] = 0
                for track, covers in tracks.items():
                    first = True
                    if strand == "-":
                        covers = reversed(covers[term["start"] - fuzzy - 2: term["end"] + fuzzy + 1])
                    elif strand == "+":
                        covers = covers[term["start"] - 2: term["end"] + 1]
                    coverage2term(covers, term, fuzzy, hl_covers, hl_poss, strand, 
                                  decrease, term_covers, track)
                if len(term_covers) != 0:
                    datas = compare_replicates(term_covers, template_texs, cond, tex_notex, replicates)
                    term_datas[cond] = datas[3]
                    detect_nums[cond] = datas[4]
                    if (diff_cover == -1) or (diff_cover < datas[1]):
                        diff_cover = datas[1]
                        diff = datas[2]
    for cond, num in detect_nums.items():
        if num >= replicates:
            if strand == "+":
                term["detect_p"] = True
            else:
                term["detect_m"] = True
    return (diff_cover, diff, term_datas, detect_nums)

def compare_term(term, terms):
    if len(terms) != 0:
        for tmp in terms:
            if term["miss"] < tmp["miss"]:
                terms = []
                terms.append(term)
                break
            elif term["miss"] == tmp["miss"]:
                if ("diff_cover" in term.keys()) and \
                   ("diff_cover" in tmp.keys()):
                    if (term["diff_cover"] > tmp["diff_cover"]):
                        terms = []
                        terms.append(term)
                    elif (term["diff_cover"] == tmp["diff_cover"]):
                        if term["ut"] > tmp["ut"]:
                            terms = []
                            terms.append(term)
                        elif term["ut"] == tmp["ut"]:
                            terms.append(term)
                break
            elif term["miss"] > tmp["miss"]:
                break
    else:
        terms.append(term)
    return terms

def first_term(strand, term, detect_terms, detect):
    if (strand == "+"):
        if (term["detect_p"]):
            detect_terms["detect"].append(term)
            detect = True
        else:
            detect_terms["undetect"].append(term)
    elif (strand == "-"):
        if (term["detect_m"]):
            detect_terms["detect"].append(term)
            detect = True
        else:
            detect_terms["undetect"].append(term)
    return detect

def get_attribute_string(num, name, parent, diff, term, coverage):
    attribute_string = ";".join(
                 ["=".join(items) for items in [("ID", "term_" + str(num)),
                  ("Name", name), ("associate", parent),
                  ("coverage_decrease", coverage),
                  ("diff_coverage", diff),
                  ("express", term["express"])]])
    return attribute_string

def print_table(term, cutoff_coverage, out_t, table_best):
    first = True
    if (term["express"] == "True") and \
       (term["diff_cover"] != -1):
        if term["diff"]["high"] >= cutoff_coverage:
            out_t.write("True\t")
        else:
            out_t.write("False\t")
        if table_best is not False:
            for cond, datas in term["datas"].items():
                for data in datas:
                    if first:
                        out_t.write("%s(diff=%s;high=%s;low=%s)" % (
                                    data["track"], data["diff"], 
                                    data["high"], data["low"]))
                        first = False
                    else:
                        out_t.write(";%s(diff=%s;high=%s;low=%s)" % (
                                    data["track"], data["diff"], 
                                    data["high"], data["low"]))
        else:
            out_t.write("%s(diff=%s;high=%s;low=%s)" % (term["diff"]["track"], 
                        term["diff_cover"], term["diff"]["high"], term["diff"]["low"]))
    elif (term["express"] == "True") and \
         (term["diff_cover"] == -1):
        out_t.write("False\t")
        out_t.write("No_coverage_decreasing")
    elif term["express"] == "False":
        out_t.write("False\t")
        out_t.write("NA")

def print2file(num, term, coverage, parent, out, out_t, method, table_best, cutoff_coverage):
    name='Term_%0*d' % (5, num)
    if ("detect_num" in term.keys()) and \
       (term["diff_cover"] != -1):
        out_t.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % \
                   (term["strain"], name, term["start"], term["end"], term["strand"],
                    ";".join(term["detect_num"].keys()), 
                    ";".join(term["detect_num"].values())))
    else:
        out_t.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % \
                   (term["strain"], name, term["start"], term["end"], term["strand"],
                    "NA", "NA"))
    if (term["express"] == "True") and (term["diff_cover"] != -1):
        if (term["diff"]["high"] >= cutoff_coverage):
            diff = ("%s(high:%s,low:%s)" % (term["diff"]["track"], 
                                            term["diff"]["high"], 
                                            term["diff"]["low"]))
            attribute_string = get_attribute_string(num, name, parent, 
                                                    diff, term, coverage)
        elif (term["diff"]["high"] < cutoff_coverage):
            attribute_string = get_attribute_string(num, name, parent, "NA", 
                                                    term, "No_coverage_decreasing")
    elif (term["express"] == "True") and (term["diff_cover"] == -1):
        attribute_string = get_attribute_string(num, name, parent, "NA", 
                                                term, "No_coverage_decreasing")
    elif (term["express"] == "False"):
        attribute_string = get_attribute_string(num, name, parent, "NA",
                                                term, "NA")
    out.write("\t".join([str(field) for field in [
              term["strain"], method, "terminator",
              str(term["start"]), str(term["end"]), ".",
              term["strand"], ".", attribute_string]]) + "\n")
    print_table(term, cutoff_coverage, out_t, table_best)
    out_t.write("\n")

def print_detect_undetect(terms, num, out, out_t, table_best, cutoff_coverage, detect):
    for term in terms:
        if term["strand"] == "+":
            print2file(num, term, detect, term["parent_p"], out, 
                       out_t, term["method"], table_best, cutoff_coverage)
            num += 1
        else:
            print2file(num, term, detect, term["parent_m"], out, 
                       out_t, term["method"], table_best, cutoff_coverage)
            num += 1
    return num

def term_validation(pre_term, term, detect, detect_terms, out, out_t, 
                    table_best, num, cutoff_coverage):
    if pre_term["name"] != term["name"]:
        if detect is True:
            num = print_detect_undetect(detect_terms["detect"], num, out, out_t, 
                                        table_best, cutoff_coverage, "True")
            detect = False
        else:
            num = print_detect_undetect(detect_terms["undetect"], num, out, out_t, 
                                        table_best, cutoff_coverage, "False")
        detect_terms["detect"] = []
        detect_terms["undetect"] = []
        detect = first_term(term["strand"], term, detect_terms, detect)
    else:
        if term["strand"] == "+":
            if (term["detect_p"]):
                detect = True
                detect_terms["detect"] = compare_term(
                                         term, detect_terms["detect"])
            else:
                if detect is False:
                    detect_terms["undetect"] = compare_term(
                                               term, detect_terms["undetect"])
        else:
            if (term["detect_m"]):
                detect = True
                detect_terms["detect"] = compare_term(
                                         term, detect_terms["detect"])
            else:
                if detect is False:
                    detect_terms["undetect"] = compare_term(
                                               term, detect_terms["undetect"])
    return (num, detect)

def print_term(terms, out, out_t, table_best, cutoff_coverage):
    first = True
    detect = False
    detect_terms = {"detect": [], "undetect": []}
    num = 0
    for term in terms:
        if "conflict" not in term.items():
            if first:
                first = False
                pre_term = term
                detect = first_term(term["strand"], term, detect_terms, detect)
            else:
                data = term_validation(pre_term, term, detect, detect_terms, out, 
                                      out_t, table_best, num, cutoff_coverage)
                num = data[0]
                detect = data[1]
                pre_term = term
    if detect is True:
        num = print_detect_undetect(detect_terms["detect"], num, out, out_t, 
                                    table_best, cutoff_coverage, "True")
    else:
        num = print_detect_undetect(detect_terms["undetect"], num, out, out_t, 
                                    table_best, cutoff_coverage, "False")

def read_data(gffs, tas, hps, seq, gff_file, tran_file, 
              tranterm_file, seq_file, term_table, fr_terms):
    gff_parser = Gff3Parser()
    for entry in gff_parser.entries(open(gff_file)):
        if (entry.feature == "CDS") or \
           (entry.feature == "tRNA") or \
           (entry.feature == "rRNA") or \
           (entry.feature == "sRNA"):
            gffs.append(entry)
    for entry in gff_parser.entries(open(tran_file)):
        tas.append(entry)
    for entry in gff_parser.entries(open(tranterm_file)):
        hps.append(entry)
    with open(seq_file, "r") as s_f:
        for line in s_f:
            line = line.strip()
            if line.startswith(">"):
                strain = line[1:]
                seq[strain] = ""
            else:
                seq[strain] = seq[strain] + line
    fh = open(term_table, "r")
    for row in csv.reader(fh, delimiter="\t"):
        fr_terms.append(import_data(row))

def input_final_terms(term, covers):
    term["diff_cover"] = covers[0]
    term["diff"] = covers[1]
    term["datas"] = covers[2]
    term["detect_num"] = {}
    for cond, num in covers[3].items():
        term["detect_num"][cond] = str(num)
    return term

def compute_wig(wig_file, libs, terms, strand, wigs, texs, fuzzy, 
                decrease, replicates, tex_notex):
    Read_wig(wigs, wig_file, strand, libs)
    for term in terms:
        if (term["strand"] == strand) and \
           (term["express"] == "True") and \
           ("conflict" not in term.keys()):
            covers = detect_coverage(term, wigs, strand, texs, fuzzy, 
                                     decrease, replicates, tex_notex)
            term = input_final_terms(term, covers)

def Detect_coverage(term_table, gff_file, tran_file, seq_file, wig_f_file, wig_r_file, 
                    fuzzy, cutoff_coverage, tranterm_file, wig_folder, input_libs, 
                    tex_notex, replicates, output_file, output_table, table_best, decrease):
    wigs_f = {}
    wigs_r = {}
    fr_terms = []
    terms = []
    hps = []
    texs = {}
    libs = []
    tas = []
    gffs = []
    seq = {}
    read_data(gffs, tas, hps, seq, gff_file, tran_file, tranterm_file, 
              seq_file, term_table, fr_terms)
    tas = sorted(tas, key = lambda x: (x.seq_id, x.start))
    gffs = sorted(gffs, key = lambda x: (x.seq_id, x.start))
    hps = sorted(hps, key = lambda x: (x.seq_id, x.start))
    compare_transtermhp(hps, fr_terms, terms)
    terms = sorted(terms, key = lambda x: (x["strain"], x["start"]))
    compare_ta(terms, tas, gffs, seq)
    Read_libs(libs, texs, input_libs, wig_folder)
    compute_wig(wig_f_file, libs, terms, "+", wigs_f, texs, fuzzy, decrease, replicates, tex_notex)
    compute_wig(wig_r_file, libs, terms, "-", wigs_r, texs, fuzzy, decrease, replicates, tex_notex)
    out = open(output_file, "w")
    out_t = open(output_table, "w")
    print_term(terms, out, out_t, table_best, cutoff_coverage)
