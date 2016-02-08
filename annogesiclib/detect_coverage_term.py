import os
import sys
import csv
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.coverage_detection import coverage_comparison, check_tex
from annogesiclib.lib_reader import read_libs, read_wig


def import_data(row):
    return{"method": "intersect_plus_minus", "strain": row[0],
           "start": int(row[1]), "end": int(row[2]), "name": row[3],
           "miss": int(row[4]), "loop": int(row[5]), "diff": [],
           "length": int(row[6]), "r_stem": int(row[7]), "strand": row[8],
           "l_stem": int(row[9]), "parent_p": row[10], "parent_m": row[11],
           "ut": int(row[12]), "print": False, "detect_p": False,
           "detect_m": False, "express": "False"}

#def compare_gff(term, gffs, seq):
#    first = True
#    detect = False
#    seq_start = -1
#    seq_end = -1
#    pre_gff = None
#    for gff in gffs:
#        if (term["strain"] == gff.seq_id) and \
#           (term["strand"] == gff.strand):
#            if first:
#                first = False
#                if (term["strand"] == "-") and \
#                   (term["end"] <= gff.end):
#                    seq_start = 1
#                    seq_end = gff.start
#                    detect = True
#            else:
#                if term["strand"] == "+":
#                    if (term["start"] >= pre_gff.start) and \
#                       (term["end"] <= gff.start):
#                        seq_start = pre_gff.end
#                        seq_end = gff.start
#                        detect = True
#                else:
#                    if (term["end"] <= gff.end) and \
#                       (term["start"] >= pre_gff.end):
#                        seq_start = pre_gff.end
#                        seq_end = gff.start
#                        detect = True
#            pre_gff = gff
#    if detect is not True:
#        if (term["strand"] == "+"):
#            for gff in reversed(gffs):
#                if gff.strand == "+":
#                    if (term["start"] >= gff.start):
#                        seq_start = gff.end
#                        seq_end = len(seq[term["strain"]])
#                    break
#    return (seq_start, seq_end)

def compare_ta(terms, tas, fuzzy):
    for term in terms:
        for ta in tas:
            start = ta.start - fuzzy
            end = ta.end + fuzzy
            if (ta.seq_id == term["strain"]) and (
                ta.strand == term["strand"]):
                if ((start <= term["start"]) and (
                     end >= term["start"]) and (
                     end <= term["end"])) or (
                    (start >= term["start"]) and (
                     end <= term["end"])) or (
                    (start >= term["start"]) and (
                     start <= term["end"]) and (
                     end >= term["end"])) or (
                    (start <= term["start"]) and (
                     end >= term["end"])):
                    term["express"] = "True"

def compare_transtermhp(hps, fr_terms):
    terms = []
    for term in fr_terms:
        detect = False
        for hp in hps:
            if (hp.seq_id == term["strain"]) and (
                hp.strand == term["strand"]):
                if ((hp.start <= term["start"]) and (
                     hp.end >= term["start"]) and (
                     hp.end <= term["end"])) or (
                    (hp.start >= term["start"]) and (
                     hp.end <= term["end"])) or (
                    (hp.start >= term["start"]) and (
                     hp.start <= term["end"]) and (
                     hp.end >= term["end"])) or (
                    (hp.start <= term["start"]) and (
                     hp.end >= term["end"])):
                    if hp.start < term["start"]:
                        term["start"] = hp.start
                    if hp.end > term["end"]:
                        term["end"] = hp.end
                    hp.attributes["print"] = True
                    detect = True
        if detect:
            term["method"] = "&".join([term["method"], "TransTermHP"])
        terms.append(term)
    for hp in hps:
        if "print" not in hp.attributes.keys():
            if hp.strand == "+":
                terms.append({"method": "TransTermHP", "strain": hp.seq_id,
                              "start": hp.start, "end": hp.end,
                              "strand": hp.strand, "name": hp.attributes["ID"],
                              "parent_p": hp.attributes["associated_gene"],
                              "print": False, "detect_p": False,
                              "detect_m": False, "express": "False",
                              "diff": []})
            else:
                terms.append({"method": "TransTermHP", "strain": hp.seq_id,
                              "start": hp.start, "end": hp.end,
                              "strand": hp.strand, "name": hp.attributes["ID"],
                              "parent_m": hp.attributes["associated_gene"],
                              "print": False, "detect_p": False,
                              "detect_m": False, "express": "False",
                              "diff": []})
    terms = sorted(terms, key=lambda x: (x["strain"], x["start"]))
    return terms

def compare_replicates(term_covers, template_texs, cond, tex_notex, replicates):
    detect_num = 0
    term_datas = []
    diff_cover = -1
    diff = []
    detect = False
    detect_num = check_tex(template_texs, term_covers, 0, term_datas, tex_notex,
                           None, None, "terminator", None, None, None, None)
    if ((detect_num >= replicates["tex"]) and \
       ("texnotex" in cond)) or \
       ((detect_num >= replicates["frag"]) and \
       ("frag" in cond)):
        detect = True
        for term in term_datas:
            if (len(diff) == 0) or (diff_cover < term["diff"]):
                diff_cover = term["diff"]
                diff = term
    return diff_cover, diff, term_datas, detect_num

def coverage2term(covers, term, fuzzy, hl_covers, hl_poss, strand,
                  decrease, term_covers, track):
    first = True
    for cover in covers:
        if (term["start"] <= cover["pos"] + fuzzy) and \
           (term["end"] >= cover["pos"] - fuzzy):
            first = coverage_comparison(cover, hl_covers, hl_poss,
                                        first, strand)
        else:
            if (strand == "+") and (cover["pos"] > term["end"] + fuzzy):
                break
            elif (strand == "-") and (cover["pos"] < term["start"] - fuzzy):
                break
        if (first is not True) and (hl_covers["high"] > 0):
            if ((hl_covers["low"] / hl_covers["high"]) < decrease) and (
                 hl_covers["low"] > -1):
                term_covers.append({"track": track, "high": hl_covers["high"],
                            "low": hl_covers["low"], "detect": "True",
                            "diff": (hl_covers["high"] - hl_covers["low"]),
                            "type": cover["type"]})
                break

def get_coverage(term, wigs, strand, template_texs, fuzzy,
                 decrease, replicates, tex_notex):
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
                    if strand == "-":
                        covers = reversed(covers[term["start"] - fuzzy - 2: \
                                          term["end"] + fuzzy + 1])
                    elif strand == "+":
                        covers = covers[term["start"] - 2: term["end"] + 1]
                    coverage2term(covers, term, fuzzy, hl_covers, hl_poss,
                                  strand, decrease, term_covers, track)
                if len(term_covers) != 0:
                    tmp_cover, tmp_diff, term_datas[cond], detect_nums[cond] = \
                            compare_replicates(term_covers, template_texs,
                                               cond, tex_notex, replicates)
                    if (diff_cover == -1) or (diff_cover < tmp_cover):
                        diff_cover = tmp_cover
                        diff = tmp_diff
    for cond, num in detect_nums.items():
        if (("texnotex" in cond) and (num >= replicates["tex"])) or (
            ("frag" in cond) and (num >= replicates["frag"])):
            if strand == "+":
                term["detect_p"] = True
            else:
                term["detect_m"] = True
    return diff_cover, diff, term_datas, detect_nums

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

def get_attribute_string(num, name, parent, diff, term, coverage, method):
    attribute_string = ";".join(
                 ["=".join(items) for items in [("ID", "term_" + str(num)),
                  ("Name", name), ("associate", parent),
                  ("coverage_decrease", coverage),
                  ("diff_coverage", diff),
                  ("express", term["express"]),
                  ("Method", method)]])
    return attribute_string

def print_table(term, cutoff_coverage, out_t, table_best):
    first = True
    if (term["express"] == "True") and \
       (term["diff_cover"] != -1):
        if term["diff"]["high"] >= cutoff_coverage:
            out_t.write("\tTrue\t")
            if not table_best:
                for datas in term["datas"].values():
                    for data in datas:
                        if first:
                            out_t.write("{0}(diff={1};high={2};low={3})".format(
                                        data["track"], data["diff"],
                                        data["high"], data["low"]))
                            first = False
                        else:
                            out_t.write(";{0}(diff={1};high={2};low={3})".format(
                                        data["track"], data["diff"],
                                        data["high"], data["low"]))
            else:
                out_t.write("{0}(diff={1};high={2};low={3})".format(
                            term["diff"]["track"], term["diff_cover"],
                            term["diff"]["high"], term["diff"]["low"]))
        else:
            out_t.write("\tFalse\t")
            out_t.write("No_coverage_decreasing")
    elif (term["express"] == "True") and \
         (term["diff_cover"] == -1):
        out_t.write("\tFalse\t")
        out_t.write("No_coverage_decreasing")
    elif term["express"] == "False":
        out_t.write("\tFalse\t")
        out_t.write("NA")

def print2file(num, term, coverage, parent, out, out_t,
               method, table_best, cutoff_coverage):
    name = 'Term_%0*d' % (5, num)
    if ("detect_num" in term.keys()) and \
       (term["diff_cover"] != -1):
        out_t.write("\t".join([term["strain"], name, str(term["start"]),
                              str(term["end"]), term["strand"], term["method"]]))
    else:
        out_t.write("\t".join([term["strain"], name, str(term["start"]),
                              str(term["end"]), term["strand"], term["method"]]))
    if (term["express"] == "True") and (term["diff_cover"] != -1):
        if (term["diff"]["high"] >= cutoff_coverage):
            diff = ("{0}(high:{1},low:{2})".format(
                    term["diff"]["track"], term["diff"]["high"],
                    term["diff"]["low"]))
            attribute_string = get_attribute_string(num, name, parent, 
                                                    diff, term, coverage, method)
        elif (term["diff"]["high"] < cutoff_coverage):
            attribute_string = get_attribute_string(num, name, parent, "NA",
                                                term, "No_coverage_decreasing", method)
    elif (term["express"] == "True") and (term["diff_cover"] == -1):
        attribute_string = get_attribute_string(num, name, parent, "NA",
                                                term, "No_coverage_decreasing", method)
    elif (term["express"] == "False"):
        attribute_string = get_attribute_string(num, name, parent, "NA",
                                                term, "NA", method)
    out.write("\t".join([str(field) for field in [
              term["strain"], "ANNOgesic", "terminator",
              str(term["start"]), str(term["end"]), ".",
              term["strand"], ".", attribute_string]]) + "\n")
    print_table(term, cutoff_coverage, out_t, table_best)
    out_t.write("\n")

def print_detect_undetect(terms, num, out, out_t, table_best,
                          cutoff_coverage, detect):
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
        if detect:
            num = print_detect_undetect(detect_terms["detect"], num, out,
                                    out_t, table_best, cutoff_coverage, "True")
            detect = False
        else:
            num = print_detect_undetect(detect_terms["undetect"], num, out,
                                    out_t, table_best, cutoff_coverage, "False")
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
                if not detect:
                    detect_terms["undetect"] = compare_term(
                                               term, detect_terms["undetect"])
        else:
            if (term["detect_m"]):
                detect = True
                detect_terms["detect"] = compare_term(
                                         term, detect_terms["detect"])
            else:
                if not detect:
                    detect_terms["undetect"] = compare_term(
                                               term, detect_terms["undetect"])
    return num, detect

def print_term(terms, out, out_t, table_best, cutoff_coverage):
    first = True
    detect = False
    detect_terms = {"detect": [], "undetect": []}
    num = 0
    for term in terms:
        if first:
            first = False
            pre_term = term
            detect = first_term(term["strand"], term, detect_terms, detect)
        else:
            num, detect = term_validation(pre_term, term, detect, detect_terms, out,
                                          out_t, table_best, num, cutoff_coverage)
            pre_term = term
    if detect:
        num = print_detect_undetect(detect_terms["detect"], num, out, out_t,
                                    table_best, cutoff_coverage, "True")
    else:
        num = print_detect_undetect(detect_terms["undetect"], num, out, out_t,
                                    table_best, cutoff_coverage, "False")


def del_repeat_term(terms):
    first = True
    new_terms = []
    for term in terms:
        detect = False
        if first:
            first = False
            pre_term = term
        else:
            if (term["strain"] == pre_term["strain"]) and (
                term["strand"] == pre_term["strand"]) and (
                term["parent_p"] == pre_term["parent_p"]) and (
                term["parent_m"] == pre_term["parent_m"]):
                if (term["start"] <= pre_term["start"]) and \
                   (term["end"] >= pre_term["start"]) and \
                   (term["end"] <= pre_term["end"]):
                    detect = True
                    pre_term["start"] = term["start"]
                elif (term["start"] <= pre_term["start"]) and \
                   (term["end"] >= pre_term["end"]):
                    detect = True
                    pre_term["start"] = term["start"]
                    pre_term["end"] = term["end"]
                elif (term["start"] >= pre_term["start"]) and \
                   (term["end"] <= pre_term["end"]):
                    detect = True
                elif (term["start"] >= pre_term["start"]) and \
                     (term["start"] <= pre_term["end"]) and \
                     (term["end"] >= pre_term["end"]):
                    detect = True
                    pre_term["end"] = term["end"]
                if detect:
                    if term["miss"] < pre_term["miss"]:
                        pre_term["miss"] = term["miss"]
                else:
                    new_terms.append(pre_term)
                    pre_term = term
            else:
                new_terms.append(pre_term)
                pre_term = term
    new_terms.append(term)
    return new_terms

def read_data(gff_file, tran_file, tranterm_file, seq_file, term_table):
    gff_parser = Gff3Parser()
    gffs = []
    tas = []
    hps = []
    fr_terms = []
    seq = {}
    for entry in gff_parser.entries(open(gff_file)):
        if (entry.feature == "CDS") or (
            entry.feature == "tRNA") or (
            entry.feature == "rRNA") or (
            entry.feature == "sRNA"):
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
    term_f = open(term_table, "r")
    for row in csv.reader(term_f, delimiter="\t"):
        fr_terms.append(import_data(row))
    new_terms = del_repeat_term(fr_terms)
    tas = sorted(tas, key=lambda x: (x.seq_id, x.start))
    gffs = sorted(gffs, key=lambda x: (x.seq_id, x.start))
    hps = sorted(hps, key=lambda x: (x.seq_id, x.start))
    return gffs, tas, hps, new_terms, seq

def compute_wig(wig_file, libs, terms, strand, texs, fuzzy,
                decrease, replicates, tex_notex):
    wigs = {}
    wigs = read_wig(wig_file, strand, libs)
    for term in terms:
        if (term["strand"] == strand) and \
           (term["express"] == "True"):
            term["diff_cover"], term["diff"], term["datas"], detect_nums = \
                get_coverage(term, wigs, strand, texs, fuzzy,
                             decrease, replicates, tex_notex)
            term["detect_num"] = {}
            for cond, num in detect_nums.items():
                term["detect_num"][cond] = str(num)

def detect_coverage(term_table, gff_file, tran_file, seq_file,
                    wig_f_file, wig_r_file, fuzzy, cutoff_coverage,
                    tranterm_file, wig_folder, input_libs, tex_notex,
                    replicates, output_file, output_table, table_best,
                    decrease):
    gffs, tas, hps, fr_terms, seq = read_data(gff_file, tran_file,
                                    tranterm_file, seq_file, term_table)
    terms = compare_transtermhp(hps, fr_terms)
    compare_ta(terms, tas, fuzzy)
    libs, texs = read_libs(input_libs, wig_folder)
    wigs = read_wig(wig_f_file, "+", libs)
    compute_wig(wig_f_file, libs, terms, "+", texs, fuzzy,
                decrease, replicates, tex_notex)
    compute_wig(wig_r_file, libs, terms, "-", texs, fuzzy,
                decrease, replicates, tex_notex)
    out = open(output_file, "w")
    out_t = open(output_table, "w")
    print_term(terms, out, out_t, table_best, cutoff_coverage)
