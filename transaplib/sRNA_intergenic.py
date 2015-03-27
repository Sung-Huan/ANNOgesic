from __future__ import division
import sys
import math
import csv
import os
from transaplib.gff3 import Gff3Parser
from transaplib.coverage_detection import Coverage_comparison, Replicate_comparison
from transaplib.lib_reader import Read_wig, Read_libs


def get_differential_cover(fuzzy_end, num, checks, cover_set, 
                           poss, cover, decrease):
    go_out = False
    if checks["detect_diff"]:
        if (num == fuzzy_end) or \
           (cover_diff == 0) or \
           ((cover["coverage"] > cover_sets["diff"]) and \
            (cover["coverage"] / cover_sets["diff"]) > (1 + decrease)):
            poss["stop_point"] = cover["pos"]
            go_out = True
        elif (cover["coverage"] <= cover_sets["diff"]):
            if (cover["coverage"] / cover_sets["diff"]) <= (decrease / 2):
                num += 1
            else:
                num = 0
            cover_sets["diff"] = cover["coverage"]
            cover_sets["low"] = cover["coverage"]
        elif (cover["coverage"] > cover_sets["diff"]) and \
             ((cover["coverage"] / cover_sets["diff"]) <= (1 + decrease)):
            num += 1
    if (checks["first"] is not True) and (cover_sets["high"] > 0):
        if ((cover_sets["low"] / cover_sets["high"]) < decrease) and \
           (cover_sets["low"] > -1):
            checks["detect_diff"] = True
            cover_sets["diff"] = cover["coverage"]
    return go_out

def get_best(wigs, strain, strand, start, end, type_, 
             decrease, cutoff_coverage, fuzzy_end):
    cover_sets = {"low": -1, "high": -1, "total": 0, "diff": 0}
    poss = {"high": 0, "low": 0, "stop_point": -1}
    srna_covers = {}
    for wig_strain, conds in wigs.items():
        if wig_strain == strain:
            for cond, tracks in conds.items():
                srna_covers[cond] = []
                for track, covers in tracks.items():
                    cover_sets["total"] = 0
                    cover_sets["diff"] = 0
                    checks = {"first": True, "detect_dff": False}
                    num = 0
                    if strand == "+":
                        covers = covers[start-2:end+1]
                    elif strand == "-":
                        covers = reversed(covers[start-2:end+1])
                    for cover in covers:
                        if (cover["strand"] == strand):
                            if (start <= cover["pos"]) and \
                               (end >= cover["pos"]):
                                cover_sets["total"] = cover_sets["total"] + cover["coverage"]
                                checks["first"] = Coverage_comparison(
                                                  cover, cover_sets, poss, checks["first"], strand)
                            else:
                                if (strand == "+") and (cover["pos"] > end):
                                    poss["stop_point"] = cover["pos"]
                                    break
                                elif (strand == "-") and (cover["pos"] < start):
                                    poss["stop_point"] = cover["pos"]
                                    break
                            if type_ == "differential":
                                go_out = get_differential_cover(fuzzy_end, num, checks, 
                                                      cover_set, poss, cover, decrease)
                                if go_out is True:
                                    break
                    avg = cover_sets["total"] / float(end - start + 1)
                    if avg > float(cutoff_coverage):
                        srna_covers[cond].append({"track": track, "high": cover_sets["high"], 
                                                  "low": cover_sets["low"], "avg": avg, 
                                                  "pos": poss["stop_point"], "type": cover["type"]})
    return srna_covers

def get_attribute_string(srna_datas, tss, num, name):
    attribute_string = ";".join(
                    ["=".join(items) for items in (["ID", "srna" + num],
                     ["Name", "sRNA_candidates_" + name])])
    with_tss = "=".join(["with_TSS", tss])
    if srna_datas is None:
        if tss != "NA":
            attribute_string = ";".join([attribute_string, with_tss])
    else:
        srna_data_string = ";".join(
                ["=".join(items) for items in (
                 ["best_avg_coverage", srna_datas["best"]],
                 ["best_high_coverage", srna_datas["high"]],
                 ["best_low_coverage", srna_datas["low"]])])
        if tss != "NA":
            attribute_string = ";".join([attribute_string, with_tss, srna_data_string])
        else:
            attribute_string = ";".join([attribute_string, srna_data_string])
    return attribute_string

def print_file(string, nums, tss, output, out_table, srna_datas, table_best):
    name='%0*d' % (5, nums["uni"])
    datas = string.split("\t")
    out_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % \
                   (datas[0], name, datas[3], datas[4], datas[6], 
                    ";".join(srna_datas["conds"].keys()), 
                    ";".join(srna_datas["conds"].values()),
                    srna_datas["best"], srna_datas["high"], srna_datas["low"]))
    get_attribute_string(srna_datas, tss, num, name)
    output.write("\t".join(string, attribute_string) + "\n")
    nums["uni"] += 1
    if (srna_datas is None):
        out_table.write(tss + "\n")
    else:
        if srna_datas["detail"] is not None:
            out_table.write(tss + "\t")
            if table_best is False:
                first = True
                for data in srna_datas["detail"]:
                    if first:
                        out_table.write("%s(avg=%s;high=%s;low=%s)" % (
                                        data["track"], data["avg"], data["high"], data["low"]))
                        first = False
                    else:
                        out_table.write(";%s(avg=%s;high=%s;low=%s)" % (
                                        data["track"], data["avg"], data["high"], data["low"]))
            else:
                out_table.write("%s(avg=%s;high=%s;low=%s)" % (
                                srna_datas["track"], srna_datas["best"], 
                                srna_datas["high"], srna_datas["low"]))
        out_table.write("\n")

def get_coverage(start, end, strain, wigs, strand, ta, nums, tss, 
                 output, template_texs, out_table, cutoff_coverage,
                 tex_notex, replicates, decrease, fuzzy_end, table_best):
    srna_covers = get_best(wigs, strain, strand, start, end, "total",
                           decrease, cutoff_coverage, fuzzy_end)
    srna_datas = Repliate_comparison(srna_covers, template_texs, strand, 
                                     cutoff_coverage, tex_notex, replicates, 
                                     "normal", None, None, None)
    string = ("\t".join([str(field) for field in [
                             ta.seq_id, ta.source, "sRNA", str(start),
                             str(end), ".", ta.strand, "."]]))
    if datas[0] != 0:
        print_file(string, nums, tss, output, out_table, srna_datas, table_best)

def detect_wig_pos(wigs, ta, start, end, nums, output, tss, template_texs, 
                   out_table, cutoff_coverage, min_len, Max_len, decrease, 
                   fuzzy_end, table_best):
    srna_covers = get_best(wigs, ta.seq_id, ta.strand, start, end, "differential",
                           decrease, cutoff_coverage, fuzzy_end)
    srna_datas = Repliate_comparison(srna_covers, template_texs, strand, 
                                     cutoff_coverage, tex_notex, replicates, 
                                     "normal", None, None, None)
    best = datas[0]
    high = datas[1]
    low = datas[2]
    pos = datas[3]
    if best > cutoff_coverage:
        if ta.strand == "+":
            if ((srna_datas["pos"] - start) >= min_len) and \
               ((srna_datas["pos"] - start) <= Max_len):
                string = ("\t".join([str(field) for field in [
                          ta.seq_id, ta.source, "sRNA", str(start),
                          str(srna_datas["pos"]), ".", ta.strand, "."]]))
                print_file(string, nums, tss, output, 
                           out_table, srna_datas, table_best)
        else:
            if ((end - srna_datas["pos"]) >= min_len) and \
               ((end - srna_datas["pos"]) <= Max_len):
                string = ("\t".join([str(field) for field in [
                          ta.seq_id, ta.source, "sRNA", str(srna_datas["pos"]),
                          str(end), ".", ta.strand, "."]]))
                print_file(string, nums, tss, output, 
                           out_table, srna_datas, table_best)

def detect_longer(tsss, ta, nums, output, wigs_f, wigs_r, template_texs, 
                  out_table, fuzzy, cutoff_coverage, tex_notex, replicates,
                  decrease, fuzzy_end, table_best):
    if tsss is not False:
        for tss in tsss:
            if (tss.strand == ta.strand) and \
               (tss.seq_id == ta.seq_id):
                if (tss.strand == "+"):
                    compare_ta_tss(tss.start, ta.start - fuzzy, ta.end, min_len, 
                                   Max_len, wigs_f, nums, ta, tss, output, out_table, 
                                   template_texs, detects, (ta.end - tss.start),
                                   cutoff_coverage, tex_notex, replicates)
                    if (tss.start >= ta.start - fuzzy) and \
                       (tss.start <= ta.end) and \
                       ((ta.end - tss.start) > Max_len):
                        if len(wigs_f) != 0:
                            detect_wig_pos(wigs_f, ta, tss.start, ta.end, nums, output, 
                                           "True", template_texs, out_table, cutoff_coverage,
                                           min_len, Max_len, decrease, fuzzy_end, table_best)
                else:
                    compare_ta_tss(tss.end, ta.start, ta.end + fuzzy, min_len, 
                                   Max_len, wigs_r, nums, ta, tss, output, out_table, 
                                   template_texs, detects, (tss.end - ta.start),
                                   cutoff_coverage, tex_notex, replicates)
                    if (tss.end >= ta.start) and \
                       (tss.end <= ta.end + fuzzy) and \
                       (tss.end - ta.start > Max_len):
                        if len(wigs_r) != 0:
                            detect_wig_pos(wigs_r, ta, ta.start, tss.start, nums, output, 
                                           "TSS_" + str(tss.start) + tss.strand, template_texs, 
                                           out_table, cutoff_coverage, min_len, Max_len,
                                           decrease, fuzzy_end, table_best)
    if tsss is False:
        if (len(wigs_f) != 0) and (len(wigs_r) != 0):
            if ta.strand == "+":
                num_uni = detect_wig_pos(wigs_f, ta, ta.start, ta.end, nums, output, "NA", 
                          template_texs, out_table, cutoff_coverage, min_len, Max_len,
                          decrease, fuzzy_end, table_best)
            else:
                num_uni = detect_wig_pos(wigs_r, ta, ta.start, ta.end, nums, output, "NA", 
                          template_texs, out_table, cutoff_coverage, min_len, Max_len,
                          decrease, fuzzy_end, table_best)

def compare_ta_tss(tss_pos, ta_start, ta_end, min_len, Max_len, wigs, nums, ta, tss, 
                   output, out_table, template_texs, detects, diff, cutoff_coverage, 
                   tex_notex, replicates, decrease, fuzzy_end, table_best):
    if (tss_pos >= ta_start) and \
       (tss_pos <= ta_end) and \
       (diff >= min_len) and \
       (diff <= Max_len):
        if tss.strand == "+":
            start = tss_pos
            end = ta_end
        else:
            start = ta_start
            end = tss_pos
        if len(wigs) != 0:
            get_coverage(start, end, ta.seq_id, wigs, tss.strand, ta,
                         nums, "TSS_" + str(tss.start) + tss.strand, 
                         output, template_texs, out_table, cutoff_coverage,
                         tex_notex, replicates, decrease, fuzzy_end, table_best)
        else:
            string = "\t".join([str(field) for field in [
                     ta.seq_id, ta.source, "sRNA", str(start),
                     str(end), ta.score, ta.strand, ta.phase]])
            print_file(string, nums, "TSS_" + str(tss.start) + tss.strand, 
                       output, out_table, None, table_best)
        if detects is not None:
            detects["uni_with_tss"] = True

def detect_include_TSS(num_uni, tsss, ta, wigs_f, wigs_r, output,out_table, 
                       template_texs, min_len, Max_len, detects, fuzzy, cutoff_coverage, 
                       tex_notex, replicates, decrease, fuzzy_end, table_best):
    detects["uni_with_tss"] = False
    for tss in tsss:
        if (tss.strand == ta.strand) and \
           (tss.seq_id == ta.seq_id):
            if (tss.strand == "+"):
                compare_ta_tss(tss.start, ta.start - fuzzy, ta.end, min_len, 
                               Max_len, wigs_f, nums, ta, tss, output, out_table, 
                               template_texs, detects, (ta.end - tss.start),
                               cutoff_coverage, tex_notex, replicates)
                if tss.start > ta.end:
                    break
            else:
                compare_ta_tss(tss.end, ta.start, ta.end + fuzzy, min_len, 
                               Max_len, wigs_r, nums, ta, tss, output, out_table, 
                               template_texs, detects, (tss.end - ta.start),
                               cutoff_coverage, tex_notex, replicates)
                if tss.end > ta.end + fuzzy:
                    break
    if detects["uni_with_tss"] is not True:
        if (ta.strand == "+") and (len(wigs_f) != 0):
            get_coverage(ta.start, ta.end, ta.seq_id, wigs_f, "+", ta, nums, 
                         "False", output, template_texs, out_table, cutoff_coverage, 
                         tex_notex, replicates, decrease, fuzzy_end, table_best)
        elif (ta.strand == "-") and (len(wigs_r) != 0):
            get_coverage(ta.start, ta.end, ta.seq_id, wigs_r, "-", ta, nums, 
                         "False", output, template_texs, out_table, cutoff_coverage, 
                         tex_notex, replicates, decrease, fuzzy_end, table_best)
        elif (len(wigs_f) == 0) and (len(wigs_r) == 0):
            print_file(ta.info_without_attributes.replace("Transcript", "sRNA"),
                       nums, "False", output, out_table, None, table_best)

def read_data(cdss, tsss, tas, wigs_f, wigs_r, gff_file, TSS_file, Tran_file):
    num_cds = 0
    num_tss = 0
    num_ta = 0
    gff_parser = gff3.Gff3Parser()
    g_f = open(gff_file, "r")
    for entry in gff_parser.entries(g_f):
        if (entry.feature == "CDS") or \
           (entry.feature == "pCDS") or \
           (entry.feature == "tRNA") or \
           (entry.feature == "rRNA"):
            cdss.append(entry)
            num_cds += 1
    if TSS_file is not False:
        tss_f = open(TSS_file, "r")
        for entry in gff_parser.entries(tss_f):
            tsss.append(entry)
            num_tss += 1
    for entry_ta in gff_parser.entries(open(Tran_file)):
        tas.append(entry_ta)
        num_ta += 1
    nums = {"cds": num_cds, "tss": num_tss, "ta": num_ta, "uni": 0}
    return nums

def compare_ta_cds(cdss, ta, detects):
    for cds in cdss:
        if (cds.strand == ta.strand) and \
           (cds.seq_id == ta.seq_id):
            if (cds.end < ta.end) and \
               (cds.end > ta.start) and \
               (cds.start <= ta.start):
                detects["overlap"] = True
                break
            elif (cds.start > ta.start) and \
                 (cds.start < ta.end) and \
                 (cds.end >= ta.end):
                detects["overlap"] = True
                break
            elif (cds.end >= ta.end) and \
                 (cds.start <= ta.start):
                detects["overlap"] = True
                break
            elif (cds.end <= ta.end) and \
                 (cds.start >= ta.start):
                detects["overlap"] = True
                break

def check_sRNA_condition(ta, min_len, Max_len, tsss, wigs_f, wigs_r, 
                         nums, output, out_table, texs, fuzzy, detects,
                         cutoff_coverage, tex_notex, replicates,
                         decrease, fuzzy_end, table_best):
    if ((ta.end - ta.start) >= min_len) and \
       ((ta.end - ta.start) <= Max_len):
        if len(tsss) != 0:
            detect_include_TSS(nums, tsss, ta, wigs_f, wigs_r, output, 
                               out_table, texs, min_len, Max_len, detects, fuzzy)
        else:
            if (ta.strand == "+") and (len(wigs_f) != 0):
                get_coverage(ta.start, ta.end, ta.seq_id, wigs_f, "+", 
                             ta, nums, "NA", output, texs, out_table,
                             cutoff_coverage, tex_notex, replicates,
                             decrease, fuzzy_end, table_best)
            elif (ta.strand == "-") and (len(wigs_r) != 0):
                get_coverage(ta.start, ta.end, ta.seq_id, wigs_r, "-", 
                             ta, nums, "NA", output, texs, out_table,
                             cutoff_coverage, tex_notex, replicates,
                             decrease,fuzzy_end, table_best)
            if (len(wigs_f) == 0) and (len(wigs_r) == 0):
                print_file(ta.info_without_attributes.replace("Transcript", "sRNA"),
                           nums, "NA", output, out_table, None, table_best)
    if ((ta.end - ta.start) > Max_len):
        if (len(tsss) != 0) and (len(wigs_f) == 0) and (len(wigs_r) == 0):
            detect_longer(tsss, ta, nums, output, False, False, texs, 
                          out_table, fuzzy, cutoff_coverage, tex_notex, 
                          replicates, decrease, fuzzy_end, table_best)
        elif (len(tsss) != 0) and (len(wigs_f) != 0) and (len(wigs_r) != 0):
            detect_longer(tsss, ta, nums, output, wigs_f, wigs_r, texs, 
                          out_table, fuzzy, cutoff_coverage, tex_notex, 
                          replicates, decrease, fuzzy_end, table_best)
        elif (len(tsss) == 0) and (len(wigs_f) != 0) and (len(wigs_r) != 0):
            detect_longer(False, ta, nums, output, wigs_f, wigs_r, texs, 
                          out_table, fuzzy, cutoff_coverage, tex_notex, 
                          replicates, decrease, fuzzy_end, table_best)

def Intergenic_sRNA(gff_file, Tran_file, TSS_file, fuzzy, Max_len, min_len,
                    wig_f_file, wig_r_file, wig_folder, input_libs, tex_notex,
                    replicates, output_file, output_table, table_best,
                    decrease, fuzzy_end, cutoff_coverage):
    wigs_f = {}
    wigs_r = {}
    cdss = []
    tas = []
    tsss = []
    libs = []
    texs = {}
    Read_libs(libs, texs, input_libs, wig_folder)  
    Read_wig(wig_f_file, wigs_f, "+", libs)
    Read_wig(wig_r_file, wigs_r, "-", libs)
    nums = read_data(cdss, tsss, tas, wigs_f, wigs_r, gff_file, TSS_file, Tran_file)
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start))
    tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start))
    tas = sorted(tas, key=lambda k: (k.seq_id, k.start))
    detects = {"overlap": False, "uni_with_tss": False}
    output = open(output_file, "w")        
    out_table = open(output_table, "w")
    output.write("##gff-version 3\n")
    for ta in tas:
        compare_ta_cds(cdss, ta, detects)
        if detects["overlap"]:
            detects["overlap"] = False
            continue
        else:
            check_sRNA_condidtion(ta, min_len, Max_len, tsss, wigs_f, wigs_r, 
                                  nums, output, out_table, texs, fuzzy, detects,
                                  cutoff_coverage, tex_notex, replicates,
                                  decrease, fuzzy_end, table_best)
    file_name = output_file.split(".")
    file_name = file_name[0] + ".stat"
    stat = open(file_name, "w")
    stat.write("number of cds = %s \n" % str(num_cds))
    stat.write("number of transcript assembly = %s \n" % str(num_ta))
    stat.write("total of sRNA candidates = %s \n" % str(num_uni))
