#!/usr/bin/python

import os	
import sys
import gff3
import csv
import argparse
import parser_wig as par
import math
from lib_reader import read_wig, read_libs
from coverage_detection import coverage_comparison, replicate_comparison

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-g","--gff_file",help="input gff file")
parser.add_argument("-a","--ta_file",help="input TA file")
parser.add_argument("-t","--tss_file",help="input TSS file")
parser.add_argument("-wf","--wig_f",help="input wig file (forward)")
parser.add_argument("-wr","--wig_r",help="input wig file (reverse)")
parser.add_argument("-p","--pro_file", default=False, help="input processing site file")
parser.add_argument("-M","--Max_len", default=500,type=int,help="Max length")
parser.add_argument("-m","--min_len",default=30, type=int,help="min length")
parser.add_argument("-f","--fuzzy",default=2, type=int,help="the fuzzy between TSS and transcript")
parser.add_argument("-fe","--fuzzy_end",default=10, type=int,help="the fuzzy between end of 5utr and cds")
parser.add_argument("-s","--seq",help="fasta file")
parser.add_argument("-b","--wig_folder", help="input wig folder")
parser.add_argument("-l","--libs", nargs="+", help="library name")
parser.add_argument("-te","--tex_notex", default=2, type=int, help="tex +/- should be both detected or just one.(1/2)")
parser.add_argument("-r","--replicates", type=int, help="replicate match")
parser.add_argument("-o","--output_file", help="output sRNA gff file")
parser.add_argument("-ot","--output_table", help="output sRNA table")
parser.add_argument("-tb","--table_best", default=False, action="store_true", help="output sRNA table only best track")
parser.add_argument("-d","--decrease", type=float, default=0.5, help="the range of coverage decreasing for 5'utr to determine the end of sRNA.")
parser.add_argument("-cf3","--utr3_coverage", default="median", help="cutoff of coverage")
parser.add_argument("-cf5","--utr5_coverage", default="median", help="cutoff of coverage")
parser.add_argument("-cfc","--interCDS_coverage", default="median", help="cutoff of coverage")
args = parser.parse_args()
err = open("err", "w")

def import_data(strand, strain, start, end, utr, type_, name, srna_cover):
    if type_ == "TSS":
         return {"strand": strand, "strain": strain, "start": start, 
                 "end": end, "utr": utr, "tss": name, "cleavage": "NA", 
                 "datas": srna_cover}
    elif type_ == "cleavage":
        return {"strand": strand, "strain": strain, "start": start, "end": end, 
                "utr": utr, "tss": "NA", "cleavage": name, "datas": srna_cover}
    else:
        return {"strand": strand, "strain": strain, "start": start, "end": end, 
                "utr": utr, "tss": "NA", "cleavage": "NA", "datas": srna_cover}

def get_terminal(cdss, inters, seq, type_):
    first_p = True
    first_m = True
    pre_strain = ""
    if type_ == "start":
        for cds in cdss:
            if cds.seq_id != pre_strain:
                pre_strain = cds.seq_id
                first_p = True
                first_m = True
            if (cds.strand == "+") and (first_p):
                first_p = False
                inters.append(import_data(cds.strand, cds.seq_id, 1, cds.start, 
                              "", "Term", "NA", None))
            elif (cds.strand == "-") and (first_m):
                first_m = False
                inters.append(import_data(cds.strand, cds.seq_id, 1, cds.start, 
                              "", "Term", "NA", None))
    elif type_ == "end":
        for cds in reversed(cdss):
            if cds.seq_id != pre_strain:
                pre_strain = cds.seq_id
                first_p = True
                first_m = True
            if (cds.strand == "+") and (first_p):
                first_p = False
                inters.append(import_data(cds.strand, cds.seq_id, cds.end, len(seq[cds.seq_id]), 
                              "", "Term", "NA", None))
            elif (cds.strand == "-") and (first_m):
                first_m = False
                inters.append(import_data(cds.strand, cds.seq_id, cds.start, len(seq[cds.seq_id]), 
                              "", "Term", "NA", None))

def set_cutoff(utr_type, coverages, median):
    if coverages[utr_type] == "median":
        cutoff = median
    else:
        cutoff = float(coverages[utr_type])
    return cutoff

def exchange_pos(tmp_poss, cover):
    if tmp_poss["start"] > cover["final_start"]:
        tmp_poss["start"] = cover["final_start"]
    if tmp_poss["end"] < cover["final_end"]:
        tmp_poss["end"] = cover["final_end"]

def check_tex(template_texs, covers, utr_type, coverages, median, 
              detect_num, tex_notex, srna_datas, tmp_poss):
    check_texs = {}
    texs = template_texs.copy()
    for key, num in texs.items():
        check_texs[key] = []
    for cover in covers:
        cutoff = set_cutoff(utr_type, coverages, median)
        if (cover["avg"] > cutoff):
            if (cover["type"] == "tex") or (cover["type"] == "notex"):
                for key, num in texs.items():
                    if cover["track"] in key:
                        texs[key] += 1
                        check_texs[key].append(cover)
                    if texs[key] == tex_notex:
                        if detect_num == 0:
                            tmp_poss["start"] = cover["final_start"]
                            tmp_poss["end"] = cover["final_end"]
                        else:
                            exchange_pos(tmp_poss, cover)
                        detect_num += 1
                        srna_datas["detail"].append(cover)
                        if tex_notex != 1:
                            srna_datas["detail"].append(check_texs[key][0])
                            exchange_pos(tmp_poss, check_texs[key][0])
            elif cover["type"] == "frag":
                if detect_num == 0:
                    tmp_poss["start"] = cover["final_start"]
                    tmp_poss["end"] = cover["final_end"]
                else:
                    exchange_pos(tmp_poss, cover)
                detect_num += 1
                srna_datas["detail"].append(cover)
    return detect_num

def compare_replicates(srna_tracks, template_texs, median, utr_type, coverages):
    srna_datas = {"best": 0, "high": 0, "low": 0, "start": -1,
                  "end": -1, "track": "", "detail": [], "conds": {}}
    tmp_poss = {"start": -1, "end": -1}
    all_poss = {"starts": [], "ends": []}
    for cond, covers in srna_tracks.items():
        detect_num = 0
        detect_num = check_tex(template_texs, covers, utr_type, coverages, median, 
                               detect_num, args.tex_notex, srna_datas, tmp_poss)
        if detect_num >= args.replicates:
            all_poss["starts"].append(tmp_poss["start"])
            all_poss["ends"].append(tmp_poss["end"])
            sort_datas = sorted(srna_datas["detail"], key=lambda k: (k["avg"]))
            avg = sort_datas[-1]["avg"]
            srna_datas["conds"][cond] = str(detect_num)
            if (avg > srna_datas["best"]):
                srna_datas["high"] = sort_datas[-1]["high"]
                srna_datas["low"] = sort_datas[-1]["low"]
                srna_datas["best"] = avg
                srna_datas["track"] = sort_datas[-1]["track"]
    if len(all_poss["starts"]) != 0:
        srna_datas["start"] = min(all_poss["starts"])
        srna_datas["end"] = max(all_poss["ends"])
    else:
        srna_datas["start"] = -1
        srna_datas["end"] = -1
    return srna_datas

def get_coverage(wigs, inter, start, end, type_, interCDS_type):
    srna_finals = ""
    srna_covers = {}
    srna_tracks = []
    cover_sets = {"high": 0, "low": 0, "best": -1}
    poss = {"low": 0, "high": 0}
    conditions = {}
    for wig_strain, conds in wigs.items():
        if wig_strain == inter["strain"]:
            for cond, tracks in conds.items():
                srna_covers[cond] = []
                for track, covers in tracks.items():
                    cover_tmp = {"5utr": 0, "total": 0}
                    checks = {"first": True, "detect_5utr": False}
                    final_poss = {"start": start, "end": end}
                    num = 0
                    if inter["strand"] == "-":
                        covers = reversed(covers[start - args.fuzzy - 2 - args.fuzzy_end: end + args.fuzzy + 1])
                    elif inter["strand"] == "+":
                        covers = covers[start - 2: end + 1 + args.fuzzy_end]
                    for cover in covers:
                        cover_tmp["total"] = cover_tmp["total"] + cover["coverage"]
                        checks["first"] = coverage_comparison(cover, cover_sets, 
                                          poss, checks["first"], inter["strand"])
                        if (checks["first"] is not True) and (cover_sets["high"] > 0):
                            if (type_ == "5utr") or ((type_ == "interCDS") and (interCDS_type == "TSS")):
                                if ((cover_sets["low"] / cover_sets["high"]) < args.decrease) and \
                                   (cover_sets["low"] > -1):
                                    checks["detect_5utr"] = True
                                    cover_tmp["5utr"] = cover["coverage"]
                        if checks["detect_5utr"]:
                            datas = get_cover_5utr(num, args.fuzzy_end, cover_sets, cover_tmp,
                                                   final_poss, args.decrease, cover, inter)
                            num = datas[0]
                            if datas[1] is True:
                                break
                    if (checks["first"] is not True) and (cover_sets["high"] > 0):
                        check_import_srna_covers(checks, final_poss, end, start, type_, cover_sets, cover_tmp,
                                                 interCDS_type, inter, srna_covers, cond, track, cover)
    return (srna_covers)

def get_cover_5utr(num, fuzzy_end, cover_sets, cover_tmp, final_poss, decrease, cover, inter):
    go_out = False
    if (num == fuzzy_end) or \
       (cover_tmp["5utr"] == 0) or \
       ((cover["coverage"] > cover_tmp["5utr"]) and \
        (cover["coverage"] / cover_tmp["5utr"]) > (1 + decrease)):
        if inter["strand"] == "+":
            final_poss["end"] = cover["pos"]
        elif inter["strand"] == "-":
            final_poss["start"] = cover["pos"]
        go_out = True
    elif (cover["coverage"] <= cover_tmp["5utr"]):
        if (cover["coverage"] / cover_tmp["5utr"]) >= (decrease / 2):
            num += 1
        else:
            num = 0
        cover_tmp["5utr"] = cover["coverage"]
        cover_sets["low"] = cover["coverage"]
    elif (cover["coverage"] > cover_tmp["5utr"]) and \
         ((cover["coverage"] / cover_tmp["5utr"]) <= (1 + decrease)):
        num += 1
    return (num, go_out)

def check_import_srna_covers(checks, final_poss, end, start, type_, cover_sets, cover_tmp,
                             interCDS_type, inter, srna_covers, cond, track, cover):
    avg = cover_tmp["total"] / float(end - start + 1)
    if ((type_ == "5utr") and (checks["detect_5utr"] is True)) or \
       ((type_ == "interCDS") and (interCDS_type == "TSS") and \
        (checks["detect_5utr"] is True)) or \
       ((type_ == "interCDS") and (interCDS_type == "two_pro")) or \
       (type_ == "3utr") or (type_ == "inter"):
        if final_poss["start"] < final_poss["end"]:
            if (inter["strand"] == "+") and (final_poss["end"] > end):
                final_poss["end"] = end
            elif (inter["strand"] == "-") and (final_poss["start"] < start):
                final_poss["start"] = start
            srna_covers[cond].append({"track": track, "high": cover_sets["high"],
                                      "low": cover_sets["low"], "avg": avg,
                                      "type": cover["type"],
                                      "final_start": final_poss["start"],
                                      "final_end": final_poss["end"]})

#def get_coverage(wigs, inter, start, end, type_, interCDS_type):
#    srna_finals = ""
#    srna_covers = {}
#    srna_tracks = []
#    cover_sets = {"high": 0, "low": 0, "best": -1, "total": 0, "5utr": 0}
#    poss = {"low": 0, "high": 0, "start": 0, "end": 0}
#    conditions = {}
#    for wig_strain, conds in wigs.items():
#        if wig_strain == inter["strain"]:
#            for cond, tracks in conds.items():
#                srna_covers[cond] = []
#                for track, covers in tracks.items():
#                    cover_sets["total"] = 0
#                    cover_sets["5utr"] = 0
#                    checks = {"first": True, "detect_5utr": False}
#                    num = 0
#                    poss["start"] = start
#                    poss["end"] = end
#                    if inter["strand"] == "-":
#                        covers = reversed(covers[start - args.fuzzy - 2 - args.fuzzy_end: end + args.fuzzy + 1])
#                    elif inter["strand"] == "+":
#                        covers = covers[start - 2: end + 1 + args.fuzzy_end]
#                    for cover in covers:
#                        cover_sets["total"] = cover_sets["total"] + cover["coverage"]
#                        checks["first"] = Coverage_comparison(cover, cover_sets, 
#                                          poss, checks["first"], inter["strand"])
#                        if (checks["first"] is not True) and (cover_sets["high"] > 0):
#                            if (type_ == "5utr") or ((type_ == "interCDS") and \
#                               (interCDS_type == "TSS")):
#                                if ((cover_sets["low"] / cover_sets["high"]) < args.decrease) and \
#                                   (cover_sets["low"] > -1):
#                                    checks["detect_5utr"] = True
#                                    cover_sets["5utr"] = cover["coverage"]
#                        if checks["detect_5utr"]:
#                            if (num == args.fuzzy_end) or \
#                               (cover_sets["5utr"] == 0) or \
#                               ((cover["coverage"] > cover_sets["5utr"]) and \
#                                (cover["coverage"] / cover_sets["5utr"]) > (1 + args.decrease)):
#                                if inter["strand"] == "+":
#                                    poss["end"] = cover["pos"]
#                                elif inter["strand"] == "-":
#                                    poss["start"] = cover["pos"]
#                                break
#                            elif (cover["coverage"] <= cover_sets["5utr"]):
#                                if (cover["coverage"] / cover_sets["5utr"]) >= (args.decrease / 2):
#                                    num += 1
#                                else:
#                                    num = 0
#                                cover_sets["5utr"] = cover["coverage"]
#                                cover_sets["low"] = cover["coverage"]
#                            elif (cover["coverage"] > cover_sets["5utr"]) and \
#                                 ((cover["coverage"] / cover_sets["5utr"]) <= (1 + args.decrease)):
#                                num += 1
#                    if (checks["first"] is not True) and \
#                       (cover_sets["high"] > 0):
#                        avg = cover_sets["total"] / float(end - start + 1)
#                        if ((type_ == "5utr") and (checks["detect_5utr"] is True)) or \
#                           ((type_ == "interCDS") and (interCDS_type == "TSS") and (checks["detect_5utr"] is True)) or \
#                           ((type_ == "interCDS") and (interCDS_type == "two_pro")) or \
#                           (type_ == "3utr") or (type_ == "inter"):
#                            if poss["start"] < poss["end"]:
#                                if (inter["strand"] == "+") and (poss["end"] > end):
#                                    poss["end"] = end
#                                elif (inter["strand"] == "-") and (poss["start"] < start):
#                                    poss["start"] = start
#                                srna_covers[cond].append({"track": track, "high": cover_sets["high"],
#                                                          "low": cover_sets["low"], "avg": avg,
#                                                          "type": cover["type"], 
#                                                          "final_start": poss["start"],
#                                                          "final_end": poss["end"]})
#    return (srna_covers)

def detect_3utr_pro(pros, inter, srnas, start, end, utrs, wigs, feature, name, utr_type):
    for pro in pros:
        if (pro.seq_id == inter["strain"]) and \
           (pro.strand == inter["strand"]):
            if (pro.start >= start) and \
               (pro.start <= end):
                if pro.strand == "+":
                    if ((end - pro.start) >= args.min_len) and \
                       ((end - pro.start) <= args.Max_len):
                        covers = get_coverage(wigs, inter, pro.start, end, utr_type, "NA")
                        utrs.append(import_data(inter["strand"], inter["strain"], pro.start, end,
                                    utr_type, feature, name, covers))
                        srnas.append(import_data(inter["strand"], inter["strain"], pro.start, end, 
                                     "3utr", "cleavage", "Cleavage_" + str(pro.start) + pro.strand, covers))
                else:
                    if ((pro.start - start) >= args.min_len) and \
                       ((pro.start - start) <= args.Max_len):
                        covers = get_coverage(wigs, inter, start, pro.start, utr_type, "NA")
                        utrs.append(import_data(inter["strand"], inter["strain"], start,
                                    pro.start, utr_type, feature, name, covers))
                        srnas.append(import_data(inter["strand"], inter["strain"], start, pro.start, 
                                     "3utr", "cleavage", "Cleavage_" + str(pro.start) + pro.strand, covers))
            if (pro.start > end + args.fuzzy):
                break

def detect_interCDS_twopro(pros, inter, srnas, start, end, utrs, wigs, feature, name, utr_type):
    poss = []
    for pro in pros:
        if (pro.seq_id == inter["strain"]) and \
           (pro.strand == inter["strand"]):
            if (pro.start >= start) and \
               (pro.start <= end):
                if len(poss) == 0:
                    poss.append(pro.start)
                else:
                    for pos in poss:
                        if ((pro.start - pos) >= args.min_len) and \
                           ((pro.start - pos) <= args.Max_len):
                            covers = get_coverage(wigs, inter, pos, pro.start, utr_type, "two_pro")
                            utrs.append(import_data(inter["strand"], inter["strain"], pos, pro.start,
                                        utr_type, feature, name, covers))
                            srnas.append(import_data(inter["strand"], inter["strain"], pos, pro.start,
                                 "interCDS", "cleavage", 
                                 "&".join(["Cleavage_" + str(pos) + pro.strand, \
                                           "Cleavage_" + str(pro.start) + pro.strand]), covers))
            if (pro.start > end + args.fuzzy):
                break

def detect_normal(diff, wigs, inter, start, end, utr_type, feature, 
                  name, tss, pros, utrs, srnas):
    if (diff >= args.min_len) and \
       (diff <= args.Max_len):
        covers = get_coverage(wigs, inter, start, end, utr_type, "TSS")
        utrs.append(import_data(inter["strand"], inter["strain"], start, end,
                    utr_type, feature, name, covers))
        srnas.append(import_data(inter["strand"], inter["strain"], start,
                     end, utr_type, "TSS", "TSS_" + str(tss.start) + tss.strand, covers))
    elif (diff > args.Max_len) and \
         (args.pro_file is not False):
        for pro in pros:
            if (pro.seq_id == inter["strain"]) and \
               (pro.strand == inter["strand"]):
                if (pro.start >= start) and \
                   (pro.start <= end):
                    if ((pro.start - start) >= args.min_len) and \
                       ((pro.start - start) <= args.Max_len):
                        if inter["strand"] == "+":
                            covers = get_coverage(wigs, inter, tss.start, pro.start, utr_type, "TSS")
                            utrs.append(import_data(inter["strand"], inter["strain"], tss.start,
                                        pro.start, utr_type, feature, name, covers))
                            srnas.append(import_data(inter["strand"], inter["strain"], tss.start,
                                         pro.start, utr_type, "TSS", "TSS_" + str(tss.start) + tss.strand, covers))
                        else:
                            covers = get_coverage(wigs, inter, pro.start, tss.start, utr_type, "TSS")
                            utrs.append(import_data(inter["strand"], inter["strain"], pro.start,
                                        tss.start, utr_type, feature, name, covers))
                            srnas.append(import_data(inter["strand"], inter["strain"], pro.start,
                                         tss.start, utr_type, "TSS", "TSS_" + str(tss.start) + tss.strand, covers))

def detect_utr(srnas, tsss, pros, start, end, inter, utr_type, utrs, wigs, feature, name):
    for tss in tsss:
        if (tss.seq_id == inter["strain"]) and \
           (tss.strand == inter["strand"]):
            if (tss.strand == "+"):
                if (utr_type == "5utr") or (utr_type == "interCDS"):
                    start_fuzzy = start - args.fuzzy
                else:
                    start_fuzzy = start
                if (tss.start >= start_fuzzy) and \
                   (tss.start <= end):
                    detect_normal((end - tss.start), wigs, inter, tss.start, end, 
                                  utr_type, feature, name, tss, pros, utrs, srnas)
                elif tss.start > end:
                    break
            else:
                if (utr_type == "5utr") or (utr_type == "interCDS"):
                    end_fuzzy = end + args.fuzzy
                else:
                    end_fuzzy = end
                if (tss.start >= start) and \
                   (tss.start <= end_fuzzy):
                    detect_normal((tss.start - start), wigs, inter, start, tss.start, 
                                  utr_type, feature, name, tss, pros, utrs, srnas)
                if tss.start > end_fuzzy:
                    break
    if (utr_type == "3utr") and (args.pro_file is not False):
        detect_3utr_pro(pros, inter, srnas, start, end, utrs, wigs, feature, name, utr_type)
    if (utr_type == "interCDS") and (args.pro_file is not False):
        detect_interCDS_twopro(pros, inter, srnas, start, end, utrs, wigs, feature, name, utr_type)

def run_utr_detection(wigs, inter, start, end, utr_type, utrs, tsss, 
                      pros, srnas, feature, name):
    if end - start >= 0:
        if (utr_type == "3utr") or (utr_type == "5utr") or (utr_type == "interCDS"):
            detect_utr(srnas, tsss, pros, start, end, inter, 
                       utr_type, utrs, wigs, feature, name)
        else:
            covers = get_coverage(wigs, inter, start, end, utr_type, "NA")
            utrs.append(import_data(inter["strand"], inter["strain"], start, end, 
                        utr_type, feature, name, covers))

def class_utr(inter, ta, utrs, srnas, tsss, pros, wig_fs, wig_rs):
    if inter["strand"] == "+":
        if (inter["start"] <= ta.end) and \
           (inter["end"] >= ta.end) and \
           (ta.start <= inter["start"]):
            run_utr_detection(wig_fs, inter, inter["start"] + 1, ta.end, "3utr", utrs, 
                              tsss, pros, srnas, "NA", "NA")
        elif (inter["start"] <= ta.start) and \
             (inter["end"] >= ta.start) and \
             (ta.end >= inter["end"]):
            run_utr_detection(wig_fs, inter, ta.start, inter["end"] - 1, "5utr", utrs, 
                              tsss, pros, srnas, "NA", "NA")
        elif (inter["start"] <= ta.start) and \
             (inter["end"] >= ta.end):
            run_utr_detection(wig_fs, inter, ta.start, ta.end, "inter", utrs, 
                              tsss, pros, srnas, "NA", "NA")
        elif (inter["start"] >= ta.start) and \
             (inter["end"] <= ta.end):
            run_utr_detection(wig_fs, inter, inter["start"] + 1, inter["end"] - 1, "interCDS", utrs,
                              tsss, pros, srnas, "NA", "NA")
    else:
        if (inter["start"] <= ta.end) and \
           (inter["end"] >= ta.end) and \
           (ta.start <= inter["start"]):
            run_utr_detection(wig_rs, inter, inter["start"] + 1, ta.end, "5utr", utrs, 
                              tsss, pros, srnas, "NA", "NA")
        elif (inter["start"] <= ta.start) and \
             (inter["end"] >= ta.start) and \
             (ta.end >= inter["end"]):
            run_utr_detection(wig_rs, inter, ta.start, inter["end"] - 1, "3utr", utrs, 
                              tsss, pros, srnas, "NA", "NA")
        elif (inter["start"] <= ta.start) and \
             (inter["end"] >= ta.end):
            run_utr_detection(wig_rs, inter, ta.start, ta.end, "inter", utrs, 
                              tsss, pros, srnas, "NA", "NA")
        elif (inter["start"] >= ta.start) and \
             (inter["end"] <= ta.end):
            run_utr_detection(wig_rs, inter, inter["start"] + 1, inter["end"] - 1, "interCDS", utrs,
                              tsss, pros, srnas, "NA", "NA")

def median_score(lst):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = (lstLen - 1) // 2
    if (lstLen % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1])/2.0

def print_file(num, out_t, out, srna, start, end, srna_datas, table_best):
    name = '%0*d' % (5, num)
    out_t.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t".format(
                srna["strain"], name, start, end, srna["strand"],
                ";".join(srna_datas["conds"].keys()),
                ";".join(srna_datas["conds"].values()),
                srna_datas["best"], srna_datas["high"], srna_datas["low"]))
    attribute_string = ";".join(
        ["=".join(items) for items in [["ID", "srna_utr" + str(num)],
        ["Name", "UTR_sRNA_" + name], ["UTR_type", srna["utr"]],
        ["best_avg_coverage", str(srna_datas["best"])],
        ["best_high_coverage", str(srna_datas["high"])],
        ["best_low_coverage", str(srna_datas["low"])],
        ["with_TSS", srna["tss"]], ["with_cleavage", srna["cleavage"]]]])
    out.write("\t".join([str(field) for field in [
              srna["strain"], srna["utr"], "UTR_sRNA", str(start),
              str(end), ".", srna["strand"], ".", attribute_string]]) + "\n")
    if table_best is False:
        first = True
        for data in srna_datas["detail"]:
            if first:
                out_t.write("{0}(avg={1};high={2};low={3})".format(
                            data["track"], data["avg"], data["high"], data["low"]))
                first = False
            else:
                out_t.write(";{0}(avg={1};high={2};low={3})".format(
                            data["track"], data["avg"], data["high"], data["low"]))
    else:
        out_t.write("{0}(avg={1};high={2};low={3})".format(
                    srna_datas["track"], srna_datas["best"],
                    srna_datas["high"], srna_datas["low"]))
    out_t.write("\n")

def detect_sRNA(srnas, median, out, out_t, template_texs, coverages):
    num = 0
    first = True
    if len(srnas) != 0:
        for srna in srnas:
            if srna["strain"] in median.keys():
#                srna_datas = compare_replicates(srna["datas"], template_texs, 
#                             median[srna["strain"]][srna["utr"]], srna["utr"], coverages)
                srna_datas = replicate_comparison(srna["datas"], template_texs, srna["strand"],
                                 None, args.tex_notex, args.replicates, "sRNA_utr_derived",
                                 median[srna["strain"]][srna["utr"]], coverages, srna["utr"])
                if srna_datas["best"] != 0:
                    if (srna["utr"] == "5utr") or \
                       (srna["utr"] == "interCDS"):
                        start = srna_datas["start"]
                        end = srna_datas["end"]
                    elif srna["utr"] == "3utr":
                        start = srna["start"]
                        end = srna["end"]
                    print_file(num, out_t, out, srna, start, end, srna_datas, args.table_best)
                    num += 1
#                if srna_datas[0] != 0:
#                    if srna["utr"] == "5utr":
#                        start = srna_datas[6]
#                        end = srna_datas[7]
#                    elif srna["utr"] == "3utr":
#                        start = srna["start"]
#                        end = srna["end"]
#                    elif srna["utr"] == "interCDS":
#                        start = srna_datas[6]
#                        end = srna_datas[7]
#                    name = '%0*d' % (5, num)
#                    out_t.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % \
#                               (srna["strain"], name, start, end,
#                                srna["strand"], ";".join(srna_datas[5].keys()),
#                                ";".join(srna_datas[5].values()), 
#                                srna_datas[0], srna_datas[1], srna_datas[2]))
#                    attribute_string = ";".join(
#                        ["=".join(items) for items in [["ID", "srna_utr" + str(num)],
#                        ["Name", "UTR_sRNA_" + name], ["UTR_type", srna["utr"]],
#                        ["best_avg_coverage", str(srna_datas[0])],
#                        ["best_high_coverage", str(srna_datas[1])],
#                        ["best_low_coverage", str(srna_datas[2])],
#                        ["with_TSS", srna["tss"]], ["with_cleavage", srna["cleavage"]]]])
#                    out.write("\t".join([str(field) for field in [
#                              srna["strain"], srna["utr"], "UTR_sRNA", str(start),
#                              str(end), ".", srna["strand"], ".", attribute_string]]) + "\n")
#                    if args.table_best is False:
#                        first = True
#                        for data in srna_datas[4]:
#                            if first:
#                                out_t.write("%s(avg=%s;high=%s;low=%s)" % (data["track"], data["avg"], data["high"], data["low"]))
#                                first = False
#                            else:
#                                out_t.write(";%s(avg=%s;high=%s;low=%s)" % (data["track"], data["avg"], data["high"], data["low"]))
#                    else:
#                        out_t.write("%s(avg=%s;high=%s;low=%s)" % (srna_datas[3], srna_datas[0], srna_datas[1], srna_datas[2]))
#                    out_t.write("\n")
#                    num += 1

def read_data(cdss, tas, tsss, pros, seq):
    gff_parser = gff3.Gff3Parser()
    for entry in gff_parser.entries(open(args.gff_file)):
        if (entry.feature == "CDS") or \
           (entry.feature == "tRNA") or \
           (entry.feature == "rRNA"):
            cdss.append(entry)
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start))
    for entry in gff_parser.entries(open(args.ta_file)):
        tas.append(entry)
    tas = sorted(tas, key=lambda k: (k.seq_id, k.start))
    for entry in gff_parser.entries(open(args.tss_file)):
        tsss.append(entry)
    tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start))
    if args.pro_file is not False:
        for entry in gff_parser.entries(open(args.pro_file)):
            pros.append(entry)
    pros = sorted(pros, key=lambda k: (k.seq_id, k.start))
    with open(args.seq, "r") as s_f:
        for line in s_f:
            line = line.strip()
            if line.startswith(">"):
                strain = line[1:]
                seq[strain] = ""
            else:
                seq[strain] = seq[strain] + line

def get_utr_coverage(utrs, covers, template_texs):
    first = True
    for utr in utrs:
        best_cover = -1
        for cond, utr_covers in utr["datas"].items():
            total = 0
            detect_num = 0
            check_texs = {}
            texs = template_texs.copy()
            for key, num in texs.items():
                check_texs[key] = []
            for cover in utr_covers:
                if (cover["type"] == "tex") or (cover["type"] == "notex"):
                    for key, num in texs.items():
                        if cover["track"] in key:
                            texs[key] += 1
                            check_texs[key].append(cover)
                        if texs[key] == args.tex_notex:
                            detect_num += 1
                            total = total + ((cover["avg"] + (check_texs[key][0]["avg"])) / float(2))
                elif (cover["type"] == "frag"):
                    detect_num += 1
                    total = total + cover["avg"]
            if detect_num != 0:
                avg = total / float(detect_num)
                if (avg > best_cover):
                    best_cover = avg
        if best_cover != -1:
            if first:
                covers[utr["strain"]] = {"3utr": [], "5utr": [], "interCDS": [], "inter": [], "total": []}
                covers[utr["strain"]][utr["utr"]].append(best_cover)
                covers[utr["strain"]]["total"].append(best_cover)
                first = False
            else:
                if pre_utr["strain"] != utr["strain"]:
                    covers[utr["strain"]] = {"3utr": [], "5utr": [], "interCDS": [], "inter": [], "total": []}
                    covers[utr["strain"]][utr["utr"]].append(best_cover)
                    covers[utr["strain"]]["total"].append(best_cover)
                else:
                    covers[utr["strain"]][utr["utr"]].append(best_cover)
                    covers[utr["strain"]]["total"].append(best_cover)
            pre_utr = utr

def get_inter(cdss, inters):
    for cds1 in cdss:
        for cds2 in cdss:
            if (cds1.seq_id == cds2.seq_id) and \
               (cds1.strand == cds2.strand):
                if (cds2.start > cds1.start):
                    if cds2.start - cds1.end > 1:
                        inters.append(import_data(cds1.strand, cds1.seq_id, cds1.end,
                                      cds2.start, "", "NA", "NA", None))
                    break

def set_median(covers, mediandict):
    for strain, cover in covers.items():
        mediandict[strain] = {"3utr": 0, "5utr": 0, "interCDS": 0, "inter": 0, "total": 0}
        mediandict[strain]["3utr"] = median_score(cover["3utr"])
        mediandict[strain]["5utr"] = median_score(cover["5utr"])
        mediandict[strain]["interCDS"] = median_score(cover["interCDS"])
        mediandict[strain]["inter"] = median_score(cover["inter"])
        mediandict[strain]["total"] = median_score(cover["total"])

def main():
    cdss = []
    inters = []
    tas = []
    tsss = []
    wig_fs = {}
    wig_rs = {}
    wigs = []
    pros = []
    seq = {}
    libs = []
    texs = {}
    read_data(cdss, tas, tsss, pros, seq)
    read_libs(libs, texs, args.libs, args.wig_folder)
    read_wig(wig_fs, args.wig_f, "+", libs)
    read_wig(wig_rs, args.wig_r, "-", libs)
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start))
    tas = sorted(tas, key=lambda k: (k.seq_id, k.start))
    tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start))
    if len(pros) != 0: 
        pros = sorted(pros, key=lambda k: (k.seq_id, k.start))
    out = open(args.output_file, "w")
    out.write("##gff-version 3\n")
    out_t = open(args.output_table, "w")
    get_terminal(cdss, inters, seq, "start")
    pre_strain = ""
    get_inter(cdss, inters)
    get_terminal(cdss, inters, seq, "end")
    inters = sorted(inters, key=lambda k: (k["strain"], k["start"]))
    utrs = []
    srnas = []
    for inter in inters:
        for ta in tas:
            if (inter["strain"] == ta.seq_id) and \
               (inter["strand"] == ta.strand):
                class_utr(inter, ta, utrs, srnas, tsss, pros, wig_fs, wig_rs)
    for utr in utrs:
        print(utr)
    covers = {}
    get_utr_coverage(utrs, covers, texs)
    mediandict = {}
    set_median(covers, mediandict)
    out_m = open("median", "w")
    for strain, types in mediandict.items():
        for type_, cover in types.items():
            out_m.write(";".join([strain, type_, str(cover)]) + "\n")
    if (len(mediandict) != 0):
        coverages = {"3utr": args.utr3_coverage, "5utr": args.utr5_coverage,
                     "interCDS": args.interCDS_coverage}
        detect_sRNA(srnas, mediandict, out, out_t, texs, coverages)
if __name__ == "__main__":
    main()
