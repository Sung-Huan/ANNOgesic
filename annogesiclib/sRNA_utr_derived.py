#!/usr/bin/python

import os	
import sys
import csv
import math
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.lib_reader import read_wig, read_libs
from annogesiclib.coverage_detection import coverage_comparison, replicate_comparison

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
                inters.append(import_data(cds.strand, cds.seq_id, cds.end, 
                              len(seq[cds.seq_id]), "", "Term", "NA", None))
            elif (cds.strand == "-") and (first_m):
                first_m = False
                inters.append(import_data(cds.strand, cds.seq_id, cds.start, 
                              len(seq[cds.seq_id]), "", "Term", "NA", None))

def check_pos(cover, check_point, checks):
    if (cover["pos"] >= min(check_point["utr_start"], check_point["utr_end"])) and \
       (cover["pos"] <= max(check_point["utr_start"], check_point["utr_end"])):
        checks["utr"] = True
    if (cover["pos"] >= min(check_point["srna_start"], check_point["srna_end"])) and \
       (cover["pos"] <= max(check_point["srna_start"], check_point["srna_end"])):
        checks["srna"] = True

def set_cover_and_point(inter, covers, ori_start, ori_end, fuzzy_end, 
                        start, end, check_point):
    if inter["strand"] == "-":
        covers = reversed(covers[ori_start - 2 - \
                          fuzzy_end: ori_end + 1])
        check_point["srna_start"] = end + 1
        check_point["srna_end"] = start - 2 - fuzzy_end
        check_point["utr_start"] = ori_end
        check_point["utr_end"] = ori_start
    elif inter["strand"] == "+":
        covers = covers[ori_start - 2: ori_end + 1 + fuzzy_end]
        check_point["srna_start"] = start - 2
        check_point["srna_end"] = end + 1 + fuzzy_end
        check_point["utr_start"] = ori_start
        check_point["utr_end"] = ori_end
    return covers

def get_coverage(wigs, inter, start, end, type_, interCDS_type, fuzzy_end, decrease,
                 ori_start, ori_end):
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
                    cover_tmp = {"5utr": 0, "total": 0, "ori_total": 0}
                    checks = {"first": True, "detect_5utr": False, "srna": False, "utr": False}
                    final_poss = {"start": start, "end": end}
                    check_point = {"srna_start": 0, "srna_end": 0, "utr_start": 0, "utr_end": 0}
                    num = 0
                    covers = set_cover_and_point(inter, covers, ori_start, ori_end, 
                                                 fuzzy_end, start, end, check_point)
                    for cover in covers:
                        checks["srna"] = False
                        checks["utr"] = False
                        check_pos(cover, check_point, checks)
                        if checks["utr"]:
                            cover_tmp["ori_total"] = cover_tmp["ori_total"] + cover["coverage"]
                        if checks["srna"]:
                            cover_tmp["total"] = cover_tmp["total"] + cover["coverage"]
                            checks["first"] = coverage_comparison(cover, cover_sets, 
                                              poss, checks["first"], inter["strand"])
                            if (checks["first"] is not True) and (cover_sets["high"] > 0):
                                if (type_ == "5utr") or \
                                   ((type_ == "interCDS") and (interCDS_type == "TSS")):
                                    if ((cover_sets["low"] / cover_sets["high"]) < decrease) and \
                                       (cover_sets["low"] > -1):
                                        checks["detect_5utr"] = True
                                        cover_tmp["5utr"] = cover["coverage"]
                            if checks["detect_5utr"]:
                                datas = get_cover_5utr(num, fuzzy_end, cover_sets, cover_tmp,
                                                       final_poss, decrease, cover, inter)
                                num = datas[0]
                                if datas[1] is True:
                                    break
                    if (checks["first"] is not True) and (cover_sets["high"] > 0):
                        check_import_srna_covers(checks, final_poss, end, start, type_, 
                                        cover_sets, cover_tmp, interCDS_type, inter, 
                                        srna_covers, cond, track, cover, ori_start, ori_end)
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

def check_import_srna_covers(checks, final_poss, end, start, type_, cover_sets, 
                             cover_tmp, interCDS_type, inter, srna_covers, cond, 
                             track, cover, ori_start, ori_end):
    avg = cover_tmp["total"] / float(end - start + 1)
    ori_avg = cover_tmp["ori_total"] / float(ori_end - ori_start + 1)
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
                                      "type": cover["type"], "ori_avg": ori_avg,
                                      "final_start": final_poss["start"],
                                      "final_end": final_poss["end"]})

def detect_3utr_pro(pros, inter, srnas, start, end, utrs, wigs, feature, name, 
                    utr_type, fuzzys, fuzzy_end, decrease, min_len, Max_len):
    for pro in pros:
        if (pro.seq_id == inter["strain"]) and \
           (pro.strand == inter["strand"]):
            if (pro.start >= start) and \
               (pro.start <= end):
                if pro.strand == "+":
                    if ((end - pro.start) >= min_len) and \
                       ((end - pro.start) <= Max_len):
                        covers = get_coverage(wigs, inter, pro.start, end, 
                                              utr_type, "NA", fuzzy_end, decrease, 
                                              start, end)
                        utrs.append(import_data(inter["strand"], inter["strain"], 
                                    pro.start, end, utr_type, feature, name, covers))
                        srnas.append(import_data(inter["strand"], inter["strain"], 
                                     pro.start, end, "3utr", "cleavage", 
                                     "Cleavage:" + "_".join([str(pro.start), pro.strand]), 
                                     covers))
                else:
                    if ((pro.start - start) >= min_len) and \
                       ((pro.start - start) <= Max_len):
                        covers = get_coverage(wigs, inter, start, pro.start, utr_type, 
                                              "NA", fuzzy_end, decrease, start, end)
                        utrs.append(import_data(inter["strand"], inter["strain"], start,
                                    pro.start, utr_type, feature, name, covers))
                        srnas.append(import_data(inter["strand"], inter["strain"], start, 
                                     pro.start, "3utr", "cleavage", 
                                     "Cleavage:" + "_".join([str(pro.start), pro.strand]), 
                                     covers))
            if (pro.start > end + fuzzys[utr_type]):
                break

def detect_interCDS_twopro(pros, inter, srnas, start, end, utrs, wigs, feature, name, 
                           utr_type, fuzzy_end, decrease, min_len, Max_len):
    poss = []
    for pro in pros:
        if (pro.seq_id == inter["strain"]) and \
           (pro.strand == inter["strand"]):
            if (pro.start >= start) and \
               (pro.start <= end):
                poss.append(pro)
    first = True
    for pos in poss:
        if first:
            first = False
        else:
            if ((pos.start - pre_pos.start) >= min_len) and \
               ((pos.start - pre_pos.start) <= Max_len):
                covers = get_coverage(wigs, inter, pre_pos.start, pos.start, utr_type, 
                                      "two_pro", fuzzy_end, decrease, start, end)
                utrs.append(import_data(inter["strand"], inter["strain"], 
                            pre_pos.start, pos.start, utr_type, feature, name, covers))
                srnas.append(import_data(inter["strand"], inter["strain"], 
                             pre_pos.start, pos.start, "interCDS", "cleavage", 
                             "&".join(["Cleavage:" + "_".join([str(pre_pos.start), pos.strand]), \
                                       "Cleavage:" + "_".join([str(pos.start), pos.strand])]), 
                                       covers))
        pre_pos = pos

def detect_normal(diff, wigs, inter, start, end, utr_type, feature, name, tss, 
                  pros, utrs, srnas, fuzzy_end, decrease, min_len, Max_len,
                  ori_start, ori_end):
    if (diff >= min_len) and \
       (diff <= Max_len):
        covers = get_coverage(wigs, inter, start, end, utr_type, 
                              "TSS", fuzzy_end, decrease, ori_start, ori_end)
        utrs.append(import_data(inter["strand"], inter["strain"], start, end,
                    utr_type, feature, name, covers))
        srnas.append(import_data(inter["strand"], inter["strain"], start,
                     end, utr_type, "TSS", 
                     "TSS:" + "_".join([str(tss.start), tss.strand]), covers))
    elif (diff > Max_len) and \
         (len(pros) != 0):
        for pro in pros:
            if (pro.seq_id == inter["strain"]) and \
               (pro.strand == inter["strand"]):
                if (pro.start >= start) and \
                   (pro.start <= end):
                    if ((pro.start - start) >= min_len) and \
                       ((pro.start - start) <= Max_len):
                        if inter["strand"] == "+":
                            covers = get_coverage(wigs, inter, tss.start, pro.start, utr_type, 
                                         "TSS", fuzzy_end, decrease, ori_start, ori_end)
                            utrs.append(import_data(inter["strand"], inter["strain"], tss.start,
                                        pro.start, utr_type, feature, name, covers))
                            srnas.append(import_data(inter["strand"], inter["strain"], tss.start,
                                         pro.start, utr_type, "TSS", 
                                         "TSS:" + "_".join([str(tss.start), tss.strand]), covers))
                        else:
                            covers = get_coverage(wigs, inter, pro.start, tss.start, utr_type, 
                                         "TSS", fuzzy_end, decrease, ori_start, ori_end)
                            utrs.append(import_data(inter["strand"], inter["strain"], pro.start,
                                        tss.start, utr_type, feature, name, covers))
                            srnas.append(import_data(inter["strand"], inter["strain"], pro.start,
                                         tss.start, utr_type, "TSS", 
                                         "TSS:" + "_".join([str(tss.start), tss.strand]), covers))

def detect_utr(srnas, tsss, pros, start, end, inter, utr_type, utrs, wigs, 
               feature, name, fuzzys, fuzzy_end, decrease, min_len, Max_len):
    for tss in tsss:
        if (tss.seq_id == inter["strain"]) and \
           (tss.strand == inter["strand"]):
            if (tss.strand == "+"):
                start_fuzzy = start - fuzzys[utr_type]
                if (tss.start >= start_fuzzy) and \
                   (tss.start <= end):
                    detect_normal((end - tss.start), wigs, inter, tss.start, end, 
                                  utr_type, feature, name, tss, pros, utrs, srnas,
                                  fuzzy_end, decrease, min_len, Max_len, start, end)
                elif tss.start > end:
                    break
            else:
                end_fuzzy = end + fuzzys[utr_type]
                if (tss.start >= start) and \
                   (tss.start <= end_fuzzy):
                    detect_normal((tss.start - start), wigs, inter, start, tss.start, 
                                  utr_type, feature, name, tss, pros, utrs, srnas,
                                  fuzzy_end, decrease, min_len, Max_len, start, end)
                if tss.start > end_fuzzy:
                    break
    if (utr_type == "3utr") and (len(pros) != 0):
        detect_3utr_pro(pros, inter, srnas, start, end, utrs, wigs, feature, name, 
                        utr_type, fuzzys, fuzzy_end, decrease, min_len, Max_len)
    if (utr_type == "interCDS") and (len(pros) != 0):
        detect_interCDS_twopro(pros, inter, srnas, start, end, utrs, wigs, feature, name,
                               utr_type, fuzzy_end, decrease, min_len, Max_len)

def run_utr_detection(wigs, inter, start, end, utr_type, utrs, tsss, 
                      pros, srnas, feature, name, fuzzys, fuzzy_end, 
                      decrease, min_len, Max_len):
    if end - start >= 0:
        if (utr_type == "3utr") or (utr_type == "5utr") or (utr_type == "interCDS"):
            detect_utr(srnas, tsss, pros, start, end, inter, utr_type, utrs, 
                       wigs, feature, name, fuzzys, fuzzy_end, decrease,
                       min_len, Max_len)
        else:
            covers = get_coverage(wigs, inter, start, end, utr_type, 
                                  "NA", fuzzy_end, decrease, start, end)
            utrs.append(import_data(inter["strand"], inter["strain"], start, end, 
                        utr_type, feature, name, covers))

def class_utr(inter, ta, utrs, srnas, tsss, pros, wig_fs, wig_rs,
              fuzzys, fuzzy_end, decrease, min_len, Max_len):
    if inter["strand"] == "+":
        if (inter["start"] <= ta.end) and \
           (inter["end"] >= ta.end) and \
           (ta.start <= inter["start"]):
            run_utr_detection(wig_fs, inter, inter["start"] + 1, ta.end, "3utr", 
                              utrs, tsss, pros, srnas, "NA", "NA", fuzzys, fuzzy_end, 
                              decrease, min_len, Max_len, )
        elif (inter["start"] <= ta.start) and \
             (inter["end"] >= ta.start) and \
             (ta.end >= inter["end"]):
            run_utr_detection(wig_fs, inter, ta.start, inter["end"] - 1, "5utr", 
                              utrs, tsss, pros, srnas, "NA", "NA", fuzzys, fuzzy_end, 
                              decrease, min_len, Max_len)
        elif (inter["start"] <= ta.start) and \
             (inter["end"] >= ta.end):
            run_utr_detection(wig_fs, inter, ta.start, ta.end, "inter", utrs,
                              tsss, pros, srnas, "NA", "NA", fuzzys, fuzzy_end,
                              decrease, min_len, Max_len)
        elif (inter["start"] >= ta.start) and \
             (inter["end"] <= ta.end):
            run_utr_detection(wig_fs, inter, inter["start"] + 1, inter["end"] - 1, "interCDS", utrs,
                              tsss, pros, srnas, "NA", "NA", fuzzys, fuzzy_end,
                              decrease, min_len, Max_len)
    else:
        if (inter["start"] <= ta.end) and \
           (inter["end"] >= ta.end) and \
           (ta.start <= inter["start"]):
            run_utr_detection(wig_rs, inter, inter["start"] + 1, ta.end, "5utr", 
                              utrs, tsss, pros, srnas, "NA", "NA", fuzzys, fuzzy_end, 
                              decrease, min_len, Max_len)
        elif (inter["start"] <= ta.start) and \
             (inter["end"] >= ta.start) and \
             (ta.end >= inter["end"]):
            run_utr_detection(wig_rs, inter, ta.start, inter["end"] - 1, "3utr", 
                              utrs, tsss, pros, srnas, "NA", "NA", fuzzys, fuzzy_end, 
                              decrease, min_len, Max_len)
        elif (inter["start"] <= ta.start) and \
             (inter["end"] >= ta.end):
            run_utr_detection(wig_rs, inter, ta.start, ta.end, "inter", utrs,
                              tsss, pros, srnas, "NA", "NA", fuzzys, fuzzy_end,
                              decrease, min_len, Max_len)
        elif (inter["start"] >= ta.start) and \
             (inter["end"] <= ta.end):
            run_utr_detection(wig_rs, inter, inter["start"] + 1, inter["end"] - 1, "interCDS", utrs,
                              tsss, pros, srnas, "NA", "NA", fuzzys, fuzzy_end,
                              decrease, min_len, Max_len)

def median_score(lst):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = (lstLen - 1) // 2
    if lstLen != 0:
        if (lstLen % 2):
            return sortedLst[index]
        else:
            return (sortedLst[index] + sortedLst[index + 1])/2.0
    else:
        return 0

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
    if not table_best:
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

def detect_sRNA(srnas, median, out, out_t, template_texs, coverages, 
                tex_notex, replicates, table_best, max_len, min_len):
    num = 0
    first = True
    if len(srnas) != 0:
        for srna in srnas:
            if srna["strain"] in median.keys():
                srna_datas = replicate_comparison(srna["datas"], template_texs, srna["strand"],
                                 None, tex_notex, replicates, "sRNA_utr_derived",
                                 median[srna["strain"]][srna["utr"]], coverages, srna["utr"])
                if srna_datas["best"] != 0:
                    if (srna["utr"] == "5utr") or \
                       (srna["utr"] == "interCDS"):
                        start = srna_datas["start"]
                        end = srna_datas["end"]
                    elif srna["utr"] == "3utr":
                        start = srna["start"]
                        end = srna["end"]
                    if (math.fabs(start - end) >= min_len) and \
                       (math.fabs(start - end) <= max_len):
                        print_file(num, out_t, out, srna, start, end, srna_datas, table_best)
                        num += 1

def read_data(cdss, tas, tsss, pros, seq, gff_file, 
              ta_file, tss_file, pro_file, seq_file):
    gff_parser = Gff3Parser()
    for entry in gff_parser.entries(open(gff_file)):
        if (entry.feature == "CDS") or \
           (entry.feature == "tRNA") or \
           (entry.feature == "rRNA"):
            cdss.append(entry)
    for entry in gff_parser.entries(open(ta_file)):
        tas.append(entry)
    for entry in gff_parser.entries(open(tss_file)):
        tsss.append(entry)
    if pro_file is not None:
        for entry in gff_parser.entries(open(pro_file)):
            pros.append(entry)
    with open(seq_file, "r") as s_f:
        for line in s_f:
            line = line.strip()
            if line.startswith(">"):
                strain = line[1:]
                seq[strain] = ""
            else:
                seq[strain] = seq[strain] + line

def get_utr_coverage(utrs, covers, template_texs, tex_notex):
    first = True
    for utr in utrs:
        best_cover = -1
        avgs = []
        if utr["strain"] not in covers.keys():
            covers[utr["strain"]] = {"3utr": {}, "5utr": {}, "interCDS": {}, "inter": {}, "total": {}}
        for cond, utr_covers in utr["datas"].items():
            for cover in utr_covers:
                if cover["track"] not in covers[utr["strain"]][utr["utr"]].keys():
                    covers[utr["strain"]][utr["utr"]][cover["track"]] = []
                if cover["track"] not in covers[utr["strain"]]["total"].keys():
                    covers[utr["strain"]]["total"][cover["track"]] = []
                covers[utr["strain"]][utr["utr"]][cover["track"]].append(cover["ori_avg"])
                covers[utr["strain"]]["total"][cover["track"]].append(cover["ori_avg"])

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

def mean_score(lst):
    total = 0
    for li in lst:
        total = total + li
    if len(lst) != 0:
        return (total / len(lst))
    else:
        return 0

def set_median(covers, mediandict):
    for strain, utrs in covers.items():
#        print(strain)
        mediandict[strain] = {"3utr": {}, "5utr": {}, "interCDS": {}, "inter": {}, "total": {}}
        for utr, tracks in utrs.items():
            for track, avgs in tracks.items():
                if track not in mediandict[strain][utr].keys():
                    mediandict[strain][utr][track] = {}
                mediandict[strain][utr][track] = {"median": median_score(avgs),
                                                  "mean": mean_score(avgs)}

def utr_derived_srna(gff_file, ta_file, tss_file, wig_f_file, wig_r_file, pro_file, 
                     Max_len, min_len, fuzzys, fuzzy_end, seq_file, wig_folder, 
                     input_libs, tex_notex, replicates, output_file, output_table, 
                     table_best, decrease, utr3_coverage, utr5_coverage, interCDS_coverage):
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
    read_data(cdss, tas, tsss, pros, seq, gff_file, 
              ta_file, tss_file, pro_file, seq_file)
    read_libs(libs, texs, input_libs, wig_folder)
    read_wig(wig_fs, wig_f_file, "+", libs)
    read_wig(wig_rs, wig_r_file, "-", libs)
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start))
    tas = sorted(tas, key=lambda k: (k.seq_id, k.start))
    tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start))
    if len(pros) != 0: 
        pros = sorted(pros, key=lambda k: (k.seq_id, k.start))
    out = open(output_file, "w")
    out.write("##gff-version 3\n")
    out_t = open(output_table, "w")
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
                class_utr(inter, ta, utrs, srnas, tsss, pros, wig_fs, wig_rs,
                          fuzzys, fuzzy_end, decrease, min_len, Max_len)
    covers = {}
    get_utr_coverage(utrs, covers, texs, tex_notex)
    mediandict = {}
    set_median(covers, mediandict)
    if (len(mediandict) != 0):
        coverages = {"3utr": utr3_coverage, "5utr": utr5_coverage,
                     "interCDS": interCDS_coverage}
        detect_sRNA(srnas, mediandict, out, out_t, texs, coverages,
                    tex_notex, replicates, table_best, Max_len, min_len)
