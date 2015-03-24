#!/usr/bin/python

import os	
import sys
import csv
import math
from transaplib.gff3 import Gff3Parser
from transaplib.coverage_detection import Coverage_comparison, Repliate_comparison
from transaplib.lib_reader import Read_wig, Read_libs


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

def get_cover_5utr(num, fuzzy_end, cover_sets, poss, decrease):
    if (num == fuzzy_end) or \
       (cover_sets["5utr"] == 0) or \
       ((cover["coverage"] > cover_sets["5utr"]) and \
        (cover["coverage"] / cover_sets["5utr"]) > (1 + decrease)):
        if inter["strand"] == "+":
            poss["end"] = cover["pos"]
        elif inter["strand"] == "-":
            poss["start"] = cover["pos"]
        break
    elif (cover["coverage"] <= cover_sets["5utr"]):
        if (cover["coverage"] / cover_sets["5utr"]) >= (decrease / 2):
            num += 1
        else:
            num = 0
        cover_sets["5utr"] = cover["coverage"]
        cover_sets["low"] = cover["coverage"]
    elif (cover["coverage"] > cover_sets["5utr"]) and \
         ((cover["coverage"] / cover_sets["5utr"]) <= (1 + decrease)):
        num += 1
    return num

def check_import_srna_covers(checks, covers, poss, end, start, type_, 
                             interCDS_type, inter, srna_covers, cond, track):
    if (checks["first"] is not True) and \
       (cover_sets["high"] > 0):
        avg = cover_sets["total"] / float(end - start + 1)
        if ((type_ == "5utr") and (checks["detect_5utr"] is True)) or \
           ((type_ == "interCDS") and (interCDS_type == "TSS") and \
            (checks["detect_5utr"] is True)) or \
           ((type_ == "interCDS") and (interCDS_type == "two_pro")) or \
           (type_ == "3utr") or (type_ == "inter"):
            if poss["start"] < poss["end"]:
                if (inter["strand"] == "+") and (poss["end"] > end):
                    poss["end"] = end
                elif (inter["strand"] == "-") and (poss["start"] < start):
                    poss["start"] = start
                srna_covers[cond].append({"track": track, "high": cover_sets["high"],
                                          "low": cover_sets["low"], "avg": avg,
                                          "type": cover["type"],
                                          "final_start": poss["start"],
                                          "final_end": poss["end"]})

def get_coverage(wigs, inter, start, end, type_, interCDS_type, fuzzy, fuzzy_end, decrease):
    srna_finals = ""
    srna_covers = {}
    srna_tracks = []
    cover_sets = {"high": 0, "low": 0, "best": -1, "total": 0, "5utr": 0}
    poss = {"low": 0, "high": 0, "start": 0, "end": 0}
    conditions = {}
    for wig_strain, conds in wigs.items():
        if wig_strain == inter["strain"]:
            for cond, tracks in conds.items():
                srna_covers[cond] = []
                for track, covers in tracks.items():
                    cover_sets["total"] = 0
                    cover_sets["5utr"] = 0
                    checks = {"first": True, "detect_5utr": False}
                    num = 0
                    poss["start"] = start
                    poss["end"] = end
                    if inter["strand"] == "-":
                        covers = reversed(covers[start - fuzzy - 2 - fuzzy_end: end + fuzzy + 1])
                    elif inter["strand"] == "+":
                        covers = covers[start - 2: end + 1 + fuzzy_end]
                    for cover in covers:
                        cover_sets["total"] = cover_sets["total"] + cover["coverage"]
                        checks["first"] = Coverage_comparison(cover, cover_sets, 
                                          poss, checks["first"]i, inter["strand"])
                        if (checks["first"] is not True) and (cover_sets["high"] > 0):
                            if (type_ == "5utr") or ((type_ == "interCDS") and \
                               (interCDS_type == "TSS")):
                                if ((cover_sets["low"] / cover_sets["high"]) < decrease) and \
                                   (cover_sets["low"] > -1):
                                    checks["detect_5utr"] = True
                                    cover_sets["5utr"] = cover["coverage"]
                        if detect_5utr:
                            num = get_cover_5utr(num, fuzzy_end, cover_sets, poss, decrease)
                    check_import_srna_covers(checks, covers, poss, end, start, type_,
                                             interCDS_type, inter, srna_covers, cond, track)
    return (srna_covers)

def detect_3utr_pro(pros, inter, srnas, start, end, utrs, wigs, 
                    utr_type, Max_len, min_len, fuzzy):
    for pro in pros:
        if (pro.seq_id == inter["strain"]) and \
           (pro.strand == inter["strand"]):
            if (pro.start >= start) and \
               (pro.start <= end):
                if pro.strand == "+":
                    if ((end - pro.start) >= min_len) and \
                       ((end - pro.start) <= Max_len):
                        covers = get_coverage(wigs, inter, pro.start, end, utr_type, "NA")
                        utrs.append(import_data(inter["strand"], inter["strain"], 
                                    pro.start, end, utr_type, "NA", "NA", covers))
                        srnas.append(import_data(inter["strand"], inter["strain"], 
                                     pro.start, end, "3utr", "cleavage", 
                                     "Cleavage_" + str(pro.start) + pro.strand, covers))
                else:
                    if ((pro.start - start) >= min_len) and \
                       ((pro.start - start) <= Max_len):
                        covers = get_coverage(wigs, inter, start, pro.start, utr_type, "NA")
                        utrs.append(import_data(inter["strand"], inter["strain"], 
                                    start, pro.start, utr_type, "NA", "NA", covers))
                        srnas.append(import_data(inter["strand"], inter["strain"], 
                                     start, pro.start, "3utr", "cleavage", 
                                     "Cleavage_" + str(pro.start) + pro.strand, covers))
            if (pro.start > end + fuzzy):
                break

def detect_interCDS_twopro(pros, inter, srnas, start, end, utrs, wigs, 
                           utr_type, Max_len, min_len, fuzzy):
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
                        if ((pro.start - pos) >= min_len) and \
                           ((pro.start - pos) <= Max_len):
                            covers = get_coverage(wigs, inter, pos, pro.start, 
                                                  utr_type, "two_pro")
                            utrs.append(import_data(inter["strand"], inter["strain"], pos, 
                                        pro.start, utr_type, "NA", "NA", covers))
                            srnas.append(import_data(inter["strand"], inter["strain"], 
                                         pos, pro.start, "interCDS", "cleavage", 
                                         "&".join(["Cleavage_" + str(pos) + pro.strand, \
                                         "Cleavage_" + str(pro.start) + pro.strand]), covers))
            if (pro.start > end + fuzzy):
                break

def detect_normal_utr(utr_type, fuzzy, inter, tss, diff, min_len, 
                      Max_len, start, end, pros, utrs, wigs):
    if (utr_type == "5utr") or (utr_type == "interCDS"):
        if tss.strand == "+":
            start = start - fuzzy
            utr_start = tss.start
            utr_end = end
            pro_utr_start = utr_start
            pro_utr_end = pro.start
        else:
            end = end + fuzzy
            utr_start = start
            utr_end = tss.end
            pro_utr_start = pro.start
            pro_utr_end = utr_start
    if (tss.start >= start) and \
       (tss.start <= end):
        if (diff >= min_len) and \
           (diff <= Max_len):
            covers = get_coverage(wigs, inter, utr_start, utr_end, utr_type, "TSS")
            utrs.append(import_data(inter["strand"], inter["strain"], 
                        utr_start, utr_end, utr_type, "NA", "NA", covers))
            srnas.append(import_data(inter["strand"], inter["strain"], 
                         utr_start, utr_end, utr_type, "TSS", 
                         "TSS_" + str(tss.start) + tss.strand, covers))
        elif (diff > Max_len) and (len(pros) != 0):
            for pro in pros:
                if (pro.seq_id == inter["strain"]) and \
                   (pro.strand == inter["strand"]):
                    if (pro.start >= utr_start) and \
                       (pro.start <= utr_end):
                        if ((pro.start - utr_start) >= min_len) and \
                           ((pro.start - utr_start) <= Max_len):
                            covers = get_coverage(wigs, inter, pro_utr_start, 
                                                  pro_utr_end, utr_type, "TSS")
                            utrs.append(import_data(inter["strand"], inter["strain"], 
                                        pro_utr_start, pro_utr_end, utr_type, 
                                        "NA", "NA", covers))
                            srnas.append(import_data(inter["strand"], inter["strain"], 
                                         pro_utr_start, pro_utr_end, utr_type, "TSS", 
                                         "TSS_" + str(tss.start) + tss.strand, covers))

def detect_utrs(srnas, tsss, pros, start, end, inter, utr_type, 
                utrs, wigs, Max_len, min_len, fuzzy):
    for tss in tsss:
        if (tss.seq_id == inter["strain"]) and \
           (tss.strand == inter["strand"]):
            if (tss.strand == "+"):
                detect_normal_utr(utr_type, fuzzy, inter, tss, diff, min_len, 
                                  Max_len, start, end, pros, utrs, wigs)
                if tss.start > end:
                    break
            else:
                detect_normal_utr(utr_type, fuzzy, inter, tss, diff, min_len, 
                                  Max_len, start, end, pros, utrs, wigs)
                if tss.start > end_fuzzy:
                    break
    if (utr_type == "3utr") and (len(pros) != 0):
        detect_3utr_pro(pros, inter, srnas, start, end, utrs, 
                        wigs, utr_type, Max_len, min_len, fuzzy)
    if (utr_type == "interCDS") and (len(pros) != 0):
        detect_interCDS_twopro(pros, inter, srnas, start, end, utrs, 
                               wigs, utr_type, Max_len, min_len, fuzzy)

def run_utr_detection(wigs, inter, start, end, utr_type, utrs, tsss, 
                      pros, srnas, Max_len, min_len, fuzzy):
    if end - start >= 0:
        if (utr_type == "3utr") or \
           (utr_type == "5utr") or \
           (utr_type == "interCDS"):
            detect_utrs(srnas, tsss, pros, start, end, inter, 
                        utr_type, utrs, wigs, Max_len, min_len, fuzzy)
        else:
            covers = get_coverage(wigs, inter, start, end, utr_type, "NA")
            utrs.append(import_data(inter["strand"], inter["strain"], start, end, 
                        utr_type, "NA", "NA", covers))

def classify_utr(inter, ta, wig_fs, wig_rs, utrs. tsss, pros, srnas,
                 Max_len, min_len, fuzzy):
    if len(pros) != 0:
        pros = sorted(pros, key=lambda k: (k.seq_id, k.start))
    if inter["strand"] == "+":
        type_1 = "3utr"
        type_2 = "5utr"
        wigs = wig_fs
    else:
        type_1 = "5utr"
        type_2 = "3utr"
        wigs = wig_rs
    if (inter["start"] <= ta.end) and \
       (inter["end"] >= ta.end) and \
       (ta.start <= inter["start"]):
        run_utr_detection(wigs, inter, inter["start"] + 1, ta.end, type_1, utrs,
                          tsss, pros, srnas, Max_len, min_len, fuzzy)
    elif (inter["start"] <= ta.start) and \
         (inter["end"] >= ta.start) and \
         (ta.end >= inter["end"]):
        run_utr_detection(wigs, inter, ta.start, inter["end"] - 1, type_2, utrs,
                          tsss, pros, srnas, Max_len, min_len, fuzzy)
    elif (inter["start"] <= ta.start) and \
         (inter["end"] >= ta.end):
        run_utr_detection(wigs, inter, ta.start, ta.end, "inter", utrs,
                          tsss, pros, srnas, Max_len, min_len, fuzzy)
    elif (inter["start"] >= ta.start) and \
         (inter["end"] <= ta.end):
        run_utr_detection(wigs, inter, inter["start"] + 1, inter["end"] - 1, "interCDS", utrs,
                          tsss, pros, srnas, Max_len, min_len, fuzzy)

def median(lst):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = (lstLen - 1) // 2
    if (lstLen % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1])/2.0

def print_file(num, out_t, out, srna, start, end, srna_datas, table_best):
    name = '%0*d' % (5, num)
    out_t.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % \
               (srna["strain"], name, start, end, srna["strand"],
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
                out_t.write("%s(avg=%s;high=%s;low=%s)" % (
                            data["track"], data["avg"], data["high"], data["low"]))
                first = False
            else:
                out_t.write(";%s(avg=%s;high=%s;low=%s)" % (
                            data["track"], data["avg"], data["high"], data["low"]))
    else:
        out_t.write("%s(avg=%s;high=%s;low=%s)" % (
                    srna_datas["track"], srna_datas["best"], 
                    srna_datas["high"], srna_datas["low"]))
    out_t.write("\n")

def detect_sRNA(srnas, median, out, out_t, template_texs, tex_notex, 
                coverages, replicates, table_best):
    num = 0
    if len(srnas) != 0:
        for srna in srnas:
            if srna["strain"] in median.keys():
                srna_datas = Repliate_comparison(srna["datas"], template_texs, srna.strand,
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
                    print_file(num, out_t, out, srna, start, end, srna_datas, table_best)
                    num += 1

def read_data(cdss, tas, tsss, pros, seq, gff_file, ta_file, 
              tss_file, pro_file, seq_file):
    for entry in Gff3Parser().entries(open(gff_file)):
        if (entry.feature == "CDS") or \
           (entry.feature == "tRNA") or \
           (entry.feature == "rRNA"):
            cdss.append(entry)
    for entry in Gff3Parser().entries(open(ta_file)):
        tas.append(entry)
    for entry in Gff3Parser().entries(open(tss_file)):
        tsss.append(entry)
    if pro_file is not False:
        for entry in Gff3Parser().entries(open(pro_file)):
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
                        if texs[key] == tex_notex:
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

def get_intergenic_region(cdss, inters):
    inters = []
    for cds1 in cdss:
        for cds2 in cdss:
            if (cds1.seq_id == cds2.seq_id) and \
               (cds1.strand == cds2.strand):
                if (cds2.start > cds1.start):
                    if cds2.start - cds1.end > 1:
                        inters.append(import_data(cds1.strand, cds1.seq_id, cds1.end,
                                      cds2.start, "", "NA", "NA", None))
                    break
    inters = sorted(inters, key=lambda k: (k["strain"], k["start"]))
    return inters

def set_mediandict(covers):
    mediandict = {}
    for strain, cover in covers.items():
        mediandict[strain] = {"3utr": 0, "5utr": 0, "interCDS": 0, "inter": 0, "total": 0}
        mediandict[strain]["3utr"] = median(cover["3utr"])
        mediandict[strain]["5utr"] = median(cover["5utr"])
        mediandict[strain]["interCDS"] = median(cover["interCDS"])
        mediandict[strain]["inter"] = median(cover["inter"])
        mediandict[strain]["total"] = median(cover["total"])
    return mediandict

def UTR_derived_sRNA(gff_file, ta_file, tss_file, wig_f_file, wig_r_file,
                     pro_file, Max_len, min_len, fuzzy, fuzzy_end, seq_file,
                     wig_folder, input_libs, tex_notex, replicates, 
                     output_file, output_table, table_best, decrease, 
                     utr3_coverage, utr5_coverage, interCDS_coverage):
    cdss = []
    tas = []
    tsss = []
    wig_fs = {}
    wig_rs = {}
    pros = []
    seq = {}
    libs = []
    texs = {}
    read_data(cdss, tas, tsss, pros, seq, gff_file, ta_file,
              tss_file, pro_file, seq_file)
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start))
    tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start))
    tas = sorted(tas, key=lambda k: (k.seq_id, k.start))
    Read_libs(libs, texs, input_libs, wig_folder)
    Read_wig(wig_f_file, wig_fs, "+", libs)
    Read_wig(wig_r_file, wig_rs, "-", libs)
    out = open(output_file, "w")
    out.write("##gff-version 3\n")
    out_t = open(output_table, "w")
    get_terminal(cdss, inters, seq, "start")
    pre_strain = ""
    get_intergenic_region(cdss, inters)
    get_terminal(cdss, inters, seq, "end")
    utrs = []
    srnas = []
    for inter in inters:
        for ta in tas:
            if (inter["strain"] == ta.seq_id) and \
               (inter["strand"] == ta.strand):
                classify_utr(inter, ta, wig_fs, wig_rs, utrs. tsss, pros, srnas,
                             Max_len, min_len, fuzzy)
    covers = {}
    inters = get_utr_coverage(utrs, covers, texs, tex_notex)
    mediandict = set_mediandict(covers)
    out_m = open("median", "w")
    for strain, types in mediandict.items():
        for type_, cover in types.items():
            out_m.write(";".join([strain, type_, str(cover)]) + "\n")
    if (len(mediandict) != 0):
        coverages = {"3utr": utr3_coverage, "5utr": utr5_coverage, 
                     "interCDS": interCDS_coverage}
        detect_sRNA(srnas, mediandict, out, out_t, texs, tex_notex, 
                    coverages, replicates, table_best)
