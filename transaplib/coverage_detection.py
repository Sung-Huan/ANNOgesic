import sys
import math
import csv
import os

def Coverage_comparison(cover, cover_sets, poss, first, strand):
    if first:
        first = False
        cover_sets["high"] = cover["coverage"]
        cover_sets["low"] = cover["coverage"]
        poss["high"] = cover["pos"]
        poss["low"] = cover["pos"]
    else:
        if cover_sets["high"] < cover["coverage"]:
            cover_sets["high"] = cover["coverage"]
            poss["high"] = cover["pos"]
            poss["low"] = cover["pos"]
        if ((strand == "+") and (poss["low"] >= poss["high"])) or \
           ((strand == "-") and (poss["low"] <= poss["high"])):
            if cover_sets["low"] > cover["coverage"]:
                cover_sets["low"] = cover["coverage"]
                poss["low"] = cover["pos"]
        elif ((strand == "+") and (poss["low"] < poss["high"])) or \
             ((strand == "-") and (poss["low"] > poss["high"])):
            poss["low"] = cover["pos"]
            cover_sets["low"] = cover["coverage"]
    return first

def define_cutoff(coverages, median, utr_type):
    if coverages[utr_type] == "median":
        cutoff = median
    else:
        cutoff = float(coverages[utr_type])
    return cutoff

def Check_Tex(template_texs, covers, cutoff, target_datas, tex_notex,
              detect_num, type_, poss, median, coverages, utr_type):
    check_texs = {}
    texs = template_texs.copy()
    run_check_tex = False
    for key, num in texs.items():
        check_texs[key] = []
    for cover in covers:
        if type_ == "sRNA_utr_derived":
            cutoff = define_cutoff(coverages, median, utr_type)
        if type_ == "terminator":
            run_check_tex = True
        else:
            if cover["avg"] > cutoff:
                run_check_tex = True
        if run_check_tex:
            if (cover["type"] == "tex") or (cover["type"] == "notex"):
                for key, num in texs.items():
                    if cover["track"] in key:
                        texs[key] += 1
                        check_texs[key].append(cover)
                    if texs[key] == tex_notex:
                        if type_ == "sRNA_utr_derived":
                            if detect_num == 0:
                                poss["start"] = cover["final_start"]
                                poss["end"] = cover["final_end"]
                            else:
                                exchange_start_end(poss, cover)
                        detect_num += 1
                        target_datas.append(cover)
                        if tex_notex != 1:
                            target_datas.append(check_texs[key][0])
                            if type_ == "sRNA_utr_derived":
                                exchange_start_end(poss, cover)
            elif cover["type"] == "frag":
                if type_ == "sRNA_utr_derived":
                    if detect_num == 0:
                        poss["start"] = cover["final_start"]
                        poss["end"] = cover["final_end"]
                    else:
                        exchange_start_end(poss, cover)
                detect_num += 1
                target_datas.append(cover)
    return detect_num

def exchange_start_end(poss, cover):
    if poss["start"] > cover["final_start"]:
        poss["start"] = cover["final_start"]
    if poss["end"] < cover["final_end"]:
        poss["end"] = cover["final_end"]

def Repliate_comparison(srna_covers, template_texs, strand, cutoff_coverage,
                        tex_notex, replicates, type_, median, coverages, utr_type):
    srna_datas = {"best": 0, "high": 0, "low": 0, "pos": -1,
                  "track": "", "detail": [], "conds": {}}
    starts = [] ### for sRNA utr derived
    ends = []
    for cond, covers in srna_covers.items():
        detect_num = 0
        tmp_poss = {"start": -1, "end": -1}
        detect_num = Check_Tex(template_texs, covers, cutoff_coverage, 
                               srna_datas["detail"], tex_notex, detect_num, 
                               type_, tmp_poss, median, coverages, utr_type)
        if detect_num >= replicates:
            if strand == "+":
                sort_datas = sorted(srna_datas["detail"], key=lambda k: (k["pos"]))
            else:
                sort_datas = sorted(srna_datas["detail"], key=lambda k:
                                    (k["pos"]), reverse=True)
            if type_ == "sRNA_utr_derived":
                starts.append(tmp_poss["start"])
                ends.append(tmp_poss["end"])
            srna_datas["pos"] = sort_datas[-1]["pos"]
            sort_datas = sorted(srna_datas["detail"], key=lambda k: (k["avg"]))
            avg = sort_datas[-1]["avg"]
            srna_datas["conds"][cond] = str(detect_num)
            if (avg > srna_datas["best"]):
                srna_datas["high"] = sort_datas[-1]["high"]
                srna_datas["low"] = sort_datas[-1]["low"]
                srna_datas["best"] = avg
                srna_datas["track"] = sort_datas[-1]["track"]
    if type_ == "sRNA_utr_derived":
        if len(starts) != 0:
            srna_datas["start"] = min(starts)
            srna_datas["end"] = max(ends)
        else:
            srna_datas["start"] = -1
            srna_datas["end"] = -1
    return srna_datas
