import sys
import math
import csv
import os

def coverage_comparison(cover, cover_sets, poss, first, strand):
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
            cover_sets["low"] = cover["coverage"]
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
    cutoffs = {}
    if coverages[utr_type] == "median":
        for track, values in median.items():
            cutoffs[track] = values["median"]
    elif coverages[utr_type] == "mean":
        for track, values in median.items():
            cutoffs[track] = values["mean"]
    else:
        for track, values in median.items():
            cutoffs[track] = float(coverages[utr_type])
    return cutoffs

def check_tex(template_texs, covers, cutoff, target_datas, tex_notex,
              type_, poss, median, coverages, utr_type):
    detect_num = 0
    check_texs = {}
    texs = template_texs.copy()
    for key, num in texs.items():
        check_texs[key] = []
    for cover in covers:
        run_check_tex = False
        if type_ == "sRNA_utr_derived":
            cutoffs = define_cutoff(coverages, median, utr_type)
            if cover["track"] in cutoffs.keys():
                if cover["avg"] > cutoffs[cover["track"]]:
                    run_check_tex = True
            else:
                run_check_tex = True
        elif type_ == "sORF":
            if cover["avg"] > coverages[cover["track"]]:
                run_check_tex = True
        elif (type_ == "terminator") or (type_ == "merge_sRNA"):
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
                                exchange_start_end(poss, check_texs[key][0])
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

def replicate_comparison(srna_covers, template_texs, strand, cutoff_coverage,
                         tex_notex, replicates, type_, median, coverages,
                         utr_type, frag_detect):
    srna_datas = {"best": 0, "high": 0, "low": 0, "start": -1,
                  "end": -1, "track": "", "detail": [], "conds": {}}
    tmp_poss = {"start": -1, "end": -1, "pos": -1,
                "all_start": [], "all_end": []}
    output = False
    for cond, covers in srna_covers.items():
        detect_num = check_tex(template_texs, covers, cutoff_coverage,
                               srna_datas["detail"], tex_notex,
                               type_, tmp_poss, median, coverages, utr_type)
        if (((detect_num >= replicates["frag"]) and \
            ("frag" in cond)) or (not (frag_detect))):
            output = True
        if ((detect_num >= replicates["tex"]) and \
           ("texnotex" in cond)) or \
           ((detect_num >= replicates["frag"]) and \
           ("frag" in cond)):
            if type_ == "sRNA_utr_derived":
                tmp_poss["all_start"].append(tmp_poss["start"])
                tmp_poss["all_end"].append(tmp_poss["end"])
            else:
                if strand == "+":
                    sort_datas = sorted(srna_datas["detail"],
                                        key=lambda k: (k["pos"]))
                else:
                    sort_datas = sorted(srna_datas["detail"],
                                        key=lambda k: (k["pos"]), reverse=True)
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
        if len(tmp_poss["all_start"]) != 0:
            srna_datas["start"] = min(tmp_poss["all_start"])
            srna_datas["end"] = max(tmp_poss["all_end"])
        else:
            srna_datas["start"] = -1
            srna_datas["end"] = -1
        if not output:
            srna_datas = None
            srna_datas = {"best": 0, "high": 0, "low": 0, "start": -1,
                          "end": -1, "track": "", "detail": [], "conds": {}}
    return srna_datas
