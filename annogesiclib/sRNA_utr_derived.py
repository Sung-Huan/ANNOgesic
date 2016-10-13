import os, gc
import math
import numpy as np
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.lib_reader import read_wig, read_libs
from annogesiclib.coverage_detection import coverage_comparison, get_repmatch
from annogesiclib.coverage_detection import replicate_comparison
from annogesiclib.args_container import ArgsContainer


def import_data(strand, strain, pos, utr, type_, name, srna_cover, pro):
    if type_ == "TSS":
        data = {"strand": strand, "strain": strain, "start": pos["start"],
                "end": pos["end"], "utr": utr, "start_tss": name,
                "end_cleavage": pro, "start_cleavage": "NA",
                "datas": srna_cover}
    elif type_ == "cleavage":
        data = {"strand": strand, "strain": strain, "start": pos["start"],
                "end": pos["end"], "utr": utr, "start_tss": "NA",
                "start_cleavage": name, "datas": srna_cover,
                "end_cleavage": pro}
    else:
        data = {"strand": strand, "strain": strain, "start": pos["start"],
                "end": pos["end"], "utr": utr, "start_tss": "NA",
                "start_cleavage": "NA", "datas": srna_cover,
                "end_cleavage": pro}
    del(srna_cover)
    return data


def import_inter(strand, strain, start, end, length):
    return {"strand": strand, "strain": strain,
            "start": start, "end": end, "len_CDS": length}


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
                inters.append(import_inter(cds.strand, cds.seq_id,
                                           1, cds.start, 0))
            elif (cds.strand == "-") and (first_m):
                first_m = False
                inters.append(import_inter(cds.strand, cds.seq_id,
                                           1, cds.start, 0))
    elif type_ == "end":
        for cds in reversed(cdss):
            if cds.seq_id != pre_strain:
                pre_strain = cds.seq_id
                first_p = True
                first_m = True
            if (cds.strand == "+") and (first_p):
                first_p = False
                inters.append(import_inter(
                    cds.strand, cds.seq_id, cds.end,
                    len(seq[cds.seq_id]), cds.end - cds.start))
            elif (cds.strand == "-") and (first_m):
                first_m = False
                inters.append(import_inter(
                    cds.strand, cds.seq_id, cds.start,
                    len(seq[cds.seq_id]), cds.end - cds.start))


def check_pos(cover, check_point, checks, cover_pos):
    if (cover_pos >= min(check_point["utr_start"],
                            check_point["utr_end"])) and (
        cover_pos <= max(check_point["utr_start"],
                            check_point["utr_end"])):
        checks["utr"] = True
    if (cover_pos >= min(check_point["srna_start"],
                            check_point["srna_end"])) and (
        cover_pos <= max(check_point["srna_start"],
                            check_point["srna_end"])):
        checks["srna"] = True


def check_start_and_end(start, end, covers, strand, fuzzy):
    if strand == "-":
        if (start - 2 - fuzzy) < 0:
            c_start = 0
        else:
            c_start = start - 2 - fuzzy
        if (end + 1) > len(covers):
            c_end = len(covers)
        else:
            c_end = end + 1
    else:
        if (start - 2) < 0:
            c_start = 0
        else:
            c_start = start - 2
        if (end + 1 + fuzzy) > len(covers):
            c_end = len(covers)
        else:
            c_end = end + 1 + fuzzy
    return c_start, c_end


def set_cover_and_point(cover_results, inter, covers, pos, fuzzy_end):
    check_point = {"srna_start": 0, "srna_end": 0,
                   "utr_start": 0, "utr_end": 0}
    if pos["ori_start"] - 2 < 0:
        ori_start = 2
    else:
        ori_start = pos["ori_start"]
    if inter["strand"] == "-":
        c_start, c_end = check_start_and_end(ori_start, pos["ori_end"],
                                             covers, "-", fuzzy_end)
        covers = covers[c_start: c_end]
        covers = covers[::-1]
        check_point["srna_start"] = pos["end"] + 1
        check_point["srna_end"] = pos["start"] - 2 - fuzzy_end
        check_point["utr_start"] = pos["ori_end"]
        check_point["utr_end"] = ori_start
        start = c_start
        end = c_end
    elif inter["strand"] == "+":
        c_start, c_end = check_start_and_end(ori_start, pos["ori_end"], 
                                             covers, "+", fuzzy_end)
        covers = covers[c_start: c_end]
        check_point["srna_start"] = pos["start"] - 2
        check_point["srna_end"] = pos["end"] + 1 + fuzzy_end
        check_point["utr_start"] = ori_start
        check_point["utr_end"] = pos["ori_end"]
        start = c_start
        end = c_end
    cover_results["check_point"] = check_point
    cover_results["covers"] = covers
    return start, end


def detect_cover_utr_srna(cover_results, pos, inter, cond, track,
                          args_srna, lib_type, start, end, strand):
    datas = {"num": 0, "cover_tmp": {"5utr": 0, "total": 0, "ori_total": 0},
             "checks": {"first": True, "detect_decrease": False,
                        "srna": False, "utr": False},
             "final_poss": {"start": pos["start"], "end": pos["end"]}}
    index_pos = 0
    for cover in cover_results["covers"]:
        if strand == "+":
            cover_pos = start + index_pos
        else:
            cover_pos = end - index_pos
        datas["checks"]["srna"] = False
        datas["checks"]["utr"] = False
        check_pos(cover, cover_results["check_point"], datas["checks"],
                  cover_pos)
        if datas["checks"]["utr"]:
            datas["cover_tmp"]["ori_total"] = \
                datas["cover_tmp"]["ori_total"] + cover
        if datas["checks"]["srna"]:
            datas["cover_tmp"]["total"] = \
                datas["cover_tmp"]["total"] + cover
            datas["checks"]["first"] = coverage_comparison(
                    cover, cover_results["cover_sets"], cover_results["pos"],
                    datas["checks"]["first"], inter["strand"], cover_pos)
            if (datas["checks"]["first"] is not True) and (
                    cover_results["cover_sets"]["high"] > 0):
                if (cover_results["type"] == "5utr") or (
                        cover_results["type"] == "3utr") or (
                        (cover_results["type"] == "interCDS") and (
                         cover_results["intercds"] == "TSS")):
                    if ((cover_results["cover_sets"]["low"] /
                            cover_results["cover_sets"]["high"]) <
                            args_srna.decrease_utr) and (
                            cover_results["cover_sets"]["low"] > -1):
                        datas["checks"]["detect_decrease"] = True
                        datas["cover_tmp"]["5utr"] = cover
            if datas["checks"]["detect_decrease"]:
                go_out = get_cover_5utr(datas, cover_results["cover_sets"],
                                        cover, inter, args_srna, cover_pos)
                if go_out is True:
                    break
        index_pos += 1
    if (datas["checks"]["first"] is not True) and (
            cover_results["cover_sets"]["high"] > 0):
        check_import_srna_covers(datas, cover_results, inter, cond, track,
                                 cover, pos, args_srna, lib_type)


def get_coverage(wigs, inter, pos, type_, intercds_type, args_srna):
    cover_results = {"srna_covers": {}, "utr_covers": {},
                     "cover_sets": {"high": 0, "low": 0, "best": -1},
                     "pos": {"low": 0, "high": 0},
                     "type": type_, "intercds": intercds_type}
    for wig_strain, conds in wigs.items():
        if wig_strain == inter["strain"]:
            for cond, tracks in conds.items():
                cover_results["srna_covers"][cond] = []
                cover_results["utr_covers"][cond] = []
                for lib_name, covers in tracks.items():
                    track = lib_name.split("|")[-3]
                    lib_strand = lib_name.split("|")[-2]
                    lib_type = lib_name.split("|")[-1]
                    start, end = set_cover_and_point(
                            cover_results, inter, covers,
                            pos, args_srna.fuzzy_utr)
                    detect_cover_utr_srna(cover_results, pos, inter,
                                          cond, track, args_srna, lib_type,
                                          start, end, lib_strand)
    return cover_results["srna_covers"], cover_results["utr_covers"]


def get_cover_5utr(datas, cover_sets, cover, inter, args_srna, cover_pos):
    go_out = False
    if (datas["num"] == args_srna.fuzzy_utr) or (
            datas["cover_tmp"]["5utr"] == 0) or (
            (cover > datas["cover_tmp"]["5utr"]) and (
            cover / datas["cover_tmp"]["5utr"]) > (
            1 + args_srna.decrease_utr)):
        if inter["strand"] == "+":
            datas["final_poss"]["end"] = cover_pos
        elif inter["strand"] == "-":
            datas["final_poss"]["start"] = cover_pos
        go_out = True
    elif (cover <= datas["cover_tmp"]["5utr"]):
        if (cover / datas["cover_tmp"]["5utr"]) >= (
                args_srna.decrease_utr / 2):
            datas["num"] += 1
        else:
            datas["num"] = 0
        datas["cover_tmp"]["5utr"] = cover
        cover_sets["low"] = cover
    elif (cover > datas["cover_tmp"]["5utr"]) and (
            (cover / datas["cover_tmp"]["5utr"]) <= (
            1 + args_srna.decrease_utr)):
        datas["num"] += 1
    return go_out


def import_cover(inter, covers, track, cover_sets, avgs, cover,
                 final_poss, pos, lib_type):
    if final_poss["start"] < final_poss["end"]:
        if (inter["strand"] == "+") and (final_poss["end"] > pos["end"]):
            final_poss["end"] = pos["end"]
        elif (inter["strand"] == "-") and (final_poss["start"] < pos["start"]):
            final_poss["start"] = pos["start"]
        covers.append({"track": track,
                       "high": cover_sets["high"],
                       "low": cover_sets["low"], "avg": avgs["avg"],
                       "type": lib_type,
                       "ori_avg": avgs["ori_avg"],
                       "final_start": final_poss["start"],
                       "final_end": final_poss["end"]})


def check_import_srna_covers(datas, cover_results, inter, cond, track,
                             cover, pos, args_srna, lib_type):
    avgs = {"avg": datas["cover_tmp"]["total"] / float(
                   pos["end"] - pos["start"] + 1),
            "ori_avg": datas["cover_tmp"]["ori_total"] / float(
                       pos["ori_end"] - pos["ori_start"] + 1)}
    if ((cover_results["type"] == "5utr") and (
            (cover_results["intercds"] == "tsspro") or (
             datas["checks"]["detect_decrease"]))) or (
        (cover_results["type"] == "3utr") and (
             cover_results["intercds"] == "two_pro")) or (
        (cover_results["type"] == "3utr") and (
             cover_results["intercds"] != "two_pro") and (
             datas["checks"]["detect_decrease"])) or (
        (cover_results["type"] == "3utr") and (
             cover_results["intercds"] != "two_pro") and (
             (pos["end"] - pos["start"]) >= args_srna.min_len) and (
             (pos["end"] - pos["start"]) <= args_srna.max_len)) or (
        (cover_results["type"] == "interCDS") and (
             cover_results["intercds"] == "TSS") and (
             datas["checks"]["detect_decrease"])) or (
        (cover_results["type"] == "interCDS") and (
            (cover_results["intercds"] == "tss_pro") or (
             cover_results["intercds"] == "two_pro"))):
        import_cover(inter, cover_results["srna_covers"][cond], track,
                     cover_results["cover_sets"], avgs, cover,
                     datas["final_poss"], pos, lib_type)
        import_cover(inter, cover_results["utr_covers"][cond], track,
                     cover_results["cover_sets"], avgs, cover,
                     datas["final_poss"], pos, lib_type)
    else:
        import_cover(inter, cover_results["utr_covers"][cond], track,
                     cover_results["cover_sets"], avgs, cover,
                     datas["final_poss"], pos, lib_type)


def detect_3utr_pro(inter, pos, wigs, utr_type, args_srna):
    '''3UTR start with processing site'''
    for pro in args_srna.pros:
        if (pro.seq_id == inter["strain"]) and (
                pro.strand == inter["strand"]):
            if (pro.start >= pos["start"]) and (pro.start <= pos["end"]):
                if pro.strand == "+":
                    if ((pos["end"] - pro.start) >= args_srna.min_len) and (
                            (pos["end"] - pro.start) <= args_srna.max_len):
                        n_pos = import_position(
                                pro.start, pos["end"],
                                pos["ori_start"], pos["ori_end"])
                        srna_covers, utr_covers = get_coverage(
                                wigs, inter, n_pos, utr_type, "pro", args_srna)
                        args_srna.utrs.append(import_data(
                            inter["strand"], inter["strain"], n_pos, utr_type,
                            "NA", "NA", utr_covers, "NA"))
                        args_srna.srnas.append(import_data(
                            inter["strand"], inter["strain"], n_pos, "3utr",
                            "cleavage", "Cleavage:" + "_".join([
                                str(pro.start), pro.strand]),
                            srna_covers, "NA"))
                    elif (pos["end"] - pro.start) > args_srna.max_len:
                        detect_twopro(inter, pos, wigs, utr_type,
                                      "3utr", args_srna)
                else:
                    if ((pro.start - pos["start"]) >= args_srna.min_len) and (
                            (pro.start - pos["start"]) <= args_srna.max_len):
                        n_pos = import_position(
                                pos["start"], pro.start,
                                pos["ori_start"], pos["ori_end"])
                        srna_covers, utr_covers = get_coverage(
                                wigs, inter, n_pos, utr_type, "pro", args_srna)
                        args_srna.utrs.append(import_data(
                            inter["strand"], inter["strain"], n_pos, utr_type,
                            "NA", "NA", utr_covers, "NA"))
                        args_srna.srnas.append(import_data(
                            inter["strand"], inter["strain"], n_pos, "3utr",
                            "cleavage", "Cleavage:" + "_".join([
                                str(pro.start), pro.strand]),
                            srna_covers, "NA"))
                    elif (pro.start - pos["start"]) > args_srna.max_len:
                        detect_twopro(inter, pos, wigs, utr_type,
                                      "3utr", args_srna)
            if (pro.start > pos["end"] + args_srna.fuzzy_tsss[utr_type]):
                break


def detect_twopro(inter, pos, wigs, utr_type, import_type, args_srna):
    '''the sRNA is associated with two processing sites'''
    pros = []
    for pro in args_srna.pros:
        if (pro.seq_id == inter["strain"]) and (
                pro.strand == inter["strand"]):
            if (pro.start >= pos["start"]) and (pro.start <= pos["end"]):
                pros.append(pro)
    first = True
    pre_pro = None
    for pro in pros:
        if first:
            first = False
        else:
            if ((pro.start - pre_pro.start) >= args_srna.min_len) and (
                    (pro.start - pre_pro.start) <= args_srna.max_len):
                n_pos = import_position(pre_pro.start, pro.start,
                                        pos["ori_start"], pos["ori_end"])
                srna_covers, utr_covers = get_coverage(
                        wigs, inter, n_pos, utr_type, "two_pro", args_srna)
                args_srna.utrs.append(import_data(
                    inter["strand"], inter["strain"], n_pos, utr_type, "NA",
                    "NA", utr_covers, "Cleavage:" + "_".join(
                        [str(pro.start), pro.strand])))
                args_srna.srnas.append(import_data(
                    inter["strand"], inter["strain"], n_pos, import_type,
                    "cleavage", "Cleavage:" + "_".join(
                        [str(pre_pro.start), pro.strand]),
                    srna_covers, "Cleavage:" + "_".join(
                        [str(pro.start), pro.strand])))
        pre_pro = pro


def decrease_pos(covers, pos, strand):
    longer = -1
    for cond, datas in covers.items():
        for data in datas:
            for key, value in data.items():
                if key == pos:
                    if longer == -1:
                        longer = value
                    elif (strand == "+") and (value > longer) and (
                            longer != -1):
                        longer = value
                    elif (strand == "-") and (value < longer) and (
                            longer != -1):
                        longer = value
    return longer


def get_decrease(inter, wigs, tss, pos, utr_type, args_srna):
    '''check the coverage decrease'''
    if inter["strand"] == "+":
        n_pos = import_position(tss.start, pos["end"], pos["ori_start"],
                                pos["ori_end"])
        srna_covers, utr_covers = get_coverage(
                wigs, inter, n_pos, utr_type, "TSS", args_srna)
        utr_pos = decrease_pos(utr_covers, "final_end", "+")
        n_pos = import_position(tss.start, utr_pos, pos["ori_start"],
                                pos["ori_end"])
        args_srna.utrs.append(import_data(
            inter["strand"], inter["strain"], n_pos, utr_type, "NA", "NA",
            utr_covers, "NA"))
        if len(srna_covers) != 0:
            srna_pos = decrease_pos(srna_covers, "final_end", "+")
            n_pos = import_position(tss.start, srna_pos, pos["ori_start"],
                                    pos["ori_end"])
            args_srna.srnas.append(import_data(
                inter["strand"], inter["strain"], n_pos, utr_type, "TSS",
                "TSS:" + "_".join([str(tss.start), tss.strand]),
                srna_covers, "NA"))
    else:
        n_pos = import_position(pos["start"], tss.start, pos["ori_start"],
                                pos["ori_end"])
        srna_covers, utr_covers = get_coverage(
                wigs, inter, n_pos, utr_type, "TSS", args_srna)
        utr_pos = decrease_pos(utr_covers, "final_start", "-")
        n_pos = import_position(utr_pos, tss.start, pos["ori_start"],
                                pos["ori_end"])
        args_srna.utrs.append(import_data(
            inter["strand"], inter["strain"], n_pos, utr_type, "NA", "NA",
            utr_covers, "NA"))
        if len(srna_covers) != 0:
            srna_pos = decrease_pos(srna_covers, "final_start", "-")
            n_pos = import_position(srna_pos, tss.start, pos["ori_start"],
                                    pos["ori_end"])
            args_srna.srnas.append(import_data(
                inter["strand"], inter["strain"], n_pos, utr_type, "TSS",
                "TSS:" + "_".join([str(tss.start), tss.strand]),
                srna_covers, "NA"))


def import_append_normal(inter, tss, pro, pos, wigs, utr_type, args_srna):
    if inter["strand"] == "+":
        n_pos = import_position(tss.start, pro.start,
                                pos["ori_start"], pos["ori_end"])
        srna_covers, utr_covers = get_coverage(
                wigs, inter, n_pos, utr_type, "tsspro", args_srna)
        args_srna.utrs.append(import_data(
            inter["strand"], inter["strain"], n_pos, utr_type, "NA", "NA",
            utr_covers, "Cleavage:" + "_".join([str(pro.start), pro.strand])))
        args_srna.srnas.append(import_data(
            inter["strand"], inter["strain"], n_pos, utr_type, "TSS",
            "TSS:" + "_".join([str(tss.start), tss.strand]), srna_covers,
            "Cleavage:" + "_".join([str(pro.start), pro.strand])))
    else:
        n_pos = import_position(pro.start, tss.start,
                                pos["ori_start"], pos["ori_end"])
        srna_covers, utr_covers = get_coverage(
                wigs, inter, n_pos, utr_type, "tsspro", args_srna)
        args_srna.utrs.append(import_data(
            inter["strand"], inter["strain"], n_pos, utr_type, "NA", "NA",
            utr_covers, "Cleavage:" + "_".join([str(pro.start), pro.strand])))
        args_srna.srnas.append(import_data(
            inter["strand"], inter["strain"], n_pos, utr_type, "TSS",
            "TSS:" + "_".join([str(tss.start), tss.strand]), srna_covers,
            "Cleavage:" + "_".join([str(pro.start), pro.strand])))


def detect_normal(diff, wigs, inter, pos, utr_type, tss, args_srna):
    '''normal case, UTR-derived sRNA with TSS'''
    if (diff >= args_srna.min_len) and (
            diff <= args_srna.max_len):
        srna_covers, utr_covers = get_coverage(
                wigs, inter, pos, utr_type, "TSS", args_srna)
        args_srna.utrs.append(import_data(
            inter["strand"], inter["strain"], pos, utr_type, "NA", "NA",
            utr_covers, "NA"))
        args_srna.srnas.append(import_data(
            inter["strand"], inter["strain"], pos, utr_type, "TSS",
            "TSS:" + "_".join([str(tss.start), tss.strand]),
            srna_covers, "NA"))
    elif (diff > args_srna.max_len) and (len(args_srna.pros) != 0):
        detect = False
        for pro in args_srna.pros:
            if (pro.seq_id == inter["strain"]) and (
                    pro.strand == inter["strand"]):
                if (pro.start >= pos["start"]) and (pro.start <= pos["end"]):
                    if ((pro.start - pos["start"]) >= args_srna.min_len) and (
                            (pro.start - pos["start"]) <= args_srna.max_len):
                        detect = True
                        import_append_normal(inter, tss, pro, pos,
                                             wigs, utr_type, args_srna)
        if not detect:
            get_decrease(inter, wigs, tss, pos, utr_type, args_srna)


def import_position(start, end, ori_start, ori_end):
    return {"start": start, "end": end,
            "ori_start": ori_start, "ori_end": ori_end}


def detect_utr(start, end, inter, utr_type, wigs, args_srna):
    ori_fuzzy = args_srna.fuzzy_tsss[utr_type]
    if "p_" in args_srna.fuzzy_tsss[utr_type]:
        per = float(args_srna.fuzzy_tsss[utr_type].split("_")[-1])
        args_srna.fuzzy_tsss[utr_type] = inter["len_CDS"]*per
    elif "n_" in args_srna.fuzzy_tsss[utr_type]:
        args_srna.fuzzy_tsss[utr_type] = float(
            args_srna.fuzzy_tsss[utr_type].split("_")[-1])
    for tss in args_srna.tsss:
        if (tss.seq_id == inter["strain"]) and (
                tss.strand == inter["strand"]):
            if tss.strand == "+":
                start_fuzzy = start - args_srna.fuzzy_tsss[utr_type]
                if (tss.start >= start_fuzzy) and (tss.start <= end):
                    n_pos = import_position(tss.start, end, start, end)
                    detect_normal(end - tss.start, wigs, inter, n_pos,
                                  utr_type, tss, args_srna)
                elif tss.start > end:
                    break
            else:
                end_fuzzy = end + args_srna.fuzzy_tsss[utr_type]
                if (tss.start >= start) and (tss.start <= end_fuzzy):
                    n_pos = import_position(start, tss.start, start, end)
                    detect_normal(tss.start - start, wigs, inter, n_pos,
                                  utr_type, tss, args_srna)
                if tss.start > end_fuzzy:
                    break
    if (utr_type == "3utr") and (len(args_srna.pros) != 0):
        pos = import_position(start, end, start, end)
        detect_3utr_pro(inter, pos, wigs, utr_type, args_srna)
    if (utr_type == "interCDS") and (len(args_srna.pros) != 0):
        pos = import_position(start, end, start, end)
        detect_twopro(inter, pos, wigs, utr_type, "interCDS", args_srna)
    args_srna.fuzzy_tsss[utr_type] = ori_fuzzy


def run_utr_detection(wigs, inter, start, end, utr_type, args_srna):
    if end - start >= 0:
        if (utr_type == "3utr") or (
                utr_type == "5utr") or (
                utr_type == "interCDS"):
            detect_utr(start, end, inter, utr_type, wigs, args_srna)
        else:
            pos = import_position(start, end, start, end)
            srna_covers, utr_covers = get_coverage(
                    wigs, inter, pos, utr_type, "NA", args_srna)
            args_srna.utrs.append(import_data(
                inter["strand"], inter["strain"],
                pos, utr_type, "NA", "NA", utr_covers, "NA"))


def class_utr(inter, ta, args_srna, wig_fs, wig_rs):
    '''classify the UTR-dervied sRNA'''
    if inter["strand"] == "+":
        if (inter["start"] <= ta.end) and (
                inter["end"] >= ta.end) and (
                ta.start <= inter["start"]):
            run_utr_detection(wig_fs, inter, inter["start"] + 1,
                              ta.end, "3utr", args_srna)
        elif (inter["start"] <= ta.start) and (
                inter["end"] >= ta.start) and (
                ta.end >= inter["end"]):
            run_utr_detection(wig_fs, inter, ta.start,
                              inter["end"] - 1, "5utr", args_srna)
        elif (inter["start"] >= ta.start) and (inter["end"] <= ta.end):
            run_utr_detection(wig_fs, inter, inter["start"] + 1,
                              inter["end"] - 1, "interCDS", args_srna)
    else:
        if (inter["start"] <= ta.end) and (
                inter["end"] >= ta.end) and (
                ta.start <= inter["start"]):
            run_utr_detection(wig_rs, inter, inter["start"] + 1,
                              ta.end, "5utr", args_srna)
        elif (inter["start"] <= ta.start) and (
                inter["end"] >= ta.start) and (
                ta.end >= inter["end"]):
            run_utr_detection(wig_rs, inter, ta.start,
                              inter["end"] - 1, "3utr", args_srna)
        elif (inter["start"] >= ta.start) and (inter["end"] <= ta.end):
            run_utr_detection(wig_rs, inter,  inter["start"] + 1,
                              inter["end"] - 1, "interCDS", args_srna)


def median_score(lst, per):
    '''if the cutoff is assigned by precentage, 
    it can get the corresponding number'''
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = int((lstLen - 1) * per)
    if lstLen != 0:
        return sortedLst[index]
    else:
        return 0


def print_file(num, srna, start, end, srna_datas, args_srna):
    name = '%0*d' % (5, num)
    args_srna.out_t.write(
            "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t".format(
                srna["strain"], name, start, end, srna["strand"],
                ";".join(srna_datas["conds"].keys()),
                ";".join(srna_datas["conds"].values()),
                srna_datas["best"], srna_datas["high"], srna_datas["low"]))
    attribute_string = ";".join(
        ["=".join(items) for items in [
            ["ID", "srna_utr" + str(num)], ["Name", "UTR_sRNA_" + name],
            ["sRNA_type", srna["utr"]],
            ["best_avg_coverage", str(srna_datas["best"])],
            ["best_high_coverage", str(srna_datas["high"])],
            ["best_low_coverage", str(srna_datas["low"])],
            ["with_TSS", srna["start_tss"]],
            ["start_cleavage", srna["start_cleavage"]],
            ["end_cleavage", srna["end_cleavage"]]]])
    args_srna.out.write("\t".join([str(field) for field in [
                        srna["strain"], "ANNOgesic", "ncRNA", str(start),
                        str(end), ".", srna["strand"], ".",
                        attribute_string]]) + "\n")
    if not args_srna.table_best:
        first = True
        for data in srna_datas["detail"]:
            if first:
                args_srna.out_t.write("{0}(avg={1};high={2};low={3})".format(
                                      data["track"], data["avg"],
                                      data["high"], data["low"]))
                first = False
            else:
                args_srna.out_t.write(";{0}(avg={1};high={2};low={3})".format(
                                      data["track"], data["avg"],
                                      data["high"], data["low"]))
    else:
        args_srna.out_t.write("{0}(avg={1};high={2};low={3})".format(
                              srna_datas["track"], srna_datas["best"],
                              srna_datas["high"], srna_datas["low"]))
    args_srna.out_t.write("\n")


def detect_srna(median, args_srna):
    '''check the sRNA candidates and print it out'''
    num = 0
    if len(args_srna.srnas) != 0:
        for srna in args_srna.srnas:
            if srna["strain"] in median.keys():
                srna_datas = replicate_comparison(
                        args_srna, srna["datas"],
                        srna["strand"], "sRNA_utr_derived",
                        median[srna["strain"]][srna["utr"]],
                        args_srna.coverages, srna["utr"], None, None,
                        args_srna.texs)
                if srna_datas["best"] != 0:
                    if (srna["utr"] == "5utr") or (
                            srna["utr"] == "interCDS"):
                        start = srna_datas["start"]
                        end = srna_datas["end"]
                    elif srna["utr"] == "3utr":
                        start = srna["start"]
                        end = srna["end"]
                    if (math.fabs(start - end) >= args_srna.min_len) and (
                            math.fabs(start - end) <= args_srna.max_len):
                        print_file(num, srna, start, end,
                                   srna_datas, args_srna)
                        num += 1


def read_data(args_srna):
    cdss = []
    tas = []
    tsss = []
    pros = []
    seq = {}
    gff_parser = Gff3Parser()
    fh = open(args_srna.gff_file, "r")
    for entry in gff_parser.entries(fh):
        if (entry.feature == "CDS") or (
                entry.feature == "tRNA") or (
                entry.feature == "rRNA"):
            if ("product" in entry.attributes.keys()) and (args_srna.hypo):
                if "hypothetical protein" not in entry.attributes["product"]:
                    cdss.append(entry)
            else:
                cdss.append(entry)
    fh.close()
    fh = open(args_srna.ta_file, "r")
    for entry in gff_parser.entries(fh):
        tas.append(entry)
    fh.close()
    fh = open(args_srna.tss_file, "r")
    for entry in gff_parser.entries(fh):
        tsss.append(entry)
    fh.close()
    if args_srna.pro_file is not None:
        fh = open(args_srna.pro_file, "r")
        for entry in gff_parser.entries(fh):
            pros.append(entry)
        fh.close()
    with open(args_srna.seq_file, "r") as s_f:
        for line in s_f:
            line = line.strip()
            if line.startswith(">"):
                strain = line[1:]
                seq[strain] = ""
            else:
                seq[strain] = seq[strain] + line
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    tas = sorted(tas, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    if len(pros) != 0:
        pros = sorted(pros, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    return cdss, tas, tsss, pros, seq


def get_utr_coverage(utrs):
    covers = {}
    for utr in utrs:
        if utr["strain"] not in covers.keys():
            covers[utr["strain"]] = {"3utr": {}, "5utr": {},
                                     "interCDS": {}}
        for cond, utr_covers in utr["datas"].items():
            for cover in utr_covers:
                if (cover["track"] not in
                        covers[utr["strain"]][utr["utr"]].keys()):
                    covers[utr["strain"]][utr["utr"]][cover["track"]] = []
                covers[utr["strain"]][utr["utr"]][cover["track"]].append(
                        cover["ori_avg"])
    return covers


def get_inter(cdss, inters):
    '''get the intergenic region'''
    for cds1 in cdss:
        for cds2 in cdss:
            if (cds1.seq_id == cds2.seq_id) and \
               (cds1.strand == cds2.strand):
                if (cds2.start > cds1.start):
                    if cds2.start - cds1.end > 1:
                        if cds1.strand == "+":
                            length = cds1.end - cds1.start
                        else:
                            length = cds2.end - cds2.start
                        inters.append(import_inter(cds1.strand, cds1.seq_id,
                                      cds1.end, cds2.start, length))
                    break


def mean_score(lst):
    total = 0
    for li in lst:
        total = total + li
    if len(lst) != 0:
        return (total / len(lst))
    else:
        return 0


def get_utr_cutoff(coverage, mediandict, avgs, strain, utr, track):
    if "n_" in coverage:
        cutoff = float(coverage.split("_")[-1])
        mediandict[strain][utr][track] = {"median": cutoff,
                                          "mean": mean_score(avgs)}
    elif "p_" in coverage:
        cutoff = float(coverage.split("_")[-1])
        mediandict[strain][utr][track] = {"median": median_score(avgs, cutoff),
                                          "mean": mean_score(avgs)}


def set_cutoff(covers, args_srna):
    '''set the cutoff based on the types of sRNA'''
    mediandict = {}
    for strain, utrs in covers.items():
        mediandict[strain] = {"3utr": {}, "5utr": {}, "interCDS": {}}
        for utr, tracks in utrs.items():
            if (utr == "3utr") or (utr == "5utr") or (utr == "interCDS"):
                for track, avgs in tracks.items():
                    if track not in mediandict[strain][utr].keys():
                        mediandict[strain][utr][track] = {}
                    if args_srna.cover_notex is not None:
                        for keys in args_srna.texs.keys():
                            tracks = keys.split("@AND@")
                            if tracks[0] == track:
                                get_utr_cutoff(args_srna.coverages[utr],
                                               mediandict,
                                               avgs, strain, utr, track)
                                break
                            elif tracks[1] == track:
                                get_utr_cutoff(args_srna.cover_notex[utr],
                                               mediandict,
                                               avgs, strain, utr, track)
                                break
                    else:
                        get_utr_cutoff(args_srna.coverages[utr], mediandict,
                                       avgs, strain, utr, track)
    return mediandict


def print_median(out_folder, mediandict):
    '''print the cutoff based on the types of sRNA'''
    out = open(os.path.join(out_folder, "tmp_median"), "a")
    for strain, utrs in mediandict.items():
        for utr, tracks in utrs.items():
            for track, value in tracks.items():
                out.write("\t".join([strain, utr, track,
                                     str(value["median"])]) + "\n")
    out.close()


def free_memory(paras):
    for data in (paras):
        del(data)
    gc.collect()


def utr_derived_srna(args_srna, libs, texs, wig_fs, wig_rs):
    inters = []
    cdss, tas, tsss, pros, seq = read_data(args_srna)
    out = open(args_srna.output_file, "w")
    out.write("##gff-version 3\n")
    out_t = open(args_srna.output_table, "w")
    get_terminal(cdss, inters, seq, "start")
    get_inter(cdss, inters)
    get_terminal(cdss, inters, seq, "end")
    inters = sorted(inters, key=lambda k: (k["strain"], k["start"],
                                           k["end"], k["strand"]))
    args_srna = ArgsContainer().extend_utr_container(
                            args_srna, cdss, tsss, pros, out,
                            out_t, texs)
    for inter in inters:
        for ta in tas:
            if (inter["strain"] == ta.seq_id) and (
                    inter["strand"] == ta.strand):
                class_utr(inter, ta, args_srna, wig_fs, wig_rs)
    covers = get_utr_coverage(args_srna.utrs)
    mediandict = set_cutoff(covers, args_srna)
    print_median(args_srna.out_folder, mediandict)
    detect_srna(mediandict, args_srna)
    args_srna.out.close()
    args_srna.out_t.close()
    paras = [args_srna.srnas, args_srna.utrs, seq, inters,
             tas, cdss, tas, tsss, pros, covers]
    free_memory(paras)
