import os
import sys
import csv
import math
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.lib_reader import read_wig, read_libs
from annogesiclib.coverage_detection import coverage_comparison, replicate_comparison

def import_data(strand, strain, start, end, utr, type_, name, srna_cover, pro):
    if type_ == "TSS":
        return {"strand": strand, "strain": strain, "start": start,
                "end": end, "utr": utr, "start_tss": name, "end_cleavage": pro,
                "start_cleavage": "NA", "datas": srna_cover}
    elif type_ == "cleavage":
        return {"strand": strand, "strain": strain, "start": start, "end": end,
                "utr": utr, "start_tss": "NA", "start_cleavage": name,
                "datas": srna_cover, "end_cleavage": pro}
    else:
        return {"strand": strand, "strain": strain, "start": start, "end": end,
                "utr": utr, "start_tss": "NA", "start_cleavage": "NA",
                "datas": srna_cover, "end_cleavage": pro}

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
                inters.append(import_inter(cds.strand, cds.seq_id, 1, cds.start, 0))
            elif (cds.strand == "-") and (first_m):
                first_m = False
                inters.append(import_inter(cds.strand, cds.seq_id, 1, cds.start, 0))
    elif type_ == "end":
        for cds in reversed(cdss):
            if cds.seq_id != pre_strain:
                pre_strain = cds.seq_id
                first_p = True
                first_m = True
            if (cds.strand == "+") and (first_p):
                first_p = False
                inters.append(import_inter(cds.strand, cds.seq_id, cds.end,
                              len(seq[cds.seq_id]), cds.end - cds.start))
            elif (cds.strand == "-") and (first_m):
                first_m = False
                inters.append(import_inter(cds.strand, cds.seq_id, cds.start,
                              len(seq[cds.seq_id]), cds.end - cds.start))

def check_pos(cover, check_point, checks):
    if (cover["pos"] >= min(check_point["utr_start"], check_point["utr_end"])) and (
        cover["pos"] <= max(check_point["utr_start"], check_point["utr_end"])):
        checks["utr"] = True
    if (cover["pos"] >= min(check_point["srna_start"], check_point["srna_end"])) and (
        cover["pos"] <= max(check_point["srna_start"], check_point["srna_end"])):
        checks["srna"] = True

def set_cover_and_point(inter, covers, ori_start, ori_end, fuzzy_end,
                        start, end):
    check_point = {"srna_start": 0, "srna_end": 0,
                   "utr_start": 0, "utr_end": 0}
    if inter["strand"] == "-":
        if ori_start - 2 < 0:
            ori_start = 2
        covers = reversed(covers[ori_start - 2 - fuzzy_end: ori_end + 1])
        check_point["srna_start"] = end + 1
        check_point["srna_end"] = start - 2 - fuzzy_end
        check_point["utr_start"] = ori_end
        check_point["utr_end"] = ori_start
    elif inter["strand"] == "+":
        if ori_start - 2 < 0:
            ori_start = 2
        covers = covers[ori_start - 2: ori_end + 1 + fuzzy_end]
        check_point["srna_start"] = start - 2
        check_point["srna_end"] = end + 1 + fuzzy_end
        check_point["utr_start"] = ori_start
        check_point["utr_end"] = ori_end
    return covers, check_point

def detect_cover_utr_srna(covers, check_point, cover_sets, start, end,
                          poss, inter, type_, intercds_type, decrease,
                          fuzzy_end, srna_covers, utr_covers, cond,
                          track, ori_start, ori_end, max_len, min_len):
    num = 0
    cover_tmp = {"5utr": 0, "total": 0, "ori_total": 0}
    checks = {"first": True, "detect_decrease": False,
              "srna": False, "utr": False}
    final_poss = {"start": start, "end": end}
    for cover in covers:
        checks["srna"] = False
        checks["utr"] = False
        check_pos(cover, check_point, checks)
        if checks["utr"]:
            cover_tmp["ori_total"] = \
                cover_tmp["ori_total"] + cover["coverage"]
        if checks["srna"]:
            cover_tmp["total"] = \
                cover_tmp["total"] + cover["coverage"]
            checks["first"] = coverage_comparison(
                              cover, cover_sets, poss,
                              checks["first"], inter["strand"])
            if (checks["first"] is not True) and (
                cover_sets["high"] > 0):
                if (type_ == "5utr") or (type_ == "3utr") or (
                    (type_ == "interCDS") and (
                     intercds_type == "TSS")):
                    if ((cover_sets["low"] / cover_sets["high"]) < decrease) and (
                        cover_sets["low"] > -1):
                        checks["detect_decrease"] = True
                        cover_tmp["5utr"] = cover["coverage"]
            if checks["detect_decrease"]:
                num, go_out = get_cover_5utr(num, fuzzy_end,
                              cover_sets, cover_tmp, final_poss,
                              decrease, cover, inter)
                if go_out is True:
                    break
    if (checks["first"] is not True) and (cover_sets["high"] > 0):
        check_import_srna_covers(checks, final_poss, end,
            start, type_, cover_sets, cover_tmp, intercds_type,
            inter, srna_covers, cond, track, cover, ori_start,
            ori_end, utr_covers, max_len, min_len)

def get_coverage(wigs, inter, start, end, type_, intercds_type, fuzzy_end,
                 decrease, ori_start, ori_end, max_len, min_len):
    srna_finals = ""
    srna_covers = {}
    utr_covers = {}
    srna_tracks = []
    cover_sets = {"high": 0, "low": 0, "best": -1}
    poss = {"low": 0, "high": 0}
    conditions = {}
    for wig_strain, conds in wigs.items():
        if wig_strain == inter["strain"]:
            for cond, tracks in conds.items():
                srna_covers[cond] = []
                utr_covers[cond] = []
                for track, covers in tracks.items():
                    covers, check_point = set_cover_and_point(inter, covers,
                                          ori_start, ori_end, fuzzy_end, start,
                                          end)
                    detect_cover_utr_srna(covers, check_point, cover_sets,
                          start, end, poss, inter, type_, intercds_type,
                          decrease, fuzzy_end, srna_covers, utr_covers, cond,
                          track, ori_start, ori_end, max_len, min_len)
    return srna_covers, utr_covers

def get_cover_5utr(num, fuzzy_end, cover_sets, cover_tmp,
                   final_poss, decrease, cover, inter):
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
    return num, go_out

def import_cover(inter, covers, track, cover_sets, avg, cover, ori_avg,
                 final_poss, start, end):
    if final_poss["start"] < final_poss["end"]:
        if (inter["strand"] == "+") and (final_poss["end"] > end):
            final_poss["end"] = end
        elif (inter["strand"] == "-") and (final_poss["start"] < start):
            final_poss["start"] = start
        covers.append({"track": track,
                       "high": cover_sets["high"],
                       "low": cover_sets["low"], "avg": avg,
                       "type": cover["type"],
                       "ori_avg": ori_avg,
                       "final_start": final_poss["start"],
                       "final_end": final_poss["end"]})

def check_import_srna_covers(checks, final_poss, end, start, type_, cover_sets,
                             cover_tmp, intercds_type, inter, srna_covers, cond,
                             track, cover, ori_start, ori_end, utr_covers,
                             max_len, min_len):
    avg = cover_tmp["total"] / float(end - start + 1)
    ori_avg = cover_tmp["ori_total"] / float(ori_end - ori_start + 1)
    if ((type_ == "5utr") and ((intercds_type == "tsspro") or (
                                checks["detect_decrease"]))) or (
        (type_ == "3utr") and (intercds_type == "two_pro")) or (
        (type_ == "3utr") and (intercds_type != "two_pro") and (
                               checks["detect_decrease"])) or (
        (type_ == "3utr") and (intercds_type != "two_pro") and (
                              (end - start) >= min_len) and (
                              (end - start) <= max_len)) or (
        (type_ == "interCDS") and (intercds_type == "TSS") and (
                                   checks["detect_decrease"])) or (
        (type_ == "interCDS") and ((intercds_type == "tss_pro") or (
                                    intercds_type == "two_pro"))):
        import_cover(inter, srna_covers[cond], track, cover_sets,
                     avg, cover, ori_avg, final_poss, start, end)
        import_cover(inter, utr_covers[cond], track, cover_sets,
                     avg, cover, ori_avg, final_poss, start, end)
    else:
        import_cover(inter, utr_covers[cond], track, cover_sets,
                     avg, cover, ori_avg, final_poss, start, end)

def detect_3utr_pro(pros, inter, utrs, srnas, start, end, wigs, feature, name,
                    utr_type, fuzzys, fuzzy_end, decrease, min_len, max_len):
    for pro in pros:
        if (pro.seq_id == inter["strain"]) and (
            pro.strand == inter["strand"]):
            if (pro.start >= start) and (pro.start <= end):
                if pro.strand == "+":
                    if ((end - pro.start) >= min_len) and (
                        (end - pro.start) <= max_len):
                        srna_covers, utr_covers = get_coverage(wigs, inter,
                                 pro.start, end, utr_type, "pro", fuzzy_end,
                                 decrease,start, end, max_len, min_len)
                        utrs.append(import_data(inter["strand"],
                             inter["strain"], pro.start, end, utr_type,
                             feature, name, utr_covers, "NA"))
                        srnas.append(import_data(inter["strand"],
                              inter["strain"], pro.start, end, "3utr",
                              "cleavage", "Cleavage:" + "_".join(
                              [str(pro.start), pro.strand]), srna_covers, "NA"))
                    elif (end - pro.start) > max_len:
                        detect_twopro(pros, inter, utrs, srnas, start,
                                      end, wigs, feature, name, utr_type,
                                      fuzzy_end, decrease, min_len,
                                      max_len, "3utr")
                else:
                    if ((pro.start - start) >= min_len) and (
                        (pro.start - start) <= max_len):
                        srna_covers, utr_covers = get_coverage(wigs, inter,
                                 start, pro.start, utr_type, "pro", fuzzy_end,
                                 decrease, start, end, max_len, min_len)
                        utrs.append(import_data(inter["strand"],
                                    inter["strain"], start, pro.start,
                                    utr_type, feature, name, utr_covers, "NA"))
                        srnas.append(import_data(inter["strand"],
                              inter["strain"], start, pro.start,
                              "3utr", "cleavage", "Cleavage:" + \
                              "_".join([str(pro.start), pro.strand]),
                              srna_covers, "NA"))
                    elif (pro.start - start) > max_len:
                        detect_twopro(pros, inter, utrs, srnas, start,
                                      end, wigs, feature, name, utr_type,
                                      fuzzy_end, decrease, min_len,
                                      max_len, "3utr")
            if (pro.start > end + fuzzys[utr_type]):
                break

def detect_twopro(pros, inter, utrs, srnas, start, end, wigs, feature, name,
                  utr_type, fuzzy_end, decrease, min_len, max_len, import_type):
    poss = []
    for pro in pros:
        if (pro.seq_id == inter["strain"]) and (
            pro.strand == inter["strand"]):
            if (pro.start >= start) and (pro.start <= end):
                poss.append(pro)
    first = True
    for pos in poss:
        if first:
            first = False
        else:
            if ((pos.start - pre_pos.start) >= min_len) and (
                (pos.start - pre_pos.start) <= max_len):
                srna_covers, utr_covers = get_coverage(wigs, inter,
                                          pre_pos.start, pos.start,
                                          utr_type, "two_pro", fuzzy_end,
                                          decrease, start, end, max_len, min_len)
                utrs.append(import_data(inter["strand"], inter["strain"],
                            pre_pos.start, pos.start, utr_type, feature,
                            name, utr_covers,
                            "Cleavage:" + "_".join([str(pos.start), pos.strand])))
                srnas.append(import_data(inter["strand"], inter["strain"],
                             pre_pos.start, pos.start, import_type, "cleavage",
                             "Cleavage:" + "_".join([str(pre_pos.start), pos.strand]),
                             srna_covers,
                             "Cleavage:" + "_".join([str(pos.start), pos.strand])))
        pre_pos = pos

def decrease_pos(covers, pos, strand):
    longer = -1
    for cond, datas in covers.items():
        for data in datas:
            for key, value in data.items():
                if key == pos:
                    if longer == -1:
                        longer = value
                    elif (strand == "+") and (value > longer) and (longer != -1):
                        longer = value
                    elif (strand == "-") and (value < longer) and (longer != -1):
                        longer = value
    return longer

def get_decrease(inter, wigs, tss, start, end, utr_type, fuzzy_end, decrease,
                 ori_start, ori_end, utrs, srnas, feature, name, max_len, min_len):
    if inter["strand"] == "+":
        srna_covers, utr_covers = get_coverage(wigs, inter,
                 tss.start, end, utr_type, "TSS",
                 fuzzy_end, decrease, ori_start, ori_end, max_len, min_len)
        utr_pos = decrease_pos(utr_covers, "final_end", "+")
        utrs.append(import_data(inter["strand"],
                    inter["strain"], tss.start,
                    utr_pos, utr_type, feature,
                    name, utr_covers, "NA"))
        if len(srna_covers) != 0:
            srna_pos = decrease_pos(srna_covers, "final_end", "+")
            srnas.append(import_data(inter["strand"],
                         inter["strain"], tss.start,
                         srna_pos, utr_type, "TSS",
                         "TSS:" + "_".join([str(tss.start),
                         tss.strand]), srna_covers, "NA"))
    else:
        srna_covers, utr_covers = get_coverage(wigs, inter,
                 start, tss.start, utr_type, "TSS",
                 fuzzy_end, decrease, ori_start, ori_end, max_len, min_len)
        utr_pos = decrease_pos(utr_covers, "final_start", "-")
        utrs.append(import_data(inter["strand"],
                    inter["strain"], utr_pos,
                    tss.start, utr_type, feature,
                    name, utr_covers, "NA"))
        if len(srna_covers) != 0:
            srna_pos = decrease_pos(srna_covers, "final_start", "-")
            srnas.append(import_data(inter["strand"],
                         inter["strain"], srna_pos,
                         tss.start, utr_type, "TSS",
                         "TSS:" + "_".join([str(tss.start),
                         tss.strand]), srna_covers, "NA"))

def detect_normal(diff, wigs, inter, start, end, utr_type, feature, name, tss,
                  pros, utrs, srnas, fuzzy_end, decrease, min_len, max_len,
                  ori_start, ori_end):
    if (diff >= min_len) and (diff <= max_len):
        srna_covers, utr_covers = get_coverage(wigs, inter, start, end,
                                  utr_type, "TSS", fuzzy_end, decrease,
                                  ori_start, ori_end, max_len, min_len)
        utrs.append(import_data(inter["strand"], inter["strain"], start, end,
                    utr_type, feature, name, utr_covers, "NA"))
        srnas.append(import_data(inter["strand"], inter["strain"], start,
                     end, utr_type, "TSS",
                     "TSS:" + "_".join([str(tss.start), tss.strand]), srna_covers, "NA"))
    elif (diff > max_len) and (len(pros) != 0):
        detect = False
        for pro in pros:
            if (pro.seq_id == inter["strain"]) and (
                pro.strand == inter["strand"]):
                if (pro.start >= start) and (pro.start <= end):
                    if ((pro.start - start) >= min_len) and (
                        (pro.start - start) <= max_len):
                        detect = True
                        if inter["strand"] == "+":
                            srna_covers, utr_covers = get_coverage(wigs, inter,
                                     tss.start, pro.start, utr_type, "tsspro",
                                     fuzzy_end, decrease, ori_start, ori_end,
                                     max_len, min_len)
                            utrs.append(import_data(inter["strand"],
                                        inter["strain"], tss.start,
                                        pro.start, utr_type, feature,
                                        name, utr_covers,
                                        "Cleavage:" + "_".join([str(pro.start), pro.strand])))
                            srnas.append(import_data(inter["strand"],
                                         inter["strain"], tss.start,
                                         pro.start, utr_type, "TSS",
                                         "TSS:" + "_".join([str(tss.start),
                                         tss.strand]), srna_covers,
                                         "Cleavage:" + "_".join([str(pro.start), pro.strand])))
                        else:
                            srna_covers, utr_covers = get_coverage(wigs, inter,
                                     pro.start, tss.start, utr_type, "tsspro",
                                     fuzzy_end, decrease, ori_start, ori_end,
                                     max_len, min_len)
                            utrs.append(import_data(inter["strand"],
                                        inter["strain"], pro.start,
                                        tss.start, utr_type, feature,
                                        name, utr_covers,
                                        "Cleavage:" + "_".join([str(pro.start), pro.strand])))
                            srnas.append(import_data(inter["strand"],
                                         inter["strain"], pro.start,
                                         tss.start, utr_type, "TSS",
                                         "TSS:" + "_".join([str(tss.start),
                                         tss.strand]), srna_covers,
                                         "Cleavage:" + "_".join([str(pro.start), pro.strand])))
        if not detect:
            get_decrease(inter, wigs, tss, start, end, utr_type, fuzzy_end, decrease,
                         ori_start, ori_end, utrs, srnas, feature, name, max_len, min_len)

def detect_utr(srnas, tsss, pros, start, end, inter, utr_type, utrs, wigs,
               feature, name, fuzzys, fuzzy_end, decrease, min_len, max_len):
    ori_fuzzy = fuzzys[utr_type]
    if "p_" in fuzzys[utr_type]:
        per = float(fuzzys[utr_type].split("_")[-1])
        fuzzys[utr_type] = inter["len_CDS"]*per
    elif "n_" in fuzzys[utr_type]:
        fuzzys[utr_type] = float(fuzzys[utr_type].split("_")[-1])
    for tss in tsss:
        if (tss.seq_id == inter["strain"]) and (
            tss.strand == inter["strand"]):
            if tss.strand == "+":
                start_fuzzy = start - fuzzys[utr_type]
                if (tss.start >= start_fuzzy) and (tss.start <= end):
                    detect_normal((end - tss.start), wigs, inter, tss.start,
                                  end, utr_type, feature, name, tss, pros,
                                  utrs, srnas, fuzzy_end, decrease, min_len,
                                  max_len, start, end)
                elif tss.start > end:
                    break
            else:
                end_fuzzy = end + fuzzys[utr_type]
                if (tss.start >= start) and (tss.start <= end_fuzzy):
                    detect_normal((tss.start - start), wigs, inter, start,
                                  tss.start, utr_type, feature, name, tss,
                                  pros, utrs, srnas, fuzzy_end, decrease,
                                  min_len, max_len, start, end)
                if tss.start > end_fuzzy:
                    break
    if (utr_type == "3utr") and (len(pros) != 0):
        detect_3utr_pro(pros, inter, utrs, srnas, start, end, wigs, feature,
                name, utr_type, fuzzys, fuzzy_end, decrease, min_len, max_len)
    if (utr_type == "interCDS") and (len(pros) != 0):
        detect_twopro(pros, inter, utrs, srnas, start, end, wigs, feature, name,
                      utr_type, fuzzy_end, decrease, min_len, max_len, "interCDS")
    fuzzys[utr_type] = ori_fuzzy

def run_utr_detection(wigs, inter, start, end, utr_type, utrs, tsss,
                      pros, srnas, feature, name, fuzzys, fuzzy_end,
                      decrease, min_len, max_len):
    if end - start >= 0:
        if (utr_type == "3utr") or (
            utr_type == "5utr") or (
            utr_type == "interCDS"):
            detect_utr(srnas, tsss, pros, start, end, inter, utr_type, utrs,
                       wigs, feature, name, fuzzys, fuzzy_end, decrease,
                       min_len, max_len)
        else:
            srna_covers, utr_covers = get_coverage(wigs, inter, start, end,
                                      utr_type, "NA", fuzzy_end, decrease,
                                      start, end, max_len, min_len)
            utrs.append(import_data(inter["strand"], inter["strain"], start,
                        end, utr_type, feature, name, utr_covers, "NA"))

def class_utr(inter, ta, utrs, srnas, tsss, pros, wig_fs, wig_rs,
              fuzzys, fuzzy_end, decrease, min_len, max_len):
    if inter["strand"] == "+":
        if (inter["start"] <= ta.end) and (
            inter["end"] >= ta.end) and (
            ta.start <= inter["start"]):
            run_utr_detection(wig_fs, inter, inter["start"] + 1, ta.end, "3utr",
                              utrs, tsss, pros, srnas, "NA", "NA", fuzzys,
                              fuzzy_end, decrease, min_len, max_len)
        elif (inter["start"] <= ta.start) and (
              inter["end"] >= ta.start) and (
              ta.end >= inter["end"]):
            run_utr_detection(wig_fs, inter, ta.start, inter["end"] - 1, "5utr",
                              utrs, tsss, pros, srnas, "NA", "NA", fuzzys,
                              fuzzy_end, decrease, min_len, max_len)
        elif (inter["start"] >= ta.start) and (inter["end"] <= ta.end):
            run_utr_detection(wig_fs, inter, inter["start"] + 1,
                              inter["end"] - 1, "interCDS", utrs,
                              tsss, pros, srnas, "NA", "NA", fuzzys, fuzzy_end,
                              decrease, min_len, max_len)
    else:
        if (inter["start"] <= ta.end) and (
            inter["end"] >= ta.end) and (
            ta.start <= inter["start"]):
            run_utr_detection(wig_rs, inter, inter["start"] + 1, ta.end, "5utr",
                              utrs, tsss, pros, srnas, "NA", "NA", fuzzys,
                              fuzzy_end, decrease, min_len, max_len)
        elif (inter["start"] <= ta.start) and (
              inter["end"] >= ta.start) and (
              ta.end >= inter["end"]):
            run_utr_detection(wig_rs, inter, ta.start, inter["end"] - 1, "3utr",
                              utrs, tsss, pros, srnas, "NA", "NA", fuzzys,
                              fuzzy_end, decrease, min_len, max_len)
        elif (inter["start"] >= ta.start) and (inter["end"] <= ta.end):
            run_utr_detection(wig_rs, inter, inter["start"] + 1,
                              inter["end"] - 1, "interCDS", utrs,
                              tsss, pros, srnas, "NA", "NA", fuzzys, fuzzy_end,
                              decrease, min_len, max_len)

def median_score(lst, per):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = int((lstLen - 1) * per)
    if lstLen != 0:
        return sortedLst[index]
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
        ["Name", "UTR_sRNA_" + name], ["sRNA_type", srna["utr"]],
        ["best_avg_coverage", str(srna_datas["best"])],
        ["best_high_coverage", str(srna_datas["high"])],
        ["best_low_coverage", str(srna_datas["low"])],
        ["with_TSS", srna["start_tss"]],
        ["start_cleavage", srna["start_cleavage"]],
        ["end_cleavage", srna["end_cleavage"]]]])
    out.write("\t".join([str(field) for field in [
              srna["strain"], "ANNOgesic", "sRNA", str(start),
              str(end), ".", srna["strand"], ".", attribute_string]]) + "\n")
    if not table_best:
        first = True
        for data in srna_datas["detail"]:
            if first:
                out_t.write("{0}(avg={1};high={2};low={3})".format(
                            data["track"], data["avg"],
                            data["high"], data["low"]))
                first = False
            else:
                out_t.write(";{0}(avg={1};high={2};low={3})".format(
                            data["track"], data["avg"],
                            data["high"], data["low"]))
    else:
        out_t.write("{0}(avg={1};high={2};low={3})".format(
                    srna_datas["track"], srna_datas["best"],
                    srna_datas["high"], srna_datas["low"]))
    out_t.write("\n")

def detect_srna(srnas, median, out, out_t, template_texs, coverages,
                tex_notex, replicates, table_best, max_len, min_len):
    num = 0
    if len(srnas) != 0:
        for srna in srnas:
            if srna["strain"] in median.keys():
                srna_datas = replicate_comparison(srna["datas"],
                             template_texs, srna["strand"], None,
                             tex_notex, replicates, "sRNA_utr_derived",
                             median[srna["strain"]][srna["utr"]],
                             coverages, srna["utr"], template_texs, None)
                if srna_datas["best"] != 0:
                    if (srna["utr"] == "5utr") or (
                        srna["utr"] == "interCDS"):
                        start = srna_datas["start"]
                        end = srna_datas["end"]
                    elif srna["utr"] == "3utr":
                        start = srna["start"]
                        end = srna["end"]
                    if (math.fabs(start - end) >= min_len) and (
                        math.fabs(start - end) <= max_len):
                        print_file(num, out_t, out, srna, start, end,
                                   srna_datas, table_best)
                        num += 1

def read_data(gff_file, ta_file, tss_file, pro_file, seq_file, hypo):
    cdss = []
    tas = []
    tsss = []
    pros = []
    seq = {}
    gff_parser = Gff3Parser()
    fh = open(gff_file, "r")
    for entry in gff_parser.entries(fh):
        if (entry.feature == "CDS") or (
            entry.feature == "tRNA") or (
            entry.feature == "rRNA"):
            if ("product" in entry.attributes.keys()) and (hypo):
                if "hypothetical protein" not in entry.attributes["product"]:
                    cdss.append(entry)
            else:
                cdss.append(entry)
    fh.close()
    fh = open(ta_file, "r")
    for entry in gff_parser.entries(fh):
        tas.append(entry)
    fh.close()
    fh = open(tss_file, "r")
    for entry in gff_parser.entries(fh):
        tsss.append(entry)
    fh.close()
    if pro_file is not None:
        fh = open(pro_file, "r")
        for entry in gff_parser.entries(fh):
            pros.append(entry)
        fh.close()
    with open(seq_file, "r") as s_f:
        for line in s_f:
            line = line.strip()
            if line.startswith(">"):
                strain = line[1:]
                seq[strain] = ""
            else:
                seq[strain] = seq[strain] + line
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start))
    tas = sorted(tas, key=lambda k: (k.seq_id, k.start))
    tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start))
    if len(pros) != 0:
        pros = sorted(pros, key=lambda k: (k.seq_id, k.start))
    fh.close()
    return cdss, tas, tsss, pros, seq

def get_utr_coverage(utrs):
    covers = {}
    first = True
    for utr in utrs:
        best_cover = -1
        avgs = []
        if utr["strain"] not in covers.keys():
            covers[utr["strain"]] = {"3utr": {}, "5utr": {},
                                     "interCDS": {}}
        for cond, utr_covers in utr["datas"].items():
            for cover in utr_covers:
                if cover["track"] not in covers[utr["strain"]][utr["utr"]].keys():
                    covers[utr["strain"]][utr["utr"]][cover["track"]] = []
                covers[utr["strain"]][utr["utr"]][cover["track"]].append(cover["ori_avg"])
    return covers

def get_inter(cdss, inters):
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

def set_cutoff(covers, coverages, cover_notex, texs):
    mediandict = {}
    for strain, utrs in covers.items():
        mediandict[strain] = {"3utr": {}, "5utr": {}, "interCDS": {}}
        for utr, tracks in utrs.items():
            if (utr == "3utr") or (utr == "5utr") or (utr == "interCDS"):
                for track, avgs in tracks.items():
                    if track not in mediandict[strain][utr].keys():
                        mediandict[strain][utr][track] = {}
                    if cover_notex is not None:
                        for keys in texs.keys():
                            tracks = keys.split("@AND@")
                            if tracks[0] == track:
                                get_utr_cutoff(coverages[utr], mediandict,
                                               avgs, strain, utr, track)
                                break
                            elif tracks[1] == track:
                                get_utr_cutoff(cover_notex[utr], mediandict,
                                               avgs, strain, utr, track)
                                break
                    else:
                        get_utr_cutoff(coverages[utr], mediandict,
                                       avgs, strain, utr, track)
    return mediandict

def print_median(out_folder, mediandict):
    out = open(os.path.join(out_folder, "tmp_median"), "a")
    for strain, utrs in mediandict.items():
        for utr, tracks in utrs.items():
            for track, value in tracks.items():
                out.write("\t".join([strain, utr, track, str(value["median"])]) + "\n")
    out.close()

def utr_derived_srna(gff_file, ta_file, tss_file, wig_f_file, wig_r_file,
        pro_file, max_len, min_len, fuzzys, fuzzy_end, seq_file, wig_folder,
        input_libs, tex_notex, replicates, output_file, output_table,
        table_best, decrease, utr_coverages, hypo, out_folder, notex):
    inters = []
    wigs = []
    cdss, tas, tsss, pros, seq = read_data(gff_file, ta_file, tss_file,
                                           pro_file, seq_file, hypo)
    libs, texs = read_libs(input_libs, wig_folder)
    wig_fs = read_wig(wig_f_file, "+", libs)
    wig_rs = read_wig(wig_r_file, "-", libs)
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
            if (inter["strain"] == ta.seq_id) and (
                inter["strand"] == ta.strand):
                class_utr(inter, ta, utrs, srnas, tsss, pros, wig_fs, wig_rs,
                          fuzzys, fuzzy_end, decrease, min_len, max_len)
    covers = get_utr_coverage(utrs)
    coverages = {"5utr": utr_coverages[0], "3utr": utr_coverages[1],
                 "interCDS": utr_coverages[2]}
    if notex is not None:
        cover_notex = {"5utr": notex[0], "3utr": notex[1],
                       "interCDS": notex[2]}
    else:
        cover_notex = None
    mediandict = set_cutoff(covers, coverages, cover_notex, texs)
    print_median(out_folder, mediandict)
    detect_srna(srnas, mediandict, out, out_t, texs, coverages,
                tex_notex, replicates, table_best, max_len, min_len)
