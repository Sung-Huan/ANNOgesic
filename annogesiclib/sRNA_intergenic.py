import sys
import os
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.coverage_detection import coverage_comparison, replicate_comparison
from annogesiclib.lib_reader import read_wig, read_libs
from annogesiclib.gen_TSS_type import compare_tss_cds, fix_primary_type


def get_differential_cover(fuzzy_end, num, checks, cover_sets,
                           poss, cover, decrease):
    go_out = False
    if checks["detect_diff"]:
        if (num == fuzzy_end) or (
            cover_sets["diff"] == 0) or (
            (cover["coverage"] > cover_sets["diff"]) and (
             cover["coverage"] / cover_sets["diff"]) > (1 + decrease)):
            poss["stop_point"] = cover["pos"]
            go_out = True
        elif (cover["coverage"] <= cover_sets["diff"]):
            if (cover["coverage"] / cover_sets["diff"]) <= (decrease / 2):
                num += 1
            else:
                num = 0
            cover_sets["diff"] = cover["coverage"]
            cover_sets["low"] = cover["coverage"]
        elif (cover["coverage"] > cover_sets["diff"]) and (
             (cover["coverage"] / cover_sets["diff"]) <= (1 + decrease)):
            num += 1
    if (not checks["first"]) and (cover_sets["high"] > 0):
        if ((cover_sets["low"] / cover_sets["high"]) < decrease) and (
            cover_sets["low"] > -1):
            checks["detect_diff"] = True
            cover_sets["diff"] = cover["coverage"]
    return go_out

def check_coverage_pos(start, end, cover, cutoff_coverage, tmps, cover_sets,
                       checks, poss, strand, tolerance):
    go_out = False
    if (start <= cover["pos"]) and (
        end >= cover["pos"]):
        if cover["coverage"] > cutoff_coverage:
            tmps["pos"] = cover["pos"]
            tmps["toler"] = 0
            cover_sets["total"] = cover_sets["total"] + \
                                  cover["coverage"]
            if tmps["total"] != 0:
                 cover_sets["total"] = cover_sets["total"] + \
                                       tmps["total"]
            tmps["total"] = 0
            tmps["toler"] = 0
            checks["first"] = coverage_comparison(
                              cover, cover_sets, poss,
                              checks["first"], strand)
        else:
            if tmps["toler"] <= tolerance:
                tmps["toler"] += 1
                tmps["total"] = cover_sets["total"] + \
                            cover["coverage"]
                checks["first"] = coverage_comparison(
                                  cover, cover_sets, poss,
                                  checks["first"], strand)
            else:
                if tmps["pos"] != 0:
                    poss["stop_point"] = tmps["pos"]
                    go_out = True
                tmps["total"] = 0
                tmps["toler"] = 0
    else:
        if (strand == "+") and (cover["pos"] > end):
            poss["stop_point"] = cover["pos"]
            tmps["total"] = 0
            tmps["toler"] = 0
            go_out = True
        elif (strand == "-") and (cover["pos"] < start):
            poss["stop_point"] = cover["pos"]
            tmps["total"] = 0
            tmps["toler"] = 0
            go_out = True
    return go_out, tmps

def get_median(medians):
    sortedLst = sorted(medians)
    lstLen = len(medians)
    index = (lstLen - 1) // 2
    if lstLen != 0:
        if (lstLen % 2):
            return sortedLst[index]
        else:
            return (sortedLst[index] + sortedLst[index + 1])/2.0
    else:
        return 0

def get_best(wigs, strain, strand, start, end, type_, decrease,
             cutoff_coverage, fuzzy_end, tolerance):
    cover_sets = {"low": -1, "high": -1, "total": 0, "diff": 0}
    poss = {"high": 0, "low": 0, "stop_point": -1}
    tmps = {"total": 0, "toler": 0, "pos": 0}
    srna_covers = {}
    for wig_strain, conds in wigs.items():
        if wig_strain == strain:
            for cond, tracks in conds.items():
                srna_covers[cond] = []
                for track, covers in tracks.items():
                    cover_sets["total"] = 0
                    cover_sets["diff"] = 0
                    checks = {"first": True, "detect_diff": False}
                    num = 0
                    if strand == "+":
                        covers = covers[start-2:end+1]
                    elif strand == "-":
                        covers = reversed(covers[start-2:end+1])
                    go_out = False
                    tmps = {"total": 0, "toler": 0, "pos": 0}
                    medians = []
                    for cover in covers:
                        if (cover["strand"] == strand):
                            medians.append(cover["coverage"])
                            go_out, tmps = check_coverage_pos(
                                           start, end, cover,
                                           cutoff_coverage, tmps,
                                           cover_sets, checks, poss,
                                           strand, tolerance)
                            if go_out:
                                break
                            if type_ == "differential":
                                go_out = get_differential_cover(
                                             fuzzy_end, num, checks, cover_sets,
                                             poss, cover, decrease)
                                if go_out:
                                    break
                    if strand == "+":
                        diff = poss["stop_point"] - start
                    else:
                        diff = end - poss["stop_point"]
                    avg = cover_sets["total"] / float(diff + 1)
                    median = get_median(medians)
                    if avg > float(cutoff_coverage):
                        srna_covers[cond].append({"track": track,
                                                  "high": cover_sets["high"],
                                                  "low": cover_sets["low"],
                                                  "avg": avg,
                                                  "pos": poss["stop_point"],
                                                  "type": cover["type"]})
    return srna_covers

def get_attribute_string(srna_datas, tss_pro, num, name, srna_type):
    attribute_string = ";".join(
                    ["=".join(items) for items in (["ID", "srna" + str(num)],
                     ["Name", "_".join(["sRNA", name])],
                     ["sRNA_type", srna_type])])
    datas = tss_pro.split(";")
    tss = ""
    pro = ""
    for data in datas:
        if "TSS" in data:
            if len(tss) == 0:
                tss = tss + data
            else:
                tss = ";".join([tss, data])
        elif "Cleavage" in data:
            if len(pro) == 0:            
                pro = pro + data
            else:
                pro = ";".join([pro, data])
    if len(tss) == 0:
        tss = "NA"
    if len(pro) == 0:
        pro = "NA"
    with_tss = "=".join(["with_TSS", tss])
    with_pro = "=".join(["end_cleavage", pro])
    if srna_datas is None:
        if (tss != "NA") and (pro != "NA"):
            attribute_string = ";".join([attribute_string, with_tss, with_pro])
        elif (tss != "NA"):
            attribute_string = ";".join([attribute_string, with_tss])
        elif (pro != "NA"):
            attribute_string = ";".join([attribute_string, with_pro])
    else:
        srna_data_string = ";".join(
                ["=".join(items) for items in (
                 ["best_avg_coverage", str(srna_datas["best"])],
                 ["best_high_coverage", str(srna_datas["high"])],
                 ["best_low_coverage", str(srna_datas["low"])])])
        if (tss != "NA") and (pro != "NA"):
            attribute_string = ";".join([attribute_string, with_tss, with_pro,
                                         srna_data_string])
        elif (tss != "NA"):
            attribute_string = ";".join([attribute_string, with_tss,
                                         srna_data_string])
        elif (pro != "NA"):
            attribute_string = ";".join([attribute_string, with_pro,
                                         srna_data_string])
        else:
            attribute_string = ";".join([attribute_string, srna_data_string])
    return attribute_string

def print_file(string, nums, tss, output, out_table, srna_datas,
               table_best, srna_type):
    name = '%0*d' % (5, nums["uni"])
    datas = string.split("\t")
    if (srna_datas is None):
        out_table.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t".format(
                        datas[0], name, datas[3], datas[4], datas[6], "NA", "NA",
                        "NA", "NA", "NA"))
        attribute_string = get_attribute_string(srna_datas, tss,
                                                nums["uni"], name, srna_type)
        output.write("\t".join([string, attribute_string]) + "\n")
        out_table.write(tss + "\n")
    else:
        out_table.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t".format(
                        datas[0], name, datas[3], datas[4], datas[6],
                        ";".join(srna_datas["conds"].keys()),
                        ";".join(srna_datas["conds"].values()),
                        srna_datas["best"], srna_datas["high"], srna_datas["low"]))
        attribute_string = get_attribute_string(srna_datas, tss,
                                                nums["uni"], name, srna_type)
        output.write("\t".join([string, attribute_string]) + "\n")
        if srna_datas["detail"] is not None:
            out_table.write(tss + "\t")
            if not table_best:
                first = True
                for data in srna_datas["detail"]:
                    if first:
                        out_table.write("{0}(avg={1};high={2};low={3})".format(
                                        data["track"], data["avg"],
                                        data["high"], data["low"]))
                        first = False
                    else:
                        out_table.write(";{0}(avg={1};high={2};low={3})".format(
                                        data["track"], data["avg"],
                                        data["high"], data["low"]))
            else:
                out_table.write("{0}(avg={1};high={2};low={3})".format(
                                srna_datas["track"], srna_datas["best"],
                                srna_datas["high"], srna_datas["low"]))
        out_table.write("\n")
    nums["uni"] += 1

def get_coverage(start, end, strain, wigs, strand, ta, nums, tss, output,
                 template_texs, out_table, cutoff_coverage, tex_notex,
                 replicates, decrease, fuzzy_end, table_best, tolerance,
                 texs, notex):
    srna_covers = get_best(wigs, strain, strand, start, end, "total",
                           decrease, cutoff_coverage, fuzzy_end, tolerance)
    srna_datas = replicate_comparison(srna_covers, template_texs, strand,
                                      cutoff_coverage, tex_notex, replicates,
                                      "normal", None, None, None, texs, notex)
    string = ("\t".join([str(field) for field in [
                             ta.seq_id, "ANNOgesic", "sRNA", str(start),
                             str(end), ".", ta.strand, "."]]))
    if srna_datas["best"] != 0:
        print_file(string, nums, tss, output, out_table, srna_datas,
                   table_best, ta.attributes["sRNA_type"])

def check_pro(pros, ta, start, end, min_len, max_len, srna_datas, type_,
              cutoff_coverage, wigs, decrease, fuzzy_end, template_texs,
              tex_notex, replicates, tolerance, texs, notex):
    pro_pos = -1
    detect_pro = "NA"
    for pro in pros:
        if (pro.seq_id == ta.seq_id) and (
            pro.strand == ta.strand):
            if ta.strand == "+":
                if (pro.start >= ta.start) and (pro.start <= ta.end) and (
                    (pro.start - start) >= min_len) and (
                    (pro.start - start) <= max_len):
                    pro_pos = pro.start
                    detect_pro = "".join(["Cleavage:", str(pro.start),
                                 "_", pro.strand])
                if pro.start > ta.end:
                    break
            if ta.strand == "-":
                if (pro.start >= ta.start) and (pro.start <= ta.end) and (
                    (end - pro.start) >= min_len) and (
                    (end - pro.start) <= max_len):
                    pro_pos = pro.start
                    detect_pro = "".join(["Cleavage:", str(pro.start),
                                 "_", pro.strand])
                    break
                if pro.start > ta.end:
                    break
    new_srna_datas = None
    if ta.strand == "+":
        if ((type_ == "within") and (srna_datas["pos"] < pro_pos)) or (
            (type_ == "longer") and (pro_pos != -1)):
            srna_covers = get_best(wigs, ta.seq_id, ta.strand, start, pro_pos,
                          "total", decrease, cutoff_coverage, fuzzy_end, tolerance)
            new_srna_datas = replicate_comparison(srna_covers, template_texs,
                            ta.strand, cutoff_coverage, tex_notex,
                            replicates, "normal", None, None, None, texs, notex)
            if new_srna_datas["best"] <= cutoff_coverage:
                new_srna_datas = None
    else:
        if ((type_ == "within") and (srna_datas["pos"] > pro_pos)) or (
            (type_ == "longer") and (pro_pos != -1)):
            srna_covers = get_best(wigs, ta.seq_id, ta.strand, pro_pos, end,
                          "total", decrease, cutoff_coverage, fuzzy_end, tolerance)
            new_srna_datas = replicate_comparison(srna_covers, template_texs,
                            ta.strand, cutoff_coverage, tex_notex,
                            replicates, "normal", None, None, None, texs, notex)
            if new_srna_datas["best"] <= cutoff_coverage:
                new_srna_datas = None
    return pro_pos, new_srna_datas, detect_pro

def exchange_to_pro(srna_datas, min_len, max_len, pros, ta, start, end, nums,
                    cutoff_coverage, output, out_table, table_best, wigs,
                    decrease, fuzzy_end, template_texs, tex_notex, replicates,
                    tolerance, texs, notex):
    detect = False
    if srna_datas["high"] != 0:
        if ((srna_datas["pos"] - start) >= min_len) and (
            (srna_datas["pos"] - start) <= max_len):
            pro_pos, pro_datas, pro = check_pro(pros, ta, start, end,
                                      min_len, max_len, srna_datas, "within",
                                      cutoff_coverage, wigs, decrease,
                                      fuzzy_end, template_texs, tex_notex,
                                      replicates, tolerance, texs, notex)
            if pro_datas is not None:
                srna_datas = pro_datas
                srna_datas["pos"] = pro_pos
                detect = True
            else:
                if srna_datas["best"] > cutoff_coverage:
                    detect = True
        else:
            pro_pos, pro_datas, pro = check_pro(pros, ta, start, end,
                                      min_len, max_len, srna_datas, "longer",
                                      cutoff_coverage, wigs, decrease,
                                      fuzzy_end, template_texs, tex_notex,
                                      replicates, tolerance, texs, notex)
            if pro_datas is not None:
                srna_datas = pro_datas
                srna_datas["pos"] = pro_pos
                detect = True
    else:
        pro = None
    return detect, srna_datas, pro

def detect_wig_pos(wigs, ta, start, end, nums, output, tss, template_texs,
                   out_table, cutoff_coverage, min_len, max_len, decrease,
                   fuzzy_end, table_best, tex_notex, replicates, pros,
                   tolerance, texs, notex):
    srna_covers = get_best(wigs, ta.seq_id, ta.strand, start, end,
                           "differential", decrease, cutoff_coverage,
                           fuzzy_end, tolerance)
    srna_datas = replicate_comparison(srna_covers, template_texs, ta.strand,
                                      cutoff_coverage, tex_notex, replicates,
                                      "normal", None, None, None, texs, notex)
    detect, srna_datas, pro = exchange_to_pro(srna_datas, min_len, max_len,
                              pros, ta, start, end, nums, cutoff_coverage,
                              output, out_table, table_best, wigs, decrease,
                              fuzzy_end, template_texs, tex_notex, replicates,
                              tolerance, texs, notex)
    if ta.strand == "+":
        if detect:
            string = ("\t".join([str(field) for field in [
                      ta.seq_id, "ANNOgesic", "sRNA", str(start),
                      str(srna_datas["pos"]), ".", ta.strand, "."]]))
            if pro != "NA":
                tss = ";".join([tss, pro])
            print_file(string, nums, tss, output, out_table, srna_datas,
                       table_best, ta.attributes["sRNA_type"])
    else:
        if detect:
            string = ("\t".join([str(field) for field in [
                      ta.seq_id, "ANNOgesic", "sRNA", str(srna_datas["pos"]),
                      str(end), ".", ta.strand, "."]]))
            if pro != "NA":
                tss = ";".join([tss, pro])
            print_file(string, nums, tss, output, out_table, srna_datas,
                       table_best, ta.attributes["sRNA_type"])

def detect_longer(tsss, ta, nums, output, wigs_f, wigs_r, template_texs,
                  out_table, fuzzy, cutoff_coverage, tex_notex, replicates,
                  decrease, fuzzy_end, table_best, min_len, max_len,
                  detects, pros, tolerance, texs, cut_notex):
    notex = None
    if len(tsss) != 0:
        for tss in tsss:
            cutoff = get_tss_type(tss, cutoff_coverage)
            if cut_notex is not None:
                notex = get_tss_type(tss, cut_notex)
            if cutoff is not None:
                if (tss.strand == ta.strand) and (
                    tss.seq_id == ta.seq_id):
                    if (tss.strand == "+"):
                        compare_ta_tss(tss.start, ta.start - fuzzy, ta.end, min_len,
                            max_len, wigs_f, nums, ta, tss, output, out_table,
                            template_texs, detects, (ta.end - tss.start),
                            cutoff, tex_notex, replicates, decrease,
                            fuzzy_end, table_best, tolerance, texs, notex)
                        if (tss.start >= ta.start - fuzzy) and (
                            tss.start <= ta.end) and (
                           (ta.end - tss.start) > max_len):
                            if len(wigs_f) != 0:
                                detect_wig_pos(wigs_f, ta, tss.start, ta.end, nums,
                                      output, "".join(["TSS:", str(tss.start),
                                                       "_", tss.strand]),
                                      template_texs, out_table, cutoff,
                                      min_len, max_len, decrease, fuzzy_end,
                                      table_best, tex_notex, replicates, pros,
                                      tolerance, texs, notex)
                    else:
                        compare_ta_tss(tss.end, ta.start, ta.end + fuzzy, min_len,
                            max_len, wigs_r, nums, ta, tss, output, out_table,
                            template_texs, detects, (tss.end - ta.start),
                            cutoff, tex_notex, replicates, decrease,
                            fuzzy_end, table_best, tolerance, texs, notex)
                        if (tss.end >= ta.start) and (
                            tss.end <= ta.end + fuzzy) and (
                            tss.end - ta.start > max_len):
                            if len(wigs_r) != 0:
                                detect_wig_pos(wigs_r, ta, ta.start, tss.end, nums,
                                      output, "".join(["TSS:", str(tss.end),
                                                       "_", tss.strand]),
                                      template_texs, out_table, cutoff,
                                      min_len, max_len, decrease, fuzzy_end,
                                      table_best, tex_notex, replicates, pros,
                                      tolerance, texs, notex)
    if len(tsss) == 0:
        if (len(wigs_f) != 0) and (len(wigs_r) != 0):
            if ta.strand == "+":
                num_uni = detect_wig_pos(wigs_f, ta, ta.start, ta.end, nums,
                          output, "NA", template_texs, out_table, cutoff,
                          min_len, max_len, decrease, fuzzy_end, table_best,
                          tex_notex, replicates, pros, tolerance, texs, notex)
            else:
                num_uni = detect_wig_pos(wigs_r, ta, ta.start, ta.end, nums,
                          output, "NA", template_texs, out_table, cutoff,
                          min_len, max_len, decrease, fuzzy_end, table_best,
                          tex_notex, replicates, pros, tolerance, texs, notex)

def get_tss_type(tss, cutoff_coverage):
    types = []
    for type_, cover in cutoff_coverage.items():
        if cover is not None:
            types.append(type_)
    cover = None
    if ("type" in tss.attributes.keys()):
        for type_ in types:
            if (type_ in tss.attributes["type"].lower()):
                if cover is None:
                    cover = cutoff_coverage[type_]
                elif cover < cutoff_coverage[type_]:
                    cover = cutoff_coverage[type_]
    else:
        cover = cutoff_coverage["no_tss"]
    return cover

def compare_ta_tss(tss_pos, ta_start, ta_end, min_len, max_len, wigs, nums,
                   ta, tss, output, out_table, template_texs, detects, diff,
                   cutoff_coverage, tex_notex, replicates, decrease,
                   fuzzy_end, table_best, tolerance, texs, notex):
    if (tss_pos >= ta_start) and (tss_pos <= ta_end) and (
        diff >= min_len) and (diff <= max_len):
        if tss.strand == "+":
            start = tss_pos
            end = ta_end
        else:
            start = ta_start
            end = tss_pos
        if len(wigs) != 0:
            get_coverage(start, end, ta.seq_id, wigs, tss.strand, ta, nums,
                         "".join(["TSS:", str(tss.start), "_", tss.strand]),
                         output, template_texs, out_table, cutoff_coverage,
                         tex_notex, replicates, decrease, fuzzy_end, table_best,
                         tolerance, texs, notex)
        else:
            string = "\t".join([str(field) for field in [
                     ta.seq_id, "ANNOgesic", "sRNA", str(start),
                     str(end), ta.score, ta.strand, ta.phase]])
            print_file(string, nums, "".join(["TSS:", str(tss.start),
                       "_", tss.strand]), output, out_table, None,
                       table_best, ta.attributes["sRNA_type"])
        if detects is not None:
            detects["uni_with_tss"] = True

def detect_include_tss(nums, tsss, ta, wigs_f, wigs_r, output, out_table,
                       template_texs, min_len, max_len, detects, fuzzy,
                       cutoff_coverage, tex_notex, replicates, decrease,
                       fuzzy_end, table_best, tolerance, texs, cut_notex):
    detects["uni_with_tss"] = False
    notex = None
    for tss in tsss:
        cutoff = get_tss_type(tss, cutoff_coverage)
        if cut_notex is not None:
            notex = get_tss_type(tss, cut_notex)
        if cutoff is not None:
            if (tss.strand == ta.strand) and (tss.seq_id == ta.seq_id):
                if (tss.strand == "+"):
                    compare_ta_tss(tss.start, ta.start - fuzzy, ta.end, min_len,
                                   max_len, wigs_f, nums, ta, tss, output,
                                   out_table, template_texs, detects,
                                   (ta.end - tss.start), cutoff,
                                   tex_notex, replicates, decrease,
                                   fuzzy_end, table_best, tolerance, texs, notex)
                    if tss.start > ta.end:
                        break
                else:
                    compare_ta_tss(tss.end, ta.start, ta.end + fuzzy, min_len,
                                   max_len, wigs_r, nums, ta, tss, output,
                                   out_table, template_texs, detects,
                                   (tss.end - ta.start), cutoff,
                                   tex_notex, replicates, decrease,
                                   fuzzy_end, table_best, tolerance, texs, notex)
                    if tss.end > ta.end + fuzzy:
                        break
    if not detects["uni_with_tss"]:
        if (ta.strand == "+") and (len(wigs_f) != 0):
            get_coverage(ta.start, ta.end, ta.seq_id, wigs_f, "+", ta, nums,
                         "False", output, template_texs, out_table,
                         cutoff_coverage["no_tss"], tex_notex, replicates, decrease,
                         fuzzy_end, table_best, tolerance, texs, notex)
        elif (ta.strand == "-") and (len(wigs_r) != 0):
            get_coverage(ta.start, ta.end, ta.seq_id, wigs_r, "-", ta, nums,
                         "False", output, template_texs, out_table,
                         cutoff_coverage["no_tss"], tex_notex, replicates, decrease,
                         fuzzy_end, table_best, tolerance, texs, notex)
        elif (len(wigs_f) == 0) and (len(wigs_r) == 0):
            print_file(ta.info_without_attributes.replace("Transcript", "sRNA"),
                       nums, "False", output, out_table, None, table_best,
                       ta.attributes["sRNA_type"])

def get_proper_tss(tss_file, cutoff_coverage):
    types = []
    gff_parser = Gff3Parser()
    for type_, cover in cutoff_coverage.items():
        if cover is not None:
            types.append(type_)
    tsss = []
    num_tss = 0
    if tss_file is not None:
        tss_f = open(tss_file, "r")
        for entry in gff_parser.entries(tss_f):
            if ("type" in entry.attributes.keys()):
                for type_ in types:
                    if (type_ in entry.attributes["type"].lower()):
                        tsss.append(entry)
                        num_tss += 1
                        break
            else:
                tsss.append(entry)
                num_tss += 1
        tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start))
        tss_f.close()
    return tsss, num_tss

def read_data(gff_file, tran_file, pro_file, hypo):
    cdss = []
    tas = []
    pros = []
    genes = []
    num_cds = 0
    num_ta = 0
    num_pro = 0
    gff_parser = Gff3Parser()
    g_f = open(gff_file, "r")
    for entry in gff_parser.entries(g_f):
        if (entry.feature == "CDS") or (
            entry.feature == "pCDS") or (
            entry.feature == "tRNA") or (
            entry.feature == "rRNA"):
            if ("product" in entry.attributes.keys()) and (hypo):
                if "hypothetical protein" not in entry.attributes["product"]:
                    cdss.append(entry)
                    num_cds += 1
            else:
                cdss.append(entry)
                num_cds += 1
        if (entry.feature == "gene"):
            genes.append(entry)
    if pro_file is not None:
        pro_f = open(pro_file, "r")
        for entry in gff_parser.entries(pro_f):
            pros.append(entry)
            num_pro += 1
        pros = sorted(pros, key=lambda k: (k.seq_id, k.start))
        pro_f.close()
    t_h = open(tran_file)
    for entry_ta in gff_parser.entries(t_h):
        tas.append(entry_ta)
        num_ta += 1
    nums = {"cds": num_cds, "ta": num_ta,
            "pro": num_pro, "uni": 0}
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start))
    tas = sorted(tas, key=lambda k: (k.seq_id, k.start))
    genes = sorted(genes, key=lambda k: (k.seq_id, k.start))
    g_f.close()
    t_h.close()
    return nums, cdss, tas, pros, genes

def read_tss(tss_file, cutoff_coverage, modify):
    if modify:
        tsss = []
        tss_f = open(tss_file, "r")
        gff_parser = Gff3Parser()
        for entry in gff_parser.entries(tss_f):
            tsss.append(entry)
        num_tss = None
        tss_f.close()
    else:
        tsss, num_tss = get_proper_tss(tss_file, cutoff_coverage)
    return tsss, num_tss

def compare_ta_cds(cdss, ta, detects):
    for cds in cdss:
        if (cds.strand == ta.strand) and (
            cds.seq_id == ta.seq_id):
            if ((cds.end < ta.end) and (
                 cds.end > ta.start) and (
                 cds.start <= ta.start)) or (
                (cds.start > ta.start) and (
                 cds.start < ta.end) and (
                 cds.end >= ta.end)) or (
                (cds.end >= ta.end) and (
                 cds.start <= ta.start)) or (
                (cds.end <= ta.end) and (
                 cds.start >= ta.start)):
                detects["overlap"] = True
                ta.attributes["sRNA_type"] = "in_CDS"
    if not detects["overlap"]:
        ta.attributes["sRNA_type"] = "intergenic"

def check_srna_condition(ta, min_len, max_len, tsss, pros, wigs_f, wigs_r,
                         nums, output, out_table, texs, fuzzy, detects,
                         cutoff_coverage, tex_notex, replicates,
                         decrease, fuzzy_end, table_best, tolerance, notex):
    if ((ta.end - ta.start) >= min_len) and (
        (ta.end - ta.start) <= max_len):
        if len(tsss) != 0:
            detect_include_tss(nums, tsss, ta, wigs_f, wigs_r, output,
                               out_table, texs, min_len, max_len, detects,
                               fuzzy, cutoff_coverage, tex_notex, replicates,
                               decrease, fuzzy_end, table_best, tolerance,
                               texs, notex)
        else:
            if (ta.strand == "+") and (len(wigs_f) != 0):
                get_coverage(ta.start, ta.end, ta.seq_id, wigs_f, "+",
                             ta, nums, "NA", output, texs, out_table,
                             cutoff_coverage["no_tss"], tex_notex, replicates,
                             decrease, fuzzy_end, table_best, tolerance, texs,
                             notex)
            elif (ta.strand == "-") and (len(wigs_r) != 0):
                get_coverage(ta.start, ta.end, ta.seq_id, wigs_r, "-",
                             ta, nums, "NA", output, texs, out_table,
                             cutoff_coverage["no_tss"], tex_notex, replicates,
                             decrease, fuzzy_end, table_best, tolerance, texs,
                             notex)
            if (len(wigs_f) == 0) and (len(wigs_r) == 0):
                print_file(ta.info_without_attributes.replace("Transcript", "sRNA"),
                           nums, "NA", output, out_table, None, table_best,
                           ta.attributes["sRNA_type"])
    if ((ta.end - ta.start) > max_len):
        detect_longer(tsss, ta, nums, output, wigs_f, wigs_r, texs,
                      out_table, fuzzy, cutoff_coverage, tex_notex,
                      replicates, decrease, fuzzy_end, table_best,
                      min_len, max_len, detects, pros, tolerance, texs, notex)

def get_cutoff(cutoffs, out_folder, file_type):
    out = open(os.path.join(out_folder, "tmp_cutoff_inter"), "a")
    coverages = {}
    num_cutoff = 0
    for cutoff in cutoffs:
        if (cutoff != "0") and (num_cutoff == 0):
            coverages["primary"] = float(cutoff)
        elif (cutoff != "0") and (num_cutoff == 1):
            coverages["secondary"] = float(cutoff)
        elif (cutoff != "0") and (num_cutoff == 2):
            coverages["internal"] = float(cutoff)
        elif (cutoff != "0") and (num_cutoff == 3):
            coverages["antisense"] = float(cutoff)
        elif (cutoff != "0") and (num_cutoff == 4):
            coverages["orphan"] = float(cutoff)
        num_cutoff += 1
    low = None
    for cover in coverages.values():
        if cover != 0:
            if low is None:
                low = cover
            elif cover < float(low):
                low = cover
    if low is None:
        print("Error:The cutoff of coverage can not be all 0...")
        sys.exit()
    coverages["no_tss"] = float(low)
    for tss, cover in coverages.items():
        out.write("\t".join([file_type, tss, str(cover)]) + "\n")
    out.close()
    return coverages

def modify_wigs_for_tss_type(wigs, strand):
    new_wigs = {}
    for strain, conds in wigs.items():
        if strain not in new_wigs.keys():
            new_wigs[strain] = {}
        for cond, tracks in conds.items():
            for track, covers in tracks.items():
                if track not in new_wigs[strain].keys():
                    new_wigs[strain][track] = []
                for cover in covers:
                    new_wigs[strain][track].append({
                                    "pos": cover["pos"],
                                    "coverage": cover["coverage"],
                                    "strand": strand})
    return new_wigs

def compute_tss_type(tss_file, cutoff_coverage, out_folder, prefix, cdss, genes,
                     wigs_f, wigs_r):
    tsss, num_tss = read_tss(tss_file, cutoff_coverage, True)
    if "TSS_class" not in os.listdir(out_folder):
        os.mkdir(os.path.join(out_folder, "TSS_class"))
    new_tss_file = os.path.join(out_folder, "TSS_class",
                                "_".join([prefix, "TSS.gff"]))
    new_tss_fh = open(new_tss_file, "w")
    num_tss = 0
    for tss in tsss:
        tss_type = compare_tss_cds(tss, cdss, genes)

        tss.attributes = tss_type[1]
        tss.attributes["ID"] = "tss" + str(num_tss)
        tss.attribute_string = "".join([tss_type[0], ";ID=tss", str(num_tss)])
        num_tss += 1
    wigs_fm = modify_wigs_for_tss_type(wigs_f, "+")
    wigs_rm = modify_wigs_for_tss_type(wigs_r, "-")
    final_tsss = fix_primary_type(tsss, wigs_fm, wigs_rm)
    for tss in final_tsss:
        name = ">" + "_".join([tss.seq_id, str(tss.start), tss.strand])
        tss.attribute_string = ";".join(
            ["=".join(items) for items in tss.attributes.items()])
        new_tss_fh.write("\t".join([str(field) for field in [
                         tss.seq_id, tss.source, tss.feature, tss.start,
                         tss.end, tss.score, tss.strand, tss.phase,
                         tss.attribute_string]]) + "\n")
    new_tss_fh.close()
    wigs_fm = {}
    wigs_rm = {}
    tss_file = new_tss_file

def intergenic_srna(gff_file, tran_file, tss_file, pro_file, fuzzy, max_len,
                    min_len, wig_f_file, wig_r_file, wig_folder, input_libs,
                    tex_notex, replicates, output_file, output_table, table_best,
                    decrease, fuzzy_end, cutoffs, tolerance, in_cds, hypo,
                    tss_source, prefix, out_folder, file_type, cut_notex):
    cutoff_coverage = get_cutoff(cutoffs, out_folder, file_type)
    notex = None
    if cut_notex is not None:
        notex = get_cutoff(cut_notex, out_folder, "notex")
    libs, texs = read_libs(input_libs, wig_folder)
    wigs_f = read_wig(wig_f_file, "+", libs)
    wigs_r = read_wig(wig_r_file, "-", libs)
    nums, cdss, tas, pros, genes = read_data(gff_file, tran_file,
                                             pro_file, hypo)
    if not tss_source:
        compute_tss_type(tss_file, cutoff_coverage, out_folder, prefix, cdss,
                         genes, wigs_f, wigs_r)
    tsss, num_tss = read_tss(tss_file, cutoff_coverage, False)
    detects = {"overlap": False, "uni_with_tss": False}
    output = open(output_file, "w")
    out_table = open(output_table, "w")
    output.write("##gff-version 3\n")
    for ta in tas:
        detects["overlap"] = False
        compare_ta_cds(cdss, ta, detects)
        if (detects["overlap"]) and (not in_cds):
            continue
        else:
            check_srna_condition(ta, min_len, max_len, tsss, pros, wigs_f, wigs_r,
                                 nums, output, out_table, texs, fuzzy, detects,
                                 cutoff_coverage, tex_notex, replicates,
                                 decrease, fuzzy_end, table_best, tolerance, notex)
    file_name = output_file.split(".")
    file_name = file_name[0] + ".stat"
    output.close()
    out_table.close()
