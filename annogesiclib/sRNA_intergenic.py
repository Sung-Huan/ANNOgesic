import sys
import os
import gc
import pandas
import numpy as np
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.coverage_detection import coverage_comparison
from annogesiclib.coverage_detection import replicate_comparison, get_repmatch
from annogesiclib.lib_reader import read_wig, read_libs
from annogesiclib.gen_TSS_type import compare_tss_cds, fix_primary_type
from annogesiclib.args_container import ArgsContainer


def get_differential_cover(num, checks, cover_sets, poss, cover,
                           args_srna, cover_pos):
    go_out = False
    if checks["detect_diff"]:
        if (num == args_srna.fuzzy_inter) or (
            cover_sets["diff"] == 0) or (
            (cover > cover_sets["diff"]) and (
             cover / cover_sets["diff"]) > (
             1 + args_srna.decrease_inter)):
            poss["stop_point"] = cover_pos
            go_out = True
        elif (cover <= cover_sets["diff"]):
            if (cover / cover_sets["diff"]) <= (
                    args_srna.decrease_inter / 2):
                num += 1
            else:
                num = 0
            cover_sets["diff"] = cover
            cover_sets["low"] = cover
        elif (cover > cover_sets["diff"]) and (
                (cover / cover_sets["diff"]) <= (
                1 + args_srna.decrease_inter)):
            num += 1
    if (not checks["first"]) and (cover_sets["high"] > 0):
        if ((cover_sets["low"] / cover_sets["high"]) <
                args_srna.decrease_inter) and (
                cover_sets["low"] > -1):
            checks["detect_diff"] = True
            cover_sets["diff"] = cover
    return go_out


def check_coverage_pos(start, end, cover, cutoff_coverage, cover_sets,
                       checks, poss, strand, cover_pos):
    go_out = False
    if (start <= cover_pos) and (
            end >= cover_pos):
        if cover > cutoff_coverage:
            cover_sets["total"] = cover_sets["total"] + cover
            checks["first"] = coverage_comparison(
                              cover, cover_sets, poss,
                              checks["first"], strand, cover_pos)
        else:
            cover_sets["total"] = cover_sets["total"] + cover
            checks["first"] = coverage_comparison(
                                cover, cover_sets, poss,
                                checks["first"], strand, cover_pos)
    else:
        if (strand == "+") and (cover_pos > end):
            poss["stop_point"] = cover_pos
            go_out = True
        elif (strand == "-") and (cover_pos < start):
            poss["stop_point"] = cover_pos
            go_out = True
    return go_out


def check_start_and_end(start, end, covers):
    if (start - 2) < 0:
        c_start = 0
    else:
        c_start = start - 2
    if (end + 2) > len(covers):
        c_end = len(covers)
    else:
        c_end = end + 2
    return c_start, c_end


def get_best(wigs, strain, strand, start, end, type_, args_srna, cutoff):
    cover_sets = {"low": -1, "high": -1, "total": 0, "diff": 0}
    poss = {"high": 0, "low": 0, "stop_point": -1}
    srna_covers = {}
    for wig_strain, conds in wigs.items():
        if wig_strain == strain:
            for cond, tracks in conds.items():
                srna_covers[cond] = []
                for lib_name, covers in tracks.items():
                    track = lib_name.split("|")[-3]
                    lib_strand = lib_name.split("|")[-2]
                    lib_type = lib_name.split("|")[-1]
                    cover_sets["total"] = 0
                    cover_sets["diff"] = 0
                    checks = {"first": True, "detect_diff": False}
                    num = 0
                    c_start, c_end = check_start_and_end(start, end, covers)
                    covers = covers[c_start: c_end]
                    if strand == "-":
                        covers = covers[::-1]
                    go_out = False
                    pos = 0
                    for cover in covers:
                        if strand == "+":
                            cover_pos = c_start + pos
                        else:
                            cover_pos = c_end - pos
                        if (lib_strand == strand):
                            go_out = check_coverage_pos(
                                    start, end, cover, cutoff,
                                    cover_sets, checks, poss,
                                    strand, cover_pos)
                            if go_out:
                                break
                            if type_ == "differential":
                                go_out = get_differential_cover(
                                        num, checks, cover_sets,
                                        poss, cover, args_srna, cover_pos)
                                if go_out:
                                    break
                        pos += 1
                    if strand == "+":
                        diff = poss["stop_point"] - start
                    else:
                        diff = end - poss["stop_point"]
                    avg = cover_sets["total"] / float(diff + 1)
                    if (avg > float(cutoff)):
                        srna_covers[cond].append({"track": track,
                                                  "high": cover_sets["high"],
                                                  "low": cover_sets["low"],
                                                  "avg": avg,
                                                  "pos": poss["stop_point"],
                                                  "type": lib_type})
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


def print_file(string, tss, srna_datas, srna_type, args_srna):
    name = '%0*d' % (5, args_srna.nums["uni"])
    datas = string.split("\t")
    if (srna_datas is None):
        args_srna.out_table.write(
            "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t".format(
                datas[0], name, datas[3], datas[4], datas[6],
                "NA", "NA", "NA", "NA", "NA"))
        attribute_string = get_attribute_string(
                srna_datas, tss, args_srna.nums["uni"], name, srna_type)
        args_srna.output.write("\t".join([string, attribute_string]) + "\n")
        args_srna.out_table.write(tss + "\n")
    else:
        args_srna.out_table.write(
            "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t".format(
                datas[0], name, datas[3], datas[4], datas[6],
                ";".join(srna_datas["conds"].keys()),
                ";".join(srna_datas["conds"].values()),
                srna_datas["best"], srna_datas["high"], srna_datas["low"]))
        attribute_string = get_attribute_string(
                srna_datas, tss, args_srna.nums["uni"], name, srna_type)
        args_srna.output.write("\t".join([string, attribute_string]) + "\n")
        if srna_datas["detail"] is not None:
            args_srna.out_table.write(tss + "\t")
            if not args_srna.table_best:
                first = True
                for data in srna_datas["detail"]:
                    if first:
                        args_srna.out_table.write(
                            "{0}(avg={1};high={2};low={3})".format(
                                data["track"], data["avg"],
                                data["high"], data["low"]))
                        first = False
                    else:
                        args_srna.out_table.write(
                            ";{0}(avg={1};high={2};low={3})".format(
                                data["track"], data["avg"],
                                data["high"], data["low"]))
            else:
                args_srna.out_table.write(
                    "{0}(avg={1};high={2};low={3})".format(
                        srna_datas["track"], srna_datas["best"],
                        srna_datas["high"], srna_datas["low"]))
        args_srna.out_table.write("\n")
    args_srna.nums["uni"] += 1


def get_coverage(start, end, strain, wigs, strand, ta, tss, cutoff_coverage,
                 notex, args_srna):
    srna_covers = get_best(wigs, strain, strand, start, end, "total",
                           args_srna, cutoff_coverage)
    srna_datas = replicate_comparison(args_srna, srna_covers, strand, "normal",
                                      None, None, None, notex, cutoff_coverage,
                                      args_srna.texs)
    string = ("\t".join([str(field) for field in [
                             ta.seq_id, "ANNOgesic", "ncRNA", str(start),
                             str(end), ".", ta.strand, "."]]))
    if srna_datas["best"] != 0:
        print_file(string, tss, srna_datas, ta.attributes["sRNA_type"],
                   args_srna)


def check_pro(ta, start, end, srna_datas, type_, cutoff,
              wigs, notex, args_srna):
    '''check the processing site for long non-coding RNA'''
    pro_pos = -1
    detect_pro = "NA"
    for pro in args_srna.pros:
        if (pro.seq_id == ta.seq_id) and (
                pro.strand == ta.strand):
            if ta.strand == "+":
                if (pro.start >= ta.start) and (pro.start <= ta.end) and (
                        (pro.start - start) >= args_srna.min_len) and (
                        (pro.start - start) <= args_srna.max_len):
                    pro_pos = pro.start
                    detect_pro = "".join(["Cleavage:", str(pro.start),
                                          "_", pro.strand])
                if pro.start > ta.end:
                    break
            if ta.strand == "-":
                if (pro.start >= ta.start) and (pro.start <= ta.end) and (
                        (end - pro.start) >= args_srna.min_len) and (
                        (end - pro.start) <= args_srna.max_len):
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
                                   "total", args_srna, cutoff)
            new_srna_datas = replicate_comparison(
                    args_srna, srna_covers, ta.strand, "normal", None, None,
                    None, notex, cutoff, args_srna.texs)
            if new_srna_datas["best"] <= cutoff:
                new_srna_datas = None
    else:
        if ((type_ == "within") and (srna_datas["pos"] > pro_pos) and (
                pro_pos != -1)) or (
                (type_ == "longer") and (pro_pos != -1)):
            srna_covers = get_best(wigs, ta.seq_id, ta.strand, pro_pos, end,
                                   "total", args_srna, cutoff)
            new_srna_datas = replicate_comparison(
                    args_srna, srna_covers, ta.strand, "normal", None, None,
                    None, notex, cutoff, args_srna.texs)
            if new_srna_datas["best"] <= cutoff:
                new_srna_datas = None
    return pro_pos, new_srna_datas, detect_pro


def exchange_to_pro(args_srna, srna_datas, ta, start,
                    end, cutoff, wigs, notex):
    detect = False
    if srna_datas["high"] != 0:
        if ((srna_datas["pos"] - start) >= args_srna.min_len) and (
                (srna_datas["pos"] - start) <= args_srna.max_len):
            pro_pos, pro_datas, pro = check_pro(
                    ta, start, end, srna_datas,
                    "within", cutoff, wigs, notex, args_srna)
            if pro_datas is not None:
                srna_datas = pro_datas
                srna_datas["pos"] = pro_pos
                detect = True
            else:
                if srna_datas["best"] > cutoff:
                    detect = True
        else:
            pro_pos, pro_datas, pro = check_pro(
                    ta, start, end, srna_datas,
                    "longer", cutoff, wigs, notex, args_srna)
            if pro_datas is not None:
                srna_datas = pro_datas
                srna_datas["pos"] = pro_pos
                detect = True
    else:
        pro = None
    return detect, srna_datas, pro


def detect_wig_pos(wigs, ta, start, end, tss, cutoff, notex, args_srna):
    '''searching the coverage decrease'''
    srna_covers = get_best(wigs, ta.seq_id, ta.strand, start, end,
                           "differential", args_srna, cutoff)
    srna_datas = replicate_comparison(
            args_srna, srna_covers, ta.strand, "normal", None, None,
            None, notex, cutoff, args_srna.texs)
    detect, srna_datas, pro = exchange_to_pro(args_srna, srna_datas, ta, start,
                                              end, cutoff, wigs, notex)
    if ta.strand == "+":
        if detect:
            string = ("\t".join([str(field) for field in [
                      ta.seq_id, "ANNOgesic", "ncRNA", str(start),
                      str(srna_datas["pos"]), ".", ta.strand, "."]]))
            if pro != "NA":
                tss = ";".join([tss, pro])
            print_file(string, tss, srna_datas, ta.attributes["sRNA_type"],
                       args_srna)
    else:
        if detect:
            string = ("\t".join([str(field) for field in [
                      ta.seq_id, "ANNOgesic", "ncRNA", str(srna_datas["pos"]),
                      str(end), ".", ta.strand, "."]]))
            if pro != "NA":
                tss = ";".join([tss, pro])
            print_file(string, tss, srna_datas, ta.attributes["sRNA_type"],
                       args_srna)


def detect_longer(ta, args_srna, cdss, wigs_f, wigs_r):
    '''deal with the long non-coding RNA'''
    notex = None
    if len(args_srna.tsss) != 0:
        for tss in args_srna.tsss:
            cutoff = get_tss_type(tss, args_srna.cutoff_coverage, ta, cdss,
                                  args_srna.file_type, args_srna.break_tran)
            if notex is not None:
                notex = get_tss_type(tss, args_srna.notex, ta, cdss, "notex",
                                     args_srna.break_tran)
            if cutoff is not None:
                if (tss.strand == ta.strand) and (
                        tss.seq_id == ta.seq_id):
                    if (tss.strand == "+"):
                        compare_ta_tss(
                            tss.start, ta.start - args_srna.fuzzy, ta.end, ta,
                            tss, ta.end - tss.start, cutoff, notex,
                            wigs_f, args_srna)
                        if (tss.start >= ta.start - args_srna.fuzzy) and (
                                tss.start <= ta.end) and (
                                (ta.end - tss.start) > args_srna.max_len):
                            if len(wigs_f) != 0:
                                detect_wig_pos(
                                    wigs_f, ta, tss.start, ta.end,
                                    "".join(["TSS:", str(tss.start),
                                             "_", tss.strand]),
                                    cutoff, notex, args_srna)
                    else:
                        compare_ta_tss(
                            tss.end, ta.start, ta.end + args_srna.fuzzy, ta,
                            tss, tss.end - ta.start, cutoff, notex,
                            wigs_r, args_srna)
                        if (tss.end >= ta.start) and (
                                tss.end <= ta.end + args_srna.fuzzy) and (
                                tss.end - ta.start > args_srna.max_len):
                            if len(wigs_r) != 0:
                                detect_wig_pos(
                                    wigs_r, ta, ta.start,
                                    tss.end, "".join(["TSS:", str(tss.end),
                                                      "_", tss.strand]),
                                    cutoff, notex, args_srna)
    if len(args_srna.tsss) == 0:
        cutoff = get_tss_type(None, args_srna.cutoff_coverage, ta, cdss,
                              args_srna.file_type, args_srna.break_tran)
        if (len(wigs_f) != 0) and (len(wigs_r) != 0):
            if ta.strand == "+":
                detect_wig_pos(wigs_f, ta, ta.start, ta.end, "NA",
                               cutoff, notex, args_srna)
            else:
                detect_wig_pos(wigs_r, ta, ta.start, ta.end, "NA",
                               cutoff, notex, args_srna)


def get_tss_type(tss, cutoff_coverage, ta, cdss, file_type, break_tran):
    '''get the cutoff based on the tss type for doing core computation
    If there are multiple TSS type, it will get the high cutoff one'''
    types = []
    for type_, cover in cutoff_coverage.items():
        if cover is not None:
            types.append(type_)
    cover = None
    break_ = False
    if tss is None:
        cover = cutoff_coverage["no_tss"]
    else:
        if ("type" in tss.attributes.keys()):
            for type_ in types:
                if (type_ in tss.attributes["type"].lower()):
                    if cover is None:
                        cover = cutoff_coverage[type_]
                        cover, break_ = check_break_tran(
                                tss, ta, cdss, cover, file_type,
                                break_tran, type_)
                    elif cover < cutoff_coverage[type_]:
                        if (break_) and ((type_.lower() == "primary") or (
                            type_.lower() == "secondary")):
                            pass
                        else:
                            cover = cutoff_coverage[type_]
        else:
            cover = cutoff_coverage["no_tss"]
    return cover


def check_break_tran(tss, ta, cdss, cover, file_type, break_tran, type_):
    '''Check the Primary or Secondary TSS which associated with a 
    non-contain CDS transcript'''
    break_ = False
    if ("type" in tss.attributes.keys()):
        if ("primary" in tss.attributes["type"].lower()) or (
                "secondary" in tss.attributes["type"].lower()):
            overlap = False
            for cds in cdss:
                if (cds.seq_id == ta.seq_id) and (
                        cds.strand == ta.strand):
                    if ((cds.start <= ta.start) and (
                            cds.end >= ta.end)) or (
                            (cds.start >= ta.start) and (
                            cds.end <= ta.end)) or (
                            (cds.start <= ta.start) and (
                            cds.end <= ta.end) and (
                            cds.end >= ta.start)) or (
                            (cds.start >= ta.start) and (
                            cds.start <= ta.end) and (
                            cds.end >= ta.end)):
                        overlap = True
                    elif (cds.start > ta.end):
                        break
            if not overlap:
                break_ = True
                if (file_type == "tex"):
                    tmp_cover = break_tran[0]
                elif (file_type == "notex"):
                    tmp_cover = break_tran[1]
                elif (file_type == "frag"):
                    tmp_cover = break_tran[2]
                if (type_.lower() == "primary") or (
                        type_.lower() == "secondary"):
                    cover = tmp_cover
                else:
                    if tmp_cover > cover:
                        cover = tmp_cover
    return cover, break_


def compare_ta_tss(tss_pos, ta_start, ta_end, ta, tss, diff, cutoff_coverage,
                   notex, wigs, args_srna):
    if (tss_pos >= ta_start) and (tss_pos <= ta_end) and (
            diff >= args_srna.min_len) and (diff <= args_srna.max_len):
        if tss.strand == "+":
            start = tss_pos
            end = ta_end
        else:
            start = ta_start
            end = tss_pos
        if len(wigs) != 0:
            get_coverage(start, end, ta.seq_id, wigs, tss.strand, ta,
                         "".join(["TSS:", str(tss.start), "_", tss.strand]),
                         cutoff_coverage, notex, args_srna)
        else:
            string = "\t".join([str(field) for field in [
                     ta.seq_id, "ANNOgesic", "ncRNA", str(start),
                     str(end), ta.score, ta.strand, ta.phase]])
            print_file(string, "".join(["TSS:", str(tss.start),
                       "_", tss.strand]), None,
                       ta.attributes["sRNA_type"], args_srna)
        if args_srna.detects is not None:
            args_srna.detects["uni_with_tss"] = True


def detect_include_tss(ta, args_srna, cdss, wigs_f, wigs_r):
    '''compare sRNA with TSS'''
    args_srna.detects["uni_with_tss"] = False
    notex = None
    for tss in args_srna.tsss:
#        print("AAA")
#        print(tss.start)
#        print(tss.attributes["type"])
        cutoff = get_tss_type(tss, args_srna.cutoff_coverage, ta, cdss,
                              args_srna.file_type, args_srna.break_tran)
#        print(cutoff)
        if args_srna.notex is not None:
            notex = get_tss_type(tss, args_srna.notex, ta, cdss, "notex",
                                 args_srna.break_tran)
        if cutoff is not None:
            if (tss.strand == ta.strand) and (tss.seq_id == ta.seq_id):
                if (tss.strand == "+"):
                    compare_ta_tss(
                            tss.start, ta.start - args_srna.fuzzy, ta.end,
                            ta, tss, ta.end - tss.start, cutoff, notex,
                            wigs_f, args_srna)
                    if tss.start > ta.end:
                        break
                else:
                    compare_ta_tss(
                            tss.end, ta.start, ta.end + args_srna.fuzzy, ta,
                            tss, tss.end - ta.start, cutoff, notex,
                            wigs_r, args_srna)
                    if (tss.end > ta.end + args_srna.fuzzy):
                        break
    if not args_srna.detects["uni_with_tss"]:
        if (ta.strand == "+") and (len(wigs_f) != 0):
            get_coverage(ta.start, ta.end, ta.seq_id, wigs_f, "+",
                         ta, "False", args_srna.cutoff_coverage["no_tss"],
                         notex, args_srna)
        elif (ta.strand == "-") and (len(wigs_r) != 0):
            get_coverage(ta.start, ta.end, ta.seq_id, wigs_r, "-",
                         ta, "False", args_srna.cutoff_coverage["no_tss"],
                         notex, args_srna)
        elif (len(wigs_f) == 0) and (len(wigs_r) == 0):
            print_file(
                ta.info_without_attributes.replace("Transcript", "ncRNA"),
                "False", None, ta.attributes["sRNA_type"], args_srna)


def get_proper_tss(tss_file, cutoff_coverage):
    '''get the associated TSS'''
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
        tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
        tss_f.close()
    return tsss, num_tss


def read_data(args_srna):
    cdss = []
    tas = []
    pros = []
    genes = []
    ncs = []
    num_cds = 0
    num_ta = 0
    num_pro = 0
    gff_parser = Gff3Parser()
    g_f = open(args_srna.gff_file, "r")
    for entry in gff_parser.entries(g_f):
        if (entry.feature == "CDS") or (
                entry.feature == "pCDS") or (
                entry.feature == "tRNA") or (
                entry.feature == "rRNA"):
            if ("product" in entry.attributes.keys()) and (args_srna.hypo):
                if "hypothetical protein" not in entry.attributes["product"]:
                    cdss.append(entry)
                    num_cds += 1
            else:
                cdss.append(entry)
                num_cds += 1
        if (entry.feature == "gene"):
            genes.append(entry)
        if (entry.feature == "ncRNA"):
            ncs.append(entry)
    if args_srna.pro_file is not None:
        pro_f = open(args_srna.pro_file, "r")
        for entry in gff_parser.entries(pro_f):
            pros.append(entry)
            num_pro += 1
        pros = sorted(pros, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
        pro_f.close()
    t_h = open(args_srna.tran_file)
    for entry_ta in gff_parser.entries(t_h):
        tas.append(entry_ta)
        num_ta += 1
    nums = {"cds": num_cds, "ta": num_ta,
            "pro": num_pro, "uni": 0}
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    tas = sorted(tas, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    genes = sorted(genes, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    g_f.close()
    t_h.close()
    return nums, cdss, tas, pros, genes, ncs


def read_tss(tss_file):
    tsss = []
    if tss_file is not None:
        tss_f = open(tss_file, "r")
        gff_parser = Gff3Parser()
        for entry in gff_parser.entries(tss_f):
            tsss.append(entry)
        tss_f.close()
    num_tss = None
    return tsss, num_tss


def check_overlap(cds, ta):
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
        return True


def compare_ta_cds(cdss, ta, detects):
    for cds in cdss:
        if (cds.strand == ta.strand) and (
                cds.seq_id == ta.seq_id):
            if check_overlap(cds, ta):
                detects["overlap"] = True
                ta.attributes["sRNA_type"] = "in_CDS"
        elif (cds.strand != ta.strand) and (
              cds.seq_id == ta.seq_id):
            if check_overlap(cds, ta):
                detects["anti"] = True
                ta.attributes["sRNA_type"] = "antisense"
    if (not detects["overlap"]) and (not detects["anti"]):
        ta.attributes["sRNA_type"] = "intergenic"


def check_srna_condition(ta, args_srna, cdss, wigs_f, wigs_r):
    '''check the long sRNA and get the coverage of sRNA and check TSS'''
    if ((ta.end - ta.start) >= args_srna.min_len) and (
            (ta.end - ta.start) <= args_srna.max_len):
        if len(args_srna.tsss) != 0:
            detect_include_tss(ta, args_srna, cdss, wigs_f, wigs_r)
        else:
            if (ta.strand == "+") and (len(wigs_f) != 0):
                get_coverage(ta.start, ta.end, ta.seq_id,
                             wigs_f, "+", ta, "NA",
                             args_srna.cutoff_coverage["no_tss"],
                             args_srna.notex["no_tss"], args_srna)
            elif (ta.strand == "-") and (len(wigs_r) != 0):
                get_coverage(ta.start, ta.end, ta.seq_id,
                             wigs_r, "-", ta, "NA",
                             args_srna.cutoff_coverage["no_tss"],
                             args_srna.notex["no_tss"], args_srna)
            if (len(wigs_f) == 0) and (len(wigs_r) == 0):
                print_file(ta.info_without_attributes.replace("Transcript",
                           "ncRNA"), "NA", None, ta.attributes["sRNA_type"],
                           args_srna)
    if ((ta.end - ta.start) > args_srna.max_len):
        detect_longer(ta, args_srna, cdss, wigs_f, wigs_r)


def get_cutoff(cutoffs, out_folder, file_type):
    '''set the cutoff of intergenic and antisense sRNA'''
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


def compute_tss_type(args_srna, cdss, genes, wigs_f, wigs_r):
    tsss, num_tss = read_tss(args_srna.tss_file)
    if "TSS_class" not in os.listdir(args_srna.out_folder):
        os.mkdir(os.path.join(args_srna.out_folder, "TSS_class"))
    new_tss_file = os.path.join(args_srna.out_folder, "TSS_class",
                                "_".join([args_srna.prefix, "TSS.gff"]))
    new_tss_fh = open(new_tss_file, "w")
    num_tss = 0
    for tss in tsss:
        tss_type = compare_tss_cds(tss, cdss, genes)
        tss.attributes = tss_type[1]
        tss.attributes["ID"] = "tss" + str(num_tss)
        tss.attribute_string = "".join([tss_type[0], ";ID=tss", str(num_tss)])
        num_tss += 1
    final_tsss = fix_primary_type(tsss, wigs_f, wigs_r)
    for tss in final_tsss:
        tss.attribute_string = ";".join(
            ["=".join(items) for items in tss.attributes.items()])
        new_tss_fh.write("\t".join([str(field) for field in [
                         tss.seq_id, tss.source, tss.feature, tss.start,
                         tss.end, tss.score, tss.strand, tss.phase,
                         tss.attribute_string]]) + "\n")
    new_tss_fh.close()


def get_intergenic_antisense_cutoff(args_srna):
    '''set the cutoff of intergenic and antisense sRNA
    also deal with the no tex library'''
    cutoff_coverage = get_cutoff(args_srna.cutoffs, args_srna.out_folder,
                                 args_srna.file_type)
    notex = None
    if args_srna.cut_notex is not None:
        notex = get_cutoff(args_srna.cut_notex, args_srna.out_folder, "notex")
    return cutoff_coverage, notex


def free_memory(paras):
    for data in paras:
        del(data)
    gc.collect()


def check_ncRNA(gene, ncs):
    detect = False
    for nc in ncs:
        if "Parent" in nc.attributes.keys():
            if gene.attributes["ID"] in nc.attributes["Parent"]:
                detect = True
        elif (nc.seq_id == gene.seq_id) and (
                nc.strand == gene.strand):
            if (gene.start <= nc.start) and (
                    gene.end >= nc.end):
                detect = True
    return detect


def check_tss_ncRNA(tsss, ncs, genes):
    for tss in tsss:
        print("VVVV")
        print(tss.start)
        print(tss.attributes["type"])
        new_types = []
        if "associated_gene" in tss.attributes.keys():
            index = 0
            types = tss.attributes["type"].split(",")
            print(types)
            for ass in tss.attributes["associated_gene"].split(","):
                detect = False
                for gene in genes:
                    if "locus_tag" in gene.attributes.keys():
                        if ass == gene.attributes["locus_tag"]:
                            detect = check_ncRNA(gene, ncs)
                    elif "gene" in gene.attributes.keys():
                        if ass == gene.attributes["gene"]:
                            detect = check_ncRNA(gene, ncs)
                    elif "Name" in gene.attributes.keys():
                        if ass == gene.attributes["Name"]:
                            detect = check_ncRNA(gene, ncs)
                if detect:
                    print("NNNN")
                    if ("Primary" == types[index]) or (
                            "Secondary" == types[index]):
                        print("BBB")
                        new_types.append("Orphan")
                else:
                    new_types.append(types[index])
                index += 1
            tss.attributes["type"] = ",".join(new_types)
        print(tss.attributes["type"])


def intergenic_srna(args_srna, libs, texs, wigs_f, wigs_r):
    '''get intergenic and antisense sRNA'''
    inter_cutoff_coverage, inter_notex = get_intergenic_antisense_cutoff(
                                         args_srna)
    anti_cutoff_coverage, anti_notex = get_intergenic_antisense_cutoff(
                                       args_srna)
    testa = pandas.DataFrame.from_dict(wigs_f)
    nums, cdss, tas, pros, genes, ncs = read_data(args_srna)
    if not args_srna.tss_source:
        print("Classifying TSS...")
        compute_tss_type(args_srna, cdss, genes, wigs_f, wigs_r)
        print("Classification of TSS has done...")
    tsss, num_tss = read_tss(args_srna.tss_file)
#    check_tss_ncRNA(tsss, ncs, genes)
    detects = {"overlap": False, "uni_with_tss": False, "anti": False}
    output = open(args_srna.output_file, "w")
    out_table = open(args_srna.output_table, "w")
    output.write("##gff-version 3\n")
    for ta in tas:
        detects["overlap"] = False
        detects["anti"] = False
        compare_ta_cds(cdss, ta, detects)
        if (detects["overlap"]) and (not args_srna.in_cds):
            continue
        else:
            if not detects["anti"]:
                cutoff_coverage = inter_cutoff_coverage
                notex = inter_notex
            else:
                cutoff_coverage = anti_cutoff_coverage
                notex = anti_notex
            args_srna = ArgsContainer().extend_inter_container(
                            args_srna, tsss, pros, nums,
                            output, out_table, texs, detects,
                            cutoff_coverage, notex)
            check_srna_condition(ta, args_srna, cdss, wigs_f, wigs_r)
    file_name = args_srna.output_file.split(".")
    file_name = file_name[0] + ".stat"
    output.close()
    out_table.close()
    paras = [tsss, tas, pros, genes, cdss]
    free_memory(paras)
