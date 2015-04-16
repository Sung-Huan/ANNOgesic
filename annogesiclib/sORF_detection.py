#!/usr/bin/python

import os	
import sys
import csv
import math
import annogesiclib.parser_wig as wig_par
from annogesiclib.helper import Helper
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.lib_reader import read_libs, read_wig
from annogesiclib.coverage_detection import replicate_comparison

def get_coverage(sorf, wigs, strand, cutoff_cover):
    high_cover = -1
    low_cover = -1
    sorf_covers = {}
    for wig_strain, conds in wigs.items():
        if wig_strain == sorf["strain"]:
            for cond, tracks in conds.items():
                sorf_covers[cond] = []
                for track, covers in tracks.items():
                    total_cover = 0
                    first = True
                    if strand == "+":
                        covers = covers[sorf["start"] - 2: sorf["end"] + 1]
                    elif strand == "-":
                        covers = reversed(covers[sorf["start"] - 2: sorf["end"] + 1])
                    for cover in covers:
                        if (cover["strand"] == strand):
                            if (sorf["start"] <= cover["pos"]) and \
                               (sorf["end"] >= cover["pos"]):
                                total_cover = total_cover + cover["coverage"]
                                if first:
                                    first = False
                                    high_cover = cover["coverage"]
                                    low_cover = cover["coverage"]
                                else:
                                    if high_cover < cover["coverage"]:
                                        high_cover = cover["coverage"]
                                    if low_cover > cover["coverage"]:
                                        low_cover = cover["coverage"]
                    avg = total_cover / float(sorf["end"] - sorf["start"] + 1)
                    if avg > float(cutoff_cover):
                        sorf_covers[cond].append({"track": track, "high": high_cover,
                                                  "low": low_cover, "avg": avg,
                                                  "type": cover["type"]})
    return sorf_covers

def import_sorf(inter, sorfs, start, end, type_, fasta):
    sorfs.append({"strain": inter.seq_id,
                  "strand": inter.strand,
                  "start": start,
                  "end": end,
                  "seq": fasta[start:end],
                  "type": type_,
                  "print": False})

def detect_start_stop(inters, seq):
    sorfs = []
    for inter in inters:
        fasta = extract_seq.extract_gene(seq[inter.seq_id], \
                inter.start, inter.end, inter.strand)
        starts = []
        stops = []
        for frame in range(0, 3):
            for index in range(frame, len(fasta), 3):
                if fasta[index:index + 3] in args.start_codon:
                    starts.append(index)
                elif fasta[index:index + 3] in args.stop_codon:
                    stops.append(index)
        for start in starts:
            for stop in stops:
                if ((stop - start) > 0) and \
                   (((stop - start) % 3) == 0) and \
                   ((stop - start) <= args.Max_len) and \
                   ((stop - start) >= args.min_len):
                    if inter.source == "intergenic":
                        if inter.strand == "+":
                            import_sorf(inter, sorfs, inter.start + start, \
                                        inter.start + stop + 2, inter.source, fasta)
                        else:
                            import_sorf(inter, sorfs, inter.end - stop - 2, \
                                        inter.end - start, inter.source, fasta)
                    elif inter.source == "UTR_derived":
                        if inter.strand == "+":
                            import_sorf(inter, sorfs, inter.start + start + 1, \
                                        inter.start + stop - 1, inter.attributes["UTR_type"], \
                                        fasta)
                        else:
                            import_sorf(inter, sorfs, inter.end - stop - 2, \
                                        inter.end - start, inter.attributes["UTR_type"], \
                                        fasta)
    return sorfs

def read_data(seq, inter_gff, TSS, sRNA_gff, fasta):
    gff_parser = gff3.Gff3Parser()
    for entry in gff_parser.entries(open(inter_gff)):
        if ((entry.source == "UTR_derived") and \
           (args.utr_detect)) or \
           (entry.source == "intergenic"):
            inters.append(entry)
    inters = sorted(inters, key=lambda k: (k.seq_id, k.start))
    if TSS is not False:
        for entry in gff_parser.entries(open(TSS)):
            tsss.append(entry)
        tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start))
    else:
        tsss = None
    if sRNA_gff is not False:
        for entry in gff_parser.entries(open(sRNA_gff)):
            srnas.append(entry)
        srnas = sorted(srnas, key=lambda k: (k.seq_id, k.start))
    else:
        srnas = None
    with open(fasta, "r") as s_f:
        for line in s_f:
            line = line.strip()
            if line.startswith(">"):
                strain = line[1:]
                seq[strain] = ""
            else:
                seq[strain] = seq[strain] + line
    return (inters, tsss, srnas)

def compare_sorf_tss(sorfs, tsss):
    if args.TSS is not False:
        for sorf in sorfs:
            detect_start = False
            sorf["with_TSS"] = []
            for tss in tsss:
                check_import = False
                if (sorf["strain"] == tss.seq_id) and \
                   (sorf["strand"] == tss.strand):
                    if sorf["strand"] == "+":
                        if (sorf["start"] - tss.start <= args.utr_fuzzy) and \
                           (sorf["start"] - tss.start >= 0):
                            sorf["start_TSS"] = "TSS_" + str(tss.start) + tss.strand
                            sorf["with_TSS"].append("TSS_" + str(tss.start) + tss.strand)
                            detect_start = True
                            check_import = True
                    else:
                        if (sorf["start"] - tss.end <= args.utr_fuzzy) and \
                           (sorf["start"] - tss.end >= 0):
                            sorf["start_TSS"] = "TSS_" + str(tss.start) + tss.strand
                            sorf["with_TSS"].append("TSS_" + str(tss.start) + tss.strand)
                            detect_start = True
                            check_import = True
                    if check_import is False:
                        if (tss.start <= sorf["start"]) and \
                           (tss.start >= sorf["end"]):
                            sorf["with_TSS"].append("TSS_" + str(tss.start) + tss.strand)
            if detect_start is False:
                sorf["start_TSS"] = "NA"
            if len(sorf["with_TSS"]) == 0:
                sorf["with_TSS"] = ["NA"]
    else:
        for sorf in sorfs:
            sorf["with_TSS"] = ["NA"]
            sorf["start_TSS"] = "NA"

def compare_sorf_sRNA(sorfs, srnas):
    if args.sRNA_gff is not False:
        for sorf in sorfs:
            sorf["srna"] = []
            for srna in srnas:
                if (sorf["strain"] == srna.seq_id) and \
                   (sorf["strand"] == srna.strand):
                    if (srna.start <= sorf["start"]) and \
                       (srna.end >= sorf["end"]):
                        sorf["srna"].append(srna.attributes["ID"] + ":" + str(srna.start) + \
                                        "-" + str(srna.end) + "_" + srna.strand)
            if len(sorf["srna"]) == 0:
                sorf["srna"] = ["NA"]
    else:
        sorf["srna"] = ["NA"]

def merge(sorfs, seq):
    finals = []
    for sorf1 in sorfs:
        final = sorf1.copy()
        first = True
        if sorf1["print"] is not True:
            sorf1["print"] = True
            for sorf2 in sorfs:
                overlap = False
                if (final["strain"] == sorf2["strain"]) and \
                   (final["strand"] == sorf2["strand"]):
                    if (final["start"] >= sorf2["start"]) and \
                       (final["end"] <= sorf2["end"]):
                        overlap = True
                    elif (final["start"] >= sorf2["start"]) and \
                         (final["start"] <= sorf2["end"]) and \
                         (final["end"] >= sorf2["end"]):
                        overlap = True
                    elif (final["start"] <= sorf2["start"]) and \
                         (final["end"] >= sorf2["start"]) and \
                         (final["end"] <= sorf2["end"]):
                        overlap = True
                    elif (final["start"] <= sorf2["start"]) and \
                         (final["end"] >= sorf2["end"]):
                        overlap = True
                    elif (sorf2["start"] > final["end"]):
                        break
                    if overlap:
                        sorf2["print"] = True
                        if final["start"] > sorf2["start"]:
                            final["start"] = sorf2["start"]
                        if final["end"] < sorf2["end"]:
                            final["end"] = sorf2["end"]
                        if first:
                            final["candidate"] = []
                            final["candidate"].append("-".join([str(sorf1["start"]), str(sorf1["end"])]) \
                                                      + "_TSS:" + sorf1["start_TSS"])
                            first = False
                        if "-".join([str(sorf2["start"]), str(sorf2["end"]) + \
                                     "_TSS:" + sorf2["start_TSS"]]) not in final["candidate"]:
                            final["candidate"].append("-".join([str(sorf2["start"]), str(sorf2["end"])]) \
                                                      + "_TSS:" + sorf1["start_TSS"])
            final["seq"] = self.helper.extract_gene(seq[final["strain"]], \
                           final["start"], final["end"], final["strand"])
            del final["print"]
            finals.append(final)
    return finals

def get_cutoff(sorf, medians, coverages):
    if sorf["type"] == "intergenic":
        cutoff_cover = float(coverages["inter"])
    elif ("5utr" in sorf["type"]) and ("3utr" in sorf["type"]):
        if coverages["3utr"] == "median":
            cutoff_utr3 = medians[sorf["strain"]]["3utr"]
        else:
            cutoff_utr3 = float(coverages["5utr"])
        if coverages["5utr"] == "median":
            cutoff_utr5 = medians[sorf["strain"]]["5utr"]
        else:
            cutoff_utr5 = float(coverages["5utr"])
        cutoff_cover = min(cutoff_utr5, cutoff_utr3)
    elif ("5utr" in sorf["type"]):
        if coverages["5utr"] == "median":
            cutoff_cover = medians[sorf["strain"]]["5utr"]
        else:
            cutoff_cover = float(coverages["5utr"])
    elif ("3utr" in sorf["type"]):
        if coverages["3utr"] == "median":
            cutoff_cover = medians[sorf["strain"]]["3utr"]
        else:
            cutoff_cover = float(coverages["3utr"])
    elif ("interCDS" in sorf["type"]):
        if coverages["interCDS"] == "median":
            cutoff_cover = medians[sorf["strain"]]["interCDS"]
        else:
            cutoff_cover = float(coverages["interCDS"])
    return cutoff_cover

def print_file(sorf, sorf_datas, num, out_g, out_t):
    name = '%0*d' % (5, num)
    if sorf["type"] == "intergenic":
        source = "intergenic"
        type_ = "Intergenic"
        attribute_string = ";".join(
            ["=".join(items) for items in (["ID", "sorf" + str(num)],
                                           ["Name", "sORF_" + name],
                                           ["start_TSS", sorf["start_TSS"]],
                                           ["with_TSS", "&".join(sorf["with_TSS"])],
                                           ["candidates", "&".join(sorf["candidate"])],
                                           ["sRNA", "&".join(sorf["srna"])])])
    else:
        source = "UTR_derived"
        if ("3utr" in sorf["type"]) and ("5utr" in sorf["type"]):
            type_ = "3'UTR_derived;5'UTR_derived"
        elif ("3utr" in sorf["type"]):
            type_ = "3'UTR_derived"
        elif ("5utr" in sorf["type"]):
            type_ = "5'UTR_derived"
        elif ("interCDS" in sorf["type"]):
            type_ = "interCDS"
        attribute_string = ";".join(
            ["=".join(items) for items in (["ID", "sorf" + str(num)],
                                           ["Name", "sORF_" + name],
                                           ["start_TSS", sorf["start_TSS"]],
                                           ["with_TSS", "&".join(sorf["with_TSS"])],
                                           ["candidates", "&".join(sorf["candidate"])],
                                           ["UTR_type", sorf["type"]],
                                           ["sRNA", "&".join(sorf["srna"])])])
    info = "\t".join([str(field) for field in [
                      sorf["strain"], source, "sORF", str(sorf["start"]),
                      str(sorf["end"]), ".", sorf["strand"], ".", attribute_string]])
    out_g.write(info + "\n")
    check_frag = False
    check_tex = False
    if ("frag" in ";".join(sorf_datas["cond"].keys())) and \
       ("tex" in ";".join(sorf_datas["cond"].keys())):
        lib_type = "TEX+/-;Fragmented"
    elif ("frag" in ";".join(sorf_datas["cond"].keys())):
        lib_type = "Fragmented"
    elif ("tex" in ";".join(sorf_datas["cond"].keys())):
        lib_type = "TEX+/-"
    print_table(out_t, sorf, name, type_, lib_type, sorf_datas, table_best)

def print_table():
    out_t.write("\t".join([sorf["strain"], "sORF_" + name, str(sorf["start"]),
                           str(sorf["end"]), sorf["strand"], type_,
                           ";".join(sorf["with_TSS"]), ";".join(sorf["candidate"]),
                           ";".join(sorf["srna"]), lib_type, str(sorf_datas["best_cover"]),
                           str(sorf_datas["high"]), str(sorf_datas["low"]),
                           ";".join(sorf["srna"]), str(sorf_datas["best_cover"]),
                           str(sorf_datas["high"]), str(sorf_datas["low"])]) + "\t")
    if table_best is False:
        first = True
        for data in sorf_datas["detail"]:
            if first:
                out_t.write("%s(avg=%s;high=%s;low=%s)" % (data["track"], data["avg"], data["high"], data["low"]))
                first = False
            else:
                out_t.write(";%s(avg=%s;high=%s;low=%s)" % (data["track"], data["avg"], data["high"], data["low"]))
    else:
        out_t.write("%s(avg=%s;high=%s;low=%s)" % (sorf_datas["best_track"], sorf_datas["best_cover"], sorf_datas["high"], sorf_datas["low"]))
    out_t.write(sorf["seq"])
    out_t.write("\n")

def get_inter_coverage(inters, inter_covers, utr_type, template_texs):
    first = True
    for inter in inters:
        best_cover = -1
        total = 0
        detect_num = 0
        check_texs = {}
        texs = template_texs.copy()
        for key, num in texs.items():
            check_texs[key] = []
        for cond, covers in inter.items():
            for cover in covers:
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
            inter_covers.append(best_cover)

def detect_utr_type(inter, utr_type, med_inters, wigs, strand):
    if inter.attributes["UTR_type"] == utr_type:
        inter_datas = {}
        inter_datas["strain"] = inter.seq_id
        inter_datas["strand"] = inter.strand
        inter_datas["start"] = inter.start
        inter_datas["end"] = inter.end
        med_inters[inter.seq_id][utr_type].append(get_coverage(inter_datas, wigs, strand, 0))

def median(lst):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = (lstLen - 1) // 2
    if (lstLen % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1])/2.0

def validate_tss(starts, ends, sorf):
    tsss = []
    start_pos = "NA"
    for tss in sorf["with_TSS"]:
        tss_start = int(tss.replace("TSS_", "")[:-1])
        if sorf["strand"] == "+":
            if (tss_start >= min(starts) - args.utr_fuzzy) and \
               (tss_start <= max(ends)):
                tsss.append(tss)
                if (tss_start >= min(starts) - args.utr_fuzzy) and \
                   (tss_start <= min(starts)):
                    start_pos = tss
        else:
            if (tss_start >= min(starts)) and \
               (tss_start <= max(ends) + args.utr_fuzzy):
                tsss.append(tss)
                if (tss_start <= min(ends) + args.utr_fuzzy) and \
                   (tss_start >= min(ends)):
                    start_pos = tss
                    break
    return (tsss, start_pos)

def validate_srna(starts, ends, sorf):
    srnas = []
    for srna in sorf["srna"]:
        if srna == "NA":
            break
        else:
            datas = srna.split(":")[1][:-2].split("-")
            start = int(datas[0])
            end = int(datas[1])
            for index in range(0, len(starts)):
                if (start <= starts[index]) and \
                   (end >= ends[index]):
                    srnas.append(srna)
                    break
    if len(srnas) == 0:
        srnas = ["NA"]
    return srnas

def get_best(sorfs, condition):
    final_sorfs = []
    for sorf in sorfs:
        if (condition == "TSS") or \
           (condition == "both"):
            if sorf["with_TSS"][0] != "NA":
                starts = []
                ends = []
                cands = []
                tmp_sorf = sorf.copy()
                for candidate in sorf["candidate"]:
                    tss = candidate.split("_TSS:")[1]
                    if tss != "NA":
                        datas = candidate.split("_TSS:")[0].split("-")
                        starts.append(int(datas[0]))
                        ends.append(int(datas[1]))
                        cands.append("-".join([str(datas[0]), str(datas[1])]) \
                                               + "_TSS:" + tss)
                tmp_sorf["start"] = min(starts)
                tmp_sorf["end"] = max(ends)
                tmp_sorf["candidate"] = cands
                tsss_datas = validate_tss(starts, ends, sorf)
                tmp_sorf["with_TSS"] = tsss_datas[0]
                tmp_sorf["start_TSS"] = tsss_datas[1]
                tmp_sorf["sRNA"] = validate_srna(starts, ends, sorf)
                if args.condition == "TSS":
                    final_sorfs.append(tmp_sorf)
                elif args.condition == "both":
                    if tmp_sorf["sRNA"][0] == "NA":
                        final_sorfs.append(tmp_sorf)
        elif (condition == "sRNA"):
            if sorf["sRNA"][0] == "NA":
                final_sorfs.append(sorf)
    return final_sorfs

def coverage_and_output(sorfs, mediandict, wigs_f, wigs_r, texs, out_g, out_t, tex_notex, replicates):
    out_g.write("##gff-version 3\n")
    out_t.write("\t".join(["strain", "Name", "start", "end", "strand", "type",
                           "TSS", "candidates", "lib_type", "sRNA_confliction",
                           "best_avg_coverage", "best_highest_coverage",
                           "best_lowest_coverage", "track_detail", "seq"]) + "\n")
    num = 0
    for sorf in sorfs:
        cutoff_cover = get_cutoff(sorf, mediandict)
        if sorf["strand"] == "+":
            sorf_covers = get_coverage(sorf, wigs_f, "+", cutoff_cover)
        else:
            sorf_covers = get_coverage(sorf, wigs_r, "-", cutoff_cover)
        if len(sorf_covers) != 0:
            sorf_info = replicate_comparison(sorf_covers, texs, sorf.strand, cutoff_cover,
                                             tex_notex, replicates, "sORF", None, None, None)
            if len(sorf_info["cond"].keys()) != 0:
                print_file(sorf, sorf_info, num, out_g, out_t)
                num += 1

def detect_inter_type(inters, wigs_f, wigs_r):
    strain = ""
    for inter in inters:
        if inter.seq_id != strain:
            strain = inter.seq_id
            med_inters[inter.seq_id] = {"5utr": [], "3utr": [], "interCDS": []}
        if (inter.source == "UTR_derived") and (inter.strand == "+"):
            detect_utr_type(inter, "5utr", med_inters, wigs_f, "+")
            detect_utr_type(inter, "3utr", med_inters, wigs_f, "+")
            detect_utr_type(inter, "interCDS", med_inters, wigs_f, "+")
        elif (inter.source == "UTR_derived") and (inter.strand == "-"):
            detect_utr_type(inter, "5utr", med_inters, wigs_r, "-")
            detect_utr_type(inter, "3utr", med_inters, wigs_r, "-")
            detect_utr_type(inter, "interCDS", med_inters, wigs_r, "-")

def sorf_detection(fasta, sRNA_gff, inter_gff, TSS, utr_fuzzy, utr_detect, input_libs, tex_notex,
                   replicates, inter_coverage, utr3_coverage, utr5_coverage, interCDS_coverage,
                   wig_f_file, wig_r_file, wig_folder, start_codon, stop_codon, table_best,
                   Max_len, min_len, condition, out_prefix):
    seq = {}
    wigs_f = {}
    wigs_r = {}
    libs = []
    texs= {}
    coverages = {"3utr": utr3_coverage, "5utr": utr5_coverage,
                 "inter": inter_coverage, "interCDS": interCDS_coverage}
    read_libs(libs, texs, input_libs, wig_folder)
    datas = read_data(seq, inter_gff, TSS, sRNA_gff, fasta)
    inters = datas[0]
    tsss = datas[1]
    srnas = datas[2]
    read_wig(wigs_f, wig_f_file, "+", libs)
    read_wig(wigs_r, wig_r_file, "-", libs)
    med_inters = {}
    detect_inter_type(inters, wigs_f, wigs_r)
    inter_covers = {}
    mediandict = {}
    if utr_detect:
        for strain, meds in med_inters.items():
            inter_covers[strain] = {"5utr": [], "3utr": [], "interCDS": []}
            mediandict[strain] = {"3utr": 0, "5utr": 0, "interCDS": 0}
            for type_, covers in meds.items():
                get_inter_coverage(covers, inter_covers[strain][type_], type_, texs)
                if len(inter_covers[strain][type_]) != 0:
                    mediandict[strain][type_] = median(inter_covers[strain][type_])
    out_ag = open("_".join([out_prefix, "all.gff"]), "w")
    out_at = open("_".join([out_prefix, "all.csv"]), "w")
    out_bg = open("_".join([out_prefix, "best.gff"]), "w")
    out_bt = open("_".join([out_prefix, "best.csv"]), "w")
    sorfs = detect_start_stop(inters, seq)
    compare_sorf_tss(sorfs, tsss)
    compare_sorf_sRNA(sorfs, srnas)
    sorfs = sorted(sorfs, key=lambda k: (k["strain"], k["start"]))
    sorfs = merge(sorfs, seq)
    final_sorfs = get_best(sorfs, condition)
    coverage_and_output(sorfs, mediandict, wigs_f, wigs_r, texs, 
                        out_ag, out_at, tex_notex, replicates)
    coverage_and_output(final_sorfs, mediandict, wigs_f, wigs_r, texs, 
                        out_bg, out_bt, tex_notex, replicates)
