#!/usr/bin/python

import os
import sys
import copy
from annogesiclib.helper import Helper
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.lib_reader import read_libs, read_wig
from annogesiclib.coverage_detection import replicate_comparison

def get_coverage(sorf, wigs, strand, coverages, medianlist, cutoffs):
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
                            if (sorf["start"] <= cover["pos"]) and (
                                sorf["end"] >= cover["pos"]):
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
                    if medianlist is not None:
                        cutoff_cover = get_cutoff(sorf, track,
                                       coverages, medianlist)
                    else:
                        cutoff_cover = coverages
                    if cutoffs is not None:
                        cutoffs[track] = cutoff_cover
                    if avg > float(cutoff_cover):
                        sorf_covers[cond].append({"track": track,
                                                  "high": high_cover,
                                                  "low": low_cover, "avg": avg,
                                                  "type": cover["type"],
                                                  "pos": sorf["start"]})
    return sorf_covers

def import_sorf(inter, sorfs, start, end, type_, fasta):
    sorfs.append({"strain": inter.seq_id,
                  "strand": inter.strand,
                  "start": start,
                  "end": end,
                  "starts": [str(start)],
                  "ends": [str(end)],
                  "seq": fasta[start:end],
                  "type": type_,
                  "print": False})

def detect_start_stop(inters, seq, start_codon, stop_codon, max_len, min_len):
    sorfs = []
    for inter in inters:
        fasta = Helper().extract_gene(seq[inter.seq_id],
                inter.start, inter.end, inter.strand)
        starts = []
        stops = []
        for frame in range(0, 3):
            for index in range(frame, len(fasta), 3):
                if fasta[index:index + 3] in start_codon:
                    starts.append(index)
                elif fasta[index:index + 3] in stop_codon:
                    stops.append(index)
        for start in starts:
            for stop in stops:
                if ((stop - start) > 0) and \
                   (((stop - start) % 3) == 0) and \
                   ((stop - start) <= max_len) and \
                   ((stop - start) >= min_len):
                    if inter.source == "intergenic":
                        if inter.strand == "+":
                            import_sorf(inter, sorfs, inter.start + start,
                                        inter.start + stop + 2,
                                        inter.source, fasta)
                        else:
                            import_sorf(inter, sorfs, inter.end - stop - 2,
                                        inter.end - start, inter.source, fasta)
                    elif inter.source == "UTR_derived":
                        if inter.strand == "+":
                            import_sorf(inter, sorfs, inter.start + start + 1,
                                        inter.start + stop - 1,
                                        inter.attributes["UTR_type"], fasta)
                        else:
                            import_sorf(inter, sorfs, inter.end - stop - 2,
                                        inter.end - start,
                                        inter.attributes["UTR_type"], fasta)
    return sorfs

def read_data(inter_gff, tss_file, srna_gff, fasta, utr_detect):
    seq = {}
    inters = []
    tsss = []
    srnas = []
    for entry in Gff3Parser().entries(open(inter_gff)):
        if ((entry.source == "UTR_derived") and \
           (utr_detect)) or \
           (entry.source == "intergenic"):
            inters.append(entry)
    inters = sorted(inters, key=lambda k: (k.seq_id, k.start))
    if tss_file is not None:
        for entry in Gff3Parser().entries(open(tss_file)):
            tsss.append(entry)
        tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start))
    else:
        tsss = None
    if srna_gff is not None:
        for entry in Gff3Parser().entries(open(srna_gff)):
            new = {}
            for key, value in entry.attributes.items():
                if "sORF" not in key:
                    new[key] = value
            entry.attributes = copy.deepcopy(new)
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
    return inters, tsss, srnas, seq

def compare_sorf_tss(sorfs, tsss, tss_file, utr_fuzzy):
    if tss_file is not None:
        for sorf in sorfs:
            detect_start = False
            sorf["with_TSS"] = []
            for tss in tsss:
                check_import = False
                if (sorf["strain"] == tss.seq_id) and (
                    sorf["strand"] == tss.strand):
                    if sorf["strand"] == "+":
                        if (sorf["start"] - tss.start <= utr_fuzzy) and (
                            sorf["start"] - tss.start >= 0):
                            sorf["start_TSS"] = "TSS_" + str(tss.start) + tss.strand
                            sorf["with_TSS"].append("TSS_" + str(tss.start) + tss.strand)
                            detect_start = True
                            check_import = True
                    else:
                        if (sorf["start"] - tss.end <= utr_fuzzy) and (
                            sorf["start"] - tss.end >= 0):
                            sorf["start_TSS"] = "TSS_" + str(tss.start) + tss.strand
                            sorf["with_TSS"].append("TSS_" + str(tss.start) + tss.strand)
                            detect_start = True
                            check_import = True
                    if not check_import:
                        if (tss.start <= sorf["start"]) and (
                            tss.start >= sorf["end"]):
                            sorf["with_TSS"].append("TSS_" + str(tss.start) + tss.strand)
            if not detect_start:
                sorf["start_TSS"] = "NA"
            if len(sorf["with_TSS"]) == 0:
                sorf["with_TSS"] = ["NA"]
    else:
        for sorf in sorfs:
            sorf["with_TSS"] = ["NA"]
            sorf["start_TSS"] = "NA"

def compare_sorf_srna(sorfs, srnas, srna_gff):
    if srna_gff is not None:
        for sorf in sorfs:
            sorf["srna"] = []
            for srna in srnas:
                if (sorf["strain"] == srna.seq_id) and (
                    sorf["strand"] == srna.strand):
                    if (srna.start <= sorf["start"]) and (
                        srna.end >= sorf["end"]):
                        strand = Helper().get_strand_name(srna.strand)
                        sorf["srna"].append(srna.attributes["ID"] + ":" + \
                                     str(srna.start) + "-" + str(srna.end) + \
                                     "_" + strand)
            if len(sorf["srna"]) == 0:
                sorf["srna"] = ["NA"]
    else:
        for sorf in sorfs:
            sorf["srna"] = ["NA"]

def import_overlap(sorf2, final, sorf1, first):
    sorf2["print"] = True
    if final["start"] > sorf2["start"]:
        final["start"] = sorf2["start"]
    if final["end"] < sorf2["end"]:
        final["end"] = sorf2["end"]
    if first:
        final["candidate"] = []
        final["candidate"].append("-".join([
                str(sorf1["start"]), str(sorf1["end"])]) + \
                "_TSS:" + sorf1["start_TSS"])
        first = False
#        final["starts"] = []
#        final["starts"].append(str(sorf1["start"]))
#        final["ends"] = []
#        final["ends"].append(str(sorf1["end"]))
    if "-".join([str(sorf2["start"]), str(sorf2["end"]) + \
                 "_TSS:" + sorf2["start_TSS"]]) not in final["candidate"]:
        final["candidate"].append("-".join([
                str(sorf2["start"]), str(sorf2["end"])]) + \
                "_TSS:" + sorf1["start_TSS"])
    if str(sorf2["start"]) not in final["starts"]:
        final["starts"].append(str(sorf2["start"]))
    if str(sorf2["end"]) not in final["ends"]:
        final["ends"].append(str(sorf2["end"]))
    return first

def merge(sorfs, seq):
    finals = []
    for sorf1 in sorfs:
        final = copy.deepcopy(sorf1)
        first = True
        if not sorf1["print"]:
            sorf1["print"] = True
            for sorf2 in sorfs:
                overlap = False
                if (final["strain"] == sorf2["strain"]) and (
                    final["strand"] == sorf2["strand"]):
                    if (final["start"] >= sorf2["start"]) and (
                        final["end"] <= sorf2["end"]):
                        overlap = True
                    elif (final["start"] >= sorf2["start"]) and (
                          final["start"] <= sorf2["end"]) and (
                          final["end"] >= sorf2["end"]):
                        overlap = True
                    elif (final["start"] <= sorf2["start"]) and (
                          final["end"] >= sorf2["start"]) and (
                          final["end"] <= sorf2["end"]):
                        overlap = True
                    elif (final["start"] <= sorf2["start"]) and (
                          final["end"] >= sorf2["end"]):
                        overlap = True
                    elif (sorf2["start"] > final["end"]):
                        break
                    if overlap:
                        first = import_overlap(sorf2, final, sorf1, first)
            final["seq"] = Helper().extract_gene(seq[final["strain"]],
                           final["start"], final["end"], final["strand"])
            new = {}
            for key, value in final.items():
                if "print" not in key:
                    new[key] = value
            final = copy.deepcopy(new)
            finals.append(final)
    return finals

def assign_utr_cutoff(coverages, utr_type, medians):
    if coverages[utr_type] == "median":
        cutoff = medians["median"]
    elif coverages[utr_type] == "mean":
        cutoff = medians["mean"]
    else:
        cutoff = float(coverages[utr_type])
    return cutoff

def get_cutoff(sorf, track, coverages, medians):
    if sorf["type"] == "intergenic":
        cutoff_cover = float(coverages["inter"])
    elif ("5utr" in sorf["type"]) and ("3utr" in sorf["type"]):
        cutoff_utr3 = assign_utr_cutoff(coverages, "3utr",
                          medians[sorf["strain"]]["3utr"][track])
        cutoff_utr5 = assign_utr_cutoff(coverages, "5utr",
                          medians[sorf["strain"]]["5utr"][track])
        cutoff_cover = min(cutoff_utr5, cutoff_utr3)
    elif ("5utr" in sorf["type"]):
        cutoff_cover = assign_utr_cutoff(coverages, "5utr",
                          medians[sorf["strain"]]["5utr"][track])
    elif ("3utr" in sorf["type"]):
        cutoff_cover = assign_utr_cutoff(coverages, "3utr",
                          medians[sorf["strain"]]["3utr"][track])
    elif ("interCDS" in sorf["type"]):
        cutoff_cover = assign_utr_cutoff(coverages, "interCDS",
                          medians[sorf["strain"]]["interCDS"][track])
    return cutoff_cover

def print_file(sorf, sorf_datas, num, out_g, out_t, table_best):
    name = '%0*d' % (5, num)
    if sorf["type"] == "intergenic":
        source = "intergenic"
        type_ = "Intergenic"
        attribute_string = ";".join(
            ["=".join(items) for items in (["ID", "sorf" + str(num)],
                                 ["Name", "sORF_" + name],
                                 ["start_TSS", sorf["start_TSS"]],
                                 ["with_TSS", "&".join(sorf["with_TSS"])],
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
                                 ["UTR_type", sorf["type"]],
                                 ["sRNA", "&".join(sorf["srna"])])])
    info = "\t".join([str(field) for field in [
                      sorf["strain"], source, "sORF", str(sorf["start"]),
                      str(sorf["end"]), ".", sorf["strand"],
                      ".", attribute_string]])
    out_g.write(info + "\n")
    check_frag = False
    check_tex = False
    if ("frag" in ";".join(sorf_datas["conds"].keys())) and (
        "tex" in ";".join(sorf_datas["conds"].keys())):
        lib_type = "TEX+/-;Fragmented"
    elif ("frag" in ";".join(sorf_datas["conds"].keys())):
        lib_type = "Fragmented"
    elif ("tex" in ";".join(sorf_datas["conds"].keys())):
        lib_type = "TEX+/-"
    print_table(out_t, sorf, name, type_, lib_type, sorf_datas, table_best)

def print_table(out_t, sorf, name, type_, lib_type, sorf_datas, table_best):
    out_t.write("\t".join([sorf["strain"], "sORF_" + name, str(sorf["start"]),
                           str(sorf["end"]), sorf["strand"], type_,
                           ";".join(sorf["with_TSS"]),
                           ";".join(sorf["starts"]), ";".join(sorf["ends"]),
                           ";".join(sorf["srna"]),
                           lib_type, str(sorf_datas["best"]),
                           str(sorf_datas["high"]), str(sorf_datas["low"]),
                           ";".join(sorf["srna"]), str(sorf_datas["best"]),
                           str(sorf_datas["high"]),
                           str(sorf_datas["low"])]) + "\t")
    if not table_best:
        first = True
        for data in sorf_datas["detail"]:
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
                    sorf_datas["track"], sorf_datas["best"],
                    sorf_datas["high"], sorf_datas["low"]))
    out_t.write("\t" + sorf["seq"])
    out_t.write("\n")

def get_inter_coverage(inters, inter_covers, utr_type,
                       template_texs, tex_notex):
    for datas in inters:
        for cond, covers in datas.items():
            for inter in covers:
                if inter["track"] not in inter_covers.keys():
                    inter_covers[inter["track"]] = []
                inter_covers[inter["track"]].append(inter["avg"])

def detect_utr_type(inter, utr_type, med_inters, wigs, strand, background):
    if inter.attributes["UTR_type"] == utr_type:
        inter_datas = {}
        inter_datas["strain"] = inter.seq_id
        inter_datas["strand"] = inter.strand
        inter_datas["start"] = inter.start
        inter_datas["end"] = inter.end
        inter_datas = get_coverage(inter_datas, wigs, strand, background, None, None)
        med_inters[inter.seq_id][utr_type].append(inter_datas)

def median_score(lst):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = (lstLen - 1) // 2
    if (lstLen % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1])/2.0

def mean_score(lst):
    total = 0
    for li in lst:
        total = total + li
    if len(lst) != 0:
        return (total / len(lst))
    else:
        return 0

def validate_tss(starts, ends, sorf, utr_fuzzy):
    tsss = []
    start_pos = "NA"
    for tss in sorf["with_TSS"]:
        tss_start = int(tss.replace("TSS_", "")[:-1])
        if sorf["strand"] == "+":
            if (tss_start >= min(starts) - utr_fuzzy) and (
                tss_start <= max(ends)):
                tsss.append(tss)
                if (tss_start >= min(starts) - utr_fuzzy) and (
                    tss_start <= min(starts)):
                    start_pos = tss
        else:
            if (tss_start >= min(starts)) and (
                tss_start <= max(ends) + utr_fuzzy):
                tsss.append(tss)
                if (tss_start <= min(ends) + utr_fuzzy) and (
                    tss_start >= min(ends)):
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
                if (start <= starts[index]) and (
                    end >= ends[index]):
                    srnas.append(srna)
                    break
    if len(srnas) == 0:
        srnas = ["NA"]
    return srnas

def get_best(sorfs, condition, utr_fuzzy):
    final_sorfs = []
    for sorf in sorfs:
        if (condition == "TSS") or (condition == "both"):
            if sorf["with_TSS"][0] != "NA":
                cands = []
                starts = []
                ends = []
                tmp_sorf = copy.deepcopy(sorf)
                for candidate in sorf["candidate"]:
                    tss = candidate.split("_TSS:")[1]
                    if tss != "NA":
                        datas = candidate.split("_TSS:")[0].split("-")
                        cands.append("-".join([
                              str(datas[0]), str(datas[1])]) + "_TSS:" + tss)
                        starts.append(int(datas[0]))
                        ends.append(int(datas[1]))
                tmp_sorf["start"] = min(starts)
                tmp_sorf["end"] = max(ends)
                tmp_sorf["starts"] = sorf["starts"]
                tmp_sorf["ends"] = sorf["ends"]
                tmp_sorf["candidate"] = cands
                tsss_datas = validate_tss(starts, ends, sorf, utr_fuzzy)
                tmp_sorf["with_TSS"] = tsss_datas[0]
                tmp_sorf["start_TSS"] = tsss_datas[1]
                tmp_sorf["sRNA"] = validate_srna(starts, ends, sorf)
                if condition == "TSS":
                    final_sorfs.append(tmp_sorf)
                elif condition == "both":
                    if tmp_sorf["sRNA"][0] == "NA":
                        final_sorfs.append(tmp_sorf)
        elif (condition == "sRNA"):
            if sorf["sRNA"][0] == "NA":
                final_sorfs.append(sorf)
    return final_sorfs

def coverage_and_output(sorfs, mediandict, wigs_f, wigs_r, texs, out_g, out_t,
                        tex_notex, replicates, coverages, table_best):
    out_g.write("##gff-version 3\n")
    out_t.write("\t".join(["strain", "Name", "start", "end", "strand", "type",
                "TSS", "start_codons", "stop_codons", "lib_type",
                "sRNA_confliction", "best_avg_coverage",
                "best_highest_coverage", "best_lowest_coverage",
                "track_detail", "seq"]) + "\n")
    num = 0
    for sorf in sorfs:
        cutoffs = {}
        if sorf["strand"] == "+":
            sorf_covers = get_coverage(sorf, wigs_f, "+", coverages, mediandict, cutoffs)
        else:
            sorf_covers = get_coverage(sorf, wigs_r, "-", coverages, mediandict, cutoffs)
        if len(sorf_covers) != 0:
            sorf_info = replicate_comparison(sorf_covers, texs,
                                sorf["strand"], None, tex_notex, replicates,
                                "sORF", None, cutoffs, None, False)
            if len(sorf_info["conds"].keys()) != 0:
                print_file(sorf, sorf_info, num, out_g, out_t, table_best)
                num += 1

def detect_inter_type(inters, wigs_f, wigs_r, background):
    med_inters = {}
    strain = ""
    for inter in inters:
        if inter.seq_id != strain:
            strain = inter.seq_id
            med_inters[inter.seq_id] = {"5utr": [], "3utr": [], "interCDS": []}
        if (inter.source == "UTR_derived") and (inter.strand == "+"):
            detect_utr_type(inter, "5utr", med_inters, wigs_f, "+", background)
            detect_utr_type(inter, "3utr", med_inters, wigs_f, "+", background)
            detect_utr_type(inter, "interCDS", med_inters,
                            wigs_f, "+", background)
        elif (inter.source == "UTR_derived") and (inter.strand == "-"):
            detect_utr_type(inter, "5utr", med_inters, wigs_r, "-", background)
            detect_utr_type(inter, "3utr", med_inters, wigs_r, "-", background)
            detect_utr_type(inter, "interCDS", med_inters,
                            wigs_r, "-", background)
    return med_inters

def set_median(covers, mediandict):
    for strain, utrs in covers.items():
        mediandict[strain] = {"3utr": {}, "5utr": {}, "interCDS": {}}
        for utr, tracks in utrs.items():
            for track, avgs in tracks.items():
                if track not in mediandict[strain][utr].keys():
                    mediandict[strain][utr][track] = {}
                mediandict[strain][utr][track] = {"median": median_score(avgs),
                                                  "mean": mean_score(avgs)}

def sorf_detection(fasta, srna_gff, inter_gff, tss_file, utr_fuzzy, utr_detect,
                   input_libs, tex_notex, replicates, inter_coverage,
                   utr3_coverage, utr5_coverage, interCDS_coverage, wig_f_file,
                   wig_r_file, wig_folder, start_codon, stop_codon, table_best,
                   max_len, min_len, condition, out_prefix, background):
    coverages = {"3utr": utr3_coverage, "5utr": utr5_coverage,
                 "inter": inter_coverage, "interCDS": interCDS_coverage}
    libs, texs = read_libs(input_libs, wig_folder)
    inters, tsss, srnas, seq = read_data(inter_gff, tss_file, srna_gff,
                                         fasta, utr_detect)
    wigs_f = read_wig(wig_f_file, "+", libs)
    wigs_r = read_wig(wig_r_file, "-", libs)
    med_inters = detect_inter_type(inters, wigs_f, wigs_r, background)
    inter_covers = {}
    mediandict = {}
    if utr_detect:
        for strain, meds in med_inters.items():
            inter_covers[strain] = {"5utr": {}, "3utr": {}, "interCDS": {}}
            for type_, covers in meds.items():
                get_inter_coverage(covers, inter_covers[strain][type_],
                                   type_, texs, tex_notex)
        set_median(inter_covers, mediandict)
    out_ag = open("_".join([out_prefix, "all.gff"]), "w")
    out_at = open("_".join([out_prefix, "all.csv"]), "w")
    out_bg = open("_".join([out_prefix, "best.gff"]), "w")
    out_bt = open("_".join([out_prefix, "best.csv"]), "w")
    sorfs = detect_start_stop(inters, seq, start_codon,
                              stop_codon, max_len, min_len)
    compare_sorf_tss(sorfs, tsss, tss_file, utr_fuzzy)
    compare_sorf_srna(sorfs, srnas, srna_gff)
    sorfs = sorted(sorfs, key=lambda k: (k["strain"], k["start"]))
    sorfs = merge(sorfs, seq)
    final_sorfs = get_best(sorfs, condition, utr_fuzzy)
    coverage_and_output(sorfs, mediandict, wigs_f, wigs_r, texs, out_ag, out_at,
                        tex_notex, replicates, coverages, table_best)
    coverage_and_output(final_sorfs, mediandict, wigs_f, wigs_r, texs, out_bg,
                        out_bt, tex_notex, replicates, coverages, table_best)
