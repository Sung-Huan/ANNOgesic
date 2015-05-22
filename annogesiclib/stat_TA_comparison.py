#!/usr/bin/python

import os        
import sys
import csv
from annogesiclib.gff3 import Gff3Parser

def assign_tss(tss, tran):
    if "ID" in tran.attributes.keys():
        tran_id = tran.attributes["ID"]
    else:
        tran_id = "".join([tran.feature, ":", str(tran.start), "-",
                           str(tran.end), "_", tran.strand])
    if "Parent_tran" not in tss.attributes.keys():
        tss.attributes["Parent_tran"] = tran_id
    else:
        tss.attributes["Parent_tran"] = "&".join([tss.attributes["Parent_tran"], tran_id])
    if "ID" in tss.attributes.keys():
        tss_name = tss.attributes["Name"]
    else:
        tss_name = "".join(["TSS:", str(tss.start), "_", tss.strand])
    if "associated_tss" not in tran.attributes.keys():
        tran.attributes["associated_tss"] = tss_name
    else:
        tran.attributes["associated_tss"] = "&".join([tran.attributes["associated_tss"], tss_name])

def del_attributes(entry, features):
    attributes = {}
    for key, value in entry.attributes.items():
        if (key not in features):
            attributes[key] = value
    return attributes

def compare_tran_tss(trans, tsss, fuzzy, stat, out):
    num_tran = 0
    for tran in trans:
        tran.attributes["ID"] = "tran" + str(num_tran)
        detect = False
        tran_type = ""
        check = [0, 0, 0]
        for tss in tsss:
            if (tss.strand == tran.strand) and \
               (tss.seq_id == tran.seq_id):
                if tss.strand == "+":
                    if (tss.start + int(fuzzy) >= tran.start) and \
                       (tss.start <= tran.end):
                        if check[0] != 1:
                            stat["with_TSS"] += 1
                        assign_tss(tss, tran)
                        check[0] = 1
                        detect = True
                        tss.attributes["detect"] = True
                else:
                    if (tss.end - int(fuzzy) <= tran.end) and \
                       (tss.end >= tran.start):
                        if check[0] != 1:
                            stat["with_TSS"] += 1
                        assign_tss(tss, tran)
                        check[0] = 1
                        detect = True
                        tss.attributes["detect"] = True
        if not detect:
            stat["no_TSS"] += 1
            tran.attributes["associated_tss"] = "NA"
            check[1] = 1
        else:
            detect = False
        tran.attributes = del_attributes(tran, ["TSS_note", "detect"])
        tran.attribute_string = ";".join(
                ["=".join(items) for items in tran.attributes.items()])
        if out is not None:
            out.write("\t".join([str(field) for field in [
                       tran.seq_id, tran.source, tran.feature, tran.start,
                       tran.end, tran.score, tran.strand, tran.phase,
                       tran.attribute_string + "\n"]]))
        num_tran += 1

def detect_tas_region(tsss, trans, out, out_tss, fuzzy):
    stat = {"with_TSS": 0, "no_TSS": 0, "TSS_no_tran": 0, "TSS_with_tran": 0}
    compare_tran_tss(trans, tsss, fuzzy, stat, out)
    for tss in tsss:
        if "Parent_tran" not in tss.attributes.keys():
            tss.attributes["Parent_tran"] = "NA"
        if ("detect" in tss.attributes.keys()):
            tss.attributes = del_attributes(tss, ["detect"])
            tss.attribute_string = ";".join(
                ["=".join(items) for items in tss.attributes.items()])
            stat["TSS_with_tran"] += 1
            if out_tss is not None:
                out_tss.write("\t".join([str(field) for field in [
                          tss.seq_id, tss.source, tss.feature, tss.start,
                          tss.end, tss.score, tss.strand, tss.phase,
                          tss.attribute_string]]) + "\n")
        else:
            stat["TSS_no_tran"] += 1
            if out_tss is not None:
                out_tss.write("\t".join([str(field) for field in [
                          tss.seq_id, tss.source, tss.feature, tss.start,
                          tss.end, tss.score, tss.strand, tss.phase,
                          tss.attribute_string]]) + "\n")
    return stat

def print_tas_stat(stat, out):
    total_tran = stat["with_TSS"] + stat["no_TSS"]
    total_TSS = stat["TSS_no_tran"] + stat["TSS_with_tran"]
    out.write("\tTranscript starts or overlap with TSS:{0} ({1})\n".format(
              stat["with_TSS"], float(stat["with_TSS"]) / float(total_tran)))
    out.write("\tTranscript has no relationship with TSS:{0} ({1})\n".format(
              stat["no_TSS"], float(stat["no_TSS"]) / float(total_tran)))
    out.write("\tTSS starts or overlap with transcript:{0} ({1})\n".format(
              stat["TSS_with_tran"], float(stat["TSS_with_tran"]) / float(total_TSS)))
    out.write("\tTSS has no relationship with transcript:{0} ({1})\n".format(
              stat["TSS_no_tran"], float(stat["TSS_no_tran"]) / float(total_TSS)))

def read_tas_file(tss_file, ta_file, tsss_uni, tsss, tas_uni, tas):
    tss_f = open(tss_file, "r")
    ta_f = open(ta_file, "r")
    pre_seq_id = ""
    for entry in Gff3Parser().entries(tss_f):
        entry.attributes = del_attributes(entry, ["Parent_tran", "tran_note"])
        if pre_seq_id != entry.seq_id:
            pre_seq_id = entry.seq_id
            tsss_uni[entry.seq_id] = []
        tsss_uni[entry.seq_id].append(entry)
        tsss.append(entry)
    tss_f.close()
    pre_seq_id = ""
    for entry in Gff3Parser().entries(ta_f):
        entry.attributes = del_attributes(entry, ["associated_tss"])
        if pre_seq_id != entry.seq_id:
            pre_seq_id = entry.seq_id
            tas_uni[entry.seq_id] = []
        tas_uni[entry.seq_id].append(entry)
        tas.append(entry)
    ta_f.close()

def stat_ta_tss(ta_file, tss_file, stat_file, out_ta_file, out_tss_file, fuzzy):
    tsss_uni = {}
    tas_uni = {}
    tsss = []
    tas = []
    read_tas_file(tss_file, ta_file, tsss_uni, tsss, tas_uni, tas)
    out_stat = open(stat_file, "w")
    out_stat.write("All strains:\n")
    out_ta = open(out_ta_file, "w")
    out_tss = open(out_tss_file, "w")
    out_ta.write("##gff-version 3\n")
    out_tss.write("##gff-version 3\n")
    sort_tas = sorted(tas, key=lambda k: (k.seq_id, k.start))
    stats = detect_tas_region(tsss, sort_tas, out_ta, out_tss, fuzzy)
    print_tas_stat(stats, out_stat)
    if (len(tsss_uni) > 1) and (len(tas_uni) > 1):
        for strain_tss in tsss_uni.keys():
            for strain_ta in tas_uni.keys():
                if strain_tss == strain_ta:
                    out_stat.write(strain_tss + ":\n")
                    sort_tas = sorted(tas_uni[strain_ta], key=lambda k: (k.seq_id, k.start))
                    stats = detect_tas_region(tsss_uni[strain_tss], sort_tas, None, None)
                    print_tas_stat(stats, out_stat)

def assign_parent(gff, tran):
    if "Parent_tran" not in gff.attributes.keys():
        gff.attributes["Parent_tran"] = tran.attributes["ID"]
    else:
        gff.attributes["Parent_tran"] = "&".join([gff.attributes["Parent_tran"], tran.attributes["ID"]])
    if "associated_cds" not in tran.attributes.keys():
        if "protein_id" in gff.attributes.keys():
            tran.attributes["associated_cds"] = gff.attributes["protein_id"]
        elif "locus_tag" in gff.attributes.keys():
            tran.attributes["associated_cds"] = gff.attributes["locus_tag"]
        elif "Name" in gff.attributes.keys():
            tran.attributes["associated_cds"] = gff.attributes["Name"]
        else:
            tran.attributes["associated_cds"] = "".join([gff.feature, ":", str(gff.start), 
                                                "-", str(gff.end), "_", gff.strand])
    else:
        if "protein_id" in gff.attributes.keys():
            tran.attributes["associated_cds"] = "&".join([tran.attributes["associated_cds"], 
                                                          gff.attributes["protein_id"]])
        elif "locus_tag" in gff.attributes.keys():
            tran.attributes["associated_cds"] = "&".join([tran.attributes["associated_cds"],
                                                          gff.attributes["locus_tag"]])
        elif "Name" in gff.attributes.keys():
            tran.attributes["associated_cds"] = "&".join([tran.attributes["associated_cds"],
                                                          gff.attributes["Name"]])
        else:
            tran.attributes["associated_cds"] = "&".join([tran.attributes["associated_cds"],
                                                "".join([gff.feature, ":", str(gff.start), 
                                                         "-", str(gff.end), "_", gff.strand])])

def compare_ta_gff(gffs, tran, check, tran_type, detect, stats):
    for gff in gffs:
        if (gff.feature == "CDS") or \
           (gff.feature == "tRNA") or \
           (gff.feature == "rRNA"):
            if (gff.strand == tran.strand) and \
               (gff.seq_id == tran.seq_id):
                if (gff.start < tran.start) and \
                   (gff.end > tran.end):
                    if check[0] != 1:
                        stats[tran.seq_id]["bsae"] += 1
                        stats["All"]["bsae"] += 1
                    tran_type.append("within_CDS")
                    assign_parent(gff, tran)
                    detect = True
                    check[0] = 1
                elif (gff.start >= tran.start) and \
                     (gff.end > tran.end) and \
                     (gff.start < tran.end):
                    if check[1] != 1:
                        stats[tran.seq_id]["asae"] += 1
                        stats["All"]["asae"] += 1
                    tran_type.append("left_shift_CDS")
                    assign_parent(gff, tran)
                    check[1] = 1
                    detect = True
                elif (gff.start < tran.start) and \
                     (gff.end <= tran.end) and \
                     (gff.end > tran.start):
                    if check[2] != 1:
                        stats[tran.seq_id]["bsbe"] += 1
                        stats["All"]["bsbe"] += 1
                    tran_type.append("right_shift_CDS")
                    assign_parent(gff, tran)
                    check[2] = 1
                    detect = True
                elif (gff.start >= tran.start) and \
                     (gff.end <= tran.end):
                    if check[3] != 1:
                        stats[tran.seq_id]["asbe"] += 1
                        stats["All"]["asbe"] += 1
                    tran_type.append("cover_CDS")
                    assign_parent(gff, tran)
                    check[3] = 1
                    detect = True
    return detect

def detect_tag_region(gffs, trans, stats, out_t, out_g):
    detect = False
    for tran in trans:
        check = [0, 0, 0, 0, 0]
        tran_type = []
        tran_type_string = ""
        detect = compare_ta_gff(gffs, tran, check, tran_type, detect, stats)
        if not detect:
            stats[tran.seq_id]["other"] += 1
            stats["All"]["other"] += 1
            check[4] = 1
            tran_type.append("not_related_CDS")
        else:
            detect = False
        tran_type_string = "&".join(tran_type)
        attribute_string = ";".join(
            ["=".join(items) for items in tran.attributes.items()])
        out_t.write("\t".join([tran.info_without_attributes, attribute_string]) + \
                        ";type=" + tran_type_string + "\n")
        pre_check = check
    for gff in gffs:
        attribute_string = ";".join(
            ["=".join(items) for items in gff.attributes.items()])
        out_g.write(gff.info_without_attributes + "\t" + attribute_string + "\n")
    return stats

def print_tag_stat(stats, out):
    total = stats["bsae"] + stats["bsbe"] + stats["asae"] + \
            stats["asbe"] + stats["other"]
    out.write("\tTranscript starts before and ends after CDS/rRNA/tRNA:%s (%s)\n" %
              (str(stats["asbe"]), str(float(stats["asbe"]) / float(total))))
    out.write("\tTranscript starts after and ends before CDS/rRNA/tRNA:%s (%s)\n" %
              (str(stats["bsae"]), str(float(stats["bsae"]) / float(total))))
    out.write("\tTranscript starts before and ends before CDS/rRNA/tRNA:%s (%s)\n" %
              (str(stats["asae"]), str(float(stats["asae"]) / float(total))))
    out.write("\tTranscript starts after and ends after CDS/rRNA/tRNA:%s (%s)\n" %
              (str(stats["bsbe"]), str(float(stats["bsbe"]) / float(total))))
    out.write("\tTranscript has no overlap of CDS/rRNA/tRNA:%s (%s)\n" %
              (str(stats["other"]), str(float(stats["other"]) / float(total))))
    out.write("\tTotal CDSs which have gene expression:%s (%s)\n" %
              (str(total), str(float(total) / float(stats["cds"]))))

def read_tag_file(gff_file, ta_file, gffs, tas, stats):
    pre_seq_id = ""
    ta_f = open(ta_file, "r")
    for entry in Gff3Parser().entries(ta_f):
        if entry.seq_id != pre_seq_id:
            pre_seq_id = entry.seq_id
            stats[entry.seq_id] = {"bsae": 0, "bsbe": 0, "asae": 0, "asbe": 0, "other": 0, "cds": 0}
        entry.attributes = del_attributes(entry, ["associated_cds", "type"])
        tas.append(entry)
    ta_f.close()
    g_f = open(gff_file, "r")
    for entry in Gff3Parser().entries(g_f):
        if (entry.feature == "CDS") or \
           (entry.feature == "rRNA") or \
           (entry.feature == "tRNA"):
            entry.attributes = del_attributes(entry, ["Parent_tran"])
            if entry.seq_id in stats.keys():
                stats[entry.seq_id]["cds"] += 1
                stats["All"]["cds"] += 1
        gffs.append(entry)
    g_f.close()

def stat_ta_gff(ta_file, gff_file, stat_file, out_ta_file, out_gff_file):
    gffs = []
    tas = []
    tsss = []
    stats = {}
    stats["All"] = {"bsae": 0, "bsbe": 0, "asae": 0, "asbe": 0, "other": 0, "cds": 0}
    o_f = open(out_ta_file, "w")
    og_f = open(out_gff_file, "w")
    o_f.write("##gff-version 3\n")
    og_f.write("##gff-version 3\n")
    read_tag_file(gff_file, ta_file, gffs, tas, stats)
    sort_tas = sorted(tas, key=lambda k: (k.seq_id, k.start))
    stats = detect_tag_region(gffs, sort_tas, stats, o_f, og_f)
    out_stat = open(stat_file, "w")
    out_stat.write("All strains:\n")
    out_stat.write("The transcriptome assembly information compares with annotation gff file:\n")
    print_tag_stat(stats["All"], out_stat)
    if len(stats) > 2:
        for strain in stats.keys():
            if strain != "All":
                out_stat.write(strain + ":\n")
                out_stat.write("The transcriptome assembly information compares with annotation gff file:\n")
                print_tag_stat(stats[strain], out_stat)
