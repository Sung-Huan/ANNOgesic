#!/usr/bin/python

import os
import sys
import csv
from copy import deepcopy
from annogesiclib.gff3 import Gff3Parser


def import_data(row, type_):
    if type_ == "gff":
        return {"strain": row[0], "name": row[1], "start": int(row[2]),
                "end": int(row[3]), "strand": row[4], "conds": row[5],
                "detect": row[6], "tss_pro": row[7], "avg": float(row[8]),
                "high": float(row[9]), "low": float(row[10]), "track": row[11]}
    elif type_ == "nr":
        if len(row) == 8:
            return {"strain": row[0], "name": row[1], "strand": row[2],
                    "start": int(row[3]), "end": int(row[4]),
                    "hits": "|".join(row[5:8])}
        elif len(row) == 6:
            return {"strain": row[0], "name": row[1], "strand": row[2],
                    "start": int(row[3]), "end": int(row[4]), "hits": row[5]}
    elif type_ == "sRNA":
        if len(row) == 7:
            return {"strain": row[0], "name": row[1], "strand": row[2],
                    "start": int(row[3]), "end": int(row[4]),
                    "hits": "|".join(row[5:6])}
        elif len(row) == 6:
            return {"strain": row[0], "name": row[1], "strand": row[2],
                    "start": int(row[3]), "end": int(row[4]), "hits": row[5]}

def merge_info(blasts):
    first = True
    finals = []
    for blast in blasts:
        if first:
            repeat = 0
            first = False
            pre_blast = deepcopy(blast)
        else:
            if (pre_blast["strain"] == blast["strain"]) and (
                pre_blast["strand"] == blast["strand"]) and (
                pre_blast["start"] == blast["start"]) and (
                pre_blast["end"] == blast["end"]):
                if (repeat < 2):
                    pre_blast["hits"] = ";".join([pre_blast["hits"], blast["hits"]])
                    repeat += 1
            else:
                repeat = 0
                finals.append(pre_blast)
                pre_blast = blast.copy()
    finals.append(pre_blast)
    return finals

def compare_srna_table(srna_tables, srna, final, min_len, max_len):
    for table in srna_tables:
        tsss = []
        pros = []
        cands = []
        if (table["strain"] == srna.seq_id) and (
            table["strand"] == srna.strand) and (
            table["start"] == srna.start) and (
            table["end"] == srna.end):
            final = dict(final, **table)
            datas = table["tss_pro"].split(";")
            tsss.append(table["start"])
            pros.append(table["end"])
            for data in datas:
                if "TSS" in data:
                    if table["start"] != int(data.split(":")[1][:-2]):
                        tsss.append(int(data.split(":")[1][:-2]))
                elif "Cleavage" in data:
                    if table["end"] != int(data.split(":")[1][:-2]):
                        pros.append(int(data.split(":")[1][:-2]))
            for tss in tsss:
                for pro in pros:
                    if ((pro - tss) >= min_len) and (
                        (pro - tss) <= max_len):
                        cands.append("-".join([str(tss), str(pro)]))
            final["candidates"] = ";".join(cands)
            if ("tex" in table["conds"]) and (
                "frag" in table["conds"]):
                final["type"] = "TEX+/-;Fragmented"
            elif ("tex" in table["conds"]):
                final["type"] = "TEX+/-"
            elif ("frag" in table["conds"]):
                final["type"] = "Fragmented"
    return final

def compare_blast(blasts, srna, final, hit):
    for blast in blasts:
        if (srna.seq_id == blast["strain"]) and (
            srna.strand == blast["strand"]) and (
            srna.start == blast["start"]) and (
            srna.end == blast["end"]):
            final[hit] = blast["hits"]
    return final

def compare(srnas, srna_tables, nr_blasts, srna_blasts, min_len, max_len):
    finals = []
    for srna in srnas:
        final = {}
        final["energy"] = srna.attributes["2d_energy"]
        final["nr_hit_num"] = srna.attributes["nr_hit"]
        final["sRNA_hit_num"] = srna.attributes["sRNA_hit"]
        if "sORF" in srna.attributes.keys():
            final["sORF"] = srna.attributes["sORF"]
        else:
            final["sORF"] = "NA"
        if srna.source == "intergenic":
            final["utr"] = "Intergenic"
        else:
            if "&" in srna.attributes["UTR_type"]:
                final["utr"] = "5'UTR_derived;3'UTR_derived"
            elif srna.attributes["UTR_type"] == "5utr":
                final["utr"] = "5'UTR_derived"
            elif srna.attributes["UTR_type"] == "3utr":
                final["utr"] = "3'UTR_derived"
            elif srna.attributes["UTR_type"] == "interCDS":
                final["utr"] = "interCDS"
        final = compare_srna_table(srna_tables, srna, final, min_len, max_len)
        final = compare_blast(nr_blasts, srna, final, "nr_hit")
        final = compare_blast(srna_blasts, srna, final, "sRNA_hit")
        finals.append(final)
    return finals

def print_file(finals, out):
    for final in finals:
        out.write("\t".join([
                  final["strain"], final["name"], str(final["start"]),
                  str(final["end"]), final["strand"], final["tss_pro"],
                  final["candidates"], final["type"], str(final["avg"]),
                  str(final["high"]), str(final["low"]), final["track"],
                  final["energy"], final["utr"], final["sORF"],
                  final["nr_hit_num"], final["sRNA_hit_num"],
                  final["nr_hit"], final["sRNA_hit"]]) + "\n")

def read_table(srna_table_file, nr_blast, srna_blast_file):
    srna_tables = []
    nr_blasts = []
    srna_blasts = []
    f_h = open(srna_table_file, "r")
    for row in csv.reader(f_h, delimiter='\t'):
        srna_tables.append(import_data(row, "gff"))
    f_h.close()
    f_h = open(nr_blast, "r")
    for row in csv.reader(f_h, delimiter='\t'):
        nr_blasts.append(import_data(row, "nr"))
    f_h.close()
    f_h = open(srna_blast_file, "r")
    for row in csv.reader(f_h, delimiter='\t'):
        srna_blasts.append(import_data(row, "sRNA"))
    f_h.close()
    return srna_tables, nr_blasts, srna_blasts

def read_gff(srna_gff):
    srnas = []
    for entry in Gff3Parser().entries(open(srna_gff)):
        srnas.append(entry)
    srnas = sorted(srnas, key=lambda k: (k.seq_id, k.start))
    return srnas

def gen_srna_table(srna_gff, srna_table_file, nr_blast, srna_blast_file,
                   max_len, min_len, out_file):
    srnas = read_gff(srna_gff)
    srna_tables, nr_blasts, srna_blasts = read_table(srna_table_file,
                                          nr_blast, srna_blast_file)
    out = open(out_file, "w")
    out.write("\t".join(["strain", "name", "start", "end", "strand",
                         "TSS/Cleavage_site", "candidates",
                         "lib_type", "best_avg_coverage",
                         "best_highest_coverage", "best_lower_coverage",
                         "track/coverage", "secondary_energy_change",
                         "UTR_derived/Intergenic", "confliction of sORF",
                         "nr_hit_number", "sRNA_hit_number",
                         "nr_hit_top3|ID|e-value", "sRNA_hit|e-value"]) + "\n")
    nr_blasts = merge_info(nr_blasts)
    srna_blasts = merge_info(srna_blasts)
    finals = compare(srnas, srna_tables, nr_blasts,
                     srna_blasts, min_len, max_len)
    sort_finals = sorted(finals, key=lambda x: (x["avg"]), reverse=True)
    print_file(sort_finals, out)

def print_best(detect, out, srna):
    no_print = False
    for value in detect.values():
        if not value:
            no_print = True
    if not no_print:
        out.write(srna.info + "\n")

def gen_best_srna(srna_file, all_srna_hit, energy, hit_nr_num,
                  compare_sorf, out_file):
    srnas = read_gff(srna_file)
    out = open(out_file, "w")
    out.write("##gff-version 3\n")
    for srna in srnas:
        detect = {"energy": False, "TSS": False, "nr_hit": False,
                  "sRNA_hit": False, "sORF": False}
        if "2d_energy" in srna.attributes.keys():
            if float(srna.attributes["2d_energy"]) < energy:
                detect["energy"] = True
        else:
            detect["energy"] = True
        if "with_TSS" in srna.attributes.keys():
            if srna.attributes["with_TSS"] != "NA":
                detect["TSS"] = True
            elif (srna.source == "UTR_derived"):
                if (("3utr" in srna.attributes["UTR_type"]) or (
                     "interCDS" in srna.attributes["UTR_type"])) and (
                     srna.attributes["with_cleavage"] != "NA"):
                    detect["TSS"] = True
        else:
            detect["TSS"] = True
        if "nr_hit" in srna.attributes.keys():
            if (srna.attributes["nr_hit"] == "NA") or (
                int(srna.attributes["nr_hit"]) <= hit_nr_num):
                detect["nr_hit"] = True
        else:
            detect["nr_hit"] = True
        if (compare_sorf):
            if ("sORF" in srna.attributes.keys()):
                if srna.attributes["sORF"] == "NA":
                    detect["sORF"] = True
        else:
            detect["sORF"] = True
        if ("sRNA_hit" in srna.attributes.keys()) and (all_srna_hit):
            if (srna.attributes["sRNA_hit"] != "NA"):
                for key in detect.keys():
                    detect[key] = True
            else:
                count = 0
                for value in detect.values():
                    if value:
                        count += 1
                if count == 4:
                    detect["sRNA_hit"] = True
        else:
            detect["sRNA_hit"] = True
        print_best(detect, out, srna)
    out.close()
