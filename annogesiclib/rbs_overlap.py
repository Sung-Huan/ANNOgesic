import os
import csv
import shutil
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.helper import Helper


def read_gff(gff_file, type_):
    cdss = []
    g_h = open(gff_file)
    for entry in Gff3Parser().entries(g_h):
        if (Helper().feature_without_notgene(entry)):
            if (type_ == "riboswitch") and (entry.feature != "riboswitch"):
                cdss.append(entry)
            elif (type_ == "thermometer") and (
                    entry.feature != "RNA_thermometer"):
                cdss.append(entry)
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    g_h.close()
    return cdss


def check_repeat(tab, strain, strand, start, end, fuzzy):
    start = start + fuzzy
    end = end - fuzzy
    if (tab["strain"] == strain) and (
            tab["strand"] == strand):
        if ((tab["start"] <= start) and (
                 tab["end"] >= end)) or (
                (tab["start"] >= start) and (
                 tab["end"] <= end)) or (
                (tab["start"] <= start) and (
                 tab["end"] <= end) and (
                 tab["end"] >= start)) or (
                (tab["start"] >= start) and (
                 tab["start"] <= end) and (
                 tab["end"] >= end)):
            return True
    return False


def rbs_overlap(table_file, gff_file, type_, fuzzy):
    tmp_tab = table_file + "_tmp"
    cdss = read_gff(gff_file, type_)
    out = open(tmp_tab, "w")
    fh = open(table_file, "r")
    tables = []
    for row in csv.reader(fh, delimiter='\t'):
        if not row[0].startswith("#"):
            tables.append({"strain": row[1], "strand": row[2],
                           "start": int(row[4]), "end": int(row[5]),
                           "info": "\t".join(row)})
    fh.close()
    for tab in tables:
        overlap = False
        for cds in cdss:
            overlap = check_repeat(tab, cds.seq_id, cds.strand,
                                   cds.start, cds.end, fuzzy)
            if overlap:
                break
        for com in tables:
            if tab != com:
                repeat = check_repeat(tab, com["strain"], com["strand"],
                                      com["start"], com["end"], 0)
                if (not overlap):
                    if ((repeat) and (
                            "print" not in tab.keys()) and (
                            "print" not in com.keys())) or (
                            not repeat):
                        overlap = False
                    else:
                        overlap = True
        if not overlap:
            tab["print"] = True
            out.write(tab["info"] + "\n")
    out.close()
    os.remove(table_file)
    shutil.move(tmp_tab, table_file)
