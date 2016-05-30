import os
import csv
import shutil
from annogesiclib.gff3 import Gff3Parser


def filter_utr(srna_gff, srna_table, min_utr):
    out = open("tmp_utr_srna.gff", "w")
    out_ta = open("tmp_utr_srna.csv", "w")
    out.write("##gff-version 3\n")
    gffs = []
    tables = []
    gff_parser = Gff3Parser()
    g_f = open(srna_gff, "r")
    for entry in gff_parser.entries(g_f):
        gffs.append(entry)
    fh = open(srna_table, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if row[0] != "rank":
            if (float(row[7]) >= min_utr):
                tables.append(row)
                out_ta.write("\t".join(row) + "\n")
    for gff in gffs:
        for table in tables:
            if (table[0] == gff.seq_id) and (
                    int(table[2]) == gff.start) and (
                    int(table[3]) == gff.end) and (
                    table[4] == gff.strand):
                out.write(gff.info + "\n")
    g_f.close()
    fh.close()
    os.remove(srna_gff)
    os.remove(srna_table)
    shutil.move("tmp_utr_srna.gff", srna_gff)
    shutil.move("tmp_utr_srna.csv", srna_table)
    out.close()
    out_ta.close()
