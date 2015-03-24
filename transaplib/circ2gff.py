#!/usr/bin/python

import os
import sys
import csv
import argparse

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-c","--circ_file",help="input circRNA table")
parser.add_argument("-d","--depth", default=5, type=int, help="depth")
parser.add_argument("-sr","--start_ratio", default=0.25, type=float, help="the ratio of (read support circ / all read) at starting point.")
parser.add_argument("-er","--end_ratio", default=0.25, type=float, help="the ratio of (read support circ / all read) at end point.")
parser.add_argument("-oa","--out_all", help="output gff of all circRNA")
parser.add_argument("-of","--out_fillter", help="output gff of circRNA after fillter out")
args = parser.parse_args()

def Convert

def import_list(row):
    return({"name": row[0], "strain": row[1],
            "strand": row[2], "start": int(row[3]),
            "end": int(row[4]), "conflict": row[5],
            "depth": int(row[6]), "per_start":float(row[7]),
            "per_end":float(row[8])})

def main():
    circs = []
    out_a = open(args.out_all, "w")
    out_f = open(args.out_fillter, "w")
    out_a.write("##gff-version 3\n")
    out_f.write("##gff-version 3\n")
    fh = open(args.circ_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if row[0] != "ID":
            circs.append(import_list(row))
    circs = sorted(circs, key=lambda k: (k["strain"], k["start"]))
    for circ in circs:
        ID = circ["name"].split("_")[1]
        attribute_string = ";".join(["=".join(items) for items in [
                                       ("ID", "circrna" + ID),
                                       ("name", circ["name"]),
                                       ("support_reads", str(circ["depth"])),
                                       ("read_at_start", str(circ["per_start"])),
                                       ("read_at_end", str(circ["per_end"])),
                                       ("confliction", circ["conflict"])]])
        out_a.write("\t".join([str(field) for field in [
                        circ["strain"], "segemehl", "circRNA", str(circ["start"]),
                        str(circ["end"]), ".",circ["strand"], ".", attribute_string]]) + "\n")
        if (circ["depth"] >= args.depth) and \
           (circ["conflict"] == "NA") and \
           (circ["per_start"] >= args.start_ratio) and \
           (circ["per_end"] >= args.end_ratio):
            out_f.write("\t".join([str(field) for field in [
                        circ["strain"], "segemehl", "circRNA", str(circ["start"]),
                        str(circ["end"]), ".",circ["strand"], ".", attribute_string]]) + "\n")
            
if __name__ == "__main__":
    main()
