#!/usr/bin/python

import os
import sys
import csv
import argparse
import math
from gff3 import Gff3Parser

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-k","--regulondb_file",help="TSS file of Salgado et. al in RegulonDB")
parser.add_argument("-p","--predict_file",help="ANNOgesic predicted TSS file")
parser.add_argument("-f","--fuzzy", type=int, help="tolerance of nts for comparison")
args = parser.parse_args()

def main():
    pros = {}
    tsss = []
    total = 0
    detect = 0
    for entry in Gff3Parser().entries(open(args.predict_file)):
        tsss.append(entry)
    refs = []
    fh = open(args.regulondb_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if not row[0].startswith("#"):
            if row[5] == "forward":
                strand = "+"
            else:
                strand = "-"
            total += 1
            refs.append({"start": int(row[0]),
                          "end": int(row[1]), "strand": strand})
    for ref in refs:
        ref["start"] = ref["start"] - args.fuzzy
        ref["end"] = ref["end"] + args.fuzzy
        for pre in tsss:
            if pre.strand == ref["strand"]:
                if ((pre.start >= ref["start"]) and (
                     pre.end <= ref["end"])) or (
                    (pre.start <= ref["start"]) and (
                     pre.end >= ref["end"])) or (
                    (pre.start >= ref["start"]) and (
                     pre.start <= ref["end"]) and (
                     pre.end >= ref["end"])) or (
                    (pre.start <= ref["start"]) and (
                     pre.end >= ref["start"]) and (
                     pre.end <= ref["end"])):
                    detect += 1
                    break
    print("the number of reported TSSs which can be detected by ANNOgesic:" + str(detect))
    print("total number of TSSs from Salgado et. al in RegulonDB:" + str(detect))
    print("detection rate:" + str(float(detect)/float(total)))

if __name__ == "__main__":
    main()
