#!/usr/bin/python

import os
import sys
import csv
import argparse
from gff3 import Gff3Parser

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-k","--regulondb_file",help="terminators in RegulonDB")
parser.add_argument("-p","--predict_file",help="ANNOgesic predicted terminators")
args = parser.parse_args()

def main():
    terms = []
    detect = 0
    total = 0
    fh = open(args.regulondb_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if row[3] == "forward":
            strand = "+"
        else:
            strand = "-"
        total += 1
        terms.append({"id": row[0], "start": int(row[1]),
                      "end": int(row[2]), "strand": strand})
        if row[3] == "both":
            terms.append({"id": row[0], "start": int(row[1]),
                          "end": int(row[2]), "strand": "+"})
            total += 1
    for pre in Gff3Parser().entries(open(args.predict_file)):
        for ref in terms:
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
    print("the number of published terminators which can be detected by ANNOgesic:" + str(detect))
    print("total number of terminators in RegulonDB:" + str(total))
    print("detection rate:" + str(float(detect)/float(total)))

if __name__ == "__main__":
    main()
