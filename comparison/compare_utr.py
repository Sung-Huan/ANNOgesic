#!/usr/bin/python

import os
import sys
import csv
import argparse
from gff3 import Gff3Parser

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-k","--regulondb_file",help="the UTRs in RegulonDB")
parser.add_argument("-p5","--predict5_file",help="ANNOgesic predicted 5'UTRs")
parser.add_argument("-p3","--predict3_file",help="ANNOgesic predicted 3'UTRs")
args = parser.parse_args()

def main():
    utr5s = {}
    utr3s = {}
    detect5 = 0
    detect3 = 0
    total5 = 0
    total3 = 0
    fh = open(args.regulondb_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if not row[0].startswith("#"):
            if row[4] == "forward":
                strand = "+"
            else:
                strand = "-"
            if len(row[9]) != 0:
                if row[0] not in utr5s.keys():
                    total5 += 1
                    start = int(row[9].split("-")[0])
                    end = int(row[9].split("-")[-1])
                    utr5s[row[0]] = {"start": start,
                                     "end": end, "strand": strand}
                else:
                    start = int(row[9].split("-")[0])
                    end = int(row[9].split("-")[-1])
                    if start < utr5s[row[0]]["start"]:
                        utr5s[row[0]]["start"] = start
                    if end < utr5s[row[0]]["end"]:
                        utr5s[row[0]]["end"] = end
            if len(row[11]) != 0:
                if row[0] not in utr3s.keys():
                    total3 += 1
                    start = int(row[11].split("-")[0])
                    end = int(row[11].split("-")[-1])
                    utr3s[row[0]] = {"start": start,
                                     "end": end, "strand": strand}
                else:
                    start = int(row[11].split("-")[0])
                    end = int(row[11].split("-")[-1])
                    if start < utr3s[row[0]]["start"]:
                        utr3s[row[0]]["start"] = start
                    if end < utr3s[row[0]]["end"]:
                        utr3s[row[0]]["end"] = end
    for pre in Gff3Parser().entries(open(args.predict5_file)):
        for ref in utr5s.values():
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
                    detect5 += 1
                    break
    for pre in Gff3Parser().entries(open(args.predict3_file)):
        for ref in utr3s.values():
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
                    detect3 += 1
                    break
    print("the number of published 5'UTRs which can be detected by ANNOgesic:" + str(detect5))
    print("the number of published 5'UTRs which can be detected by ANNOgesic:" + str(detect3))
    print("total number of 5'UTRs in RegulonDB:" + str(total5))
    print("total number of 5'UTRs in RegulonDB:" + str(total3))
    print("detection rate:" + str(float(detect5)/float(total5)))
    print("detection rate:" + str(float(detect3)/float(total3)))

if __name__ == "__main__":
    main()
