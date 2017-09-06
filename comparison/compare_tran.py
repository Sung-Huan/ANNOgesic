#!/usr/bin/python

import os
import sys
import csv
import argparse
from gff3 import Gff3Parser

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-k","--ecocyc_file",help="the transcripts from EcoCyc")
parser.add_argument("-p","--predict_file",help="ANNOgesic predicted transcripts")
args = parser.parse_args()

def main():
    trans = {}
    pres = []
    total = 0
    detect = 0
    for entry in Gff3Parser().entries(open(args.predict_file)):
        pres.append(entry)
    fh = open(args.ecocyc_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if row[0] not in trans.keys():
            total += 1
            trans[row[0]] = {"start": int(row[1]), "end": int(row[2])}
        else:
            if int(row[1]) < trans[row[0]]["start"]:
                trans[row[0]]["start"] = int(row[1])
            if int(row[2]) > trans[row[0]]["end"]:
                trans[row[0]]["end"] = int(row[2])
    for ref in trans.values():
        for pre in pres:
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
    print("the number of published transcripts which can be detected by ANNOgesic:" + str(detect))
    print("total number of transcripts in EcoCyc:" + str(total))
    print("detection rate:" + str(float(detect)/float(total)))

if __name__ == "__main__":
    main()
