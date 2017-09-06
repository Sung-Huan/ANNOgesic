#!/usr/bin/python

import os
import sys
import csv
import argparse
from gff3 import Gff3Parser

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-k","--ecocyc_file",help="terminators of EcoCyc")
parser.add_argument("-p","--predict_file",help="ANNOgesic predicted terminator file")
args = parser.parse_args()

def main():
    terms = []
    detect = 0
    total = 0
    fh = open(args.ecocyc_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if len(row) >= 4:
            total += 1
            terms.append({"id": row[0], "start": int(row[1]),
                          "end": int(row[2])})
    tot_term = 0
    for pre in Gff3Parser().entries(open(args.predict_file)):
        tot_term += 1
        for ref in terms:
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
    print("the number of published terminators can be detected by ANNOgesic:" + str(detect))
    print("total number of terminators in EcoCyc:" + str(total))
    print("detection rate:" + str(float(detect)/float(total)))

if __name__ == "__main__":
    main()
