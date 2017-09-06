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
parser.add_argument("-k","--regulondb_file",help="RegulonDB promoter file")
parser.add_argument("-p","--predict_file",help="ANNOgesic promoter table")
parser.add_argument("-f","--fuzzy", type=int, help="the tolerance nts for comparison")
args = parser.parse_args()

def main():
    pros = {}
    pres = []
    total = 0
    detect = 0
    ph = open(args.predict_file, "r")
    for row in csv.reader(ph, delimiter='\t'): 
        if row[0] != "strain":      
            pres.append({"start": int(row[1]), "strand": row[2]})
    fh = open(args.regulondb_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if (not row[0].startswith("#")) and (row[-1].lower() != "weak"):
            if row[2] == "forward":
                strand = "+"
            else:
                strand = "-"
            if int(row[3]) != 0:
                total += 1
                pros[row[0]] = {"start": int(row[3]), "strand": strand}
    for ref in pros.values():
        for pre in pres:
            if pre["strand"] == ref["strand"]:
                if (math.fabs(ref["start"] - pre["start"]) <= args.fuzzy):
                    detect += 1
                    break
    print("the number of published promoters which can be found by ANNOgesic:" + str(detect))
    print("total number of promoters in RegulonDB:" + str(total))
    print("detection rate:" + str(float(detect)/float(total)))

if __name__ == "__main__":
    main()
