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
parser.add_argument("-k","--regulondb_file",help="TSS of Mendoza-Vargas in RegulonDB")
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
    fh = open(args.regulondb_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if (not row[0].startswith("#")) and (row[-1] != "weak"):
            total += 1
            if row[5] == "forward":
                strand = "+"
            else:
                strand = "-"
            pros[row[1]] = {"start": int(row[3]), "strand": strand}
    for ref in pros.values():
        for pre in tsss:
            if pre.strand == ref["strand"]:
                if (math.fabs(ref["start"] - pre.start) <= args.fuzzy):
                    detect += 1
                    break
    print("the number of published TSSs which can be detected by ANNOgesic:" + str(detect))
    print("the total number of TSSs from Mendoza-Vargas in Regulon DB" + str(total))
    print("detection rate:" + str(float(detect)/float(total)))

if __name__ == "__main__":
    main()
