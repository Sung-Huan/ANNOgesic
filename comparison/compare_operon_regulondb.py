#!/usr/bin/python

import os
import sys
import csv
import argparse

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-d","--regulondb_file",help="RegulonDB file")
parser.add_argument("-p","--predict_file",help="ANNOgesic operon table")
args = parser.parse_args()

def main():
    pre_op = ""
    operons = []
    nums = {"detect": 0, "total": 0}
    fh = open(args.regulondb_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if (not row[0].startswith("#")) and (row[-1] != "Weak"):
            nums["total"] += 1
            if row[3] == "forward":
                row[3] = "+"
            else:
                row[3] = "-"
            operons.append({"start": int(row[1]), "end": int(row[2]), "strand": row[3]})
    sh = open(args.predict_file, "r")
    uniqs = []
    for row in csv.reader(sh, delimiter='\t'):
        if row[0] != "Operon_ID":
            start = int(row[2].split("-")[0])
            end = int(row[2].split("-")[-1])
            for operon in operons:
                if operon["strand"] == row[3]:
                    if ((operon["start"] <= start) and (
                            operon["end"] >= end)) or (
                            (operon["start"] >= start) and (
                            operon["end"] >= end)) or (
                            (operon["start"] >= start) and (
                            operon["start"] <= end) and (
                            operon["end"] >= end)) or (
                            (operon["start"] <= start) and (
                            operon["end"] >= start) and (
                            operon["end"] <= end)):
                        if operon not in uniqs :
                            nums["detect"] += 1
                            uniqs.append(operon)
                            operon["detect"] = True
                            break
            pre_op = {"start": start, "end": end, "strand": row[3]}
    print("detected operons by ANNOgesic:" + str(nums["detect"]))
    print("total operon in RegulonDB:" + str(nums["total"]))
    print("detection rate:" + str(float(nums["detect"]/nums["total"])))

if __name__ == "__main__":
    main()
