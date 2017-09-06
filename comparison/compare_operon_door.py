#!/usr/bin/python

import os
import sys
import csv
import argparse

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-d","--door2_file",help="door2 data")
parser.add_argument("-p","--predict_file",help="ANNOgesic gff file")
args = parser.parse_args()

def main():
    pre_op = ""
    operons = []
    nums = {"detect": 0, "total": 0}
    fh = open(args.door2_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if row[0] != "OperonID":
            if row[0] != pre_op:
                if len(pre_op) != 0:
                    nums["total"] += 1
                    operons.append({"start": start, "end": end, "strand": strand})
                    start = int(row[3])
                    end = int(row[4])
                    strand = row[5]
                else:
                    start = int(row[3])
                    end = int(row[4])
                    strand = row[5]
                pre_op = row[0]
            else:
                if start > int(row[3]):
                    start = int(row[3])
                if end < int(row[4]):
                    end = int(row[4])
    sh = open(args.predict_file, "r")
    total_p = 0
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
    print("detected operon by ANNOgesic:" + str(nums["detect"]))
    print("detection rate:" + str(float(nums["detect"]/nums["total"])))
    print("total number of DOOR2:" + str(nums["total"]))

if __name__ == "__main__":
    main()
