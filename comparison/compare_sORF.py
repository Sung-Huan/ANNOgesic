#!/usr/bin/python

import os
import sys
import csv
import argparse
from gff3 import Gff3Parser

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-k","--benchmark_file",help="the benchmarking set of sORF")
parser.add_argument("-p","--predict_file",help="ANNOgesic predicted sORF file")
args = parser.parse_args()

def main():
    sorfs = []
    pres = []
    num_ref = 0
    detect = 0
    for sorf in Gff3Parser().entries(open(args.benchmark_file)):
        num_ref += 1
        sorfs.append(sorf)       
    for pre in Gff3Parser().entries(open(args.predict_file)):
        pres.append(pre)
    for sorf in sorfs:
        for pre in pres:
            if pre.strand == sorf.strand:
                if ((pre.start >= sorf.start) and (
                     pre.end <= sorf.end)) or (
                    (pre.start <= sorf.start) and (
                     pre.end >= sorf.end)) or (
                    (pre.start >= sorf.start) and (
                     pre.start <= sorf.end) and (
                     pre.end >= sorf.end)) or (
                    (pre.start <= sorf.start) and (
                     pre.end >= sorf.start) and (
                     pre.end <= sorf.end)):
                    detect += 1
                    sorf.attributes["detect"] = True
                    break
    print("the number of known sORFs which can be detected by ANNOgesic:" + str(detect))
    print("the total number of known sORFs:" + str(num_ref))
    print("the detection rate:" str(float(detect) / float(num_ref)))


if __name__ == "__main__":
    main()
