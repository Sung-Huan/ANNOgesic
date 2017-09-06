#!/usr/bin/python

import os
import sys
import csv
import argparse
from gff3 import Gff3Parser


__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-p","--predict_file",help="ANNOgesic predicted sRNA file")
parser.add_argument("-r","--refseq_file",help="RefSeq gff file")
args = parser.parse_args()

def main():
    pres = []
    for entry in Gff3Parser().entries(open(args.predict_file)):    
        pres.append(entry)
    refs = []
    num_ref = 0
#    fh = open(args.refseq_file, "r")
#    for row in csv.reader(fh, delimiter='\t'):
#        num_ref += 1
#        refs.append({"start": int(row[0]), "end": int(row[1]), "strand": row[2]})
    for entry in Gff3Parser().entries(open(args.refseq_file)):
        if entry.feature == "ncRNA":
            num_ref += 1
            refs.append(entry)
    detect = 0
    for ref in refs:
        for pre in pres:
#            if pre.strand == ref["strand"]:
#                if ((pre.start >= ref["start"]) and (
#                     pre.end <= ref["end"])) or (
#                    (pre.start <= ref["start"]) and (
#                     pre.end >= ref["end"])) or (
#                    (pre.start >= ref["start"]) and (
#                     pre.start <= ref["end"]) and (
#                     pre.end >= ref["end"])) or (
#                    (pre.start <= ref["start"]) and (
#                     pre.end >= ref["start"]) and (
#                     pre.end <= ref["end"])):
            if pre.strand == ref.strand:
                if ((pre.start >= ref.start) and (
                     pre.end <= ref.end)) or (
                    (pre.start <= ref.start) and (
                     pre.end >= ref.end)) or (
                    (pre.start >= ref.start) and (
                     pre.start <= ref.end) and (
                     pre.end >= ref.end)) or (
                    (pre.start <= ref.start) and (
                     pre.end >= ref.start) and (
                     pre.end <= ref.end)):
                    detect += 1
                    break
    print("the number of published sRNAs which can be detected by ANNOgesic:" + str(detect))
    print("total number of sRNAs in RefSeq:" + str(num_ref))
    print("detection rate:" + str(float(detect) / float(num_ref)))

if __name__ == "__main__":
    main()
