#!/usr/bin/python

import os
import sys
import csv
import argparse

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file",help="input file")
parser.add_argument("-o","--output_file",help="output file")
args = parser.parse_args()

def main():
    num = 1
    with open(args.input_file) as f_h:
        for line in f_h:
            line = line.strip()
            if line.startswith(">"):
                datas = line.split("|")
                if datas[0][1:] == "NA":
                    datas[0] = ">srn_" + str(num)
                    num += 1
                print("|".join(datas[:3]))
            else:
                print(line)


if __name__ == "__main__":
    main()
