#!/usr/bin/python

import os
import sys
import csv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_files", nargs="+", help="input file")
parser.add_argument("-o","--output_file",help="output file")
args = parser.parse_args()

def main():
    out = open(args.output_file, "w")
    out.write("Name\tStmts\tMiss\tCover\n----------------------------------------------------\n")
    sts = 0
    miss = 0
    for input_file in args.input_files:
        with open(input_file) as fh:
            for line in fh:
                line = line.strip()
                datas = line.split(" ")
                covers = []
                if datas[0].split("/")[-1] == input_file.split("/")[-1].replace("unitest_test_", ""):
                    for data in datas:
                        if len(data):
                            covers.append(data)
                    sts = sts + int(covers[1])
                    miss = miss + int(covers[2])
                    out.write("\t".join(covers))
                    out.write("\n")
    out.write("----------------------------------------------------\n")
    out.write("Total = " + str(100 - (100*(float(miss) / float(sts)))))
    out.close()
if __name__ == "__main__":
    main()
