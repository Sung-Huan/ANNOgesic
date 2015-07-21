#!/usr/bin/python

import os
import sys
import csv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i","--table",help="input file")
parser.add_argument("-o","--output_all", default=False, action="store_true", help="output file")
args = parser.parse_args()


def import_data(row):
    datas = row[0].split("|")
    return{"strain": datas[1], "strand": datas[2],
           "associate": datas[3], "start_seq": int(datas[4]),
           "end_seq": int(datas[5]), "rfam": row[1], "e": row[2],
           "start_align": int(row[3]), "end_align": int(row[4]),
           "info": row[0], "ID": datas[0]}

def main():
    first = True
    rbss = []
    out = open("tmp.csv", "w")
    out.write("#ID|strain|strand|associated_CDS|start|end\tRfam\te_value\tstart_align\tend_align\n")
    if args.output_all:
        with open(args.table) as fh:
            for line in fh:
                line = line.strip()
                if first:
                    first = False
                    rbss.append(line)
                    out.write(line + "\n")
                else:
                    if line not in rbss:
                        rbss.append(line)
                        out.write(line + "\n")
    else:
        fh = open(args.table, "r")
        for row in csv.reader(fh, delimiter='\t'):
            rbss.append(import_data(row))
        for rbs1 in rbss:
            print(rbs1)
            repeat = False
            if "print" not in rbs1.keys():
                rbs1["print"] = True
                for rbs2 in rbss:
                    print("AAA")
                    print(rbs2)
                    print("BBB")
                    if (rbs1["strain"] == rbs2["strain"]) and \
                       (rbs1["strand"] == rbs2["strand"]) and \
                       (rbs1["ID"] == rbs2["ID"]):
                        if "print" not in rbs2.keys():
                            rbs2["print"] = True
                            repeat = True
                print("CCC")
                print(repeat)
                if not repeat:
                    out.write("\t".join([rbs1["info"], rbs1["rfam"], rbs1["e"],
                                         str(rbs1["start_align"]),
                                         str(rbs1["end_align"])]) + "\n")
    out.close()

if __name__ == "__main__":
    main()
