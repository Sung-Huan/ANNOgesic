#!/usr/bin/python

import os
import sys
import csv

def import_data(row):
    datas = row[0].split("|")
    return{"strain": datas[1], "strand": datas[2],
           "associate": datas[3], "start_seq": int(datas[4]),
           "end_seq": int(datas[5]), "rfam": row[1], "e": row[2],
           "start_align": int(row[3]), "end_align": int(row[4]),
           "info": row[0], "ID": datas[0]}

def modify_table(table, output_all):
    first = True
    rbss = []
    out = open("tmp.csv", "w")
    out.write("#ID|strain|strand|associated_CDS|start|end\tRfam\te_value\tstart_align\tend_align\n")
    if output_all:
        with open(table) as fh:
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
        fh = open(table, "r")
        for row in csv.reader(fh, delimiter='\t'):
            rbss.append(import_data(row))
        for rbs1 in rbss:
            repeat = False
            if "print" not in rbs1.keys():
                rbs1["print"] = True
                for rbs2 in rbss:
                    if (rbs1["strain"] == rbs2["strain"]) and \
                       (rbs1["strand"] == rbs2["strand"]) and \
                       (rbs1["ID"] == rbs2["ID"]):
                        if "print" not in rbs2.keys():
                            rbs2["print"] = True
                            repeat = True
                if not repeat:
                    out.write("\t".join([rbs1["info"], rbs1["rfam"], rbs1["e"],
                                         str(rbs1["start_align"]),
                                         str(rbs1["end_align"])]) + "\n")
    os.remove(table)
    os.rename("tmp.csv", table)
    out.close()

