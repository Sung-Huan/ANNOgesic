#!/usr/bin/python

import os
import sys
import csv
from transaplib.gff3 import Gff3Parser

def uni(tas, genes, out):
    start_tmp=0
    stop_tmp=0
    check=0
    for ta in tas:
        for gene in genes:
            if (ta.strand == gene.strand) and \
               (ta.seq_id == gene.seq_id):
                if ((ta.start < gene.start) and (ta.end > gene.start) and (ta.end < gene.end)) or \
                   ((ta.start > gene.start) and (ta.end < gene.end)) or \
                   ((ta.start > gene.start) and (ta.start < gene.end) and (ta.end > gene.end)):
                    check = 1
        if (check == 0) and (start_tmp != ta.start) and (stop_tmp != ta.end):
            out.write(ta.info + "\n")
            start_tmp = ta.start
            stop_tmp = ta.end
        check = 0

def overlap(tas, genes, print_list, out):
    start_tmp=0
    stop_tmp=0
    printed = 0
    check=0
    combine = False
    for gene in genes:
        start_tmp = 0
        stop_tmp = 0
        for ta in tas:
            if (ta.strand == gene.strand) and \
               (ta.seq_id == gene.seq_id):
                if ((ta.start < gene.start) and (ta.end > gene.start) and (ta.end < gene.end)) or \
                   ((ta.start > gene.start) and (ta.end < gene.end)) or \
                   ((ta.start > gene.start) and (ta.start < gene.end) and (ta.end > gene.end)):
                    check = 1
                    if ta in print_list:
                        printed = 1
                    else:
                        print_list.append(ta)
                    tmp_ta = ta
                    if start_tmp == 0:
                        start_tmp = ta.start
                        stop_tmp = ta.end
                    else:
                        combine = True
                        if stop_tmp < ta.end:
                            stop_tmp = ta.end
                if (ta.start > gene.end) and (start_tmp != 0):
                    check = 0
                    if combine or (printed != 1):
                        out.write("\t".join([str(field) for field in [ \
                                  ta.seq_id, ta.source, ta.feature, start_tmp, stop_tmp, \
                                  ta.score, ta.strand, ta.phase, ta.attribute_string]]) + "\n")
                    combine = False
                    printed = 0
                    break
        if (start_tmp != 0) and (check != 0):
            if combine or (printed != True):
                out.write('\t'.join([str(field) for field in [ \
                          ta.seq_id, ta.source, ta.feature, start_tmp, stop_tmp, \
                          ta.score, ta.strand, ta.phase, ta.attribute_string]]) + "\n")
def Fill_gap(gff_file, ta_file, type_, output):
    tas = []
    genes = []
    print_list = []
    ta_f = open(ta_file, "r");
    gff_f = open(gff_file, "r");
    for entry in Gff3Parser().entries(ta_f):
        tas.append(entry)
    ta_f.close()
    for entry in Gff3Parser().entries(gff_f):
        if entry.feature == "gene":
            genes.append(entry)
    gff_f.close()
    out = open(output, "w")
    if type_ == "overlap":
        overlap(tas, genes, print_list, out)
    elif type_ == "uni":
        uni(tas, genes, out)

def check_confirm(entry, confirm, check):
    if (entry.end >= confirm.end) or \
       ((entry.end >= confirm.start) and \
       (entry.end <= confirm.end)):
        if entry.start < confirm.start:
            confirm.start = entry.start
        if entry.end > confirm.end:
            confirm.end = entry.end
        check = True
    return check

def Longer_TA(TA_file, length, out_file):
    tas = []
    confirms = []
    check = False
    for entry in Gff3Parser().entries(open(TA_file)):
        tas.append(entry)
    tas = sorted(tas, key=lambda k: (k.seq_id, k.start))
    for entry in tas:
        if len(confirms) == 0:
            confirms.append(entry)
        else:
            for confirm in confirms:
                if (entry.strand == confirm.strand) and \
                   (entry.seq_id == confirm.seq_id):
                    if entry.start <= confirm.start:
                        check = check_confirm(entry, confirm, check)
                    elif (entry.start > confirm.start) and \
                         (entry.start < confirm.end):
                        check = check_confirm(entry, confirm, check)
                    if check:
                        break
            if check:
                check = False
            else:
                confirms.append(entry)
    out = open(out_file, "w")
    num = 0
    for confirm in confirms:
        if (confirm.end - confirm.start) >= length:
            confirm.attributes["ID"] = "tran" + str(num)
            confirm.attributes["Name"] = "Transcript_" + ('%0*d' % (5, num))
            attribute_string = ";".join(
                ["=".join(items) for items in confirm.attributes.items()])
            out.write("\t".join([confirm.info_without_attributes, attribute_string]) + "\n")
            num += 1
