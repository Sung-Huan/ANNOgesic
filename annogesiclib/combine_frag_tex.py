#!/usr/bin/python

import os        
import sys
import csv
from annogesiclib.gff3 import Gff3Parser

def modify_position(frag, norm):
    if frag.end < norm.end:
        frag.end = norm.end
    if frag.start > norm.start:
        frag.start = norm.start
    norm.attributes["print"] = True
    frag.attributes["print"] = True

def print_file(data, out, name, num):
    attributes = {}
    attributes["ID"] = "tran" + str(num)
    attributes["Name"] = "Tran_" + name
    attribute_string = ";".join(["=".join(items) for items in attributes.items()])
    out.write("\t".join([str(field) for field in [
                        data.seq_id, data.source, data.feature, data.start,
                        data.end, data.score, data.strand, data.phase,
                        attribute_string]]) + "\n")

def store(data, source, finals):
    data.source = source
    data.attributes["print"] = False
    finals.append(data)

def compare(data1, data2, overlap, tolerance):
    if (data1.seq_id == data2.seq_id) and \
       (data1.strand == data2.strand):
        if (data1.start <= (data2.end + tolerance)) and \
           (data1.start >= data2.start):
            modify_position(data1, data2)
            overlap = True
        elif (data1.end >= (data2.start - tolerance)) and \
             (data1.end <= data2.end):
            modify_position(data1, data2)
            overlap = True
        elif (data1.start <= data2.start) and \
             (data1.end >= data2.end):
            modify_position(data1, data2)
            overlap = True
    return overlap

def combine(frag_file, tex_file, tolerance, output_file):
    frags = []
    norms = []
    finals = []
    out = open(output_file, "w")
    fh = open(frag_file, "r")
    for entry in Gff3Parser().entries(fh):
        entry.attributes["print"] = False
        frags.append(entry)
    fh.close()
    nh = open(tex_file, "r")
    for entry in Gff3Parser().entries(nh):
        entry.attributes["print"] = False
        norms.append(entry)
    nh.close()
    sort_frags = sorted(frags, key=lambda k: (k.seq_id, k.start))
    sort_norms = sorted(norms, key=lambda k: (k.seq_id, k.start))
    for frag in sort_frags:
        overlap = False
        for norm in sort_norms:
            overlap = compare(frag, norm, overlap, tolerance)
        if overlap:
            store(frag, "fragmented_and_normal", finals)
        else:
            store(frag, "fragmented", finals)
    for norm in sort_norms:
        if norm.attributes["print"] is False:
            store(norm, "normal", finals)
    sort_finals = sorted(finals, key=lambda k: (k.seq_id, k.start))
    num = 0
    for tar in sort_finals:
        if tar.attributes["print"] is True:
            continue
        overlap = False
        for ref in sort_finals:
             overlap = compare(tar, ref, overlap, tolerance)
        name='%0*d' % (5, num)
        print_file(tar, out, name, num)
        num += 1
