#!/usr/bin/python

import os        
import sys
import csv
import itertools
from collections import defaultdict
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
from annogesiclib.gff3 import Gff3Parser

def plot(pri, sec, anti, inter, orph, total, total_more, name, feature_name, file_type):
    tsss = [pri, sec, anti, inter, orph]
    ind = np.arange(5)  # the x locations for the groups
    width = 0.5       # the width of the bars
    fig, ax = plt.subplots()
    if feature_name == "processing site":
        plt.text(0.85, 0.95, "Total processing sites", ha='center', va='center', transform=ax.transAxes)
        plt.text(0.85, 0.9, str(total), ha='center', va='center', transform=ax.transAxes)
    elif feature_name == "TSS":
        plt.text(0.9, 0.95, "Total TSSs", ha='center', va='center', transform=ax.transAxes)
        plt.text(0.9, 0.9, str(total), ha='center', va='center', transform=ax.transAxes)
    rects = ax.bar(ind, tsss, width, color='#9999FF')
    ax.set_ylabel("the number of " + feature_name)
    ax.set_xticks(ind + width/2)
    ax.set_xticklabels(('Primary', 'Secondary', 'Antisense', 'Internal', 'Orphan'))
    ax.set_xlabel("The type of " + feature_name)
    for rect in rects:
        height = rect.get_height()
        plt.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                ha='center', va='bottom')
    plt.savefig(file_type + "_class_" + name + ".png")
    

def stat(tsss, strain, feature_name, out_stat, file_type):
    tss_type = {"Primary": [], "Secondary": [], "Internal": [], "Antisense": [], "Orphan": []}
    num_tss = 0
    num_tss_more = 0
    for entry in tsss:
        num_tss += 1
        if entry.attributes["type"].find("Primary") != -1: 
            tss_type["Primary"].append(num_tss)
        if entry.attributes["type"].find("Secondary") != -1: 
            tss_type["Secondary"].append(num_tss)
        if entry.attributes["type"].find("Antisense") != -1: 
            tss_type["Antisense"].append(num_tss)
        if entry.attributes["type"].find("Internal") != -1: 
            tss_type["Internal"].append(num_tss)
        if entry.attributes["type"].find("Orphan") != -1: 
            tss_type["Orphan"].append(num_tss)
    for key in tss_type.keys():
        num_tss_more = num_tss_more + len(tss_type[key])
    plot(len(tss_type["Primary"]), len(tss_type["Secondary"]), 
         len(tss_type["Antisense"]), len(tss_type["Internal"]), 
         len(tss_type["Orphan"]), num_tss, num_tss_more, 
         strain, feature_name, file_type)
    out_stat.write(strain + ":\n")
    out_stat.write("total number of {0} (if one {1} belong to two class, it count two times) = {2}\n".format(
                   feature_name, feature_name, num_tss_more))
    out_stat.write("total number of unique {0} (if one {1} belong to two class, it count only one time) = {2}\n".format(
                   feature_name, feature_name, num_tss))
    for it in range(1,5):
        for tss in itertools.combinations(tss_type.keys(), it):
            union = []
            for key in tss:
                union = list(set(tss_type[key]) | set(union))
            out_stat.write("{0} = {1} ({2})\n".format(
                           '-'.join(tss), len(union), float(len(union))/float(num_tss)))
    out_stat.write("\n")

def stat_tsspredator(tss_file, file_type, stat_file):
    if file_type == "processing":
        feature_name = "processing site"
    else:
        feature_name = "TSS"
    num_tss = 0
    num_tss_more = 0
    union = []
    tss_type = {"Primary": [], "Secondary": [], "Internal": [], "Antisense": [], "Orphan": []}
    tsss = []
    tsss_strain = {}
    pre_seq_id = ""
    out_stat = open(stat_file, "w")
    gff_parser = Gff3Parser()
    for entry in gff_parser.entries(open(tss_file)):    
        if entry.seq_id != pre_seq_id:
            pre_seq_id = entry.seq_id
            tsss_strain[entry.seq_id] = []
        tsss_strain[entry.seq_id].append(entry)
        tsss.append(entry)
    tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start))
    if len(tsss_strain) > 1:
        stat(tsss, "All_strains", feature_name, out_stat, file_type)
    for strain in tsss_strain.keys():
        stat(tsss_strain[strain], strain, feature_name, out_stat, file_type)
