#!/usr/bin/python

import os	
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from transaplib.gff3 import Gff3Parser

def ellipse(x, y, angle, face, al, plt):
    ellipse= mpl.patches.Ellipse(xy=(x,y), width=0.8, height=0.3, 
                                 angle=angle, facecolor=face, alpha=al)
    plt.gca().add_artist(ellipse)
    return plt
    
def line(x, y, angle, plt):
    line = mpl.patches.Ellipse(xy=(x,y), width=0.8, height=0.3, 
                               angle=angle, facecolor="none", 
                               edgecolor="#000000", linewidth=3)
    plt.gca().add_artist(line)
    return plt

def text_total(xy, tss_type, num, plt):
    plt.text(xy[0], xy[1], tss_type, ha="center", 
             va="center", fontsize=15, fontweight='bold')
    if tss_type != "Orphan":
        plt.text(xy[0], xy[1] - 0.05, str(num), ha="center", 
                 va="center", fontsize=15, fontweight='bold')

def text(xy, tss_type, num, plt):
    if (tss_type == "Primary") or \
       (tss_type == "Antisense") or \
       (tss_type == "Antisense_Primary"):
        plt.text(xy[0], xy[1], str(num), fontsize=16, ha="center",
                 va="center", fontweight='bold', color="white")
    else:
        plt.text(xy[0], xy[1], str(num), fontsize=16, ha="center",
                 va="center", fontweight='bold')

def check_tss_class(total_types, strain, tss, tss_type):
    if tss_type not in total_types[strain].keys():
        total_types[strain][tss_type] = 0
    if tss_type in tss.attributes["type"]:
        total_types[strain][tss_type] += 1

def import_types(tsss, types, total_types):
    for strain, datas in tsss.items():
        if strain not in types.keys():
            types[strain] = {}
            total_types[strain] = {}
        for tss in datas:
            check_tss_class(total_types, strain, tss, "Primary")
            check_tss_class(total_types, strain, tss, "Secondary")
            check_tss_class(total_types, strain, tss, "Internal")
            check_tss_class(total_types, strain, tss, "Antisense")
            check_tss_class(total_types, strain, tss, "Orphan")
            sorted_types = sorted(tss.attributes["type"].split(" "))
            ty = None
            for tss_type in sorted_types:
                if ty is None:
                    ty = tss_type
                else:
                    if tss_type not in ty:
                        ty = ty + "_" + tss_type
            if ty not in types[strain].keys():
                types[strain][ty] = 0
            types[strain][ty] += 1


def read_gff(tss_file, tsss, tss_num):
    pre_strain = ""
    gff_parser = Gff3Parser()
    for entry in gff_parser.entries(open(tss_file)):
        if pre_strain != entry.seq_id:
            tsss[entry.seq_id] = []
            tss_num[entry.seq_id] = 0
            pre_strain = entry.seq_id
        tsss[entry.seq_id].append(entry)
        tsss["all"].append(entry)
        tss_num[entry.seq_id] += 1
        tss_num["all"] += 1
    for strain in tsss.keys():
        tsss[strain] = sorted(tsss[strain], key=lambda k: (k.seq_id, k.start))

def plot(types, file_type, feature_name, total_types, tss_num):
    for strain, tss_types in types.items():
        if len(types.keys()) <= 2:
            if strain == "all":
                continue
        plt.figure(figsize=(12,6))
        coordinate_total = {"Primary": (0.05, 0.85), "Secondary": (0.2, 0.95),
                            "Internal": (0.575, 0.95), "Antisense": (0.7, 0.85),
                            "Orphan":(0.8, 0.3)}
        if feature_name == "processing site":
            plt.text(0.05, 0.05, "Total processing sites", ha="center",
                     va="center", fontsize=15, fontweight='bold')
            plt.text(0.05, 0, str(tss_num[strain]), ha="center",
                     va="center", fontsize=15, fontweight='bold')
        elif feature_name == "TSS":
            plt.text(0.025, 0.05, "Total TSSs", ha="center",
                     va="center", fontsize=15, fontweight='bold')
            plt.text(0.025, 0, str(tss_num[strain]), ha="center",
                     va="center", fontsize=15, fontweight='bold')
        for tss_type, num in total_types[strain].items():
            text_total(coordinate_total[tss_type], tss_type, num, plt)
        ellipse(0.5, 0.4, 70, "#E83241", 1.0, plt)
        ellipse(0.25, 0.4, -70, "#6648DC", 0.8, plt)
        ellipse(0.38, 0.495, 70, "#13C139", 0.5, plt)
        ellipse(0.37, 0.495, -70, "#E8D632", 0.4, plt)
        circ = mpl.patches.Ellipse(xy=(0.8, 0.2), width=0.09, height=0.15,
                                   facecolor='none', edgecolor="#000000", linewidth=3)
        plt.gca().add_artist(circ)
        line(0.25, 0.4, -70, plt)
        line(0.37, 0.495, -70, plt)
        line(0.5, 0.4, 70, plt)
        line(0.38, 0.495, 70, plt)
        coordinates = {"Primary": (0.15, 0.5), "Secondary": (0.275, 0.75),
                       "Internal": (0.476, 0.75), "Antisense": (0.625, 0.5),
                       "Primary_Secondary": (0.225, 0.625),
                       "Internal_Primary": (0.25, 0.225),
                       "Antisense_Primary": (0.375, 0.075),
                       "Internal_Secondary": (0.375, 0.625),
                       "Antisense_Secondary": (0.5, 0.225),
                       "Antisense_Internal": (0.525, 0.625),
                       "Internal_Primary_Secondary": (0.3, 0.45),
                       "Antisense_Primary_Secondary": (0.42, 0.18),
                       "Antisense_Internal_Primary": (0.335, 0.18),
                       "Antisense_Internal_Secondary": (0.45, 0.45),
                       "Antisense_Internal_Primary_Secondary": (0.375, 0.3),
                       "Orphan": (0.8, 0.19)}
        for tss_type, xy in coordinates.items():
            if tss_type not in tss_types.keys():
                tss_types[tss_type] = 0
            text(xy, tss_type, tss_types[tss_type], plt)
        plt.axis('off')
        plt.savefig(file_type + "_venn_" + strain + ".png")

def Plot_Venn(tss_file, file_type):
    if file_type == "processing":
        feature_name = "processing site"
    else:
        feature_name = "TSS"
    tsss = {"all": []}
    types = {"all": {}}
    tss_num = {"all": 0}
    total_types = {"all": {}}
    read_gff(tss_file, tsss, tss_num)
    import_types(tsss, types, total_types)
    plot(types, file_type, feature_name, total_types, tss_num)
