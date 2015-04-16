#!/usr/bin/python

import os	
import sys
import random
import csv
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from annogesiclib.gff3 import Gff3Parser

def get_num_class_color(subs, nums, classes, explode, total, unknown, ini_color,
                        color_elements, color, colors):
    count = 0
    for local, num in subs.items():
        if local == "Unknown":
            tmp_unknown = [local, num]
        else:
            nums["all"].append((float(num) / float(total)) * 100)
            nums["no_unknown"].append((float(num) / float(total - unknown)) * 100)
            classes["all"].append(("{0} {1} ({2}{3})").format(
                    local, num, round((float(num) / float(total)) * 100, 2), "%"))
            classes["no_unknown"].append(("{0} {1} ({2}{3})").format(
                    local, num, round((float(num) / float(total - unknown)) * 100, 2), "%"))
            explode["no_unknown"].append(0)
        explode["all"].append(0)
        if count <= 5:
            colors.append("#" + ini_color[count])
        else:
            while True:
                color_index = random.randint(0, 5)
                if count % 3 == 0:
                    color = color[:4] + color_elements[color_index]
                elif count % 3 == 1:
                    color = color_elements[color_index] + color[2:]
                elif count % 3 == 2:
                    color = color[:2] + color_elements[color_index] + color[4:]
                if (color != "000000") and (color != "FFFFFF") and \
                   ("#" + color not in colors):
                    colors.append("#" + color)
                    break
        count += 1
    nums.append((float(tmp_unknown[1]) / float(total)) * 100)
    classes.append(("{0} {1} ({2}{3})").format(
                    tmp_unknown[0], tmp_unknown[1],
                    round((float(tmp_unknown[1]) / float(total)) * 100, 2), "%"))

def plot(subs, total, unknown, strain, prefix_name):
    color_elements = ["FF", "CC", "99", "66", "33", "00"]
    colors = []
    nums = {"all": [], "no_unknown": []}
    classes = {"all": [], "no_unknown": []}
    explode = {"all": [], "no_unknown": []}
    color = "FFFFFF"
    ini_color = ["FFCCCC", "FFFF99", "99FFCC", "CCFF66", "FFCC33", "CC99FF"]
    get_num_class_color(subs, nums, classes, explode, total, unknown, ini_color,
                        color_elements, color, colors)
    fig = plt.figure(figsize=(12, 10))
    plt.subplot(221)
    plt.pie(nums, explode=explode, colors=colors, startangle=90)
    plt.axis('equal')
    plt.title("Subcellualr localization with Unknown\n", fontsize="16")
    plt.legend(classes, bbox_to_anchor=(1.7, 0.5), loc=10, shadow=True, prop={'size':16})
    plt.subplot(223)
    plt.pie(nums_no_unknown, explode=explode_no_unknown, colors=colors, startangle=90)
    plt.axis('equal')
    plt.title("Subcellualr localization without Unknown\n", fontsize="16")
    plt.legend(classes_no_unknown, bbox_to_anchor=(1.7, 0.5), loc=10, shadow=True, prop={'size':16})
    plt.savefig("_".join([prefix_name, strain + ".png"]))
    plt.clf()
    plt.cla()

def read_table(psortb_file, subs, total_nums, unknown_nums):
    pre_strain = ""
    fh = open(psortb_file, "r")
    for row in csv.reader(fh, delimiter="\t"):
        if pre_strain != row[0]:
            subs[row[0]] = {}
            pre_strain = row[0]
            total_nums[row[0]] = 0
            unknown_nums[row[0]] = 0
        if row[5] not in subs[row[0]].keys():
            subs[row[0]][row[5]] = 1
        else:
            if row[5] == "Unknown":
                unknown_nums[row[0]] += 1
            subs[row[0]][row[5]] += 1
            total_nums[row[0]] += 1
        if row[5] not in subs["all_strain"].keys():
            subs["all_strain"][row[5]] = 1
        else:
            if row[5] == "Unknown":
                unknown_nums["all_strain"] += 1
            subs["all_strain"][row[5]] += 1
            total_nums["all_strain"] += 1

def print_file_and_plot(sub, total_nums, unknown_nums, strain, out_stat):
    plot(subs[strain], total_nums[strain], unknown_nums[strain], strain)
    out_stat.write(strain + ":\n")
    out_stat.write("Total with Unknown is {0}; Total_wihout_Unknown is {1}\n".format(
                   total_nums[strain], total_nums[strain] - unknown_nums[strain]))
    for local, num in sub.items():
        if local != "Unknown":
            out_stat.write("\t{0}\t{1}(include Unknown {2}; exclude Unknonwn {3})\n".format(
            local, num, float(num) / float(total_nums[strain]),
            float(num) / (float(total_nums[strain]) - float(unknown_nums[strain]))))
        else:
            out_stat.write("\t{0}\t{1}(include Unknown {2})\n".format(
            local, num, float(num) / float(total_nums[strain])))

def stat_sublocal(psortb_file, prefix_name, stat_file):
    subs = {}
    subs["all_strain"] = {}
    total_nums = {}
    total_nums["all_strain"] = 0
    unknown_nums = {}
    unknown_nums["all_strain"] = 0
    read_table(psortb_file, subs, total_nums, unknown_nums)
    out_stat = open(stat_file, "w")
    if len(subs) > 2:
        print_file_and_plot(sub["all_strain"], total_nums, unknown_nums, "all_strain", out_stat)
    for strain, sub in subs.items():
        if strain != "all_strain":
            print_file_and_plot(sub, total_nums, unknown_nums, strain, out_stat)
