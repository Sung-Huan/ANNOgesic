#!/usr/bin/python

import os	
import sys
import csv
from annogesiclib.gff3 import Gff3Parser

def plus_num(nums, strain, feature):
    nums[strain][feature] += 1
    nums["total"][feature] += 1

def print_file(nums, out, strain):
    out.write(strain + ":")
    out.write("Total terminators = %s\n" % nums["total"])
    out.write("Total terminators which located in gene expression region = %s (percentage of total terminators = %s)\n" % \
              (nums["express"], float(nums["express"]) / float(nums["total"])))
    out.write("Total terminators which are detected by method_1 = %s\n" % str(nums["fr"] + nums["frhp"]))
    out.write("\t(percentage of total terminators = %s, percentage of total terminators which have gene expression = %s)\n" % \
              (float(nums["fr"] + nums["frhp"]) / float(nums["total"]), float(nums["fr"] + nums["frhp"]) / float(nums["express"])))
    out.write("Total terminators which are detected by method_2 = %s\n" % str(nums["hp"] + nums["frhp"]))
    out.write("\t(percentage of total terminators = %s, percentage of total terminators which have gene expression = %s)\n" % \
              (float(nums["hp"] + nums["frhp"]) / float(nums["total"]), float(nums["hp"] + nums["frhp"]) / float(nums["express"])))
    out.write("Total terminators which intersect betweeen these two methods = %s\n" % nums["frhp"])
    out.write("\t(percentage of total terminators = %s, percentage of total terminators which have gene expression = %s)\n" % \
              (float(nums["frhp"]) / float(nums["total"]), float(nums["frhp"]) / float(nums["express"])))
    out.write("Total terminators which are only detected by method_1 = %s\n" % str(nums["fr"]))
    out.write("\t(percentage of total terminators = %s)\n" % str(float(nums["fr"]) / float(nums["total"])))
    out.write("\t(percentage of total terminators of method_1 = %s)\n" % str(float(nums["fr"]) / float(nums["fr"] + nums["frhp"])))
    out.write("\t(percentage of terminators which have gene expression = %s)\n" % str(float(nums["fr"]) / float(nums["express"])))
    out.write("Total terminators which are only detected by method_2 = %s\n" % str(nums["hp"]))
    out.write("\t(percentage of total terminators = %s)\n" % str(float(nums["hp"]) / float(nums["total"])))
    out.write("\t(percentage of total terminators of method_2 = %s)\n" % str(float(nums["hp"]) / float(nums["hp"] + nums["frhp"])))
    out.write("\t(percentage of terminators which have gene expression = %s)\n" % str(float(nums["hp"]) / float(nums["express"])))
    out.write("Total terminators which have dramatic coverage decreasing = %s\n" % str(nums["total_de"]))
    out.write("\t(percentage of total terminators = %s)\n" % str(float(nums["total_de"]) / float(nums["total"])))
    out.write("\t(percentage of terminators which have gene expression = %s)\n" % str(float(nums["total_de"]) / float(nums["express"])))
    out.write("Total terminators which have dramatic coverage decreasing in method_1 = %s\n" % str(nums["de_fr"]))
    out.write("\t(percentage of total terminators = %s)\n" % str(float(nums["de_fr"]) / float(nums["total"])))
    out.write("\t(percentage of total terminators of method_1 = %s)\n" % str(float(nums["de_fr"]) / (float(nums["fr"]) + float(nums["frhp"]))))
    out.write("\t(percentage of terminators which have gene expression = %s)\n" % str(float(nums["de_fr"]) / float(nums["express"])))
    out.write("Total terminators which have dramatic coverage decreasing in method_2 = %s\n" % nums["de_hp"])
    out.write("\t(percentage of total terminators = %s)\n" % str(float(nums["de_hp"]) / float(nums["total"])))
    out.write("\t(percentage of total terminators of method_2 = %s)\n" % str(float(nums["de_hp"]) / (float(nums["hp"]) + float(nums["frhp"]))))
    out.write("\t(percentage of terminators which have gene expression = %s)\n" % str(float(nums["de_hp"]) / float(nums["express"])))
    out.write("Total terminators which have dramatic coverage decreasing in intersection of two methods = %s\n" % str(nums["de_frhp"]))
    out.write("\t(percentage of total terminators = %s)\n" % str(float(nums["de_frhp"]) / float(nums["total"])))
    out.write("\t(percentage of total terminators of method_1 = %s)\n" % str(float(nums["de_frhp"]) / (float(nums["fr"]) + float(nums["frhp"]))))
    out.write("\t(percentage of total terminators of method_2 = %s)\n" % str(float(nums["de_frhp"]) / (float(nums["hp"]) + float(nums["frhp"]))))
    out.write("\t(percentage of total terminators of intersection = %s)\n" % str(float(nums["de_frhp"]) / float(nums["frhp"])))
    out.write("\t(percentage of terminators which have gene expression = %s)\n" % str(float(nums["de_frhp"]) / float(nums["express"])))

def classify_terms(terms, nums, out_d, out_e, pre_strain):
    for term in terms:
        if term.seq_id != pre_strain:
            pre_strain = term.seq_id
            strain = term.seq_id
            nums[strain] = {"fr": 0, "hp": 0, "frhp": 0, "re_fr": 0, "re_hp": 0,
                            "de_fr": 0, "de_hp": 0, "de_frhp": 0, "total": 0,
                            "total_de": 0, "only_fr": 0, "only_hp": 0, "express": 0}
        if term.attributes["coverage_decrease"] == "True":
            out_d.write(term.info + "\n")
        if term.attributes["express"] == "True":
            out_e.write(term.info + "\n")
            plus_num(nums, strain, "express")
        if term.source == "forward_reverse":
            plus_num(nums, strain, "total")
            plus_num(nums, strain, "fr")
            plus_num(nums, strain, "re_fr")
            if term.attributes["coverage_decrease"] == "True":
                plus_num(nums, strain, "de_fr")
                plus_num(nums, strain, "only_fr")
                plus_num(nums, strain, "total_de")
        elif term.source == "TransTermHP":
            plus_num(nums, strain,  "total")
            plus_num(nums, strain,  "hp")
            plus_num(nums, strain,  "re_hp")
            if term.attributes["coverage_decrease"] == "True":
                plus_num(nums, strain, "de_hp")
                plus_num(nums, strain, "only_hp")
                plus_num(nums, strain, "total_de")
        elif term.source == "forward_reverse&TransTermHP":
            plus_num(nums, strain, "total")
            plus_num(nums, strain, "frhp")
            if term.attributes["coverage_decrease"] == "True":
                plus_num(nums, strain, "de_frhp")
                plus_num(nums, strain, "de_fr")
                plus_num(nums, strain, "de_hp")
                plus_num(nums, strain, "total_de")

def stat_term(term_gff, term_table, stat, output_decrease, output_expression):
    terms = []
    nums = {}
    nums["total"] = {"fr": 0, "hp": 0, "frhp": 0, "re_fr": 0, "re_hp": 0,
                     "de_fr": 0, "de_hp": 0, "de_frhp": 0, "total": 0,
                     "total_de": 0, "only_fr": 0, "only_hp": 0, "express": 0}
    pre_strain = ""
    out_te = open(output_expression + ".csv", "w")
    out_td = open(output_decrease + ".csv", "w")
    fh = open(term_table, "r");
    for row in csv.reader(fh, delimiter="\t"):
        if (row[-1] != "NA") and (row[-1] != "No_coverage_decreasing"):
            out_td.write("\t".join(row) + "\n")
            out_te.write("\t".join(row) + "\n")
        if (row[-1] == "No_coverage_decreasing"):
            out_te.write("\t".join(row) + "\n")
    for entry in Gff3Parser().entries(open(term_gff)):
        terms.append(entry) 
    out = open(stat, "w")   
    out_e = open(output_expression + ".gff", "w")
    out_d = open(output_decrease + ".gff", "w")
    out_e.write("##gff-version 3\n")
    out_d.write("##gff-version 3\n")
    classify_terms(terms, nums, out_d, out_e, pre_strain)
    out.write("method_1 is searching the intergenic region between forward strand and reverse strand.\n")
    out.write("method_2 is TransTermH.\n")
    if len(nums) > 2:
        print_file(nums["total"], out, "All strain")
    else:
        for strain, datas in nums.items():
            if strain != "total":
                print_file(datas, out, strain)

