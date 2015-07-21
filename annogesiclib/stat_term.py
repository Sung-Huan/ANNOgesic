import os	
import sys
import csv
from annogesiclib.gff3 import Gff3Parser

def plus_num(nums, strain, feature):
    nums[strain][feature] += 1
    nums["total"][feature] += 1

def print_method(nums, method_name, method, express, detect, only, out):
    out.write(method_name + "\n")
    out.write("\tTotal {0} terminators = {1}\n".format(
              method_name, nums[method] + nums["frhp"]))
    out.write("\t\t(percentage of total terminator = {0})\n".format(
              float(nums[method] + nums["frhp"]) / float(nums["total"])))
    out.write("\tTotal terminators which only can be detected in {0} = {1}\n".format(
              method_name, nums[method]))
    out.write("\t\t(percentage of total terminators = {0})\n".format(
              float(nums[method]) / float(nums["total"])))
    out.write("\t\t(percentage of total {0} terminators = {1})\n".format(
              method_name, float(nums[method]) / float(nums[method] + nums["frhp"])))
    out.write("\tTotal {0} terminators which located in gene expression region = {1}\n".format(
              method_name, nums[express]))
    out.write("\t\t(percentage of total terminators = {0})\n".format(
              float(nums[express]) / float(nums["total"])))
    out.write("\t\t(percentage of total {0} terminators = {1})\n".format(
              method_name, float(nums[express]) / float(nums[method] + nums["frhp"])))
    out.write("\tTotal {0} terminators which have dramatic coverage decreasing = {1}\n".format(
              method_name, nums[detect]))
    out.write("\t\t(percentage of total terminators = {0})\n".format(
              float(nums[detect]) / float(nums["total"])))
    out.write("\t\t(percentage of total {0} terminators = {1})\n".format(
              method_name, float(nums[detect]) / float(nums[method] + nums["frhp"])))
    out.write("\t\t(percentage of terminators which have gene expression = {0})\n".format(
              float(nums[detect]) / float(nums["total_ex"])))
    out.write("\t\t(percentage of {0} terminators which have gene expression = {1})\n".format(
              method_name, float(nums[detect]) / float(nums[express])))
    out.write("\tTotal terminators which have dramatic coverage decreasing(unique in {0}) = {1}\n".format(
              method_name, nums[only]))
    out.write("\t\t(percentage of total terminator which have dramatic coverage decreasing = {0})\n".format(
              float(nums[only]) / float(nums["total_de"])))
    out.write("\t\t(percentage of total {0} terminator which have dramatic coverage decreasing = {1})\n\n".format(
              method_name, float(nums[only]) / float(nums[detect])))

def print_intersection_number(out, nums, type_):
    out.write("\t\t(percentage of total terminator = {0})\n".format(
              float(nums[type_]) / float(nums["total"])))
    out.write("\t\t(percentage of total method_1 terminator = {0})\n".format(
              float(nums[type_]) / float(nums["fr"])))
    out.write("\t\t(percentage of total method_2 terminator = {0})\n".format(
              float(nums[type_]) / float(nums["hp"])))

def print_intersection_express(out, nums, type_):
    out.write("\t\t(percentage of total terminator which have gene expression = {0})\n".format(
              float(nums[type_]) / float(nums["total_ex"])))
    out.write("\t\t(percentage of total method_1 terminator which have gene expression = {0})\n".format(
              float(nums[type_]) / float(nums["ex_fr"])))
    out.write("\t\t(percentage of total method_2 terminator which have gene expression = {0})\n".format(
              float(nums[type_]) / float(nums["ex_hp"])))

def print_file(nums, out, strain):
    out.write(strain + ":\n")
    out.write("Combine two methods:\n")
    out.write("\tTotal terminators = {0}\n".format(nums["total"]))
    out.write("\tTotal terminators which located in gene expression region = {0}\n".format(nums["total_ex"]))
    out.write("\t\t(percentage of total terminator = {0})\n".format(
              float(nums["total_ex"]) / float(nums["total"])))
    out.write("\tTotal terminators which have dramatic coverage decreasing = {0}\n".format(nums["total_de"]))
    out.write("\t\t(percentage of total terminators = {0})\n".format(
              float(nums["total_de"]) / float(nums["total"])))
    out.write("\t\t(percentage of terminators which have gene expression = {0})\n\n".format(
              float(nums["total_de"]) / float(nums["total_ex"])))
    print_method(nums, "method_1", "fr", "ex_fr", "de_fr", "only_de_fr", out)
    print_method(nums, "method_2", "hp", "ex_hp", "de_hp", "only_de_hp", out)
    out.write("intersection two methods:\n")
    out.write("\tTotal terminators which overlap with two methods = {0}\n".format(nums["frhp"]))
    print_intersection_number(out, nums, "frhp")
    out.write("\tTotal overlaped terminators which located in gene expression region = {0}\n".format(
              nums["ex_frhp"]))
    print_intersection_number(out, nums, "ex_frhp")
    print_intersection_express(out, nums, "ex_frhp")
    out.write("\tTotal overlaped terminators which have dramatic coverage decreasing = {0}\n".format(
              nums["de_frhp"]))
    print_intersection_number(out, nums, "de_frhp")
    print_intersection_express(out, nums, "de_frhp")
    out.write("\t\t(percentage of total terminator which have coverage dramatic decreasing = {0})\n".format(
              float(nums["de_frhp"]) / float(nums["total_de"])))
    out.write("\t\t(percentage of total method_1 terminator which have coverage dramatic decreasing = {0})\n".format(
              float(nums["de_frhp"]) / float(nums["de_fr"])))
    out.write("\t\t(percentage of total method_2 terminator which have coverage dramatic decreasing = {0})\n".format(
              float(nums["de_frhp"]) / float(nums["de_hp"])))

def classify_terms(terms, nums, out_d, out_e, pre_strain):
    for term in terms:
        if term.seq_id != pre_strain:
            pre_strain = term.seq_id
            strain = term.seq_id
            nums[strain] = {"fr": 0, "hp": 0, "frhp": 0, "ex_fr": 0, "ex_hp": 0, "ex_frhp": 0,
                            "de_fr": 0, "de_hp": 0, "de_frhp": 0, "total": 0,
                            "total_de": 0, "total_ex": 0, "only_de_fr": 0, "only_de_hp": 0, 
                            "only_ex_fr": 0, "only_ex_hp": 0, "de_frhp": 0, "ex_frhp": 0}
        if term.attributes["coverage_decrease"] == "True":
            out_d.write(term.info + "\n")
        if term.attributes["express"] == "True":
            out_e.write(term.info + "\n")
        if term.source == "forward_reverse":
            plus_num(nums, strain, "total")
            plus_num(nums, strain, "fr")
            if term.attributes["coverage_decrease"] == "True":
                plus_num(nums, strain, "de_fr")
                plus_num(nums, strain, "only_de_fr")
                plus_num(nums, strain, "total_de")
            if term.attributes["express"] == "True":
                plus_num(nums, strain, "ex_fr")
                plus_num(nums, strain, "only_ex_fr")
                plus_num(nums, strain, "total_ex")
        elif term.source == "TransTermHP":
            plus_num(nums, strain,  "total")
            plus_num(nums, strain,  "hp")
            if term.attributes["coverage_decrease"] == "True":
                plus_num(nums, strain, "de_hp")
                plus_num(nums, strain, "only_de_hp")
                plus_num(nums, strain, "total_de")
            if term.attributes["express"] == "True":
                plus_num(nums, strain, "ex_hp")
                plus_num(nums, strain, "only_ex_hp")
                plus_num(nums, strain, "total_ex")
        elif term.source == "forward_reverse&TransTermHP":
            plus_num(nums, strain, "total")
            plus_num(nums, strain, "frhp")
            if term.attributes["coverage_decrease"] == "True":
                plus_num(nums, strain, "de_frhp")
                plus_num(nums, strain, "de_fr")
                plus_num(nums, strain, "de_hp")
                plus_num(nums, strain, "total_de")
            if term.attributes["express"] == "True":
                plus_num(nums, strain, "ex_frhp")
                plus_num(nums, strain, "ex_fr")
                plus_num(nums, strain, "ex_hp")
                plus_num(nums, strain, "total_ex")

def stat_term(term_gff, term_table, stat, output_decrease, output_expression):
    terms = []
    nums = {}
    nums["total"] = {"fr": 0, "hp": 0, "frhp": 0, "ex_fr": 0, "ex_hp": 0, "ex_frhp": 0,
                     "de_fr": 0, "de_hp": 0, "de_frhp": 0, "total": 0,
                     "total_de": 0, "total_ex": 0, "only_de_fr": 0, "only_de_hp": 0,
                     "only_ex_fr": 0, "only_ex_hp": 0, "de_frhp": 0, "ex_frhp": 0}
    pre_strain = ""
    out_te = open(output_expression + ".csv", "w")
    out_td = open(output_decrease + ".csv", "w")
    fh = open(term_table, "r");
#    print(term_gff)
#    with open(term_gff) as fh:
#        for line in fh:
#            print(line)
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
    out.write("method_2 is TransTermHP.\n")
    if len(nums) > 2:
        print_file(nums["total"], out, "All strain")
    else:
        for strain, datas in nums.items():
            if strain != "total":
                print_file(datas, out, strain)

