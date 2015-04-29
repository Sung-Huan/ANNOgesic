#!/usr/bin/python

import os        
import sys
import math
import csv
from annogesiclib.gff3 import Gff3Parser

def import_to_operon(start, end, strand):
    return {"start": start, "end": end, "strand": strand}

def get_term_feature(ta, data, term_fuzzy, features, datas,
                     ta_check_point, data_check_start, data_check_end):
    jump = False
    if (ta.strand == data.strand) and \
       (ta.seq_id == data.seq_id) and \
       ((math.fabs(data.start - ta_check_point) <= term_fuzzy) or \
       (math.fabs(data.end - ta_check_point) <= term_fuzzy) or \
       ((ta_check_point >= data.start) and (ta_check_point <= data.end))):
        features["detect"] = True
    if (ta.strand == data.strand) and \
       (ta.seq_id == data.seq_id) and \
       (ta.start <= data_check_start) and \
       (ta.end >= data_check_end):
        features["num"] += 1
        datas.append(data)
    elif (ta.seq_id == data.seq_id) and \
         (data.start - term_fuzzy > ta.end):
        jump = True
    return jump

def get_tss_feature(ta, data, features, tss_fuzzy, datas, ta_check_point,
                    data_check_start, data_check_end):
    jump = False
    if (ta.strand == data.strand) and \
       (ta.seq_id == data.seq_id) and \
       (math.fabs(ta_check_point - data.start) <= tss_fuzzy):
        features["detect"] = True
    if (ta.strand == data.strand) and \
       (ta.seq_id == data.seq_id) and \
       (ta.start <= data_check_start) and \
       (ta.end >= data_check_end):
        features["num"] += 1
        datas.append(data)
    elif (ta.seq_id == data.seq_id) and \
         (data_check_end > ta.end):
        jump = True
    return jump

def detect_features(ta, inputs, feature, term_fuzzy, tss_fuzzy):
    features = {"num": 0, "detect": False}
    datas = []
    for data in inputs:
        if (feature == "term"):
            if ta.strand == "+":
                jump_term = get_term_feature(ta, data, term_fuzzy, features, datas,
                                        ta.end, data.start, data.start - term_fuzzy)
            elif ta.strand == "-":
                jump_term = get_term_feature(ta, data, term_fuzzy, features, datas,
                                        ta.start, data.end + term_fuzzy, data.end)
            if jump_term:
                break
        elif (feature == "tss"):
            if ta.strand == "+":
                jump_tss = get_tss_feature(ta, data, features, tss_fuzzy, datas,
                                           ta.start, data.start + tss_fuzzy, data.end)
            elif ta.strand == "-":
                jump_tss = get_tss_feature(ta, data, features, tss_fuzzy, datas,
                                           ta.end, data.end, data.start - tss_fuzzy)
            if jump_tss:
                break
        else:
            if feature == "gene":
                if (ta.strand == data.strand) and \
                   (ta.seq_id == data.seq_id) and \
                   (ta.start <= data.start) and \
                   (ta.end >= data.end):
                    features["num"] += 1
                    features["detect"] = True
                    datas.append(data)
                elif (ta.seq_id == data.seq_id) and \
                     (data.start > ta.end):
                    break
    return {"data_list": datas, "num_feature": features["num"], "with_feature": features["detect"]}

def check_gene(genes, pos, strand):
    conflict = False
    for gene in genes["data_list"]:
        if (gene.strand == strand) and \
           (gene.start < pos) and \
           (gene.end >= pos):
            conflict = True
            break
    return conflict

def sub_operon(operons, strand, tsss, ta_pos, end, genes, min_length):
    first = True
    operon_pos = ta_pos
    if tsss["with_feature"]:
        if tsss["num_feature"] == 1:
            pass
        else:
            for tss in tsss["data_list"]:
                conflict = check_gene(genes, tss.start, strand)
                if conflict is False:
                    datas = _sub_operon(operons, strand, tss, ta_pos, first, operon_pos, min_length)
                    first = datas[1]
                    operon_pos = datas[0]
    else:
        for tss in tsss["data_list"]:
            conflict = check_gene(genes, tss.start, strand)
            if conflict is False:
                datas = _sub_operon(operons, strand, tss, ta_pos, first, operon_pos, min_length)
                first = datas[1]
                operon_pos = datas[0]
    if (operon_pos != ta_pos) and \
       (math.fabs(end - operon_pos) > min_length):
        if strand == "+":
            operons.append(import_to_operon(operon_pos, end, strand))
        else:
            operons.append(import_to_operon(end, operon_pos, strand))
    return operons

def _sub_operon(operons, strand, tss, ta_pos, first, operon_pos, min_length):
    if first:
        operon_pos = ta_pos
        first = False
    else:
        if math.fabs(tss.start - operon_pos) > min_length:
            if strand == "+":
                operons.append(import_to_operon(operon_pos, tss.start - 1, strand))
                operon_pos = tss.start
            else:
                operons.append(import_to_operon(tss.start + 1, operon_pos, strand))
                operon_pos = tss.start
    return (operon_pos, first)

def read_gff(TA_file, GFF_file, TSS_file, terminator_file, 
             tas, gffs, tss_gffs, term_gffs):
    gff_parser = Gff3Parser()
    for ta in gff_parser.entries(open(TA_file)):
        tas.append(ta)
    for entry in gff_parser.entries(open(GFF_file)):
        gffs.append(entry)
    for entry in gff_parser.entries(open(TSS_file)):
        tss_gffs.append(entry)
    if terminator_file is not False:
        for entry in gff_parser.entries(open(terminator_file)):
            term_gffs.append(entry)

def print_file(ta, operons, out, Operon_ID, whole_operon, tsss, terms, genes, whole_gene):
    if len(operons) == 0:
        out.write("\t".join([Operon_ID, ta.seq_id, "-".join([str(whole_operon.start), str(whole_operon.end)]),
                             whole_operon.strand, str(len(operons)), "NA", str(tsss["with_feature"]),
                             str(tsss["num_feature"]), str(terms["with_feature"]), str(terms["num_feature"]),
                             "NA", str(genes["num_feature"]), "NA", ", ".join(whole_gene)]) + "\n")
    else:
        for sub in operons:
            sub_gene = []
            num_sub_gene = 0
            for gene in genes["data_list"]:
                if (sub["strand"] == gene.strand) and \
                   (sub["start"] <= gene.start) and \
                   (sub["end"] >= gene.end):
                    sub_gene.append(gene.attributes["locus_tag"])
                    num_sub_gene += 1
            if num_sub_gene == 0:
                sub_gene.append("NA")
            out.write("\t".join([Operon_ID, ta.seq_id, "-".join([str(whole_operon.start), str(whole_operon.end)]),
                                 whole_operon.strand, str(len(operons)), "-".join([str(sub["start"]), str(sub["end"])]),
                                 str(tsss["with_feature"]), str(tsss["num_feature"]), str(terms["with_feature"]),
                                 str(terms["num_feature"]), str(num_sub_gene), str(genes["num_feature"]),
                                 ", ".join(sub_gene), ", ".join(whole_gene)]) + "\n")

def operon(TA_file, TSS_file, GFF_file, terminator_file, tss_fuzzy, term_fuzzy, min_length, out_file):
    out = open(out_file, "w")
    out.write("Operon_ID\tstrain\tOperon_position\tstrand\tNumber_of_suboperon\t")
    out.write("Position_of_suboperon\tStart_with_TSS\tNumber_of_TSS\t")
    out.write("Terminated_with_terminator\tNumber_of_terminator\t")
    out.write("Number_of_gene_associated_suboperon\tNumber_of_gene_associated_operon\t")
    out.write("Associated_genes_with_suboperon\tAssociated_genes_with_whole_operon\n")
    num_operon = 0
    genes = []
    tas = []
    tsss = []
    gffs = []
    tss_gffs = []
    term_gffs = []
    read_gff(TA_file, GFF_file, TSS_file, terminator_file, 
             tas, gffs, tss_gffs, term_gffs)
    tas = sorted(tas, key=lambda k: (k.seq_id, k.start))
    gffs = sorted(gffs, key=lambda k: (k.seq_id, k.start))
    tss_gffs = sorted(tss_gffs, key=lambda k: (k.seq_id, k.start))
    if terminator_file is not False:
        term_gffs = sorted(term_gffs, key=lambda k: (k.seq_id, k.start))
    for ta in tas:
        operons = []
        whole_gene = []
        check_operon = False
        if (math.fabs(ta.start - ta.end) >= min_length):
            whole_operon = ta
            check_operon = True
            num_operon += 1
        genes = detect_features(ta, gffs, "gene", term_fuzzy, tss_fuzzy)
        tsss = detect_features(ta, tss_gffs, "tss", term_fuzzy, tss_fuzzy)
        if terminator_file is False:
            terms = {"with_feature": "NA", "num_feature": "NA"}
        else:
            terms = detect_features(ta, term_gffs, "term", term_fuzzy, tss_fuzzy)
        if ta.strand == "+":
            operons = sub_operon(operons, ta.strand, tsss, ta.start, ta.end, genes, min_length)
        else:
            operons = sub_operon(operons, ta.strand, tsss, ta.end, ta.start, genes, min_length)
        Operon_ID = "Operon" + str(num_operon)
        if genes["num_feature"] != 0:
            for gene in genes["data_list"]:
                whole_gene.append(gene.attributes["locus_tag"])
        else:
            whole_gene.append("NA")
        if check_operon:
            print_file(ta, operons, out, Operon_ID, whole_operon, tsss, terms, genes, whole_gene)
