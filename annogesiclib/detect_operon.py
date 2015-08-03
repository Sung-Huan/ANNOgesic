import os
import sys
import math
import csv
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.helper import Helper

def import_to_operon(start, end, strand):
    return {"start": start, "end": end, "strand": strand}

def get_gene_info(cds):
    if "locus_tag" in cds.attributes.keys():
        feature = cds.attributes["locus_tag"]
#    elif "protein_id" in cds.attributes.keys():
#        feature = cds.attributes["protein_id"]
    else:
        strand = Helper().get_strand_name(cds.strand)
        feature = "".join([cds.feature, ":", str(cds.start),
                           "-", str(cds.end), "_", strand])
    return feature

def get_term_feature(ta, data, term_fuzzy, features, datas,
                     ta_check_point, data_check_start, data_check_end):
    jump = False
    if (ta.strand == data.strand) and (
        ta.seq_id == data.seq_id) and (
       (math.fabs(data.start - ta_check_point) <= term_fuzzy) or (
        math.fabs(data.end - ta_check_point) <= term_fuzzy) or (
       (ta_check_point >= data.start) and (
        ta_check_point <= data.end))):
        features["detect"] = True
    if (ta.strand == data.strand) and (
        ta.seq_id == data.seq_id):
        if (ta.start <= data_check_start) and (
            ta.end >= data_check_end):
            features["num"] += 1
            datas.append(data)
        elif (ta_check_point >= data.start) and (
              ta_check_point <= data.end):
            features["num"] += 1
            datas.append(data)
    if (ta.seq_id == data.seq_id) and (
        data.start - term_fuzzy > ta.end):
        jump = True
    return jump

def get_tss_feature(ta, data, features, tss_fuzzy, datas, ta_check_point,
                    data_check_start, data_check_end):
    jump = False
    if (ta.strand == data.strand) and (
        ta.seq_id == data.seq_id) and (
        math.fabs(ta_check_point - data.start) <= tss_fuzzy):
        features["detect"] = True
    if (ta.strand == data.strand) and (
        ta.seq_id == data.seq_id) and (
        ta.start <= data_check_start) and (
        ta.end >= data_check_end):
        features["num"] += 1
        datas.append(data)
    if (ta.seq_id == data.seq_id) and (
        data_check_end > ta.end):
        jump = True
    return jump

def detect_features(ta, inputs, feature, term_fuzzy, tss_fuzzy):
    features = {"num": 0, "detect": False}
    datas = []
    for data in inputs:
        if (feature == "term"):
            if ta.strand == "+":
                jump_term = get_term_feature(ta, data, term_fuzzy, features,
                                    datas, ta.end, data.start,
                                    data.start - term_fuzzy)
            elif ta.strand == "-":
                jump_term = get_term_feature(ta, data, term_fuzzy, features,
                                    datas, ta.start, data.end + term_fuzzy,
                                    data.end)
            if jump_term:
                break
        elif (feature == "tss"):
            if ta.strand == "+":
                jump_tss = get_tss_feature(ta, data, features, tss_fuzzy,
                                   datas, ta.start, data.start + tss_fuzzy,
                                   data.end)
            elif ta.strand == "-":
                jump_tss = get_tss_feature(ta, data, features, tss_fuzzy,
                                   datas, ta.end, data.end,
                                   data.start - tss_fuzzy)
            if jump_tss:
                break
        else:
            if feature == "gene":
                if (ta.strand == data.strand) and (
                    ta.seq_id == data.seq_id) and (
                    ta.start <= data.start) and (
                    ta.end >= data.end) and (
                    (data.feature == "CDS") or (
                     data.feature == "tRNA") or (
                     data.feature == "rRNA") or (
                     data.feature == "sRNA")):
                    features["num"] += 1
                    features["detect"] = True
                    datas.append(data)
                elif (ta.seq_id == data.seq_id) and (
                      data.start > ta.end):
                    break
    return {"data_list": datas, "num_feature": features["num"],
            "with_feature": features["detect"]}

def check_gene(genes, pos, strand):
    conflict = False
    for gene in genes["data_list"]:
        if (gene.strand == strand) and (
            gene.start < pos) and (
            gene.end >= pos):
            conflict = True
            break
    return conflict

def sub_operon(strand, tsss, ta_pos, end, genes, min_length):
    first = True
    operons = []
    operon_pos = ta_pos
    if tsss["with_feature"]:
        if tsss["num_feature"] == 1:
            pass
        else:
            for tss in tsss["data_list"]:
                conflict = check_gene(genes, tss.start, strand)
                if not conflict:
                    operon_pos, first = compute_sub_operon(strand, tss, ta_pos,
                                        first, min_length, end, operons, operon_pos)
    else:
        for tss in tsss["data_list"]:
            conflict = check_gene(genes, tss.start, strand)
            if not conflict:
                operon_pos, first = compute_sub_operon(strand, tss, ta_pos,
                                    first, min_length, end, operons, operon_pos)
    return operons

def compute_sub_operon(strand, tss, ta_pos, first,
                       min_length, end, operons, operon_pos):
    if first:
        operon_pos = ta_pos
        first = False
    else:
        if math.fabs(tss.start - operon_pos) > min_length:
            if strand == "+":
                operons.append(import_to_operon(operon_pos,
                               tss.start - 1, strand))
                operon_pos = tss.start
            else:
                operons.append(import_to_operon(tss.start + 1,
                               operon_pos, strand))
                operon_pos = tss.start
    if (operon_pos != ta_pos) and \
       (math.fabs(end - operon_pos) > min_length):
        if strand == "+":
            operons.append(import_to_operon(operon_pos, end, strand))
        else:
            operons.append(import_to_operon(end, operon_pos, strand))
    return operon_pos, first

def read_gff(ta_file, gff_file, tss_file, terminator_file):
    tas = []
    gffs = []
    tss_gffs = []
    term_gffs = []
    gff_parser = Gff3Parser()
    for ta in gff_parser.entries(open(ta_file)):
        tas.append(ta)
    for entry in gff_parser.entries(open(gff_file)):
        gffs.append(entry)
    for entry in gff_parser.entries(open(tss_file)):
        tss_gffs.append(entry)
    if terminator_file is not False:
        for entry in gff_parser.entries(open(terminator_file)):
            term_gffs.append(entry)
        term_gffs = sorted(term_gffs, key=lambda k: (k.seq_id, k.start))
    tas = sorted(tas, key=lambda k: (k.seq_id, k.start))
    gffs = sorted(gffs, key=lambda k: (k.seq_id, k.start))
    tss_gffs = sorted(tss_gffs, key=lambda k: (k.seq_id, k.start))
    return tas, gffs, tss_gffs, term_gffs

def print_file(ta, operons, out, operon_id, whole_operon, tsss,
               terms, genes, whole_gene):
    if len(operons) == 0:
        out.write("\t".join([operon_id, ta.seq_id,
                  "-".join([str(whole_operon.start), str(whole_operon.end)]),
                  whole_operon.strand, str(len(operons)), "NA",
                  str(tsss["with_feature"]), str(tsss["num_feature"]),
                  str(terms["with_feature"]), str(terms["num_feature"]), "NA",
                  str(genes["num_feature"]), "NA",
                  ", ".join(whole_gene)]) + "\n")
    else:
        for sub in operons:
            sub_gene = []
            num_sub_gene = 0
            for gene in genes["data_list"]:
                if (sub["strand"] == gene.strand) and (
                    sub["start"] <= gene.start) and (
                    sub["end"] >= gene.end):
                    if "locus_tag" in gene.attributes.keys():
                        sub_gene.append(gene.attributes["locus_tag"])
                    else:
                        sub_gene.append("".join([gene.feature, ":",
                                        str(gene.start), "-", str(gene.end),
                                        "_", gene.strand]))
                    num_sub_gene += 1
            if num_sub_gene == 0:
                sub_gene.append("NA")
            out.write("\t".join([operon_id, ta.seq_id,
                      "-".join([str(whole_operon.start),
                                str(whole_operon.end)]),
                      whole_operon.strand, str(len(operons)),
                      "-".join([str(sub["start"]), str(sub["end"])]),
                      str(tsss["with_feature"]), str(tsss["num_feature"]),
                      str(terms["with_feature"]), str(terms["num_feature"]),
                      str(num_sub_gene), str(genes["num_feature"]),
                      ", ".join(sub_gene), ", ".join(whole_gene)]) + "\n")

def operon(ta_file, tss_file, gff_file, terminator_file, tss_fuzzy,
           term_fuzzy, min_length, out_file):
    out = open(out_file, "w")
    out.write("Operon_ID\tstrain\tOperon_position\tstrand\t")
    out.write("Number_of_suboperon\tPosition_of_suboperon\tStart_with_TSS\t")
    out.write("Number_of_TSS\tTerminated_with_terminator\t")
    out.write("Number_of_terminator\tNumber_of_gene_associated_suboperon\t")
    out.write("Number_of_gene_associated_operon\t")
    out.write("Associated_genes_with_suboperon\t")
    out.write("Associated_genes_with_whole_operon\n")
    num_operon = 0
    tas, gffs, tss_gffs, term_gffs = read_gff(ta_file, gff_file, tss_file,
                                              terminator_file)
    for ta in tas:
        whole_gene = []
        check_operon = False
        if (math.fabs(ta.start - ta.end) >= min_length):
            whole_operon = ta
            check_operon = True
            num_operon += 1
        genes = detect_features(ta, gffs, "gene", term_fuzzy, tss_fuzzy)
        tsss = detect_features(ta, tss_gffs, "tss", term_fuzzy, tss_fuzzy)
        if terminator_file is None:
            terms = {"with_feature": "NA", "num_feature": "NA"}
        else:
            terms = detect_features(ta, term_gffs, "term",
                                    term_fuzzy, tss_fuzzy)
        if ta.strand == "+":
            operons = sub_operon(ta.strand, tsss, ta.start,
                                 ta.end, genes, min_length)
        else:
            operons = sub_operon(ta.strand, tsss, ta.end,
                                 ta.start, genes, min_length)
        operon_id = "Operon" + str(num_operon)
        if genes["num_feature"] != 0:
            for gene in genes["data_list"]:
                whole_gene.append(get_gene_info(gene))
        else:
            whole_gene.append("NA")
        if check_operon:
            print_file(ta, operons, out, operon_id, whole_operon,
                       tsss, terms, genes, whole_gene)
    out.close()
