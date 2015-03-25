#!/usr/bin/python

import sys
import csv
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
from transaplib.gff3 import Gff3Parser


def plot(utr, utr_pri, utr_sec, filename, source, utr_type):
    bin_num = np.arange(0, 300, 5)
    if utr_type == "5utr":
        if source:
            n, bins, hist1 = plt.hist(utr, bin_num, color="#FF9999", label='secondary')
            n, bins, hist2 = plt.hist(utr_pri, bin_num, color="#9999FF", label='primary')
            plt.legend((hist2[0],hist1[0]), ("Primary TSSs", "Secondary TSSs"))
        else:
            n, bins, hist1 = plt.hist(utr, bin_num, color="#FF9999")
        plt.xlabel("5'UTR_length")
    elif utr_type == "3utr":
        n, bins, hist = plt.hist(utr, bin_num, color="#9999FF", label='3\'UTR')
        plt.legend([hist[0]], ["3'UTR"])
        plt.xlabel("3'UTR_length")
    plt.savefig(filename)

def get_feature(cds, cds_name):
    if "protein_id" in cds.attributes.keys():
        cds_name = cds.attributes["protein_id"]
    elif "locus_tag" in cds.attributes.keys():
        cds_name = cds.attributes["locus_tag"]
    else:
        cds_name = "".join([cds.feature, ":", str(cds.start), 
                            "-", str(cds.end), "_", cds.strand])
    return cds_name

def check_ta(tas, utr_start, utr_end, seq_id, strand):
    detect = False
    for ta in tas:
        if (ta.seq_id == seq_id) and \
           (ta.strand == strand):
            if (ta.start <= utr_start) and \
               (ta.end >= utr_end):
                detect = True
                break
    return detect

def import_utr(source, tss, utr_strain, utr_all, start, end, tas, length):
    if source:
        if "Primary" in tss.attributes["type"]:
            utr_strain["pri"][tss.seq_id].append(length)
            utr_all["pri"].append(length)
        elif "Secondary" in tss.attributes["type"]:
            utr_strain["sec"][tss.seq_id].append(length)
            utr_all["sec"].append(length)
    utr_strain["all"][tss.seq_id].append(length)
    utr_all["all"].append(length)
    detect = check_ta(tas, start, end, tss.seq_id, tss.strand)
    return detect

def get_5utr(tss, near_cds, utr_strain, utr_all, tas, 
            num_utr, cds_name, locus_tag, source, out):
    if tss.strand == "+":
        start = tss.start
        end = near_cds.start
        length = end - start
        detect = import_utr(source, tss, utr_strain, utr_all, 
                            start, end, tas, length)
    else:
        start = near_cds.end
        end = tss.end
        length = end - start
        detect = import_utr(source, tss, utr_strain, utr_all, 
                            start, end, tas, length)
    if detect:
        name_utr='%0*d' % (5, num_utr)
        if source is True:
            attribute_string = ";".join(
                 ["=".join(items) for items in [("ID", "_".join(["utr5", str(num_utr)])),
                  ("Name", "_".join(["5'UTR", name_utr])), ("length", str(length)),
                  ("TSS_type", tss.attributes["type"]),
                  ("associated_cds", cds_name),
                  ("associated_gene", locus_tag),
                  ("associated_tss", tss.attributes["Name"])]]) 
            out.write("{0}\t{1}\t5UTR\t{2}\t{3}\t.\t{4}\t.\t{5}\n".format(
                  tss.seq_id, tss.source, start, end, tss.strand, attribute_string))
        else:
            attribute_string = ";".join(
                 ["=".join(items) for items in [("ID", "_".join(["utr5", str(num_utr)])),
                  ("Name", "_".joinn(["5'UTR", name_utr])), ("length", length),
                  ("associated_cds", cds_name),
                  ("associated_gene", locus_tag),
                  ("associated_tss","".join(["TSS:", str(tss.start), "_", tss.strand]))]])
            out.write("{0}\t{1}\t5UTR\t{2}\t{3}\t.\t{4}\t.\t{5}\n".format(
                  tss.seq_id, tss.source, start, end, tss.strand, attribute_string))
        num_utr += 1
    return num_utr

def detect_cds(cdss, gene):
    detect = False
    for cds in cdss:
        if "Parent" in cds.attributes.keys():
            if gene.attributes["ID"] == cds.attributes["Parent"]:
                detect = True
                near_cds = cds
                check_utr = True
                cds_name = get_feature(cds, cds_name)
        else:
            if "locus_tag" in cds.attributes.keys():
                if gene.attributes["locus_tag"] == cds.attributes["locus_tag"]:
                    near_cds = cds
                    cds_name = cds.attributes["locus_tag"]
                    check_utr = True
                    detect = True
            elif (gene.seq_id == cds.seq_id) and \
                 (gene.strand == cds.strand):
                if (gene.start >= cds.start) and \
                   (gene.end <= cds.end):
                    near_cds = cds
                    cds_name = get_feature(cds, cds_name)
                    check_utr = True
                    detect = True
    if detect is False:
        return (None, None, False)
    else:
        return (near_cds, cds_name, check_utr)

def read_file(TSS_file, tsss, GFF_file, cdss, genes, 
              TA_file, tas, term_file, terms):
    gff_f = open(GFF_file, "r")
    for entry in Gff3Parser().entries(gff_f):
        if (entry.feature == "CDS") or \
           (entry.feature == "tRNA") or \
           (entry.feature == "rRNA"):
            cdss.append(entry)
        elif (entry.feature == "gene"):
            genes.append(entry)
    genes = sorted(genes, key=lambda k: (k.seq_id, k.start))
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start))
    for entry in Gff3Parser().entries(open(TA_file)):
        tas.append(entry)
    tas = sorted(tas, key=lambda k: (k.seq_id, k.start))
    if term_file is not False:
        for entry_term in Gff3Parser().entries(open(term_file)):
            terms.append(entry_term)
        terms = sorted(terms, key=lambda k: (k.seq_id, k.start))
    if TSS_file is not False:
        for entry in Gff3Parser().entries(open(TSS_file)):
            tsss.append(entry)
        tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start))


def get_5utr_from_TSSpredator(tss, genes, cdss, check_utr, cds_name):
    print(tss)
    if ("Primary" in tss.attributes["type"]) or \
       ("Secondary" in tss.attributes["type"]):
        for gene in genes:
            if (tss.seq_id == gene.seq_id) and \
               (tss.strand == gene.strand):
                ass_gene = tss.attributes["associated_gene"].split(" ")
                tss_type = tss.attributes["type"].split(" ")
                for index in range(len(tss_type)):
                    if (tss_type[index] == "Primary") or \
                       (tss_type[index] == "Secondary"):
                        locus_tag = ass_gene[index]
                if gene.attributes["locus_tag"] == locus_tag:
                    print("aaa")
                    for cds in cdss:
                        if "Parent" in cds.attributes.keys():
                            if gene.attributes["ID"] == cds.attributes["Parent"]:
                                near_cds = cds
                                check_utr = True
                                cds_name = get_feature(cds, cds_name)
                        else:
                            if "locus_tag" in cds.attributes.keys():
                                if gene.attributes["locus_tag"] == cds.attributes["locus_tag"]:
                                    near_cds = cds
                                    cds_name = cds.attributes["locus_tag"]
                                    check_utr = True
                            elif (gene.seq_id == cds.seq_id) and \
                                 (gene.strand == cds.strand):
                                if (gene.start >= cds.start) and \
                                   (gene.end <= cds.end):
                                    near_cds = cds
                                    cds_name = get_feature(cds, cds_name)
                                    check_utr = True
    utr_datas = {"check": check_utr, "cds_name": cds_name, "near_cds": near_cds, "locus": locus_tag}
    return (utr_datas)

def get_5utr_from_other(tss, genes, cdss, check_utr, cds_name):
    for gene in genes:
        if (tss.seq_id == gene.seq_id) and \
           (tss.strand == gene.strand):
            if (tss.strand == "+") and \
               ((gene.start - tss.start) <= 300) and \
               (tss.start <= gene.start):
                if "locus_tag" in gene.attributes.keys():
                    locus_tag = gene.attributes["locus_tag"]
                else:
                    locus_tag = gene.attributes["ID"]
                datas = detect_cds(cdss, gene)
                near_cds = datas[0]
                cds_name = datas[1]
                check_utr = datas[2]
                break
            elif (tss.strand == "-") and \
                 ((tss.start - gene.end) <= 300) and \
                 (tss.start >= gene.end):
                if "locus_tag" in gene.attributes.keys():
                    locus_tag = gene.attributes["locus_tag"]
                else:
                    locus_tag = gene.attributes["ID"]
                datas = detect_cds(cdss, gene)
                near_cds = datas[0]
                cds_name = datas[1]
                check_utr = datas[2]
                break
    utr_datas = {"check": check_utr, "cds_name": cds_name, "near_cds": near_cds, "locus": locus_tag}
    return (utr_datas)

def Detect_5UTR(TSS_file, GFF_file, TA_file, source, out_file):
    num_utr = 0
    genes = []
    tas = []
    cdss = []
    tsss = []
    utr_all = {"all": [], "pri": [], "sec": []}
    utr_strain = {"all": {}, "pri": {}, "sec": {}}
    pre_seq_id = ""
    out = open(out_file, "w")
    out.write("##gff-version 3\n")
    read_file(TSS_file, tsss, GFF_file, cdss, genes, TA_file, tas, False, False)
    for tss in tsss:
        check_utr = False
        cds_name = "NA"
        if (source is True):
            utr_datas = get_5utr_from_TSSpredator(tss, genes, cdss, check_utr, cds_name)
        else:
            utr_datas = get_5utr_from_other(tss, genes, cdss, check_utr, cds_name)
        if utr_datas["check"]:
            if tss.seq_id != pre_seq_id:
                pre_seq_id = tss.seq_id
                utr_strain["pri"][tss.seq_id] = []
                utr_strain["sec"][tss.seq_id] = []
                utr_strain["all"][tss.seq_id] = []
            num_utr = get_5utr(tss, utr_datas["near_cds"], utr_strain, utr_all, tas, num_utr, 
                              utr_datas["cds_name"], utr_datas["locus"], source, out)
    name = (GFF_file.split("/"))[-1].replace(".gff", "")
    plot(utr_all["all"], utr_all["pri"], utr_all["sec"], 
         "_".join([name, "all_5utr_length.png"]), source, "5utr")
    if len(utr_strain["all"]) > 1:
        for strain in utr_strain["all"].keys():
            plot(utr_strain["all"][strain], utr_strain["pri"][strain],
                 utr_strain["sec"][strain], "_".join([strain, "5utr_length.png"]),
                 source, "5utr")

def compare_term(ta, terms, fuzzy):
    for term in terms:
        if ta.strand == term.strand:
            if term.strand == "+":
                if (math.fabs(ta.end - term.start) <= fuzzy) or \
                   ((ta.end >= term.start) and (ta.end <= term.end)):
                    return term
            else:
                if (math.fabs(ta.start - term.end) <= fuzzy) or \
                   ((ta.start >= term.start) and (ta.start <= term.end)):
                    return term

def get_3utr(ta, near_cds, utr_all, utr_strain, attributes, num_utr, out):
    if ta.strand == "+":
        start = near_cds.end
        end = ta.end
        length = ta.end - near_cds.end
        utr_all.append(length)
        utr_strain[ta.seq_id].append(length)
    else:
        start = ta.start
        end = near_cds.start
        length = near_cds.start - ta.start
        utr_all.append(length)
        utr_strain[ta.seq_id].append(length)
    attributes.append("=".join(["lenght", str(length)]))
    attributes.append("=".join(["associated_tran", ta.attributes["ID"]]))
    attribute = ";".join(attributes)
    if (length <= 300) and (length > 0):
        name_utr='%0*d' % (5, num_utr)
        name = "=".join(["Name", "_".join(["3'UTR", ta.attributes["ID"]])])
        id_ = "ID=utr3_" + str(num_utr)
        num_utr += 1
        attribute_string = ";".join([id_, name, attribute])
        out.write("\t".join([ta.seq_id, "Transcript", "3UTR", str(start), str(end),
                             ta.score, ta.strand, ta.phase, attribute_string]) + "\n")
    return num_utr

def get_near_cds(cdss, genes, ta, attributes):
    first = True
    detect = False
    for cds in cdss:
        if (ta.seq_id == cds.seq_id) and \
           (ta.strand == cds.strand):
            if first:
                near_cds = cds
                first = False
            if ta.strand == "+":
                if (cds.end <= ta.end) and \
                   (cds.start >= ta.start) and \
                   (cds.end > near_cds.end):
                    near_cds = cds
                    detect = True
            else:
                if (cds.start >= ta.start) and \
                   (cds.end <= ta.end):
                    near_cds = cds
                    detect = True
                    break
    if detect:
        for gene in genes:
            if ("Parent" in near_cds.attributes.keys()):
                if near_cds.attributes["Parent"] == gene.attributes["ID"]:
                    attributes.append("=".join(["associated_gene", gene.attributes["locus_tag"]]))
        if "protein_id" in near_cds.attributes.keys():
            attributes.append("=".join(["associated_cds", near_cds.attributes["protein_id"]]))
        elif "locus_tag" in near_cds.attributes.keys():
            attributes.append("=".join(["associated_cds", near_cds.attributes["locus_tag"]]))
        else:
            attributes.append("=".join(["associated_cds", "_".join([near_cds.feature,
                              str(near_cds.start), str(near_cds.end), near_cds.strand])]))
    return near_cds

def Detect_3UTR(TA_file, GFF_file, term_file, fuzzy, out_file):
    tas = []
    cdss = []
    genes = []
    terms = []
    num_utr = 0
    utr_all = []
    utr_strain = {}
    pre_seq_id = ""
    out = open(out_file, "w")
    out.write("##gff-version 3\n")
    read_file(False, False, GFF_file, cdss, genes, TA_file, tas, term_file, terms)
    for ta in tas:
        if term_file is not False:
            term = compare_term(ta, terms, fuzzy)
        else:
            term = None
        attributes = []
        if term is not None:
            attributes.append("=".join(["associated_term", term.attributes["ID"]]))
        near_cds = get_near_cds(cdss, genes, ta, attributes)
        if ta.seq_id != pre_seq_id:
            pre_seq_id = ta.seq_id
            utr_strain[ta.seq_id] = []
        num_utr = get_3utr(ta, near_cds, utr_all, utr_strain, attributes, num_utr, out)
    name = (GFF_file.split("/"))[-1].replace(".gff", "")
    plot(utr_all, None, None, "_".join([name, "all_3utr_length.png"]), None, "3utr")
    if len(utr_strain) > 1:
        for strain in utr_strain.keys():
            plot(utr_strain[strain], None, None, "_".join([strain, "3utr_length.png"]), None, "3utr")
