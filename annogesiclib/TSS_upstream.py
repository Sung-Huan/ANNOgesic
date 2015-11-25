import os
import sys
import csv
import math
import copy
import shutil
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.TSSpredator import TSSPredatorReader
from annogesiclib.helper import Helper
from annogesiclib.parser_wig import WigParser
from annogesiclib.gen_TSS_type import compare_tss_cds, fix_primary_type

def get_upstream(seq, tss, out, name):
    if tss.strand == "+":
        fasta = Helper().extract_gene(seq, tss.start - 50,
                                      tss.start, tss.strand)
    else:
        fasta = Helper().extract_gene(seq, tss.start,
                                      tss.start + 50, tss.strand)
    out.write("{0}\n{1}\n".format(name, fasta))

def print_fasta(seq, tss, files, name):
    for key in seq.keys():
        if tss.seq_id == key:
            if "Primary" in tss.attributes["type"]:
                get_upstream(seq[key], tss, files["pri"], name)
            if "Secondary" in tss.attributes["type"]:
                get_upstream(seq[key], tss, files["sec"], name)
            if "Internal" in tss.attributes["type"]:
                get_upstream(seq[key], tss, files["inter"], name)
            if "Antisense" in tss.attributes["type"]:
                get_upstream(seq[key], tss, files["anti"], name)
            if "Orphan" in tss.attributes["type"]:
                get_upstream(seq[key], tss, files["orph"], name)

def read_wig(filename, strand):
    wigs = {}
    wig_parser = WigParser()
    wig_fh = open(filename)
    if filename:
        for entry in wig_parser.parser(wig_fh, strand):
            if entry.strain not in wigs.keys():
                strain = entry.strain
                wigs[strain] = {}
            if entry.track not in wigs[strain].keys():
                wigs[strain][entry.track] = []
            wigs[strain][entry.track].append({
                 "pos": entry.pos, "coverage": entry.coverage,
                 "strand": entry.strand})
    wig_fh.close()
    return wigs

def read_data(tss_file, fasta_file):
    seq = {}
    tsss = []
    t_f = open(tss_file, "r")
    for entry in Gff3Parser().entries(t_f):
        tsss.append(entry)
    with open(fasta_file, "r") as f_h:
        for line in f_h:
            line = line.strip()
            if line[0] == ">":
                seq[line[1:]] = ""
                seq_id = line[1:]
            else:
                seq[seq_id] = seq[seq_id] + line
    tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start))
    return tsss, seq

def read_gff(gff_file):
    cdss = []
    genes = []
    g_f = open(gff_file, "r")
    for entry in Gff3Parser().entries(g_f):
        if (entry.feature == "CDS") or \
           (entry.feature == "rRNA") or \
           (entry.feature == "tRNA"):
            cdss.append(entry)
        if entry.feature == "gene":
            genes.append(entry)
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start))
    genes = sorted(genes, key=lambda k: (k.seq_id, k.start))
    return cdss, genes

def read_libs(input_libs, wig_folder):
    libs = {}
    if "merge_forward.wig" in os.listdir(os.path.join(os.getcwd(), "tmp")):
        os.remove("tmp/merge_forward.wig")
    if "merge_reverse.wig" in os.listdir(os.path.join(os.getcwd(), "tmp")):
        os.remove("tmp/merge_reverse.wig")
    for lib in input_libs:
        datas = lib.split(":")
        if (datas[1] == "tex") and (datas[4] == "+"):
            Helper().merge_file(os.path.join(wig_folder, datas[0]),
                                os.path.join("tmp", "merge_forward.wig"))
        elif (datas[1] == "tex") and (datas[4] == "-"):
            Helper().merge_file(os.path.join(wig_folder, datas[0]),
                                os.path.join("tmp", "merge_reverse.wig"))
    return libs

def upstream(tss_file, fasta_file, gff_file, source, wig_folder,
             input_libs, out_class):
    files = {"pri": open("tmp/primary.fa", "w"),
             "sec": open("tmp/secondary.fa", "w"),
             "inter": open("tmp/internal.fa", "w"),
             "anti": open("tmp/antisense.fa", "w"),
             "orph": open("tmp/orphan.fa", "w")}
    tsss, seq = read_data(tss_file, fasta_file)
    num_tss = 0
    if not source:
        out = open(out_class, "w")
        out.write("##gff-version 3\n")
        cdss, genes = read_gff(gff_file)
    for tss in tsss:
        if source is True:
            name = ">" + "_".join([tss.seq_id, str(tss.start), tss.strand])
            print_fasta(seq, tss, files, name)
        else:
            tss_type = compare_tss_cds(tss, cdss, genes)
            tss.attributes = tss_type[1]
            tss.attributes["ID"] = "tss" + str(num_tss)
            tss.attribute_string = "".join([tss_type[0], ";ID=tss", str(num_tss)])
            num_tss += 1
    if not source:
        libs = read_libs(input_libs, wig_folder)
        wigs_f = read_wig("tmp/merge_forward.wig", "+")
        wigs_r = read_wig("tmp/merge_reverse.wig", "-")
        sort_tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start))
        final_tsss = fix_primary_type(sort_tsss, wigs_f, wigs_r)
        for tss in final_tsss:
            name = ">" + "_".join([tss.seq_id, str(tss.start), tss.strand])
            tss.attribute_string = ";".join(
                ["=".join(items) for items in tss.attributes.items()])
            out.write("\t".join([str(field) for field in [
                            tss.seq_id, tss.source, tss.feature, tss.start,
                            tss.end, tss.score, tss.strand, tss.phase,
                            tss.attribute_string]]) + "\n")
            print_fasta(seq, tss, files, name)

def del_repeat_fasta(input_file, out_file):
    data = {}
    seq = ""
    check_same = False
    first_file = True
    pre_line = None
    out = open(out_file, "w")
    with open(input_file, "r") as f_h:
        for line in f_h:
            line = line.strip()
            if line[0] == ">":
                if check_same:
                    check_same = False
                if first_file:
                    seq_id = line[1:]
                    first_file = False
                    data[seq_id] = ""
                else:
                    if line[1:] in data.keys():
                        check_same = True
                    else:
                        seq_id = line[1:]
                        data[seq_id] = ""
            else:
                if check_same:
                    pass
                else:
                    data[seq_id] = data[seq_id] + line
    for strain, fasta in data.items():
        out.write(">" + strain + "\n")
        out.write(fasta + "\n")
