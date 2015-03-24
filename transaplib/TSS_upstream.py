#!/usr/bin/python

import os        
import sys
import csv
import math
mport transaplib.parser_wig as par_wig
from transaplib.gff3 import Gff3Parser
from transaplib.TSSpredator import TSSPredatorReader
from transaplib.helper import extract_gene


def get_upstream(seq, tss, out, name):
    if tss.strand == "+":
        fasta = extract_seq.extract_gene(seq, tss.start - 50, tss.start, tss.strand)
    else:
        fasta = extract_seq.extract_gene(seq, tss.start, tss.start + 50, tss.strand)
    out.write("%s\n%s\n" % (name, fasta))

def get_primary_locus_tag(tss):
    tsss = []
    tss_types = tss.attributes["type"].split(" ")
    tss_locus_tags = tss.attributes["associated_gene"].split(" ")
    tss_utr_lengths = tss.attributes["UTR_length"].split(" ")
    index = 0
    for tss_type in tss_types:
        if "Primary" in tss_type:
            tsss.append({"locus": tss_locus_tags[index],
                         "utr": int(tss_utr_lengths[index].split("_")[1]),
                         "type": tss_type})
        index += 1
    return tsss

def import_to_tss(tss_type, cds_pos, tss, locus_tag, tss_entry):
    if cds_pos == "NA":
        utr = tss_type + "_" + "NA"
    else:
        utr = tss_type + "_" + str(int(math.fabs(cds_pos - tss.start)))
    if len(tss_entry) != 0:
        tss_dict = tss_entry[1]
        tss_dict_types = tss_dict["type"].split(" ")
        tss_dict_utrs = tss_dict["UTR_length"].split(" ")
        tss_dict_tags = tss_dict["associated_gene"].split(" ")
        if tss_type == "Primary" and ("Primary" in tss_dict["type"]):
            index = 0
            for tss_dict_type in tss_dict_types:
                if "Primary" in tss_dict_type:
                    utr_length = tss_dict_utrs[index].split("_")
                    if math.fabs(cds_pos - tss.start) < int(utr_length[1]):
                        tss_dict_utrs[index] = utr
                        tss_dict_tags[index] = locus_tag
                index += 1
        else:
            tss_dict_types.append(tss_type)
            tss_dict_utrs.append(utr)
            tss_dict_tags.append(locus_tag)
        tss_dict = {"Name": "TSS_" + str(tss.start) + tss.strand,
                    "type": " ".join(tss_dict_types),
                    "UTR_length": " ".join(tss_dict_utrs),
                    "associated_gene": " ".join(tss_dict_tags)}
    else:
        tss_dict = {"Name": "TSS_" + str(tss.start) + tss.strand,
                    "type": tss_type,
                    "UTR_length": utr,
                    "associated_gene": locus_tag}
    tss_string = ";".join(["=".join(["UTR_length", tss_dict["UTR_length"]]),
                           "=".join(["associated_gene", tss_dict["associated_gene"]]),
                           "=".join(["type", tss_dict["type"]]),
                           "=".join(["Name", tss_dict["Name"]])])
    return (tss_string, tss_dict)

def detect_coverage(wigs, tss, ref):
    for strain, tracks in wigs.items():
        if strain == tss.seq_id:
            tss_cover = 0
            ref_cover = 0
            for track, wig in tracks.items():
                if tss.start <= len(wig):
                    if tss.strand == "+":
                        diff_t = (wig[tss.start - 1]["coverage"] - \
                                wig[tss.start - 2]["coverage"])
                        diff_r = (wig[ref.start - 1]["coverage"] - \
                                wig[ref.start - 2]["coverage"])
                    else:
                        diff_t = (wig[tss.start - 1]["coverage"] - \
                                wig[tss.start]["coverage"])
                        diff_r = (wig[ref.start - 1]["coverage"] - \
                                wig[ref.start]["coverage"])
                    tss_cover = tss_cover + diff_t
                    ref_cover = ref_cover + diff_r
    return (tss_cover, ref_cover)

def fix_attributes(tss, tss_entry):
    index = 0
    genes = tss.attributes["associated_gene"].split(" ")
    utrs = tss.attributes["UTR_length"].split(" ")
    types = tss.attributes["type"].split(" ")
    for gene in genes:
        if gene == tss_entry["locus"]:
            utrs[index] = utrs[index].replace("Primary", "Secondary")
            types[index] = types[index].replace("Primary", "Secondary")
        index += 1
    tss.attributes["UTR_length"] = " ".join(utrs)
    tss.attributes["type"] = " ".join(types)

def get_real_pri_sec(type_, detects, indexs, utrs, real_utrs):
    if (type_ == "Primary") and (detects["pri"] is False):
        detects["pri"] = True
        real_utrs["pri"] = int(utrs[indexs["index_num"]].split("_")[1])
        indexs["real_pri"] = indexs["index_num"]
    elif (type_ == "Primary") and (detects["pri"] is True):
        compare_utr = int(utrs[indexs["index_num"]].split("_")[1])
        if compare_utr < real_utrs["pri"]:
            real_utrs["pri"] = compare_utr
            indexs["real_pri"] = indexs["index_num"]
    elif (type_ == "Secondary") and (detects["sec"] is False):
        detects["sec"] = True
        real_utrs["sec"] = int(utrs[indexs["index_num"]].split("_")[1])
        indexs["real_sec"] = indexs["index_num"]
    elif (type_ == "Secondary") and (detects["sec"] is True):
        compare_utr = int(utrs[indexs["index_num"]].split("_")[1])
        if compare_utr < real_utrs["sec"]:
            real_utrs["sec"] = compare_utr
            indexs["real_sec"] = indexs["index_num"]

def del_repeat(tsss):
    for tss in tsss:
        types = tss.attributes["type"].split(" ")
        utrs = tss.attributes["UTR_length"].split(" ")
        genes = tss.attributes["associated_gene"].split(" ")
        detects = {"pri": False, "sec": False}
        indexs = {"index_num": 0, "real_pri": 0, "real_sec": 0}
        final_types = []
        final_utrs = []
        final_genes = []
        real_utrs = {"pri": -1, "sec": -1}
        for type_ in types:
            get_real_pri_sec(type_, detects, indexs, utrs, real_utrs):
            elif (type_ == "Antisense") or \
                 (type_ == "Internal") or \
                 (type_ == "Orphan"):
                final_types.append(types[index])
                final_utrs.append(utrs[index])
                final_genes.append(genes[index])
            index += 1
        if detect_pri is True:    
            final_types.append(types[real_index])
            final_utrs.append(utrs[real_index])
            final_genes.append(genes[real_index])
        else:
            if detect_sec is True:
                final_types.append(types[real_index2])
                final_utrs.append(utrs[real_index2])
                final_genes.append(genes[real_index2])
        tss.attributes["type"] = " ".join(final_types)
        tss.attributes["UTR_length"] = " ".join(final_utrs)
        tss.attributes["associated_gene"] = " ".join(final_genes)

def fix_primary_type(tsss, wigs_f, wigs_r):
    num_man = 0
    for tss in tsss:
        num_auto = 0
        if ("Primary" in tss.attributes["type"]):
            tss_entrys = get_primary_locus_tag(tss)
            for ref in tsss:
                if (ref.seq_id == tss.seq_id) and \
                   (ref.strand == tss.strand) and \
                   (ref.start == tss.start):
                    pass
                else:
                    if ("Primary" in ref.attributes["type"]):
                        ref_entrys = get_primary_locus_tag(ref)
                        for tss_entry in tss_entrys:
                            for ref_entry in ref_entrys:
                                if (tss_entry["locus"] == ref_entry["locus"]) and \
                                   (tss_entry["type"] == "Primary") and \
                                   (ref_entry["type"] == "Primary") and \
                                   (tss.seq_id == ref.seq_id):
                                    if tss.strand == "+":
                                        covers = detect_coverage(wigs_f, tss, ref)
                                    else:
                                        covers = detect_coverage(wigs_r, tss, ref)
                                    tss_cover = covers[0]
                                    ref_cover = covers[1]
                                    if tss_cover < ref_cover:
                                        fix_attributes(tss, tss_entry)
                                    elif tss_cover > ref_cover:
                                        fix_attributes(ref, ref_entry)
                                    elif tss_cover == ref_cover:
                                        if (tss_entry["utr"] < ref_entry["utr"]):
                                            fix_attributes(ref, ref_entry)
                                        elif (tss_entry["utr"] > ref_entry["utr"]):
                                            fix_attributes(tss, tss_entry)
    del_repeat(tsss)
    return tsss

def remove_primary(tss, tss_entry):
    final_types = []
    final_utrs = []
    final_genes = []
    tss_dict = tss_entry[1]
    types = tss_dict["type"].split(" ")
    utrs = tss_dict["UTR_length"].split(" ")
    genes = tss_dict["associated_gene"].split(" ")
    index = 0
    for type_ in types:
        if type_ != "Primary":
            final_types.append(type_)
            final_utrs.append(utrs[index])
            final_genes.append(genes[index])
        index += 1
    tss_dict = {"Name": "TSS_" + str(tss.start) + tss.strand,
                "type": " ".join(final_types),
                "UTR_length": " ".join(final_utrs),
                "associated_gene": " ".join(final_genes)}
    tss_string = ";".join(["=".join(["UTR_length", tss_dict["UTR_length"]]),
                           "=".join(["associated_gene", tss_dict["associated_gene"]]),
                           "=".join(["type", tss_dict["type"]]),
                           "=".join(["Name", tss_dict["Name"]])])
    return [tss_string, tss_dict]

def check_pri_inter(gene, tss, anti_ends, checks, gene_ends, tss_entry):
    if _is_primary(gene.start, gene.end, tss.start, tss.strand):
        locus_tag = gene.attributes["locus_tag"]
        if tss.strand == "+":
            if ((anti_ends["reverse"] != -1) and (anti_ends["reverse"] - gene.start) > 0) or \
               (anti_ends["reverse"] == -1):
                tss_entry = import_to_tss("Primary", gene.start, tss, locus_tag, tss_entry)
                checks["orphan"] = False
                gene_ends["forward"] = gene.start
            elif (anti_ends["reverse"] != -1) and \
                 ((anti_ends["reverse"] - gene.start) < 0):
                if (checks["int_anti"] is True) or (tss.start - anti_ends["reverse"]) > 0:
                    tss_entry = import_to_tss("Primary", gene.start, tss, locus_tag, tss_entry)
                    checks["orphan"] = False
                    gene_ends["forward"] = gene.start
        else:
            if ((anti_ends["forward"] != -1) and (gene.end - anti_ends["forward"]) > 0) or \
               (anti_ends["forward"] == -1):
                tss_entry = import_to_tss("Primary", gene.end, tss, locus_tag, tss_entry)
                checks["orphan"] = False
                gene_ends["reverse"] = gene.end
    if _is_internal(gene.start, gene.end, tss.start, tss.strand):
        locus_tag = gene.attributes["locus_tag"]
        tss_entry = import_to_tss("Internal", "NA", tss, locus_tag, tss_entry)
        checks["orphan"] = False
    return tss_entry

def check_anti_enter(gene, tss, anti_ends, checks, gene_ends, tss_entry):
    if _is_antisense(gene.start, gene.end, tss.start, tss.strand):
        checks["int_anti"] = False
        if tss.strand == "-":
            anti_ends["forward"] = gene.start
            if (gene_ends["reverse"] != -1) and (gene.start - gene_ends["reverse"]) > 0:
                if _is_internal(gene.start, gene.end, tss.start, tss.strand):
                    pass
                else:
                    if (tss.start - gene.end) > 0:
                        tss_entry = remove_primary(tss, tss_entry)
        else:
            anti_ends["reverse"] = gene.end
            if _is_internal(gene.start, gene.end, tss.start, tss.strand):
                checks["int_anti"] = True
            if (gene_ends["forward"] != -1) and (gene.start - gene_ends["forward"]) > 0:
                if (checks["int_anti"] is not True) and \
                   (gene.start - tss.start) > 0:
                    tss_entry = remove_primary(tss, tss_entry)
        locus_tag = gene.attributes["locus_tag"]
        tss_entry = import_to_tss("Antisense", "NA", tss, locus_tag, tss_entry)
        checks["orphan"] = False

def compare_tss_cds(tss, cdss, genes):
    tss_entry = []
    checks = {"operon": True, "inter_anti": None}
    gene_ends = {"forward": -1, "reverse": -1}
    anti_ends = {"forward": -1, "reverse": -1}
    for gene in genes:
        if gene.strand == tss.strand:
            tss_entry = check_pri_inter(gene, tss, anti_ends, checks, gene_ends, tss_entry)
        else:
            tss_entry = check_anti_inter(gene, tss, anti_ends, checks, gene_ends, tss_entry)
    if checks["orphan"] == True:
        tss_entry = import_to_tss("Orphan", "NA", tss, "NA", tss_entry)
    return tss_entry

def _is_primary(cds_start, cds_end, tss_pos, strand):
    if strand == "+":
        if (_is_utr(cds_start, tss_pos, 300) and (cds_start >= tss_pos)):
            return True
    else:
        if (_is_utr(tss_pos, cds_end, 300) and (cds_end <= tss_pos)):
            return True

def _is_internal(cds_start, cds_end, tss_pos, strand):
    if ((cds_start < tss_pos) and (cds_end > tss_pos)) or \
       ((strand == "+") and (tss_pos == cds_end)) or \
       ((strand == "-") and (tss_pos == cds_start)):
        return True

def _is_antisense(cds_start, cds_end, tss_pos, strand):
    if ((_is_utr(cds_start, tss_pos, 100)) and (cds_start >= tss_pos)) or \
       ((_is_utr(tss_pos, cds_end, 100)) and (cds_end <= tss_pos)) or \
       (_is_internal(cds_start, cds_end, tss_pos, strand)):
        return True

def _is_utr(pos1, pos2, length):
    if (pos1 - pos2 <= length):
        return True

def print_fasta(seq, tss, pri, sec, inter, anti, orph, name):
    for key in seq.keys():
        if tss.seq_id == key:
            if "Primary" in tss.attributes["type"]:
                get_upstream(seq[key], tss, pri, name)
            if "Secondary" in tss.attributes["type"]:
                get_upstream(seq[key], tss, sec, name)
            if "Internal" in tss.attributes["type"]:
                get_upstream(seq[key], tss, inter, name)
            if "Antisense" in tss.attributes["type"]:
                get_upstream(seq[key], tss, anti, name)
            if "Orphan" in tss.attributes["type"]:
                get_upstream(seq[key], tss, orph, name)

def read_wig(wigs, filename, strand):
    wig_parser = par_wig.parser_wig()
    if filename is not False:
        for entry in wig_parser.parser(filename, strand):
            if entry.strain not in wigs.keys():
                strain = entry.strain
                wigs[strain] = {}
            if entry.track not in wigs[strain].keys():
                wigs[strain][entry.track] = []
            wigs[strain][entry.track].append({
                 "pos": entry.pos, "coverage": entry.coverage,
                 "strand": entry.strand})

def read_data(tsss, seq, TSS_file, fasta_file):
    t_f = open(TSS_file, "r")
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

def read_gff(cdss, genes):
    gff_parser = gff3.Gff3Parser()
    g_f = open(args.gff_file, "r")
    for entry in gff_parser.entries(g_f):
        if (entry.feature == "CDS") or \
           (entry.feature == "rRNA") or \
           (entry.feature == "tRNA"):
            cdss.append(entry)
        if entry.feature == "gene":
            genes.append(entry)

def read_libs(input_libs, wig_folder):
    if "merge_forward.wig" in os.listdir(os.getcwd()):
        os.system("rm merge_forward.wig")
    if "merge_reverse.wig" in os.listdir(os.getcwd()):
        os.system("rm merge_reverse.wig")
    for lib in input_libs:
        datas = lib.split(":")
        if (datas[1] == "tex") and (datas[4] == "+"):
            os.system("cat " + wig_folder + "/" + datas[0] + " >> merge_forward.wig")
        elif (datas[1] == "tex") and (datas[4] == "-"):
            os.system("cat " + wig_folder + "/" + datas[0] + " >> merge_reverse.wig")

def TSS_upstream(TSS_file, fasta_file, gff_file, source, input_libs, wig_folder, no_class_output):
    pri = open("primary.fa", "w")
    sec = open("secondary.fa", "w")
    inter = open("internal.fa", "w")
    anti = open("antisense.fa", "w")
    orph = open("orphan.fa", "w")
    tsss = []
    seq = {}
    cdss = []
    genes = []
    wigs_f = {}
    wigs_r = {}
    read_data(tsss, seq, TSS_file, fasta_file)
    num_tss = 0
    if source is False:
        out = open(no_class_output)
        out.write("##gff-version 3\n")
        read_gff(cdss, genes)
        cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start))
        genes = sorted(genes, key=lambda k: (k.seq_id, k.start))
    for tss in tsss:
        if source is True:
            name = ">" + "_".join([str(tss.start), tss.strand, tss.seq_id])
            print_fasta(seq, tss, pri, sec, inter, anti, orph, name)
        else:
            tss_type = compare_tss_cds(tss, cdss, genes)
            tss.attributes = tss_type[1]
            tss.attribute_string = tss_type[0] + ";ID=tss" + str(num_tss)
            num_tss += 1
    if source is False:
        read_libs(input_libs, wig_folder)
        read_wig(wigs_f, "merge_forward.wig", "+")
        read_wig(wigs_r, "merge_reverse.wig", "-")
        sort_tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start))
        final_tsss = fix_primary_type(sort_tsss, wigs_f, wigs_r)
        for tss in final_tsss:
            name = ">" + "_".join([str(tss.start), tss.strand, tss.seq_id])
            tss.attribute_string = ";".join(
                ["=".join(items) for items in tss.attributes.items()])
            out.write("\t".join([str(field) for field in [
                            tss.seq_id, tss.source, tss.feature, tss.start,
                            tss.end, tss.score, tss.strand, tss.phase,
                            tss.attribute_string]]) + "\n")
            print_fasta(seq, tss, pri, sec, inter, anti, orph, name)
        os.system("rm merge_forward.wig")
        os.system("rm merge_reverse.wig")           
