import shutil
import math
import csv
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.helper import Helper


def read_file(seq_file, srna_table):
    seq = ""
    with open(seq_file) as fh:
        for line in fh:
            if not line.startswith(">"):
                line = line.strip()
                seq = seq + line
    tabs = []
    sh = open(srna_table, "r")
    for row in csv.reader(sh, delimiter='\t'):
        tabs.append({"info": row, "seq_id": row[0],
                     "start": int(row[2]), "end": int(row[3]),
                     "strand": row[4]})
    return seq, tabs

def get_table_entry(tabs, srna):
    for tab in tabs:
        if (srna.seq_id == tab["seq_id"]) and (
                srna.start == tab["start"]) and (
                srna.end == tab["end"]) and (
                srna.strand == tab["strand"]):
            return tab

def backward_t(seq, start, end, strand, mut_u):
    no_ut = 0
    ext = 0
    bord = end - start
    while 1:
        if strand == "+":
            nt = Helper().extract_gene(seq, end - ext, end - ext, strand)
        else:
            nt = Helper().extract_gene(seq, start + ext,
                                       start + ext, strand)
        if (nt == "U") or (nt == "T"):
            pass
        else:
            no_ut += 1
            if no_ut > mut_u:
                break
        ext += 1
        if ext >= bord:
            break
    return ext    

def forward_t(seq, start, end, strand, mut_u):
    no_ut = 0
    ext = 0
    bord = end - start
    while 1:
        if strand == "+":
            nt = Helper().extract_gene(seq, end + ext, end + ext, strand)
        else:
            nt = Helper().extract_gene(seq, start - ext,
                                       start - ext, strand)
        if (nt == "U") or (nt == "T"):
            pass
        else:
            no_ut += 1
            if no_ut > mut_u:
                break
        ext += 1
        if ext >= bord:
            break
    return ext

def iterate_seq(seq_u, args_srna):
    pos = 0
    first = True
    nts = {"ut": 0, "no_ut": 0}
    while 1:
        if (len(seq_u) - pos) < args_srna.num_u:
            break
        for nt in reversed(seq_u[:(len(seq_u) - pos)]):
            if first:
                first = False
                if (nt != "U") and (nt != "T"):
                    break
            if (nt == "U") or (nt == "T"):
                nts["ut"] += 1
            else:
                nts["no_ut"] += 1
                if nts["no_ut"] > args_srna.mut_u:
                    break
        if nts["ut"] < args_srna.num_u:
            nts = {"ut": 0, "no_ut": 0}
            first = True
        else:
            break
        pos += 1
    return pos


def search_t(seq, start, end, strand, ext_b, ext_f, args_srna):
    if strand == "+":
        seq_end = end + args_srna.len_u + ext_f + 1
        if seq_end > len(seq):
            seq_end = len(seq)
        seq_u = Helper().extract_gene(
                    seq, end - ext_b - 1, seq_end, strand)
    else:
        seq_start = start - args_srna.len_u - ext_f - 1
        if (seq_start) < 1:
            seq_start = 1
        seq_u = Helper().extract_gene(seq, seq_start,
                                      start + ext_b + 1, strand)
    pos = iterate_seq(seq_u, args_srna)
    if strand == "+":
        final_end = (seq_end - pos)
        if (final_end - end) <= 0:
            final_end = end
        elif (final_end - start) >= args_srna.max_len:
            diff = final_end - start - args_srna.max_len
            pos = iterate_seq(seq_u[:(final_end - diff)], args_srna)
            final_end = (final_end - diff - pos)
            if (final_end - end) <= 0:
                final_end = end
        final_start = start
    else:
        final_start = (seq_start + pos)
        if (start - final_start) <= 0:
            final_start = start
        elif (end - final_start) >= args_srna.max_len:
            diff = end - final_start - args_srna.max_len
            pos = iterate_seq(seq_u[(final_start + diff):], args_srna)
            final_start = (final_start + diff + pos)
            if (start - final_start) <= 0:
                final_start = start
        final_end = end
    return final_start, final_end

def check_term(srna, tab, seq, len_u, out, out_t):
    if "with_term" in srna.attributes.keys():
        feature = srna.attributes["with_term"].split(":")[0]
        info = srna.attributes["with_term"].split(":")[-1]
        start = int(info.split("-")[0])
        end = int(info.split("-")[-1].split("_")[0])
        strand = info.split("_")[-1]

def get_srna_poly_u(srna_file, seq_file, srna_table, args_srna):
    seq, tabs = read_file(seq_file, srna_table)
    out = open(srna_file + "_ext_polyu", "w")
    out_t = open(srna_table + "_ext_polyu", "w")
    gff_f = open(srna_file, "r")
    for entry in Gff3Parser().entries(gff_f):
        tab = get_table_entry(tabs, entry)
        ext_b = backward_t(seq, entry.start, entry.end,
                           entry.strand, args_srna.mut_u)
        ext_f = forward_t(seq, entry.start - args_srna.len_u,
                          entry.end + args_srna.len_u, entry.strand,
                          args_srna.mut_u)
        final_start, final_end = search_t(seq, entry.start, entry.end,
                                          entry.strand, ext_b, ext_f,
                                          args_srna)
        tab["info"][2] = str(final_start)
        tab["info"][3] = str(final_end)
        info_without_attributes = "\t".join([str(field) for field in [
                    entry.seq_id, entry.source, entry.feature, str(final_start),
                    str(final_end), entry.score, entry.strand, entry.phase]])
        out.write("\t".join([info_without_attributes,
                             entry.attribute_string]) + "\n")
        out_t.write("\t".join(tab["info"]) + "\n")
    gff_f.close()
    out.close()
    out_t.close()
    shutil.move(srna_file + "_ext_polyu", srna_file)
    shutil.move(srna_table + "_ext_polyu", srna_table)
