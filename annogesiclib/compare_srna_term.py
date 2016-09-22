import os
import csv
import shutil
from annogesiclib.gff3 import Gff3Parser


def read_gff(gff_file):
    datas = []
    for entry in Gff3Parser().entries(open(gff_file)):
        datas.append(entry)
    datas = sorted(datas, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    return datas


def print_table(srna_table, out_t, srnas):
    fh = open(srna_table, "r")
    for row in csv.reader(fh, delimiter='\t'):
        for srna in srnas:
            if (row[0] == srna.seq_id) and (
                int(row[2]) == srna.start) and (
                int(row[3]) == srna.end) and (
                    row[4] == srna.strand):
                if "with_term" in srna.attributes.keys():
                    with_term = [srna.attributes["with_term"]]
                else:
                    with_term = ["NA"]
                out_t.write("\t".join(row + with_term) + "\n")


def compare_srna_term(srna_gff, srna_table, term_file, fuzzy_b, fuzzy_a):
    '''Comparison of sRNA and terminator. 
    It can search the sRNA which is associated with terminator'''
    srnas = read_gff(srna_gff)
    terms = read_gff(term_file)
    out_g = open("tmp_srna.gff", "w")
    out_t = open("tmp_srna.csv", "w")
    out_g.write("##gff-version 3\n")
    for srna in srnas:
        detect = False
        for term in terms:
            if (srna.seq_id == term.seq_id) and (
                    srna.strand == term.strand):
                if (srna.strand == "+"):
                    if (((srna.end - term.end) <= fuzzy_b) and (
                            srna.end >= term.end) and (
                            srna.start < term.start)) or (
                            ((term.start - srna.end) <= fuzzy_a) and (
                             term.start >= srna.end)) or (
                            (srna.end > term.start) and (
                             srna.end < term.end) and (
                             srna.start < term.start)):
                        term_string = (term.feature + ":" + str(term.start) +
                                       "-" + str(term.end) + "_" + term.strand)
                        srna.attributes["with_term"] = term_string
                        detect = True
                        break
                else:
                    if (((term.start - srna.start) <= fuzzy_b) and (
                            term.start >= srna.start) and (
                            term.end < srna.end)) or (
                            ((srna.start - term.end) <= fuzzy_a) and (
                             srna.start >= term.end)) or (
                            (srna.start > term.start) and (
                             srna.start < term.end) and (
                             srna.end > term.end)):
                        term_string = (term.feature + ":" + str(term.start) +
                                       "-" + str(term.end) + "_" + term.strand)
                        srna.attributes["with_term"] = term_string
                        detect = True
                        break
        if "end_cleavage" in srna.attributes.keys():
            if (srna.attributes["end_cleavage"] != "NA") and (
                    "with_term" not in srna.attributes.keys()):
                srna.attributes["with_term"] = srna.attributes["end_cleavage"]
            elif (srna.attributes["end_cleavage"] != "NA") and (
                    "with_term" in srna.attributes.keys()):
                srna.attributes["with_term"] = ",".join([
                    srna.attributes["with_term"],
                    srna.attributes["end_cleavage"]])
        if detect:
            out_g.write(srna.info + ";with_term=" +
                        srna.attributes["with_term"] + "\n")
        else:
            out_g.write(srna.info + ";with_term=NA" + "\n")
    print_table(srna_table, out_t, srnas)
    os.remove(srna_gff)
    os.remove(srna_table)
    out_t.close()
    out_g.close()
    shutil.move("tmp_srna.gff", srna_gff)
    shutil.move("tmp_srna.csv", srna_table)
