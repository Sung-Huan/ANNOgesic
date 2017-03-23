import os
import csv
import shutil
from annogesiclib.gff3 import Gff3Parser


def read_file(gff_file, args_srna):
    srnas = []
    for entry in Gff3Parser().entries(open(gff_file)):
        attributes = {}
        for key, value in entry.attributes.items():
            if "promoter" not in key:
                attributes[key] = value
        entry.attributes = attributes
        srnas.append(entry)
    srnas = sorted(srnas, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    fh = open(args_srna.promoter_table, "r")
    pros = []
    for row in csv.reader(fh, delimiter='\t'):
        if (row[0] != "Genome") and (
                row[3] in args_srna.promoter_name):
            pros.append({"strain": row[0], "pos": row[1],
                         "strand": row[2], "name": row[3]})
    fh.close()
    return srnas, pros


def print_table(srna_table, out_t, srnas):
    fh = open(srna_table, "r")
    for row in csv.reader(fh, delimiter='\t'):
        for srna in srnas:
            if (row[0] == srna.seq_id) and (
                    int(row[2]) == srna.start) and (
                    int(row[3]) == srna.end) and (
                    row[4] == srna.strand):
                if "promoter" in srna.attributes.keys():
                    promoter = [srna.attributes["promoter"]]
                else:
                    promoter = ["NA"]
                out_t.write("\t".join(row + promoter) + "\n")


def compare_srna_promoter(srna_gff, srna_table, args_srna):
    '''compare sRNA and promoter to find the sRNA 
    which is associated with a promoter.
    it is for the ranking of sRNA'''
    srnas, pros = read_file(srna_gff, args_srna)
    out_g = open("tmp_srna.gff", "w")
    out_t = open("tmp_srna.csv", "w")
    out_g.write("##gff-version 3\n")
    for srna in srnas:
        tsss = []
        detect = False
        if "with_TSS" in srna.attributes.keys():
            if srna.attributes["with_TSS"] != "NA":
                datas = srna.attributes["with_TSS"].split(",")
                for data in datas:
                    info = data.split(":")[-1]
                    tss = info.split("_")
                    tsss.append({"pos": tss[0], "strand": tss[-1]})
        if len(tsss) != 0:
            for tss in tsss:
                for pro in pros:
                    if (srna.seq_id == pro["strain"]) and (
                            tss["strand"] == pro["strand"]) and (
                            tss["pos"] == pro["pos"]):
                        detect = True
                        if "promoter" not in srna.attributes.keys():
                            srna.attributes["promoter"] = pro["name"]
                        else:
                            srna.attributes["promoter"] = ",".join([
                                srna.attributes["promoter"],
                                pro["name"]])
        if detect:
            out_g.write(srna.info + ";promoter=" +
                        srna.attributes["promoter"] + "\n")
        else:
            out_g.write(srna.info + ";promoter=NA" + "\n")
    print_table(srna_table, out_t, srnas)
    os.remove(srna_gff)
    os.remove(srna_table)
    out_t.close()
    out_g.close()
    shutil.move("tmp_srna.gff", srna_gff)
    shutil.move("tmp_srna.csv", srna_table)
