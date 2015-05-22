import os
import sys
from annogesiclib.helper import Helper
from annogesiclib.gff3 import Gff3Parser

def print_fasta(entry, seq, out):
    try:
        if entry.attributes["locus_tag"] == entry.attributes["Name"]:
            out.write(">{0}\n{1}\n".format(
                      entry.attributes["locus_tag"], seq))
        else:
            out.write(">" + "|".join([entry.attributes["locus_tag"], entry.attributes["Name"]]) + "\n")
            out.write(seq + "\n")
    except KeyError:
        out.write(">{0}:{1}-{2}_{3}\n{4}\n".format(
                  entry.feature, entry.start, entry.end, entry.strand, seq))

def read_file(seq_file, gff_file, cdss_f, cdss_r, genes):
    fastas = []
    with open (seq_file, "r") as seq_f:
        for line in seq_f:
            if line.startswith(">"):
                continue
            else:
                line=line.strip()
                fastas.append(line)
    fasta = "".join(fastas)
    for entry in Gff3Parser().entries(open(gff_file)):
        if (entry.feature == "CDS") and (entry.strand == "+"):
            cdss_f.append(entry)
        elif (entry.feature == "CDS") and (entry.strand == "-"):
            cdss_r.append(entry)
        elif entry.feature == "gene":
            genes.append(entry)
    return fasta

def deal_cds_forward(cdss_f, target_folder, fasta, genes):
    pre_id = ""
    for cds in cdss_f:
        if cds.seq_id != pre_id:
            out = open(os.path.join(target_folder, 
                       "_".join([cds.seq_id, "target.fa"])), "w")
            pre_id = cds.seq_id
        if (cds.start > 350):
            start = cds.start - 350
        else:
            start = 1
        if ((cds.start + 300) < len(fasta)) and ((cds.end - cds.start) >= 300):
            end = cds.start + 299
        elif cds.start + 300 >= len(fasta):
            end = len(fasta)
        elif (cds.end - cds.start) < 300:
            end = cds.end
        pre_cds_f = cds
        seq = Helper().extract_gene(fasta, start, end, cds.strand)
        taraget = cds
        for gene in genes:
            if "Parent" in cds.attributes.keys():
                if cds.attributes["Parent"] == gene.attributes["ID"]:
                    target = gene
                    break
        print_fasta(target, seq, out)

def deal_cds_reverse(cdss_r, target_folder, fasta, genes):
    pre_id = ""
    for cds in cdss_r:
        if cds.seq_id != pre_id:
            out = open(os.path.join(target_folder, 
                       "_".join([cds.seq_id, "target.fa"])), "a")
            pre_id = cds.seq_id
        if (len(fasta) - cds.end > 350):
            end = cds.end + 350
        else:
            end = len(fasta)
        if ((cds.end - 300) > 1) and ((cds.end - cds.start) >= 300):
            start = cds.end - 299
        elif cds.end - 300 < 1:
            start = 1
        elif (cds.end - cds.start) < 300:
            start = cds.start
        pre_cds_r = cds
        seq = Helper().extract_gene(fasta, start, end, cds.strand)
        taraget = cds
        for gene in genes:
            if "Parent" in cds.attributes.keys():
                if cds.attributes["Parent"] == gene.attributes["ID"]:
                    target = gene
                    break
        print_fasta(target, seq, out)

def potential_target(gff_file, seq_file, target_folder):
    cdss_f = []
    cdss_r = []
    genes = []
    fasta = read_file(seq_file, gff_file, cdss_f, cdss_r, genes)
    sort_cdss_f = sorted(cdss_f, key=lambda k: (k.seq_id, k.start))
    deal_cds_forward(sort_cdss_f, target_folder, fasta, genes)
    sort_cdss_r = sorted(cdss_r, reverse=True, key=lambda k: (k.seq_id, k.start))
    deal_cds_reverse(sort_cdss_r, target_folder, fasta, genes)
