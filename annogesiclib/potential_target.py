import os
from annogesiclib.helper import Helper
from annogesiclib.gff3 import Gff3Parser


def assign_name(entry):
    if entry.attributes["locus_tag"] == entry.attributes["Name"]:
        return None
    else:
        return entry.attributes["Name"]

def print_fasta(entry, seq, out, gene):
    if gene is not None:
        if ("locus_tag" in gene.attributes.keys()):
            locus = gene.attributes["locus_tag"]
        else:
            locus = "NA"
    else:
        locus = "NA"
    if ("ID" in entry.attributes.keys()):
        out.write(">{0}_{1}-{2}_{3}\n{4}\n".format(
                  "_".join([locus,
                            entry.attributes["ID"]]),
                  entry.start, entry.end,
                  entry.strand, seq))
    else:
        out.write(">{0}_{1}-{2}_{3}\n{4}\n".format(
                  "_".join([locus, "NA"]) , entry.start, entry.end,
                  entry.strand, seq))


def read_file(seq_file, gff_file, target_folder, features):
    fastas = []
    cdss_f = []
    cdss_r = []
    genes = []
    with open(seq_file, "r") as seq_f:
        for line in seq_f:
            if line.startswith(">"):
                continue
            else:
                line = line.strip()
                fastas.append(line)
    fasta = "".join(fastas)
    g_h = open(gff_file)
    for entry in Gff3Parser().entries(g_h):
        if os.path.exists(os.path.join(target_folder,
                          "_".join([entry.seq_id, "target.fa"]))):
            os.remove(os.path.join(target_folder,
                      "_".join([entry.seq_id, "target.fa"])))
        for feature in features:
            if (entry.feature == feature) and (entry.strand == "+"):
                cdss_f.append(entry)
            elif (entry.feature == feature) and (entry.strand == "-"):
                cdss_r.append(entry)
        if entry.feature == "gene":
            genes.append(entry)
    g_h.close()
    cdss_f = sorted(cdss_f, key=lambda k: (k.seq_id, k.start,
                                           k.end, k.strand))
    cdss_r = sorted(cdss_r, key=lambda k: (k.seq_id, k.start,
                                           k.end, k.strand))
    genes = sorted(genes, key=lambda k: (k.seq_id, k.start,
                                         k.end, k.strand))
    return fasta, cdss_f, cdss_r, genes


def check_parent_gene(cds, genes):
    target_gene = None
    for gene in genes:
        if (gene.seq_id == cds.seq_id) and (
                gene.strand == cds.strand) and (
                gene.start > cds.end):
            break
        elif "Parent" in cds.attributes.keys():
            if (gene.attributes["ID"] in
                    cds.attributes["Parent"].split(",")):
                target_gene = gene
    if target_gene is None:
        for gene in genes:
            if (gene.seq_id == cds.seq_id) and (
                    gene.strand == cds.strand):
                if ((cds.start <= gene.start) and (
                        cds.end >= gene.end)) or (
                        (cds.start >= gene.start) and (
                        cds.end <= gene.end)) or (
                        (cds.start <= gene.start) and (
                        cds.end <= gene.end) and (
                        cds.end >= gene.start)) or (
                        (cds.start >= gene.start) and (
                        cds.start <= gene.end) and (
                        cds.end >= gene.end)):
                    target_gene = gene
                if (cds.start == gene.start) and (
                        cds.end == gene.end):
                    target_gene = gene
                    break
    return target_gene


def deal_cds_forward(cdss_f, target_folder, fasta, genes, tar_start, tar_end):
    '''for forward strand'''
    pre_id = ""
    out = None
    for cds in cdss_f:
        if cds.seq_id != pre_id:
            out = open(os.path.join(target_folder,
                       "_".join([cds.seq_id, "target.fa"])), "w")
            pre_id = cds.seq_id
        if (cds.start > tar_start):
            start = cds.start - tar_start
        else:
            start = 1
        if ((cds.start + tar_end) < len(fasta)) and (
                (cds.end - cds.start) >= tar_end):
            end = cds.start + tar_end - 1
        elif cds.start + tar_end >= len(fasta):
            end = len(fasta)
        elif (cds.end - cds.start) < tar_end:
            end = cds.end
        seq = Helper().extract_gene(fasta, start, end, cds.strand)
        target = cds
        target_gene = check_parent_gene(cds, genes)
        print_fasta(target, seq, out, target_gene)
    if out is not None:
        out.close()


def deal_cds_reverse(cdss_r, target_folder, fasta, genes, tar_start, tar_end):
    '''for the reverse strand'''
    pre_id = ""
    out = None
    for cds in cdss_r:
        if cds.seq_id != pre_id:
            out = open(os.path.join(target_folder,
                       "_".join([cds.seq_id, "target.fa"])), "a")
            pre_id = cds.seq_id
        if (len(fasta) - cds.end > tar_start):
            end = cds.end + tar_start
        else:
            end = len(fasta)
        if ((cds.end - tar_end) > 1) and ((cds.end - cds.start) >= tar_end):
            start = cds.end - tar_end - 1
        elif cds.end - tar_end < 1:
            start = 1
        elif (cds.end - cds.start) < tar_end:
            start = cds.start
        seq = Helper().extract_gene(fasta, start, end, cds.strand)
        target = cds
        target_gene = check_parent_gene(cds, genes)
        print_fasta(target, seq, out, target_gene)
    if out is not None:
        out.close()


def potential_target(gff_file, seq_file, target_folder, args_tar):
    '''get the sequence of the potential target of sRNA'''
    fasta, cdss_f, cdss_r, genes = read_file(seq_file, gff_file,
                                             target_folder, args_tar.features)
    sort_cdss_f = sorted(cdss_f, key=lambda k: (k.seq_id, k.start,
                                                k.end, k.strand))
    deal_cds_forward(sort_cdss_f, target_folder, fasta, genes,
                     args_tar.tar_start, args_tar.tar_end)
    sort_cdss_r = sorted(cdss_r, reverse=True,
                         key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    deal_cds_reverse(sort_cdss_r, target_folder, fasta, genes,
                     args_tar.tar_start, args_tar.tar_end)
