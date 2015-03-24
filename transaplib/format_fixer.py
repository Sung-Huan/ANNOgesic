#!/usr/bin/python

import os
import sys
import csv
import argparse
from transaplib.gff3 import Gff3Parser

class Format_Fixer(object):

    def _read_gff(self, gff_file, genes, datas, strain):
        gene_num = 0
        gff_parser = Gff3Parser()
        for entry in gff_parser.entries(open(gff_file, "r")):
            entry.seq_id = strain
            entry.info_without_attributes = "\t".join([str(field) for field in [
                        entry.seq_id, entry.source, entry.feature, entry.start,
                        entry.end, entry.score, entry.strand, entry.phase]])
            datas.append(entry)
            if entry.feature == "gene":
                if "locus_tag" in entry.attributes.keys():
                    name = entry.attributes["locus_tag"]
                if "gene" in entry.attributes.keys():
                    name = entry.attributes["gene"]
                entry.attribute_string = ";".join(["ID=gene" + str(gene_num),
                                         "Name=" + name, entry.attribute_string])
                gene_id = "gene" + str(gene_num)
                entry.attributes["ID"] = gene_id
                genes.append(entry)
                gene_num += 1

    def Fix_ratt(self, gff_file, strain, out_file):
        out = open(out_file, "w")
        out.write("##gff-version 3\n")
        cds_num = 0
        rna_num = 0
        gene_num = 0
        genes = []
        datas = []
        check_parent = False
        self._read_gff(gff_file, genes, datas, strain)
        check_parent = False
        for data in datas:
            if data.feature == "gene":
                data = genes[gene_num]
                gene_num += 1
            elif (data.feature == "rRNA") or \
                 (data.feature == "tRNA"):
                name = data.attributes["locus_tag"]
                data.attribute_string = ";".join(["ID=rna" + str(rna_num),
                                        "Name=" + name, data.attribute_string])
                rna_num += 1
            elif data.feature == "CDS":
                if "protein_id" in data.attributes.keys():
                    name = data.attributes["protein_id"]
                for gene in genes:
                    if ((gene.start <= data.start) and \
                       (gene.end >= data.end)) or \
                       (gene.attributes["locus_tag"] == data.attributes["locus_tag"]):
                        data.attribute_string = ";".join(["ID=cds" + str(cds_num), 
                                                "Name=" + name, "Parent=" + gene.attributes["ID"],
                                                data.attribute_string])
                        check_parent = True
                        break
                if check_parent:
                    check_parent = False
                    pass
                else:
                    data.attribute_string = ";".join(["ID=cds" + str(cds_num), "Name=" + name,
                                                      data.attribute_string])
                cds_num += 1
            if "group" in data.attributes.keys():
                ref_f = open(gff_file, "r")
                for ref in gff_parser.entries(ref_f):
                    if "group" in ref.attributes.keys():
                        if (data.attributes["group"] == ref.attributes["group"]):
                            if (data.strand != ref.strand):
                                data.strand = ref.strand
                            break
                ref_f.close()
            out.write("\t".join([data.info_without_attributes, data.attribute_string]) + "\n")

    def Fix_rnaplex(self, rnaplex_file, out_file):
        out = open(out_file, "w")
        with open(rnaplex_file, "r") as f_h:
            for line in f_h:
                line = line.strip()
                if line != "Error during initialization of the duplex in duplexfold_XS":
                    out.write(line + "\n")
