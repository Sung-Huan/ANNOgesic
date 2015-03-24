#!/usr/bin/python

import os	
import sys
import csv
from transaplib.gff3 import Gff3Parser


def _import_uniprot_data(entry, name_list, feature):
    ref_name = entry.attributes[feature]
    if ref_name not in name_list:
        name_list.add(ref_name)

def Retrieve_Uniprot(database_file, gff_file, out_file):
    name_list = set()
    gffs = []
    out = open(out_file, "w")
    out.write("\t".join(["strain", "strand", "start", "end", 
                         "protein_id", "Go_term"]) + "\n")
    for entry in Gff3Parser.entries(open(gff_file)):
	if entry.feature == "CDS":
	    if ("Name" in entry.attributes.keys()) and \
               ("protein_id" in entry.attributes.keys()):
                if entry.attributes["Name"] == entry.attributes["protein_id"]:
                    _import_uniprot_data(entry, name_list, "Name")
                else:
                    _import_uniprot_data(entry, name_list, "Name")
                    _import_uniprot_data(entry, name_list, "protein_id")
            elif ("Name" in entry.attributes.keys()):
                _import_uniprot_data(entry, name_list, "Name")
            elif ("protein_id" in entry.attributes.keys()):
                _import_uniprot_data(entry, name_list, "protein_id")
            else:
                ref_name = ''
            gffs.append(entry)
    idmapping = open(database_file,"r")
    for uni_id in idmapping:
        uni_line = uni_id.rstrip("\n")
        uni_lines = uni_line.split("\t")
        if uni_lines[3] in name_list:
            for gff in gffs:
                if (uni_lines[3] == gff.attributes["Name"]) or \
                   (uni_lines[3] == gff.attributes["protein_id"]):
                    out.write("\t".join([gff.seq_id, gff.strand, str(gff.start), 
                                   str(gff.end), uni_lines[3], uni_lines[6]]) + "\n")
