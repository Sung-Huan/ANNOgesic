import sys
import os
import csv
import re
import itertools
import math
import shutil
from Bio import SeqIO
from subprocess import call
from collections import defaultdict
from annogesiclib.gff3 import Gff3Parser, Gff3Entry
from annogesiclib.TSSpredator import TSSPredatorReader

class Converter(object):

    def __init__(self):
        self.gff3parser = Gff3Parser()
        self.tssparser = TSSPredatorReader()

    def _print_rntptt_file(self, out, entrys, genes):
        for entry in entrys:
            gene_tag = "-"
            locus_tag = "-"
            location = "..".join([str(entry.start), str(entry.end)])
            length = str(entry.end - entry.start + 1)
            if entry.feature == "CDS":
                if "protein_id" in entry.attributes.keys():
                    pid = entry.attributes["protein_id"]
                else:
                    pid = "-"
                if "Parent" in entry.attributes.keys():
                    for gene in genes:
                        if gene.attributes["ID"] == entry.attributes["Parent"]:
                            if "gene" in gene.attributes.keys():
                                gene_tag = gene.attributes["gene"]
                            locus_tag = gene.attributes["locus_tag"]
            else:
                pid = "-"
                gene_tag = "-"
                if "locus_tag" in entry.attributes.keys():
                    locus_tag = entry.attributes["locus_tag"]
                elif "Parent" in entry.attributes.keys():
                    for gene in genes:
                        if gene.attributes["ID"] == entry.attributes["Parent"]:
                            if "gene" in gene.attributes.keys():
                                gene_tag = gene.attributes["gene"]
                            locus_tag = gene.attributes["locus_tag"]
            product = entry.attributes["product"]
            out.write("\t".join([location, entry.strand, length,
                                 pid, gene_tag, locus_tag, "-", "-",
                                 product]) + "\n")
    
    def _print_rntptt_title(self, out, num, seq_id, length):
        out.write(seq_id + " - 1.." + length + "\n")
        out.write(num + " proteins\n")
        out.write("\t".join(["Location", "Strand", "Length", "PID", 
                  "Gene", "Synonym Code", "COG", "Product"]) + "\n")    
    
    def _read_file(self, gff_file, fasta_file, rnas, cdss, genes):
        num_cds = 0
        num_rna = 0
        seq = ""
        g_f = open(gff_file, "r")
        for entry in self.gff3parser.entries(g_f):
            if (entry.feature == "rRNA") or \
               (entry.feature == "tRNA"):
                num_rna += 1
                rnas.append(entry)
            elif entry.feature == "CDS":
                num_cds += 1
                cdss.append(entry)
            elif entry.feature == "gene":
                genes.append(entry)
        g_f.close()
        with open(fasta_file, "r") as f_f:
            for line in f_f:
                line = line.strip()
                if line[0] != ">":
                    seq = seq + line    
        return (num_cds, num_rna, seq)

    def _srna2rntptt(self, sRNA_input_file, sRNA_output_file, srnas, length):
        num_srna = 0
        r_s = open(sRNA_input_file, "r")
        for entry in Gff3Parser().entries(r_s):
            num_srna += 1
            srnas.append(entry)
        srnas = sorted(srnas, key=lambda k: (k.seq_id, k.start))
        r_s.close()
        out_s = open(sRNA_output_file, "w")
        self._print_rntptt_title(out_s, str(num_srna), srnas[0].seq_id, str(length))
        num_srna = 0
        for srna in srnas:
            num_srna += 1
            name = '%0*d' % (5, num_srna)
            gene_tag = "-"
            locus_tag = "ncRNA_" + name
            pid = "ncRNA_" + name
            product = "sRNA"
            location = "..".join([str(srna.start), str(srna.end)])
            length = str(srna.end - srna.start + 1)
            out_s.write("\t".join([location, srna.strand, length,
                                 pid, gene_tag, locus_tag, "-", "-",
                                 product]) + "\n")
        out_s.close()

    def _deal_embl_join(self, info):
        info = info.replace("(", "")
        info = info.replace(")", "")
        info = info.replace("join", "")
        joins = info.split(",")
        return joins
    
    def _multi_embl_pos(self, row, mark, id_name):
        poss = []
        num_join = 1
        if row[21:31] == "complement":
            comple = row[32:-1]
            if comple.find("join") != -1:
                joins = self._deal_embl_join(comple)
                for join in joins:
                    pos = join.split("..")
                    if len(pos) < 2:
                        return "Wrong"
                    poss.append({"start": pos[0], "end": pos[1]})
            else:
                pos = comple.split("..")
                if len(pos) < 2:
                    return "Wrong"
                poss.append({"start": pos[0], "end": pos[1]})
            strand = "-"
        else:
            if row[21:].find("join") != -1:
                joins = self._deal_embl_join(row[21:])
                for join in joins:
                    pos = join.split("..")
                    if len(pos) < 2:
                        return "Wrong"
                    poss.append({"start": pos[0], "end": pos[1]})
            else:
                pos = row[21:].split("..")
                if len(pos) < 2:
                    return "Wrong"
                poss.append({"start": pos[0], "end": pos[1]})
            strand ="+"
        source = row[5:21].rstrip()
        return {"pos": poss, "strand": strand, "source": source}

    def _parser_embl_data(self, embl_file, out):
        first = True
        line = ""
        note_name = ""
        info = "Wrong"
        with open(embl_file, "r") as f_h:
            for row in f_h:
                row = row.strip()
                if row[0:2] == "SQ":
                    break
                if row[0:2] == "ID":
                    name = row.split(";")
                    id_name = name[0][21:]
                if (row[0:2] == "FT"):
                    if row[5] != " ":
                        note_name = row[5:9]
                    if row[5:11] == "source":
                        info = self._multi_embl_pos(row, ".", id_name)
                    if (note_name != "misc")&(row[5] == " ")&(row[21] == "/"):
                        if first:
                            first = False
                        else:
                            line = line + ";"
                        data = row[22:].replace(";", ",")
                        data = data.split("=")
                        try:
                            note = data[1].replace("\"", "")
                            line = line + data[0]+"="+note
                        except:
                            line = line + data[0]+"="+"True"
                    if (note_name != "misc")&(row[5] == " ")&(row[21] != "/"):
                        note = row[21:].replace("\"", "")
                        note = note.replace(";", ",")
                        line = line + " "+note
                    if (note_name != "misc")&(row[5] != " ")&(row[5:11] != "source"):
                        first = True
                        if info != "Wrong":
                            for pos in info["pos"]:
                                out.write(("{0}\tRefseq\t{1}\t{2}\t{3}\t.\t{4}\t.\t{5}\n").format(
                                          id_name, info["source"], pos["start"], 
                                          pos["end"], info["strand"], line))
                        if (row[5:8] != "CDS") and (row[5:9] != "misc"):
                            info = self._multi_embl_pos(row, ".", id_name)
                        elif (row[5:8] == "CDS"):
                            info = self._multi_embl_pos(row, "0", id_name)
                        line = ""

    def _assign_tss_type(self, tss, utr_pri, utr_sec):
        if tss.is_primary:
            tss_type = "Primary"
            utr_pri.append(int(tss.utr_length))
        elif tss.is_secondary:
            tss_type = "Secondary"
            utr_sec.append(int(tss.utr_length))
        elif tss.is_internal:
            tss_type = "Internal"
        elif tss.is_antisense:
            tss_type = "Antisense"
        else:
            tss_type = "Orphan"
        return tss_type
    
    def _multi_tss_class(self, tss, tss_index, tss_features, nums, utrs):
        tss_type = self._assign_tss_type(tss, utrs["pri"], utrs["sec"])
        if (tss_type not in tss_features["tss_types"]) or \
           (tss.locus_tag not in tss_features["locus_tags"]):
            if (tss_type not in tss_features["tss_types"]):
                tss_index[tss_type] += 1
                nums["tss"] += 1
            if (nums["class"] == 1):
                tss_features["tss_types"].append(tss_type)
                tss_features["utr_lengths"].append(tss_type+"_"+tss.utr_length)
                tss_features["locus_tags"].append(tss.locus_tag)
            else:
                if tss_type not in tss_features["tss_types"]:
                    tss_features["tss_types"].append(tss_type)
                    tss_features["utr_lengths"].append(tss_type+"_"+tss.utr_length)
                    tss_features["locus_tags"].append(tss.locus_tag)
            nums["class"] += 1

    def _uni_tss_class(self, tss, utrs, tss_index, tss_features, nums):
        tss_type = self._assign_tss_type(tss, utrs["pri"], utrs["sec"])
        tss_index[tss_type] += 1
        tss_features["tss_types"].append(tss_type)
        tss_features["utr_lengths"].append(tss_type+"_"+tss.utr_length)
        tss_features["locus_tags"].append(tss.locus_tag)
        nums["tss"] += 1

    def _print_tssfile(self, nums, tss_features, tss, tss_pro, 
                       strain, method, out):
        tss_pro = tss_pro[0].upper() + tss_pro[1:]
        tss_merge_type = " ".join(tss_features["tss_types"])
        utr_length = " ".join(tss_features["utr_lengths"])
        merge_locus_tag = " ".join(tss_features["locus_tags"])
        attribute_string = ";".join(
                          ["=".join(items) for items in (["Name", "".join([tss_pro, ":",
                                                          str(tss.super_pos), "_", tss.super_strand])],
                                                         ["ID", tss_pro.lower() + str(nums["tss_uni"])],
                                                         ["type", tss_merge_type],
                                                         ["UTR_length", str(utr_length)],
                                                         ["associated_gene", merge_locus_tag])])
        out.write("\t".join([strain, method, tss_pro, str(tss.super_pos),
                             str(tss.super_pos), ".", tss.super_strand, ".",
                             attribute_string]) + "\n")
    
    def convert_gff2rntptt(self, gff_file, fasta_file, ptt_file, rnt_file, 
                           sRNA_input_file, sRNA_output_file):
        genes = []
        rnas = []
        cdss = []
        srnas = []
        datas = self._read_file(gff_file, fasta_file, rnas, cdss, genes)
        rnas = sorted(rnas, key=lambda k: (k.seq_id, k.start))
        cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start))
        genes = sorted(genes, key=lambda k: (k.seq_id, k.start))
        num_cds = datas[0]
        num_rna = datas[1]
        seq = datas[2]
        out_p = open(ptt_file, "w")
        out_r = open(rnt_file, "w")
        self._print_rntptt_title(out_p, str(num_cds), cdss[0].seq_id, str(len(seq)))
        self._print_rntptt_file(out_p, cdss, genes)
        if len(rnas) != 0:
            self._print_rntptt_title(out_r, str(num_rna), rnas[0].seq_id, str(len(seq)))
            self._print_rntptt_file(out_r, rnas, genes)
        out_p.close()
        out_r.close()
        if (sRNA_input_file is not None) and \
           (sRNA_output_file is not None):
            self._srna2rntptt(sRNA_input_file, sRNA_output_file, srnas, str(len(seq)))
        elif (sRNA_input_file is None) and \
             (sRNA_output_file is None):
            pass
        else:
            print("Error: lack sRNA input gff file or the name sRNA output rnt file\n")

    def convert_embl2gff(self, embl_file, gff_file):
        info = "Wrong"
        out = open(gff_file, "w")
        out.write("##gff-version 3\n")
        self._parser_embl_data(embl_file, out)
        if info != "Wrong":
            for pos in info["pos"]:
                out.write(("{0}\tRefseq\t{1}\t{2}\t{3}\t.\t{4}\t.\t{5}\n").format(
                          id_name, info["source"], pos["start"], pos["end"], info["strand"], line))
    def convert_mastertable2gff(self, tss_file, method, tss_pro, strain, gff_file):
        temps = {"tss": 0, "strand": "#"}
        nums = {"tss": 0, "tss_uni": 0, "class": 1}
        check_print = False
        utrs = {"total": [], "pri": [], "sec": []}
        tss_features = {"tss_types": [], "locus_tags": [], "utr_lengths": []}
        tss_index = defaultdict(lambda: 0)
        tss_fh = open(tss_file, "r");
        out = open(gff_file, "w")
        out.write("##gff-version 3\n")
        for tss in self.tssparser.entries(tss_fh):
            if ((tss.super_pos == temps["tss"])) and \
                (temps["strand"] == tss.super_strand) and \
                (tss.class_count == 1):
                pass
            else:
                if ((tss.super_pos != temps["tss"])) or \
                    (temps["strand"] != tss.super_strand):
                    check_print = False
                    nums["class"] = 1
                    if tss.utr_length != "NA":
                        utrs["total"].append(int(tss.utr_length))
                temps["tss"] = tss.super_pos
                temps["strand"] = tss.super_strand
                if (tss.class_count != 1) and (nums["class"] <= tss.class_count):
                    self._multi_tss_class(tss, tss_index, tss_features, nums, utrs)
                if (tss.class_count == 1) or \
                   (nums["class"] > tss.class_count):
                    if (tss.class_count == 1):
                        self._uni_tss_class(tss, utrs, tss_index, tss_features, nums)
                    if (check_print is False):
                        self._print_tssfile(nums, tss_features, tss, tss_pro, strain, method, out)
                        check_print = True
                        nums["tss_uni"] += 1
                    tss_features = {"tss_types": [], "locus_tags": [], "utr_lengths": []}
        tss_fh.close()
    def convert_transtermhp2gff(self, transterm_file, gff_file):
        out = open(gff_file, "w")
        out.write("##gff-version 3\n")
        num = 0
        terms = []
        for line in open(transterm_file):
            row = line[:-1].split()
            if len(row) < 10:
                continue
            if len(row) == 14:
                start = row[1]
                end = row[3]
                strand = row[4]
                gene = row[0]
            else:
                start = row[0]
                end = row[2]
                strand = row[3]
                gene = "missing"
            entry = Gff3Entry({
                "seq_id" : transterm_file.split("/")[-1][:-1 * len("_best_terminator_after_gene.bag")],
                "source" : "TransTermHP",
                "feature" : "terminator",
                "start" : start,
                "end" : end,
                "score" : ".",
                "strand" : strand,
                "phase" : ".",
                "attributes" : "associated_gene=%s" % (gene)
                })
            terms.append(entry)
        sort_terms = sorted(terms, key=lambda k: (k.seq_id, k.start))
        for term in sort_terms:
            out.write("\t".join([str(field) for field in [
                                term.seq_id, term.source, term.feature, term.start,
                                term.end, term.score, term.strand, term.phase,
                                term.attribute_string]]))
            name = '%0*d' % (5, num)
            out.write(";ID=term{0};Name=terminator_{1}\n".format(num, name))
            num += 1
    
    def convert_circ2gff(self, circ_file, depth, start_ratio, end_ratio, out_all, out_filter):
        circs = []
        out_a = open(out_all, "w")
        out_f = open(out_filter, "w")
        out_a.write("##gff-version 3\n")
        out_f.write("##gff-version 3\n")
        fh = open(circ_file, "r")
        for row in csv.reader(fh, delimiter='\t'):
            if row[0] != "ID":
                circs.append({"name": row[0], "strain": row[1],
                              "strand": row[2], "start": int(row[3]),
                              "end": int(row[4]), "conflict": row[5],
                              "depth": int(row[6]), "per_start":float(row[7]),
                              "per_end":float(row[8])})
        circs = sorted(circs, key=lambda k: (k["strain"], k["start"]))
        for circ in circs:
            ID = circ["name"].split("_")[1]
            attribute_string = ";".join(["=".join(items) for items in [
                                           ("ID", "circrna" + ID),
                                           ("name", circ["name"]),
                                           ("support_reads", str(circ["depth"])),
                                           ("read_at_start", str(circ["per_start"])),
                                           ("read_at_end", str(circ["per_end"])),
                                           ("confliction", circ["conflict"])]])
            out_a.write("\t".join([str(field) for field in [
                            circ["strain"], "segemehl", "circRNA", str(circ["start"]),
                            str(circ["end"]), ".",circ["strand"], ".", attribute_string]]) + "\n")
            if (circ["depth"] >= depth) and \
               (circ["conflict"] == "NA") and \
               (circ["per_start"] >= start_ratio) and \
               (circ["per_end"] >= end_ratio):
                out_f.write("\t".join([str(field) for field in [
                            circ["strain"], "segemehl", "circRNA", str(circ["start"]),
                            str(circ["end"]), ".",circ["strand"], ".", attribute_string]]) + "\n")

    def convert_gbk2embl(self, input_folder):
        """Convert gbk to embl."""
        print("Convert gbk to embl...")
        for annotation_file in os.listdir(input_folder):
            if annotation_file[-3:] == "gbk":
                gbk_file = annotation_file
                embl_file = gbk_file[0:-3] + "embl"
                gbk_entry = SeqIO.parse(os.path.join(input_folder, gbk_file), "genbank")
                count = SeqIO.write(gbk_entry, os.path.join(input_folder, embl_file), "embl")
                print("Converted %i records" % count)
