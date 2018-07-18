import os
import csv
from Bio import SeqIO
from collections import defaultdict
from annogesiclib.gff3 import Gff3Parser, Gff3Entry
from annogesiclib.TSSpredator_parser import TSSPredatorReader
from annogesiclib.helper import Helper


class Converter(object):
    '''Converting from one format to another format'''

    def __init__(self):
        self.gff3parser = Gff3Parser()
        self.tssparser = TSSPredatorReader()

    def _check_locus_tag(self, entry, genes):
        gene_tag = "-"
        locus_tag = "-"
        if "locus_tag" in entry.attributes.keys():
            locus_tag = entry.attributes["locus_tag"]
        elif "Parent" in entry.attributes.keys():
            for gene in genes:
                if (gene.attributes["ID"] in
                        entry.attributes["Parent"].split(",")):
                    if "gene" in gene.attributes.keys():
                        gene_tag = gene.attributes["gene"]
                    if "locus_tag" in gene.attributes.keys():
                        locus_tag = gene.attributes["locus_tag"]
        if locus_tag == "-":
            locus_tag = "".join([
                entry.feature, ":", str(entry.start), "-",
                str(entry.end), "_", entry.strand])
        return locus_tag, gene_tag


    def _print_rntptt_file(self, out, entrys, genes):
        '''output to rnt and ptt file'''
        for entry in entrys:
            location = "..".join([str(entry.start), str(entry.end)])
            length = str(entry.end - entry.start + 1)
            if entry.feature == "CDS":
                if "protein_id" in entry.attributes.keys():
                    pid = entry.attributes["protein_id"]
                else:
                    pid = "-"
                gene_tag, locus_tag = self._check_locus_tag(entry, genes)
            else:
                pid = "-"
                gene_tag = "-"
                gene_tag, locus_tag = self._check_locus_tag(entry, genes)
            if "product" in entry.attributes.keys():
                product = entry.attributes["product"]
            else:
                product = "-"
            out.write("\t".join([location, entry.strand, length,
                                 pid, gene_tag, locus_tag, "-", "-",
                                 product]) + "\n")

    def _print_rntptt_title(self, out, num, seq_id, length):
        '''print the title of rnt and ptt file'''
        out.write(seq_id + " - 1.." + length + "\n")
        out.write(num + " proteins\n")
        out.write("\t".join(["Location", "Strand", "Length", "PID",
                  "Gene", "Synonym", "Code", "COG", "Product"]) + "\n")

    def _read_file(self, gff_file, fasta_file, rnas, cdss, genes):
        num_cds = 0
        num_rna = 0
        seq = ""
        g_f = open(gff_file, "r")
        for entry in self.gff3parser.entries(g_f):
            if (entry.feature == "rRNA") or (entry.feature == "tRNA"):
                num_rna += 1
                rnas.append(entry)
            elif entry.feature == "CDS":
                num_cds += 1
                cdss.append(entry)
            elif entry.feature == "gene":
                genes.append(entry)
        g_f.close()
        if fasta_file == "0":
            seq = "-1"
        else:
            with open(fasta_file, "r") as f_f:
                for line in f_f:
                    line = line.strip()
                    if len(line) != 0:
                        if line[0] != ">":
                            seq = seq + line
        return (num_cds, num_rna, seq)

    def _srna2rntptt(self, srna_input_file, srna_output_file, srnas, length):
        '''convert the sRNA gff file to rnt file'''
        num_srna = 0
        r_s = open(srna_input_file, "r")
        for entry in Gff3Parser().entries(r_s):
            num_srna += 1
            srnas.append(entry)
        srnas = sorted(srnas, key=lambda k: (k.seq_id, k.start,
                                             k.end, k.strand))
        r_s.close()
        out_s = open(srna_output_file, "w")
        self._print_rntptt_title(out_s, str(num_srna),
                                 srnas[0].seq_id, str(length))
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
        '''deal with the embl file which contain join'''
        info = info.replace("(", "")
        info = info.replace(")", "")
        info = info.replace("join", "")
        joins = info.split(",")
        return joins

    def _multi_embl_pos(self, row):
        '''deal with the feature which has multiple positions'''
        poss = []
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
            strand = "+"
        source = row[5:21].rstrip()
        return {"pos": poss, "strand": strand, "source": source}

    def _parser_embl_data(self, embl_file, out):
        '''Parser of embl file for converting to other format'''
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
                    name[0] = name[0].replace("ID", "")
                    if "SV" in name[1]:
                        version = name[1].split(" ")[-1]
                        id_name = ".".join([name[0].strip(), version.strip()])
                    else:
                        id_name = name[0].strip()
                if (row[0:2] == "FT"):
                    if row[5] != " ":
                        note_name = row[5:9]
                    if row[5:11] == "source":
                        info = self._multi_embl_pos(row)
                    if (note_name != "misc") and (row[5] == " ") and (
                            row[21] == "/"):
                        if first:
                            first = False
                        else:
                            line = line + ";"
                        data = row[22:].replace(";", ",")
                        data = data.split("=")
                        try:
                            note = data[1].replace("\"", "")
                            line = line + data[0] + "=" + note
                        except:
                            line = line + data[0] + "=" + "True"
                    if (note_name != "misc") and (row[5] == " ") and (
                            row[21] != "/"):
                        note = row[21:].replace("\"", "")
                        note = note.replace(";", ",")
                        line = line + " " + note
                    if (note_name != "misc") and (row[5] != " ") and (
                            row[5:11] != "source"):
                        first = True
                        if info != "Wrong":
                            for pos in info["pos"]:
                                out.write(("{0}\tRefseq\t{1}\t{2}\t{3}"
                                           "\t.\t{4}\t.\t{5}\n").format(
                                          id_name, info["source"],
                                          pos["start"], pos["end"],
                                          info["strand"], line))
                        if (row[5:8] != "CDS") and (row[5:9] != "misc"):
                            info = self._multi_embl_pos(row)
                        elif (row[5:8] == "CDS"):
                            info = self._multi_embl_pos(row)
                        line = ""
        return (id_name, info, line)

    def _assign_tss_type(self, tss, utr_pri, utr_sec):
        '''Assigning the TSS types'''
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
        '''deal with the TSS which has multiple TSS types'''
        tss_type = self._assign_tss_type(tss, utrs["pri"], utrs["sec"])
        if (tss_type not in tss_features["tss_types"]) or (
                tss.locus_tag not in tss_features["locus_tags"]):
            if (tss_type not in tss_features["tss_types"]):
                tss_index[tss_type] += 1
                nums["tss"] += 1
            if (nums["class"] == 1):
                tss_features["tss_types"].append(tss_type)
                tss_features["utr_lengths"].append(
                        tss_type + "_" + tss.utr_length)
                tss_features["locus_tags"].append(tss.locus_tag)
            else:
                if tss_type not in tss_features["tss_types"]:
                    tss_features["tss_types"].append(tss_type)
                    tss_features["utr_lengths"].append(
                            tss_type + "_" + tss.utr_length)
                    tss_features["locus_tags"].append(tss.locus_tag)
            nums["class"] += 1

    def _uni_tss_class(self, tss, utrs, tss_index, tss_features, nums):
        '''It is for TSS which has only one type'''
        tss_type = self._assign_tss_type(tss, utrs["pri"], utrs["sec"])
        tss_index[tss_type] += 1
        tss_features["tss_types"].append(tss_type)
        tss_features["utr_lengths"].append(tss_type+"_"+tss.utr_length)
        tss_features["locus_tags"].append(tss.locus_tag)
        nums["tss"] += 1

    def _print_tssfile(self, nums, tss_features, tss, tss_pro,
                       strain, method, out, tss_libs):
        '''print gff file of TSS'''
        tss_merge_type = ",".join(tss_features["tss_types"])
        utr_length = ",".join(tss_features["utr_lengths"])
        merge_locus_tag = ",".join(tss_features["locus_tags"])
        libs = ",".join(tss_libs)
        strand = Helper().get_strand_name(tss.super_strand)
        attribute_string = ";".join(
                          ["=".join(items) for items in (
                              ["Name", "".join([tss_pro, ":",
                               str(tss.super_pos), "_", strand])],
                              ["ID", strain + "_" + tss_pro.lower() 
                               + str(nums["tss_uni"])],
                              ["type", tss_merge_type],
                              ["utr_length", str(utr_length)],
                              ["associated_gene", merge_locus_tag],
                              ["libs", libs], ["method", "TSSpredator"])])
        out.write("\t".join([strain, method, tss_pro, str(tss.super_pos),
                             str(tss.super_pos), ".", tss.super_strand, ".",
                             attribute_string]) + "\n")

    def convert_gff2rntptt(self, gff_file, fasta_file, ptt_file, rnt_file,
                           srna_input_file, srna_output_file):
        '''Convert gff format to rnt and ptt format'''
        genes = []
        rnas = []
        cdss = []
        srnas = []
        datas = self._read_file(gff_file, fasta_file, rnas, cdss, genes)
        rnas = sorted(rnas, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
        cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
        genes = sorted(genes, key=lambda k: (k.seq_id, k.start,
                                             k.end, k.strand))
        num_cds = datas[0]
        num_rna = datas[1]
        seq = datas[2]
        out_p = open(ptt_file, "w")
        out_r = open(rnt_file, "w")
        if len(cdss) != 0:
            self._print_rntptt_title(out_p, str(num_cds),
                                     cdss[0].seq_id, str(len(seq)))
            self._print_rntptt_file(out_p, cdss, genes)
        if len(rnas) != 0:
            self._print_rntptt_title(out_r, str(num_rna),
                                     rnas[0].seq_id, str(len(seq)))
            self._print_rntptt_file(out_r, rnas, genes)
        out_p.close()
        out_r.close()
        if (srna_input_file is not None) and \
           (srna_output_file is not None):
            self._srna2rntptt(srna_input_file, srna_output_file,
                              srnas, str(len(seq)))
        elif (srna_input_file is None) and \
             (srna_output_file is None):
            pass
        else:
            print("Error: Lack sRNA input gff files or "
                  "the name sRNA output rnt files\n")

    def convert_embl2gff(self, embl_file, gff_file):
        '''Convert embl format to gff format'''
        info = "Wrong"
        out = open(gff_file, "w")
        out.write("##gff-version 3\n")
        datas = self._parser_embl_data(embl_file, out)
        id_name = datas[0]
        info = datas[1]
        line = datas[2]
        if info != "Wrong":
            for pos in info["pos"]:
                out.write(("{0}\tRefseq\t{1}\t{2}\t{3}"
                           "\t.\t{4}\t.\t{5}\n").format(
                          id_name, info["source"], pos["start"],
                          pos["end"], info["strand"], line))
        out.close()

    def _get_libs(self, tss_file):
        '''Get the library which can detect this specific TSS'''
        tss_libs = {}
        tss_fh = open(tss_file, "r")
        for tss in self.tssparser.entries(tss_fh):
            key = "_".join([str(tss.super_pos), tss.super_strand])
            if key not in tss_libs.keys():
                tss_libs[key] = []
            if (tss.is_detected) and (tss.genome not in tss_libs[key]):
                tss_libs[key].append(tss.genome)
        tss_fh.close()
        return tss_libs

    def convert_mastertable2gff(self, tss_file, method, tss_pro,
                                strain, out_gff):
        '''Convert MasterTable to gff format'''
        temps = {"tss": 0, "strand": "#"}
        nums = {"tss": 0, "tss_uni": 0, "class": 1}
        check_print = False
        utrs = {"total": [], "pri": [], "sec": []}
        tss_features = {"tss_types": [], "locus_tags": [], "utr_lengths": []}
        tss_index = defaultdict(lambda: 0)
        tss_fh = open(tss_file, "r")
        out = open(out_gff, "w")
        out.write("##gff-version 3\n")
        tss_libs = self._get_libs(tss_file)
        detect_run = False
        for tss in self.tssparser.entries(tss_fh):
            detect_run = True
            key = "_".join([str(tss.super_pos), tss.super_strand])
            if ((tss.super_pos == temps["tss"])) and (
                    temps["strand"] == tss.super_strand) and (
                        tss.class_count == 1):
                pass
            else:
                if ((tss.super_pos != temps["tss"])) or (
                        temps["strand"] != tss.super_strand):
                    check_print = False
                    nums["class"] = 1
                    if tss.utr_length != "NA":
                        utrs["total"].append(int(tss.utr_length))
                temps["tss"] = tss.super_pos
                temps["strand"] = tss.super_strand
                if (tss.class_count != 1) and (
                        nums["class"] <= tss.class_count):
                    self._multi_tss_class(tss, tss_index, tss_features,
                                          nums, utrs)
                if (tss.class_count == 1) or (
                        nums["class"] > tss.class_count):
                    if (tss.class_count == 1):
                        self._uni_tss_class(tss, utrs, tss_index,
                                            tss_features, nums)
                    if (check_print is False):
                        self._print_tssfile(nums, tss_features, tss, tss_pro,
                                            strain, method, out, tss_libs[key])
                        check_print = True
                        nums["tss_uni"] += 1
                    tss_features = {"tss_types": [], "locus_tags": [],
                                    "utr_lengths": []}
        if (check_print is False) and detect_run:
            self._print_tssfile(nums, tss_features, tss, tss_pro,
                                strain, method, out, tss_libs[key])
        tss_fh.close()
        out.close()

    def convert_transtermhp2gff(self, transterm_file, gff_file):
        '''Convert the output of TransTermHP to gff format'''
        out = open(gff_file, "w")
        out.write("##gff-version 3\n")
        terms = []
        with open(transterm_file) as t_h:
            for line in t_h:
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
                    "seq_id": transterm_file.split("/")[-1][:-1 * len(
                               "_best_terminator_after_gene.bag")],
                    "source": "TransTermHP",
                    "feature": "terminator",
                    "start": start,
                    "end": end,
                    "score": ".",
                    "strand": strand,
                    "phase": ".",
                    "attributes": "associated_gene=%s" % (gene)
                    })
                terms.append(entry)
        sort_terms = sorted(terms, key=lambda k: (k.seq_id, k.start,
                                                  k.end, k.strand))
        num = 0
        for term in sort_terms:
            out.write("\t".join([str(field) for field in [
                                term.seq_id, term.source, term.feature,
                                term.start, term.end, term.score,
                                term.strand, term.phase,
                                term.attribute_string]]))
            name = '%0*d' % (5, num)
            out.write(";ID={0}_terminator{1};Name=terminator_{2}\n".format(
                      term.seq_id, num, name))
            num += 1
        out.close()

    def convert_circ2gff(self, circ_file, args_circ, out_all, out_filter):
        '''Convert the circRNA output of segemehl to gff format'''
        circs = []
        out_a = open(out_all, "w")
        out_f = open(out_filter, "w")
        out_a.write("##gff-version 3\n")
        out_f.write("##gff-version 3\n")
        f_h = open(circ_file, "r")
        for row in csv.reader(f_h, delimiter='\t'):
            if row[0] != "Genome":
                circs.append({"strain": row[0],
                              "strand": row[1], "start": int(row[2]),
                              "end": int(row[3]), "conflict": row[4],
                              "depth": int(row[5]), "per_start": float(row[6]),
                              "per_end": float(row[7])})
        circs = sorted(circs, key=lambda k: (k["strain"], k["start"],
                                             k["end"], k["strand"]))
        id_ = 0
        for circ in circs:
            attribute_string = ";".join(["=".join(items) for items in [
                                           ("ID", circ["strain"] +
                                            "_circrna" + str(id_)),
                                           ("name", "circRNA_" + str(id_)),
                                           ("support_reads",
                                            str(circ["depth"])),
                                           ("read_at_start",
                                            str(circ["per_start"])),
                                           ("read_at_end",
                                            str(circ["per_end"])),
                                           ("conflict", circ["conflict"]),
                                           ("method", "segemehl")]])
            out_a.write("\t".join([str(field) for field in [
                            circ["strain"], "ANNOgesic", "circRNA",
                            str(circ["start"]), str(circ["end"]),
                            ".", circ["strand"], ".",
                            attribute_string]]) + "\n")
            if (circ["depth"] >= args_circ.support) and (
                    circ["conflict"] == "NA") and (
                    circ["per_start"] >= args_circ.start_ratio) and (
                    circ["per_end"] >= args_circ.end_ratio):
                out_f.write("\t".join([str(field) for field in [
                            circ["strain"], "ANNOgesic", "circRNA",
                            str(circ["start"]), str(circ["end"]),
                            ".", circ["strand"], ".",
                            attribute_string]]) + "\n")
            id_ += 1
        f_h.close()
        out_a.close()
        out_f.close()

    def convert_gbk2embl(self, input_folder):
        """Convert gbk to embl."""
        print("Converting gbk files to embl files")
        for annotation_file in os.listdir(input_folder):
            if annotation_file[-3:] == "gbk":
                gbk_file = annotation_file
                embl_file = gbk_file[0:-3] + "embl"
                gbk_entry = SeqIO.parse(os.path.join(
                                        input_folder, gbk_file), "genbank")
                count = SeqIO.write(gbk_entry, os.path.join(
                                    input_folder, embl_file), "embl")
                print("Converted %i records" % count)
