import os
import sys
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.helper import Helper
from annogesiclib.parser_wig import WigParser
from annogesiclib.gen_TSS_type import compare_tss_cds, fix_primary_type
from annogesiclib.lib_reader import read_libs, read_wig


def get_upstream(seq, tss, out, name, nt_before):
    if tss.strand == "+":
        if (tss.start - nt_before + 1) <= 0:
            start = 1
        else:
            start = tss.start - nt_before + 1
        fasta = Helper().extract_gene(seq, start,
                                      tss.start, tss.strand)
    else:
        if (tss.start + nt_before - 1) > len(seq):
            end = len(seq)
        else:
            end = tss.start + nt_before - 1
        fasta = Helper().extract_gene(seq, tss.start,
                                      end, tss.strand)
    if len(fasta) >= nt_before:
        out.write("{0}\n{1}\n".format(name, fasta))


def print_fasta(seq, tss, files, name, nt_before):
    for key in seq.keys():
        if tss.seq_id == key:
            if "Primary" in tss.attributes["type"]:
                get_upstream(seq[key], tss, files["pri"], name, nt_before)
            if "Secondary" in tss.attributes["type"]:
                get_upstream(seq[key], tss, files["sec"], name, nt_before)
            if "Internal" in tss.attributes["type"]:
                get_upstream(seq[key], tss, files["inter"], name, nt_before)
            if "Antisense" in tss.attributes["type"]:
                get_upstream(seq[key], tss, files["anti"], name, nt_before)
            if "Orphan" in tss.attributes["type"]:
                get_upstream(seq[key], tss, files["orph"], name, nt_before)


def read_data(tss_file, fasta_file):
    seq = {}
    tsss = []
    t_f = open(tss_file, "r")
    for entry in Gff3Parser().entries(t_f):
        tsss.append(entry)
    if fasta_file is not None:
        with open(fasta_file, "r") as f_h:
            for line in f_h:
                line = line.strip()
                if line[0] == ">":
                    seq[line[1:]] = ""
                    seq_id = line[1:]
                else:
                    seq[seq_id] = seq[seq_id] + line
    tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    return tsss, seq


def read_gff(gff_file):
    cdss = []
    genes = []
    g_f = open(gff_file, "r")
    for entry in Gff3Parser().entries(g_f):
        if (Helper().feature_without_notgene(entry)):
            cdss.append(entry)
        if entry.feature == "gene":
            genes.append(entry)
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    genes = sorted(genes, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    return cdss, genes


def upstream(tss_file, fasta_file, gff_file, out_class, args_pro, prefix):
    '''get the upstream sequence of TSS'''
    if fasta_file is not None:
        files = {"pri": open("tmp/primary.fa", "w"),
                 "sec": open("tmp/secondary.fa", "w"),
                 "inter": open("tmp/internal.fa", "w"),
                 "anti": open("tmp/antisense.fa", "w"),
                 "orph": open("tmp/orphan.fa", "w")}
    tsss, seq = read_data(tss_file, fasta_file)
    num_tss = 0
    if not args_pro.source:
        out = open(out_class, "w")
        out.write("##gff-version 3\n")
        cdss, genes = read_gff(gff_file)
    for tss in tsss:
        if ("type" not in tss.attributes.keys()) and (args_pro.source):
            print("Error: The TSS gff file may not generated from ANNOgesic."
                  "Please run with --tss_source!")
            sys.exit()
        if args_pro.source:
            name = ">" + "_".join([str(tss.start), tss.strand, tss.seq_id])
            print_fasta(seq, tss, files, name, args_pro.nt_before)
        else:
            tss_type = compare_tss_cds(tss, cdss, genes)
            tss.attributes = tss_type[1]
            tss.attributes["ID"] = tss.seq_id + "_tss" + str(num_tss)
            tss.attribute_string = "".join([
                tss_type[0], ";ID=", tss.seq_id, "_tss", str(num_tss)])
            num_tss += 1
    if not args_pro.source:
        if args_pro.tex_wigs is not None:
            libs, texs = read_libs(args_pro.input_libs, args_pro.tex_wigs)
            wigs_f = read_wig(os.path.join(
                args_pro.wig_path, prefix + "_forward.wig"), "+", libs)
            wigs_r = read_wig(os.path.join(
                args_pro.wig_path, prefix + "_reverse.wig"), "+", libs)
        else:
            wigs_f = None
            wigs_r = None
        sort_tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start,
                                                k.end, k.strand))
        final_tsss = fix_primary_type(sort_tsss, wigs_f, wigs_r)
        for tss in final_tsss:
            name = ">" + "_".join([str(tss.start), tss.strand, tss.seq_id])
            tss.attribute_string = ";".join(
                ["=".join(items) for items in tss.attributes.items()])
            out.write("\t".join([str(field) for field in [
                            tss.seq_id, tss.source, tss.feature, tss.start,
                            tss.end, tss.score, tss.strand, tss.phase,
                            tss.attribute_string]]) + "\n")
            if fasta_file is not None:
                print_fasta(seq, tss, files, name, args_pro.nt_before)


def del_repeat_fasta(input_file, out_file):
    data = {}
    check_same = False
    first_file = True
    out = open(out_file, "w")
    with open(input_file, "r") as f_h:
        for line in f_h:
            line = line.strip()
            if line[0] == ">":
                if check_same:
                    check_same = False
                if first_file:
                    seq_id = line[1:]
                    first_file = False
                    data[seq_id] = ""
                else:
                    if line[1:] in data.keys():
                        check_same = True
                    else:
                        seq_id = line[1:]
                        data[seq_id] = ""
            else:
                if check_same:
                    pass
                else:
                    data[seq_id] = data[seq_id] + line
    for strain, fasta in data.items():
        out.write(">" + strain + "\n")
        out.write(fasta + "\n")
