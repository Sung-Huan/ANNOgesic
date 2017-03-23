import os
import csv
import shutil
from annogesiclib.gff3 import Gff3Parser


def del_attributes(entry, features):
    attributes = {}
    for key, value in entry.attributes.items():
        if (key not in features):
            attributes[key] = value
    return attributes


def comparing(ta, ter, fuzzy_down_ta, fuzzy_up_ta, stats):
    '''main part for comparing terminator and transcript'''
    if (ta.seq_id == ter.seq_id) and (
            ta.strand == ter.strand):
        if ta.strand == "+":
            if ((ta.end >= ter.start) and (
                     ta.end <= ter.end)) or (
                    (ta.end <= ter.start) and (
                     (ter.start - ta.end) <= fuzzy_down_ta)) or (
                    (ta.end >= ter.end) and (
                     (ta.end - ter.end) <= fuzzy_up_ta)) or (
                    (ta.start <= ter.start) and (
                     ta.end >= ter.end)):
                if ter.attributes["Parent"] == "NA":
                    stats[ta.seq_id]["overlap"] += 1
                    ta.attributes["associated_term"] = (
                        "terminator:" + str(ter.start) + "-" +
                        str(ter.end) + "_" + ter.strand)
                    if "ID" in ta.attributes.keys():
                        ter.attributes["Parent"] = ta.attributes["ID"]
                    else:
                        ter.attributes["Parent"] = (
                            "transcript:" + str(ta.start) + "-" +
                            str(ta.end) + "_" + ta.strand)
        else:
            if ((ta.start >= ter.start) and (
                     ta.start <= ter.end)) or (
                    (ta.start <= ter.start) and (
                     (ter.start - ta.start) <= fuzzy_up_ta)) or (
                    (ta.start >= ter.end) and (
                     (ta.start - ter.end) <= fuzzy_down_ta)) or (
                    (ta.start <= ter.start) and (
                     ta.end >= ter.end)):
                if ter.attributes["Parent"] == "NA":
                    stats[ta.seq_id]["overlap"] += 1
                ta.attributes["associated_term"] = (
                    "terminator:" + str(ter.start) + "-" +
                    str(ter.end) + "_" + ter.strand)
                if "ID" in ta.attributes.keys():
                    ter.attributes["Parent"] = ta.attributes["ID"]
                else:
                    ter.attributes["Parent"] = (
                        "transcript:" + str(ta.start) + "-" +
                        str(ta.end) + "_" + ta.strand)


def output_term(ters, term_file, type_, term_outfolder):
    out = open(term_file + "tmp", "w")
    out.write("##gff-version 3\n")
    for ter in ters:
        attribute_string = ";".join(
            ["=".join(items) for items in ter.attributes.items()])
        out.write("\t".join([ter.info_without_attributes,
                             attribute_string]) + "\n")
    out.close()
    os.remove(term_file)
    filename = term_file.split("/")[-1]
    if filename in os.listdir(term_outfolder):
        os.remove(os.path.join(term_outfolder, filename))
    shutil.copy(term_file + "tmp", term_file)
    shutil.move(term_file + "tmp", os.path.join(term_outfolder, filename))
    if type_ == "terminator":
        table_file = term_file.replace("/gffs/", "/tables/")
        table_file = table_file.replace(".gff", ".csv")
        out_t = open(table_file + "tmp", "w")
        out_t.write("\t".join(["Genome", "Name", "Start", "End", "Strand",
                               "Detect", "Associated_gene",
                               "Associated_transcript",
                               "Coverage_decrease", "Coverage_detail"]) + "\n")
        fh = open(table_file, "r")
        for row in csv.reader(fh, delimiter='\t'):
            if row[0] != "genome":
                for ter in ters:
                    if (row[0] == ter.seq_id) and (
                            row[2] == str(ter.start)) and (
                            row[3] == str(ter.end)) and (
                            row[4] == ter.strand):
                        out_t.write("\t".join([row[0], row[1], row[2], row[3],
                                    row[4], row[5], row[6],
                                    ter.attributes["Parent"],
                                    row[7], row[8]]) + "\n")
                        break
        fh.close()
        out_t.close()
        os.remove(table_file)
        shutil.move(table_file + "tmp", table_file)


def read_gff(filename, index):
    gf = open(filename, "r")
    gff_parser = Gff3Parser()
    datas = []
    for entry in gff_parser.entries(gf):
        entry.attributes[index] = "NA"
        datas.append(entry)
        datas = sorted(datas, key=lambda k: (k.seq_id, k.start,
                                             k.end, k.strand))
    gf.close()
    return datas


def compare_term_tran(trans, terms, fuzzy_up_ta, fuzzy_down_ta,
                      out_folder, type_, term_outfolder, tran_outfolder):
    '''Comparison of terminator and transcript. It can realise the 
    relationship of terminator and transcript'''
    for tran in os.listdir(trans):
        if tran.endswith("_transcript.gff"):
            prefix = tran.replace("_transcript.gff", "")
            out_g = open(os.path.join(trans, tran) + "tmp", "w")
            out_g.write("##gff-version 3\n")
            tas = read_gff(os.path.join(trans, tran), "associated_term")
            ters = read_gff(os.path.join(terms, prefix + "_term.gff"),
                            "Parent")
            stats = {}
            pre_seq = ""
            for ta in tas:
                if ta.seq_id != pre_seq:
                    stats[ta.seq_id] = {"all_tran": 0, "all_term": 0,
                                        "overlap": 0}
                    pre_seq = ta.seq_id
                    new_term = True
                stats[ta.seq_id]["all_tran"] += 1
                for ter in ters:
                    if new_term:
                        stats[ta.seq_id]["all_term"] += 1
                    comparing(ta, ter, fuzzy_down_ta, fuzzy_up_ta, stats)
                if new_term:
                    new_term = False
                attribute_string = ";".join(
                    ["=".join(items) for items in ta.attributes.items()])
                out_g.write("\t".join([ta.info_without_attributes,
                                       attribute_string]) + "\n")
            os.remove(os.path.join(trans, tran))
            if tran in os.listdir(tran_outfolder):
                os.remove(os.path.join(tran_outfolder, tran))
            shutil.copy(os.path.join(trans, tran) + "tmp",
                        os.path.join(tran_outfolder, tran))
            shutil.move(os.path.join(trans, tran) + "tmp",
                        os.path.join(trans, tran))
            output_term(ters, os.path.join(terms, prefix + "_term.gff"), type_,
                        term_outfolder)
            out = open(os.path.join(out_folder,
                       "statistics/stat_compare_transcript_terminator_" + prefix + ".csv"), "w")
            for strain, stat in stats.items():
                out.write(strain + ":\n")
                out.write("\tThe overlap between transcripts "
                          "and terminators are {0}\n".format(
                           stat["overlap"]))
                out.write("\tThe overlap percentage of transcripts are {0}\n".format(
                          float(stat["overlap"])/float(stat["all_tran"])))
                out.write("\tThe overlap percentage of terminators are {0}\n".format(
                          float(stat["overlap"])/float(stat["all_term"])))
            out.close()
