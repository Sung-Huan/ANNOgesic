import os
import csv
import sys
import shutil
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.helper import Helper


def del_attributes(entry, features):
    attributes = {}
    for key, value in entry.attributes.items():
        if (key not in features):
            attributes[key] = value
    entry.attributes = attributes


def read_gffs(gff_files, feature):
    gffs = {}
    if feature == "transcript":
        gffs["transcript"] = []
        gff_f = open(gff_files, "r")
        for entry in Gff3Parser().entries(gff_f):
            gffs["transcript"].append(entry)
        gff_f.close()
        gffs["transcript"] = sorted(
                gffs["transcript"], key=lambda x: (x.seq_id, x.start,
                                                   x.end, x.strand))
    else:
        num = 0
        for gff_file in gff_files:
            gffs[num] = []
            gff_f = open(gff_file, "r")
            for entry in Gff3Parser().entries(gff_f):
                parent = None
                if (entry.feature == "CDS") or (
                        entry.feature == "exon") or (
                        entry.feature == "repeat_unit") or (
                        entry.feature == "tRNA") or (
                        entry.feature == "rRNA") or (
                        entry.feature == "ncRNA"):
                    if "Parent" in entry.attributes.keys():
                        parent = entry.attributes["Parent"]
                del_attributes(entry, ["associated_tran", "parent_tran",
                                       "Parent", "Parent"])
                if parent is not None:
                    entry.attributes["Parent"] = parent
                entry.attributes["print"] = False
                gffs[num].append(entry)
            gff_f.close()
            gffs[num] = sorted(
                gffs[num], key=lambda x: (x.seq_id, x.start,
                                          x.end, x.strand))
            num += 1
    return gffs


def assign_parent(other, tran):
    '''assign the parent transcript to all features'''
    if "Parent" not in other.attributes.keys():
        other.attributes["Parent"] = tran.attributes["ID"]
    else:
        other.attributes["Parent"] = ",".join([
            other.attributes["Parent"], tran.attributes["ID"]])


def compare_tran(tran_gffs, other_gffs, fuzzy_tss, fuzzy_term):
    for tran in tran_gffs["transcript"]:
        for num, others in other_gffs.items():
            for other in others:
                if (tran.seq_id == other.seq_id) and (
                        tran.strand == other.strand) and (
                        other.feature != "source") and (
                        other.feature != "region"):
                    if other.feature == "TSS":
                        if (other.start >= tran.start - fuzzy_tss) and (
                                other.end <= tran.end + fuzzy_tss):
                            assign_parent(other, tran)
                    elif other.feature.lower() == "terminator":
                        if ta.strand == "+":
                            if ((ta.end >= other.start) and (
                                     ta.end <= other.end)) or (
                                    (ta.end <= other.start) and (
                                     (other.start - ta.end) <= fuzzy_term)) or (
                                    (ta.end >= other.end) and (
                                     (ta.end - other.end) <= fuzzy_term)) or (
                                    (ta.start <= other.start) and (
                                     ta.end >= other.end)):
                                assign_parent(other, tran)
                        else:
                            if ((ta.start >= other.start) and (
                                     ta.start <= other.end)) or (
                                    (ta.start <= other.start) and (
                                     (other.start - ta.start) <= fuzzy_term)) or (
                                    (ta.start >= other.end) and (
                                     (ta.start - other.end) <= fuzzy_term)) or (
                                    (ta.start <= other.start) and (
                                     ta.end >= other.end)):
                                assign_parent(other, tran)
                    else:
                        if ((tran.start <= other.start) and (
                                tran.end >= other.end)) or (
                                (tran.start >= other.start) and (
                                tran.end <= other.end)) or (
                                (tran.start <= other.start) and (
                                tran.end >= other.start) and (
                                tran.end <= other.end)) or (
                                (tran.start >= other.start) and (
                                tran.start <= other.end) and (
                                tran.end >= other.end)):
                            assign_parent(other, tran)


def combine_gffs(tran_gffs, other_gffs):
    gffs = []
    o_gffs = []
    s_gffs = []
    if tran_gffs is not None:
        for tran in tran_gffs["transcript"]:
            gffs.append(tran)
            for num, others in other_gffs.items():
                for other in others:
                    if tran.seq_id == other.seq_id:
                        if "Parent" in other.attributes.keys():
                            attributes = {}
                            for key, value in other.attributes.items():
                                if key != "print":
                                    attributes[key] = value
                            if (tran.attributes["ID"] in 
                                    other.attributes["Parent"].split(",")):
                                other.attribute_string = ";".join(
                                    ["=".join(items) for items in 
                                        attributes.items()])
                                other.attributes["print"] = True
                                gffs.append(other)
    for num, others in other_gffs.items():
        for other in others:
            if (other.feature == "source") or (
                    other.feature == "region"):
                s_gffs.append(other)
            if not other.attributes["print"]:
                attributes = {}
                for key, value in other.attributes.items():
                    if key != "print":
                        attributes[key] = value
                other.attribute_string = ";".join(
                    ["=".join(items) for items in attributes.items()])
                o_gffs.append(other)
    return gffs, o_gffs, s_gffs


def print_gff(gffs, o_gffs, s_gffs, output):
    sort_others = sorted(o_gffs, key=lambda x: (x.seq_id, x.start,
                                                x.end, x.strand))
    out = open(output, "w")
    if len(gffs) != 0:
        pre_strain = None
        for gff in gffs:
            if (pre_strain is not None) and (pre_strain != gff.seq_id):
                for other in sort_others:
                    if other.seq_id == pre_strain:
                        if (not other.attributes["print"]):
                            out.write("\t".join([other.info_without_attributes,
                                                 other.attribute_string]) + "\n")
                            other.attributes["print"] = True
            for source in s_gffs:
                if (source.seq_id == gff.seq_id) and (
                        not source.attributes["print"]):
                    out.write("\t".join([source.info_without_attributes,
                                         source.attribute_string]) + "\n")
                    source.attributes["print"] = True
            out.write("\t".join([gff.info_without_attributes,
                                 gff.attribute_string]) + "\n")
            pre_strain = gff.seq_id
        for other in sort_others:
            if other.seq_id == gff.seq_id:
                if (not other.attributes["print"]):
                    out.write("\t".join([other.info_without_attributes,
                                         other.attribute_string]) + "\n")
                    other.attributes["print"] = True
    else:
        for other in sort_others:
            out.write("\t".join([gff.info_without_attributes,
                                 gff.attribute_string]) + "\n")
    out.close()


def run_merge(out_folder, tran, others, fuzzy_term, fuzzy_tss, strain):
    '''merge all features to be one gff file'''
    output = "_".join([strain, "merge_features.gff"])
    if tran is None and others is None:
        print("Error: There is no input file...")
        sys.exit()
    elif (tran is not None) and (others is None):
        shutil.copy(tran, os.path.join(out_folder, output))
    elif others is not None:
        others = others.split(",")
        if (tran is not None):
            tran_gffs = read_gffs(tran, "transcript")
            other_gffs = read_gffs(others, "others")
            compare_tran(tran_gffs, other_gffs, fuzzy_tss, fuzzy_term)
        else:
            other_gffs = read_gffs(others, "others")
        gffs, o_gffs, s_gffs = combine_gffs(tran_gffs, other_gffs)
        print_gff(gffs, o_gffs, s_gffs, output)
