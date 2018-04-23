import os
import csv
import sys
import shutil
from glob import glob
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
        for files in gff_files:
            for gff_file in glob(files):
                gffs[num] = []
                gff_f = open(gff_file, "r")
                for entry in Gff3Parser().entries(gff_f):
                    parent = None
                    if (entry.feature != "gene") and (
                            entry.feature != "transcript") and (
                            entry.feature != "source") and (
                            entry.feature != "region") and (
                            entry.feature != "remark"):
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
        if tran.attributes["ID"] not in other.attributes["Parent"].split(","):
            other.attributes["Parent"] = ",".join([
                other.attributes["Parent"], tran.attributes["ID"]])


def compare_tran(tran_gffs, other_gffs, fuzzy_tss, fuzzy_term):
    for tran in tran_gffs["transcript"]:
        for num, others in other_gffs.items():
            for other in others:
                if (tran.seq_id == other.seq_id) and (
                        tran.strand == other.strand) and (
                        other.feature != "source") and (
                        other.feature != "region") and (
                        other.feature != "remark"):
                    if other.feature == "TSS":
                        if (other.start >= tran.start - fuzzy_tss) and (
                                other.end <= tran.end + fuzzy_tss):
                            assign_parent(other, tran)
                    elif other.feature.lower() == "terminator":
                        if tran.strand == "+":
                            if ((tran.end >= other.start) and (
                                     tran.end <= other.end)) or (
                                    (tran.end <= other.start) and (
                                     (other.start - tran.end) <= fuzzy_term)) or (
                                    (tran.end >= other.end) and (
                                     (tran.end - other.end) <= fuzzy_term)) or (
                                    (tran.start <= other.start) and (
                                     tran.end >= other.end)):
                                assign_parent(other, tran)
                        else:
                            if ((tran.start >= other.start) and (
                                     tran.start <= other.end)) or (
                                    (tran.start <= other.start) and (
                                     (other.start - tran.start) <= fuzzy_term)) or (
                                    (tran.start >= other.end) and (
                                     (tran.start - other.end) <= fuzzy_term)) or (
                                    (tran.start <= other.start) and (
                                     tran.end >= other.end)):
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
                    if not other.attributes["print"]:
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
                    other.feature == "region") or (
                    other.feature == "remark"):
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
    out.write("##gff-version 3\n")
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


def run_merge(out_folder, tran, others, fuzzy_term, fuzzy_tss, strain, log):
    '''merge all features to be one gff file'''
    output = "_".join([strain, "merge_features.gff"])
    if tran is None and others is None:
        lgo.write("No input files are found.\n")
        print("Error: There is no input file...")
        sys.exit()
    elif (tran is not None) and (others is None):
        shutil.copy(tran, os.path.join(out_folder, output))
    elif others is not None:
        if (tran is not None):
            tran_gffs = read_gffs(tran, "transcript")
            other_gffs = read_gffs(others, "others")
            log.write("Comparing transripts and other features to get "
                      "parental transcripts.\n")
            compare_tran(tran_gffs, other_gffs, fuzzy_tss, fuzzy_term)
        else:
            other_gffs = read_gffs(others, "others")
        log.write("Combining all the gff files and merge the features.\n")
        gffs, o_gffs, s_gffs = combine_gffs(tran_gffs, other_gffs)
        print_gff(gffs, o_gffs, s_gffs, output)
        log.write("\t" + output + " is generated.\n")
