import os
import csv
import math
import copy
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.coverage_detection import coverage_comparison, replicate_comparison
from annogesiclib.lib_reader import read_wig, read_libs

def modify_attributes(pre_srna, srna, srna_type, input_type):
    if srna_type == "UTR":
        if pre_srna.attributes["UTR_type"] != srna.attributes["UTR_type"]:
            if input_type == "pre":
                if "&" not in pre_srna.attributes["UTR_type"]:
                    pre_srna.attributes["UTR_type"] = \
                           "&".join([srna.attributes["UTR_type"],
                                     pre_srna.attributes["UTR_type"]])
            else:
                if "&" not in pre_srna.attributes["UTR_type"]:
                    srna.attributes["UTR_type"] = \
                           "&".join([srna.attributes["UTR_type"],
                                     pre_srna.attributes["UTR_type"]])
                else:
                    srna.attributes["UTR_type"] = pre_srna.attributes["UTR_type"]

def del_attributes(feature, entry):
    attributes = {}
    for key, value in entry.attributes.items():
        if feature not in key:
            attributes[key] = value
    return attributes

def detect_overlap(srna, pre_srna, srna_type, overlap):
    if (srna.seq_id == pre_srna.seq_id) and (
        srna.strand == pre_srna.strand):
        if (pre_srna.start >= srna.start) and (
            pre_srna.end <= srna.end):
            modify_attributes(pre_srna, srna, srna_type, None)
            overlap = True
        elif (pre_srna.start >= srna.start) and (
              pre_srna.start <= srna.end) and (
              pre_srna.end >= srna.end):
            modify_attributes(pre_srna, srna, srna_type, None)
            overlap = True
        elif (pre_srna.start <= srna.start) and (
              pre_srna.end >= srna.start) and (
              pre_srna.end <= srna.end):
            modify_attributes(pre_srna, srna, srna_type, None)
            overlap = True
        elif (pre_srna.start <= srna.start) and \
             (pre_srna.end >= srna.end):
            overlap = True
            modify_attributes(pre_srna, srna, srna_type, "pre")
    return overlap

def merge_tss_pro(pre_srna, srna, feature):
    if (feature not in pre_srna.attributes.keys()) and (
        feature in srna.attributes.keys()):
        if srna.attributes[feature] != "NA":
            pre_srna.attributes[feature] = srna.attributes[feature]
    elif (feature in pre_srna.attributes.keys()) and (
        feature in srna.attributes.keys()):
        if (pre_srna.attributes[feature] == "NA") and (
            srna.attributes[feature] != "NA"):
            pre_srna.attributes[feature] = srna.attributes[feature]
        elif (srna.attributes[feature] not in pre_srna.attributes[feature]) and (
              srna.attributes[feature] != "NA"):
            pre_srna.attributes[feature] = "&".join(
                                              [pre_srna.attributes[feature],
                                               srna.attributes[feature]])

def modify_overlap(pre_srna, srna):
    merge_tss_pro(pre_srna, srna, "with_TSS")
    merge_tss_pro(pre_srna, srna, "end_cleavage")
    if "UTR_type" in srna.attributes.keys():
        merge_tss_pro(pre_srna, srna, "start_cleavage")
    if (srna.start < pre_srna.start):
        pre_srna.start = srna.start
    if (srna.end > pre_srna.end):
        pre_srna.end = srna.end
    return pre_srna

def merge_srna(srnas, srna_type):
    final_srnas = []
    first = True
    for srna in srnas:
        if srna_type == "UTR":
            srna.source = "UTR_derived"
            srna.feature = "sRNA"
        else:
            if srna.source != "in_CDS":
                srna.source = "intergenic"
            if "with_TSS" in srna.attributes.keys():
                if srna.attributes["with_TSS"] == "False":
                    srna.attributes["with_TSS"] = "NA"
            else:
                srna.attributes["with_TSS"] = "NA"
            if "end_cleavage" in srna.attributes.keys():
                if srna.attributes["end_cleavage"] == "False":
                    srna.attributes["end_cleavage"] = "NA"
            else:
                srna.attributes["end_cleavage"] = "NA"
        overlap = False
        if first:
            first = False
            pre_srna = srna
        else:
            overlap = detect_overlap(srna, pre_srna, srna_type, overlap)
            if overlap:
                pre_srna = modify_overlap(pre_srna, srna)
            if not overlap:
                final_srnas.append(pre_srna)
                pre_srna = srna
    if not overlap:
        final_srnas.append(srna)
    return final_srnas

def read_gff(gff_file, type_):
    datas = []
    if os.path.exists(gff_file):
        for entry in Gff3Parser().entries(open(gff_file)):
            if type_ == "sRNA":
                datas.append(entry)
            elif type_ == "tss":
                datas.append(entry)
            else:
                if (entry.feature == "CDS") or (
                    entry.feature == "tRNA") or (
                    entry.feature == "rRNA"):
                    datas.append(entry)
        datas = sorted(datas, key=lambda k: (k.seq_id, k.start))
    return datas

def read_table(table_file, file_type):
    datas = []
    if os.path.exists(table_file):
        f_h = open(table_file, "r")
        for row in csv.reader(f_h, delimiter='\t'):
            datas.append(import_data(row, file_type))
    return datas

def merge_incds_utr(utrs, inters):
    new_inters = []
    for inter in inters:
        remove = False
        for utr in utrs:
            if inter.source == "in_CDS":
                if (inter.seq_id == utr.seq_id) and (
                    inter.strand == utr.strand):
                    if ((inter.end < utr.end) and (
                         inter.end > utr.start) and (
                         inter.start <= utr.start)) or (
                        (inter.start > utr.start) and (
                         inter.start < utr.end) and (
                         inter.end >= utr.end)) or (
                        (inter.end >= utr.end) and (
                         inter.start <= utr.start)) or (
                        (inter.end <= utr.end) and (
                         inter.start >= utr.start)):
                        utr.start = min(inter.start, utr.start)
                        utr.end = max(inter.end, utr.end)
                        remove = True
        if not remove:
            new_inters.append(inter)
    return new_inters

def compare_srna_cds(srna, cdss, cutoff_overlap):
    detect = False
    overlap = False
    for cds in cdss:
        if (srna.seq_id == cds.seq_id) and (
            srna.strand == cds.strand):
            if ((srna.end < cds.end) and (
                 srna.end > cds.start) and (
                 srna.start <= cds.start)) or (
                (srna.start > cds.start) and (
                 srna.start < cds.end) and (
                 srna.end >= cds.end)) or (
                (srna.end >= cds.end) and (
                 srna.start <= cds.start)) or (
                (srna.end <= cds.end) and (
                 srna.start >= cds.start)):
                overlap = True
                per_c = float(min(srna.end, cds.end) - max(
                        srna.start, cds.start)) / float(cds.end - cds.start)
                if per_c <= cutoff_overlap:
                    if "product" in cds.attributes.keys():
                        cds_name = "".join([cds.feature, ":", str(cds.start),
                                            "-", str(cds.end), "_", cds.strand,
                                            "(", cds.attributes["product"], ")"])
                    else:
                        cds_name = "".join([cds.feature, ":", str(cds.start),
                                            "-", str(cds.end), "_", cds.strand])
                    if "overlap_cds" not in srna.attributes.keys():
                        srna.attributes["overlap_cds"] = cds_name
                        srna.attributes["overlap_percent"] = str(per_c)
                    else:
                        srna.attributes["overlap_cds"] = \
                            "&".join([srna.attributes["overlap_cds"], cds_name])
                        srna.attributes["overlap_percent"] = \
                            "&".join([srna.attributes["overlap_percent"], str(per_c)])
                    detect = True
    if not overlap:
        srna.attributes["overlap_cds"] = "NA"
        srna.attributes["overlap_percent"] = "NA"
        return srna
    elif overlap and detect:
        return srna
    else:
        return None
                        
def merge_srna_gff(srna_utr, srna_inter, out_file, in_cds, cutoff_overlap,
                   gff_file):
    out = open(out_file, "w")
    out.write("##gff-version 3\n")
    utrs = read_gff(srna_utr, "sRNA")
    inters = read_gff(srna_inter, "sRNA")
    cdss = read_gff(gff_file, "CDS")
    num_srna = 0
    srnas = None
    if (in_cds) and (len(utrs) != 0) and (len(inters) != 0):
        inters = merge_incds_utr(utrs, inters)
    if len(utrs) != 0:
        srnas = merge_srna(utrs, "UTR")
    if len(inters) != 0:
        if srnas is not None:
            srnas = srnas + merge_srna(inters, "inter")
        else:
            srnas = merge_srna(inters, "inter")
    sort_srnas = sorted(srnas, key=lambda x: (x.seq_id, x.start))
    for srna in sort_srnas:
        new_srna = compare_srna_cds(srna, cdss, cutoff_overlap)
        if new_srna:
            new_srna.attributes["ID"] = "srna" + str(num_srna)
            name = '%0*d' % (5, num_srna)
            new_srna.attributes["Name"] = "sRNA_candidate_" + str(name)
            new_srna.attributes = del_attributes("best_high_coverage", new_srna)
            new_srna.attributes = del_attributes("best_low_coverage", new_srna)
            new_srna.attributes = del_attributes("best_avg_coverage", new_srna)
            attribute_string = ";".join(
                              ["=".join(items) for items in new_srna.attributes.items()])
            new_srna.info_without_attributes = "\t".join([str(field) for field in [
                            new_srna.seq_id, new_srna.source, new_srna.feature,
                            new_srna.start, new_srna.end, new_srna.score,
                            new_srna.strand, new_srna.phase]])
            out.write(srna.info_without_attributes + "\t" + attribute_string + "\n")
            num_srna += 1
    out.close()

def import_data(row, type_):
    if type_ == "inter":
        return {"strain": row[0], "name": row[1],
                "start": int(row[2]), "end": int(row[3]),
                "strand": row[4], "libs": row[5],
                "detect": row[6], "avg": row[7],
                "high": row[8], "low": row[9],
                "detail": row[11], "tss": row[10]}
    if type_ == "utrr":
        return {"strain": row[0], "name": row[1],
                "start": int(row[2]), "end": int(row[3]),
                "strand": row[4], "libs": row[5],
                "detect": row[6], "avg": row[7],
                "high": row[8], "low": row[9],
                "detail": row[10]}

def check_real_cut(inter_cuts, tss_type, cut):
    for tss, value in inter_cuts.items():
        if tss in tss_type.lower():
            if cut is None:
                cut = inter_cuts[tss]
            else:
                if cut > inter_cuts[tss]:
                    cut = inter_cuts[tss]
    if cut is None:
        if "no_tss" not in inter_cuts.keys():
            cut = 0
        else:
            cut = inter_cuts["no_tss"]
    return cut

def get_cutoff(srna, tsss, cutoff_tex, cutoff_frag, type_, fuzzy_tss,
               tables, out_folder):
    if type_ == "inter":
        tss_type = None
        inter_cuts = {"frag": {}, "tex": {}, "notex": {}}
        fh = open(os.path.join(out_folder, "tmp_cutoff_inter"), "r")
        for row in csv.reader(fh, delimiter='\t'):
            inter_cuts[row[0]][row[1]] = float(row[2])
        if tsss is not None:
            for tss in tsss:
                if (srna.seq_id == tss.seq_id) and (
                    srna.strand == tss.strand):
                    if srna.strand == "+":
                        if math.fabs(srna.start - tss.start) <= fuzzy_tss:
                            tss_type = tss.attributes["type"]
                            if srna.start == tss.start:
                                break
                    else:
                        if math.fabs(srna.end - tss.start) <= fuzzy_tss:
                            tss_type = tss.attributes["type"]
                            if srna.end == tss.start:
                                break
            cut = {"frag": None, "tex": None, "notex": None}
            if tss_type is None:
                tss_type = "no_tss"
            for key, types in inter_cuts.items():
                cut[key] = check_real_cut(types, tss_type, cut[key])
        else:
            cut = {}
            for key, types in inter_cuts.items():
                cut[key] = inter_cuts[key]["no_tss"]
    elif type_ == "utr":
        cut = {}
        fh = open(os.path.join(out_folder, "tmp_median"), "r")
        for row in csv.reader(fh, delimiter='\t'):
            if (row[0] == srna.seq_id) and (
                row[1] == srna.attributes["UTR_type"]):
                if row[1] not in cut.keys():
                    cut[row[1]] = {}
                cut[row[1]][row[2]] = {"median": float(row[3])}
    fh.close()
    return cut

def devide_covers(covers):
    frag_covers = {}
    tex_covers = {}
    for cond, tracks in covers.items():
        if "frag" in cond:
            frag_covers[cond] = tracks
        elif "tex" in cond:
            tex_covers[cond] = tracks
    return frag_covers, tex_covers

def merge_srna_datas(srna_datas_tex, srna_datas_frag):
    if (len(srna_datas_tex["conds"]) != 0) and (
        len(srna_datas_frag["conds"]) != 0):
        srna_datas = copy.deepcopy(srna_datas_tex)
        for key, values in srna_datas_frag.items():
            if key == "conds":
                srna_datas["conds"] = dict(srna_datas["conds"], **values)
            elif key == "best":
                if srna_datas["best"] < values:
                    srna_datas["best"] = values
                    srna_datas["high"] = srna_datas_frag["high"]
                    srna_datas["low"] = srna_datas_frag["low"]
                    srna_datas["track"] = srna_datas_frag["track"]
            elif key == "detail":
                srna_datas["detail"] = srna_datas["detail"] + values
    elif len(srna_datas_tex["conds"]) != 0:
        srna_datas = copy.deepcopy(srna_datas_tex)
    elif len(srna_datas_frag["conds"]) != 0:
        srna_datas = copy.deepcopy(srna_datas_frag)
    else:
        srna_datas = copy.deepcopy(srna_datas_tex)
    return srna_datas

def compare_table(srna, tables, type_, wigs_f, wigs_r, template_texs, 
                  table_best, out, tex_notex, replicates, cutoff_tex,
                  cutoff_notex, cutoff_frag, tsss, fuzzy_tss, out_folder):
    detect = False
    tss_pro, end_pro = get_tss_pro(type_, srna)
    if not detect:
        if type_ == "inter":
            if srna.strand == "+":
                covers = get_coverage(wigs_f, srna.seq_id, srna.strand,
                                      srna.start, srna.end)
            else:
                covers = get_coverage(wigs_r, srna.seq_id, srna.strand,
                                      srna.start, srna.end)
            cut = get_cutoff(srna, tsss, cutoff_tex, cutoff_frag, type_,
                             fuzzy_tss, tables, out_folder)
            frag_covers, tex_covers = devide_covers(covers)
            srna_datas_tex = replicate_comparison(tex_covers, template_texs, srna.strand,
                                  cut["tex"], tex_notex, replicates, "normal",
                                  None, None, None, template_texs, cut["notex"])
            srna_datas_frag = replicate_comparison(frag_covers, template_texs, srna.strand,
                                  cut["frag"], tex_notex, replicates, "normal",
                                  None, None, None, template_texs, None)
            srna_datas = merge_srna_datas(srna_datas_tex, srna_datas_frag)
        elif type_ == "utr":
            if srna.strand == "+":
                covers = get_coverage(wigs_f, srna.seq_id, srna.strand,
                                      srna.start, srna.end)
            else:
                covers = get_coverage(wigs_r, srna.seq_id, srna.strand,
                                      srna.start, srna.end)
            cut = get_cutoff(srna, tsss, cutoff_tex, cutoff_frag, type_,
                             fuzzy_tss, tables, out_folder)
            srna_datas = replicate_comparison(covers, template_texs, srna.strand,
                                      cut, tex_notex, replicates,
                                      "sRNA_utr_derived", cut[srna.attributes["UTR_type"]],
                                      cut, srna.attributes["UTR_type"], template_texs, None)
        if len(srna_datas["conds"]) != 0:
            out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t".format(
                       srna.seq_id, srna.attributes["Name"], srna.start,
                       srna.end, srna.strand,
                       ";".join(srna_datas["conds"].keys()),
                       ";".join(srna_datas["conds"].values()), tss_pro, end_pro,
                       srna_datas["best"], srna_datas["high"], srna_datas["low"]))
            if not table_best:
                first = True
                for data in srna_datas["detail"]:
                    if first:
                        out.write("{0}(avg={1};high={2};low={3})".format(
                                  data["track"], data["avg"], data["high"],
                                  data["low"]))
                        first = False
                    else:
                        out.write(";{0}(avg={1};high={2};low={3})".format(
                                  data["track"], data["avg"], data["high"],
                                  data["low"]))
            else:
                out.write("{0}(avg={1};high={2};low={3})".format(
                          srna_datas["track"], srna_datas["best"],
                          srna_datas["high"], srna_datas["low"]))
            out.write("\t{0}\t{1}\n".format(
                      srna.attributes["overlap_cds"].replace("&", ";"),
                      srna.attributes["overlap_percent"].replace("&", ";")))

def get_coverage(wigs, strain, strand, start, end):
    cover_sets = {"high": -1, "low": -1, "total": 0, "diff": 0}
    poss = {"high": 0, "low": 0, "pos": 0}
    srna_covers = {}
    for wig_strain, conds in wigs.items():
        if wig_strain == strain:
            for cond, tracks in conds.items():
                srna_covers[cond] = []
                for track, covers in tracks.items():
                    cover_sets["total"] = 0
                    cover_sets["diff"] = 0
                    first = True
                    if strand == "+":
                        covers = covers[start-2:end+1]
                    elif strand == "-":
                        covers = reversed(covers[start-2:end+1])
                    for cover in covers:
                        if (cover["strand"] == strand):
                            if (start <= cover["pos"]) and (
                                end >= cover["pos"]):
                                cover_sets["total"] = cover_sets["total"] + cover["coverage"]
                                first = coverage_comparison(
                                        cover, cover_sets, poss, first, strand)
                            else:
                                if (strand == "+") and (cover["pos"] > end):
                                    cover_sets["pos"] = cover["pos"]
                                    break
                                elif (strand == "-") and (cover["pos"] < start):
                                    cover_sets["pos"] = cover["pos"]
                                    break
                    avg = cover_sets["total"] / float(end - start + 1)
                    srna_covers[cond].append({"track": track,
                                              "high": cover_sets["high"],
                                              "low": cover_sets["low"],
                                              "avg": avg,
                                              "pos": poss["pos"],
                                              "type": cover["type"],
                                              "final_start": start,
                                              "final_end": end})
    return srna_covers

def get_tss_pro(type_, srna):
    if type_ == "utr":
        if (srna.attributes["with_TSS"] != "NA") and (
            srna.attributes["start_cleavage"] != "NA"):
            tss_pro = ";".join([srna.attributes["with_TSS"],
                                srna.attributes["start_cleavage"]])
        elif (srna.attributes["with_TSS"] != "NA"):
            tss_pro = srna.attributes["with_TSS"]
        elif srna.attributes["start_cleavage"] != "NA":
            tss_pro = srna.attributes["start_cleavage"]
        else:
            tss_pro = "NA"
        if (srna.attributes["end_cleavage"] != "NA"):
            end_pro = srna.attributes["end_cleavage"]
        else:
            end_pro = "NA"
        tss_pro = tss_pro.replace("&", ";")
        end_pro = end_pro.replace("&", ";")
    elif type_ == "inter":
        tss_pro = ""
        end_pro = ""
        if (srna.attributes["with_TSS"] != "NA"):
            tss_pro = srna.attributes["with_TSS"].replace("&", ";")
        else:
            tss_pro = "NA"
        if (srna.attributes["end_cleavage"] != "NA"):
            end_pro = srna.attributes["end_cleavage"].replace("&", ";")
        else:
            end_pro = "NA"
    return tss_pro, end_pro

def merge_srna_table(srna_file, inter_table, utr_table, wig_f_file,
                     wig_r_file, wig_folder, input_libs, tex_notex, replicates,
                     table_best, out_file, in_cds, tss_file, cutoff_tex, cutoff_notex,
                     cutoff_frag, fuzzy_tss, out_folder, utr_tex, utr_notex, utr_frag):
    libs, texs = read_libs(input_libs, wig_folder)
    wigs_f = read_wig(wig_f_file, "+", libs)
    wigs_r = read_wig(wig_r_file, "-", libs)
    srnas = read_gff(srna_file, "sRNA")
    if tss_file is not None:
        tsss = read_gff(tss_file, "tss")
    else:
        tsss = None
    inters = read_table(inter_table, "inter")
    utrs = read_table(utr_table, "utr")
    out = open(out_file, "w")
    for srna in srnas:
        if srna.source == "UTR_derived":
            compare_table(srna, utrs, "utr", wigs_f, wigs_r, texs,
                          table_best, out, tex_notex, replicates, cutoff_tex,
                          cutoff_notex, cutoff_frag, tsss, fuzzy_tss, out_folder)
        elif (srna.source == "intergenic") or (srna.source == "in_CDS"):
            compare_table(srna, inters, "inter", wigs_f, wigs_r, texs,
                          table_best, out, tex_notex, replicates, utr_tex,
                          utr_notex, utr_frag, tsss, fuzzy_tss, out_folder)
