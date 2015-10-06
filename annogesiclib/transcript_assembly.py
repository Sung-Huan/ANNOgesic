import os
import sys
import csv
from annogesiclib.parser_wig import WigParser

def check_tex_conds(tracks, libs, texs, check_tex, conds, tex_notex):
    for track in tracks:
        for lib in libs:
            if lib["name"] == track:
                index = "_".join([lib["cond"], lib["type"]])
                if len(texs) != 0:
                    for key, num in texs.items():
                        if track in key:
                            if (texs[key] >= tex_notex) and (
                                key not in check_tex):
                                check_tex.append(key)
                                if index not in conds.keys():
                                    conds[index] = 1
                                else:
                                    conds[index] += 1
                else:
                    if index not in conds.keys():
                        conds[index] = 1
                    else:
                        conds[index] += 1

def detect_hight_toler(cover, height, tmp_covers, tracks):
    if cover["coverage"] > height:
        if tmp_covers["best"] < cover["coverage"]:
            tmp_covers["best"] = cover["coverage"]
        tracks.append(cover["track"])
    else:
        if cover["coverage"] > tmp_covers["toler"]:
            tmp_covers["toler"] = cover["coverage"]

def elongation(covers, height, template_texs, libs, reps,
               tex_notex, strand, trans, strain, tolers):
    first = True
    pre_pos = -1
    check_tex = []
    tracks = []
    conds = {}
    texs = template_texs.copy()
    tmp_covers = {"best": 0, "toler": -1}
    for cover in covers:
        if (pre_pos == cover["pos"]) or first:
            first = False
            detect_hight_toler(cover, height, tmp_covers, tracks)
        else:
            for track in tracks:
                if len(texs) != 0:
                    for key, num in texs.items():
                        if track in key:
                            texs[key] += 1
            check_tex_conds(tracks, libs, texs, check_tex, conds, tex_notex)
            for cond, num in conds.items():
                if ((num >= reps["tex"]) and ("tex" in cond)) or (
                    (num >= reps["frag"]) and ("frag" in cond)):
                    trans[strain].append({
                        "strand": strand, "pos": pre_wig["pos"],
                        "coverage": tmp_covers["best"], "cond": num})
            if (tmp_covers["toler"] != -1):
                tolers.append(tmp_covers["toler"])
            else:
                tolers.append(height + 10)
            tmp_covers = {"best": 0, "toler": -1}
            tracks = []
            conds = {}
            check_tex = []
            texs = template_texs.copy()
            detect_hight_toler(cover, height, tmp_covers, tracks)
        pre_wig = cover
        pre_pos = cover["pos"]
    return tmp_covers["best"], conds, tracks, texs, pre_pos

def transfer_to_tran(wigs, height, libs, template_texs, strand,
                     reps, tex_notex):
    tolers = {}
    trans = {}
    for strain, covers in wigs.items():
        if strain not in trans:
            trans[strain] = []
            tolers[strain] = []
        sort_covers = sorted(covers, key=lambda k: (k["pos"], k["cond"]))
        best_cover, conds, tracks, texs, pos = elongation(sort_covers, height,
                                               template_texs, libs, reps,
                                               tex_notex, strand, trans,
                                               strain, tolers[strain])
    for track in tracks:
        if len(texs) != 0:
            for key, num in texs.items():
                if track in key:
                    texs[key] += 1
    for track in tracks:
        for lib in libs:
            if lib["name"] == track:
                index = "_".join([lib["cond"], lib["type"]])
                if len(texs) != 0:
                    for key, num in texs.items():
                        if (texs[key] >= tex_notex) and (
                            lib["name"] in key):
                            if index not in conds.keys():
                                conds[index] = 1
                            else:
                                conds[index] += 1
                else:
                    if index not in conds.keys():
                        conds[index] = 1
                    else:
                        conds[index] += 1
    for cond, num in conds.items():
        if ((num >= reps["tex"]) and ("tex" in cond)) or (
            (num >= reps["frag"]) and ("frag" in cond)):
            trans[strain].append({"strand": strand, "pos": pos,
                                  "coverage": best_cover, "cond": num})
    return tolers, trans

def print_transctipt(start, end, width, num, high_cover,
                     low_cover, out, strain, strand):
    if (start != -1) and (end != -1) and (
        (end - start) >= width):
        name = '%0*d' % (5, num)
        attribute = ";".join(
                    ["=".join(items) for items in ([("ID", "tran_" + str(num)),
                     ("Name", "Transcript_" + name),
                     ("high_coverage", str(high_cover)),
                     ("low_coverage", str(low_cover))])])
        out.write("\t".join([str(field) for field in [
                  strain, "Transcript", "Transcript", str(start),
                  str(end), ".", strand, ".", attribute]]) + "\n")

def fill_gap_and_print(trans, strand, tolerance, width, out, low_cutoff, tolers):
    for strain, datas in tolers.items():
        num = 0
        for data in datas:
            num += 1
    for strain, covers in trans.items():
        sort_covers = sorted(covers, key=lambda k: k["pos"])
        first = True
        start = -1
        end = -1
        num = 0
        for cover in covers:
            fit = True
            if first:
                first = False
                start = cover["pos"]
                high_cover = cover["coverage"]
                low_cover = cover["coverage"]
            else:
                if (cover["pos"] - pre_cover["pos"]) <= tolerance:
                    if cover["pos"] - pre_cover["pos"] > 1:
                        for toler_strain, toler_datas in tolers.items():
                            if toler_strain == strain:
                                toler_covers = toler_datas[(pre_cover["pos"] - 1): cover["pos"]]
                                for toler_cover in toler_covers:
                                    if (toler_cover < low_cutoff):
                                        fit = False
                                        break
                    if fit:
                        end = cover["pos"]
                        if high_cover < cover["coverage"]:
                            high_cover = cover["coverage"]
                        if low_cover > cover["coverage"]:
                            low_cover = cover["coverage"]
                if ((cover["pos"] - pre_cover["pos"]) > tolerance) or (not fit):
                    if (start != -1) and (end != -1) and (
                        (end - start) >= width):
                        name = '%0*d' % (5, num)
                        attribute = ";".join(
                                    ["=".join(items) for items in ([
                                     ("ID", "tran_" + str(num)),
                                     ("Name", "Transcript_" + name),
                                     ("high_coverage", str(high_cover)),
                                     ("low_coverage", str(low_cover))])])
                        out.write("\t".join([str(field) for field in [
                                  strain, "Transcript", "Transcript", str(start),
                                  str(end), ".", strand, ".", attribute]]) + "\n")
                        num += 1
                    start = cover["pos"]
                    end = -1
                    high_cover = cover["coverage"]
                    low_cover = cover["coverage"]
            pre_cover = cover
        if len(covers) != 0:
            print_transctipt(start, end, width, num, high_cover,
                             low_cover, out, strain, strand)

def read_wig(filename, libs, strand):
    wigs = {}
    tracks = []
    wig_parser = WigParser()
    wig_fh = open(filename)
    for entry in wig_parser.parser(wig_fh, strand):
        if entry.track not in tracks:
            tracks.append(entry.track)
        if entry.strain not in wigs.keys():
            strain = entry.strain
            wigs[entry.strain] = []
        for lib in libs:
            if strand == lib["strand"]:
                for key, value in lib.items():
                    if (key == "name") and (value == entry.track):
                        cond = lib["cond"]
                        track = value
                        wigs[strain].append({
                              "pos": entry.pos, "coverage": entry.coverage,
                              "strand": entry.strand, "track": track,
                              "cond": cond})
                        break
    wig_fh.close()
    return wigs

def assembly(wig_f_file, wig_r_file, height, width, tolerance, low_cutoff,
             wig_folder, tex_notex, input_lib, replicates, out_file):
    libs = []
    conds = []
    out = open(out_file, "w")
    out.write("##gff-version 3\n")
    for lib in input_lib:
        datas = lib.split(":")
        if datas[2] not in conds:
            conds.append(datas[2])
        for wig in os.listdir(wig_folder):
            if wig == datas[0]:
                with open(wig_folder + "/" + wig, "r") as w_h:
                    for line in w_h:
                        line = line.strip()
                        if line.startswith("track"):
                            name = line.split("=")[-1][1:-1]
                            break
        libs.append({"name": name, "type": datas[1], 
                     "cond": datas[2], "rep": datas[3], "strand": datas[4]})
    texs = {}
    for lib1 in libs:
        if lib1["type"] == "frag":
            type_ = "frag"
        elif (lib1["type"] == "tex") or (lib1["type"] == "notex"):
            for lib2 in libs:
                if (lib1["cond"] == lib2["cond"]) and (
                    lib1["rep"] == lib2["rep"]) and (
                    lib1["type"] == "tex") and (
                    lib2["type"] == "notex") and (
                    lib1["strand"] == lib2["strand"]):
                    texs[lib1["name"] + "_" + lib2["name"]] = 0
            type_ = "tex"
        else:
            print("Error:not correct library type (fragmented or Tex treated??)")
            sys.exit()
    wig_fs = read_wig(wig_f_file, libs, "+")
    wig_rs = read_wig(wig_r_file, libs, "-")
    tolers_f, tran_fs = transfer_to_tran(wig_fs, height,
                        libs, texs, "+", replicates, tex_notex)
    tolers_r, tran_rs = transfer_to_tran(wig_rs, height,
                        libs, texs, "-", replicates, tex_notex)
    fill_gap_and_print(tran_fs, "+", tolerance, width, out,
                       low_cutoff, tolers_f)
    fill_gap_and_print(tran_rs, "-", tolerance, width, out,
                       low_cutoff, tolers_r)
    out.close()
