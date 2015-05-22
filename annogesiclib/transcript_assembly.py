#!/usr/bin/python

import os        
import sys
import csv
import annogesiclib.parser_wig as par_wig

def check_tex_conds(tracks, libs, texs, check_tex, conds, tex_notex):
    for track in tracks:
        for lib in libs:
            if lib["name"] == track:
                index = "_".join([lib["cond"], lib["type"]])
                if len(texs) != 0:
                    for key, num in texs.items():
                        if track in key:
                            if (texs[key] >= tex_notex) and \
                               (key not in check_tex):
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
#        print(tmp_covers["toler"])

def elongation(covers, height, template_texs, libs, reps, 
               tex_notex, strand, trans, strain, tolers):
    first = True
    pre_pos = -1
    check_tex = []
    tracks = []
    conds = {}
    texs = template_texs.copy()
#    conti = False
    tmp_covers = {"best": 0, "toler": -1}
    for cover in covers:
#        print(cover)
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
                if ((num >= reps["tex"]) and ("tex" in cond)) or \
                   ((num >= reps["frag"]) and ("frag" in cond)):
                    trans[strain].append({"strand": strand, "pos": pre_wig["pos"], 
                                          "coverage": tmp_covers["best"], "cond": num})
            if (tmp_covers["toler"] != -1):
                tolers.append(tmp_covers["toler"])
            else:
                tolers.append(height + 10)
#                    conti = False
#            print(conti)
#            if (tmp_covers["toler"] != -1) and (not conti):
#                tolers[pre_pos] = {"coverage": tmp_covers["toler"]}
#                tmps = {"pos": pre_pos, "cover": tmp_covers["toler"]}
#                conti = True
#            elif (tmp_covers["toler"] != -1) and (conti):
#                if tmps["cover"] < tmp_covers["toler"]:
#                    tolers[tmps["pos"]] = {"coverage": tmp_covers["toler"]}
#            print(tmps["pos"])
#            print(tolers[tmps["pos"]])
            tmp_covers = {"best": 0, "toler": -1}
            tracks = []
            conds = {}
            check_tex = []
            texs = template_texs.copy()
            detect_hight_toler(cover, height, tmp_covers, tracks)
        pre_wig = cover
        pre_pos = cover["pos"]
    return (tmp_covers["best"], conds, tracks, texs, pre_pos)

def transfer_to_tran(wigs, trans, height, libs, template_texs, strand,
                     reps, tex_notex):
    tolers = {}
    for strain, covers in wigs.items():
        if strain not in trans:
            trans[strain] = []
            tolers[strain] = []
        sort_covers = sorted(covers, key=lambda k: (k["pos"], k["cond"]))
        datas = elongation(sort_covers, height, template_texs, libs, 
                           reps, tex_notex, strand, trans, strain, tolers[strain])
        best_cover = datas[0]
        conds = datas[1]
        tracks = datas[2]
        texs = datas[3]
        pos = datas[4]
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
                        if (texs[key] >= tex_notex) and \
                           (lib["name"] in key):
                            if index not in conds.keys():
                                conds[index] = 1
                            else:
                                conds[index] += 1
    for cond, num in conds.items():
        if ((num >= reps["tex"]) and ("tex" in cond)) or \
           ((num >= reps["frag"]) and ("frag" in cond)):
            trans[strain].append({"strand": strand, "pos": pos, "coverage": best_cover, "cond": num})
    return tolers

def print_transctipt(start, end, width, num, high_cover, low_cover, out, strain, strand):
    if (start != -1) and (end != -1) and \
       ((end - start) >= width):
        name='%0*d' % (5, num)
        attribute = ";".join(
                    ["=".join(items) for items in ([("ID", "tran_" + str(num)),
                     ("Name", "Transcript_" + name), ("high_coverage", str(high_cover)),
                     ("low_coverage", str(low_cover))])])
        out.write("\t".join([str(field) for field in [
                  strain, "Transcript", "Transcript", str(start),
                  str(end), ".", strand, ".", attribute]]) + "\n")

def fill_gap_and_print(trans, strand, tolerance, width, out, low_cutoff, num_track, tolers):
#    print("XXXXXXXXXx")
    for strain, datas in tolers.items():
        num = 0
        for data in datas:
            print(num)
            print(data)
            num += 1
    print("WWW")
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
                        print(cover)
                        print(pre_cover)
                        for toler_strain, toler_datas in tolers.items():
                            if toler_strain == strain:
                                toler_covers = toler_datas[(pre_cover["pos"] - 1): cover["pos"]]
                                for toler_cover in toler_covers:
                                    print(toler_cover)
                                    if (toler_cover < low_cutoff):
                                        print("AAAAA")
                                        fit = False
                                        break
#                                if toler_datas[pre_cover["pos"]]["coverage"] >= low_cutoff:
#                                    fit = True
                    if fit:
                        end = cover["pos"]
                        if high_cover < cover["coverage"]:
                            high_cover = cover["coverage"]
                        if low_cover > cover["coverage"]:
                            low_cover = cover["coverage"]
                if ((cover["pos"] - pre_cover["pos"]) > tolerance) or (not fit):
                    print("BBBB")
                    if (start != -1) and (end != -1) and ((end - start) >= width):
                        name='%0*d' % (5, num)
                        attribute = ";".join(
                                    ["=".join(items) for items in ([("ID", "tran_" + str(num)), 
                                     ("Name", "Transcript_" + name), ("high_coverage", str(high_cover)),
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
        print_transctipt(start, end, width, num, high_cover, low_cover, out, strain, strand)

def read_wig(wigs, filename, libs, strand):
    tracks = []
    wig_parser = par_wig.parser_wig()
    for entry in wig_parser.parser(filename, strand):
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
                              "strand": entry.strand, "track": track, "cond": cond})
                        break
    return len(tracks)

def assembly(wig_f_file, wig_r_file, height, width, tolerance, low_cutoff,
             wig_folder, tex_notex, input_lib, replicates, out_file):
    wig_fs = {}
    wig_rs = {}
    tran_fs = {}
    tran_rs = {}
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
                if (lib1["cond"] == lib2["cond"]) and \
                   (lib1["rep"] == lib2["rep"]) and \
                   (lib1["type"] == "tex") and \
                   (lib2["type"] == "notex") and \
                   (lib1["strand"] == lib2["strand"]):
                    texs[lib1["name"] + "_" + lib2["name"]] = 0
            type_ = "tex"
        else:
            print("Error:not correct library type (fragmented or Tex treated??)")
            sys.exit()
    num_track = read_wig(wig_fs, wig_f_file, libs, "+")
    num_track = read_wig(wig_rs, wig_r_file, libs, "-")
    tolers_f = transfer_to_tran(wig_fs, tran_fs, height, libs, texs, "+", replicates, tex_notex)
    tolers_r = transfer_to_tran(wig_rs, tran_rs, height, libs, texs, "-", replicates, tex_notex)
    fill_gap_and_print(tran_fs, "+", tolerance, width, out, low_cutoff, num_track, tolers_f)
    fill_gap_and_print(tran_rs, "-", tolerance, width, out, low_cutoff, num_track, tolers_r)
