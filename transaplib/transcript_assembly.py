#!/usr/bin/python

import os        
import sys
import csv
import transaplib.parser_wig as par_wig

def check_tex_conds(tracks, libs, texs, check_tex, conds, tex_notex):
    for track in tracks:
        for lib in libs:
            if lib["name"] == track:
                if len(texs) != 0:
                    for key, num in texs.items():
                        if track in key:
                            if (texs[key] >= tex_notex) and \
                               (key not in check_tex):
                                check_tex.append(key)
                                if lib["cond"] not in conds.keys():
                                    conds[lib["cond"]] = 1
                                else:
                                    conds[lib["cond"]] += 1
                else:
                    if lib["cond"] not in conds.keys():
                        conds[lib["cond"]] = 1
                    else:
                        conds[lib["cond"]] += 1

def elongation(covers, height, template_texs, libs, rep, 
               tex_notex, strand, trans, strain):
    first = True
    best_cover = 0
    pre_pos = -1
    check_tex = []
    tracks = []
    conds = {}
    texs = template_texs.copy()
    for cover in covers:
        if (pre_pos == cover["pos"]) or first:
            first = False
            if cover["coverage"] > height:
                if best_cover < cover["coverage"]:
                    best_cover = cover["coverage"]
                tracks.append(cover["track"])
        else:
            for track in tracks:
                if len(texs) != 0:
                    for key, num in texs.items():
                        if track in key:
                            texs[key] += 1
            check_tex_conds(tracks, libs, texs, check_tex, conds, tex_notex)
            for cond, num in conds.items():
                if num >= rep:
                    trans[strain].append({"strand": strand, "pos": pre_wig["pos"], 
                                          "coverage": best_cover, "cond": num})
            best_cover = 0
            tracks = []
            conds = {}
            check_tex = []
            texs = template_texs.copy()
            if cover["coverage"] > height:
                if best_cover < cover["coverage"]:
                    best_cover = cover["coverage"]
                tracks.append(cover["track"])
        pre_wig = cover
        pre_pos = cover["pos"]
    return (best_cover, conds, tracks, texs)

def transfer_to_tran(wigs, trans, height, libs, template_texs, strand,
                     replicate, tex_notex):
    for strain, covers in wigs.items():
        if strain not in trans:
            trans[strain] = []
        sort_covers = sorted(covers, key=lambda k: (k["pos"], k["cond"]))
        datas = elongation(sort_covers, height, template_texs, libs, 
                           replicate, tex_notex, strand, trans, strain)
        best_cover = datas[0]
        conds = datas[1]
        tracks = datas[2]
        texs = datas[3]
    for track in tracks:
        if len(texs) != 0:
            for key, num in texs:
                if track in key:
                    texs[key] += 1
    for track in tracks:
        for lib in libs:
            if lib["name"] == track:
                if len(texs) != 0:
                    if texs[lib["name"]] >= tex_notex:
                        if lib["cond"] not in conds.keys():
                            conds[lib["cond"]] = 1
                        else:
                            conds[lib["cond"]] += 1
    for cond, num in conds.items():
        if num >= replicate_supported:
            trans[strain].append({"strand": strand, "pos": wig["pos"], "coverage": best_cover, "cond": num})

def print_transcript(trans, strand, tolerance, width, out):
    for strain, covers in trans.items():
        sort_covers = sorted(covers, key=lambda k: k["pos"])
        first = True
        start = -1
        end = -1
        num = 0
        for cover in covers:
            if first:
                first = False
                start = cover["pos"]
                high_cover = cover["coverage"]
                low_cover = cover["coverage"]
            else:
                if (cover["pos"] - pre_cover["pos"]) <= tolerance:
                    end = cover["pos"]
                    if high_cover < cover["coverage"]:
                        high_cover = cover["coverage"]
                    if low_cover > cover["coverage"]:
                        low_cover = cover["coverage"]
                else:
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
                        num += 1
                    start = cover["pos"]
                    end = -1
                    high_cover = cover["coverage"]
                    low_cover = cover["coverage"]
            pre_cover = cover
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

def read_wig(wigs, filename, libs, strand):
    wig_parser = par_wig.parser_wig()
    for entry in wig_parser.parser(filename, strand):
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

def Assembly(wig_f_file, wig_r_file, height, width, tolerance,
             wig_folder, tex_notex, input_lib, replicates, out_file):
    wig_fs = {}
    wig_rs = {}
    tran_fs = {}
    tran_rs = {}
    libs = []
    conds = []
    out = open(out_file, "w")
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
    read_wig(wig_fs, wig_f_file, libs, "+")
    read_wig(wig_rs, wig_r_file, libs, "-")
    transfer_to_tran(wig_fs, tran_fs, height, libs, texs, "+", replicates, tex_notex)
    transfer_to_tran(wig_rs, tran_rs, height, libs, texs, "-", replicates, tex_notex)
    print_transcript(tran_fs, "+", tolerance, width, out)
    print_transcript(tran_rs, "-", tolerance, width, out)
