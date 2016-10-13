import os
import sys
import numpy as np
from annogesiclib.lib_reader import read_wig, read_libs
from annogesiclib.coverage_detection import get_repmatch


def check_tex_conds(tracks, libs, texs, check_tex, conds, tex_notex):
    for track in tracks:
        for lib in libs:
            if lib["name"] == track:
                if "tex" in lib["type"]:
                    type_ = "tex"
                else:
                    type_ = "frag"
                index = "_".join([lib["cond"]])
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


def detect_hight_toler(cover, height, tmp_covers, tracks, lib_track):
    if cover > height:
        if tmp_covers["best"] < cover:
            tmp_covers["best"] = cover
        tracks.append(lib_track)
    else:
        if cover > tmp_covers["toler"]:
            tmp_covers["toler"] = cover


def get_repmatch(replicates, cond):
    if "all" in replicates:
        rep = int(replicates.split("_")[-1])
    else:
        for match in replicates.split(","):
            if cond.split("_")[0] == match.split("_")[0]:
                rep = int(match.split("_")[-1])
    return rep

def elongation(lib_conds, template_texs, libs, strand, trans,
               args_tran, strain, tolers):
    '''check coverage and replicate match to form transcript'''
    first = True
    pre_pos = -1
    check_tex = []
    tracks = []
    conds = {}
    pre_wig = None
    detect = False
    texs = template_texs.copy()
    tmp_covers = {"best": 0, "toler": -1}
    for cond, lib_tracks in lib_conds.items():
        for lib_name, covers in lib_tracks.items():
            index_pos = 0
            for cover in covers:
                for cond, lib_tracks in lib_conds.items():
                    for lib_track in lib_tracks.keys():
                        real_track = lib_track.split("|")[-3]
                        if index_pos < len(lib_tracks[lib_track]):
                            compare_cover = lib_tracks[lib_track][index_pos]
                        else:
                            compare_cover = 0
                        detect_hight_toler(
                                compare_cover, args_tran.height,
                                tmp_covers, tracks, real_track)
                for track in tracks:
                    if len(texs) != 0:
                        for key, num in texs.items():
                            if track in key:
                                texs[key] += 1
                check_tex_conds(tracks, libs, texs, check_tex,
                                conds, args_tran.tex)
                for cond, detect_num in conds.items():
                    if ("tex" in cond):
                        tex_rep = get_repmatch(args_tran.replicates["tex"], cond)
                        if detect_num >= tex_rep:
                            detect = True
                    elif ("frag" in cond):
                        frag_rep = get_repmatch(args_tran.replicates["frag"], cond)
                        if detect_num >= frag_rep:
                            detect = True
                if detect:
                    detect = False
                    trans[strain].append(tmp_covers["best"])
                else:
                    trans[strain].append(-1)
                if (tmp_covers["toler"] != -1):
                    tolers.append(tmp_covers["toler"])
                else:
                    tolers.append(args_tran.height + 10)
                tmp_covers = {"best": 0, "toler": -1}
                tracks = []
                conds = {}
                check_tex = []
                texs = template_texs.copy()
                index_pos += 1
            break
        break


def transfer_to_tran(wigs, libs, template_texs, strand, args_tran):
    '''check coverage and replicate match to form transcript'''
    tolers = {}
    trans = {}
    detect = False
    for strain, lib_conds in wigs.items():
        if strain not in trans:
            trans[strain] = []
            tolers[strain] = []
        elongation(lib_conds, template_texs, libs, strand, trans,
                   args_tran, strain, tolers[strain])
    return tolers, trans


def print_transctipt(start, end, width, num, high_cover, wig_type,
                     low_cover, out, strain, strand):
    if (start != -1) and (end != -1) and (
            (end - start) >= width):
        name = '%0*d' % (5, num)
        attribute = gen_attribute_string(num, name, high_cover,
                                         low_cover, wig_type)
        out.write("\t".join([str(field) for field in [
                  strain, "ANNOgesic", "transcript", str(start),
                  str(end), ".", strand, ".", attribute]]) + "\n")


def gen_attribute_string(num, name, high_cover, low_cover, wig_type):
    attribute = ";".join(
                ["=".join(items) for items in ([
                    ("ID", "tran_" + str(num)),
                    ("Name", "transcript_" + name),
                    ("high_coverage", str(high_cover)),
                    ("low_coverage", str(low_cover)),
                    ("detect_lib", wig_type)])])
    return attribute


def fill_gap_and_print(trans, strand, out, tolers, wig_type, args_tran):
    '''compare transcript with CDS to modify transcript(merge mutliple 
    transcript based on overlap with the same CDS)'''
    for strain, datas in tolers.items():
        num = 0
        for data in datas:
            num += 1
    for strain, covers in trans.items():
        first = True
        start = -1
        end = -1
        num = 0
        pre_cover = None
        cover_pos = 1
        for cover in covers:
            fit = True
            if cover != -1:
                if first:
                    first = False
                    start = cover_pos
                    high_cover = cover
                    low_cover = cover
                else:
                    if (cover_pos - pre_pos) <= args_tran.tolerance:
                        if cover_pos - pre_pos > 1:
                            for toler_strain, toler_datas in tolers.items():
                                if toler_strain == strain:
                                    toler_covers = toler_datas[
                                            (pre_pos - 1): cover_pos]
                                    for toler_cover in toler_covers:
                                        if (toler_cover < args_tran.low_cutoff):
                                            fit = False
                                            break
                        if fit:
                            end = cover_pos
                            if high_cover < cover:
                                high_cover = cover
                            if low_cover > cover:
                                low_cover = cover
                    if ((cover_pos - pre_pos) >
                            args_tran.tolerance) or (not fit):
                        if (start != -1) and (end != -1) and (
                                (end - start) >= args_tran.width):
                            name = '%0*d' % (5, num)
                            attribute = gen_attribute_string(num, name, high_cover,
                                                             low_cover, wig_type)
                            out.write("\t".join([str(field) for field in [
                                strain, "ANNOgesic", "transcript", str(start),
                                str(end), ".", strand, ".", attribute]]) + "\n")
                            num += 1
                        start = cover_pos
                        end = -1
                        high_cover = cover
                        low_cover = cover
                pre_cover = cover
                pre_pos = cover_pos
            cover_pos += 1
        if (len(covers) != 0) and (not first) and (
                (start != -1) and (end != -1) and (
                (end - start) >= args_tran.width)):
            print_transctipt(start, end, args_tran.width, num, high_cover,
                             wig_type, low_cover, out, strain, strand)


def assembly(wig_f_file, wig_r_file, wig_folder, input_lib,
             out_file, wig_type, args_tran):
    out = open(out_file, "w")
    out.write("##gff-version 3\n")
    libs, texs = read_libs(input_lib, wig_folder)
    wig_fs = read_wig(wig_f_file, "+", libs)
    wig_rs = read_wig(wig_r_file, "-", libs)
    tolers_f, tran_fs = transfer_to_tran(wig_fs, libs, texs, "+", args_tran)
    tolers_r, tran_rs = transfer_to_tran(wig_rs, libs, texs, "-", args_tran)
    fill_gap_and_print(tran_fs, "+", out, tolers_f, wig_type, args_tran)
    fill_gap_and_print(tran_rs, "-", out, tolers_r, wig_type, args_tran)
    out.close()
    del wig_fs
    del wig_rs
