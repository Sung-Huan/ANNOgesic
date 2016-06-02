import os
import sys
import copy
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.lib_reader import read_wig, read_libs
from annogesiclib.plot_coverage_table import plot_table


def read_data(gff, features):
    gffs = {}
    stats = {}
    outs = {}
    for entry in Gff3Parser().entries(open(gff)):
        for feature in features:
            if feature not in gffs.keys():
                gffs[feature] = []
                stats[feature] = {}
                outs[feature] = {"all": [], "least_one": [], "none": []}
            if entry.feature == feature:
                gffs[feature].append(entry)
    for feature in gffs.keys():
        gffs[feature] = sorted(gffs[feature], key=lambda k: (
            k.seq_id, k.start, k.end, k.strand))
    return gffs, stats, outs


def set_cutoff(cond, percent_tex, percent_frag, detects, gff):
    if ("tex" in cond) or ("notex" in cond):
        cutoffs = percent_tex.split("_")
    elif ("frag" in cond):
        cutoffs = percent_frag.split("_")
    if cutoffs[0] == "p":
        cutoff_percent = float(cutoffs[1])
        diff = detects["express"] / (gff.end - gff.start + 1)
    elif cutoffs[0] == "n":
        cutoff_percent = float(cutoffs[1])
        diff = detects["express"]
    elif cutoffs[0] == "all":
        cutoff_percent = 0
        diff = detects["express"]
    else:
        print("Error: Please assign the valid cutoff_overlap!!")
    return diff, cutoff_percent


def detect_express(wigs, gff, cutoff_coverage, detects, percent_tex,
                   percent_frag, texs, cond, tex_notex, track, plots,
                   cover_type, name):
    total = 0
    high = 0
    for wig in wigs[(gff.start - 1): gff.end]:
        total = wig["coverage"] + total
        if wig["coverage"] > high:
            high = wig["coverage"]
        if wig["coverage"] >= cutoff_coverage:
            detects["express"] += 1
    if cover_type == "average":
        plots[cond][name] = float(total) / float(gff.end - gff.start + 1)
    elif cover_type == "high":
        plots[cond][name] = high
    else:
        print("Error: The coverage_type is not correct!!!")
        sys.exit()
    diff, cutoff_percent = set_cutoff(cond, percent_tex, percent_frag,
                                      detects, gff)
    if (diff >= float(cutoff_percent)) and (diff != 0):
        if ("tex" in cond) or ("notex" in cond):
            for key, num in texs.items():
                if track in key:
                    texs[key] += 1
                if texs[key] == tex_notex:
                    detects["track"] += 1
        elif "frag" in cond:
            detects["track"] += 1


def compare_wigs(wigs, gff, tex_notex, template_texs, replicates, stats,
                 outs, plots, cover_type, cutoff_coverage, percent_tex,
                 percent_frag):
    detects = {"cond": 0, "track": 0, "import": False, "express": 0}
    texs = copy.deepcopy(template_texs)
    for strain, conds in wigs.items():
        if gff.seq_id == strain:
            detects["cond"] = 0
            num_conds = 0
            for cond, tracks in conds.items():
                num_conds += 1
                plots[cond] = {}
                if cond not in stats[strain].keys():
                    stats[strain][cond] = 0
                if cond not in stats["total"].keys():
                    stats["total"][cond] = 0
                if cond not in outs.keys():
                    outs[cond] = []
                detects["track"] = 0
                for track, wigs in tracks.items():
                    name = track
                    plots[cond][name] = 0
                    detects["express"] = 0
                    detect_express(wigs, gff, cutoff_coverage, detects,
                                   percent_tex, percent_frag, texs, cond,
                                   tex_notex, track, plots, cover_type, name)
                if ("tex" in cond) or ("notex" in cond):
                    if detects["track"] >= replicates["tex"]:
                        detects["import"] = True
                elif ("frag" in cond):
                    if detects["track"] >= replicates["frag"]:
                        detects["import"] = True
                if detects["import"]:
                    detects["import"] = False
                    stats["total"][cond] += 1
                    stats[strain][cond] += 1
                    outs[cond].append(gff)
                    detects["cond"] += 1
            if detects["cond"] == 0:
                stats["total"]["none"] += 1
                stats[strain]["none"] += 1
                outs["none"].append(gff)
            if detects["cond"] == num_conds:
                stats["total"]["all"] += 1
                stats[strain]["all"] += 1
                outs["all"].append(gff)
            if (detects["cond"] <= num_conds) and (
                  detects["cond"] > 0):
                stats["total"]["least_one"] += 1
                stats[strain]["least_one"] += 1
                outs["least_one"].append(gff)


def print_stat(out, stats):
    out.write("\t".join(["total input:", str(stats["total"])]) + "\n")
    for cond, num in stats.items():
        if cond == "least_one":
            tag = "expression at lease one condition:"
        elif cond == "all":
            tag = "expression at all conditions:"
        elif cond == "none":
            tag = "no expression:"
        else:
            tag = "condition " + cond + ":"
        if cond != "total":
            per = "(" + str(float(num) / float(stats["total"])) + ")"
            out.write("\t".join([tag, " ".join([str(num), per])]) + "\n")


def output_stat(stats, stat_folder, prefix):
    for feature, strains in stats.items():
        out = open(os.path.join(stat_folder,
                   "_".join([prefix, feature + ".csv"])), "w")
        if len(strains.keys()) > 2:
            out.write("All strain:\n")
            print_stat(out, strains["total"])
        for strain, stat in strains.items():
            if strain != "total":
                out.write(strain + ":\n")
                print_stat(out, stat)
        out.close()


def output_gff(outs, out_gff_folder, prefix):
    for feature, conds in outs.items():
        for cond, gffs in conds.items():
            if cond == "least_one":
                out = open(os.path.join(
                    out_gff_folder, "_".join([
                        prefix, feature, "at_least_one_lib.gff"])), "w")
            elif cond == "all":
                out = open(os.path.join(
                    out_gff_folder, "_".join([
                        prefix, feature, "all_libs.gff"])), "w")
            elif cond == "none":
                out = open(os.path.join(
                    out_gff_folder, "_".join([
                        prefix, feature, "no_express.gff"])), "w")
            else:
                out = open(os.path.join(
                    out_gff_folder, "_".join([
                        prefix, feature, cond + ".gff"])), "w")
            out.write("##gff-version 3\n")
            for gff in gffs:
                out.write(gff.info + "\n")
            out.close()


def deal_repeat_tag(gff, plots, feature, repeat, tag, tags):
    if (gff.attributes[tag] in tags) and (
            gff.attributes[tag] not in repeat.keys()):
        plots[feature].append({gff.attributes[tag] + "_2": {}})
        repeat[gff.attributes[tag]] = 2
        name = gff.attributes[tag] + "_2"
    elif (gff.attributes[tag] in tags) and (
            gff.attributes[tag] in repeat.keys()):
            plots[feature].append({"_".join([gff.attributes[tag],
                                   str(repeat[gff.attributes[tag]] + 1)]): {}})
            name = "_".join([gff.attributes[tag],
                             str(repeat[gff.attributes[tag]] + 1)])
            repeat[gff.attributes[tag]] += 1
    else:
        plots[feature].append({gff.attributes[tag]: {}})
        name = gff.attributes[tag]
    return name


def get_name(plots, gff, feature, repeat, tags):
    name = "".join([gff.feature, ":", str(gff.start),
                    "-", str(gff.end), "_", gff.strand])
    if feature == "gene":
        if "locus_tag" in gff.attributes.keys():
            name = deal_repeat_tag(gff, plots, feature,
                                   repeat, "locus_tag", tags)
        else:
            plots[feature].append({name: {}})
    elif feature == "CDS":
        if "locus_tag" in gff.attributes.keys():
            name = deal_repeat_tag(gff, plots, feature,
                                   repeat, "locus_tag", tags)
        elif "protein_id" in gff.attributes.keys():
            name = deal_repeat_tag(gff, plots, feature,
                                   repeat, "protein_id", tags)
        else:
            plots[feature].append({name: {}})
    else:
        plots[feature].append({name: {}})
    tags.append(name)
    return name


def plot(plots, stat_folder, max_color, min_color, cover_type):
    for feature in plots:
        plot_table(plots[feature], max_color, min_color,
                   os.path.join(stat_folder, "_".join([
                       feature, cover_type, "express_analysis.png"])))


def gene_expression(input_libs, gff_folder, percent_tex, percent_frag,
                    wig_f_file, wig_r_file, features, wigs, cutoff_coverage,
                    tex_notex, replicates, stat_folder, out_gff_folder,
                    cover_type, max_color, min_color):
    print("Loading wiggle file...")
    libs, texs = read_libs(input_libs, wigs)
    wig_fs = read_wig(wig_f_file, "+", libs)
    wig_rs = read_wig(wig_r_file, "-", libs)
    plots = {}
    repeat = {}
    for gff in os.listdir(gff_folder):
        if gff.endswith(".gff"):
            prefix = gff.replace(".gff", "")
            print("Computing " + prefix)
            gff_list, stats, outs = read_data(os.path.join(gff_folder, gff),
                                              features)
            for feature, gffs in gff_list.items():
                plots[feature] = []
                repeat[feature] = {}
                tags = []
                stats[feature]["total"] = {"total": 0, "least_one": 0,
                                           "all": 0, "none": 0}
                num = 0
                for gff in gffs:
                    if gff.seq_id not in stats[feature].keys():
                        stats[feature][gff.seq_id] = {
                                "total": 0, "least_one": 0,
                                "all": 0, "none": 0}
                    stats[feature]["total"]["total"] += 1
                    stats[feature][gff.seq_id]["total"] += 1
                    name = get_name(plots, gff, feature, repeat[feature], tags)
                    if gff.strand == "+":
                        compare_wigs(
                                wig_fs, gff, tex_notex, texs, replicates,
                                stats[feature], outs[feature],
                                plots[feature][num][name], cover_type,
                                cutoff_coverage, percent_tex, percent_frag)
                    elif gff.strand == "-":
                        compare_wigs(
                                wig_rs, gff, tex_notex, texs, replicates,
                                stats[feature], outs[feature],
                                plots[feature][num][name], cover_type,
                                cutoff_coverage, percent_tex, percent_frag)
                    num += 1
            output_stat(stats, stat_folder, prefix)
            output_gff(outs, out_gff_folder, prefix)
    plot(plots, stat_folder, max_color, min_color, cover_type)
