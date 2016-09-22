import math
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.lib_reader import read_libs, read_wig


def read_gff(input_file):
    datas = []
    gff_parser = Gff3Parser()
    f_h = open(input_file, "r")
    for entry in gff_parser.entries(f_h):
        entry.attributes["print"] = False
        datas.append(entry)
    datas = sorted(datas, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    return datas


def get_coverage(tar, wigs):
    '''get coverage'''
    coverage = 0
    for strain, conds in wigs.items():
        if tar.seq_id == strain:
            for tracks in conds.values():
                for wigs in tracks.values():
                    if coverage < wigs[tar.start - 1]["coverage"]:
                        coverage = wigs[tar.start - 1]["coverage"]
    return coverage


def compare_wig(tars, wig_fs, wig_rs):
    '''get the coverage of TSS for comparison'''
    for tar in tars:
        if tar.strand == "+":
            tar.attributes["coverage"] = get_coverage(tar, wig_fs)
        elif tar.strand == "-":
            tar.attributes["coverage"] = get_coverage(tar, wig_rs)


def stat(tars, refs, cutoff, gene_length, cluster):
    '''do statistics and print it out'''
    stats = {"tp": 0, "fp": 0, "miss": 0, "fp_rate": 0,
             "tp_rate": 0, "miss_rate": 0}
    num_ref = 0
    for ref in refs:
        num_ref += 1
        detect = False
        for tar in tars:
            if (ref.seq_id == tar.seq_id) and (
                    ref.strand == tar.strand) and (
                    float(tar.attributes["coverage"]) >= cutoff) and (
                    tar.start <= int(gene_length)):
                if math.fabs(ref.start - tar.start) <= cluster:
                    stats["tp"] += 1
                    tar.attributes["print"] = True
                    detect = True
        if not detect:
            stats["miss"] += 1
    for tar in tars:
        if (not tar.attributes["print"]) and (
                float(tar.attributes["coverage"]) >= cutoff) and (
                tar.start <= int(gene_length)):
            stats["fp"] += 1
    stats["fp_rate"] = float(stats["fp"]) / float(int(gene_length) - num_ref)
    stats["tp_rate"] = float(stats["tp"]) / float(num_ref)
    stats["miss_rate"] = float(stats["miss"]) / float(num_ref)
    return stats, num_ref


def print_file(tars, cutoff, out_file):
    out = open(out_file, "w")
    for tar in tars:
        if tar.attributes["coverage"] >= cutoff:
            out.write(tar.info + "\n")


def change_best(num_ref, best, stat_value):
    '''scoring function for evaluate the change of TSS candidates'''
    change = False
    if num_ref > 100:
        if best["tp_rate"] - stat_value["tp_rate"] >= 0.1:
            change = False
        else:
            if (best["tp_rate"] <= stat_value["tp_rate"]) and (
                    best["fp_rate"] >= stat_value["fp_rate"]):
                best = stat_value.copy()
                change = True
            elif (stat_value["tp_rate"] - best["tp_rate"] >= 0.01) and (
                    stat_value["fp_rate"] - best["fp_rate"] <= 0.00005):
                best = stat_value.copy()
                change = True
            elif (best["tp_rate"] - stat_value["tp_rate"] <= 0.01) and (
                    best["fp_rate"] - stat_value["fp_rate"] >= 0.00005):
                best = stat_value.copy()
                change = True
    else:
        if best["tp"] - stat_value["tp"] >= 5:
            change = False
        else:
            if (best["tp"] <= stat_value["tp"]) and (
                    best["fp"] >= stat_value["fp"]):
                best = stat_value.copy()
                change = True
            tp_diff = float(best["tp"] - stat_value["tp"])
            if tp_diff > 0:
                if float(best["fp"] - stat_value["fp"]) >= 5 * tp_diff:
                    best = stat_value.copy()
                    change = True
            elif tp_diff < 0:
                tp_diff = tp_diff * -1
                if float(stat_value["fp"] - best["fp"]) <= 5 * tp_diff:
                    best = stat_value.copy()
                    change = True
    return best, change


def filter_low_expression(gff_file, args_tss, wig_f_file,
                          wig_r_file, out_file):
    '''filter the low expressed TSS'''
    tars = read_gff(gff_file)
    refs = read_gff(args_tss.manual_file)
    libs, texs = read_libs(args_tss.input_lib, args_tss.wig_folder)
    wig_fs = read_wig(wig_f_file, "+", args_tss.libs)
    wig_rs = read_wig(wig_r_file, "-", args_tss.libs)
    compare_wig(tars, wig_fs, wig_rs)
    cutoff = 1
    first = True
    while True:
        stat_value, num_ref = stat(tars, refs, cutoff,
                                   args_tss.gene_length, args_tss.cluster)
        if first:
            first = False
            best = stat_value.copy()
            continue
        else:
            best, change = change_best(num_ref, best, stat_value)
            if not change:
                break
        cutoff = cutoff + 0.1
    print_file(tars, cutoff, out_file)
    return cutoff
