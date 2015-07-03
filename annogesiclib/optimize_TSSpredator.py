#!/usr/bin/python

import os
import sys
import random
import csv
from subprocess import call, Popen
from multiprocessing import Pool
import math
import time
import signal
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.converter import Converter
import copy

def compute_stat(stat_value, best, best_para, cores, list_num, out_path, indexs):
    if indexs["change"]:
        indexs["change"] = False
        best = stat_value
        best_para = copy.deepcopy(list_num[-1 * cores + indexs["count"]])
    print("_".join(["Current Parameter:step={0}", "height={1}",
                    "reduction_height={2}", "factor={3}", "reduction_factor={4}",
                    "base_height={5}", "enrichment_factor={6}",
                    "processing_factor={7}"]).format(
          indexs["step"] - cores + 1  + indexs["count"],
          list_num[-1 * cores + indexs["count"]]["height"],
          list_num[-1 * cores + indexs["count"]]["re_height"],
          list_num[-1 * cores + indexs["count"]]["factor"],
          list_num[-1 * cores + indexs["count"]]["re_factor"],
          list_num[-1 * cores + indexs["count"]]["base_height"],
          list_num[-1 * cores + indexs["count"]]["enrichment"],
          list_num[-1 * cores + indexs["count"]]["processing"]))
    print("Current:TP={0}\tTP_rate={1}\tFP={2}\tFP_rate={3}\tFN={4}\tMissing_ratio={5}".format(
           stat_value["tp"], stat_value["tp_rate"], stat_value["fp"],
           stat_value["fp_rate"], stat_value["fn"], stat_value["missing_ratio"]))
    print("\t".join(["Best Parameter:height={0}", "reduction_height={1}",
                     "factor={2}", "reduction_factor={3}", "base_height={4}",
                     "enrichment_factor={5}", "processing_factor={6}"]).format(
          best_para["height"],
          best_para["re_height"],
          best_para["factor"],
          best_para["re_factor"],
          best_para["base_height"],
          best_para["enrichment"],
          best_para["processing"]))
    print("Best:TP={0}\tTP_rate={1}\tFP={2}\tFP_rate={3}\tFN={4}\tMissing_ratio={5}".format(
          best["tp"], best["tp_rate"], best["fp"], best["fp_rate"],
          best["fn"], best["missing_ratio"]))
    best_out = open(out_path + "/best.csv", "w")
    para_line = "_".join(["he", str(best_para["height"]),
                          "rh", str(best_para["re_height"]),
                          "fa", str(best_para["factor"]),
                          "rf", str(best_para["re_factor"]),
                          "bh", str(best_para["base_height"]),
                          "ef", str(best_para["enrichment"]),
                          "pf", str(best_para["processing"])])
    best_out.write("{0}\tTP={1}\tTP_rate={2}\tFP={3}\tFP_rate={4}\tFN={5}\tMissing_ratio={6}\t".format(
                   para_line, best["tp"], best["tp_rate"], best["fp"],
                   best["fp_rate"], best["fn"], best["missing_ratio"]))
    best_out.close()
    indexs["count"] += 1
    return (best_para, best)

def scoring_function(best, stat_value, indexs):
    indexs["change"] = False
    if (stat_value["tp_rate"] == best["tp_rate"]) and (
        stat_value["fp_rate"] == best["fp_rate"]):
        pass
    else:
        if (stat_value["tp_rate"] - best["tp_rate"]) >= 0.1:
            indexs["change"] = True
        elif (best["tp_rate"] - stat_value["tp_rate"]) <= 0.1:
            if (best["tp_rate"] <= stat_value["tp_rate"]) and (
                best["fp_rate"] >= stat_value["fp_rate"]):
                indexs["change"] = True
            elif (stat_value["tp_rate"] - best["tp_rate"] >= 0.01) and (
                  stat_value["fp_rate"] - best["fp_rate"] <= 0.00005):
                indexs["change"] = True
            elif (best["tp_rate"] - stat_value["tp_rate"] <= 0.01) and (
                  best["fp_rate"] - stat_value["fp_rate"] >= 0.00005):
                indexs["change"] = True
            else:
                tp_diff = float(best["tp"] - stat_value["tp"])
                if tp_diff > 0:
                    if float(best["fp"] - stat_value["fp"]) >= 5 * tp_diff:
                        indexs["change"] = True
                elif tp_diff < 0:
                    tp_diff = tp_diff * -1
                    if float(stat_value["fp"] - best["fp"]) <= 5 * tp_diff:
                        indexs["change"] = True

def check_overlap(overlap, pre_tss, nums, length, manual, predict, pre_pos):
    if overlap:
        if pre_tss:
            pre_tss.attributes["print"] = True
            tss = pre_tss
        else:
            tss = predict
        if (tss.start <= int(length)):
            if (pre_pos != -1):
                if (tss.start - pre_pos != 0):
                    nums["overlap"] += 1
                else:
                    nums["overlap"] += 1
            else:
                nums["overlap"] += 1
        overlap = False
        pre_pos = tss.start
    else:
        if (manual.start <= int(length)):
            nums["manual"] += 1
    return (overlap, pre_pos)

def comparison(manuals, predicts, nums, length, cluster):
    overlap = False
    pre_pos = -1
    for tss_m in manuals:
        pre_tss = None
        for tss_p in predicts:
            if (tss_p.strand == tss_m.strand) and (
                (tss_p.seq_id == tss_m.seq_id) or (
                (tss_p.seq_id == tss_m.seq_id[:-2]) and (
                 tss_m.seq_id[-2] == "."))):
                if (tss_p.start == tss_m.start):
                    tss_p.attributes["print"] = True
                    overlap = True
                    pre_tss = None
                    break
                elif (math.fabs(tss_p.start - tss_m.start) <= cluster):
                    overlap = True
                    pre_tss = tss_p
        datas = check_overlap(overlap, pre_tss, nums, length,
                              tss_m, tss_p, pre_pos)
        overlap = datas[0]
        pre_pos = datas[1]
    for tss_p in predicts:
        if tss_p.attributes["print"] is False:
            if (tss_p.start <= int(length)):
                nums["predict"] += 1

def read_predict_manual_gff(gff_file, gene_length, gffs, cluster):
    num = 0
    f_h = open(gff_file, "r")
    for entry in Gff3Parser().entries(f_h):
        if (entry.start <= int(gene_length)):
            num += 1
            entry.attributes["print"] = False
            gffs.append(entry)
    return num

def compare_manual_predict(total_step, para_list, gff_files, out_path,
                           out, manual, cores, gene_length, cluster):
    manuals = []
    manual_fh = open(manual, "r")
    stats = []
    gff_parser = Gff3Parser()
    count = 0
    total_step = total_step - int(cores) + 1
    num_manual = read_predict_manual_gff(manual, gene_length, manuals, cluster)
    for gff_file in gff_files:
        nums = {"overlap": 0, "predict": 0, "manual": 0}
        predicts = []
        para = "_".join(["he", str(para_list[count]["height"]),
                         "rh", str(para_list[count]["re_height"]),
                         "fa", str(para_list[count]["factor"]),
                         "rf", str(para_list[count]["re_factor"]),
                         "bh", str(para_list[count]["base_height"]),
                         "ef", str(para_list[count]["enrichment"]),
                         "pf", str(para_list[count]["processing"])])
        num_predict = read_predict_manual_gff(gff_file,
                      gene_length, predicts, cluster)
        comparison(manuals, predicts, nums, gene_length, cluster)
        out.write("{0}\t{1}\tTP\t{2}\tTP_rate\t{3}\t".format(
                  total_step, para, nums["overlap"],
                  float(nums["overlap"]) / float(num_manual)))
        out.write("FP\t{0}\tFP_rate\t{1}\tFN\t{2}\tmissing_ratio\t{3}\n".format(
                  nums["predict"],
                  float(nums["predict"]) / float(int(gene_length) - num_manual),
                  nums["manual"],
                  float(nums["manual"]) / float(num_manual)))
        stats.append({"tp": nums["overlap"],
                      "tp_rate": float(nums["overlap"]) / float(num_manual),
                      "fp": nums["predict"],
                      "fp_rate": float(nums["predict"]) / float(gene_length - num_manual),
                      "fn": nums["manual"],
                      "missing_ratio": float(nums["manual"]) / float(num_manual)})
        total_step += 1
        count += 1
    manual_fh.close()
    return stats

def convert2gff(cores, out_path, strain, gff_files, tss_pro):
    for core in range(1, cores+1):
        output_folder = os.path.join(out_path,
                        "_".join(["MasterTable", str(core)]))
        gff_file = os.path.join(output_folder,
                   "_".join(["TSSpredator", str(core) + ".gff"]))
        Converter().convert_mastertable2gff(
                    os.path.join(output_folder, "MasterTable.tsv"),
                    "TSSpredator", tss_pro, strain, gff_file)
        gff_files.append(gff_file)

def run_TSSpredator(tsspredator_path, config_file):
    folders = config_file.split("/")
    out_path = "/".join(folders[:-1])
    out = open(os.path.join(out_path, "log.txt"), "w")
    p = Popen(["java", "-Xmx2G", "-jar",
           tsspredator_path, config_file], stdout=out)
    return p

def run_TSSpredator_paralle(config_files, tsspredator_path, processes):
    for config_file in config_files:
        process = run_TSSpredator(tsspredator_path, config_file)
        processes.append(process)
    for p in processes:
        p.wait()
        if p.stdout:
            p.stdout.close()
        if p.stdin:
            p.stdin.close()
        if p.stderr:
            p.stderr.close()
        try:
            p.kill()
        except OSError:
            pass
    time.sleep(5)

def print_lib(lib_num, lib_list, out, wig_folder, prefix):
    for num_id in range(1, lib_num+1):
        cond_list = []
        for lib in lib_list:
            if num_id == lib["condition"]:
                cond_list.append(lib)
        cond_sort_list = sorted(cond_list, key=lambda k: k['replicate'])
        for cond in cond_sort_list:
            out.write("{0}_{1}{2} = {3}/{4}\n".format(
                      prefix, cond["condition"], cond["replicate"],
                      wig_folder, cond["wig"]))

def assign_dict(lib_datas):
    return {"wig": lib_datas[0],
            "tex": lib_datas[1],
            "condition": int(lib_datas[2]),
            "replicate": lib_datas[3],
            "strand": lib_datas[4]}

def import_lib(libs, wig_folder, project_strain_name, rep_set, lib_dict,
               out, gff, program, list_num_id, fasta):
    lib_num = 0
    for lib in libs:
        lib_datas = lib.split(":")
        if lib_datas[0].endswith(".wig") is not True:
            print("Error:Exist a not proper wig files!!")
            sys.exit()
        for wig in os.listdir(wig_folder):
            filename = wig.split("_STRAIN_")
            if (filename[0] == lib_datas[0][:-4]) and \
               (filename[1][:-4] == project_strain_name):
                lib_datas[0] = wig
            elif (filename[0] == lib_datas[0][:-4]) and \
                 ("." == filename[1][-6]) and \
                 (filename[1][:-6] == project_strain_name):
                lib_datas[0] = wig
        if int(lib_datas[2]) > lib_num:
            lib_num = int(lib_datas[2])
        if lib_datas[3] not in rep_set:
            rep_set.add(lib_datas[3])
        if (lib_datas[1] == "tex") and \
           (lib_datas[4] == "+"):
            lib_dict["fp"].append(assign_dict(lib_datas))
        elif (lib_datas[1] == "tex") and \
             (lib_datas[4] == "-"):
            lib_dict["fm"].append(assign_dict(lib_datas))
        elif (lib_datas[1] == "notex") and \
             (lib_datas[4] == "+"):
            lib_dict["np"].append(assign_dict(lib_datas))
        elif (lib_datas[1] == "notex") and \
             (lib_datas[4] == "-"):
            lib_dict["nm"].append(assign_dict(lib_datas))
    for num_id in range(1, lib_num+1):
        out.write("annotation_{0} = {1}\n".format(num_id, gff))
    if program.lower() == "tss":
        print_lib(lib_num, lib_dict["fm"], out, wig_folder, "fivePrimeMinus")
        print_lib(lib_num, lib_dict["fp"], out, wig_folder, "fivePrimePlus")
    elif program.lower() == "processing_site":
        print_lib(lib_num, lib_dict["nm"], out, wig_folder, "fivePrimeMinus")
        print_lib(lib_num, lib_dict["np"], out, wig_folder, "fivePrimePlus")
    else:
        print("Error:the program name is wrong!!")
        sys.exit()
    for num_id in range(1, lib_num+1):
        out.write("genome_%s = %s\n" % (str(num_id), fasta))
    for num_id in range(1, lib_num+1):
        list_num_id.append(str(num_id))
    return lib_num

def gen_config(para_list, out_path, core, libs, wig, project_strain, utr,
               fasta, output_prefix, gff, program, cluster, replicate_num,):
    files = os.listdir(out_path)
    if "MasterTable_" + str(core) not in files:
        os.mkdir(os.path.join(out_path, "MasterTable_" + str(core)))
    lib_dict = {"fp": [], "fm": [], "nm": [], "np": []}
    rep_set = set()
    list_num_id = []
    filename = os.path.join(out_path, "config_" + str(core) + ".ini")
    out = open(filename, "w")
    out.write("TSSinClusterSelectionMethod = HIGHEST\n")
    out.write("allowedCompareShift = 1\n")
    out.write("allowedRepCompareShift = 1\n")
    lib_num = import_lib(libs, wig, project_strain, rep_set, lib_dict,
               out, gff, program, list_num_id, fasta)
    out.write("idList = ")
    out.write(",".join(list_num_id) + "\n")
    out.write("maxASutrLength = 100\n")
    out.write("maxGapLengthInGene = 500\n")
    out.write("maxNormalTo5primeFactor = {0}\n".format(para_list["processing"]))
    out.write("maxTSSinClusterDistance = {0}\n".format(cluster + 1))
    out.write("maxUTRlength = {0}\n".format(utr))
    out.write("min5primeToNormalFactor = {0}\n".format(para_list["enrichment"]))
    out.write("minCliffFactor = {0}\n".format(para_list["factor"]))
    out.write("minCliffFactorDiscount = {0}\n".format(para_list["re_factor"]))
    out.write("minCliffHeight = {0}\n".format(para_list["height"]))
    out.write("minCliffHeightDiscount = {0}\n".format(para_list["re_height"]))
    out.write("minNormalHeight = {0}\n".format(para_list["base_height"]))
    out.write("minNumRepMatches = {0}\n".format(replicate_num))
    out.write("minPlateauLength = 0\n")
    out.write("mode = cond\n")
    out.write("normPercentile = 0.9\n")
    if (program.lower() == "tss"):
        print_lib(lib_num, lib_dict["nm"], out, wig, "normalMinus")
        print_lib(lib_num, lib_dict["np"], out, wig, "normalPlus")
    elif (program.lower() == "processing_site"):
        print_lib(lib_num, lib_dict["fm"], out, wig, "normalMinus")
        print_lib(lib_num, lib_dict["fp"], out, wig, "normalPlus")
    out.write("numReplicates = {0}\n".format(len(rep_set)))
    out.write("numberOfDatasets = {0}\n".format(lib_num))
    out.write("outputDirectory = {0}\n".format(
              os.path.join(out_path, "_".join(["MasterTable", str(core)]))))
    for prefix_id in range(len(output_prefix)):
        out.write("outputPrefix_{0} = {1}\n".format(
                  prefix_id + 1, output_prefix[prefix_id]))
    out.write("projectName = {0}\n".format(project_strain))
    out.write("superGraphCompatibility = igb\n")
    out.write("texNormPercentile = 0.5\n")
    out.write("writeGraphs = 0\n")
    out.write("writeNocornacFiles = 0\n")
    out.close()
    return filename

def run_tss_and_stat(indexs, steps, cores, list_num, seeds, diff_h, diff_f,
                     out_path, tsspredator_path, stat_out, best_para,
                     current_para, gene_length, wig, project_strain, fasta,
                     output_prefix, gff, program, libs, manual, best, cluster,
                     utr, replicate_num):
    if indexs["step"] >= steps + int(cores):
        return (True, best_para)
    elif len(list_num) == indexs["length"]:
        indexs["step"] = indexs["step"] - 1
        seeds["pre_seed"].append(seeds["seed"])
    elif (diff_h <= 0) or \
         (diff_f <= 0):
        indexs["step"] = indexs["step"] - 1
        list_num = list_num[:-1]
    else:
        indexs["num"] += 1
        seeds["pre_seed"] = []
        if indexs["num"] == cores:
            index = 0
            config_files = []
            gff_files = []
            for para in list_num[-1 * cores:]:
                index += 1
                print("_".join([str(para["height"]), str(para["re_height"]),
                      str(para["factor"]), str(para["re_factor"]),
                      str(para["base_height"]), str(para["enrichment"]),
                      str(para["processing"])]))
                config_files.append(gen_config(para, out_path, index, libs, wig,
                                    project_strain, utr, fasta, output_prefix,
                                    gff, program, cluster, replicate_num))
            indexs["count"] = 0
            processes = []
            run_TSSpredator_paralle(config_files, tsspredator_path, processes)###
            convert2gff(cores, out_path, project_strain, gff_files, program)
            stat_values = compare_manual_predict(
                              indexs["step"], list_num[-1 * cores:],
                              gff_files, out_path, stat_out, manual,
                              cores, gene_length, cluster)
            for stat_value in stat_values:
                if indexs["first"]:
                    indexs["first"] = False
                    best = stat_value
                    best_para = copy.deepcopy(list_num[-1 * cores + indexs["count"]])
                else:
                    scoring_function(best, stat_value, indexs)
                datas = compute_stat(stat_value, best, best_para,
                                     cores, list_num, out_path, indexs)
                best_para = datas[0]
                best = datas[1]
            indexs["switch"] += 1
            stat_values = []
            indexs["num"] = 0
        current_para = copy.deepcopy(best_para)
    return (False, best_para, best)

def minus_process(num_type, new_para, max_num, best_num,
                  actions, list_num, compare):
    if num_type == "base_height":
        new_para[num_type] = new_para[num_type] - 0.001
        new_para[num_type] = float('%.3f' % new_para[num_type])
        while True:
            if new_para[num_type] < 0.0:
                new_para[num_type] = 0.0
                if new_para in list_num:
                    new_para[num_type] = best_num
                    actions["minus"] = True
                    actions["in_or_de"] = 2
                    break
                else:
                    list_num.append(new_para)
                    return new_para[num_type]
            elif (new_para in list_num):
                new_para[num_type] = new_para[num_type] - 0.001
                new_para[num_type] = float('%.3f' % new_para[num_type])
                continue
            else:
                list_num.append(copy.deepcopy(new_para))
                return new_para[num_type]
    else:
        new_para[num_type] = new_para[num_type] - 0.1
        new_para[num_type] = float('%.1f' % new_para[num_type])
        while True:
            if new_para[num_type] <= 0.0:
                new_para[num_type] = best_num
                actions["in_or_de"] = 2
                actions["minus"] = True
                break
            elif (new_para in list_num):
                new_para[num_type] = new_para[num_type] - 0.1
                new_para[num_type] = float('%.1f' % new_para[num_type])
                continue
            elif ((num_type == "factor") or \
                 (num_type == "height")) and \
                 (new_para[num_type] <= compare):
                new_para[num_type] = best_num
                actions["in_or_de"] = 2
                actions["minus"] = True
                break
            else:
                list_num.append(copy.deepcopy(new_para))
                return new_para[num_type]
    return None

def plus_process(num_type, new_para, max_num,
                 best_num, actions, list_num, compare):
    if num_type == "base_height":
        new_para[num_type] = new_para[num_type] + 0.001
        new_para[num_type] = float('%.3f' % new_para[num_type])
        while True:
            if new_para[num_type] >= max_num:
                new_para[num_type] = best_num
                actions["in_or_de"] = 1
                actions["plus"] = True
                break
            elif (new_para in list_num):
                new_para[num_type] = new_para[num_type] + 0.001
                new_para[num_type] = float('%.3f' % new_para[num_type])
                continue
            else:
                list_num.append(copy.deepcopy(new_para))
                return new_para[num_type]
    else:
        new_para[num_type] = new_para[num_type] + 0.1
        new_para[num_type] = float('%.1f' % new_para[num_type])
        while True:
            if new_para[num_type] >= max_num:
                new_para[num_type] = best_num
                actions["in_or_de"] = 1
                actions["plus"] = True
                break
            elif (new_para in list_num):
                new_para[num_type] = new_para[num_type] + 0.1
                new_para[num_type] = float('%.1f' % new_para[num_type])
                continue
            elif ((num_type == "re_factor") or \
                 (num_type == "re_height")) and \
                 (new_para[num_type] >= compare):
                new_para[num_type] = best_num
                actions["in_or_de"] = 1
                actions["plus"] = True
                break
            else:
                list_num.append(copy.deepcopy(new_para))
                return new_para[num_type]
    return None

def small_change(max_num, num_type, compare, list_num, best_num, best_para):
    new_para = copy.deepcopy(best_para)
    actions = {"plus": False, "minus": False}
    step = 0
    if new_para[num_type] >= max_num:
        actions["in_or_de"] = 1
    elif new_para[num_type] <= 0:
        actions["in_or_de"] = 2
    else:
        actions["in_or_de"] = random.randint(0, 9)
    while True:
        step += 1
        if step >= 1000:
            return best_num
        if (actions["plus"] is True) and (actions["minus"] is True):
            new_para[num_type] = best_num
            return new_para[num_type]
        if actions["in_or_de"] % 2 == 0:
            tmp_para = plus_process(num_type, new_para, max_num, best_num,
                                    actions, list_num, compare)
        if actions["in_or_de"] % 2 == 1:
            tmp_para = minus_process(num_type, new_para, max_num, best_num,
                                     actions, list_num, compare)
        if tmp_para is not None:
            return tmp_para

def run_small_change_part(seeds, features, indexs, current_para,
                          best_para, list_num, max_num):
    while True:
        seeds["seed"] = random.randint(0, 6)
        if seeds["seed"] in seeds["pre_seed"]:
            if len(seeds["pre_seed"]) == 7:
                indexs["switch"] += 1
                features["pre_feature"] = features["feature"]
                break
            else:
                continue
        else:
            break
    if seeds["seed"] == 0:
        current_para["height"] = small_change(max_num["height"], "height",
                                 best_para["re_height"], list_num,
                                 best_para["height"], best_para)
    elif seeds["seed"] == 1:
        current_para["re_height"] = small_change(max_num["re_height"],
                                    "re_height", best_para["height"],
                                    list_num, best_para["re_height"], best_para)
    elif seeds["seed"] == 2:
        current_para["factor"] = small_change(max_num["factor"], "factor",
                                 best_para["re_factor"], list_num,
                                 best_para["factor"], best_para)
    elif seeds["seed"] == 3:
        current_para["re_factor"] = small_change(max_num["re_factor"],
                                    "re_factor", best_para["factor"], list_num,
                                    best_para["re_factor"], best_para)
    elif seeds["seed"] == 4:
        current_para["base_height"] = small_change(max_num["base_height"],
                                      "base_height", best_para["base_height"],
                                      list_num, best_para["base_height"],
                                      best_para)
    elif seeds["seed"] == 5:
        current_para["enrichment"] = small_change(max_num["enrichment"],
                                     "enrichment", best_para["enrichment"],
                                     list_num, best_para["enrichment"],
                                     best_para)
    elif seeds["seed"] == 6:
        current_para["processing"] = small_change(max_num["processing"],
                                     "processing", best_para["processing"],
                                     list_num, best_para["processing"],
                                     best_para)
    return current_para

def check_fit(num_type, index_large, seed, number, number_par, compare):
    if (((num_type == "height") and (index_large[seed] == "re_height")) or (
        (num_type == "factor") and (index_large[seed] == "re_factor"))) and (
        number <= number_par):
        print("A")
        return False
    elif (((num_type == "re_factor") and (index_large[seed] == "factor")) or (
          (num_type == "re_height") and (index_large[seed] == "height"))) and (
          number >= number_par):
        print("B")
        return False
    elif ((num_type == "factor") or (num_type == "height")) and (
          number <= float(compare)):
        print("C")
        return False
    elif ((num_type == "re_factor") or (num_type == "re_height")) and (
          number >= float(compare)):
        print("D")
        return False
    elif ((index_large[seed] == "factor") or (
          index_large[seed] == "height")) and (
          number_par <= float(compare)):
        print("E")
        return False
    elif ((index_large[seed] == "re_factor") or (
          index_large[seed] == "re_height")) and (
          number_par >= float(compare)):
        print("F")
        return False
    return True

def gen_large_random(max_num, num_type, compare, list_num, origin_num,
                     best_para, index_large, indexs):
    new_para = copy.deepcopy(best_para)
    step = 0
    while True:
        step += 1
        if step >= 1000000:
            return best_para
        seed = random.randint(0, 6)
        if num_type == index_large[seed]:
            continue
        if num_type == "base_height":
            number = round(random.uniform(0.001, max_num[num_type]), 3)
            number = '%.3f' % number
            number = float(number)
        else:
            number = round(random.uniform(0.1, max_num[num_type]), 1)
            number = '%.1f' % number
            number = float(number)
        if index_large[seed] == "base_height":
            number_par = round(random.uniform(0.001,
                               max_num[index_large[seed]]), 3)
            number_par = '%.3f' % number_par
            number_par = float(number_par)
        else:
            number_par = round(random.uniform(0.1,
                               max_num[index_large[seed]]), 1)
            number_par = '%.1f' % number_par
            number_par = float(number_par)
#        fit = check_fit(num_type, index_large, seed,
#                        number, number_par, compare)
#        if fit is False:
#            continue
#        else:
        new_para[num_type] = number
        new_para[index_large[seed]] = number_par
        if new_para in list_num:
            continue
        else:
            if (new_para["height"] <= new_para["re_height"]) or \
               (new_para["factor"] <= new_para["re_factor"]):
                continue
            else:
                list_num.append(copy.deepcopy(new_para))
                return new_para

def run_large_change_part(seeds, features, indexs, current_para, max_num,
                          best_para, list_num):
    index_large = {0: "height", 1: "re_height", 2: "factor", 3: "re_factor",
                   4:"base_height", 5: "enrichment", 6: "processing"}
    while True:
        seeds["seed"] = random.randint(0, 6)
        if seeds["seed"] in seeds["pre_seed"]:
            if len(seeds["pre_seed"]) == 7:
                features["pre_feature"] = features["feature"]
                indexs["switch"] += 1
                break
            else:
                continue
        else:
            break
    if seeds["seed"] == 0:
        current_para = gen_large_random(max_num, "height",
                       best_para["re_height"], list_num, best_para["height"],
                       best_para, index_large, indexs)
    elif seeds["seed"] == 1:
        current_para = gen_large_random(max_num, "re_height",
                       best_para["height"], list_num, best_para["re_height"],
                       best_para, index_large, indexs)
    elif seeds["seed"] == 2:
        current_para = gen_large_random(max_num, "factor",
                       best_para["re_factor"], list_num, best_para["factor"],
                       best_para, index_large, indexs)
    elif seeds["seed"] == 3:
        current_para = gen_large_random(max_num, "re_factor",
                       best_para["factor"], list_num, best_para["re_factor"],
                       best_para, index_large, indexs)
    elif seeds["seed"] == 4:
        current_para = gen_large_random(max_num, "base_height",
                       best_para["base_height"], list_num,
                       best_para["base_height"], best_para, index_large, indexs)
    elif seeds["seed"] == 5:
        current_para = gen_large_random(max_num, "enrichment",
                       best_para["enrichment"], list_num, best_para["enrichment"],
                       best_para, index_large, indexs)
    elif seeds["seed"] == 6:
        current_para = gen_large_random(max_num, "processing",
                       best_para["processing"], list_num, best_para["processing"],
                       best_para, index_large, indexs)
    return current_para

def run_random_part(current_para, list_num, max_num, steps, indexs):
    tmp_random_step = 0
    while True:
        current_para["height"] = round(random.uniform(
                                       0.1, max_num["height"]), 1)
        current_para["re_height"] = round(random.uniform(
                                          0.1, max_num["re_height"]), 1)
        current_para["factor"] = round(random.uniform(
                                       0.1, max_num["factor"]), 1)
        current_para["re_factor"] = round(random.uniform(
                                          0.1, max_num["re_factor"]), 1)
        current_para["enrichment"] = round(random.uniform(
                                           0.1, max_num["enrichment"]), 1)
        current_para["processing"] = round(random.uniform(
                                           0.1, max_num["processing"]), 1)
        current_para["base_height"] = round(random.uniform(
                                            0.001, max_num["base_height"]), 3)
        if (current_para["height"] > current_para["re_height"]) and (
            current_para["factor"] > current_para["re_factor"]) and (
            current_para not in list_num):
            list_num.append(copy.deepcopy(current_para))
            break
        tmp_random_step += 1
        if tmp_random_step >= steps:
            indexs["switch"] += 1
            return None
    return current_para

def optimization_process(indexs, current_para, list_num, max_num, best_para,
        steps, cores, out_path, tsspredator_path, stat_out, best, libs, wig,
        project_strain, fasta, output_prefix, gff, program, gene_length,
        manual, cluster, utr, replicate_num):
    features = {"pre_feature": "", "feature": ""}
    seeds = {"pre_seed": [], "seed" : 0}
    tests = {"test1": [], "test2": ""}
    tmp_step = 0
    while True:
        if indexs["exist"] is False:
            indexs["exist"] = True
            features["feature"] = ""
        elif (indexs["switch"] % 3 == 0):
            features["feature"] = "r"
            if features["feature"] != features["pre_feature"]:
                seeds["pre_seed "] = []
            current_para = run_random_part(current_para, list_num,
                                           max_num, steps, indexs)
            if current_para is None:
                tmp_step += 1
        elif (indexs["switch"] % 3 == 1):
            features["feature"] = "l"
            if features["feature"] != features["pre_feature"]:
                seeds["pre_seed"] = []
            current_para = run_large_change_part(seeds, features, indexs,
                               current_para, max_num, best_para, list_num)
        else:
            features["feature"] = "s"
            if features["feature"] != features["pre_feature"]:
                seeds["pre_seed"] = []
            current_para = run_small_change_part(seeds, features, indexs,
                               current_para, best_para, list_num, max_num)
        diff_h = float(current_para["height"]) - float(current_para["re_height"])
        diff_f = float(current_para["factor"]) - float(current_para["re_factor"])
        if current_para is not None:
            datas = run_tss_and_stat(indexs, steps, cores, list_num, seeds,
                        diff_h, diff_f, out_path, tsspredator_path, stat_out,
                        best_para, current_para, gene_length, wig,
                        project_strain, fasta, output_prefix, gff, program,
                        libs, manual, best, cluster, utr, replicate_num)
            tmp_step = 0
        if tmp_step >= 2:
            print("The number of steps may be enough..., it may not be able to find more parameters...\n")
            sys.exit()
        best_para = datas[1]
        if datas[0]:
            break
        else:
            best = datas[2]
        indexs["length"] = len(list_num)
        features["pre_feature"] = features["feature"]
        indexs["step"] += 1
        if indexs["step"] >= steps:
            break

def start_data(out_path, current_para, best_para, indexs, max_num):
    indexs["step"] = 0
    list_num = []
    while True:
        current_para["height"] = round(random.uniform(
                                       0.1, max_num["height"]), 1)
        current_para["re_height"] = round(random.uniform(
                                          0.1, max_num["re_height"]), 1)
        current_para["factor"] = round(random.uniform(
                                       0.1, max_num["factor"]), 1)
        current_para["re_factor"] = round(random.uniform(
                                          0.1, max_num["re_factor"]), 1)
        current_para["enrichment"] = round(random.uniform(
                                           0.1, max_num["enrichment"]), 1)
        current_para["processing"] = round(random.uniform(
                                           0.1, max_num["processing"]), 1)
        current_para["base_height"] = round(random.uniform(
                                            0.001, max_num["base_height"]), 3)
        best_para = copy.deepcopy(current_para)
        if (current_para["height"] > current_para["re_height"]) and \
           (current_para["factor"] > current_para["re_factor"]):
            break
    list_num = [{"height": current_para["height"],
                 "re_height": current_para["re_height"],
                 "factor": current_para["factor"],
                 "re_factor": current_para["re_factor"],
                 "base_height": current_para["base_height"],
                  "enrichment": current_para["enrichment"],
                 "processing": current_para["processing"]}]
    return list_num

def extend_data(out_path, best, best_para, step):
    print("extend step from {0}".format(step))
    print("\t".join(["Best Parameter:height={0}", "reduction_height={1}",
                     "factor={2}", "reduction_factor={3}", "base_height={4}",
                     "enrichment_factor={5}", "processing_factor={6}"]).format(
          best_para["height"],
          best_para["re_height"],
          best_para["factor"],
          best_para["re_factor"],
          best_para["base_height"],
          best_para["enrichment"],
          best_para["processing"]))
    print("Best:TP={0}\tTP_rate={1}\tFP={2}\tFP_rate={3}\tFN={4}\tMissing_ratio={5}".format(
          best["tp"], best["tp_rate"], best["fp"], best["fp_rate"],
          best["fn"], best["missing_ratio"]))
    current_para = copy.deepcopy(best_para)
    return current_para

def load_stat_csv(out_path, list_num, best, best_para, indexs):
    f_h = open(os.path.join(out_path, "stat.csv"), "r")
    first_line = True
    line_num = 0
    for row in csv.reader(f_h, delimiter="\t"):
        line_num += 1
        paras = row[1].split("_")
        if len(row) == 14:
            prev_stat = {"tp": int(row[3]), "tp_rate": float(row[5]),
                         "fp": int(row[7]), "fp_rate": float(row[9])}
            list_num.append({"height": float(paras[1]),
                             "re_height": float(paras[3]),
                             "factor": float(paras[5]),
                             "re_factor": float(paras[7]),
                             "base_height": float(paras[9]),
                             "enrichment": float(paras[11]),
                             "processing": float(paras[13])})
            if first_line:
                first_line = False
                indexs["change"] = True
            else:
                scoring_function(best, prev_stat, indexs)
            if indexs["change"]:
                best_para = {"height": float(paras[1]),
                             "re_height": float(paras[3]),
                             "factor": float(paras[5]),
                             "re_factor": float(paras[7]),
                             "base_height": float(paras[9]),
                             "enrichment": float(paras[11]),
                             "processing": float(paras[13])}
                best["tp"] = float(row[3])
                best["tp_rate"] = float(row[5])
                best["fp"] = float(row[7])
                best["fp_rate"] = float(row[9])
                best["fn"] = float(row[11])
                best["missing_ratio"] = float(row[13])
            indexs["step"] = int(row[0]) + 1
    return (line_num, best, best_para)

def reload_data(out_path, list_num, best, best_para, indexs):
    first_line = True
    change_para = False
    indexs["switch"] = 1
    indexs["exist"] = True
    datas = load_stat_csv(out_path, list_num, best, best_para, indexs)
    line_num = datas[0]
    best = datas[1]
    best_para = datas[2]
    if len(list_num) > 0:
        indexs["extend"] = True
    else:
        print("Error: the stat.csv has something wrong, please check it!!!")
        sys.exit()
    new_line = 0
    new_stat = open("tmp.csv", "w")
    with open(os.path.join(out_path, "stat.csv"), "r") as fh:
        for line in fh:
            new_line += 1
            line = line.strip()
            if new_line >= line_num:
                break
            else:
                new_stat.write(line + "\n")
    os.rename("tmp.csv", os.path.join(out_path, "stat.csv"))
    return (best_para, best)

def get_gene_length(fasta, strain):
    seq = ""
    detect = False
    with open(fasta, "r") as f_h:
        for line in f_h:
            line = line.strip()
            if line.startswith(">"):
                detect = False
                if line[1:] == strain:
                    detect = True
            else:
                if detect:
                    seq = seq + line
    return len(seq)

def initiate(height, reduction_height, factor, reduction_factor, base_height,
             enrichment, processing):
    max_num = {"height": height, "re_height": reduction_height,
               "factor": factor, "re_factor": reduction_factor,
               "base_height": base_height, "enrichment": enrichment,
               "processing": processing}
    best_para = {"height": 0, "re_height": 0,
                 "factor": 0, "re_factor": 0,
                 "base_height": 0, "enrichment": 0,
                 "processing": 0}
    current_para = {"height": 0, "re_height": 0,
                    "factor": 0, "re_factor": 0,
                    "base_height": 0, "enrichment": 0,
                    "processing": 0}
    indexs = {"switch": 0, "extend": False, "exist": False, "step": 0,
              "first": True, "num": 0, "length": 0, "change": False,
              "count": 0}
    return max_num, best_para, current_para, indexs

def optimization(tsspredator_path, height, reduction_height, factor, 
                 reduction_factor, base_height, enrichment, processing, 
                 output_folder, cores, wig, project_strain, fasta,
                 output_prefix, steps, gff, program, manual, libs,
                 gene_length, cluster, utr, replicate_num):
    pre_seed = []
    best = {}
    pre_feature = ""
    max_num, best_para, current_para, indexs = initiate(height,
                                               reduction_height, factor,
                                               reduction_factor, 
                                               base_height, enrichment,
                                               processing)
    out_path = os.path.join(output_folder, "optimized_TSSpredator")
    files = os.listdir(output_folder)
    stat_file = os.path.join(out_path, "stat.csv")
    if gene_length is None:
        gene_length = get_gene_length(fasta, project_strain)
    else:
        gene_length = int(gene_length)
    if "optimized_TSSpredator" not in files:
        os.mkdir(out_path)
        list_num = []
        stat_out = open(stat_file, "w")
    else:
        if ("stat.csv" in os.listdir(out_path)):
            list_num = []
            datas = reload_data(out_path, list_num, best, best_para, indexs)
            best_para = datas[0]
            best = datas[1]
            current_para = extend_data(out_path, best, best_para, indexs["step"])
            stat_out = open(stat_file, "a")
            indexs["first"] = False
        else:
            list_num = []
            stat_out = open(stat_file, "w")
    optimization_process(indexs, current_para, list_num, max_num, best_para,
                         steps, cores, out_path, tsspredator_path, stat_out,
                         best, libs, wig, project_strain, fasta, output_prefix,
                         gff, program, gene_length, manual, cluster, utr,
                         replicate_num)
