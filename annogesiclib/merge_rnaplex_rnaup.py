from annogesiclib.gff3 import Gff3Parser


def detect_energy(line, srna):
    duplex = line.split(" ")[0]
    if ("(" in duplex) or (")" in duplex):
        energy = line.split(":")[1].split("(")[1].split("=")[0]
        energy = float(energy.strip())
        if energy < srna["energy"]:
            srna["energy"] = energy


def mod_srna_tar_pos(gff, pos, type_, pre_target, suf_target):
    start = int(pos.split(",")[0])
    end = int(pos.split(",")[-1])
    if (gff.strand == "+"):
        if type_ == "srna":
            g_start = gff.start + start - 1
            g_end = gff.start + end - 1
        else:
            g_start = gff.start - pre_target + start - 1
            g_end = gff.start - pre_target + end - 1
    else:
        if type_ == "srna":
            g_end = gff.end - start + 1
            g_start = gff.end - end + 1
        else:
            g_start = gff.start + suf_target - start + 1
            g_end = gff.start + suf_target - end + 1
    return g_start, g_end

def print_rank_one(srnas, out, feature, gffs, srna_gffs, args_tar):
    out.write("\t".join(["sRNA", "strain", "sRNA_position",
                         "sRNA_interacted_position_" + feature,
                         "sRNA_strand",
                         "target", "target_position",
                         "target_interacted_position_" + feature,
                         "target_strand", "energy_" + feature,
                         "rank_" + feature]) + "\n")
    for method, srna_datas in srnas.items():
        for srna_id, targets in srna_datas.items():
            rank = 0
            sort_targets = sorted(targets, key=lambda k: (k["energy"]))
            for target in sort_targets:
                if target["energy"] < 0:
                    rank += 1
                    target["rank"] = rank
                    if rank <= args_tar.top:
                        if method == feature:
                            srna_infos = get_srna_name(srna_gffs, srna_id)
                            name = srna_infos[0]
                            srna_info = srna_infos[1]
                            target_info = get_target_info(
                                gffs, target["target"])
                            s_start, s_end = mod_srna_tar_pos(
                                    srna_info, target["srna_pos"], "srna",
                                    args_tar.tar_start, args_tar.tar_end)
                            t_start, t_end = mod_srna_tar_pos(
                                    target_info, target["tar_pos"], "tar",
                                    args_tar.tar_start, args_tar.tar_end)
                            if srna_info is not None:
                                out.write("\t".join([
                                    name, str(srna_info.seq_id),
                                    "-".join([str(srna_info.start),
                                              str(srna_info.end)]),
                                    "-".join([str(s_start), str(s_end)]),
                                    srna_info.strand, target["target"],
                                    "-".join([str(target_info.start),
                                              str(target_info.end)]),
                                    "-".join([str(t_start), str(t_end)]),
                                    target_info.strand, str(target["energy"]),
                                    str(target["rank"])]) + "\n")

def extract_pos(line, srnas, method):
    if (line.startswith("(")) or (
            line.startswith(")")) or (
            line.startswith(".")):
        tar_pos = line.split(":")[0].strip().split(" ")[-1]
        srna_pos = line.split(":")[-1].strip().split(" ")[0]
        srnas["tar_pos"] = tar_pos
        srnas["srna_pos"] = srna_pos

def read_table(gffs, rnaplex, rnaup):
    srnas = {"RNAup": {}, "RNAplex": {}}
    start = False
    count_seq = 0
    srna_names = set()
    for gff in gffs:
        if gff.attributes["ID"] not in srna_names:
            srna_names.add(gff.attributes["ID"])
    if rnaplex is not None:
        with open(rnaplex, "r") as p_h:
            for line in p_h:
                line = line.strip()
                if line.startswith(">"):
                    start = True
                    count_seq += 1
                    if count_seq == 1:
                        target = line[1:]
                    elif count_seq == 2:
                        srna = line[1:]
                        if srna not in srnas["RNAplex"].keys():
                            srnas["RNAplex"][srna] = []
                        srnas["RNAplex"][srna].append({"target": target,
                                                       "energy": 0})
                        count_seq = 0
                else:
                    if start:
                        detect_energy(line, srnas["RNAplex"][srna][-1])
                        extract_pos(line, srnas["RNAplex"][srna][-1], "RNAplex")
    if rnaup is not None:
        with open(rnaup, "r") as u_h:
            for line in u_h:
                line = line.strip()
                if line.startswith(">"):
                    if line[1:] in srna_names:
                        srna = line[1:]
                    else:
                        if srna in srnas["RNAup"].keys():
                            srnas["RNAup"][srna].append({"target": line[1:],
                                                         "energy": 0})
                        else:
                            srnas["RNAup"][srna] = []
                            srnas["RNAup"][srna].append({"target": line[1:],
                                                         "energy": 0})
                else:
                    detect_energy(line, srnas["RNAup"][srna][-1])
                    extract_pos(line, srnas["RNAup"][srna][-1], "RNAplex")
    return srnas


def import_merge(merges, name, srna_info, srna_plex, srna_up, target_info,
                 pre_target, suf_target):
    ps_start, ps_end = mod_srna_tar_pos(srna_info, srna_plex["srna_pos"],
                                        "srna", pre_target, suf_target)
    pt_start, pt_end = mod_srna_tar_pos(target_info, srna_plex["tar_pos"],
                                        "tar", pre_target, suf_target)
    us_start, us_end = mod_srna_tar_pos(srna_info, srna_up["srna_pos"],
                                        "srna", pre_target, suf_target)
    ut_start, ut_end = mod_srna_tar_pos(target_info, srna_up["tar_pos"],
                                        "tar", pre_target, suf_target)
    merges.append([name, srna_info.seq_id,
                   "-".join([str(srna_info.start), str(srna_info.end)]),
                   "-".join([str(ps_start), str(ps_end)]),
                   "-".join([str(us_start), str(us_end)]),
                   srna_info.strand, srna_plex["target"],
                   "-".join([str(target_info.start), str(target_info.end)]),
                   "-".join([str(pt_start), str(pt_end)]),
                   "-".join([str(ut_start), str(ut_end)]),
                   target_info.strand,
                   str(srna_plex["energy"]), str(srna_plex["rank"]),
                   str(srna_up["energy"]), str(srna_up["rank"])])


def get_srna_name(gffs, srna):
    detect_name = False
    srna_info = None
    for gff in gffs:
        if (gff.attributes["ID"] == srna) and \
           ("Name" in gff.attributes.keys()):
            name = gff.attributes["Name"]
            detect_name = True
            srna_info = gff
            break
        elif (gff.attributes["ID"] == srna):
            srna_info = gff
            break
    if not detect_name:
        name = srna
    return (name, srna_info)


def get_target_info(gffs, target):
    name = target.split("|")
    for gff in gffs:
        if "locus_tag" in gff.attributes.keys():
            if gff.attributes["locus_tag"] == name[0]:
                target_info = gff
                return target_info
        elif (gff.feature + ":" + str(gff.start) + "-" +
              str(gff.end) + "_" + gff.strand) == name[0]:
            target_info = gff
            return target_info
        elif "ID" in gff.attributes.keys():
            if gff.attributes["ID"] == name[0]:
                target_info = gff
                return target_info


def print_file(merges, out):
    merges = sorted(merges, key=lambda k: (k[0], int(k[12])))
    for merge in merges:
        out.write("\t".join(merge) + "\n")


def read_gff(filename):
    gffs = []
    for entry in Gff3Parser().entries(open(filename)):
        gffs.append(entry)
    gffs = sorted(gffs, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    return gffs


def print_title(out):
    out.write("\t".join(["sRNA", "strain", "sRNA_position",
                         "sRNA_interacted_position_RNAplex",
                         "sRNA_interacted_position_RNAup", "sRNA_strand",
                         "target", "target_position",
                         "target_interacted_position_RNAplex",
                         "target_interacted_position_RNAup", "target_strand",
                         "energy_RNAplex", "rank_RNAplex",
                         "energy_RNAup", "rank_RNAup"]) + "\n")


def merge_base_rnaplex(srnas, srna_gffs, args_tar, gffs, merges):
    '''merge the results based on the ranking of RNAplex'''
    overlaps = []
    for srna, srna_plexs in srnas["RNAplex"].items():
        srna_datas = get_srna_name(srna_gffs, srna)
        name = srna_datas[0]
        srna_info = srna_datas[1]
        for srna_plex in srna_plexs:
            if "rank" in srna_plex.keys():
                if ("print" not in srna_plex.keys()) and (
                        srna_plex["rank"] <= args_tar.top):
                    for srna_up in srnas["RNAup"][srna]:
                        if "rank" in srna_up.keys():
                            if (srna_plex["target"] == srna_up["target"]) and (
                                (srna_plex["rank"] <= args_tar.top) and (
                                 srna_up["rank"] <= args_tar.top)) and (
                                 "print" not in srna_up.keys()):
                                target_info = get_target_info(
                                    gffs, srna_plex["target"])
                                import_merge(
                                    overlaps, name, srna_info, srna_plex,
                                    srna_up, target_info,
                                    args_tar.tar_start, args_tar.tar_end)
                                srna_plex["print"] = True
                                srna_up["print"] = True
                                import_merge(
                                    merges, name, srna_info, srna_plex,
                                    srna_up, target_info,
                                    args_tar.tar_start, args_tar.tar_end)
                            elif (srna_plex["target"] ==
                                  srna_up["target"]) and (
                                      "print" not in srna_up.keys()):
                                target_info = get_target_info(
                                    gffs, srna_plex["target"])
                                import_merge(
                                    merges, name, srna_info, srna_plex,
                                    srna_up, target_info,
                                    args_tar.tar_start, args_tar.tar_end)
                                srna_plex["print"] = True
    return overlaps


def merge_base_rnaup(srnas, srna_gffs, args_tar, gffs, merges):
    '''merge the results based on the ranking of RNAup'''
    for srna, srna_ups in srnas["RNAup"].items():
        srna_datas = get_srna_name(srna_gffs, srna)
        name = srna_datas[0]
        srna_info = srna_datas[1]
        for srna_up in srna_ups:
            if "rank" in srna_up.keys():
                if ("print" not in srna_up.keys()) and (
                        srna_up["rank"] <= args_tar.top) and (
                        srna in srnas["RNAplex"].keys()):
                    for srna_plex in srnas["RNAplex"][srna]:
                        if "rank" in srna_plex.keys():
                            if (srna_plex["target"] == srna_up["target"]) and (
                                    "print" not in srna_plex.keys()):
                                target_info = get_target_info(
                                    gffs, srna_up["target"])
                                import_merge(
                                    merges, name, srna_info, srna_plex,
                                    srna_up, target_info,
                                    args_tar.tar_start, args_tar.tar_end)
                                srna_up["print"] = True


def merge_srna_target(rnaplex, rnaup, args_tar, out_rnaplex, out_rnaup, output,
                      out_overlap, srna_gff_file, annotation_gff):
    '''merge the results of RNAup and RNAplex'''
    merges = []
    srna_gffs = read_gff(srna_gff_file)
    gffs = read_gff(annotation_gff)
    srnas = read_table(srna_gffs, rnaplex, rnaup)
    if out_rnaplex is not None:
        out_p = open(out_rnaplex, "w")
        print_rank_one(srnas, out_p, "RNAplex", gffs, srna_gffs, args_tar)
    if out_rnaup is not None:
        out_u = open(out_rnaup, "w")
        print_rank_one(srnas, out_u, "RNAup", gffs, srna_gffs, args_tar)
    if (out_rnaup is not None) and (out_rnaplex is not None):
        out_m = open(output, "w")
        out_o = open(out_overlap, "w")
        print_title(out_m)
        print_title(out_o)
        print("Merging now...")
        overlaps = merge_base_rnaplex(srnas, srna_gffs, args_tar, gffs, merges)
        merge_base_rnaup(srnas, srna_gffs, args_tar, gffs, merges)
        print_file(merges, out_m)
        print_file(overlaps, out_o)
