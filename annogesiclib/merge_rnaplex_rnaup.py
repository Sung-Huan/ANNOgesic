from annogesiclib.gff3 import Gff3Parser

def detect_energy(line, srna):
    duplex = line.split(" ")[0]
    if ("(" in duplex) or (")" in duplex):
        energy = line.split(":")[1].split("(")[1].split("=")[0]
        energy = float(energy.strip())
        if energy < srna["energy"]:
            srna["energy"] = energy

def print_rank_one(srnas, out, feature, gffs, srna_gffs, top):
    out.write("\t".join(["sRNA", "strain", "sRNA_start", "sRNA_end",
                         "sRNA_strand", "target", "target_start",
                         "target_end", "target_strand",
                         "method", "energy", "rank"]) + "\n")
    for method, srna_datas in srnas.items():
        for srna_name, targets in srna_datas.items():
            rank = 0
            sort_targets = sorted(targets, key=lambda k: (k["energy"]))
            for target in sort_targets:
                if target["energy"] < 0:
                    rank += 1
                    target["rank"] = rank
                    if rank <= top:
                        if method == feature:
                            srna_infos = get_srna_name(srna_gffs, srna_name)
                            name = srna_infos[0]
                            srna_info = srna_infos[1]
                            target_info = get_target_info(gffs, target["target"])
                            if srna_info is not None:
                                out.write("\t".join([name, str(srna_info.seq_id),
                                          str(srna_info.start), str(srna_info.end),
                                          srna_info.strand, target["target"],
                                          str(target_info.start),
                                          str(target_info.end),
                                          target_info.strand, method,
                                          str(target["energy"]),
                                          str(target["rank"])]) + "\n")

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
                                                         "energy" : 0})
                        else:
                            srnas["RNAup"][srna] = []
                            srnas["RNAup"][srna].append({"target": line[1:],
                                                         "energy" : 0})
                else:
                    detect_energy(line, srnas["RNAup"][srna][-1])
    return srnas

def import_merge(merges, name, srna_info, srna_plex, srna_up, target_info):
    merges.append([name, srna_info.seq_id,
                   str(srna_info.start),
                   str(srna_info.end),
                   srna_info.strand, srna_plex["target"],
                   str(target_info.start),
                   str(target_info.end), target_info.strand,
                   "RNAplex", str(srna_plex["energy"]),
                   str(srna_plex["rank"]),
                   "RNAup", str(srna_up["energy"]),
                   str(srna_up["rank"])])

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
        elif "ID" in gff.attributes.keys():
            if gff.attributes["ID"] == name[0]:
                target_info = gff
                return target_info

def print_file(merges, out):
    merges = sorted(merges, key=lambda k: (k[0], int(k[11])))
    for merge in merges:
        out.write("\t".join(merge) + "\n")

def read_gff(filename):
    gffs = []
    for entry in Gff3Parser().entries(open(filename)):
        gffs.append(entry)
    gffs = sorted(gffs, key=lambda k: (k.seq_id, k.start))
    return gffs

def print_title(out):
    out.write("\t".join(["sRNA", "strain", "srna_start", "srna_end",
                         "srna_strand", "target", "target_start",
                         "target_end", "target_strand", "method",
                         "energy", "rank", "methon", "energy", "rank"]) + "\n")

def merge_base_rnaplex(srnas, srna_gffs, top, gffs, merges):
    overlaps = []
    for srna, srna_plexs in srnas["RNAplex"].items():
        srna_datas = get_srna_name(srna_gffs, srna)
        name = srna_datas[0]
        srna_info = srna_datas[1]
        for srna_plex in srna_plexs:
            if "rank" in srna_plex.keys():
                if ("print" not in srna_plex.keys()) and (
                    srna_plex["rank"] <= top):
                    for srna_up in srnas["RNAup"][srna]:
                        if "rank" in srna_up.keys():
                            if (srna_plex["target"] == srna_up["target"]) and (
                                (srna_plex["rank"] <= top) and (
                                 srna_up["rank"] <= top)) and (
                                 "print" not in srna_up.keys()):
                                target_info = get_target_info(gffs, srna_plex["target"])
                                import_merge(overlaps, name, srna_info,
                                             srna_plex,
                                             srna_up, target_info)
                                srna_plex["print"] = True
                                srna_up["print"] = True
                                import_merge(merges, name, srna_info,
                                             srna_plex,
                                             srna_up, target_info)
                            elif (srna_plex["target"] == srna_up["target"]) and (
                                  "print" not in srna_up.keys()):
                                target_info = get_target_info(gffs, srna_plex["target"])
                                import_merge(merges, name, srna_info,
                                             srna_plex,
                                             srna_up, target_info)
                                srna_plex["print"] = True
    return overlaps

def merge_base_rnaup(srnas, srna_gffs, top, gffs, merges):
    for srna, srna_ups in srnas["RNAup"].items():
        srna_datas = get_srna_name(srna_gffs, srna)
        name = srna_datas[0]
        srna_info = srna_datas[1]
        for srna_up in srna_ups:
            if "rank" in srna_up.keys():
                if ("print" not in srna_up.keys()) and (
                    srna_up["rank"] <= top) and (
                    srna in srnas["RNAplex"].keys()):
                    for srna_plex in srnas["RNAplex"][srna]:
                        if "rank" in srna_plex.keys():
                            if (srna_plex["target"] == srna_up["target"]) and (
                                "print" not in srna_plex.keys()):
                                target_info = get_target_info(gffs, srna_up["target"])
                                import_merge(merges, name, srna_info,
                                             srna_plex, srna_up, target_info)
                                srna_up["print"] = True

def merge_srna_target(rnaplex, rnaup, top, out_rnaplex, out_rnaup, output,
                      out_overlap, srna_gff_file, annotation_gff):
    merges = []
    srna_gffs = read_gff(srna_gff_file)
    gffs = read_gff(annotation_gff)
    srnas = read_table(srna_gffs, rnaplex, rnaup)
    if out_rnaplex is not None:
        out_p = open(out_rnaplex, "w")
        print_rank_one(srnas, out_p, "RNAplex", gffs, srna_gffs, top)
    if out_rnaup is not None:
        out_u = open(out_rnaup, "w")
        print_rank_one(srnas, out_u, "RNAup", gffs, srna_gffs, top)
    if (out_rnaup is not None) and (out_rnaplex is not None):
        out_m = open(output, "w")
        out_o = open(out_overlap, "w")
        print_title(out_m)
        print_title(out_o)
        print("Merging now...")
        overlaps = merge_base_rnaplex(srnas, srna_gffs, top, gffs, merges)
        merge_base_rnaup(srnas, srna_gffs, top, gffs, merges)
        print_file(merges, out_m)
        print_file(overlaps, out_o)
