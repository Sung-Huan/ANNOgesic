from annogesiclib.gff3 import Gff3Parser


def detect_energy(line, srna):
    duplex = line.split(" ")[0]
    if ("(" in duplex) or (")" in duplex):
        energy = line.split(":")[1].split("(")[1].split("=")[0]
        energy = float(energy.strip())
        if energy < srna["energy"]:
            srna["energy"] = energy


def mod_srna_tar_pos(gff, pos, type_, pre_target, suf_target, length):
    start = int(pos.split(",")[0])
    end = int(pos.split(",")[-1])
    if (gff.strand == "+"):
        if type_ == "srna":
            g_start = gff.start + start - 1
            g_end = gff.start + end - 1
        else:
            if (gff.start - pre_target) <= 0:
                g_start = start
                g_end = end
            else:
                g_start = gff.start - pre_target + start - 1
                g_end = gff.start - pre_target + end - 1
            if (gff.start - pre_target + end - 1) > length:
                g_end = length
    else:
        if type_ == "srna":
            g_end = gff.end - start + 1
            g_start = gff.end - end + 1
        else:
            if (gff.end + pre_target) > length:
                g_start = length - end + 1
                g_end = length - start + 1
            else:
                g_start = gff.end + pre_target - end + 1
                g_end = gff.end + pre_target - start + 1
            if (gff.end + pre_target - end + 1) <= 0:
                g_start = 1
    return g_start, g_end

def print_rank_one(srnas, out, feature, gffs, srna_gffs, args_tar, length):
    out.write("\t".join(["sRNA", "strain", "sRNA_position",
                         "sRNA_interacted_position_" + feature,
                         "sRNA_strand", "target_gene_ID", "target_ID",
                         "target_locus_tag", "target_position",
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
                            target_info = get_target_info(gffs, target)
                            s_start, s_end = mod_srna_tar_pos(
                                    srna_info, target["srna_pos"], "srna",
                                    args_tar.tar_start, args_tar.tar_end,
                                    length)
                            t_start, t_end = mod_srna_tar_pos(
                                    target_info, target["tar_pos"], "tar",
                                    args_tar.tar_start, args_tar.tar_end,
                                    length)
                            if srna_info is not None:
                                out.write("\t".join([
                                    name, str(srna_info.seq_id),
                                    "-".join([str(srna_info.start),
                                              str(srna_info.end)]),
                                    "-".join([str(s_start), str(s_end)]),
                                    srna_info.strand, target["gene_id"],
                                    target["target_id"],
                                    target["target_locus"],
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

def read_table(gffs, rnaplex, rnaup, genes, genomes, features):
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
                        tags = line[1:].split("|")
                        target_locus = tags[0]
                        target_id = tags[1]
                        detail = tags[-1]
                        gene_id = get_gene_id(detail, target_id,
                                              genes, genomes, features)
                    elif count_seq == 2:
                        srna = line[1:]
                        if srna not in srnas["RNAplex"].keys():
                            srnas["RNAplex"][srna] = []
                        srnas["RNAplex"][srna].append({
                            "target_id": target_id, "gene_id": gene_id,
                            "target_locus": target_locus,
                            "detail": detail, "energy": 0})
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
                        tags = line[1:].split("|")
                        gene_id = get_gene_id(tags[-1], tags[1], genes, genomes,
                                              features)
                        if srna in srnas["RNAup"].keys():
                            srnas["RNAup"][srna].append({
                                "target_id": tags[1], "target_locus": tags[0],
                                "detail": tags[-1], "energy": 0,
                                "gene_id": gene_id})
                        else:
                            srnas["RNAup"][srna] = []
                            srnas["RNAup"][srna].append({
                                "target_id": tags[1], "target_locus": tags[0],
                                "detail": tags[-1], "energy": 0,
                                "gene_id": gene_id})
                else:
                    detect_energy(line, srnas["RNAup"][srna][-1])
                    extract_pos(line, srnas["RNAup"][srna][-1], "RNAplex")
    return srnas


def get_gene_id(detail, tar_id, genes, gffs, features):
    for gff in gffs:
        if gff.feature in features:
            if tar_id != "NA":
                if "ID" in gff.attributes.keys():
                    if gff.attributes["ID"] == tar_id:
                        tar = gff
                        break
            else:
                start = int(detail.split("-")[0])
                end = int(detail.split("-")[-1].split("_")[0])
                strand = int(detail.split("_")[-1])
                if (gff.start == start) and (gff.end == end) and (
                        gff.strand == strand):
                    tar = gff
                    break
    gene_id = "NA"
    for gene in genes:
        if "Parent" in tar.attributes.keys():
            if "ID" in gene.attributes.keys():
                if (gene.attributes["ID"] in
                        tar.attributes["Parent"].split(",")):
                    gene_id = gene.attributes["ID"]
                    return gene_id
        else:
            if (gene.seq_id == tar.seq_id) and (
                    gene.strand == tar.strand):
                if ((tar.start <= gene.start) and (
                        tar.end >= gene.end)) or (
                        (tar.start >= gene.start) and (
                        tar.end <= gene.end)) or (
                        (tar.start <= gene.start) and (
                        tar.end <= gene.end) and (
                        tar.end >= gene.start)) or (
                        (tar.start >= gene.start) and (
                        tar.start <= gene.end) and (
                        tar.end >= gene.end)):
                    if "ID" in gene.attributes.keys():
                        gene_id = gene.attributes["ID"]
                        return gene_id
    return gene_id


def import_merge(merges, name, srna_info, srna_plex, srna_up, target_info,
                 pre_target, suf_target, length):
    ps_start, ps_end = mod_srna_tar_pos(
            srna_info, srna_plex["srna_pos"], "srna", pre_target,
            suf_target, length)
    pt_start, pt_end = mod_srna_tar_pos(
            target_info, srna_plex["tar_pos"], "tar", pre_target,
            suf_target, length)
    us_start, us_end = mod_srna_tar_pos(
            srna_info, srna_up["srna_pos"], "srna", pre_target,
            suf_target, length)
    ut_start, ut_end = mod_srna_tar_pos(
            target_info, srna_up["tar_pos"], "tar", pre_target,
            suf_target, length)
    merges.append([name, srna_info.seq_id,
                   "-".join([str(srna_info.start), str(srna_info.end)]),
                   "-".join([str(ps_start), str(ps_end)]),
                   "-".join([str(us_start), str(us_end)]),
                   srna_info.strand, srna_plex["gene_id"],
                   srna_plex["target_id"], srna_plex["target_locus"],
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
    for gff in gffs:
        if (str(gff.start) + "-" +
            str(gff.end) + "_" + gff.strand) == target["detail"]:
            target_info = gff
            if ("locus_tag" in gff.attributes.keys()) and (
                    "gene" in gff.attributes.keys()):
                if (gff.attributes["gene"] not in target["target_locus"]):
                    target["target_locus"] = "|".join([target["target_locus"],
                                                       gff.attributes["gene"]])
            elif ("locus_tag" in gff.attributes.keys()) and (
                    "Name" in gff.attributes.keys()):
                if (gff.attributes["Name"] not in target["target_locus"]):
                    target["target_locus"] = "|".join([target["target_locus"],
                                                       gff.attributes["Name"]])
            return target_info


def print_file(merges, out):
    merges = sorted(merges, key=lambda k: (k[0], int(k[14])))
    for merge in merges:
        out.write("\t".join(merge) + "\n")


def read_gff(filename):
    gffs = []
    genes = []
    for entry in Gff3Parser().entries(open(filename)):
        if entry.feature == "gene":
            genes.append(entry)
        gffs.append(entry)
    gffs = sorted(gffs, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    if len(genes) != 0:
        genes = sorted(genes, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    return gffs, genes


def print_title(out):
    out.write("\t".join(["sRNA", "strain", "sRNA_position",
                         "sRNA_interacted_position_RNAplex",
                         "sRNA_interacted_position_RNAup", "sRNA_strand",
                         "target_gene_ID", "target_ID", "target_locus_tag", "target_position",
                         "target_interacted_position_RNAplex",
                         "target_interacted_position_RNAup", "target_strand",
                         "energy_RNAplex", "rank_RNAplex",
                         "energy_RNAup", "rank_RNAup"]) + "\n")


def merge_base_rnaplex(srnas, srna_gffs, args_tar, gffs, merges, length):
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
                        if (srna_plex["detail"] == srna_up["detail"]) and (
                                "rank" in srna_up.keys()):
                            if (srna_up["rank"] <= args_tar.top) and (
                                    "print" not in srna_up.keys()):
                                target_info = get_target_info(
                                    gffs, srna_plex)
                                import_merge(
                                    overlaps, name, srna_info, srna_plex,
                                    srna_up, target_info, args_tar.tar_start,
                                    args_tar.tar_end, length)
                                srna_plex["print"] = True
                                srna_up["print"] = True
                                import_merge(
                                    merges, name, srna_info, srna_plex,
                                    srna_up, target_info, args_tar.tar_start,
                                    args_tar.tar_end, length)
                            elif ("print" not in srna_up.keys()):
                                target_info = get_target_info(
                                    gffs, srna_plex)
                                import_merge(
                                    merges, name, srna_info, srna_plex,
                                    srna_up, target_info, args_tar.tar_start,
                                    args_tar.tar_end, length)
                                srna_plex["print"] = True
                        elif (srna_plex["detail"] == srna_up["detail"]) and (
                                "rank" not in srna_up.keys()):
                            srna_up["rank"] = "NA"
                            srna_up["energy"] = "NA"
                            if ("print" not in srna_up.keys()):
                                target_info = get_target_info(
                                    gffs, srna_plex)
                                import_merge(
                                    merges, name, srna_info, srna_plex,
                                    srna_up, target_info, args_tar.tar_start,
                                    args_tar.tar_end, length)
                                srna_plex["print"] = True
    return overlaps


def merge_base_rnaup(srnas, srna_gffs, args_tar, gffs, merges, length):
    '''merge the results based on the ranking of RNAup'''
    for srna, srna_ups in srnas["RNAup"].items():
        srna_datas = get_srna_name(srna_gffs, srna)
        name = srna_datas[0]
        srna_info = srna_datas[1]
        for srna_up in srna_ups:
            if "rank" in srna_up.keys():
                if srna_up["rank"] != "NA":
                    if ("print" not in srna_up.keys()) and (
                            srna_up["rank"] <= args_tar.top) and (
                            srna in srnas["RNAplex"].keys()):
                        for srna_plex in srnas["RNAplex"][srna]:
                            if (srna_plex["detail"] == srna_up["detail"]):
                                if ("rank" in srna_plex.keys()) and (
                                        "print" not in srna_plex.keys()):
                                    target_info = get_target_info(
                                        gffs, srna_up)
                                    import_merge(
                                      merges, name, srna_info, srna_plex,
                                      srna_up, target_info, args_tar.tar_start,
                                      args_tar.tar_end, length)
                                    srna_up["print"] = True
                                elif ("rank" not in srna_plex.keys()) and (
                                        "print" not in srna_plex.keys()):
                                    srna_plex["rank"] = "NA"
                                    srna_plex["energy"] = "NA"
                                    import_merge(
                                      merges, name, srna_info, srna_plex,
                                      srna_up, target_info, args_tar.tar_start,
                                      args_tar.tar_end, length)
                                    srna_up["print"] = True


def read_fasta(seq_file):
    length = 0
    with open(seq_file) as fh:
        for line in fh:
            line = line.strip()
            if not line.startswith(">"):
                length = length + len(line)
    return length


def merge_srna_target(rnaplex, rnaup, args_tar, out_rnaplex, out_rnaup,
                      seq_file, output, out_overlap, srna_gff_file,
                      annotation_gff):
    '''merge the results of RNAup and RNAplex'''
    length = read_fasta(seq_file)
    merges = []
    srna_gffs, NA = read_gff(srna_gff_file)
    gffs, genes = read_gff(annotation_gff)
    srnas = read_table(srna_gffs, rnaplex, rnaup, genes, gffs,
                       args_tar.features)
    if out_rnaplex is not None:
        out_p = open(out_rnaplex, "w")
        print_rank_one(srnas, out_p, "RNAplex", gffs, srna_gffs, args_tar,
                       length)
    if out_rnaup is not None:
        out_u = open(out_rnaup, "w")
        print_rank_one(srnas, out_u, "RNAup", gffs, srna_gffs, args_tar,
                       length)
    if (out_rnaup is not None) and (out_rnaplex is not None):
        out_m = open(output, "w")
        out_o = open(out_overlap, "w")
        print_title(out_m)
        print_title(out_o)
        print("Merging now...")
        overlaps = merge_base_rnaplex(srnas, srna_gffs, args_tar, gffs,
                                      merges, length)
        merge_base_rnaup(srnas, srna_gffs, args_tar, gffs, merges, length)
        print_file(merges, out_m)
        print_file(overlaps, out_o)
