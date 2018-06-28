from annogesiclib.gff3 import Gff3Parser


def detect_energy(line, srna):
    duplex = line.split(" ")[0]
    if ("(" in duplex) or (")" in duplex):
        energy = line.split(":")[1].split("(")[1].split("=")[0]
        energy = float(energy.strip())
        if energy < srna["energy"]:
            srna["energy"] = energy


def get_locus(gene):
    if gene is not None:
        if ("locus_tag" in gene.attributes.keys()):
            locus = gene.attributes["locus_tag"]
        else:
            locus = "NA"
    else:
        locus = "NA"
    return locus


def check_parent_gene(cds, genes):
    target_gene = None
    if target_gene is None:
        for gene in genes:
            if (gene.seq_id == cds.seq_id) and (
                    gene.strand == cds.strand):
                if ((cds.start <= gene.start) and (
                        cds.end >= gene.end)) or (
                        (cds.start >= gene.start) and (
                        cds.end <= gene.end)) or (
                        (cds.start <= gene.start) and (
                        cds.end <= gene.end) and (
                        cds.end >= gene.start)) or (
                        (cds.start >= gene.start) and (
                        cds.start <= gene.end) and (
                        cds.end >= gene.end)):
                    target_gene = gene
                if (cds.start == gene.start) and (
                        cds.end == gene.end):
                    target_gene = gene
                    break
    return target_gene


def mod_srna_tar_pos(gff, pos, type_, pre_target, suf_target, length):
    if "NA" not in pos:
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
    else:
        g_start = "NA"
        g_end = "NA"
    return g_start, g_end

def print_rank_one(srnas, out, feature, gffs, srna_gffs, args_tar, length):
    out.write("\t".join(["sRNA", "Genome", "sRNA_position",
                         "sRNA_interacted_position_" + feature,
                         "sRNA_strand", "Target_gene_ID", "Target_ID",
                         "Target_locus_tag", "Target_position",
                         "Target_interacted_position_" + feature,
                         "Target_strand", "Energy_" + feature,
                         "Rank_" + feature]) + "\n")
    for method, srna_datas in srnas.items():
        for srna_id, targets in srna_datas.items():
            rank = 0
            sort_targets = sorted(targets, key=lambda k: (k["energy"]))
            for target in sort_targets:
                if (target["tar_pos"] != "NA"):
                    if target["energy"] < 0:
                        rank += 1
                        target["rank"] = rank
                        if (rank <= args_tar.top) and (method == feature):
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
        tar_pos = line.split(" : ")[0].strip().split(" ")[-1]
        srna_pos = line.split(" : ")[-1].strip().split(" ")[0]
        srnas["tar_pos"] = tar_pos
        srnas["srna_pos"] = srna_pos

def read_rnaplex(rnaplex, genes, genomes, features, srnas):
    start = False
    count_seq = 0
    with open(rnaplex, "r") as p_h:
        for line in p_h:
            line = line.strip()
            if line.startswith(">"):
                start = True
                count_seq += 1
                if count_seq == 1:
                    tags = line[1:].split("_")
                    target_locus = "_".join(tags[:-3])
                    target_id = tags[-3]
                    detail = "_".join(tags[-2:])
                    gene_id, tar = get_gene_id(detail, target_id,
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

def read_rnaup(rnaup, srna_names, srnas, genes, genomes, features):
    with open(rnaup, "r") as u_h:
        for line in u_h:
            line = line.strip()
            if line.startswith(">"):
                if line[1:] in srna_names:
                    srna = line[1:]
                else:
                    tags = line[1:].split("_")
                    target_locus = "_".join(tags[:-3])
                    target_id = tags[-3]
                    detail = "_".join(tags[-2:])
                    gene_id, tar = get_gene_id(detail, target_id, genes, genomes,
                                               features)
                    if srna in srnas["RNAup"].keys():
                        srnas["RNAup"][srna].append({
                            "target_id": target_id, "target_locus": target_locus,
                            "detail": detail, "energy": 0,
                            "gene_id": gene_id})
                    else:
                        srnas["RNAup"][srna] = []
                        srnas["RNAup"][srna].append({
                            "target_id": target_id, "target_locus": target_locus,
                            "detail": detail, "energy": 0,
                            "gene_id": gene_id})
            else:
                detect_energy(line, srnas["RNAup"][srna][-1])
                extract_pos(line, srnas["RNAup"][srna][-1], "RNAplex")

def read_intarna(intarna, srnas, genes, genomes, features):
    with open(intarna, "r") as i_h:
        for line in i_h:
            inter = line.strip().split(";")
            if inter[0] != "id1":
                if len(inter) == 9:
                    srna = inter[3]
                    tags = inter[0].split("_")
                    if len(tags) >= 4:
                        target_locus = "_".join(tags[:-3])
                        target_id = tags[-3]
                        detail = "_".join(tags[-2:])
                        if (len(tags[0])) != 0:
                            gene_id, tar = get_gene_id(detail, target_id, genes,
                                                       genomes, features)
                            if tar is not None:
                                if (srna not in srnas["IntaRNA"].keys()):
                                    srnas["IntaRNA"][srna] = []
                                srnas["IntaRNA"][srna].append({
                                    "target_id": target_id, "target_locus": target_locus,
                                    "detail": detail, "energy": float(inter[-1]),
                                    "gene_id": gene_id,
                                    "tar_pos": ",".join(inter[1:3]),
                                    "srna_pos": ",".join(inter[4:6])})

def read_table(gffs, rnaplex, rnaup, intarna, genes, genomes, features):
    srnas = {"RNAup": {}, "RNAplex": {}, "IntaRNA": {}}
    srna_names = set()
    for gff in gffs:
        if gff.attributes["ID"] not in srna_names:
            srna_names.add(gff.attributes["ID"])
    if rnaplex is not None:
        read_rnaplex(rnaplex, genes, genomes, features, srnas)
    if rnaup is not None:
        read_rnaup(rnaup, srna_names, srnas, genes, genomes, features)
    if intarna is not None:
        read_intarna(intarna, srnas, genes, genomes, features)
    return srnas


def get_gene_id(detail, tar_id, genes, gffs, features):
    tar = None
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
                strand = detail.split("_")[-1]
                if (gff.start == start) and (gff.end == end) and (
                        gff.strand == strand):
                    tar = gff
                    break
    gene_id = "NA"
    if tar is not None:
        for gene in genes:
            if "Parent" in tar.attributes.keys():
                if "ID" in gene.attributes.keys():
                    if (gene.attributes["ID"] in
                            tar.attributes["Parent"].split(",")):
                        gene_id = gene.attributes["ID"]
                        return gene_id, tar
            if gene_id == "NA":
                if (gene.seq_id == tar.seq_id) and (
                        gene.strand == tar.strand) and (
                        (tar.start == gene.start) and (
                        tar.end == gene.end)):
                    if "ID" in gene.attributes.keys():
                        gene_id = gene.attributes["ID"]
                        return gene_id, tar
        if gene_id == "NA":
            for gene in genes:
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
                            return gene_id, tar
    return gene_id, tar


def append_merge_three_methods(name, srna_info, ps_pos, pt_pos, u_data, i_data,
                               srna_m1, target_info, merges):
    merges.append([name, srna_info.seq_id,
                   "-".join([str(srna_info.start), str(srna_info.end)]),
                   ps_pos, u_data["s_pos"], i_data["s_pos"],
                   srna_info.strand, srna_m1["gene_id"],
                   srna_m1["target_id"], srna_m1["target_locus"],
                   "-".join([str(target_info.start), str(target_info.end)]),
                   pt_pos, u_data["t_pos"], i_data["t_pos"],
                   target_info.strand,
                   str(srna_m1["energy"]), str(srna_m1["rank"]),
                   str(u_data["energy"]), str(u_data["rank"]),
                   str(i_data["energy"]), str(i_data["rank"])])


def import_merge(merges, name, srna_info, srna_m1, srna_m2, srna_m3,
                 target_info, pre_target, suf_target, length, num_method):
    ps_start, ps_end = mod_srna_tar_pos(
            srna_info, srna_m1["srna_pos"], "srna", pre_target,
            suf_target, length)
    pt_start, pt_end = mod_srna_tar_pos(
            target_info, srna_m1["tar_pos"], "tar", pre_target,
            suf_target, length)
    ps_pos = "-".join([str(ps_start), str(ps_end)])
    pt_pos = "-".join([str(pt_start), str(pt_end)])
    u_data = {"s_pos": "NA", "t_pos": "NA", "energy": 1000, "rank": "NA"}
    if (srna_m1["detail"] == srna_m2["detail"]):
        us_start, us_end = mod_srna_tar_pos(
                srna_info, srna_m2["srna_pos"], "srna", pre_target,
                suf_target, length)
        ut_start, ut_end = mod_srna_tar_pos(
                target_info, srna_m2["tar_pos"], "tar", pre_target,
                suf_target, length)
        u_data["s_pos"] = "-".join([str(us_start), str(us_end)])
        u_data["t_pos"] = "-".join([str(ut_start), str(ut_end)])
        u_data["energy"] = srna_m2["energy"]
        u_data["rank"] = srna_m2["rank"]
    i_data = {"s_pos": "NA", "t_pos": "NA", "energy": 1000, "rank": "NA"}
    if srna_m3 is not None:
        if (srna_m1["detail"] == srna_m3["detail"]):
            is_start, is_end = mod_srna_tar_pos(
                    srna_info, srna_m3["srna_pos"], "srna", pre_target,
                    suf_target, length)
            it_start, it_end = mod_srna_tar_pos(
                    target_info, srna_m3["tar_pos"], "tar", pre_target,
                    suf_target, length)
            i_data["s_pos"] = "-".join([str(is_start), str(is_end)])
            i_data["t_pos"] = "-".join([str(it_start), str(it_end)])
            i_data["energy"] = srna_m3["energy"]
            i_data["rank"] = srna_m3["rank"]
        append_merge_three_methods(name, srna_info, ps_pos, pt_pos, u_data,
                                   i_data, srna_m1, target_info, merges)
    else:
        if num_method == 2:
            merges.append([name, srna_info.seq_id,
                           "-".join([str(srna_info.start), str(srna_info.end)]),
                           ps_pos, u_data["s_pos"],
                           srna_info.strand, srna_m1["gene_id"],
                           srna_m1["target_id"], srna_m1["target_locus"],
                           "-".join([str(target_info.start), str(target_info.end)]),
                           pt_pos, u_data["t_pos"],
                           target_info.strand,
                           str(srna_m1["energy"]), str(srna_m1["rank"]),
                           str(u_data["energy"]), str(u_data["rank"])])
        elif num_method == 3:
            append_merge_three_methods(name, srna_info, ps_pos, pt_pos, u_data,
                                       i_data, srna_m1, target_info, merges)


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
    tmp_gff = None
    for gff in gffs:
        if (str(gff.start) + "-" +
            str(gff.end) + "_" + gff.strand) == target["detail"]:
            if gff.feature == "gene":
                target_info = gff
                if ("locus_tag" in gff.attributes.keys()) and (
                        "gene" in gff.attributes.keys()):
                    if (gff.attributes["gene"] not in target["target_locus"]):
                        target["target_locus"] = "|".join([
                            target["target_locus"], gff.attributes["gene"]])
                elif ("locus_tag" in gff.attributes.keys()) and (
                        "Name" in gff.attributes.keys()):
                    if (gff.attributes["Name"] not in target["target_locus"]):
                        target["target_locus"] = "|".join([
                            target["target_locus"], gff.attributes["Name"]])
                return target_info
            else:
                tmp_gff = gff
    if tmp_gff is not None:
        return tmp_gff


def remove_no_rank(merges, index):
    new_merges = []
    for merge in merges:
        if merge[index] != "NA":
            new_merges.append(merge)
    return new_merges


def print_file(merges, out, num_method):
    if num_method == 2:
        merges = remove_no_rank(merges, 14)
        merges = sorted(merges, key=lambda k: (k[0], int(k[14])))
        for merge in merges:
            if float(merge[13]) == 1000:
                merge[13] = "NA"
            if float(merge[15]) == 1000:
                merge[15] = "NA"
            out.write("\t".join(merge) + "\n")
    elif num_method == 3:
        merges = remove_no_rank(merges, 16)
        merges = sorted(merges, key=lambda k: (k[0], int(k[16])))
        for merge in merges:
            if float(merge[15]) == 1000:
                merge[15] = "NA"
            if float(merge[17]) == 1000:
                merge[17] = "NA"
            if float(merge[19]) == 1000:
                merge[19] = "NA"
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


def print_title(out, methods):
    if len(methods) == 2: 
        out.write("\t".join(["sRNA", "Genome", "sRNA_position",
                             "sRNA_interacted_position_" + methods[0],
                             "sRNA_interacted_position_" + methods[1], "sRNA_strand",
                             "Target_gene_ID", "Target_ID", "Target_locus_tag", "Target_position",
                             "Target_interacted_position_" + methods[0],
                             "Target_interacted_position_" + methods[1], "Target_strand",
                             "Energy_" + methods[0], "Rank_" + methods[1],
                             "Energy_" + methods[0], "Rank_" + methods[1]]) + "\n")
    if len(methods) == 3:
        out.write("\t".join(["sRNA", "Genome", "sRNA_position",
                             "sRNA_interacted_position_" + methods[0],
                             "sRNA_interacted_position_" + methods[1],
                             "sRNA_interacted_position_" + methods[2], "sRNA_strand",
                             "Target_gene_ID", "Target_ID", "Target_locus_tag", "Target_position",
                             "Target_interacted_position_" + methods[0],
                             "Target_interacted_position_" + methods[1],
                             "Target_interacted_position_" + methods[2], "Target_strand",
                             "Energy_" + methods[0], "Rank_" + methods[0],
                             "Energy_" + methods[1], "Rank_" + methods[1],
                             "Energy_" + methods[2], "Rank_" + methods[2]]) + "\n")



def merge_three(srnas, method3, srna, detect, srna_m1, top):
    srna_m3 = None
    if method3 is not None:
        for srna_m3 in srnas[method3][srna]:
            if (srna_m1["detail"] == srna_m3["detail"]) and (
                    "rank" in srna_m3.keys()):
                if (srna_m3["rank"] <= top) and (
                        "print" not in srna_m3.keys()) and (
                        srna_m3["tar_pos"] != "NA"):
                    srna_m3["print"] = True
                    detect["3"] = True
                return srna_m3
    else:
        return srna_m3


def check_method3(srna_m3, srna_m1, method3, srna, srnas):
    if (srna_m3 is not None):
        if "rank" not in srna_m3.keys():
            srna_m3["rank"] = "NA"
            srna_m3["energy"] = 1000
        return srna_m3
    else:
        if method3 is not None:
            for srna_m3 in srnas[method3][srna]:
                if (srna_m1["detail"] == srna_m3["detail"]):
                    return srna_m3


def check_non_overlap(detect, methods, srna_m1, srna_m2, srna_m3, gffs,
                      merges, name, srna_info, args_tar, length,
                      method3, srna, srnas):
    if (not detect["2"] and len(methods) == 2) or (
            ((not detect["2"]) or (not detect["3"]))
            and len(methods) == 3):
        if "rank" not in srna_m2.keys():
            srna_m2["rank"] = "NA"
            srna_m2["energy"] = 1000
        srna_m3 = check_method3(srna_m3, srna_m1, method3, srna, srnas)
        target_info = get_target_info(
            gffs, srna_m1)
        import_merge(
            merges, name, srna_info, srna_m1,
            srna_m2, srna_m3, target_info,
            args_tar.tar_start,
            args_tar.tar_end, length, len(methods))
        srna_m1["print"] = True

def merge_result(srnas, srna_gffs, args_tar, gffs, merges, length, methods):
    '''merge the results based on the ranking of RNAplex'''
    overlaps = []
    method1 = methods[0]
    method2 = methods[1]
    method3 = None
    if len(methods) == 3:
        method3 = methods[2]
    for srna, srna_m1s in srnas[method1].items():
        srna_datas = get_srna_name(srna_gffs, srna)
        name = srna_datas[0]
        srna_info = srna_datas[1]
        for srna_m1 in srna_m1s:
            if "rank" in srna_m1.keys():
                if ("print" not in srna_m1.keys()) and (
                        srna_m1["rank"] <= args_tar.top) and (
                        srna_m1["tar_pos"] != "NA"):
                    detect = {"2": False, "3": False}
                    srna_m3 = None
                    for srna_m2 in srnas[method2][srna]:
                        if (srna_m1["detail"] == srna_m2["detail"]):
                            if ("rank" in srna_m2.keys()):
                                if (srna_m2["rank"] <= args_tar.top) and (
                                        "print" not in srna_m2.keys()) and (
                                         srna_m2["tar_pos"] != "NA"):
                                    detect["2"] = True
                                    srna_m3 = merge_three(srnas, method3, srna,
                                                          detect, srna_m1,
                                                          args_tar.top)
                                    if (len(methods) == 2) or (
                                            (len(methods) == 3) and detect["3"]):
                                        target_info = get_target_info(
                                            gffs, srna_m1)
                                        import_merge(
                                            overlaps, name, srna_info, srna_m1,
                                            srna_m2, srna_m3, target_info,
                                            args_tar.tar_start,
                                            args_tar.tar_end, length,
                                            len(methods))
                                        srna_m1["print"] = True
                                        srna_m2["print"] = True
                                        import_merge(
                                            merges, name, srna_info, srna_m1,
                                            srna_m2, srna_m3, target_info,
                                            args_tar.tar_start,
                                            args_tar.tar_end, length,
                                            len(methods))
                            break
                    check_non_overlap(detect, methods, srna_m1, srna_m2,
                                      srna_m3, gffs, merges, name, srna_info,
                                      args_tar, length, method3, srna, srnas)
    return overlaps


def compare_rest(srnas, rest, srna_last, srna):
    srna_m3 = None
    for srna_m3 in srnas[rest][srna]:
        if (srna_last["detail"] == srna_m3["detail"]) and (
                "print" not in srna_m3.keys()):
            if ("rank" not in srna_m3.keys()):
                srna_m3["rank"] = "NA"
                srna_m3["energy"] = 1000
            return srna_m3
    return srna_m3


def merge_last(srnas, srna_gffs, args_tar, gffs, merges, length,
               method, ref, num_method, rest, switch):
    '''merge the results based on the ranking of RNAup'''
    for srna, srnas_last in srnas[method].items():
        srna_datas = get_srna_name(srna_gffs, srna)
        name = srna_datas[0]
        srna_info = srna_datas[1]
        for srna_last in srnas_last:
            if "rank" in srna_last.keys():
                if srna_last["rank"] != "NA":
                    if ("print" not in srna_last.keys()) and (
                            srna_last["rank"] <= args_tar.top) and (
                            srna in srnas[ref].keys()) and (
                            srna_last["tar_pos"] != "NA"):
                        for srna_ref in srnas[ref][srna]:
                            if (srna_ref["detail"] == srna_last["detail"]) and (
                                    "print" not in srna_ref.keys()):
                                if ("rank" not in srna_ref.keys()):
                                    srna_ref["rank"] = "NA"
                                    srna_ref["energy"] = 1000
                                else:
                                    target_info = get_target_info(
                                        gffs, srna_last)
                                    if num_method == 2:
                                        import_merge(
                                            merges, name, srna_info, srna_ref,
                                            srna_last, None, target_info,
                                            args_tar.tar_start,
                                            args_tar.tar_end, length, num_method)
                                    elif num_method == 3:
                                        srna_m3 = compare_rest(
                                            srnas, rest, srna_last, srna)
                                        if switch:
                                            import_merge(
                                                merges, name, srna_info,
                                                srna_ref, srna_m3, srna_last,
                                                target_info, args_tar.tar_start,
                                                args_tar.tar_end, length,
                                                num_method)
                                        else:
                                            import_merge(
                                                merges, name, srna_info,
                                                srna_ref, srna_last, srna_m3,
                                                target_info, args_tar.tar_start,
                                                args_tar.tar_end, length,
                                                num_method)
                                        if srna_m3 is not None:
                                            srna_m3["print"] = True
                                    srna_last["print"] = True
                                    srna_ref["print"] = True


def read_fasta(seq_file):
    length = 0
    with open(seq_file) as fh:
        for line in fh:
            line = line.strip()
            if not line.startswith(">"):
                length = length + len(line)
    return length


def merge_srna_target(rnaplex, rnaup, intarna, args_tar, out_rnaplex,
                      out_rnaup, out_intarna, seq_file, output,
                      out_overlap, srna_gff_file, annotation_gff):
    '''merge the results of RNAup and RNAplex'''
    length = read_fasta(seq_file)
    merges = []
    methods = []
    srna_gffs, NA = read_gff(srna_gff_file)
    gffs, genes = read_gff(annotation_gff)
    srnas = read_table(srna_gffs, rnaplex, rnaup, intarna, genes, gffs,
                       args_tar.features)
    if out_rnaplex is not None:
        print("Ranking for RNAplex")
        methods.append("RNAplex")
        out_p = open(out_rnaplex, "w")
        print_rank_one(srnas, out_p, "RNAplex", gffs, srna_gffs, args_tar,
                       length)
    if out_rnaup is not None:
        print("Ranking for RNAup")
        methods.append("RNAup")
        out_u = open(out_rnaup, "w")
        print_rank_one(srnas, out_u, "RNAup", gffs, srna_gffs, args_tar,
                       length)
    if out_intarna is not None:
        print("Ranking for IntaRNA")
        methods.append("IntaRNA")
        out_i = open(out_intarna, "w")
        print_rank_one(srnas, out_i, "IntaRNA", gffs, srna_gffs, args_tar,
                       length)
    if (len(args_tar.program) >= 2):
        out_m = open(output, "w")
        out_o = open(out_overlap, "w")
        print_title(out_m, methods)
        print_title(out_o, methods)
        print("Merging now...")
        overlaps = merge_result(srnas, srna_gffs, args_tar, gffs,
                                merges, length, methods)
        if len(methods) == 2:
            merge_last(srnas, srna_gffs, args_tar, gffs, merges, length,
                       methods[1], methods[0], 2, None, False)
        elif len(methods) == 3:
            merge_last(srnas, srna_gffs, args_tar, gffs, merges, length,
                       methods[1], methods[0], 3, methods[2], False)
            merge_last(srnas, srna_gffs, args_tar, gffs, merges, length,
                       methods[2], methods[0], 3, methods[1], True)
        print_file(merges, out_m, len(methods))
        print_file(overlaps, out_o, len(methods))
