from annogesiclib.helper import Helper
from annogesiclib.gff3 import Gff3Parser


def get_feature(gene, file_type):
    if "Name" in gene.attributes.keys():
        feature = gene.attributes["Name"]
    elif "locus_tag" in gene.attributes.keys():
        feature = gene.attributes["locus_tag"]
    else:
        if file_type != "tran":
            file_type = gene.feature
        strand = Helper().get_strand_name(gene.strand)
        feature = "".join([file_type, ":", str(gene.start),
                           "-", str(gene.end), "_", strand])
    pos = str(gene.start) + "-" + str(gene.end)
    return feature, pos


def import_data(f_1, f_2, start, end, file_type):
    if f_1 != "terminal":
        feature_1, pos_1 = get_feature(f_1, file_type)
        strain = f_1.seq_id
    else:
        feature_1 = "NA"
        pos_1 = "NA"
    if f_2 != "terminal":
        feature_2, pos_2 = get_feature(f_2, file_type)
        strain = f_2.seq_id
    else:
        feature_2 = "NA"
        pos_2 = "NA"
    return {"strain": strain, "start": start, "end": end,
            "parent_p": feature_1, "parent_m": feature_2, "print": False,
            "p_pos": pos_1, "m_pos": pos_2}


def get_terminal(genes, inters, gene_len, type_, file_type):
    '''deal with the intergenic region which located at two ends'''
    if type_ == "start":
        for gene in genes:
            if (gene.strand == "-"):
                inters.append(import_data(gene.strand, gene.seq_id, 1,
                                          gene.start, file_type))
                break
    elif type_ == "end":
        for gene in reversed(genes):
            if (gene.strand == "+"):
                inters.append(import_data(gene.strand, gene.seq_id, gene.end,
                                          gene_len, file_type))
                break


def get_inter(features, seq, file_type):
    '''compare features and get the intergenic region'''
    inters = []
    first = True
    pre_strain = ""
    pre_feature1 = None
    sort_features = sorted(features, key=lambda x: (x.seq_id, x.start,
                                                    x.end, x.strand))
    for feature1 in sort_features:
        if pre_strain != feature1.seq_id:
            if feature1.strand == "-":
                inters.append(import_data("terminal", feature1, 1,
                                          feature1.start, file_type))
            if not first:
                if pre_feature1.strand == "+":
                    inters.append(import_data(pre_feature1, "terminal",
                                  pre_feature1.start,
                                  len(seq[pre_feature1.seq_id]), file_type))
            if first:
                first = False
        for feature2 in sort_features:
            if (feature1.seq_id == feature2.seq_id) and (
                    feature1.end < feature2.end) and (
                    feature1.end >= feature2.start) and (
                    feature1.strand == feature2.strand):
                break
            elif (feature1.seq_id == feature2.seq_id) and (
                    feature1.end <= feature2.start) and (
                    feature1.strand == "+"):
                if feature1.strand == feature2.strand:
                    break
                else:
                    inters.append(import_data(feature1, feature2, feature1.end,
                                              feature2.start, file_type))
                    break
        pre_feature1 = feature1
        pre_strain = feature1.seq_id
    return inters


def import_merge(id_, strain, start, end, parent_p, parent_m, p_pos, m_pos):
    return {"ID": "_".join(["inter_" + str(id_)]), "strain": strain,
            "start": start, "end": end, "parent_p": parent_p,
            "parent_m": parent_m, "print": False, "p_pos": p_pos,
            "m_pos": m_pos}


def get_overlap_inters(inter1, inter2, merges, id_):
    if (inter1["end"] < inter2["end"]) and (
            inter1["end"] > inter2["start"]) and (
            inter1["start"] <= inter2["start"]):
        merges.append(import_merge(id_, inter1["strain"], inter2["start"],
                      inter1["end"], inter2["parent_p"], inter1["parent_m"],
                      inter2["p_pos"], inter1["m_pos"]))
        inter1["print"] = True
        inter2["print"] = True
        id_ += 1
    elif (inter1["start"] > inter2["start"]) and (
            inter1["start"] < inter2["end"]) and (
            inter1["end"] >= inter2["end"]):
        merges.append(import_merge(id_, inter1["strain"], inter1["start"],
                      inter2["end"], inter1["parent_p"], inter2["parent_m"],
                      inter1["p_pos"], inter2["m_pos"]))
        inter1["print"] = True
        inter2["print"] = True
        id_ += 1
    elif (inter1["end"] >= inter2["end"]) and (
            inter1["start"] <= inter2["start"]):
        merges.append(import_merge(id_, inter2["strain"], inter2["start"],
                      inter2["end"], inter2["parent_p"], inter2["parent_m"],
                      inter2["p_pos"], inter2["m_pos"]))
        inter1["print"] = True
        inter2["print"] = True
        id_ += 1
    elif (inter1["end"] <= inter2["end"]) and (
            inter1["start"] >= inter2["start"]):
        merges.append(import_merge(id_, inter1["strain"], inter1["start"],
                      inter1["end"], inter1["parent_p"], inter1["parent_m"],
                      inter1["p_pos"], inter1["m_pos"]))
        inter1["print"] = True
        inter2["print"] = True
        id_ += 1
    return id_


def merge_inter(inters1, inters2):
    '''check and merge the overlap intergenic seq'''
    merges = []
    id_ = 0
    for inter1 in inters1:
        for inter2 in inters2:
            if (inter1["strain"] == inter2["strain"]):
                id_ = get_overlap_inters(inter1, inter2, merges, id_)
        if not inter1["print"]:
            merges.append(import_merge(id_, inter1["strain"], inter1["start"],
                                       inter1["end"], inter1["parent_p"],
                                       inter1["parent_m"], inter1["p_pos"],
                                       inter1["m_pos"]))
            inter1["print"] = True
            id_ += 1
    for inter2 in inters2:
        if not inter2["print"]:
            merges.append(import_merge(id_, inter2["strain"], inter2["start"],
                                       inter2["end"], inter2["parent_p"],
                                       inter2["parent_m"], inter2["p_pos"],
                                       inter2["m_pos"]))
            inter2["print"] = True
            id_ += 1
    sort_merges = sorted(merges, key=lambda x: (x["strain"],
                                                x["start"], x["end"]))
    return sort_merges


def detect_confliction(gc, genes, seq):
    '''check the intergenic region 
    (add the region for folding proper sec str.) 
    will overlap with gene or not'''
    corr_merges = []
    overlap = False
    tmp_start = gc["start"]
    tmp_end = gc["end"]
    if "tran" in gc["parent_p"]:
        if gc["start"] > 80:
            for gene in genes:
                if ((gc["start"] - 80) > gene.start) and (
                        (gc["start"] - 80) < gene.end) and (
                        gene.strand == "+"):
                    tmp_start = gene.end - 30
                    overlap = True
            if not overlap:
                tmp_start = gc["start"] - 80
        else:
            tmp_start = 1
    else:
        if gc["start"] > 30:
            tmp_start = gc["start"] - 30
        else:
            tmp_start = 1
    corr_merges.append(import_merge(gc["ID"], gc["strain"], tmp_start,
                       gc["end"], gc["parent_p"], gc["parent_m"],
                       gc["p_pos"], gc["m_pos"]))
    corr_merges[-1]["strand"] = "+"
    if "tran" in gc["parent_m"]:
        if gc["end"] < (len(seq[gc["strain"]]) - 80):
            for gene in genes:
                if ((gc["end"] + 80) > gene.start) and (
                        (gc["end"] + 80) < gene.end) and (
                        gene.strand == "-"):
                    tmp_end = gene.start + 30
                    overlap = True
            if not overlap:
                tmp_end = gc["end"] + 80
        else:
            tmp_end = len(seq[gc["strain"]])
    else:
        if gc["end"] < len(seq[gc["strain"]]) - 30:
            tmp_end = gc["end"] + 30
        else:
            tmp_end = len(seq[gc["strain"]])
    corr_merges.append(import_merge(gc["ID"], gc["strain"], gc["start"],
                       tmp_end, gc["parent_p"], gc["parent_m"],
                       gc["p_pos"], gc["m_pos"]))
    corr_merges[-1]["strand"] = "-"
    return corr_merges


def read_file(seq_file, tran_file, gff_file):
    seq = {}
    tas = []
    genes = []
    merges = []
    with open(seq_file, "r") as f_h:
        for line in f_h:
            line = line.strip()
            if line.startswith(">"):
                strain = line[1:]
                seq[strain] = ""
            else:
                seq[strain] = seq[strain] + line
    ta_fh = open(tran_file, "r")
    for entry in Gff3Parser().entries(ta_fh):
        tas.append(entry)
        merges.append(entry)
    for entry in Gff3Parser().entries(open(gff_file)):
        if (entry.feature == "gene"):
            genes.append(entry)
            merges.append(entry)
    tas = sorted(tas, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    genes = sorted(genes, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    return seq, tas, merges, genes


def intergenic_seq(seq_file, tran_file, gff_file, out_file):
    '''get intergenic seq'''
    out = open(out_file, "w")
    seq, tas, merges, genes = read_file(seq_file, tran_file, gff_file)
    inter_tas = get_inter(tas, seq, "tran")
    inter_genes = get_inter(genes, seq, "gene")
    merges = merge_inter(inter_tas, inter_genes)
    num = 0
    for tmp_merge in merges:
        corr_merges = detect_confliction(tmp_merge, genes, seq)
        for merge in corr_merges:
            if merge["start"] < merge["end"]:
                if merge["strand"] == "+":
                    inter_seq = Helper().extract_gene(
                        seq[merge["strain"]], merge["start"],
                        merge["end"], "+")
                    out.write(">" + "|".join([
                        "inter_" + str(num), str(merge["start"]),
                        str(merge["end"]), merge["strain"], merge["parent_p"],
                        merge["parent_m"], merge["p_pos"], merge["m_pos"],
                        "+"]) + "\n")
                    out.write(inter_seq + "\n")
                    num += 1
                else:
                    inter_seq = Helper().extract_gene(
                            seq[merge["strain"]],
                            merge["start"], merge["end"], "-")
                    out.write(">" + "|".join([
                        "inter_" + str(num), str(merge["start"]),
                        str(merge["end"]), merge["strain"], merge["parent_p"],
                        merge["parent_m"], merge["p_pos"], merge["m_pos"],
                        "-"]) + "\n")
                    out.write(inter_seq + "\n")
                    num += 1
