import math
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.helper import Helper


def import_candidate(cands, term_features, strain, start, end, ut, name,
                     total_length, strand, parent_p, parent_m, p_pos, m_pos):
    cands.append({"strain": strain, "start": start, "end": end, "print": False,
                  "ut": ut, "name": name, "miss": term_features["real_miss"],
                  "loop": term_features["loop"], "length": total_length,
                  "r_stem": term_features["r_stem"], "strand": strand,
                  "l_stem": term_features["l_stem"], "parent_p": parent_p,
                  "parent_m": parent_m, "detect_p": False, "detect_m": False,
                  "p_pos": p_pos, "m_pos": m_pos})


def get_feature(gene):
    if "Name" in gene.attributes.keys():
        feature = gene.attributes["Name"]
    elif "locus_tag" in gene.attributes.keys():
        feature = gene.attributes["locus_tag"]
    else:
        strand = Helper().get_strand_name(gene.strand)
        feature = "".join([gene.feature, ":", str(gene.start),
                           "-", str(gene.end), "_", strand])
    return feature


def check_miss(cand1, cand2, cutoff_miss):
    stem_len = (cand2["r_stem"] + cand2["l_stem"] - cand2["miss"])
    if (float(cand2["miss"]) / float(stem_len)) <= cutoff_miss:
        cand1["miss"] = cand2["miss"]
        cand1["r_stem"] = cand2["r_stem"]
        cand1["l_stem"] = cand2["l_stem"]


def filter_term(cands, terms, miss_rate):
    '''remove the low possibilty terminator'''
    cutoff_miss = miss_rate
    for cand1 in cands:
        stem_len = (cand1["r_stem"] + cand1["l_stem"] - cand1["miss"])
        if (not cand1["print"]) and \
           ((float(cand1["miss"]) / float(stem_len)) <= cutoff_miss):
            tmp_term = cand1.copy()
            for cand2 in cands:
                if (tmp_term["strain"] == cand2["strain"]) and (
                        tmp_term["miss"] >= cand2["miss"]):
                    if (tmp_term["start"] >= cand2["start"]) and (
                            tmp_term["start"] < cand2["end"]) and (
                            tmp_term["end"] > cand2["end"]):
                        tmp_term["start"] = cand2["start"]
                        check_miss(tmp_term, cand2, cutoff_miss)
                        cand2["print"] = True
                    elif (cand2["start"] > tmp_term["start"]) and (
                            cand2["start"] < tmp_term["end"]) and (
                            cand2["end"] >= tmp_term["end"]):
                        tmp_term["end"] = cand2["end"]
                        check_miss(tmp_term, cand2, cutoff_miss)
                        cand2["print"] = True
                    elif (tmp_term["start"] >= cand2["start"]) and (
                            tmp_term["end"] <= cand2["end"]):
                        tmp_term["start"] = cand2["start"]
                        tmp_term["end"] = cand2["end"]
                        check_miss(tmp_term, cand2, cutoff_miss)
                        cand2["print"] = True
                    elif (cand2["start"] >= tmp_term["start"]) and (
                            cand2["end"] <= tmp_term["end"]):
                        cand2["print"] = True
                        check_miss(tmp_term, cand2, cutoff_miss)
            terms.append(tmp_term)


def check_sec(sec, nts):
    '''check the criteria of sec str of terminator'''
    term_features = {"st_pos": 0, "rights": 0, "lefts": 0,
                     "tmp_miss": 0, "real_miss": 0, "loop": 0,
                     "r_stem": 0, "l_stem": 0}
    detects = {"detect_r": False, "detect_l": False,
               "conflict": False}
    for s_t in reversed(sec[0:nts]):
        term_features["st_pos"] += 1
        if s_t == ")":
            if not detects["detect_l"]:
                term_features["rights"] += 1
                term_features["real_miss"] = term_features["tmp_miss"]
                detects["detect_r"] = True
            else:
                detects["conflict"] = True
                break
        elif s_t == ".":
            if detects["detect_r"]:
                term_features["tmp_miss"] += 1
        elif s_t == "(":
            term_features["lefts"] += 1
            if not detects["detect_l"]:
                term_features["loop"] = (
                    term_features["tmp_miss"] - term_features["real_miss"])
                term_features["tmp_miss"] = term_features["real_miss"]
                term_features["r_stem"] = (
                    term_features["rights"] + term_features["real_miss"])
            else:
                term_features["real_miss"] = term_features["tmp_miss"]
            detects["detect_l"] = True
            if term_features["lefts"] == term_features["rights"]:
                break
    return term_features, detects


def detect_candidates(seq, sec, name, strain, start, end, parent_p, parent_m,
                      strand, args_term, p_pos, m_pos):
    '''check the criteria of sec str of terminator'''
    term_len = 2 * args_term.max_stem + 2 * (
               args_term.max_stem * args_term.miss_rate) + args_term.max_loop
    cands = []
    for nts in range(0, len(seq) - 6):
        ut = 0
        for nt in seq[nts:nts + args_term.range_u]:
            if (nt == "U") or (nt == "T"):
                ut += 1
        if (ut >= args_term.at_tail) and (nts > 10):
            if sec[nts - 1] == ")":
                term_features = {"st_pos": 0, "rights": 0, "lefts": 0,
                                 "tmp_miss": 0, "real_miss": 0, "loop": 0,
                                 "r_stem": 0, "l_stem": 0}
                detects = {"detect_r": False, "detect_l": False,
                           "conflict": False}
                term_features, detects = check_sec(sec, nts)
                if detects["conflict"] is False:
                    total_length = (
                        (nts) - (nts - term_features["st_pos"] + 1) + 1)
                    term_features["l_stem"] = (
                        total_length - term_features["r_stem"] -
                        term_features["loop"])
                    if (total_length <= term_len) and (
                            term_features["loop"] <= args_term.max_loop) and (
                            term_features["loop"] >= args_term.min_loop) and (
                            ((term_features["r_stem"] +
                              term_features["l_stem"] -
                              term_features["real_miss"]) / 2) >=
                            args_term.min_stem) and (
                            ((term_features["r_stem"] +
                              term_features["l_stem"] -
                              term_features["real_miss"]) / 2) <=
                            args_term.max_stem):
                        if strand == "+":
                            import_candidate(
                                cands, term_features, strain,
                                start + (nts - term_features["st_pos"]) - 10,
                                start + nts - 1 + 10, ut, name, total_length,
                                strand, parent_p, parent_m, p_pos, m_pos)
                        else:
                            import_candidate(
                                cands, term_features, strain,
                                end - (nts - 1) - 10,
                                end - (nts - term_features["st_pos"]) + 10,
                                ut, name, total_length, strand,
                                parent_p, parent_m, p_pos, m_pos)
    return cands


def check_parent(genes, term, detects, strand, fuzzy_up, fuzzy_down, type_):
    tmp = None
    for gene in genes:
        if (term["strain"] == gene.seq_id) and (
                gene.strand == strand):
            if type_ == "parent_p":
                if ((term["start"] - fuzzy_down) <= gene.end) and (
                        term["start"] >= gene.end):
                    detects[type_] = True
                    tmp = get_feature(gene)
                elif ((gene.end - term["end"]) <= fuzzy_up) and (
                        (gene.end - term["end"]) >= 0):
                    detects[type_] = True
                    tmp = get_feature(gene)
                elif ((gene.end - term["start"]) > fuzzy_up) and (
                        (gene.end - term["start"]) >= 0):
                    break
            elif type_ == "parent_m":
                if ((term["end"] + fuzzy_down) >= gene.start) and (
                        term["end"] <= gene.start):
                    detects[type_] = True
                    tmp = get_feature(gene)
                elif ((term["start"] - gene.start) <= fuzzy_up) and (
                        (term["start"] - gene.start) >= 0):
                    detects[type_] = True
                    tmp = get_feature(gene)
                elif (gene.start - term["end"] > fuzzy_down):
                    break
    return tmp


def parents(terms, genes, args_term):
    '''assign the associated gene to terminator'''
    for term in terms:
        detects = {"parent_p": False, "parent_m": False}
        if "tran" in term["parent_p"]:
            tmp_p = check_parent(genes, term, detects, "+",
                                 args_term.fuzzy_up_gene,
                                 args_term.fuzzy_down_gene, "parent_p")
#            pos = term["parent_p"].split(":")[-1].split("_")[0].split("-")[-1]
            pos = term["p_pos"].split("-")[-1]
            if ((term["start"] - int(pos) <= args_term.fuzzy_down_ta) and (
                 term["start"] - int(pos) >= 0)) or (
                (int(pos) - term["end"] <= args_term.fuzzy_up_ta) and (
                 int(pos) - term["end"] >= 0)):
                pass
            else:
                term["parent_p"] = ""
        if "tran" in term["parent_m"]:
            tmp_m = check_parent(genes, term, detects, "-",
                                 args_term.fuzzy_up_gene,
                                 args_term.fuzzy_down_gene, "parent_m")
            pos = term["m_pos"].split("-")[0]
#            pos = term["parent_m"].split(":")[-1].split("_")[0].split("-")[0]
            if ((int(pos) - term["end"] <= args_term.fuzzy_down_ta) and (
                 int(pos) - term["end"] >= 0)) or (
                (term["start"] - int(pos) <= args_term.fuzzy_up_ta) and (
                 term["start"] - int(pos) >= 0)):
                pass
            else:
                term["parent_m"] = ""
        if detects["parent_p"]:
            if len(term["parent_p"]) == 0:
                term["parent_p"] = tmp_p
            else:
                term["parent_p"] = ",".join([term["parent_p"], tmp_p])
        if detects["parent_m"]:
            if len(term["parent_m"]) == 0:
                term["parent_m"] = tmp_m
            else:
                term["parent_m"] = ",".join([term["parent_m"], tmp_m])


def read_gff(seq_file, gff_file, tran_file):
    genome = {}
    genes = []
    trans = []
    for entry in Gff3Parser().entries(open(gff_file)):
        if (entry.feature == "gene"):
            genes.append(entry)
    for entry in Gff3Parser().entries(open(tran_file)):
        trans.append(entry)
    with open(seq_file, "r") as q_h:
        for line in q_h:
            line = line.strip()
            if line.startswith(">"):
                strain = line[1:]
                genome[strain] = ""
            else:
                genome[strain] = genome[strain] + line
    genes = sorted(genes, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    trans = sorted(trans, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    return genes, genome, trans


def compare_anno(gffs, cands, fuzzy_up, fuzzy_down):
    '''compare the terminator with CDS'''
    detect = False
    new_cands = []
    for cand in cands:
        for gff in gffs:
            if (gff.seq_id == cand["strain"]) and (
                    gff.strand == cand["strand"]):
                if cand["strand"] == "+":
                    if (gff.start <= cand["start"]) and (
                            gff.end >= cand["start"]) and (
                            gff.end <= cand["end"]):
                        detect = True
                        break
                    elif (math.fabs(gff.end - cand["end"]) <= fuzzy_up) and (
                            gff.end >= cand["end"]):
                        detect = True
                        break
                    elif (math.fabs(gff.end - cand["start"]) <=
                            fuzzy_down) and (gff.end <= cand["start"]):
                        detect = True
                        break
                else:
                    if (gff.start >= cand["start"]) and (
                            gff.start <= cand["end"]) and (
                            gff.end >= cand["end"]):
                        detect = True
                        break
                    elif (math.fabs(gff.start - cand["end"]) <=
                            fuzzy_down) and (gff.start >= cand["end"]):
                        detect = True
                        break
                    elif (math.fabs(gff.start - cand["start"]) <=
                            fuzzy_up) and (cand["start"] >= gff.start):
                        detect = True
                        break
        if detect:
            detect = False
            new_cands.append(cand)
    return new_cands


def merge_cands(new_cands_gene, new_cands_ta):
    new_cands = []
    for cand_gene in new_cands_gene:
        new_cands.append(cand_gene)
    for cand_ta in new_cands_ta:
        if cand_ta not in new_cands:
            new_cands.append(cand_ta)
    new_cands = sorted(new_cands, key=lambda k: (k["strain"], k["start"],
                                                 k["end"], k["strand"]))
    return new_cands


def get_seq_sec(s_h, sec_seq):
    '''extract the secondary structure information'''
    for line in s_h:
        if ("(" in line) or ("." in line) or (")" in line):
            line = line.split(" ")
            sec_seq["sec"] = line[0]
            break
        else:
            sec_seq["seq"] = line


def poly_t(seq_file, sec_file, gff_file, tran_file, out_file, args_term):
    '''detect the sec str of terminator'''
    terms = []
    genes, genome, trans = read_gff(seq_file, gff_file, tran_file)
    out = open(out_file, "w")
    with open(sec_file, "r") as s_h:
        for line in s_h:
            line = line.strip()
            if line.startswith(">"):
                line = line[1:]
                name = line.split("|")[0]
                start = int(line.split("|")[1])
                end = int(line.split("|")[2])
                strain = line.split("|")[3]
                parent_p = line.split("|")[4]
                parent_m = line.split("|")[5]
                p_pos = line.split("|")[6]
                m_pos = line.split("|")[7]
                strand = line.split("|")[-1]
                sec_seq = {"sec": "", "seq": ""}
                get_seq_sec(s_h, sec_seq)
                if len(sec_seq["seq"]) <= 6:
                    continue
                else:
                    cands = detect_candidates(
                        sec_seq["seq"], sec_seq["sec"], name, strain, start,
                        end, parent_p, parent_m, strand, args_term,
                        p_pos, m_pos)
                cands = sorted(cands, key=lambda x: (x["miss"], x["start"]))
                new_cands_gene = compare_anno(genes, cands,
                                              args_term.fuzzy_up_gene,
                                              args_term.fuzzy_down_gene)
                new_cands_ta = compare_anno(trans, cands,
                                            args_term.fuzzy_up_ta,
                                            args_term.fuzzy_down_ta)
                new_cands = merge_cands(new_cands_gene, new_cands_ta)
                filter_term(new_cands, terms, args_term.miss_rate)
    parents(terms, genes, args_term)
    for term in terms:
        print_ = False
        if (term["strand"] == "+") and (len(term["parent_p"]) != 0):
            print_ = True
        elif (term["strand"] == "-") and (len(term["parent_m"]) != 0):
            print_ = True
        if print_:
            out.write("\t".join([term["strain"], str(term["start"]),
                      str(term["end"]), term["name"], str(term["miss"]),
                      str(term["loop"]), str(term["length"]),
                      str(term["r_stem"]), term["strand"], str(term["l_stem"]),
                      term["parent_p"], term["parent_m"],
                      str(term["ut"])]) + "\n")
