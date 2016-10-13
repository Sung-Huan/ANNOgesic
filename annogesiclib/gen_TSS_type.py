import math
import copy
from annogesiclib.helper import Helper


def remove_primary(tss, tss_entry):
    final_types = []
    final_utrs = []
    final_genes = []
    tss_dict = tss_entry[1]
    types = tss_dict["type"].split("&")
    utrs = tss_dict["UTR_length"].split("&")
    genes = tss_dict["associated_gene"].split("&")
    index = 0
    for type_ in types:
        if type_ != "Primary":
            final_types.append(type_)
            final_utrs.append(utrs[index])
            final_genes.append(genes[index])
        index += 1
    strand = Helper().get_strand_name(tss.strand)
    tss_dict = {"Name": "_".join(["TSS:" + str(tss.start), strand]),
                "type": "&".join(final_types),
                "UTR_length": "&".join(final_utrs),
                "associated_gene": "&".join(final_genes)}
    tss_string = ";".join(["=".join(["UTR_length", tss_dict["UTR_length"]]),
                           "=".join(["associated_gene",
                                     tss_dict["associated_gene"]]),
                           "=".join(["type", tss_dict["type"]]),
                           "=".join(["Name", tss_dict["Name"]])])
    return (tss_string, tss_dict)


def import_to_tss(tss_type, cds_pos, tss, locus_tag, tss_entry):
    if cds_pos == "NA":
        utr = "_".join([tss_type, "NA"])
    else:
        utr = "_".join([tss_type, str(int(math.fabs(cds_pos - tss.start)))])
    if len(tss_entry) != 0:
        tss_dict = tss_entry[1]
        tss_dict_types = tss_dict["type"].split("&")
        tss_dict_utrs = tss_dict["UTR_length"].split("&")
        tss_dict_tags = tss_dict["associated_gene"].split("&")
        if tss_type == "Primary" and ("Primary" in tss_dict["type"]):
            index = 0
            for tss_dict_type in tss_dict_types:
                if "Primary" in tss_dict_type:
                    utr_length = tss_dict_utrs[index].split("_")
                    if math.fabs(cds_pos - tss.start) < int(utr_length[1]):
                        tss_dict_utrs[index] = utr
                        tss_dict_tags[index] = locus_tag
                index += 1
        else:
            tss_dict_types.append(tss_type)
            tss_dict_utrs.append(utr)
            tss_dict_tags.append(locus_tag)
        strand = Helper().get_strand_name(tss.strand)
        tss_dict = {"Name": "_".join(["TSS:" + str(tss.start), strand]),
                    "type": "&".join(tss_dict_types),
                    "UTR_length": "&".join(tss_dict_utrs),
                    "associated_gene": "&".join(tss_dict_tags)}
    else:
        strand = Helper().get_strand_name(tss.strand)
        tss_dict = {"Name": "_".join(["TSS:" + str(tss.start), strand]),
                    "type": tss_type,
                    "UTR_length": utr,
                    "associated_gene": locus_tag}
    tss_string = ";".join(["=".join(["UTR_length", tss_dict["UTR_length"]]),
                           "=".join(["associated_gene",
                                     tss_dict["associated_gene"]]),
                           "=".join(["type", tss_dict["type"]]),
                           "=".join(["Name", tss_dict["Name"]])])
    return (tss_string, tss_dict)


def is_primary(cds_start, cds_end, tss_pos, strand):
    if strand == "+":
        if (is_utr(cds_start, tss_pos, 300) and (cds_start >= tss_pos)):
            return True
    else:
        if (is_utr(tss_pos, cds_end, 300) and (cds_end <= tss_pos)):
            return True


def is_internal(cds_start, cds_end, tss_pos, strand):
    if ((cds_start < tss_pos) and (cds_end > tss_pos)) or (
            (strand == "+") and (tss_pos == cds_end)) or (
            (strand == "-") and (tss_pos == cds_start)):
        return True


def is_antisense(cds_start, cds_end, tss_pos, strand):
    if ((is_utr(cds_start, tss_pos, 100)) and (cds_start >= tss_pos)) or (
            (is_utr(tss_pos, cds_end, 100)) and (cds_end <= tss_pos)) or (
             is_internal(cds_start, cds_end, tss_pos, strand)):
        return True


def is_utr(pos1, pos2, length):
    if (pos1 - pos2 <= length):
        return True


def same_strand_tss_gene(gene, tss, anti_ends, gene_ends, checks, tss_entry):
    '''deal with the the gene and TSS which located at the same strands'''
    if is_primary(gene.start, gene.end, tss.start, tss.strand):
        ori_entry = copy.deepcopy(tss_entry)
        if "locus_tag" in gene.attributes.keys():
            locus_tag = gene.attributes["locus_tag"]
        else:
            locus_tag = "".join([gene.feature, ":", str(gene.start), "-",
                                 str(gene.end), "_", gene.strand])
        if tss.strand == "+":
            if ((anti_ends["reverse"] != -1) and (
                    anti_ends["reverse"] - gene.start) > 0) or (
                    anti_ends["reverse"] == -1):
                tss_entry = import_to_tss("Primary", gene.start, tss,
                                          locus_tag, ori_entry)
                checks["orphan"] = False
                gene_ends["forward"] = gene.start
            elif (anti_ends["reverse"] != -1) and (
                    (anti_ends["reverse"] - gene.start) < 0):
                if (checks["int_anti"]) or (
                        (tss.start - anti_ends["reverse"]) > 0):
                    tss_entry = import_to_tss("Primary", gene.start, tss,
                                              locus_tag, ori_entry)
                    checks["orphan"] = False
                    gene_ends["forward"] = gene.start
        else:
            if ((anti_ends["forward"] != -1) and (
                    gene.end - anti_ends["forward"]) > 0) or (
                    anti_ends["forward"] == -1):
                tss_entry = import_to_tss("Primary", gene.end, tss,
                                          locus_tag, ori_entry)
                checks["orphan"] = False
                gene_ends["reverse"] = gene.end
    if is_internal(gene.start, gene.end, tss.start, tss.strand):
        ori_entry = copy.deepcopy(tss_entry)
        if "locus_tag" in gene.attributes.keys():
            locus_tag = gene.attributes["locus_tag"]
        else:
            locus_tag = "".join([gene.feature, ":", str(gene.start), "-",
                                 str(gene.end), "_", gene.strand])
        tss_entry = import_to_tss("Internal", "NA", tss, locus_tag, ori_entry)
        checks["orphan"] = False
    return tss_entry


def diff_strand_tss_gene(gene, tss, anti_ends, gene_ends, checks, tss_entry):
    '''deal with the the gene and TSS which located at different strands'''
    if is_antisense(gene.start, gene.end, tss.start, tss.strand):
        checks["int_anti"] = False
        if tss.strand == "-":
            anti_ends["forward"] = gene.start
            if (gene_ends["reverse"] != -1) and (
                    gene.start - gene_ends["reverse"]) > 0:
                if is_internal(gene.start, gene.end, tss.start, tss.strand):
                    pass
        else:
            anti_ends["reverse"] = gene.end
            if is_internal(gene.start, gene.end, tss.start, tss.strand):
                checks["int_anti"] = True
        if "locus_tag" in gene.attributes.keys():
            locus_tag = gene.attributes["locus_tag"]
        else:
            locus_tag = "".join([gene.feature, ":", str(gene.start), "-",
                                 str(gene.end), "_", gene.strand])
        ori_entry = copy.deepcopy(tss_entry)
        tss_entry = import_to_tss("Antisense", "NA", tss, locus_tag, ori_entry)
        checks["orphan"] = False
    return tss_entry


def compare_tss_cds(tss, cdss, genes):
    '''compare TSS and CDS to classify the TSSs'''
    tss_entry = []
    gene_ends = {"forward": -1, "reverse": -1}
    anti_ends = {"forward": -1, "reverse": -1}
    checks = {"orphan": True, "int_anti": None}
    if (len(genes) == 0):
        datas = copy.deepcopy(cdss)
    else:
        datas = copy.deepcopy(genes)
    for data in datas:
        ori_entry = copy.deepcopy(tss_entry)
        if data.strand == tss.strand:
            tss_entry = same_strand_tss_gene(data, tss, anti_ends,
                                             gene_ends, checks, ori_entry)
        else:
            tss_entry = diff_strand_tss_gene(data, tss, anti_ends,
                                             gene_ends, checks, ori_entry)
    if checks["orphan"]:
        ori_entry = copy.deepcopy(tss_entry)
        tss_entry = import_to_tss("Orphan", "NA", tss, "NA", ori_entry)
    return tss_entry


def fix_attributes(tss, tss_entry):
    '''change the primary TSS to secondary TSS'''
    index = 0
    genes = tss.attributes["associated_gene"].split("&")
    utrs = tss.attributes["UTR_length"].split("&")
    types = tss.attributes["type"].split("&")
    for gene in genes:
        if gene == tss_entry["locus"]:
            utrs[index] = utrs[index].replace("Primary", "Secondary")
            types[index] = types[index].replace("Primary", "Secondary")
        index += 1
    tss.attributes["UTR_length"] = "&".join(utrs)
    tss.attributes["type"] = "&".join(types)


def detect_coverage(wigs, tss, ref):
    tss_cover = -1
    ref_cover = -1
    for strain, conds in wigs.items():
        if strain == tss.seq_id:
            tss_cover = 0
            ref_cover = 0
            for cond, tracks in conds.items():
                for lib_name, covers in tracks.items():
                    if ((tss.start + 1) <= len(covers)) and (
                            (ref.start + 1) <= len(covers)):
                        if tss.strand == "+":
                            diff_t = (covers[tss.start - 1] -
                                      covers[tss.start - 2])
                            diff_r = (covers[ref.start - 1] -
                                      covers[ref.start - 2])
                        else:
                            diff_t = (covers[tss.start - 1] -
                                      covers[tss.start])
                            diff_r = (covers[ref.start - 1] -
                                      covers[ref.start])
                        tss_cover = tss_cover + diff_t
                        ref_cover = ref_cover + diff_r
    return (tss_cover, ref_cover)


def del_repeat(tsss):
    '''delete the repeat TSSs'''
    for tss in tsss:
        types = tss.attributes["type"].split("&")
        utrs = tss.attributes["UTR_length"].split("&")
        genes = tss.attributes["associated_gene"].split("&")
        detect_pri = False
        detect_sec = False
        index = 0
        final_types = []
        final_utrs = []
        final_genes = []
        for type_ in types:
            if (type_ == "Primary") and (not detect_pri):
                detect_pri = True
                pri_utr = int(utrs[index].split("_")[1])
                real_index = index
            elif (type_ == "Primary") and (detect_pri):
                compare_utr = int(utrs[index].split("_")[1])
                if compare_utr < pri_utr:
                    pri_utr = compare_utr
                    real_index = index
            elif (type_ == "Secondary") and (not detect_sec):
                detect_sec = True
                sec_utr = int(utrs[index].split("_")[1])
                real_index2 = index
            elif (type_ == "Secondary") and (detect_sec):
                compare_utr = int(utrs[index].split("_")[1])
                if compare_utr < sec_utr:
                    sec_utr = compare_utr
                    real_index2 = index
            elif (type_ == "Antisense") or \
                 (type_ == "Internal") or \
                 (type_ == "Orphan"):
                final_types.append(types[index])
                final_utrs.append(utrs[index])
                final_genes.append(genes[index])
            index += 1
        if detect_pri:
            final_types.append(types[real_index])
            final_utrs.append(utrs[real_index])
            final_genes.append(genes[real_index])
        else:
            if detect_sec:
                final_types.append(types[real_index2])
                final_utrs.append(utrs[real_index2])
                final_genes.append(genes[real_index2])
        tss.attributes["type"] = "&".join(final_types)
        tss.attributes["UTR_length"] = "&".join(final_utrs)
        tss.attributes["associated_gene"] = "&".join(final_genes)


def get_primary_locus_tag(tss):
    tsss = []
    tss_types = tss.attributes["type"].split("&")
    tss_locus_tags = tss.attributes["associated_gene"].split("&")
    tss_utr_lengths = tss.attributes["UTR_length"].split("&")
    index = 0
    for tss_type in tss_types:
        if "Primary" in tss_type:
            tsss.append({"locus": tss_locus_tags[index],
                         "utr": int(tss_utr_lengths[index].split("_")[1]),
                         "type": tss_type})
        index += 1
    return tsss


def fix_primary_type(tsss, wigs_f, wigs_r):
    '''If one gene is associated with multiple primary TSSs, 
    it will assing the low expressed one to be secondary TSS'''
    for tss in tsss:
        if ("Primary" in tss.attributes["type"]):
            tss_entrys = get_primary_locus_tag(tss)
            for ref in tsss:
                if (ref.seq_id == tss.seq_id) and (
                        ref.strand == tss.strand) and (
                        ref.start == tss.start):
                    pass
                else:
                    if ("Primary" in ref.attributes["type"]):
                        ref_entrys = get_primary_locus_tag(ref)
                        for tss_entry in tss_entrys:
                            for ref_entry in ref_entrys:
                                if (tss_entry["locus"] ==
                                    ref_entry["locus"]) and (
                                        tss_entry["type"] == "Primary") and (
                                        ref_entry["type"] == "Primary") and (
                                        tss.seq_id == ref.seq_id):
                                    if tss.strand == "+":
                                        covers = detect_coverage(
                                            wigs_f, tss, ref)
                                    else:
                                        covers = detect_coverage(
                                            wigs_r, tss, ref)
                                    tss_cover = covers[0]
                                    ref_cover = covers[1]
                                    if tss_cover < ref_cover:
                                        fix_attributes(tss, tss_entry)
                                    elif tss_cover > ref_cover:
                                        fix_attributes(ref, ref_entry)
                                    elif tss_cover == ref_cover:
                                        if (tss_entry["utr"] <
                                                ref_entry["utr"]):
                                            fix_attributes(ref, ref_entry)
                                        elif (tss_entry["utr"] >
                                              ref_entry["utr"]):
                                            fix_attributes(tss, tss_entry)
    del_repeat(tsss)
    return tsss
