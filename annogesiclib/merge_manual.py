import os
import math
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.parser_wig import WigParser
from annogesiclib.helper import Helper


def get_primary_locus_tag(tss):
    tsss = []
    tss_types = tss.attributes["type"].split(",")
    tss_locus_tags = tss.attributes["associated_gene"].split(",")
    tss_utr_lengths = tss.attributes["utr_length"].split(",")
    index = 0
    for tss_type in tss_types:
        if "Primary" in tss_type:
            tsss.append({"locus": tss_locus_tags[index],
                         "utr": int(tss_utr_lengths[index].split("_")[1]),
                         "type": tss_type})
        index += 1
    return tsss


def detect_coverage(wigs, tss, ref):
    for strain, tracks in wigs.items():
        if strain == tss.seq_id:
            tss_cover = 0
            ref_cover = 0
            for wig in tracks.values():
                if ((tss.start + 1) <= len(wig)) and (
                        (ref.start + 1) <= len(wig)):
                    if tss.strand == "+":
                        diff_t = (wig[tss.start - 1]["coverage"] -
                                  wig[tss.start - 2]["coverage"])
                        diff_r = (wig[ref.start - 1]["coverage"] -
                                  wig[ref.start - 2]["coverage"])
                    else:
                        diff_t = (wig[tss.start - 1]["coverage"] -
                                  wig[tss.start]["coverage"])
                        diff_r = (wig[ref.start - 1]["coverage"] -
                                  wig[ref.start]["coverage"])
                    tss_cover = tss_cover + diff_t
                    ref_cover = ref_cover + diff_r
    return tss_cover, ref_cover


def fix_attributes(tss, tss_entry):
    '''change the primary TSS to secondary TSS'''
    index = 0
    genes = tss.attributes["associated_gene"].split(",")
    utrs = tss.attributes["utr_length"].split(",")
    types = tss.attributes["type"].split(",")
    for gene in genes:
        if gene == tss_entry["locus"]:
            utrs[index] = utrs[index].replace("Primary", "Secondary")
            types[index] = types[index].replace("Primary", "Secondary")
        index += 1
    tss.attributes["utr_length"] = ",".join(utrs)
    tss.attributes["type"] = ",".join(types)


def del_repeat(tsss):
    '''delete the repeat TSS'''
    for tss in tsss:
        types = tss.attributes["type"].split(",")
        utrs = tss.attributes["utr_length"].split(",")
        genes = tss.attributes["associated_gene"].split(",")
        detect = {"pri": False, "sec": False}
        index = 0
        finals = {"types": [], "utrs": [], "genes": []}
        for type_ in types:
            if (type_ == "Primary") and (detect["pri"] is False):
                detect["pri"] = True
                pri_utr = int(utrs[index].split("_")[1])
                real_index = index
            elif (type_ == "Primary") and (detect["pri"] is True):
                compare_utr = int(utrs[index].split("_")[1])
                if compare_utr < pri_utr:
                    pri_utr = compare_utr
                    real_index = index
            elif (type_ == "Secondary") and (detect["sec"] is False):
                detect["sec"] = True
                sec_utr = int(utrs[index].split("_")[1])
                real_index2 = index
            elif (type_ == "Secondary") and (detect["sec"] is True):
                compare_utr = int(utrs[index].split("_")[1])
                if compare_utr < sec_utr:
                    sec_utr = compare_utr
                    real_index2 = index
            elif (type_ == "Antisense") or \
                 (type_ == "Internal") or \
                 (type_ == "Orphan"):
                finals["types"].append(types[index])
                finals["utrs"].append(utrs[index])
                finals["genes"].append(genes[index])
            index += 1
        if detect["pri"] is True:
            finals["types"].append(types[real_index])
            finals["utrs"].append(utrs[real_index])
            finals["genes"].append(genes[real_index])
        else:
            if detect["sec"] is True:
                finals["types"].append(types[real_index2])
                finals["utrs"].append(utrs[real_index2])
                finals["genes"].append(genes[real_index2])
        tss.attributes["type"] = ",".join(finals["types"])
        tss.attributes["utr_length"] = ",".join(finals["utrs"])
        tss.attributes["associated_gene"] = ",".join(finals["genes"])


def fix_primary_type(tsss, wigs_f, wigs_r):
    '''if one gene is associated with multiple TSSs, it will 
    check the coverage and assign the low expressed TSS to be 
    secondary TSS'''
    for tss in tsss:
        if ("Primary" in tss.attributes["type"]):
            tss_entrys = get_primary_locus_tag(tss)
            for ref in tsss:
                if (ref.seq_id == tss.seq_id) and \
                   (ref.strand == tss.strand) and \
                   (ref.start == tss.start):
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
                                        tss_cover, ref_cover = detect_coverage(
                                            wigs_f, tss, ref)
                                    else:
                                        tss_cover, ref_cover = detect_coverage(
                                            wigs_r, tss, ref)
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


def define_attributes(tss):
    string = []
    for key, value in tss.attributes.items():
        if key != "print":
            if key != "ID":
                string.append("=".join([key, value]))
            elif key == "Name":
                string.append("=".join([key, str(tss.start) + tss.strand]))
    return ";".join(string)


def remove_primary(tss, tss_entry):
    final_types = []
    final_utrs = []
    final_genes = []
    tss_dict = tss_entry[1]
    types = tss_dict["type"].split(",")
    utrs = tss_dict["utr_length"].split(",")
    genes = tss_dict["associated_gene"].split(",")
    index = 0
    for type_ in types:
        if type_ != "Primary":
            final_types.append(type_)
            final_utrs.append(utrs[index])
            final_genes.append(genes[index])
        index += 1
    tss_dict = {"Name": "TSS_" + str(tss.start) + tss.strand,
                "type": ",".join(final_types),
                "utr_length": ",".join(final_utrs),
                "associated_gene": ",".join(final_genes)}
    tss_string = ";".join(["=".join(["utr_length", tss_dict["utr_length"]]),
                           "=".join(["associated_gene",
                                     tss_dict["associated_gene"]]),
                           "=".join(["type", tss_dict["type"]]),
                           "=".join(["Name", tss_dict["Name"]])])
    return [tss_string, tss_dict]


def import_to_tss(tss_type, cds_pos, tss, locus_tag, tss_entry):
    if cds_pos == "NA":
        utr = "_".join([tss_type, "NA"])
    else:
        utr = "_".join([tss_type, str(int(math.fabs(cds_pos - tss.start)))])
    if len(tss_entry) != 0:
        tss_dict = tss_entry[1]
        tss_dict_types = tss_dict["type"].split(",")
        tss_dict_utrs = tss_dict["utr_length"].split(",")
        tss_dict_tags = tss_dict["associated_gene"].split(",")
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
        tss_dict = {"Name": "TSS_" + str(tss.start) + tss.strand,
                    "type": ",".join(tss_dict_types),
                    "utr_length": ",".join(tss_dict_utrs),
                    "associated_gene": ",".join(tss_dict_tags)}
    else:
        tss_dict = {"Name": "TSS_" + str(tss.start) + tss.strand,
                    "type": tss_type,
                    "utr_length": utr,
                    "associated_gene": locus_tag}
    tss_string = ";".join(["=".join(["utr_length", tss_dict["utr_length"]]),
                           "=".join(["associated_gene",
                                     tss_dict["associated_gene"]]),
                           "=".join(["type", tss_dict["type"]]),
                           "=".join(["Name", tss_dict["Name"]])])
    return (tss_string, tss_dict)


def same_strand_tss_gene(gene, tss, anti_ends, gene_ends, checks, tss_entry):
    '''check the TSS and gene which are at the same strand'''
    if is_primary(gene.start, gene.end, tss.start, tss.strand):
        locus_tag = gene.attributes["locus_tag"]
        if tss.strand == "+":
            if ((anti_ends["reverse"] != -1) and (
                    anti_ends["reverse"] - gene.start) > 0) or (
                    anti_ends["reverse"] == -1):
                tss_entry = import_to_tss("Primary", gene.start, tss,
                                          locus_tag, tss_entry)
                checks["orphan"] = False
                gene_ends["forward"] = gene.start
            elif (anti_ends["reverse"] != -1) and (
                    (anti_ends["reverse"] - gene.start) < 0):
                if (checks["int_anti"] is True) or (
                        (tss.start - anti_ends["reverse"]) > 0):
                    tss_entry = import_to_tss("Primary", gene.start, tss,
                                              locus_tag, tss_entry)
                    checks["orphan"] = False
                    gene_ends["forward"] = gene.start
        else:
            if ((anti_ends["forward"] != -1) and (
                    gene.end - anti_ends["forward"]) > 0) or (
                    anti_ends["forward"] == -1):
                tss_entry = import_to_tss("Primary", gene.end, tss,
                                          locus_tag, tss_entry)
                checks["orphan"] = False
                gene_ends["reverse"] = gene.end
    if is_internal(gene.start, gene.end, tss.start, tss.strand):
        locus_tag = gene.attributes["locus_tag"]
        tss_entry = import_to_tss("Internal", "NA", tss, locus_tag, tss_entry)
        checks["orphan"] = False
    return tss_entry


def diff_strand_tss_gene(gene, tss, anti_ends, gene_ends, checks, tss_entry):
    '''check the TSS and gene which are at the same strand'''
    if is_antisense(gene.start, gene.end, tss.start, tss.strand):
        checks["int_anti"] = False
        if tss.strand == "-":
            anti_ends["forward"] = gene.start
            if (gene_ends["reverse"] != -1) and (
                    (gene.start - gene_ends["reverse"]) > 0):
                if is_internal(gene.start, gene.end, tss.start, tss.strand):
                    pass
                else:
                    if (tss.start - gene.end) > 0:
                        tss_entry = remove_primary(tss, tss_entry)
        else:
            anti_ends["reverse"] = gene.end
            if is_internal(gene.start, gene.end, tss.start, tss.strand):
                checks["int_anti"] = True
            if (gene_ends["forward"] != -1) and (
                    (gene.start - gene_ends["forward"]) > 0):
                if (gene.start - tss.start) > 0:
                    tss_entry = remove_primary(tss, tss_entry)
        locus_tag = gene.attributes["locus_tag"]
        tss_entry = import_to_tss("Antisense", "NA", tss, locus_tag, tss_entry)
        checks["orphan"] = False
    return tss_entry


def compare_tss_gene(tss, genes):
    '''compare TSS and gene to find the relation'''
    tss_entry = []
    gene_ends = {"forward": -1, "reverse": -1}
    anti_ends = {"forward": -1, "reverse": -1}
    checks = {"orphan": True, "int_anti": None}
    for gene in genes:
        if gene.strand == tss.strand:
            tss_entry = same_strand_tss_gene(gene, tss, anti_ends,
                                             gene_ends, checks, tss_entry)
        else:
            tss_entry = diff_strand_tss_gene(gene, tss, anti_ends,
                                             gene_ends, checks, tss_entry)
    if checks["orphan"]:
        tss_entry = import_to_tss("Orphan", "NA", tss, "NA", tss_entry)
    return tss_entry


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


def print_all_unique(out, overlap_num, nums):
    if ((nums["tss_p"] != 0) or (overlap_num != 0)) and (
            (nums["tss_m"] != 0) or (overlap_num != 0)):
        out.write("the number of overlap between "
                  "TSSpredator and manual = {0} ".format(
                      overlap_num))
        out.write("(overlap of all TSSpredator = {0}, ".format(
                  float(overlap_num) /
                  (float(nums["tss_p"]) + float(overlap_num))))
        out.write("overlap of all manual = {0})\n".format(
                  float(overlap_num) /
                  (float(nums["tss_m"]) + float(overlap_num))))
        out.write("the number of unique in TSSpredator = {0} ({1})\n".format(
                  nums["tss_p"], float(nums["tss_p"]) /
                  (float(nums["tss_p"]) + float(overlap_num))))
        out.write("the number of unique in manual = {0} ({1})\n".format(
                  nums["tss_m"], float(nums["tss_m"]) /
                  (float(nums["tss_m"]) + float(overlap_num))))
    else:
        out.write("No TSS candidates which be predicted by TSSpredator.")


def print_stat(num_strain, stat_file, overlap_num, nums):
    if len(num_strain) != 0:
        out = open(stat_file, "w")
        if len(num_strain.keys()) == 1:
            print_all_unique(out, overlap_num, nums)
        else:
            out.write("All genomes: \n")
            print_all_unique(out, overlap_num, nums)
            for strain in num_strain.keys():
                if (num_strain[strain]["tsspredator"] == 0) and \
                   (num_strain[strain]["overlap"] == 0):
                    perc_tsspredator = "NA"
                    perc_tsspredator_uni = "NA"
                else:
                    perc_tsspredator = str(
                        float(num_strain[strain]["overlap"]) / (
                            float(num_strain[strain]["tsspredator"]) +
                            float(num_strain[strain]["overlap"])))
                    perc_tsspredator_uni = str(
                        float(num_strain[strain]["tsspredator"]) / (
                            float(num_strain[strain]["tsspredator"]) +
                            float(num_strain[strain]["overlap"])))
                if (num_strain[strain]["manual"] == 0) and (
                        num_strain[strain]["overlap"] == 0):
                    perc_manual = "NA"
                    perc_manual_uni = "NA"
                else:
                    perc_manual = str(
                        float(num_strain[strain]["overlap"]) / (
                            float(num_strain[strain]["manual"]) +
                            float(num_strain[strain]["overlap"])))
                    perc_manual_uni = str(
                        float(num_strain[strain]["manual"]) / (
                            float(num_strain[strain]["manual"]) +
                            float(num_strain[strain]["overlap"])))
                out.write(strain + ": \n")
                out.write("the number of overlap between "
                          "TSSpredator and manual = {0} ".format(
                              num_strain[strain]["overlap"]))
                out.write("(overlap of all TSSpredator = {0}, ".format(
                          perc_tsspredator))
                out.write("overlap of all manual = {0})\n".format(perc_manual))
                out.write("the number of unique in "
                          "TSSpredator = {0} ({1})\n".format(
                              num_strain[strain]["tsspredator"],
                              perc_tsspredator_uni))
                out.write("the number of unique in manual = {0} ({1})\n".format(
                          num_strain[strain]["manual"], perc_manual_uni))


def read_wig(filename, strand):
    wigs = {}
    wig_parser = WigParser()
    if filename:
        wig_fh = open(filename)
        for entry in wig_parser.parser(wig_fh, strand):
            if entry.strain not in wigs.keys():
                strain = entry.strain
                wigs[strain] = {}
            if entry.track not in wigs[strain].keys():
                wigs[strain][entry.track] = []
            wigs[strain][entry.track].append({
                 "pos": entry.pos, "coverage": entry.coverage,
                 "strand": entry.strand})
        wig_fh.close()
    return wigs


def read_gff(tss_predict_file, tss_manual_file, gff_file, lengths):
    tsss = {"tsss_p": [], "tsss_m": [], "merge": []}
    cdss = []
    genes = []
    gff_parser = Gff3Parser()
    tssp_fh = open(tss_predict_file, "r")
    tssm_fh = open(tss_manual_file, "r")
    g_f = open(gff_file, "r")
    for entry in gff_parser.entries(tssp_fh):
        entry.attributes["print"] = False
        tsss["tsss_p"].append(entry)
    tssp_fh.close()
    tsss["tsss_p"] = sorted(tsss["tsss_p"], key=lambda k: (k.seq_id, k.start,
                                                           k.end, k.strand))
    for entry in gff_parser.entries(tssm_fh):
        if (entry.seq_id in lengths.keys()) or ("all" in lengths.keys()):
            entry.attributes["print"] = False
            entry.attributes["libs"] = "manual"
            entry.attributes["method"] = "manual"
            tsss["tsss_m"].append(entry)
    tssm_fh.close()
    tsss["tsss_m"] = sorted(tsss["tsss_m"], key=lambda k: (k.seq_id, k.start,
                                                           k.end, k.strand))
    for entry in gff_parser.entries(g_f):
        if (Helper().feature_without_notgene(entry)):
            cdss.append(entry)
        if entry.feature == "gene":
            genes.append(entry)
    g_f.close()
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    genes = sorted(genes, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    return tsss, cdss, genes


def merge_libs(input_libs, wig_folder, program):
    if "merge_forward.wig" in os.listdir(os.getcwd()):
        os.remove("merge_forward.wig")
    if "merge_reverse.wig" in os.listdir(os.getcwd()):
        os.remove("merge_reverse.wig")
    if program == "TSS":
        type_ = "tex"
    elif program == "processing":
        type_ = "notex"
    for lib in input_libs:
        datas = lib.split(":")
        if (datas[1] == type_) and (datas[4] == "+"):
            Helper().merge_file(os.path.join(wig_folder, datas[0]),
                                os.path.join(os.getcwd(), "merge_forward.wig"))
        elif (datas[1] == type_) and (datas[4] == "-"):
            Helper().merge_file(os.path.join(wig_folder, datas[0]),
                                os.path.join(os.getcwd(), "merge_reverse.wig"))


def check_overlap(overlap, pre_tss, nums, length, num_strain, overlap_num,
                  tss_m, tss_p, tsss, pre_pos, cdss, genes):
    '''find the TSS which be detected in manual detection and TSSpredator'''
    if overlap:
        if pre_tss:
            pre_tss.attributes["print"] = True
            tss = pre_tss
        else:
            tss = tss_p
        tss.attribute_string = define_attributes(tss)
        tss.attributes["method"] = "TSSpredator,manual"
        if (not length) or \
           (tss.start <= int(length)):
            num_strain[tss.seq_id]["overlap"] += 1
            if (pre_pos != -1):
                if (tss.start - pre_pos != 0):
                    tsss["merge"].append(tss)
                    nums["tss"] += 1
                    overlap_num += 1
                else:
                    overlap_num += 1
            else:
                tsss["merge"].append(tss)
                nums["tss"] += 1
                overlap_num += 1
        overlap = False
        pre_pos = tss.start
    else:
        if tss_m.seq_id == tss_p.seq_id:
            tss_entry = compare_tss_gene(tss_m, genes)
            tss_m.attributes = tss_entry[1]
            tss_m.attribute_string = tss_entry[0]
            tss_m.attributes["method"] = "manual"
            tsss["merge"].append(tss_m)
            if (not length) or \
               (tss_m.start <= int(length)):
                num_strain[tss_m.seq_id]["manual"] += 1
                nums["tss_m"] += 1
                nums["tss"] += 1
    return (overlap, pre_pos, overlap_num)


def intersection(tsss, cluster, nums, lengths, cdss, genes, seqs):
    '''compare the predicted TSS and manual TSS'''
    num_strain = {}
    overlap = False
    overlap_num = 0
    pre_pos = -1
    start = False
    length = None
    for tss_m in tsss["tsss_m"]:
        pre_tss = None
        start = False
        if "all" in lengths.keys():
            length = seqs[tss_m.seq_id]
        else:
            if lengths[tss_m.seq_id] == "all":
                length = seqs[tss_m.seq_id]
            else:
                length = lengths[tss_m.seq_id]
        for tss_p in tsss["tsss_p"]:
            start = True
            if (tss_p.strand == tss_m.strand) and \
               (tss_p.seq_id == tss_m.seq_id):
                if (tss_p.start == tss_m.start):
                    tss_p.attributes["print"] = True
                    overlap = True
                    pre_tss = None
                    break
                elif (math.fabs(tss_p.start - tss_m.start) <= cluster):
                    overlap = True
                    pre_tss = tss_p
        if start:
            if tss_p.seq_id not in num_strain.keys():
                num_strain[tss_p.seq_id] = {"overlap": 0, "tsspredator": 0,
                                            "manual": 0}
            datas = check_overlap(overlap, pre_tss, nums, length, num_strain,
                                  overlap_num, tss_m, tss_p, tsss, pre_pos,
                                  cdss, genes)
            overlap = datas[0]
            pre_pos = datas[1]
            overlap_num = datas[2]
    if (start) or (len(tsss["tsss_m"]) == 0):
        for tss_p in tsss["tsss_p"]:
            run = False
            if not tss_p.attributes["print"]:
                tss_p.attribute_string = define_attributes(tss_p)
                tsss["merge"].append(tss_p)
                if (length is None):
                    run = True
                else:
                   if (tss_p.start <= int(length)):
                       run = True
                if run and (tss_p.seq_id in num_strain):
                    num_strain[tss_p.seq_id]["tsspredator"] += 1
                    nums["tss"] += 1
                    nums["tss_p"] += 1
    return overlap_num, num_strain


def print_file(final_tsss, program, out_gff):
    num_final = 0
    out = open(out_gff, "w")
    out.write("##gff-version 3\n")
    for tss in final_tsss:
        if "print" in tss.attributes.keys():
            del tss.attributes["print"]
        tss.attributes["ID"] = "_".join([
            tss.seq_id, program.lower() + str(num_final)])
        num_final += 1
        if program == "TSS":
            strand = Helper().get_strand_name(tss.strand)
            tss.attributes["Name"] = "TSS:" + "_".join(
                                              [str(tss.start), strand])
        else:
            strand = Helper().get_strand_name(tss.strand)
            tss.attributes["Name"] = "processing:" + "_".join(
                                              [str(tss.start), strand])
        tss.attribute_string = ";".join(
            ["=".join(items) for items in tss.attributes.items()])
        out.write("\t".join([str(field) for field in [
                             tss.seq_id, "ANNOgesic", tss.feature, tss.start,
                             tss.end, tss.score, tss.strand, tss.phase,
                             tss.attribute_string]]) + "\n")


def read_seq(seq_file):
    seqs = {}
    with open(seq_file) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                strain = line[1:]
                seqs[strain] = 0
            else:
                seqs[strain] = seqs[strain] + len(line)
    return seqs

def merge_manual_predict_tss(tss_predict_file, stat_file, out_gff,
                             gff_file, args_tss, manual, seq_file):
    '''merge the manual detected TSS and TSSpredator predicted TSS'''
    nums = {"tss_p": 0, "tss_m": 0, "tss": 0}
    merge_libs(args_tss.libs, args_tss.wig_folder, args_tss.program)
    wigs_f = read_wig("merge_forward.wig", "+")
    wigs_r = read_wig("merge_reverse.wig", "-")
    seqs = read_seq(seq_file)
    tsss, cdss, genes, = read_gff(tss_predict_file, manual,
                                  gff_file, args_tss.strain_lengths)
    overlap_num, num_strain = intersection(tsss, args_tss.cluster, nums,
                                           args_tss.strain_lengths,
                                           cdss, genes, seqs)
    sort_tsss = sorted(tsss["merge"], key=lambda k: (k.seq_id, k.start,
                                                     k.end, k.strand))
    final_tsss = fix_primary_type(sort_tsss, wigs_f, wigs_r)
    print_file(final_tsss, args_tss.program, out_gff)
    print_stat(num_strain, stat_file, overlap_num, nums)
