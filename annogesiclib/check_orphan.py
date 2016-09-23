import os
from annogesiclib.helper import Helper
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.parser_wig import WigParser


def read_gff(gff_file, features):
    gffs = []
    if not os.path.isfile(gff_file):
        filename = gff_file.split(".")
        gff_file = ".".join(filename[0:-2]) + ".gff"
    g_f = open(gff_file, "r")
    for entry in Gff3Parser().entries(g_f):
        if entry.feature in features:
            gffs.append(entry)
    gffs = sorted(gffs, key=lambda k: (k.seq_id, k.start))
    return gffs


def is_primary(cds_start, cds_end, tss_pos, strand):
    '''check primary TSS'''
    if strand == "+":
        if is_utr(cds_start, tss_pos, 300) and (cds_start >= tss_pos):
            return True
    else:
        if is_utr(tss_pos, cds_end, 300) and (cds_end <= tss_pos):
            return True


def is_internal(cds_start, cds_end, tss_pos, strand):
    '''check internal TSS'''
    if ((cds_start < tss_pos) and (cds_end > tss_pos)) or (
            (strand == "+") and (tss_pos == cds_end)) or (
            (strand == "-") and (tss_pos == cds_start)):
        return True


def is_antisense(cds_start, cds_end, tss_pos, strand):
    '''check antisense TSS'''
    if ((is_utr(cds_start, tss_pos, 100)) and (cds_start >= tss_pos)) or (
            (is_utr(tss_pos, cds_end, 100)) and (cds_end <= tss_pos)) or (
            is_internal(cds_start, cds_end, tss_pos, strand)):
        return True


def is_utr(pos1, pos2, length):
    '''check the utr'''
    if pos1 - pos2 <= length:
        return True


def get_attributes(tss, cds):
    if tss.attributes["associated_gene"] == "orphan":
        if "locus_tag" in cds.attributes.keys():
            tss.attributes["associated_gene"] = cds.attributes["locus_tag"]
        else:
            strand = Helper().get_strand_name(cds.strand)
            tss.attributes["associated_gene"] = cds.feature + ":" + \
                str(cds.start) + "-" + str(cds.end) + "_" + strand
    else:
        if "locus_tag" in cds.attributes.keys():
            tss.attributes["associated_gene"] = "&".join([
                tss.attributes["associated_gene"],
                cds.attributes["locus_tag"]])
        else:
            strand = Helper().get_strand_name(cds.strand)
            tss.attributes["associated_gene"] = "&".join([
                tss.attributes["associated_gene"],
                cds.feature + ":" + str(cds.start) + "-" +
                str(cds.end) + "_" + strand])


def detect_coverage(wigs, tss, ref):
    '''compare the coverage of TSS in order to get
    proper primary TSS'''
    tss_cover = -1
    ref_cover = -1
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


def del_repeat(tsss):
    '''delete redundant assigned types of TSS'''
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
            elif (type_ == "Antisense") or (type_ == "Internal") or (
                    type_ == "Orphan"):
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


def fix_attributes(tss, tss_entry):
    '''change primary TSS to secondary TSS'''
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
    '''Deal with the multiple primary TSS of one gene.
    change the low expressed one to be secondary TSS'''
    for tss in tsss:
        if "Primary" in tss.attributes["type"]:
            tss_entrys = get_primary_locus_tag(tss)
            for ref in tsss:
                if (ref.seq_id == tss.seq_id) and (
                    ref.strand == tss.strand) and (
                        ref.start == tss.start):
                    pass
                else:
                    if "Primary" in ref.attributes["type"]:
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
                                        if tss_entry["utr"] < ref_entry["utr"]:
                                            fix_attributes(ref, ref_entry)
                                        elif (tss_entry["utr"] >
                                              ref_entry["utr"]):
                                            fix_attributes(tss, tss_entry)
    del_repeat(tsss)
    return tsss


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


def get_attributes_int_anti(tss, cds, type_):
    '''import useful information to attributes'''
    if tss.attributes["type"] != "Orphan":
        tss.attributes["type"] = "&".join(
                [tss.attributes["type"], type_])
        tss.attributes["UTR_length"] = "&".join(
                [tss.attributes["UTR_length"],
                 type_ + "_NA"])
    else:
        tss.attributes["type"] = type_
        tss.attributes["UTR_length"] = type_ + "_NA"
    get_attributes(tss, cds)


def compare_cds_check_orphan(tsss, cdss):
    '''main part of checking all orphan TSS'''
    for tss in tsss:
        if tss.attributes["type"] == "Orphan":
            for cds in cdss:
                if (tss.seq_id == cds.seq_id) and \
                   (tss.strand == cds.strand):
                    if is_primary(cds.start, cds.end, tss.start, tss.strand):
                        if tss.attributes["type"] != "Orphan":
                            tss.attributes["type"] = "&".join(
                                    [tss.attributes["type"], "Primary"])
                            if tss.strand == "+":
                                tss.attributes["UTR_length"] = "&".join([
                                    tss.attributes["UTR_length"],
                                    "Primary_" + str(cds.start - tss.start)])
                            else:
                                tss.attributes["UTR_length"] = "&".join([
                                    tss.attributes["UTR_length"],
                                    "Primary_" + str(tss.start - cds.end)])
                        else:
                            tss.attributes["type"] = "Primary"
                            if tss.strand == "+":
                                tss.attributes["UTR_length"] = (
                                    "Primary_" + str(cds.start - tss.start))
                            else:
                                tss.attributes["UTR_length"] = (
                                    "Primary_" + str(tss.start - cds.end))
                        get_attributes(tss, cds)
                    if is_internal(cds.start, cds.end, tss.start, tss.strand):
                        if "locus_tag" in cds.attributes.keys():
                            if (cds.attributes["locus_tag"] not in
                                    tss.attributes["associated_gene"]):
                                get_attributes_int_anti(tss, cds, "Internal")
                        else:
                            strand = Helper().get_strand_name(cds.strand)
                            if ("".join([cds.feature, ":", str(cds.start),
                                "-", str(cds.end), "_", strand]) not in
                                    tss.attributes["associated_gene"]):
                                get_attributes_int_anti(tss, cds, "Internal")
                    if is_antisense(cds.start, cds.end, tss.start, tss.strand):
                        if "locus_tag" in cds.attributes.keys():
                            if (cds.attributes["locus_tag"] not in
                                    tss.attributes["associated_gene"]):
                                get_attributes_int_anti(tss, cds, "Antisense")
                        else:
                            strand = Helper().get_strand_name(cds.strand)
                            if ("".join([cds.feature, ":", str(cds.start),
                                "-", str(cds.end), "_", strand]) not in
                                    tss.attributes["associated_gene"]):
                                get_attributes_int_anti(tss, cds, "Antisense")


def check_orphan(tss_file, gff_file, wig_f_file, wig_r_file, out_gff):
    '''If the genome annotation gff file has no locus tag, TSSpredator
    will classify all TSS into orphan. It is for fixing this mistake.
    It will compare the TSS and gene to classify the TSS.'''
    cdss = read_gff(gff_file, ["CDS", "tRNA", "rRNA"])
    tsss = read_gff(tss_file, ["TSS"])
    wigs_f = read_wig(wig_f_file, "+")
    wigs_r = read_wig(wig_r_file, "-")
    out = open(out_gff, "w")
    out.write("##gff-version 3\n")
    compare_cds_check_orphan(tsss, cdss)
    final_tsss = fix_primary_type(tsss, wigs_f, wigs_r)
    for tss in final_tsss:
        tss.attribute_string = ";".join(
            ["=".join(items) for items in tss.attributes.items()])
        out.write("\t".join([str(field) for field in [
                        tss.seq_id, tss.source, tss.feature, tss.start,
                        tss.end, tss.score, tss.strand, tss.phase,
                        tss.attribute_string]]) + "\n")
