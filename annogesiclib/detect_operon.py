import math
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.helper import Helper


def import_to_operon(start, end, strand):
    return {"start": start, "end": end, "strand": strand}


def get_gene_info(cds):
    if "locus_tag" in cds.attributes.keys():
        feature = cds.attributes["locus_tag"]
    else:
        strand = Helper().get_strand_name(cds.strand)
        feature = "".join([cds.feature, ":", str(cds.start),
                           "-", str(cds.end), "_", strand])
    return feature


def get_term_feature(ta, data, term_fuzzy, features, datas,
                     ta_check_point, data_check_start, data_check_end):
    '''verify and get proper terminator to operon'''
    jump = False
    if (ta.strand == data.strand) and (
            ta.seq_id == data.seq_id) and (
            (math.fabs(data.start - ta_check_point) <= term_fuzzy) or (
            math.fabs(data.end - ta_check_point) <= term_fuzzy) or (
            (ta_check_point >= data.start) and (
            ta_check_point <= data.end))):
        features["detect"] = True
    if (ta.strand == data.strand) and (
            ta.seq_id == data.seq_id):
        if (ta.start <= data_check_start) and (
                ta.end >= data_check_end):
            features["num"] += 1
            datas.append(data)
        elif (ta_check_point >= data.start) and (
                ta_check_point <= data.end):
            features["num"] += 1
            datas.append(data)
    if (ta.seq_id == data.seq_id) and (
            data.start - term_fuzzy > ta.end):
        jump = True
    return jump


def get_tss_feature(ta, data, features, tss_fuzzy, datas, ta_check_point,
                    data_check_start, data_check_end):
    '''verify and get the proper TSS for operon'''
    jump = False
    if (ta.strand == data.strand) and (
            ta.seq_id == data.seq_id) and (
            math.fabs(ta_check_point - data.start) <= tss_fuzzy):
        features["detect"] = True
    if (ta.strand == data.strand) and (
            ta.seq_id == data.seq_id) and (
            ta.start <= data_check_start) and (
            ta.end >= data_check_end):
        features["num"] += 1
        datas.append(data)
    if (ta.seq_id == data.seq_id) and (
            data_check_end > ta.end):
        jump = True
    return jump


def detect_features(ta, inputs, feature, term_fuzzy, tss_fuzzy):
    '''Detect the feature which should group as a operon'''
    features = {"num": 0, "detect": False}
    datas = []
    for data in inputs:
        if (feature == "term"):
            if ta.strand == "+":
                jump_term = get_term_feature(ta, data, term_fuzzy, features,
                                             datas, ta.end, data.start,
                                             data.start - term_fuzzy)
            elif ta.strand == "-":
                jump_term = get_term_feature(ta, data, term_fuzzy, features,
                                             datas, ta.start,
                                             data.end + term_fuzzy, data.end)
            if jump_term:
                break
        elif (feature == "tss"):
            if ta.strand == "+":
                jump_tss = get_tss_feature(ta, data, features, tss_fuzzy,
                                           datas, ta.start,
                                           data.start + tss_fuzzy, data.end)
            elif ta.strand == "-":
                jump_tss = get_tss_feature(ta, data, features, tss_fuzzy,
                                           datas, ta.end, data.end,
                                           data.start - tss_fuzzy)
            if jump_tss:
                break
        else:
            if feature == "gene":
                if (ta.strand == data.strand) and (
                        ta.seq_id == data.seq_id) and (
                        data.feature == "gene"):
                    if ((ta.start <= data.start) and (
                             ta.end >= data.end)) or (
                            (ta.start >= data.start) and (
                             ta.end <= data.end)) or (
                            (ta.start >= data.start) and (
                             ta.start <= data.end) and (
                             ta.end >= data.end)) or (
                            (ta.start <= data.start) and (
                             ta.end <= data.end) and (
                             ta.end >= data.start)):
                        features["num"] += 1
                        features["detect"] = True
                        datas.append(data)
    return {"data_list": datas, "num_feature": features["num"],
            "with_feature": features["detect"]}


def check_conflict(genes, pos, strand):
    '''check TSS which is not primary or secondary TSS'''
    conflict = False
    for gene in genes["data_list"]:
        if (gene.strand == strand):
            if (gene.start < pos) and (
                    gene.end >= pos):
                conflict = True
                break
    return conflict


def check_gene(tsss, genes, strand, ta_pos, first, min_length, end,
               operons, operon_pos):
    '''Check TSS and annotated feature. It can group the feature and TSS to 
    be operon or sub-operon'''
    no_count_tsss = []
    for tss in tsss:
        if tss not in no_count_tsss:
            end_points = [ta_pos]
            for pos in tsss:
                if (pos not in no_count_tsss) and (
                        tss.start != pos.start):
                    end_points.append(pos.start)
            end_points.append(end)
            if tss.strand == "+":
                end_points.sort()
            else:
                end_points.sort(reverse=True)
            for point in end_points:
                detect_pos = False
                if tss.strand == "+":
                    for gene in genes["data_list"]:
                        if (gene.seq_id == tss.seq_id) and (
                                gene.strand == tss.strand):
                            if (gene.start >= tss.start) and (
                                    gene.end <= point):
                                detect_pos = True
                                break
                else:
                    for gene in genes["data_list"]:
                        if (gene.seq_id == tss.seq_id) and (
                                gene.strand == tss.strand):
                            if (gene.start >= point) and (
                                    gene.end <= tss.start):
                                detect_pos = True
                                break
                if not detect_pos:
                    no_count_tsss.append(tss)
                else:
                    operon_pos, first = compute_sub_operon(
                        strand, point, ta_pos, first, min_length,
                        end, operons, operon_pos)
                    break


def sub_operon_gene_conflict(tsss, strand, genes, ta_pos, first, min_length,
                             end, operons, operon_pos):
    '''remove the TSS which is not primary or secondary TSS of gene
    This TSS can not form sub-operon'''
    new_tsss = []
    for tss in tsss["data_list"]:
        conflict = check_conflict(genes, tss.start, strand)
        if not conflict:
            new_tsss.append(tss)
    check_gene(new_tsss, genes, strand, ta_pos, first,
               min_length, end, operons, operon_pos)


def sub_operon(strand, tsss, ta_pos, end, genes, min_length):
    '''verify the sub-operon'''
    first = True
    operons = []
    operon_pos = ta_pos
    if tsss["with_feature"]:
        if tsss["num_feature"] == 1:
            pass
        else:
            sub_operon_gene_conflict(
                tsss, strand, genes, ta_pos, first, min_length,
                end, operons, operon_pos)
    else:
        sub_operon_gene_conflict(
            tsss, strand, genes, ta_pos, first,
            min_length, end, operons, operon_pos)
    return operons


def compute_sub_operon(strand, point, ta_pos, first,
                       min_length, end, operons, operon_pos):
    '''For computting and import the sub-operon'''
    if first:
        operon_pos = ta_pos
        first = False
    if math.fabs(point - operon_pos) > min_length:
        if strand == "+":
            operons.append(import_to_operon(operon_pos,
                           point - 1, strand))
            operon_pos = point
        else:
            operons.append(import_to_operon(point + 1,
                           operon_pos, strand))
            operon_pos = point
    return operon_pos, first


def read_gff(ta_file, gff_file, tss_file, terminator_file):
    tas = []
    gffs = []
    tss_gffs = []
    term_gffs = []
    gff_parser = Gff3Parser()
    for ta in gff_parser.entries(open(ta_file)):
        tas.append(ta)
    for entry in gff_parser.entries(open(gff_file)):
        gffs.append(entry)
    if tss_file is not False:
        for entry in gff_parser.entries(open(tss_file)):
            tss_gffs.append(entry)
        tss_gffs = sorted(tss_gffs, key=lambda k: (k.seq_id, k.start,
                                                   k.end, k.strand))
    if terminator_file is not False:
        for entry in gff_parser.entries(open(terminator_file)):
            term_gffs.append(entry)
        term_gffs = sorted(term_gffs, key=lambda k: (k.seq_id, k.start,
                                                     k.end, k.strand))
    tas = sorted(tas, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    gffs = sorted(gffs, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    return tas, gffs, tss_gffs, term_gffs


def print_file(ta, operons, out, operon_id, whole_operon, tsss,
               terms, genes, whole_gene, out_g):
    attribute_string = ";".join(
             ["=".join(items) for items in [
              ("ID", "_".join([ta.seq_id, operon_id.replace("_", "")])),
              ("Name", operon_id),
              ("associated_gene", ",".join(whole_gene))]])
    out_g.write("{0}\tANNOgesic\toperon\t{1}"
                "\t{2}\t.\t{3}\t.\t{4}\n".format(
                    ta.seq_id, str(whole_operon.start), str(whole_operon.end),
                    whole_operon.strand, attribute_string))
    if len(operons) <= 1:
        out.write("\t".join([operon_id, ta.seq_id,
                  "-".join([str(whole_operon.start), str(whole_operon.end)]),
                  whole_operon.strand, "0", "NA",
                  str(tsss["with_feature"]), str(tsss["num_feature"]),
                  str(terms["with_feature"]), str(terms["num_feature"]), "NA",
                  str(genes["num_feature"]), "NA",
                  ", ".join(whole_gene)]) + "\n")
    else:
        for sub in operons:
            sub_gene = []
            num_sub_gene = 0
            for gene in genes["data_list"]:
                if (sub["strand"] == gene.strand) and (
                        sub["start"] <= gene.start) and (
                        sub["end"] >= gene.end):
                    if "locus_tag" in gene.attributes.keys():
                        sub_gene.append(gene.attributes["locus_tag"])
                    else:
                        sub_gene.append("".join([gene.feature, ":",
                                        str(gene.start), "-", str(gene.end),
                                        "_", gene.strand]))
                    num_sub_gene += 1
            if num_sub_gene == 0:
                sub_gene.append("NA")
            out.write("\t".join([operon_id, ta.seq_id,
                      "-".join([str(whole_operon.start),
                                str(whole_operon.end)]),
                      whole_operon.strand, str(len(operons)),
                      "-".join([str(sub["start"]), str(sub["end"])]),
                      str(tsss["with_feature"]), str(tsss["num_feature"]),
                      str(terms["with_feature"]), str(terms["num_feature"]),
                      str(num_sub_gene), str(genes["num_feature"]),
                      ", ".join(sub_gene), ", ".join(whole_gene)]) + "\n")


def operon(ta_file, tss_file, gff_file, terminator_file, tss_fuzzy,
           term_fuzzy, min_length, out_file, out_gff):
    '''main part for detection of operon'''
    out = open(out_file, "w")
    out_g = open(out_gff, "w")
    out_g.write("##gff-version 3\n")
    out.write("Operon_ID\tGenome\tOperon_position\tStrand\t")
    out.write("Number_of_suboperon\tPosition_of_suboperon\tStart_with_TSS\t")
    out.write("Number_of_TSS\tTerminated_with_terminator\t")
    out.write("Number_of_terminator\tNumber_of_gene_associated_suboperon\t")
    out.write("Number_of_gene_associated_operon\t")
    out.write("Associated_genes_with_suboperon\t")
    out.write("Associated_genes_with_whole_operon\n")
    num_operon = 0
    tas, gffs, tss_gffs, term_gffs = read_gff(ta_file, gff_file, tss_file,
                                              terminator_file)
    for ta in tas:
        whole_gene = []
        check_operon = False
        if (math.fabs(ta.start - ta.end) >= min_length):
            whole_operon = ta
            check_operon = True
        genes = detect_features(ta, gffs, "gene", term_fuzzy, tss_fuzzy)
        if len(tss_gffs) != 0:
            tsss = detect_features(ta, tss_gffs, "tss", term_fuzzy, tss_fuzzy)
        else:
            tsss = {"with_feature": "NA", "num_feature": "NA"}
        if terminator_file is None:
            terms = {"with_feature": "NA", "num_feature": "NA"}
        else:
            terms = detect_features(ta, term_gffs, "term",
                                    term_fuzzy, tss_fuzzy)
        if len(tss_gffs) != 0:
            if ta.strand == "+":
                operons = sub_operon(ta.strand, tsss, ta.start,
                                     ta.end, genes, min_length)
            else:
                operons = sub_operon(ta.strand, tsss, ta.end,
                                     ta.start, genes, min_length)
        else:
            operons = [{"start": ta.start, "end": ta.end, "strand": ta.strand}]
        if genes["num_feature"] != 0:
            for gene in genes["data_list"]:
                whole_gene.append(get_gene_info(gene))
        else:
            whole_gene.append("NA")
        if check_operon:
            if whole_gene != ["NA"]:
                operon_id = "Operon" + str(num_operon)
                num_operon += 1
                if len(tss_gffs) != 0:
                    print_file(ta, operons, out, operon_id, whole_operon,
                               tsss, terms, genes, whole_gene, out_g)
                else:
                    print_file(ta, operons, out, operon_id, ta,
                               tsss, terms, genes, whole_gene, out_g)
    out.close()
    out_g.close()
