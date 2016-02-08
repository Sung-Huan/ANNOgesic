import sys
import csv
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
plt.style.use('ggplot')
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.helper import Helper

def plot(utr, utr_pri, utr_sec, filename, source, utr_type, base_5utr):
    bin_num = np.arange(0, 300, 5)
    if utr_type == "5utr":
        if source and (base_5utr != "transcript"):
            n, bins, hist1 = plt.hist(utr, bin_num,
                             color="#FF9999", label='secondary')
            n, bins, hist2 = plt.hist(utr_pri, bin_num,
                             color="#9999FF", label='primary')
            plt.legend((hist2[0], hist1[0]), ("Primary TSSs", "Secondary TSSs"))
        else:
            n, bins, hist1 = plt.hist(utr, bin_num, color="#FF9999")
        plt.xlabel("5'UTR_length")
    elif utr_type == "3utr":
        n, bins, hist = plt.hist(utr, bin_num, color="#9999FF", label='3\'UTR')
        plt.legend([hist[0]], ["3'UTR"])
        plt.xlabel("3'UTR_length")
    plt.ylabel("Amount")
    plt.savefig(filename)
    plt.clf()

def get_feature(cds):
    if "protein_id" in cds.attributes.keys():
        cds_name = cds.attributes["protein_id"]
    elif "locus_tag" in cds.attributes.keys():
        cds_name = cds.attributes["locus_tag"]
    else:
        strand = Helper().get_strand_name(cds.strand)
        cds_name = "".join([cds.feature, ":", str(cds.start),
                            "-", str(cds.end), "_", strand])
    return cds_name

def check_ta(tas, utr_start, utr_end, seq_id, strand):
    detect = False
    for ta in tas:
        if (ta.seq_id == seq_id) and \
           (ta.strand == strand):
            if (ta.start <= utr_start) and \
               (ta.end >= utr_end):
                detect = True
                break
    return detect

def import_utr(source, tss, utr_strain, utr_all,
               start, end, tas, length, base_5utr):
    if source:
        if "Primary" in tss.attributes["type"]:
            if base_5utr.lower() == "both":
                detect = check_ta(tas, start, end, tss.seq_id, tss.strand)
            elif base_5utr.lower() == "tss":
                detect = True
            if detect:
                utr_strain["pri"][tss.seq_id].append(length)
                utr_all["pri"].append(length)
        elif "Secondary" in tss.attributes["type"]:
            if base_5utr.lower() == "both":
                detect = check_ta(tas, start, end, tss.seq_id, tss.strand)
            elif base_5utr.lower() == "tss":
                detect = True
            if detect:
                utr_strain["sec"][tss.seq_id].append(length)
                utr_all["sec"].append(length)
    else:
        if base_5utr.lower() == "both":
            detect = check_ta(tas, start, end, tss.seq_id, tss.strand)
        elif base_5utr.lower() == "tss":
            detect = True
    if detect:
        utr_strain["all"][tss.seq_id].append(length)
        utr_all["all"].append(length)
    return detect

def get_5utr(tss, near_cds, utr_strain, utr_all, tas, num_utr,
             cds_name, locus_tag, source, out, base_5utr):
    if tss.strand == "+":
        start = tss.start
        end = near_cds.start
        length = end - start
        detect = import_utr(source, tss, utr_strain, utr_all,
                            start, end, tas, length, base_5utr)
    else:
        start = near_cds.end
        end = tss.end
        length = end - start
        detect = import_utr(source, tss, utr_strain, utr_all,
                            start, end, tas, length, base_5utr)
    if detect:
        name_utr = '%0*d' % (5, num_utr)
        if source:
            attribute_string = ";".join(
                 ["=".join(items) for items in [
                  ("ID", "_".join(["utr5", str(num_utr)])),
                  ("Name", "_".join(["5'UTR", name_utr])),
                  ("length", str(length)),
                  ("TSS_type", tss.attributes["type"]),
                  ("associated_cds", cds_name),
                  ("associated_gene", locus_tag),
                  ("associated_tss", tss.attributes["Name"])]])
            out.write("{0}\tANNOgesic\t5UTR\t{1}\t{2}\t.\t{3}\t.\t{4}\n".format(
                  tss.seq_id, start, end,
                  tss.strand, attribute_string))
        else:
            attribute_string = ";".join(
                 ["=".join(items) for items in [
                  ("ID", "_".join(["utr5", str(num_utr)])),
                  ("Name", "_".join(["5'UTR", name_utr])),
                  ("length", str(length)),
                  ("associated_cds", cds_name),
                  ("associated_gene", locus_tag),
                  ("associated_tss",
                  "".join(["TSS:", str(tss.start), "_", tss.strand]))]])
            out.write("{0}\tANNOgesic\t5UTR\t{1}\t{2}\t.\t{3}\t.\t{4}\n".format(
                  tss.seq_id, start, end,
                  tss.strand, attribute_string))
        num_utr += 1
    return num_utr

def detect_cds(cdss, gene):
    detect = False
    for cds in cdss:
        if "Parent" in cds.attributes.keys():
            if gene.attributes["ID"] == cds.attributes["Parent"]:
                detect = True
                near_cds = cds
                check_utr = True
                cds_name = get_feature(cds)
        else:
            if "locus_tag" in cds.attributes.keys():
                if gene.attributes["locus_tag"] == cds.attributes["locus_tag"]:
                    near_cds = cds
                    cds_name = cds.attributes["locus_tag"]
                    check_utr = True
                    detect = True
            elif (gene.seq_id == cds.seq_id) and \
                 (gene.strand == cds.strand):
                if (gene.start >= cds.start) and \
                   (gene.end <= cds.end):
                    near_cds = cds
                    cds_name = get_feature(cds)
                    check_utr = True
                    detect = True
    if not detect:
        return None, None, False
    else:
        return near_cds, cds_name, check_utr

def read_file(tss_file, gff_file, ta_file, term_file):
    genes = []
    cdss = []
    terms = []
    tsss = []
    tas = []
    gff_f = open(gff_file, "r")
    for entry in Gff3Parser().entries(gff_f):
        if (entry.feature == "CDS") or \
           (entry.feature == "tRNA") or \
           (entry.feature == "rRNA"):
            cdss.append(entry)
        elif (entry.feature == "gene"):
            genes.append(entry)
    for entry in Gff3Parser().entries(open(ta_file)):
        tas.append(entry)
    if term_file is not None:
        for entry_term in Gff3Parser().entries(open(term_file)):
            terms.append(entry_term)
        terms = sorted(terms, key=lambda k: (k.seq_id, k.start))
    if tss_file is not None:
        for entry in Gff3Parser().entries(open(tss_file)):
            tsss.append(entry)
    genes = sorted(genes, key=lambda k: (k.seq_id, k.start))
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start))
    tas = sorted(tas, key=lambda k: (k.seq_id, k.start))
    tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start))
    return genes, cdss, terms, tsss, tas

def check_associated_TSSpredator(genes, tss, cdss, check_utr, cds_name, locus):
    near_cds = None
    for gene in genes:
        if (tss.seq_id == gene.seq_id) and (
            tss.strand == gene.strand):
            if "locus_tag" in gene.attributes.keys():
                if gene.attributes["locus_tag"] == locus:
                    for cds in cdss:
                        if "Parent" in cds.attributes.keys():
                            if gene.attributes["ID"] == \
                               cds.attributes["Parent"]:
                                near_cds = cds
                                check_utr = True
                                cds_name = get_feature(cds)
                        else:
                            if "locus_tag" in cds.attributes.keys():
                                if gene.attributes["locus_tag"] == \
                                   cds.attributes["locus_tag"]:
                                    near_cds = cds
                                    cds_name = cds.attributes["locus_tag"]
                                    check_utr = True
                            elif (gene.seq_id == cds.seq_id) and \
                                 (gene.strand == cds.strand):
                                if (gene.start >= cds.start) and \
                                   (gene.end <= cds.end):
                                    near_cds = cds
                                    cds_name = get_feature(cds)
                                    check_utr = True
    return check_utr, cds_name, near_cds

def get_5utr_from_TSSpredator(tss, genes, cdss):
    check_utr = False
    cds_name = "NA"
    if ("Primary" in tss.attributes["type"]) or \
       ("Secondary" in tss.attributes["type"]):
        ass_gene = tss.attributes["associated_gene"].split("&")
        tss_type = tss.attributes["type"].split("&")
        for index in range(len(tss_type)):
            if (tss_type[index] == "Primary") or \
               (tss_type[index] == "Secondary"):
                locus_tag = ass_gene[index]
                check_utr, cds_name, near_cds = check_associated_TSSpredator(
                                                genes, tss, cdss, check_utr,
                                                cds_name, locus_tag)
                if not check_utr:
                    for cds in cdss:
                        if (cds.seq_id == tss.seq_id) and \
                           (cds.strand == tss.strand):
                            strand = Helper().get_strand_name(cds.strand)
                            if locus_tag == (cds.feature + ":" + \
                                             str(cds.start) + \
                                             "-" + str(cds.end) + \
                                             "_" + strand):
                                near_cds = cds
                                cds_name = get_feature(cds)
                                check_utr = True
        if check_utr:
            utr_datas = {"check": check_utr, "cds_name": cds_name,
                         "near_cds": near_cds, "locus": locus_tag}
        else:
            utr_datas = {"check": False, "cds_name": None,
                         "near_cds": None, "locus": None}
    else:
        utr_datas = {"check": False, "cds_name": None,
                     "near_cds": None, "locus": None}
    return utr_datas

def get_5utr_from_other(tss, genes, cdss):
    check_utr = False
    cds_name = "NA"
    for gene in genes:
        if (tss.seq_id == gene.seq_id) and \
           (tss.strand == gene.strand):
            if (tss.strand == "+") and \
               ((gene.start - tss.start) <= 300) and \
               (tss.start <= gene.start):
                if "locus_tag" in gene.attributes.keys():
                    locus_tag = gene.attributes["locus_tag"]
                else:
                    locus_tag = gene.attributes["ID"]
                near_cds, cds_name, check_utr = detect_cds(cdss, gene)
                break
            elif (tss.strand == "-") and \
                 ((tss.start - gene.end) <= 300) and \
                 (tss.start >= gene.end):
                if "locus_tag" in gene.attributes.keys():
                    locus_tag = gene.attributes["locus_tag"]
                else:
                    locus_tag = gene.attributes["ID"]
                near_cds, cds_name, check_utr = detect_cds(cdss, gene)
                break
    if check_utr:
        utr_datas = {"check": check_utr, "cds_name": cds_name,
                     "near_cds": near_cds, "locus": locus_tag}
    else:
        utr_datas = {"check": False, "cds_name": None,
                     "near_cds": None, "locus": None}
    return utr_datas

def compare_ta_tss(tsss, ta, source):
    tss_name = "NA"
    for tss in tsss:
        if (tss.seq_id == ta.seq_id) and \
           (tss.strand == ta.strand):
            if ta.strand == "+":
                if math.fabs(tss.start - ta.start) <= 2:
                    if source:
                        tss_name = tss.attributes["Name"]
                    else:
                        tss_name = "".join(["TSS:", str(tss.start),
                                            "_", tss.strand])
                    break
            else:
                if math.fabs(tss.start - ta.end) <= 2:
                    if source:
                        tss_name = tss.attributes["Name"]
                    else:
                        tss_name = "".join(["TSS:", str(tss.start),
                                            "_", tss.strand])
                    break
    return tss_name

def get_gene_info(genes, cds):
    gene_name = "NA"
    for gene in genes:
        if (cds.seq_id == gene.seq_id) and \
           (cds.strand == gene.strand) and \
           (cds.start >= gene.start) and \
           (cds.end <= gene.end):
            if "locus_tag" in gene.attributes.keys():
                gene_name = gene.attributes["locus_tag"]
            else:
                strand = Helper().get_strand_name(gene.strand)
                gene_name = "".join([gene.feature, ":", str(gene.start),
                                     "-", str(gene.end), "_", strand])
            break
    return gene_name

def compare_ta_cds_gene(cds, ta, tsss, genes, utr_start,
                        utr_end, pre_pos, source):
    detect = False
    if ((utr_end - utr_start) >= 0) and \
       ((utr_end - utr_start) <= 300):
        cds_name = get_feature(cds)
        gene_name = "NA"
        length = utr_end - utr_start
        detect = True
        start = utr_start
        end = utr_end
        gene_name = get_gene_info(genes, cds)
    elif ((utr_end - utr_start) >= 300):
        cds_name = get_feature(cds)
        for tss in tsss:
            length = 400
            if (tss.seq_id == ta.seq_id) and \
               (tss.strand == ta.strand):
                if ((tss.start >= utr_start) and (tss.start <= utr_end)):
                    length = utr_end - tss.start
                    utr_start = tss.start
                elif (tss.strand == "+") and \
                     (tss.start >= pre_pos) and \
                     (tss.start <= utr_end) and \
                     (pre_pos is not None):
                    length = utr_end - pre_pos
                    utr_start = pre_pos
                elif (tss.strand == "-") and \
                     (tss.start <= pre_pos) and \
                     (tss.start >= utr_start) and \
                     (pre_pos is not None):
                    length = pre_pos - utr_start
                    utr_end = pre_pos
                if length <= 300:
                    detect = True
                    if source:
                        tss_name = tss.attributes["Name"]
                    else:
                        tss_name = "".join(["TSS:", str(tss.start),
                                            "_", tss.strand])
                    break
        start = utr_start
        end = utr_end
        gene_name = get_gene_info(genes, cds)
    return (cds_name, gene_name, detect, length, start, end)

def get_feature_data(ta, cdss, tsss, genes, source):
    datas = None
    pre_pos = None
    if ta.strand == "+":
        for cds in cdss:
            if (ta.seq_id == cds.seq_id) and (
                ta.strand == cds.strand):
                if (cds.start >= ta.start) and (
                    cds.start <= ta.end):
                    datas = compare_ta_cds_gene(cds, ta, tsss, genes,
                                    ta.start, cds.start, pre_pos, source)
                    if datas[2]:
                        break
                pre_pos = cds.end
    else:
        for cds in reversed(cdss):
            if (ta.seq_id == cds.seq_id) and (
                ta.strand == cds.strand):
                if (cds.end <= ta.end) and (
                    cds.end >= ta.start):
                    datas = compare_ta_cds_gene(cds, ta, tsss, genes,
                                    cds.end, ta.end, pre_pos, source)
                    if datas[2]:
                        break
                pre_pos = cds.start
    return datas

def compare_ta(tas, tsss, genes, cdss, utr_strain, utr_all, out, source):
    num_utr = 0
    for ta in tas:
        if ta.seq_id not in utr_strain["all"].keys():
            utr_strain["all"][ta.seq_id] = []
        if ta.seq_id not in utr_strain["pri"].keys():
            utr_strain["pri"][ta.seq_id] = []
        if ta.seq_id not in utr_strain["sec"].keys():
            utr_strain["sec"][ta.seq_id] = []
        tss_name = compare_ta_tss(tsss, ta, source)
        datas = get_feature_data(ta, cdss, tsss, genes, source)
        if datas is not None:
            cds_name = datas[0]
            gene_name = datas[1]
            detect = datas[2]
            length = datas[3]
            start = datas[4]
            end = datas[5]
        else:
            detect = False
        if detect:
            utr_strain["all"][ta.seq_id].append(length)
            utr_all["all"].append(length)
            name_utr = '%0*d' % (5, num_utr)
            attribute_string = ";".join(
                 ["=".join(items) for items in [
                  ("ID", "_".join(["utr5", str(num_utr)])),
                  ("Name", "_".join(["5'UTR", name_utr])),
                  ("length", str(length)),
                  ("associated_cds", cds_name),
                  ("associated_gene", gene_name),
                  ("associated_tss", tss_name)]])
            out.write("{0}\tANNOgesic\t5UTR\t{1}\t{2}\t.\t{3}\t.\t{4}\n".format(
                  ta.seq_id, start, end,
                  ta.strand, attribute_string))
            num_utr += 1

def detect_5utr(tss_file, gff_file, ta_file, source, base_5utr, out_file):
    num_utr = 0
    utr_all = {"all": [], "pri": [], "sec": []}
    utr_strain = {"all": {}, "pri": {}, "sec": {}}
    pre_seq_id = ""
    out = open(out_file, "w")
    out.write("##gff-version 3\n")
    genes, cdss, terms, tsss, tas = read_file(tss_file, gff_file,
                                              ta_file, None)
    if (base_5utr.upper() == "TSS") or (base_5utr.lower() == "both"):
        for tss in tsss:
            check_utr = False
            cds_name = "NA"
            if source:
                utr_datas = get_5utr_from_TSSpredator(tss, genes, cdss)
            else:
                utr_datas = get_5utr_from_other(tss, genes, cdss)
            if utr_datas["check"]:
                if tss.seq_id != pre_seq_id:
                    pre_seq_id = tss.seq_id
                    utr_strain["pri"][tss.seq_id] = []
                    utr_strain["sec"][tss.seq_id] = []
                    utr_strain["all"][tss.seq_id] = []
                num_utr = get_5utr(tss, utr_datas["near_cds"], utr_strain,
                                   utr_all, tas, num_utr,
                                   utr_datas["cds_name"], utr_datas["locus"],
                                   source, out, base_5utr)
    if base_5utr.lower() == "transcript":
        compare_ta(tas, tsss, genes, cdss, utr_strain, utr_all, out, source)
    name = (gff_file.split("/"))[-1].replace(".gff", "")
    plot(utr_all["all"], utr_all["pri"], utr_all["sec"],
         "_".join([name, "all_5utr_length.png"]), source, "5utr", base_5utr)
    if len(utr_strain["all"]) > 1:
        for strain in utr_strain["all"].keys():
            plot(utr_strain["all"][strain], utr_strain["pri"][strain],
                 utr_strain["sec"][strain],
                 "_".join([strain, "5utr_length.png"]),
                 source, "5utr", base_5utr)
    out.close()

def compare_term(ta, terms, fuzzy):
    for term in terms:
        if ta.strand == term.strand:
            if term.strand == "+":
                if (math.fabs(ta.end - term.start) <= fuzzy) or \
                   ((ta.end >= term.start) and (ta.end <= term.end)):
                    return term
            else:
                if (math.fabs(ta.start - term.end) <= fuzzy) or \
                   ((ta.start >= term.start) and (ta.start <= term.end)):
                    return term

def get_3utr(ta, near_cds, utr_all, utr_strain, attributes, num_utr, out):
    if ta.strand == "+":
        start = near_cds.end
        end = ta.end
        length = ta.end - near_cds.end
        utr_all.append(length)
        utr_strain[ta.seq_id].append(length)
    else:
        start = ta.start
        end = near_cds.start
        length = near_cds.start - ta.start
        utr_all.append(length)
        utr_strain[ta.seq_id].append(length)
    attributes.append("=".join(["length", str(length)]))
    attributes.append("=".join(["associated_tran", ta.attributes["ID"]]))
    attribute = ";".join(attributes)
    if (length <= 300) and (length > 0):
        name_utr = '%0*d' % (5, num_utr)
        name = "=".join(["Name", "_".join(["3'UTR", ta.attributes["ID"]])])
        id_ = "ID=utr3_" + str(num_utr)
        num_utr += 1
        attribute_string = ";".join([id_, name, attribute])
        out.write("\t".join([ta.seq_id, "ANNOgesic", "3UTR",
                             str(start), str(end),
                             ta.score, ta.strand, ta.phase,
                             attribute_string]) + "\n")
    return num_utr

def get_near_cds(cdss, genes, ta, attributes):
    first = True
    detect = False
    for cds in cdss:
        if (ta.seq_id == cds.seq_id) and \
           (ta.strand == cds.strand):
            if first:
                near_cds = cds
                first = False
            if ta.strand == "+":
                if (cds.end <= ta.end) and \
                   (cds.start >= ta.start) and \
                   (cds.end > near_cds.end):
                    near_cds = cds
                    detect = True
            else:
                if (cds.start >= ta.start) and \
                   (cds.end <= ta.end):
                    near_cds = cds
                    detect = True
                    break
    if detect:
        for gene in genes:
            if ("Parent" in near_cds.attributes.keys()):
                if near_cds.attributes["Parent"] == gene.attributes["ID"]:
                    attributes.append("=".join(["associated_gene",
                                      gene.attributes["locus_tag"]]))
        if "protein_id" in near_cds.attributes.keys():
            attributes.append("=".join(["associated_cds",
                              near_cds.attributes["protein_id"]]))
        elif "locus_tag" in near_cds.attributes.keys():
            attributes.append("=".join(["associated_cds",
                              near_cds.attributes["locus_tag"]]))
        else:
            attributes.append("=".join(["associated_cds",
                              "_".join([near_cds.feature,
                              str(near_cds.start), str(near_cds.end),
                              near_cds.strand])]))
    return near_cds

def detect_3utr(ta_file, gff_file, term_file, fuzzy, out_file):
    num_utr = 0
    utr_all = []
    utr_strain = {}
    pre_seq_id = ""
    out = open(out_file, "w")
    out.write("##gff-version 3\n")
    genes, cdss, terms, tsss, tas = read_file(None, gff_file,
                                              ta_file, term_file)
    for ta in tas:
        if term_file is not None:
            term = compare_term(ta, terms, fuzzy)
        else:
            term = None
        attributes = []
        if term is not None:
            attributes.append("=".join(["associated_term",
                              term.attributes["ID"]]))
        near_cds = get_near_cds(cdss, genes, ta, attributes)
        if ta.seq_id != pre_seq_id:
            pre_seq_id = ta.seq_id
            utr_strain[ta.seq_id] = []
        num_utr = get_3utr(ta, near_cds, utr_all, utr_strain,
                           attributes, num_utr, out)
    name = (gff_file.split("/"))[-1].replace(".gff", "")
    plot(utr_all, None, None, "_".join([name, "all_3utr_length.png"]),
                                        None, "3utr", "3utr")
    out.close()
    if len(utr_strain) > 1:
        for strain in utr_strain.keys():
            plot(utr_strain[strain], None, None,
                 "_".join([strain, "3utr_length.png"]), None, "3utr", "3utr")
