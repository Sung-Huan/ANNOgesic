import math
import shutil
import matplotlib as mpl
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.helper import Helper
import numpy as np
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')


def plot(utr, utr_pri, utr_sec, filename, source, utr_type, base_5utr):
    '''plot the distribution of length of UTRs'''
    bin_num = np.arange(0, 300, 5)
    if utr_type == "5utr":
        if source and (base_5utr != "transcript"):
            n, bins, hist1 = plt.hist(utr, bin_num,
                                      color="#FF9999", label='secondary',
                                      edgecolor='black', linewidth=1)
            n, bins, hist2 = plt.hist(utr_pri, bin_num,
                                      color="#9999FF", label='primary',
                                      edgecolor='black', linewidth=1)
            plt.legend((hist2[0], hist1[0]),
                       ("Primary TSSs", "Secondary TSSs"))
        else:
            n, bins, hist1 = plt.hist(utr, bin_num, color="#FF9999")
        plt.xlabel("5'UTR_length")
    elif utr_type == "3utr":
        n, bins, hist = plt.hist(utr, bin_num, color="#9999FF", label='3\'UTR',
                                 edgecolor='black', linewidth=1)
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
    '''chekc transcript to verify the UTR'''
    detect = False
    for ta in tas:
        if (ta.seq_id == seq_id) and \
           (ta.strand == strand):
            if (ta.start <= utr_start) and \
               (ta.end >= utr_end):
                detect = True
                break
    return detect, ta


def import_utr(tss, utr_strain, utr_all, start, end, tas, length, args_utr):
    if args_utr.source:
        if "Primary" in tss.attributes["type"]:
            if args_utr.base_5utr.lower() == "both":
                detect, ta = check_ta(tas, start, end, tss.seq_id, tss.strand)
            elif args_utr.base_5utr.lower() == "tss":
                detect = True
            if detect:
                utr_strain["pri"][tss.seq_id].append(length)
                utr_all["pri"].append(length)
        elif "Secondary" in tss.attributes["type"]:
            if args_utr.base_5utr.lower() == "both":
                detect, ta = check_ta(tas, start, end, tss.seq_id, tss.strand)
            elif args_utr.base_5utr.lower() == "tss":
                detect = True
            if detect:
                utr_strain["sec"][tss.seq_id].append(length)
                utr_all["sec"].append(length)
    else:
        if args_utr.base_5utr.lower() == "both":
            detect, ta = check_ta(tas, start, end, tss.seq_id, tss.strand)
        elif args_utr.base_5utr.lower() == "tss":
            detect = True
    if detect:
        utr_strain["all"][tss.seq_id].append(length)
        utr_all["all"].append(length)
    if args_utr.base_5utr.lower() == "tss":
        ta = None
    return detect, ta


def get_print_string_5utr(num_utr, name_utr, length, tss, cds_name,
                          locus_tag, ta, source, out, start, end, utrs_tss):
    if "Name" not in tss.attributes.keys():
        tss.attributes["Name"] = (tss.feature + ":" + str(tss.start) + "-" + 
                                  str(tss.end) + "_" + tss.strand)
    attribute_string = ";".join(
             ["=".join(items) for items in [
              ("ID", "_".join([tss.seq_id, "utr5", str(num_utr)])),
              ("Name", "_".join(["5'UTR", name_utr])),
              ("length", str(length)),
              ("associated_cds", cds_name),
              ("associated_gene", locus_tag),
              ("associated_tss", tss.attributes["Name"])]])
    if source:
        attribute_string = ";".join([
            attribute_string, "=".join(["tss_type", tss.attributes["type"]])])
    if ta is not None:
        if "ID" in ta.attributes.keys():
            attribute_string = ";".join([attribute_string, "Parent=" + ta.attributes["ID"]])
        else:
            attribute_string = ";".join([
                attribute_string, "=".join(["Parent", "Transcript:" +
                                            str(ta.start) + "-" + str(ta.end) +
                                            "_" + ta.strand])])
    out.write("{0}\tANNOgesic\t5UTR\t{1}\t{2}\t.\t{3}\t.\t{4}\n".format(
              tss.seq_id, start, end,
              tss.strand, attribute_string))
    utrs_tss.append({"strain": tss.seq_id, "start": start, "end": end,
                     "strand": tss.strand})


def get_5utr(tss, near_cds, utr_strain, utr_all, tas, num_utr,
             num, num_tss, cds_name, locus_tag, out, args_utr,
             utrs_tss, check_cdss, pres):
    '''print and import the 5UTR information'''
    detect = False
    if tss.strand == "+":
        start = tss.start
        end = near_cds.start - 1
        length = end - start + 1
        if str(near_cds.start) + "+" not in check_cdss:
            detect, ta = import_utr(tss, utr_strain, utr_all,
                                    start, end, tas, length, args_utr)
            check_cdss.append(str(near_cds.start) + "+")
    else:
        start = near_cds.end + 1
        end = tss.end
        length = end - start + 1
        if (str(near_cds.end) + "-" not in check_cdss) or [num == num_tss]:
            check_cdss.append(str(near_cds.end) + "-")
            if pres["tss"] is not None:
                detect, ta = import_utr(pres["tss"], utr_strain, utr_all,
                                        pres["start"], pres["end"], tas,
                                        pres["len"], args_utr)
        pres["start"] = start
        pres["end"] = end
        pres["len"] = length
        pres["tss"] = tss
    if detect:
        name_utr = '%0*d' % (5, num_utr)
        if length >= 0:
            get_print_string_5utr(num_utr, name_utr, length, tss, cds_name,
                                  locus_tag, ta, args_utr.source, out,
                                  start, end, utrs_tss)
            num_utr += 1
    return num_utr


def detect_cds(cdss, gene):
    '''get the information of CDS'''
    detect = False
    for cds in cdss:
        if "Parent" in cds.attributes.keys():
            if gene.attributes["ID"] in cds.attributes["Parent"].split(","):
                detect = True
                near_cds = cds
                check_utr = True
                cds_name = get_feature(cds)
        if not detect:
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
    source = False
    gff_f = open(gff_file, "r")
    for entry in Gff3Parser().entries(gff_f):
        if (entry.feature == "CDS"):
            cdss.append(entry)
        elif (entry.feature == "gene"):
            genes.append(entry)
    if ta_file is not None:
        for entry in Gff3Parser().entries(open(ta_file)):
            tas.append(entry)
        tas = sorted(tas, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    if term_file is not None:
        for entry_term in Gff3Parser().entries(open(term_file)):
            terms.append(entry_term)
        terms = sorted(terms, key=lambda k: (k.seq_id, k.start,
                                             k.end, k.strand))
    if tss_file is not None:
        for entry in Gff3Parser().entries(open(tss_file)):
            if "type" in entry.attributes.keys():
                if (entry.attributes["type"] != "Orphan"):
                    source = True
            tsss.append(entry)
        tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    genes = sorted(genes, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    return genes, cdss, terms, tsss, tas, source


def check_associated_TSSpredator(genes, tss, cdss, check_utr, cds_name, locus):
    '''get the associated TSS which is generated from TSSpredator'''
    near_cds = None
    for gene in genes:
        if (tss.seq_id == gene.seq_id) and (
                tss.strand == gene.strand):
            if "locus_tag" in gene.attributes.keys():
                if gene.attributes["locus_tag"] == locus:
                    for cds in cdss:
                        if "Parent" in cds.attributes.keys():
                            if (gene.attributes["ID"] in
                                    cds.attributes["Parent"].split(",")):
                                near_cds = cds
                                check_utr = True
                                cds_name = get_feature(cds)
                        if not check_utr:
                            if "locus_tag" in cds.attributes.keys():
                                if (gene.attributes["locus_tag"] ==
                                        cds.attributes["locus_tag"]):
                                    near_cds = cds
                                    cds_name = cds.attributes["locus_tag"]
                                    check_utr = True
                            elif (gene.seq_id == cds.seq_id) and (
                                    gene.strand == cds.strand):
                                if (gene.start >= cds.start) and (
                                        gene.end <= cds.end):
                                    near_cds = cds
                                    cds_name = get_feature(cds)
                                    check_utr = True
    return check_utr, cds_name, near_cds


def get_5utr_from_TSSpredator(tss, genes, cdss):
    '''It is for TSS file which is generated from ANNOgesic'''
    check_utr = False
    cds_name = "NA"
    if ("Primary" in tss.attributes["type"]) or (
            "Secondary" in tss.attributes["type"]):
        ass_gene = tss.attributes["associated_gene"].split(",")
        tss_type = tss.attributes["type"].split(",")
        for index in range(len(tss_type)):
            if (tss_type[index] == "Primary") or (
                    tss_type[index] == "Secondary"):
                locus_tag = ass_gene[index]
                check_utr, cds_name, near_cds = check_associated_TSSpredator(
                                                genes, tss, cdss, check_utr,
                                                cds_name, locus_tag)
                if not check_utr:
                    for cds in cdss:
                        if (cds.seq_id == tss.seq_id) and (
                                cds.strand == tss.strand):
                            strand = Helper().get_strand_name(cds.strand)
                            if locus_tag == (cds.feature + ":" +
                                             str(cds.start) +
                                             "-" + str(cds.end) +
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


def detect_feature_5utr(feas, tss, cdss, length, check_cds):
    locus_tag = None
    near_cds = None
    cds_name = None
    check_utr = False
    for fea in feas:
        if (tss.seq_id == fea.seq_id) and (
                tss.strand == fea.strand):
            if (tss.strand == "+") and (
                    (fea.start - tss.start) <= length) and (
                    tss.start <= fea.start):
                if "locus_tag" in fea.attributes.keys():
                    locus_tag = fea.attributes["locus_tag"]
                else:
                    locus_tag = fea.attributes["ID"]
                if check_cds:
                    near_cds, cds_name, check_utr = detect_cds(cdss, fea)
                else:
                    near_cds = fea
                    cds_name = get_feature(fea)
                    check_utr = True
                break
            elif (tss.strand == "-") and (
                    (tss.start - fea.end) <= length) and (
                    tss.start >= fea.end):
                if "locus_tag" in fea.attributes.keys():
                    locus_tag = fea.attributes["locus_tag"]
                else:
                    locus_tag = fea.attributes["ID"]
                if check_cds:
                    near_cds, cds_name, check_utr = detect_cds(cdss, fea)
                else:
                    near_cds = fea
                    cds_name = get_feature(fea)
                    check_utr = True
                break
    return near_cds, cds_name, check_utr, locus_tag

def get_5utr_from_other(tss, genes, cdss, length):
    '''It is for TSS file which is not generated from ANNOgesic'''
    check_utr = False
    cds_name = "NA"
    if len(genes) != 0:
        near_cds, cds_name, check_utr, locus_tag = detect_feature_5utr(
                genes, tss, cdss, length, True)
    else:
        near_cds, cds_name, check_utr, locus_tag = detect_feature_5utr(
                cdss, tss, cdss, length, False)
    if check_utr:
        utr_datas = {"check": check_utr, "cds_name": cds_name,
                     "near_cds": near_cds, "locus": locus_tag}
    else:
        utr_datas = {"check": False, "cds_name": None,
                     "near_cds": None, "locus": None}
    return utr_datas


def get_attribute_string(num_utr, length, cds, gene_name, ta, id_name, name,
                         feature, feature_name):
    name_utr = '%0*d' % (5, num_utr)
    cds_name = get_feature(cds)
    attribute_string = ";".join(
                 ["=".join(items) for items in [
                  ("ID", "_".join([ta.seq_id, id_name, str(num_utr)])),
                  ("Name", "_".join([name, name_utr])),
                  ("length", str(length)),
                  ("associated_cds", cds_name),
                  ("associated_gene", gene_name)]])
    if feature == "Parent":
        if "ID" in ta.attributes.keys():
            attribute_string = ";".join([attribute_string,
                                         "Parent=" + ta.attributes["ID"]])
        else:
            attribute_string = ";".join(
                [attribute_string,
                 "Parent=" + (feature, feature_name + str(ta.start) + "-" +
                 str(ta.end) + "_" + ta.strand)])
    elif feature == "associated_term":
            attribute_string = ";".join(
                [attribute_string,
                 "associated_term=" + (feature_name + str(ta.start) + "-" +
                 str(ta.end) + "_" + ta.strand)])
    return attribute_string


def get_gene_name(genes, cds):
    gene_name = "NA"
    for gene in genes:
        if (cds.seq_id == gene.seq_id) and (
                cds.strand == gene.strand):
            if ("Parent" in cds.attributes.keys()) and (
                    "ID" in gene.attributes.keys()):
                if gene.attributes["ID"] in cds.attributes["Parent"].split(","):
                    if "locus_tag" in gene.attributes.keys():
                        gene_name = gene.attributes["locus_tag"]
                    else:
                        gene_name = "Gene:" + str(gene.start) + "-" + \
                                    str(gene.end) + "_" + gene.strand
                    break
            if ((cds.start >= gene.start) and (
                    cds.end <= gene.end)):
                if "locus_tag" in gene.attributes.keys():
                    gene_name = gene.attributes["locus_tag"]
                else:
                    gene_name = ("Gene:" + str(gene.start) + "-" +
                                 str(gene.end) + "_" + gene.strand)
                break
    return gene_name


def set_utr_strain(ta, type_, utr_strain):
    if ta.seq_id not in utr_strain[type_].keys():
        utr_strain[type_][ta.seq_id] = []


def check_repeat(start, end, strain, strand, utrs_tss):
    for utr in utrs_tss:
        if (utr["strain"] == strain) and (
                utr["strand"] == strand):
            if ((utr["start"] >= start) and (
                utr["end"] <= end)) or (
                (utr["start"] <= start) and (
                utr["end"] >= end)) or (
                (utr["start"] >= start) and (
                utr["start"] <= end) and (
                utr["end"] >= end)) or (
                (utr["start"] <= start) and (
                utr["end"] >= start) and (
                utr["end"] <= end)):
                return True
    return False

def compare_ta(tas, genes, cdss, utr_strain, utr_all, out,
               args_utr, utrs_tss, num_utr):
    '''Comparing CDS and trancript to find the 5UTR'''
    for ta in tas:
        detect = False
        set_utr_strain(ta, "all", utr_strain)
        set_utr_strain(ta, "pri", utr_strain)
        set_utr_strain(ta, "sec", utr_strain)
        for cds in cdss:
            if (ta.seq_id == cds.seq_id) and (
                    ta.strand == cds.strand):
                if ta.strand == "+":
                    if ((cds.start - ta.start) <= args_utr.length) and (
                            (cds.start - ta.start) >= 0):
                        if ((ta.start <= cds.start) and (
                                ta.end >= cds.end)) or (
                                (ta.start <= cds.start) and (
                                ta.end <= cds.end) and (
                                ta.end >= cds.start)):
                            if (not check_repeat(ta.start, cds.start, ta.seq_id,
                                                 ta.strand, utrs_tss)):
                                length = cds.start - ta.start
                                utr_strain["all"][ta.seq_id].append(length)
                                utr_all["all"].append(length)
                                gene_name = get_gene_name(genes, cds)
                                string = get_attribute_string(
                                    num_utr, length, cds, gene_name, ta, "utr5",
                                    "5'UTR", "Parent", "Transcript:")
                                detect = True
                                start = ta.start
                                end = cds.start - 1
                                break
                else:
                    if ((ta.end - cds.end) <= args_utr.length) and (
                            (ta.end - cds.end) >= 0):
                        if ((ta.start <= cds.start) and (
                                ta.end >= cds.end)) or (
                                (ta.start >= cds.start) and (
                                ta.start <= cds.end) and (
                                ta.end >= cds.end)):
                            if (not check_repeat(cds.end, ta.end, ta.seq_id,
                                                 ta.strand, utrs_tss)):
                                near_cds = cds
                                detect = True
        if (ta.strand == "-") and (detect):
            length = ta.end - near_cds.end
            utr_strain["all"][ta.seq_id].append(length)
            utr_all["all"].append(length)
            gene_name = get_gene_name(genes, near_cds)
            string = get_attribute_string(
                num_utr, length, near_cds, gene_name,
                ta, "utr5", "5'UTR", "Parent", "Transcript:")
            start = near_cds.end + 1
            end = ta.end
        if detect:
            if end - start > 0:
                string = ";".join([string, "associated_tss=NA", "tss_type=NA"])
                out.write("{0}\tANNOgesic\t5UTR\t{1}"
                          "\t{2}\t.\t{3}\t.\t{4}\n".format(
                              ta.seq_id, start, end, ta.strand, string))
                num_utr += 1
    return num_utr

def detect_5utr(tss_file, gff_file, ta_file, out_file, args_utr):
    '''detection of 5UTR'''
    num_utr = 0
    utr_all = {"all": [], "pri": [], "sec": []}
    utr_strain = {"all": {}, "pri": {}, "sec": {}}
    pre_seq_id = ""
    out = open(out_file, "w")
    out.write("##gff-version 3\n")
    genes, cdss, terms, tsss, tas, source = read_file(
            tss_file, gff_file, ta_file, None)
    utrs_tss = []
    check_cdss = []
    pres = {"check": None, "tss": None, "start": -1, "end": -1, "len": -1}
    if (args_utr.source) and (not source):
        args_utr.source = False
    if (args_utr.base_5utr.upper() == "TSS") or (
            args_utr.base_5utr.lower() == "both"):
        num = 0
        for tss in tsss:
            num = num + 1
            if args_utr.source:
                utr_datas = get_5utr_from_TSSpredator(tss, genes, cdss)
            else:
                utr_datas = get_5utr_from_other(tss, genes, cdss,
                                                args_utr.length)
            if utr_datas["check"]:
                if tss.seq_id != pre_seq_id:
                    pre_seq_id = tss.seq_id
                    utr_strain["pri"][tss.seq_id] = []
                    utr_strain["sec"][tss.seq_id] = []
                    utr_strain["all"][tss.seq_id] = []
                num_utr = get_5utr(tss, utr_datas["near_cds"], utr_strain,
                                   utr_all, tas, num_utr, num, len(tsss),
                                   utr_datas["cds_name"], utr_datas["locus"],
                                   out, args_utr, utrs_tss, check_cdss, pres)
    if (args_utr.base_5utr.lower() == "transcript") or (
            args_utr.base_5utr.lower() == "both"):
        num_utr = compare_ta(tas, genes, cdss, utr_strain,
                             utr_all, out, args_utr, utrs_tss, num_utr)
    out.close()
    Helper().sort_gff(out_file, out_file + "sort")
    shutil.move(out_file + "sort", out_file)
    name = (gff_file.split("/"))[-1].replace(".gff", "")
    plot(utr_all["all"], utr_all["pri"], utr_all["sec"], "_".join([name,
         "all_5utr_length.png"]), args_utr.source, "5utr", args_utr.base_5utr)
    if len(utr_strain["all"]) > 1:
        for strain in utr_strain["all"].keys():
            plot(utr_strain["all"][strain], utr_strain["pri"][strain],
                 utr_strain["sec"][strain],
                 "_".join([strain, "5utr_length.png"]), args_utr.source,
                 "5utr", args_utr.base_5utr)


def compare_term(ta, terms, fuzzy):
    '''Comparing of transcript and terminator to get the 
    terminator which is associated with transcript'''
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


def get_3utr(ta, near_cds, utr_all, utr_strain,
             attributes, num_utr, out, args_utr, utrs_ta):
    '''print the 3UTR'''
    if ta.strand == "+":
        start = near_cds.end + 1
        end = ta.end
        length = ta.end - near_cds.end
        utr_all.append(length)
        utr_strain[ta.seq_id].append(length)
    else:
        start = ta.start
        end = near_cds.start - 1
        length = near_cds.start - ta.start
        utr_all.append(length)
        utr_strain[ta.seq_id].append(length)
    attributes.append("=".join(["length", str(length)]))
    if "ID" not in ta.attributes.keys():
        attributes.append("=".join([
            "Parent", "Transcript:" + str(ta.start) +
            "-" + str(ta.end) + "_" + ta.strand]))
    else:
        attributes.append("=".join([
            "Parent", ta.attributes["ID"]]))
    attribute = ";".join(attributes)
    if (length <= args_utr.length) and (length > 0):
        name_utr = '%0*d' % (5, num_utr)
        name = "=".join(["Name", "_".join(["3'UTR", name_utr])])
        id_ = "ID=" + ta.seq_id + "_utr3_" + str(num_utr)
        num_utr += 1
        attribute_string = ";".join([id_, name, attribute])
        out.write("\t".join([ta.seq_id, "ANNOgesic", "3UTR",
                             str(start), str(end),
                             ta.score, ta.strand, ta.phase,
                             attribute_string]) + "\n")
        utrs_ta.append({"strain": ta.seq_id, "start": start, "end": end,
                        "strand": ta.strand})
    return num_utr

def get_gene_string(gene, attributes):
    if "locus_tag" in gene.attributes.keys():
        attributes.append("=".join(["associated_gene",
                          gene.attributes["locus_tag"]]))
    else:
        gene_string = (gene.feature + ":" + str(gene.start) +
                       "-" + str(gene.end) + "_" + gene.strand)
        attributes.append("=".join(["associated_gene",
                         gene_string]))


def get_near_cds(cdss, genes, ta, attributes, utr_length):
    '''Get the associated CDS of terminator'''
    detect = False
    near_cds = None
    for cds in cdss:
        if (ta.seq_id == cds.seq_id) and (
                ta.strand == cds.strand) and (
                cds.feature == "CDS"):
            if ta.strand == "+":
                if (cds.end <= ta.end) and (
                        (ta.end - cds.end) <= utr_length) and (
                        (ta.end - cds.end) > 0):
                    if (near_cds == None):
                        near_cds = cds
                    else:
                        if (cds.end > near_cds.end):
                            near_cds = cds
                    detect = True
            else:
                if (cds.start >= ta.start) and (
                        (cds.start - ta.start) <= utr_length) and (
                        (cds.start - ta.start) > 0):
                    if (near_cds == None):
                        near_cds = cds
                    else:
                        if (cds.start < near_cds.start):
                            near_cds = cds
                    detect = True
    if detect:
        check_gene = False
        for gene in genes:
            if ("Parent" in near_cds.attributes.keys()):
                if gene.attributes["ID"] in near_cds.attributes["Parent"].split(","):
                    get_gene_string(gene, attributes)
                    check_gene = True
                    break
            else:
                if (gene.seq_id == near_cds.seq_id) and (
                    gene.strand == near_cds.strand):
                    if ((gene.start >= near_cds.start) and (
                        gene.end <= near_cds.end)) or (
                        (gene.start <= near_cds.start) and (
                        gene.end >= near_cds.end)) or (
                        (gene.start >= near_cds.start) and (
                        gene.start <= near_cds.end) and (
                        gene.end >= near_cds.end)) or (
                        (gene.start <= near_cds.start) and (
                        gene.end >= near_cds.start) and (
                        gene.end <= near_cds.end)):
                        get_gene_string(gene, attributes)
                        check_gene = True
                        break
        if not check_gene:
            attributes.append("assoicated_gene=NA")
        if "protein_id" in near_cds.attributes.keys():
            attributes.append("=".join(["associated_cds",
                              near_cds.attributes["protein_id"]]))
        elif "locus_tag" in near_cds.attributes.keys():
            attributes.append("=".join(["associated_cds",
                              near_cds.attributes["locus_tag"]]))
        else:
            cds_string = (near_cds.feature + ":" + str(near_cds.start) + 
                          "-" + str(near_cds.end) + "_" + near_cds.strand)
            attributes.append("=".join(["associated_cds", cds_string]))
    else:
        near_cds = None
    return near_cds


def compare_term_3utr(terms, cdss, genes, utr_all, utr_strain, args_utr,
                      out, utrs_ta, num_utr):
    '''Comparing of terminator and 3UTR to get the relationship'''
    for term in terms:
        detect = False
        if term.seq_id not in utr_strain.keys():
            utr_strain[term.seq_id] = []
        for cds in cdss:
            if (cds.seq_id == term.seq_id) and (
                    cds.strand == term.strand):
                if term.strand == "+":
                    if (term.end >= cds.end) and (
                            ((term.end - cds.end) <= args_utr.length) or (
                            ((term.start - cds.end) <= args_utr.length) and (
                            term.start >= cds.end))):
                        if (not check_repeat(cds.end, term.end, cds.seq_id, 
                            cds.strand, utrs_ta)):
                            detect = True
                            near_cds = cds
                else:
                    if (term.start <= cds.start) and (
                            ((cds.start - term.start) <= args_utr.length) or (
                            ((cds.start - term.end) <= args_utr.length) and (
                            cds.start >= term.end))):
                        if (not check_repeat(term.start, cds.start, cds.seq_id,
                            cds.strand, utrs_ta)):
                            length = term.start - cds.start
                            utr_strain[term.seq_id].append(length)
                            utr_all.append(length)
                            gene_name = get_gene_name(genes, cds)
                            string = get_attribute_string(
                                num_utr, length, cds, gene_name, term, "utr3",
                                "3'UTR", "associated_term", "Terminator:")
                            detect = True
                            start = term.start
                            end = cds.start - 1
                            break
            if (term.strand == "+") and detect:
                length = term.end - near_cds.end
                utr_strain[term.seq_id].append(length)
                utr_all.append(length)
                gene_name = get_gene_name(genes, near_cds)
                string = get_attribute_string(
                    num_utr, length, near_cds, gene_name, term, "utr3",
                    "3'UTR", "associated_term", "Terminator:")
                detect = True
                start = near_cds.end + 1
                end = term.end
        if detect:
            if end - start > 0:
                out.write("{0}\tANNOgesic\t3UTR\t{1}\t"
                          "{2}\t.\t{3}\t.\t{4}\n".format(
                              term.seq_id, start, end, term.strand, string))
                num_utr += 1
    return num_utr


def detect_3utr(ta_file, gff_file, term_file, out_file, args_utr):
    '''For detection of 3UTR'''
    num_utr = 0
    utr_all = []
    utr_strain = {}
    pre_seq_id = ""
    out = open(out_file, "w")
    out.write("##gff-version 3\n")
    genes, cdss, terms, tsss, tas, args_utr.source = read_file(
            None, gff_file, ta_file, term_file)
    utrs_ta = []
    if (args_utr.base_3utr == "transcript") or (
            args_utr.base_3utr == "both"):
        for ta in tas:
            if term_file is not None:
                term = compare_term(ta, terms, args_utr.fuzzy)
            else:
                term = None
            attributes = []
            if term_file is not None:
                if term is not None:
                    attributes.append(
                        "=".join(["associated_term", "Terminator:" +
                                  str(term.start) + "-" +
                                  str(term.end) + "_" + term.strand]))
                else:
                    attributes.append("=".join(["associated_term", "NA"]))
            near_cds = get_near_cds(cdss, genes, ta, attributes, args_utr.length)
            if ta.seq_id != pre_seq_id:
                pre_seq_id = ta.seq_id
                utr_strain[ta.seq_id] = []
            if near_cds is not None:
                num_utr = get_3utr(ta, near_cds, utr_all, utr_strain,
                                   attributes, num_utr, out, args_utr, utrs_ta)
    if (args_utr.base_3utr == "terminator") or (
            args_utr.base_3utr == "both"):
        num_utr = compare_term_3utr(terms, cdss, genes, utr_all, utr_strain,
                                    args_utr, out, utrs_ta, num_utr)
    out.close()
    Helper().sort_gff(out_file, out_file + "sort")
    shutil.move(out_file + "sort", out_file)
    name = (gff_file.split("/"))[-1].replace(".gff", "")
    plot(utr_all, None, None, "_".join([name, "all_3utr_length.png"]),
         None, "3utr", "3utr")
    if len(utr_strain) > 1:
        for strain in utr_strain.keys():
            plot(utr_strain[strain], None, None,
                 "_".join([strain, "3utr_length.png"]), None, "3utr", "3utr")
