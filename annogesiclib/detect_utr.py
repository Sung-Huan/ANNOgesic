import math
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
                                      color="#FF9999", label='secondary')
            n, bins, hist2 = plt.hist(utr_pri, bin_num,
                                      color="#9999FF", label='primary')
            plt.legend((hist2[0], hist1[0]),
                       ("Primary TSSs", "Secondary TSSs"))
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
                          locus_tag, ta, source, out, start, end):
    attribute_string = ";".join(
             ["=".join(items) for items in [
              ("ID", "_".join(["utr5", str(num_utr)])),
              ("Name", "_".join(["5'UTR", name_utr])),
              ("length", str(length)),
              ("associated_cds", cds_name),
              ("associated_gene", locus_tag),
              ("associated_tss", tss.attributes["Name"])]])
    if source:
        attribute_string = ";".join([
            attribute_string, "=".join(["tss_type", tss.attributes["type"]])])
    if ta is not None:
        attribute_string = ";".join([
            attribute_string, "=".join(["Parent", "Transcript:" +
                                        str(ta.start) + "-" + str(ta.end) +
                                        "_" + ta.strand])])
    out.write("{0}\tANNOgesic\t5UTR\t{1}\t{2}\t.\t{3}\t.\t{4}\n".format(
              tss.seq_id, start, end,
              tss.strand, attribute_string))


def get_5utr(tss, near_cds, utr_strain, utr_all, tas, num_utr,
             cds_name, locus_tag, out, args_utr):
    '''print and import the 5UTR information'''
    if tss.strand == "+":
        start = tss.start
        end = near_cds.start
        length = end - start
        detect, ta = import_utr(tss, utr_strain, utr_all,
                                start, end, tas, length, args_utr)
    else:
        start = near_cds.end
        end = tss.end
        length = end - start
        detect, ta = import_utr(tss, utr_strain, utr_all,
                                start, end, tas, length, args_utr)
    if detect:
        name_utr = '%0*d' % (5, num_utr)
        get_print_string_5utr(num_utr, name_utr, length, tss, cds_name,
                              locus_tag, ta, args_utr.source, out, start, end)
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
            tsss.append(entry)
        tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    genes = sorted(genes, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    return genes, cdss, terms, tsss, tas


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
                        else:
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


def get_5utr_from_other(tss, genes, cdss, length):
    '''It is for TSS file which is not generated from ANNOgesic'''
    check_utr = False
    cds_name = "NA"
    for gene in genes:
        if (tss.seq_id == gene.seq_id) and (
                tss.strand == gene.strand):
            if (tss.strand == "+") and (
                    (gene.start - tss.start) <= length) and (
                    tss.start <= gene.start):
                if "locus_tag" in gene.attributes.keys():
                    locus_tag = gene.attributes["locus_tag"]
                else:
                    locus_tag = gene.attributes["ID"]
                near_cds, cds_name, check_utr = detect_cds(cdss, gene)
                break
            elif (tss.strand == "-") and (
                    (tss.start - gene.end) <= length) and (
                    tss.start >= gene.end):
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


def get_attribute_string(num_utr, length, cds, gene_name, ta, id_name, name,
                         feature, feature_name):
    name_utr = '%0*d' % (5, num_utr)
    cds_name = get_feature(cds)
    attribute_string = ";".join(
                 ["=".join(items) for items in [
                  ("ID", "_".join([id_name, str(num_utr)])),
                  ("Name", "_".join([name, name_utr])),
                  ("length", str(length)),
                  ("associated_cds", cds_name),
                  ("associated_gene", gene_name),
                  (feature, feature_name + str(ta.start) + "-" +
                   str(ta.end) + "_" + ta.strand)]])
    return attribute_string


def get_gene_name(genes, cds):
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


def compare_ta(tas, genes, cdss, utr_strain, utr_all, out, args_utr):
    '''Comparing CDS and trancript to find the 5UTR'''
    num_utr = 0
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
                            length = cds.start - ta.start
                            utr_strain["all"][ta.seq_id].append(length)
                            utr_all["all"].append(length)
                            gene_name = get_gene_name(genes, cds)
                            string = get_attribute_string(
                                num_utr, length, cds, gene_name, ta, "utr5",
                                "5'UTR", "Parent", "Transcript:")
                            detect = True
                            start = ta.start
                            end = cds.start
                            break
                else:
                    if ((ta.end - cds.end) <= args_utr.length) and (
                            (ta.end - cds.end) >= 0):
                        if ((ta.start <= cds.start) and (
                                ta.end >= cds.end)) or (
                                (ta.start >= cds.start) and (
                                ta.start <= cds.end) and (
                                ta.end >= cds.end)):
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
            start = near_cds.end
            end = ta.end
        if detect:
            if end - start > 0:
                out.write("{0}\tANNOgesic\t5UTR\t{1}"
                          "\t{2}\t.\t{3}\t.\t{4}\n".format(
                              ta.seq_id, start, end, ta.strand, string))
                num_utr += 1


def detect_5utr(tss_file, gff_file, ta_file, out_file, args_utr):
    '''detection of 5UTR'''
    num_utr = 0
    utr_all = {"all": [], "pri": [], "sec": []}
    utr_strain = {"all": {}, "pri": {}, "sec": {}}
    pre_seq_id = ""
    out = open(out_file, "w")
    out.write("##gff-version 3\n")
    genes, cdss, terms, tsss, tas = read_file(tss_file, gff_file,
                                              ta_file, None)
    if (args_utr.base_5utr.upper() == "TSS") or (
            args_utr.base_5utr.lower() == "both"):
        for tss in tsss:
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
                                   utr_all, tas, num_utr,
                                   utr_datas["cds_name"], utr_datas["locus"],
                                   out, args_utr)
    if args_utr.base_5utr.lower() == "transcript":
        compare_ta(tas, genes, cdss, utr_strain,
                   utr_all, out, args_utr)
    name = (gff_file.split("/"))[-1].replace(".gff", "")
    plot(utr_all["all"], utr_all["pri"], utr_all["sec"], "_".join([name,
         "all_5utr_length.png"]), args_utr.source, "5utr", args_utr.base_5utr)
    if len(utr_strain["all"]) > 1:
        for strain in utr_strain["all"].keys():
            plot(utr_strain["all"][strain], utr_strain["pri"][strain],
                 utr_strain["sec"][strain],
                 "_".join([strain, "5utr_length.png"]), args_utr.source,
                 "5utr", args_utr.base_5utr)
    out.close()


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
             attributes, num_utr, out, args_utr):
    '''print the 3UTR'''
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
    attributes.append("=".join([
        "Parent", "Transcript:" + str(ta.start) +
        "-" + str(ta.end) + "_" + ta.strand]))
    attribute = ";".join(attributes)
    if (length <= args_utr.length) and (length > 0):
        name_utr = '%0*d' % (5, num_utr)
        name = "=".join(["Name", "_".join(["3'UTR", name_utr])])
        id_ = "ID=utr3_" + str(num_utr)
        num_utr += 1
        attribute_string = ";".join([id_, name, attribute])
        if args_utr.base_3utr == "transcript":
            out.write("\t".join([ta.seq_id, "ANNOgesic", "3UTR",
                                 str(start), str(end),
                                 ta.score, ta.strand, ta.phase,
                                 attribute_string]) + "\n")
        elif args_utr.base_3utr == "both":
            if "associated_term=NA" not in attributes:
                out.write("\t".join([ta.seq_id, "ANNOgesic", "3UTR",
                                     str(start), str(end),
                                     ta.score, ta.strand, ta.phase,
                                     attribute_string]) + "\n")
    return num_utr


def get_near_cds(cdss, genes, ta, attributes):
    '''Get the associated CDS of terminator'''
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
                if gene.attributes["ID"] in near_cds.attributes["Parent"].split(","):
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


def compare_term_3utr(terms, cdss, genes, utr_all, utr_strain, args_utr, out):
    '''Comparing of terminator and 3UTR to get the relationship'''
    num_utr = 0
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
                        detect = True
                        near_cds = cds
                else:
                    if (term.start <= cds.start) and (
                            ((cds.start - term.start) <= args_utr.length) or (
                            ((cds.start - term.end) <= args_utr.length) and (
                            cds.start >= term.end))):
                        length = term.start - cds.start
                        utr_strain[term.seq_id].append(length)
                        utr_all.append(length)
                        gene_name = get_gene_name(genes, cds)
                        string = get_attribute_string(
                            num_utr, length, cds, gene_name, term, "utr3",
                            "3'UTR", "associated_term", "Terminator:")
                        detect = True
                        start = term.start
                        end = cds.start
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
                start = near_cds.end
                end = term.end
        if detect:
            if end - start > 0:
                out.write("{0}\tANNOgesic\t5UTR\t{1}\t"
                          "{2}\t.\t{3}\t.\t{4}\n".format(
                              term.seq_id, start, end, term.strand, string))
                num_utr += 1


def detect_3utr(ta_file, gff_file, term_file, out_file, args_utr):
    '''For detection of 3UTR'''
    num_utr = 0
    utr_all = []
    utr_strain = {}
    pre_seq_id = ""
    out = open(out_file, "w")
    out.write("##gff-version 3\n")
    genes, cdss, terms, tsss, tas = read_file(None, gff_file,
                                              ta_file, term_file)
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
            near_cds = get_near_cds(cdss, genes, ta, attributes)
            if ta.seq_id != pre_seq_id:
                pre_seq_id = ta.seq_id
                utr_strain[ta.seq_id] = []
            num_utr = get_3utr(ta, near_cds, utr_all, utr_strain,
                               attributes, num_utr, out, args_utr)
    else:
        compare_term_3utr(terms, cdss, genes, utr_all,
                          utr_strain, args_utr, out)
    name = (gff_file.split("/"))[-1].replace(".gff", "")
    plot(utr_all, None, None, "_".join([name, "all_3utr_length.png"]),
         None, "3utr", "3utr")
    out.close()
    if len(utr_strain) > 1:
        for strain in utr_strain.keys():
            plot(utr_strain[strain], None, None,
                 "_".join([strain, "3utr_length.png"]), None, "3utr", "3utr")
