from annogesiclib.gff3 import Gff3Parser
from annogesiclib.helper import Helper


def get_feature(cds):
    if "locus_tag" in cds.attributes.keys():
        feature = cds.attributes["locus_tag"]
    elif "protein_id" in cds.attributes.keys():
        feature = cds.attributes["protein_id"]
    else:
        feature = cds.attributes["ID"]
    return feature


def import_data(seq, cds, start, end):
    feature = get_feature(cds)
    return {"seq": seq, "strain": cds.seq_id, "strand": cds.strand,
            "protein": feature, "start": start, "end": end}


def detect_site(inters, args_ribo):
    '''Detection of ribosome binding site'''
    rbss = []
    for inter in inters:
        if args_ribo.without_rbs:
            rbss.append(inter)
        else:
            for ribo_seq in args_ribo.rbs_seq:
                detect = False
                for nts in range(0, (len(inter["seq"]) - len(ribo_seq))):
                    miss = 0
                    for index in range(len(ribo_seq)):
                        if miss > args_ribo.fuzzy_rbs:
                            break
                        else:
                            if inter["seq"][nts:(nts + len(ribo_seq))][index] != ribo_seq[index]:
                                miss += 1
                    if (miss <= args_ribo.fuzzy_rbs) and (
                            len(inter["seq"][nts:(nts + len(ribo_seq))]) >= len(ribo_seq)):
                        rbss.append(inter)
                        detect = True
                        break
                if detect:
                    break
    return rbss


def read_file(seq_file, gff_file, tss_file, tran_file):
    cdss = []
    tsss = []
    trans = []
    seq = {}
    with open(seq_file, "r") as f_h:
        for line in f_h:
            line = line.strip()
            if line.startswith(">"):
                strain = line[1:]
                seq[strain] = ""
            else:
                seq[strain] = seq[strain] + line
    g_h = open(gff_file)
    for entry in Gff3Parser().entries(g_h):
        if (entry.feature == "CDS"):
            cdss.append(entry)
    if tss_file is not None:
        t_h = open(tss_file)
        for entry in Gff3Parser().entries(t_h):
            tsss.append(entry)
        tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
        t_h.close()
    a_h = open(tran_file)
    for entry in Gff3Parser().entries(a_h):
        trans.append(entry)
    cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    trans = sorted(trans, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    g_h.close()
    a_h.close()
    return cdss, seq, tsss, trans


def extract_inter_seq(inter, cds, seq, fuzzy, inters):
    '''extract the sequence of pre-CDS region'''
    helper = Helper()
    start = inter["start"] - fuzzy
    end = inter["end"] + fuzzy
    if inter["start"] - fuzzy <= 0:
        start = 1
    if inter["end"] + fuzzy >= len(seq[cds.seq_id]):
        end = len(seq)
    if cds.strand == "+":
        inter_seq = helper.extract_gene(seq[cds.seq_id], start,
                                        end, "+")
    else:
        inter_seq = helper.extract_gene(seq[cds.seq_id], start,
                                        end, "-")
    inters.append(import_data(inter_seq, cds, inter["start"], inter["end"]))


def compare_tss(tsss, cds, inters, fuzzy, seq, utr):
    '''Compare CDS with TSS to get the pre-CDS region'''
    for tss in tsss:
        if (cds.seq_id == tss.seq_id) and (
                cds.strand == tss.strand):
            if tss.strand == "+":
                if ((cds.start - tss.start) <= utr) and (
                        (cds.start - tss.start) > 0):
                    inter = {"start": tss.start, "end": cds.start,
                             "strain": cds.seq_id, "strand": cds.strand}
                    extract_inter_seq(inter, cds, seq, fuzzy, inters)
            else:
                if ((tss.end - cds.end) <= utr) and (
                        (tss.end - cds.end) > 0):
                    inter = {"start": cds.end, "end": tss.end,
                             "strain": cds.seq_id, "strand": cds.strand}
                    extract_inter_seq(inter, cds, seq, fuzzy, inters)


def compare_pre_cds(first, cdss, cds, seq):
    '''Search the front position CDS of the query one 
    to get the pre-CDS region'''
    detect_cds = False
    start = None
    end = None
    for pre_cds in cdss:
        if (pre_cds.seq_id == cds.seq_id) and (
                pre_cds.strand == cds.strand):
            if pre_cds.strand == "+":
                if first:
                    start = 1
                    end = cds.start
                    detect_cds = True
                    first = False
                    break
                elif pre_cds.end < cds.start:
                    start = pre_cds.end
                    end = cds.start
                    detect_cds = True
                elif pre_cds.end >= cds.start:
                    break
            else:
                if pre_cds.start > cds.end:
                    start = cds.end
                    end = pre_cds.start
                    detect_cds = True
                    break
    if (not detect_cds) and (cds.strand == "-"):
        start = cds.end
        end = len(seq[cds.seq_id])
    return first, start, end


def compare_tran(cds, trans, seq, inters, fuzzy, start, end):
    '''For detect the expressed region of candidates'''
    detect = False
    for tran in trans:
        if (tran.seq_id == cds.seq_id) and (
                tran.strand == cds.strand):
            if tran.strand == "+":
                if (cds.start > tran.start) and (
                        cds.start <= tran.end):
                    if start < tran.start:
                        start = tran.start
                        detect = True
                elif (tran.start > cds.start):
                    break
            else:
                if (cds.end > tran.start) and (
                        cds.end <= tran.end):
                    if end > tran.start:
                        end = tran.start
                        detect = True
                elif (tran.start > cds.end):
                    break
    if detect:
        inter = {"start": tran.start, "end": cds.start,
                 "strain": cds.seq_id, "strand": cds.strand}
        extract_inter_seq(inter, cds, seq, fuzzy, inters)
    else:
        for tran in trans:
            if (tran.seq_id == cds.seq_id) and (
                    tran.strand == cds.strand):
                if ((start <= tran.start) and (
                         end >= tran.end)) or (
                        (start >= tran.start) and (
                         end <= tran.end)) or (
                        (start <= tran.start) and (
                         end <= tran.end) and (
                         end >= tran.start)) or (
                        (start >= tran.start) and (
                         start <= tran.end) and (
                         end >= tran.end)):
                    inter = {"start": tran.start, "end": cds.start,
                             "strain": cds.seq_id, "strand": cds.strand}
                    extract_inter_seq(inter, cds, seq, fuzzy, inters)
                    break


def extract_seq(cdss, seq, tsss, trans, fuzzy, utr):
    '''extract the sequence for searching the riboswitch or RNA thermometer
    by comparing with TSS, transcript and CDS'''
    first = True
    inters = []
    for cds in cdss:
        if len(tsss) != 0:
            compare_tss(tsss, cds, inters, fuzzy, seq, utr)
        first, start, end = compare_pre_cds(first, cdss, cds, seq)
        if (start is not None) and (end is not None):
            compare_tran(cds, trans, seq, inters, fuzzy, start, end)
    return inters


def extract_potential_rbs(seq_file, gff_file, tss_file, tran_file,
                          out_file, args_ribo, feature):
    '''Get the potential riboswitch or RNA-thermometer'''
    out = open(out_file, "w")
    cdss, seq, tsss, trans = read_file(seq_file, gff_file, tss_file, tran_file)
    inters = extract_seq(cdss, seq, tsss, trans,
                         args_ribo.fuzzy, args_ribo.utr)
    rbss = detect_site(inters, args_ribo)
    num = 0
    for rbs in rbss:
        out.write(">" + feature + "_{0}\n".format(
                  "|".join([str(num), rbs["strain"], rbs["strand"],
                            rbs["protein"], str(rbs["start"]),
                            str(rbs["end"])])))
        out.write(rbs["seq"] + "\n")
        num += 1
    out.close()
