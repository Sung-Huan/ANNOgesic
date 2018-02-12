from annogesiclib.gff3 import Gff3Parser
from annogesiclib.helper import Helper


def uni(tas, genes, out):
    '''This is for the transcript which is not overlap with annotation'''
    start_tmp = 0
    stop_tmp = 0
    strand_tmp = ""
    detect = False
    for ta in tas:
        for gene in genes:
            if (ta.strand == gene.strand) and (
                    ta.seq_id == gene.seq_id):
                if ((ta.start < gene.start) and (
                         ta.end > gene.start) and (
                         ta.end < gene.end)) or (
                        (ta.start > gene.start) and (
                         ta.end < gene.end)) or (
                        (ta.start > gene.start) and (
                         ta.start < gene.end) and (
                         ta.end > gene.end)):
                    detect = True
        if ((not detect) and (start_tmp != ta.start) and (
                stop_tmp != ta.end)) or (
                (not detect) and (((start_tmp == ta.start) or (
                    stop_tmp == ta.end)) and (
                    strand_tmp != ta.strand))):
            out.write(ta.info + "\n")
            start_tmp = ta.start
            stop_tmp = ta.end
            strand_tmp = ta.strand
        detect = False


def check_modify(start_tmp, stop_tmp, gene, ta, modify):
    if (ta.start >= gene.start) and (
            ta.end <= gene.end):
        if "within_extend_ends" in modify:
            start_tmp = gene.start
            stop_tmp = gene.end
        else:
            start_tmp = ta.start
            stop_tmp = ta.end
    elif ((ta.start <= gene.start) and (
            ta.end >= gene.start) and (
            ta.end <= gene.end)):
        if (ta.strand == "+") and ("extend_3end" in modify):
            start_tmp = ta.start
            stop_tmp = gene.end
        elif (ta.strand == "-") and ("extend_5end" in modify):
            start_tmp = ta.start
            stop_tmp = gene.end
        else:
            start_tmp = ta.start
            stop_tmp = ta.end
    elif ((ta.start >= gene.start) and (
            ta.start <= gene.end) and (
            ta.end >= gene.end)):
        if (ta.strand == "+") and ("extend_5end" in modify):
            start_tmp = gene.start
            stop_tmp = ta.end
        elif (ta.strand == "-") and ("extend_3end" in modify):
            start_tmp = gene.start
            stop_tmp = ta.end
        else:
            start_tmp = ta.start
            stop_tmp = ta.end
    return start_tmp, stop_tmp


def overlap(tas, genes, out, modify):
    '''Check the overlap of annotation and transcript'''
    check = False
    for gene in genes:
        start_tmp = 0
        stop_tmp = 0
        start = 0
        stop = 0
        for ta in tas:
            if (ta.strand == gene.strand) and (
                    ta.seq_id == gene.seq_id):
                if ((ta.start <= gene.start) and (
                         ta.end >= gene.start) and (
                         ta.end <= gene.end)) or (
                        (ta.start >= gene.start) and (
                         ta.end <= gene.end)) or (
                        (ta.start >= gene.start) and (
                         ta.start <= gene.end) and (
                         ta.end >= gene.end)):
                    check = True
                    tmp_ta = ta
                    if start_tmp == 0:
                        start_tmp, stop_tmp = check_modify(
                            start_tmp, stop_tmp, gene, ta, modify)
                        start = start_tmp
                        stop = stop_tmp
                    else:
                        start_tmp, stop_tmp = check_modify(
                            start_tmp, stop_tmp, gene, ta, modify)
                        if "merge_overlap" in modify:
                            if stop < stop_tmp:
                                stop = stop_tmp
                        else:
                            if (start_tmp != 0):
                                out.write("\t".join([str(field) for field in [
                                    tmp_ta.seq_id, tmp_ta.source,
                                    tmp_ta.feature,
                                    start, stop, tmp_ta.score,
                                    tmp_ta.strand, tmp_ta.phase,
                                    tmp_ta.attribute_string]]) + "\n")
                                start = start_tmp
                                stop = stop_tmp
                if (ta.start > gene.end) and (start != 0) and (check):
                    check = False
                    out.write("\t".join([str(field) for field in [
                              tmp_ta.seq_id, tmp_ta.source, tmp_ta.feature,
                              start, stop, tmp_ta.score,
                              tmp_ta.strand, tmp_ta.phase,
                              tmp_ta.attribute_string]]) + "\n")
                    break
        if (start != 0) and (check):
            out.write('\t'.join([str(field) for field in [
                      tmp_ta.seq_id, tmp_ta.source, tmp_ta.feature,
                      start, stop, tmp_ta.score, tmp_ta.strand,
                      tmp_ta.phase, tmp_ta.attribute_string]]) + "\n")


def fill_gap(gff_file, ta_file, type_, output, modify):
    '''compare transcript with genome annotation to modify the transcript'''
    tas = []
    genes = []
    ta_f = open(ta_file, "r")
    gff_f = open(gff_file, "r")
    for entry in Gff3Parser().entries(ta_f):
        tas.append(entry)
    ta_f.close()
    tas = sorted(tas, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    for entry in Gff3Parser().entries(gff_f):
        if (entry.feature != "source") and (
                entry.feature != "region") and (
                entry.feature != "repeat_region") and (
                entry.feature != "STS") and (
                entry.feature != "remark"):
            genes.append(entry)
    gff_f.close()
    genes = sorted(genes, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    out = open(output, "w")
    out.write("##gff-version 3\n")
    if type_ == "overlap":
        overlap(tas, genes, out, modify)
    elif type_ == "uni":
        uni(tas, genes, out)


def print_file(ta, num, out):
    ta.attributes["ID"] = ta.seq_id + "_transcript" + str(num)
    ta.attributes["Name"] = "transcript_" + ('%0*d' % (5, num))
    attribute_string = ";".join(
        ["=".join(items) for items in ta.attributes.items()])
    out.write("\t".join([str(field) for field in [
              ta.seq_id, ta.source, ta.feature, ta.start,
              ta.end, ta.score, ta.strand, ta.phase,
              attribute_string]]) + "\n")


def longer_ta(ta_file, length, out_file):
    '''merge overlaped transcript to for a complete transcript'''
    tas = []
    for entry in Gff3Parser().entries(open(ta_file)):
        tas.append(entry)
    tas = sorted(tas, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    for ta_1 in tas:
        for ta_2 in tas:
            if (ta_1.seq_id == ta_2.seq_id) and (
                    ta_1.strand == ta_2.strand):
                if (ta_1.start <= ta_2.start) and (
                        ta_1.end >= ta_2.start) and (
                        ta_1.end <= ta_2.end):
                    ta_1.end = ta_2.end
                elif (ta_1.start >= ta_2.start) and (
                        ta_1.start <= ta_2.end) and (
                        ta_1.end >= ta_2.end):
                    ta_1.start = ta_2.start
                elif (ta_1.start <= ta_2.start) and (
                        ta_1.end >= ta_2.end):
                    pass
                elif (ta_1.start >= ta_2.start) and (
                        ta_1.end <= ta_2.end):
                    ta_1.start = ta_2.start
                    ta_1.end = ta_2.end
    first = True
    out = open(out_file, "w")
    out.write("##gff-version 3\n")
    num = 0
    pre_ta = None
    tas = sorted(tas, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    for ta in tas:
        if (ta.end - ta.start) >= length:
            if first:
                first = False
                print_file(ta, num, out)
                num += 1
            else:
                if (ta.seq_id == pre_ta.seq_id) and (
                        ta.strand == pre_ta.strand) and (
                        ta.start == pre_ta.start) and (
                        ta.end == pre_ta.end):
                    pass
                else:
                    print_file(ta, num, out)
                    num += 1
        pre_ta = ta
    out.close()
