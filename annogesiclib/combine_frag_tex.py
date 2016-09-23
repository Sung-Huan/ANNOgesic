from annogesiclib.gff3 import Gff3Parser


def modify_position(frag, norm):
    '''get proper position, we choose the long one'''
    if frag.end < norm.end:
        frag.end = norm.end
    if frag.start > norm.start:
        frag.start = norm.start
    norm.attributes["print"] = True
    frag.attributes["print"] = True


def print_file(data, out, name, num):
    attributes = {}
    attributes["ID"] = "tran" + str(num)
    attributes["Name"] = "Tran_" + name
    attributes["detect_lib"] = data.attributes["detect_lib"]
    attribute_string = ";".join(["=".join(items)
                                 for items in attributes.items()])
    out.write("\t".join([str(field) for field in [
                        data.seq_id, data.source, data.feature, data.start,
                        data.end, data.score, data.strand, data.phase,
                        attribute_string]]) + "\n")


def store(data, source, finals):
    data.attributes["detect_lib"] = source
    data.attributes["print"] = False
    finals.append(data)


def compare(data1, data2, overlap, tolerance):
    '''search the sRNA which can be detected in frag and tex libs.
    Then, try to merge them to be a longer one'''
    if (data1.seq_id == data2.seq_id) and (data1.strand == data2.strand):
        if (data1.start <= (data2.end + tolerance)) and (
                data1.start >= data2.start):
            modify_position(data1, data2)
            overlap = True
        elif (data1.end >= (data2.start - tolerance)) and (
                data1.end <= data2.end):
            modify_position(data1, data2)
            overlap = True
        elif (data1.start <= data2.start) and (
                data1.end >= data2.end):
            modify_position(data1, data2)
            overlap = True
    return overlap


def combine(frag_file, tex_file, tolerance, output_file):
    '''merge the results of sRNA which detected by fragmented and dRNA'''
    frags = []
    norms = []
    finals = []
    out = open(output_file, "w")
    out.write("##gff-version 3\n")
    f_h = open(frag_file, "r")
    for entry in Gff3Parser().entries(f_h):
        entry.attributes["print"] = False
        frags.append(entry)
    f_h.close()
    n_h = open(tex_file, "r")
    for entry in Gff3Parser().entries(n_h):
        entry.attributes["print"] = False
        norms.append(entry)
    n_h.close()
    sort_frags = sorted(frags, key=lambda k: (k.seq_id, k.start,
                                              k.end, k.strand))
    sort_norms = sorted(norms, key=lambda k: (k.seq_id, k.start,
                                              k.end, k.strand))
    for frag in sort_frags:
        overlap = False
        for norm in sort_norms:
            overlap = compare(frag, norm, overlap, tolerance)
        if overlap:
            store(frag, "fragmented,tex_notex", finals)
        else:
            store(frag, "fragmented", finals)
    for norm in sort_norms:
        if norm.attributes["print"] is False:
            store(norm, "tex_notex", finals)
    sort_finals = sorted(finals, key=lambda k: (k.seq_id, k.start,
                                                k.end, k.strand))
    num = 0
    for tar in sort_finals:
        if tar.attributes["print"] is True:
            continue
        overlap = False
        for ref in sort_finals:
            overlap = compare(tar, ref, overlap, tolerance)
        name = '%0*d' % (5, num)
        print_file(tar, out, name, num)
        num += 1
    out.close()
