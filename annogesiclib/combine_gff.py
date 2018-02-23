from annogesiclib.gff3 import Gff3Parser


def read_file(filename):
    datas = []
    f_h = open(filename, "r")
    for entry in Gff3Parser().entries(f_h):
        datas.append(entry)
    datas = sorted(datas, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    f_h.close()
    return datas


def del_attributes(entry):
    '''delete the feature (besides genome annotation) 
    which already has Parent'''
    if (entry.feature == "CDS") or (
            entry.feature == "tRNA") or (
            entry.feature == "rRNA") or (
            entry.feature == "exon"):
        pass
    else:
        if "Parent" in entry.attributes.keys():
            del entry.attributes["Parent"]


def print_file(entry, tran, out):
    if "Parent" in entry.attributes.keys():
        if str(tran.attributes["ID"]) not in entry.attributes["Parent"]:
            entry.attributes["Parent"] = ",".join([
                entry.attributes["Parent"], str(tran.attributes["ID"])])
    else:
        entry.attributes["Parent"] = str(tran.attributes["ID"])
    attributes = {}
    for key, value in entry.attributes.items():
        if (key != "print") and (key != "parent_tran"):
            attributes[key] = value
    attribute_string = ";".join(
        ["=".join(items) for items in attributes.items()])
    out.write("".join([entry.info_without_attributes, "\t",
                       attribute_string, "\n"]))
    entry.attributes["print"] = True


def compare_tran(datas, tran, out):
    '''comare transcript and 5UTR, 3UTR, gene/CDS for merging'''
    for data in datas:
        del_attributes(data)
        if (data.seq_id == tran.seq_id) and (
                data.strand == tran.strand):
            if (data.start >= tran.start) and (
                    data.end <= tran.end):
                print_file(data, tran, out)


def print_rest(datas, out):
    '''print the rest data which is not related to operon'''
    for data in datas:
        if "print" not in data.attributes.keys():
            out.write(data.info + "\n")


def compare_tran_term(term, tran, out, fuzzy_term):
    '''compare transcript and terminator for merging'''
    if (term.seq_id == tran.seq_id) and (
            term.strand == tran.strand):
        if (term.start >= tran.start) and (
                term.end <= tran.end):
            print_file(term, tran, out)
        else:
            if term.strand == "+":
                if ((term.start - fuzzy_term) <= tran.end) and (
                        term.end + fuzzy_term >= tran.end):
                    print_file(term, tran, out)
                elif (term.start <= tran.end) and (
                        term.end >= tran.end):
                    print_file(term, tran, out)
            else:
                if (term.end + fuzzy_term >= tran.start) and (
                        term.start - fuzzy_term <= tran.start):
                    print_file(term, tran, out)
                elif (term.start <= tran.start) and (
                        term.end >= tran.start):
                    print_file(term, tran, out)


def combine_gff(gff_file, ta_file, tss_file, utr5_file, utr3_file,
                term_file, fuzzy_tss, fuzzy_term, out_file):
    '''combine the features which related to operon to 
    form a operon gff file'''
    gffs = read_file(gff_file)
    trans = read_file(ta_file)
    if tss_file is not None:
        tsss = read_file(tss_file)
    if utr5_file is not None:
        utr5s = read_file(utr5_file)
    if utr3_file is not None:
        utr3s = read_file(utr3_file)
    out = open(out_file, "w")
    out.write("##gff-version 3\n")
    if term_file is not None:
        terms = read_file(term_file)
    for tran in trans:
        out.write(tran.info + "\n")
        if tss_file is not None:
            for tss in tsss:
                del_attributes(tss)
                if (tss.seq_id == tran.seq_id) and (tss.strand == tran.strand):
                    if tss.strand == "+":
                        if ((tss.start + fuzzy_tss) >= tran.start) and (
                                tss.start <= tran.end):
                            print_file(tss, tran, out)
                    else:
                        if (tss.start >= tran.start) and (
                                tss.end - fuzzy_tss <= tran.end):
                            print_file(tss, tran, out)
        if utr5_file is not None:
            compare_tran(utr5s, tran, out)
        compare_tran(gffs, tran, out)
        if utr3_file is not None:
            compare_tran(utr3s, tran, out)
        if term_file is not None:
            for term in terms:
                del_attributes(term)
                compare_tran_term(term, tran, out, fuzzy_term)
    if tss_file is not None:
        print_rest(tsss, out)
    if utr5_file is not None:
        print_rest(utr5s, out)
    print_rest(gffs, out)
    if utr3_file is not None:
        print_rest(utr3s, out)
    if term_file is not None:
        print_rest(terms, out)
    out.close()
