from annogesiclib.gff3 import Gff3Parser
from annogesiclib.helper import Helper


def print_file(datas, out, feature):
    for data in datas:
        if feature not in data.attributes.keys():
            data.attributes[feature] = "NA"
        else:
            data.attributes[feature] = ",".join(data.attributes[feature])
        data.attribute_string = ";".join(
            ["=".join(items) for items in data.attributes.items()])
        out.write("\t".join([data.info_without_attributes,
                  data.attribute_string]) + "\n")


def del_attributes(feature, entry):
    '''Remove to the useless attributes'''
    attributes = {}
    for key, value in entry.attributes.items():
        if feature not in key:
            attributes[key] = value
    return attributes


def srna_sorf_comparison(sRNA_file, sORF_file, sRNA_out, sORF_out):
    '''Comparison of sRNA and sORF. It can be a filter of sRNA detection'''
    sorfs = []
    srnas = []
    out_r = open(sRNA_out, "w")
    out_o = open(sORF_out, "w")
    out_r.write("##gff-version 3\n")
    out_o.write("##gff-version 3\n")
    for entry in Gff3Parser().entries(open(sRNA_file)):
        entry.attributes = del_attributes("sORF", entry)
        srnas.append(entry)
    srnas = sorted(srnas, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    for entry in Gff3Parser().entries(open(sORF_file)):
        entry.attributes = del_attributes("sRNA", entry)
        sorfs.append(entry)
    sorfs = sorted(sorfs, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    for srna in srnas:
        for sorf in sorfs:
            if (srna.seq_id == sorf.seq_id) and (srna.strand == sorf.strand):
                if ((srna.start <= sorf.start) and (
                        srna.end >= sorf.end)) or (
                        (srna.start >= sorf.start) and (
                         srna.end <= sorf.end)) or (
                        (srna.start <= sorf.start) and (
                         srna.end >= sorf.start) and (
                         srna.end <= sorf.end)) or (
                        (srna.start >= sorf.start) and (
                         srna.start <= sorf.end) and (
                         srna.end >= sorf.end)):
                    if "sORF" not in srna.attributes.keys():
                        srna.attributes["sORF"] = []
                        strand = Helper().get_strand_name(sorf.strand)
                    srna.attributes["sORF"].append("".join([
                                               "sORF:",
                                               str(sorf.start), "-",
                                               str(sorf.end),
                                               "_", strand]))
                    if "sRNA" not in sorf.attributes.keys():
                        sorf.attributes["sRNA"] = []
                        strand = Helper().get_strand_name(srna.strand)
                    sorf.attributes["sRNA"].append("".join([
                                               "sRNA:",
                                               str(srna.start), "-",
                                               str(srna.end),
                                               "_", strand]))
    print_file(sorfs, out_o, "sRNA")
    print_file(srnas, out_r, "sORF")
    out_r.close()
    out_o.close()
