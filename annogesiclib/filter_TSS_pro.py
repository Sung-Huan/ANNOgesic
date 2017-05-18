import os
import shutil
import math
from annogesiclib.gff3 import Gff3Parser


def read_gff(input_file):
    datas = []
    gff_parser = Gff3Parser()
    f_h = open(input_file, "r")
    for entry in gff_parser.entries(f_h):
        datas.append(entry)
    datas = sorted(datas, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    return datas


def compare_tss_pro(tars, refs, out, cluster):
    '''compare between TSS and processing site'''
    for tar in tars:
        for ref in refs:
            if (tar.seq_id == ref.seq_id) and (
                    tar.strand == ref.strand):
                if math.fabs(tar.start - ref.start) <= cluster:
                    break
                elif (ref.start - tar.start) > cluster:
                    out.write(tar.info + "\n")
                    break


def filter_tss_pro(tss_file, pro_file, feature, cluster):
    '''deal with the overlap of TSS and processing site'''
    tsss = read_gff(tss_file)
    pros = read_gff(pro_file)
    out = open("tmp_filter", "w")
    out.write("##gff-version 3\n")
    if feature.lower() == "tss":
        compare_tss_pro(pros, tsss, out, cluster)
        os.remove(pro_file)
        shutil.move("tmp_filter", pro_file)
    elif feature.lower() == "processing":
        compare_tss_pro(tsss, pros, out, cluster)
        os.remove(tss_file)
        shutil.move("tmp_filter", tss_file)
