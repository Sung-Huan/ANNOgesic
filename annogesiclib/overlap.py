import os
import shutil
from annogesiclib.gff3 import Gff3Parser

def get_overlap(anno, source, finals, overlaps, detect, out):
    if (anno.source in source) and (
            anno not in overlaps):
        finals.append(anno)
        detect = True
    return detect

def deal_overlap(out_folder, source):
    gffs = {}
    num = 0
    for gff_file in os.listdir(out_folder):
        if gff_file.endswith(".gff"):
            gff_f = open(os.path.join(out_folder, gff_file), "r")
            for entry in Gff3Parser().entries(gff_f):
                if entry.feature not in gffs.keys():
                    gffs[entry.feature] = []
                gffs[entry.feature].append(entry)
            gff_f.close()
            out = open(os.path.join(out_folder, gff_file + "tmp"), "w")
            finals = []
            overlaps = []
            for feature, annos in gffs.items():
                for anno1 in annos:
                    detect = False
                    for anno2 in annos:
                        if (anno1.seq_id == anno2.seq_id) and (
                            anno1.strand == anno2.strand) and (
                            anno1 != anno2) and (
                            anno1.feature == anno2.feature) and (
                            anno1.source != anno2.source):
                            if ((anno1.start <= anno2.start) and (
                                    anno1.end >= anno2.end)) or (
                                    (anno1.start >= anno2.start) and (
                                    anno1.end <= anno2.end)) or (
                                    (anno1.start <= anno2.start) and (
                                    anno1.end <= anno2.end) and (
                                    anno1.end >= anno2.start)) or (
                                    (anno1.start >= anno2.start) and (
                                    anno1.start <= anno2.end) and (
                                    anno1.end >= anno2.end)):
                                detect = get_overlap(anno1, source, finals,
                                                     overlaps, detect, out)
                                detect = get_overlap(anno2, source, finals,
                                                     overlaps, detect, out)
                                if detect:
                                    overlaps.append(anno1)
                                    overlaps.append(anno2)
                    if (not detect) and (anno1 not in overlaps):
                        finals.append(anno1)
            finals = sorted(finals, key=lambda x: (x.seq_id, x.start,
                                                   x.end, x.strand))
            for final in finals:
                if (final.feature == "region") or (
                        final.feature == "source") or (
                        final.feature == "remark"):
                    out.write(final.info + "\n")
                    break
            for final in finals:
                if (final.feature != "region") and (
                        final.feature != "source"):
                    out.write(final.info + "\n")
            out.close()
            shutil.move(os.path.join(out_folder, gff_file + "tmp"),
                        os.path.join(out_folder, gff_file))
