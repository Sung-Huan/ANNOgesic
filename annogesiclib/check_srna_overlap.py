import csv
import shutil
from annogesiclib.helper import Helper
from annogesiclib.gff3 import Gff3Parser


def import_cds(gff):
    if "Name" in gff.attributes.keys():
        return gff.attributes["Name"]
    elif "ID" in gff.attributes.keys():
        return gff.attributes["ID"]
    else:
        name = "".join([gff.feature, ":", str(gff.start), "-", str(gff.end),
                        "_", gff.strand])
        return name


def check_overlap(table_file, gff_file):
    out = open(table_file + "tmp", "w")
    gffs = []
    gff_f = open(gff_file, "r")
    for entry in Gff3Parser().entries(gff_f):
        if Helper().feature_without_notgene(entry):
            gffs.append(entry)
    fh = open(table_file, "r")
    out.write("\t".join([
        "Rank", "Genome", "Name", "Start", "End", "Strand",
        "Start_with_TSS/Cleavage_site", "End_with_cleavage", "Candidates",
        "Lib_type", "Best_avg_coverage", "Track/Coverage",
        "Normalized_secondary_energy_change(by_length)", "sRNA_types",
        "Conflict_sORF", "nr_hit_number", "sRNA_hit_number",
        "nr_hit_top3|ID|e-value|score", "sRNA_hit|e-value|score", "Overlap_CDS_forward",
        "Overlap_nts_forward", "Overlap_CDS_reverse",
        "Overlap_nts_reverse","End_with_terminator",
        "Associated_promoter", "sRNA_length"]) + "\n")
    for row in csv.reader(fh, delimiter='\t'):
        if row[3] != "Start":
            overlaps = {"forward": [], "reverse": [],
                        "CDS_f": [], "CDS_r": []}
            start = int(row[3])
            end = int(row[4])
            for gff in gffs:
                if ((gff.end < end) and (
                         gff.end > start) and (
                         gff.start <= start)) or (
                        (gff.start > start) and (
                         gff.start < end) and (
                         gff.end >= end)) or (
                        (gff.end >= end) and (
                         gff.start <= start)) or (
                        (gff.end <= end) and (
                         gff.start >= start)):
                    overlap = min(gff.end, end) - max(gff.start, start) + 1
                    percent = "{0:.0f}%".format((float(overlap) / float(end - start + 1)) * 100)
                    if gff.strand == "+":
                        overlaps["forward"].append(str(overlap) + "(" + str(percent) + ")")
                        overlaps["CDS_f"].append(import_cds(gff))
                    else:
                        overlaps["reverse"].append(str(overlap) + "(" + str(percent) + ")")
                        overlaps["CDS_r"].append(import_cds(gff))
            if len(overlaps["forward"]) == 0:
                overlaps["forward"] = ["NA"]
                overlaps["CDS_f"] = ["NA"]
            if len(overlaps["reverse"]) == 0:
                overlaps["reverse"] = ["NA"]
                overlaps["CDS_r"] = ["NA"]
            out.write("\t".join(row[0:19] + [";".join(overlaps["CDS_f"]), ";".join(overlaps["forward"]),
                                             ";".join(overlaps["CDS_r"]), ";".join(overlaps["reverse"])] +
                                row[21:]) + "\n")
    shutil.move(table_file + "tmp", table_file)
