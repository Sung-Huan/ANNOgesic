import os
import csv
from annogesiclib.gff3 import Gff3Parser


def output_coverage(table_file, gff_file, cutoff_cover, stat_file, out_folder):
    out = open(os.path.join(out_folder, "tmp_srna_table"), "w")
    out_g = open(os.path.join(out_folder, "tmp_srna_gff"), "w")
    out.write("\t".join([
        "rank", "strain", "name", "start", "end", "strand",
        "start_with_TSS/Cleavage_site", "end_with_cleavage", "candidates",
        "lib_type", "best_avg_coverage", "best_highest_coverage",
        "best_lower_coverage", "track/coverage",
        "normalized_secondary_energy_change(by_length)",
        "UTR_derived/Intergenic", "confliction_of_sORF", "nr_hit_number",
        "sRNA_hit_number", "nr_hit_top3|ID|e-value", "sRNA_hit|e-value",
        "overlap_CDS", "overlap_percent", "end_with_terminator"]) + "\n")
    out_g.write("##gff-version 3\n")
    stat_out = open(stat_file, "w")
    nums = {5: 0}
    for i in range(10, 100, 10):
        nums[i] = 0
    for i in range(100, 1000, 100):
        nums[i] = 0
    for i in range(1000, 5000, 500):
        nums[i] = 0
    gffs = []
    gh = open(gff_file, "r")
    for entry in Gff3Parser().entries(gh):
        gffs.append(entry)
    fh = open(table_file, "r")
    rank = 1
    new_gffs = []
    for row in csv.reader(fh, delimiter='\t'):
        if row[0] != "rank":
            for cutoff in nums.keys():
                if float(row[10]) >= cutoff:
                    nums[cutoff] += 1
            if float(row[10]) >= cutoff_cover:
                row[0] = str(rank)
                out.write("\t".join(row) + "\n")
                rank += 1
                for gff in gffs:
                    if (row[1] == gff.seq_id) and (
                            row[3] == str(gff.start)) and (
                            row[4] == str(gff.end)) and (
                            row[5] == gff.strand):
                        new_gffs.append(gff)
    sort_gffs = sorted(new_gffs, key=lambda k: (k.seq_id, k.start,
                                                k.end, k.strand))
    for gff in sort_gffs:
        out_g.write(gff.info + "\n")
    coverlist = sorted(nums, key=lambda key: nums[key])
    stat_out.write("coverage\tfrequency\n")
    for cover in coverlist:
        stat_out.write("\t".join([str(cover), str(nums[cover])]) + "\n")
