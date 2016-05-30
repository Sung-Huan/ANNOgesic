import os
import csv
import shutil


def print_rank_all(all_table, best_table):
    out = open("tmp_rank_table", "w")
    fh = open(best_table, "r")
    rank = 0
    bests = []
    for row in csv.reader(fh, delimiter='\t'):
        if row[0] != "rank":
            bests.append(row)
            rank = int(row[0])
        out.write("\t".join(row) + "\n")
    fh.close()
    fh = open(all_table, "r")
    for row in csv.reader(fh, delimiter='\t'):
        detect = False
        if row[0] != "rank":
            for best in bests:
                if (row[1] == best[1]) and (
                        row[3] == best[3]) and (
                        row[4] == best[4]) and (
                        row[5] == best[5]):
                    detect = True
                    break
            if not detect:
                rank += 1
                row[0] = str(rank)
                out.write("\t".join(row) + "\n")
    os.remove(all_table)
    shutil.move("tmp_rank_table", all_table)
