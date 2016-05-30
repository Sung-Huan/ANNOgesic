import os
import shutil
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.lib_reader import read_wig, read_libs


def detect_coverage(wigs, tran, infos):
    for strain, conds in wigs.items():
        if strain == tran.seq_id:
            for cond, tracks in conds.items():
                for track, covers in tracks.items():
                    infos[track] = {"avg": -1, "high": -1, "low": -1}
                    total = 0
                    for cover in covers[tran.start-2:tran.end+1]:
                        if (cover["pos"] >= tran.start) and (
                                cover["pos"] <= tran.end):
                            total = cover["coverage"] + total
                            if cover["coverage"] > infos[track]["high"]:
                                infos[track]["high"] = cover["coverage"]
                            if (cover["coverage"] < infos[track]["low"]) or (
                                    infos[track]["low"] == -1):
                                infos[track]["low"] = cover["coverage"]
                    infos[track]["avg"] = (float(total) /
                                           float(tran.end - tran.start + 1))

def print_associate(associate, tran, out):
    if associate in tran.attributes.keys():
        out.write("\t" + tran.attributes[associate])
    else:
        out.write("\tNA")

def print_coverage(trans, out, out_gff, wigs_f, wigs_r, table_best):
    for tran in trans:
        infos = {}
        tran.attributes["detect_lib"] = tran.attributes["detect_lib"].replace(
                                                        "tex_notex", "TEX+/-")
        out.write("\t".join([tran.seq_id, tran.attributes["Name"],
                             str(tran.start), str(tran.end), tran.strand,
                             tran.attributes["detect_lib"]]))
        print_associate("associated_gene", tran, out)
        print_associate("associated_tss", tran, out)
        print_associate("associated_term", tran, out)
        if tran.strand == "+":
            detect_coverage(wigs_f, tran, infos)
        else:
            detect_coverage(wigs_r, tran, infos)
        out.write("\t")
        best = -1
        best_track = ""
        best_cover = {}
        for track, cover in infos.items():
            if not table_best:
                if best != -1:
                    out.write(";")
                out.write("{0}(avg={1})".format(
                          track, str(cover["avg"])))
            if cover["avg"] > best:
                best = cover["avg"]
                best_track = track
                best_cover = cover
        if table_best:
            out.write("{0}(avg={1})".format(
                      best_track, str(best_cover["avg"])))
        out.write("\n")
        new_attrs = {}
        for key, value in tran.attributes.items():
            if ("high_coverage" not in key) and (
                    "low_coverage" not in key):
                new_attrs[key] = value
        new_attrs["best_avg_coverage"] = str(best_cover["avg"])
        attribute_string = ";".join(
            ["=".join(items) for items in new_attrs.items()])
        out_gff.write("\t".join([tran.info_without_attributes,
                                 attribute_string]) + "\n")

def gen_table_transcript(gff_folder, args_tran):
    libs, texs = read_libs(args_tran.libs, args_tran.merge_wigs)
    for gff in os.listdir(gff_folder):
        if os.path.isfile(os.path.join(gff_folder, gff)):
            wigs_f = read_wig(os.path.join(args_tran.wig_path, "_".join([
                              gff.replace("_transcript.gff", ""),
                              "forward.wig"])), "+", libs)
            wigs_r = read_wig(os.path.join(args_tran.wig_path, "_".join([
                              gff.replace("_transcript.gff", ""),
                              "reverse.wig"])), "-", libs)
            th = open(os.path.join(gff_folder, gff), "r")
            trans = []
            out = open(os.path.join(args_tran.out_folder, "tables",
                       gff.replace(".gff", ".csv")), "w")
            out_gff = open(os.path.join(args_tran.out_folder, "tmp_gff"), "w")
            out_gff.write("##gff-version 3\n")
            out.write("\t".join(["strain", "Name", "start", "end", "strand",
                                 "detect_lib_type", "associated_gene",
                                 "associated_tss", "associated_term",
                                 "coverage_details"]) + "\n")
            gff_parser = Gff3Parser()
            for entry in gff_parser.entries(th):
                trans.append(entry)
            print_coverage(trans, out, out_gff, wigs_f, wigs_r,
                           args_tran.table_best)
            out.close()
            out_gff.close()
            shutil.move(os.path.join(args_tran.out_folder, "tmp_gff"),
                        os.path.join(gff_folder, gff))
    if os.path.exists(os.path.join(args_tran.out_folder, "merge_wigs")):
        shutil.rmtree(os.path.join(args_tran.out_folder, "merge_wigs"))
