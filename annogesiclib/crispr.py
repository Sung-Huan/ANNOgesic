import os
import sys
import csv
import shutil
from subprocess import call
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.gff3 import Gff3Parser


class Crispr(object):
    '''Detection of CRISPR'''

    def __init__(self, args_cris):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.gff_parser = Gff3Parser()
        self.gff_path = os.path.join(args_cris.gffs, "tmp")
        self.fasta_path = os.path.join(args_cris.fastas, "tmp")
        self.stat_folder = os.path.join(args_cris.out_folder, "statistics")
        self.gff_out = os.path.join(args_cris.out_folder, "gffs")
        self.all_out = os.path.join(args_cris.out_folder,
                                    "gffs", "all_candidates")
        self.best_out = os.path.join(args_cris.out_folder,
                                     "gffs", "best")
        self.helper.check_make_folder(self.all_out)
        self.helper.check_make_folder(self.best_out)
        self.data_folder = os.path.join(args_cris.out_folder, "CRT_output")
        self.helper.check_make_folder(self.data_folder)
        self.helper.check_make_folder(self.stat_folder)

    def _run_crt(self, args_cris):
        '''Running CRT'''
        for seq in os.listdir(self.fasta_path):
            prefix = ".".join(seq.split(".")[:-1])
            call(["java", "-cp", args_cris.crt_path, "crt", "-minNR",
                  str(args_cris.min_num_r), "-minRL",
                  str(args_cris.min_len_r), "-maxRL",
                  str(args_cris.max_len_r), "-minSL",
                  str(args_cris.min_len_s), "-maxSL",
                  str(args_cris.max_len_s), "-searchWL",
                  str(args_cris.win_size), os.path.join(self.fasta_path, seq),
                  os.path.join(self.data_folder, prefix + ".txt")])

    def _read_gff(self, txt):
        gffs = []
        gh = open(os.path.join(self.gff_path,
                               txt.replace(".txt", ".gff")), "r")
        for entry in Gff3Parser().entries(gh):
            if (entry.feature == "gene") or (
                    entry.feature == "CDS") or (
                    entry.feature == "tRNA") or (
                    entry.feature == "rRNA"):
                gffs.append(entry)
        gh.close()
        return gffs

    def _compare_gff(self, strain, start, end, gffs, bh, indexs, ignore_hypo):
        '''Compare CRISPR and genome annotation to 
        remove the false positives'''
        overlap = False
        id_ = None
        for gff in gffs:
            if (gff.seq_id == strain):
                if ((gff.start <= start) and (gff.end >= end)) or (
                        (gff.start >= start) and (gff.end <= end)) or (
                        (gff.start <= start) and (gff.end > start) and (
                            gff.end <= end)) or (
                        (gff.start >= start) and (gff.start < end) and (
                            gff.end >= end)):
                    if "product" in gff.attributes.keys():
                        if ((not ignore_hypo) and ("hypothetical protein" in
                                gff.attributes["product"])) or (
                                "hypothetical protein" not in 
                                gff.attributes["product"]):
                            overlap = True
        if not overlap:
            id_ = "CRISPR_" + str(indexs["best"])
            attribute = ";".join(["ID=" + id_,
                                  "method=CRT"])
            bh.write("\t".join([strain, "ANNOgesic", "CRISPR", str(start),
                                str(end), ".", ".", ".", attribute]) + "\n")
            indexs["best"] += 1
        return overlap, id_

    def _print_repeat(self, row, strain, file_h, indexs, id_, best):
        '''Print the repeat units'''
        if best:
            num = indexs["re_best"]
        else:
            num = indexs["re_all"]
        if (not row[0].startswith("-")) and (
                not row[0].startswith("Repeats:")) and (
                not row[0].startswith("CRISPR")) and (
                not row[0].startswith("POSITION")):
            start = row[0].strip()
            end = str(int(start) + len(row[2].strip()) - 1)
            attribute = ";".join(["ID=Repeat_" + str(num),
                                  "method=CRT", "Parent=" + id_])
            file_h.write("\t".join([strain, "ANNOgesic", "repeat_unit",
                                    start, end, ".", ".", ".",
                                    attribute]) + "\n")
            num += 1
        if row[0].startswith("Repeats:"):
            indexs["run"] = False
        return num

    def _convert_gff(self, ignore_hypo):
        '''Convert the final CRT output to gff format'''
        for txt in os.listdir(self.data_folder):
            gffs = self._read_gff(txt)
            fh = open(os.path.join(self.data_folder, txt), "r")
            oh = open(os.path.join(
                self.all_out, txt.replace(".txt", "_CRISPR.gff")), "w")
            bh = open(os.path.join(
                self.best_out, txt.replace(".txt", "_CRISPR.gff")), "w")
            indexs = {"all": 0, "re_all": 0, "best": 0,
                      "re_best": 0, "run": False}
            for row in csv.reader(fh, delimiter='\t'):
                if len(row) != 0:
                    if row[0].startswith("ORGANISM:"):
                        strain = row[0].split(" ")[-1]
                    elif row[0].startswith("CRISPR"):
                        end = row[0].split("-")[-1].strip()
                        start = row[0].split("-")[0].split(":")[-1].strip()
                        id_ = "CRISPR_" + str(indexs["all"])
                        attribute = ";".join(["ID=" + id_,
                                              "method=CRT"])
                        oh.write("\t".join([
                            strain, "ANNOgesic", "CRISPR", start,
                            end, ".", ".", ".", attribute]) + "\n")
                        overlap, over_id = self._compare_gff(
                                strain, int(start), int(end), gffs, bh,
                                indexs, ignore_hypo)
                        indexs["all"] += 1
                        indexs["run"] = True
                    if indexs["run"]:
                        indexs["re_all"] = self._print_repeat(
                                row, strain, oh, indexs, id_, False)
                        if not overlap:
                            indexs["re_best"] = self._print_repeat(
                                    row, strain, bh, indexs, over_id, True)
            fh.close()
            oh.close()
            bh.close()

    def _stat_and_correct(self, stats, folder):
        '''do statistics and print the final gff file'''
        for gff in os.listdir(folder):
            prefix = gff.replace("_CRISPR.gff", "")
            stats[prefix] = {"all": {"cri": 0, "re": {}}}
            gh = open(os.path.join(folder, gff), "r")
            oh = open("tmp_cri.gff", "w")
            oh.write("##gff-version 3\n")
            cr_num = 0
            re_num = 0
            first = True
            for entry in Gff3Parser().entries(gh):
                if entry.seq_id not in stats[prefix].keys():
                    stats[prefix][entry.seq_id] = {"cri": 0, "re": {}}
                if entry.feature == "CRISPR":
                    id_ = "CRISPR_" + str(cr_num)
                    attribute = ";".join(["ID=" + id_,
                                          "method=CRT"])
                    cr_num += 1
                    if first:
                        first = False
                    else:
                        if repeat not in stats[prefix][entry.seq_id]["re"].keys():
                            stats[prefix][entry.seq_id]["re"][repeat] = 1
                        else:
                            stats[prefix][entry.seq_id]["re"][repeat] += 1
                        if repeat not in stats[prefix]["all"]["re"].keys():
                            stats[prefix]["all"]["re"][repeat] = 1
                        else:
                            stats[prefix]["all"]["re"][repeat] += 1
                    repeat = 0
                    stats[prefix][entry.seq_id]["cri"] += 1
                    stats[prefix]["all"]["cri"] += 1
                elif entry.feature == "repeat_unit":
                    attribute = ";".join(["ID=Repeat_" + str(re_num),
                                          "method=CRT", "Parent=" + id_])
                    re_num += 1
                    repeat += 1
                oh.write("\t".join([entry.info_without_attributes,
                                    attribute]) + "\n")
            if not first:
                if repeat not in stats[prefix][entry.seq_id]["re"].keys():
                    stats[prefix][entry.seq_id]["re"][repeat] = 1
                else:
                    stats[prefix][entry.seq_id]["re"][repeat] += 1
                if repeat not in stats[prefix]["all"]["re"].keys():
                    stats[prefix]["all"]["re"][repeat] = 1
                else:
                    stats[prefix]["all"]["re"][repeat] += 1
            gh.close()
            oh.close()
            os.remove(os.path.join(folder, gff))
            shutil.move("tmp_cri.gff", os.path.join(folder, gff))

    def _print_file(self, sh, cri_res_all, cri_res_best):
        sh.write("\tthe number of CRISPR - {0}\n".format(
                 cri_res_all["cri"]))
        for index, num in cri_res_all["re"].items():
            sh.write("\t\tCRISPR with {0} repeat units - {1}\n".format(
                         index, num))
        sh.write("\tthe number of CRISPR which not overlap "
                 "with genome annotation - {0}\n".format(
                    cri_res_best["cri"]))
        for index, num in cri_res_best["re"].items():
            sh.write("\t\tCRISPR with {0} repeat units - {1}\n".format(
                        index, num))

    def _print_stat(self, stats):
        '''print the statistics file'''
        for prefix, strains in stats["all"].items():
            sh = open(os.path.join(self.stat_folder, prefix + ".csv"), "w")
            if len(strains) == 1:
                sh.write("No CRISPR can be detected")
            elif len(strains) <= 2:
                for strain, cri_res in strains.items():
                    if strain != "all":
                        sh.write(strain + ":\n")
                        self._print_file(sh, cri_res,
                                         stats["best"][prefix][strain])
            else:
                sh.write("All strains:\n")
                self._print_file(sh, stats["all"][prefix]["all"],
                                 stats["best"][prefix]["all"])
                for strain, cri_res in strains.items():
                    if strain != "all":
                        sh.write(strain + ":\n")
                        if strain not in stats["best"][prefix].keys():
                            stats["best"][prefix][strain] = {"cri": 0,
                                                             "re": {}}
                        self._print_file(sh, cri_res,
                                         stats["best"][prefix][strain])
            sh.close()

    def run_crispr(self, args_cris):
        '''detection of CRISPR'''
        self.multiparser.parser_fasta(args_cris.fastas)
        self.multiparser.parser_gff(args_cris.gffs, None)
        self._run_crt(args_cris)
        self._convert_gff(args_cris.ignore_hypo)
        self.multiparser.combine_gff(args_cris.gffs, self.all_out,
                                     None, "CRISPR")
        self.multiparser.combine_gff(args_cris.gffs, self.best_out,
                                     None, "CRISPR")
        stats = {"all": {}, "best": {}}
        self._stat_and_correct(stats["all"], self.all_out)
        self._stat_and_correct(stats["best"], self.best_out)
        self._print_stat(stats)
        self.helper.remove_tmp(args_cris.gffs)
        self.helper.remove_tmp(args_cris.fastas)
