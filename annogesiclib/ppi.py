#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
import csv
import time
from subprocess import call
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.plot_PPI import plot_ppi


class PPI_network(object):

    def __init__(self, out_folder):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.tmp_id = os.path.join(out_folder, "tmp_id_list")
        self.all_result = os.path.join(out_folder, "all_results")
        self.best_result = os.path.join(out_folder, "best_results")
        self.fig = os.path.join(out_folder, "figures")
        self.with_strain = "with_strain"
        self.without_strain = "without_strain"
        self.tmp_log = "tmp.log"
        self.tmp_action = "tmp_action.log"
        self.tmp_pubmed = "tmp_pubmed.log"
        self.tmp_specific = os.path.join(out_folder, "tmp_specific")
        self.tmp_nospecific = os.path.join(out_folder, "tmp_nospecific")
        self.tmp_wget_action = os.path.join(out_folder, "tmp_action")

    def _make_folder_no_exist(self, path, folder):
        if folder not in os.listdir(path):
            os.mkdir(os.path.join(path, folder))

    def _make_subfolder(self, path, strain, ptt):
        os.mkdir(os.path.join(path, strain))
        os.mkdir(os.path.join(path, strain, ptt))

    def _wget_id(self, genes, strain_id, files):
        detect_id = False
        for gene in genes:
            if gene["strain"] == strain_id["ptt"]:
                print("Retrieving STRING ID for {0} of {1} -- {2}".format(
                       gene["locus_tag"], strain_id["string"], strain_id["file"]))
                id_source = "http://string-db.org/api/tsv/resolve?identifier={0}&species={1}".format(
                             gene["locus_tag"], strain_id["string"])
                call(["wget", id_source,
                      "-O", os.path.join(files["id_list"], gene["locus_tag"])],
                      stderr=files["id_log"])
                time.sleep(3)
                detect_id = True
        return detect_id

    def _retrieve_id(self, strain_id, query_id, genes, files):
        if strain_id["protein"] == "all":
            detect_id = self._wget_id(genes, strain_id, files)
        else:
            detect_id = self._wget_id(genes, strain_id, files)
            if detect_id is False:
                print("Error:there is no {0} in {1}".format(strain["protein"], strain["file"]))

    def _get_prefer_name(self, row_a, strain_id, files, genes):
        prefername = ""
        filename = row_a.split(".")
        if (filename[1] not in os.listdir(files["id_list"])) and \
           (strain_id["protein"] != "all"):
            self._retrieve_id(strain_id, filename[1], genes, files)
        if filename[1] in os.listdir(files["id_list"]):
            id_h = open(os.path.join(files["id_list"], filename[1]), "r")
            for row_i in csv.reader(id_h, delimiter="\t"):
                if row_a == row_i[0]:
                    prefername = row_i[3]
            id_h.close()
        return prefername

    def _print_title(self, out, id_file, id_folder):
        id_h = open(os.path.join(id_folder, id_file), "r")
        prefername = id_file
        for row_i in csv.reader(id_h, delimiter="\t"):
            prefername = row_i[3]
        id_h.close()
        out.write("Interaction of {0} | {1}\n".format(id_file, prefername))
        out.write("strain\titem_id_a\titem_id_b\tmode\taction\ta_is_acting\t"
                   "STRING_action_score\tpubmed_id\tpubmed_score\n")

    def _get_pubmed(self, row, out_folder, strain_id, mode, actor, score, 
                    id_file, first_output, no_specific, ptt, genes, files, paths):
        prefer1 = self._get_prefer_name(row[0], strain_id, files, genes)
        prefer2 = self._get_prefer_name(row[1], strain_id, files, genes)
        if (len(prefer1) > 0) and (len(prefer2) > 0): ## retrieve from PIE
            if no_specific:
                pubmed_source = \
                "http://www.ncbi.nlm.nih.gov/CBBresearch/Wilbur/IRET/PIE/getppi.cgi?term={0}+{1}".format(
                 prefer1, prefer2)
                call(["wget", pubmed_source,
                      "-O", self.tmp_nospecific],
                      stderr=files["pubmed_log"])
                time.sleep(3)
            strain_id["pie"] = "+".join(strain_id["pie"].split(" "))
            pubmed_source = \
            "http://www.ncbi.nlm.nih.gov/CBBresearch/Wilbur/IRET/PIE/getppi.cgi?term={0}+{1}+{2}".format(
             prefer1, prefer2, strain_id["pie"])
            call(["wget", pubmed_source,
                  "-O", self.tmp_specific],
                  stderr=files["pubmed_log"])
            time.sleep(3)
            row[2] = mode ### merge information
            row[4] = actor
            row[0] = prefer1
            row[1] = prefer2
            self._merge_information(first_output, self.tmp_specific, 
                                    files["all_specific"], files["best_specific"], row, score, 
                                    id_file, files["id_list"], "specific",
                                    os.path.join(paths["all"], self.with_strain), 
                                    os.path.join(paths["best"], self.with_strain), ptt)
            if no_specific:
                self._merge_information(first_output, self.tmp_nospecific, 
                                        files["all_nospecific"], files["best_nospecific"], row, 
                                        score, id_file, files["id_list"], "nospecific",
                                        os.path.join(paths["all"], self.without_strain), 
                                        os.path.join(paths["best"], self.without_strain), ptt)

    def _print_single_file(self, out_single, row_a, ptt, row):
        if row == "NA":
            out_single.write("\t".join([ptt, "\t".join(row_a), "NA", "NA"]) + "\n")
        else:
            out_single.write("\t".join([ptt, "\t".join(row_a), "\t".join(row)]) + "\n")

    def _merge_information(self, first_output, filename, out_all, out_best, row_a, score, 
                           id_file, id_folder, file_type, all_folder, best_folder, ptt):
        if os.path.getsize(filename) != 0:
            f_h = open(filename, "r")
            out_all_single = open(os.path.join(all_folder, ptt, 
                                  "_".join([row_a[0], row_a[1] + ".csv"])), "w")
            out_best_single = open(os.path.join(best_folder, ptt, 
                                   "_".join([row_a[0], row_a[1] + ".csv"])), "w")
            self._print_title(out_all_single, id_file, id_folder)
            self._print_title(out_best_single, id_file, id_folder)
            detect = False
            for row in csv.reader(f_h, delimiter="\t"):
                self._print_single_file(out_all_single, row_a, ptt, row)
                if first_output["_".join([file_type, "all"])]:
                    first_output["_".join([file_type, "all"])] = False
                    self._print_title(out_all, id_file, id_folder)
                out_all.write("\t".join([ptt, "\t".join(row_a), "\t".join(row)]) + "\n")
                if (float(row[1]) >= score):
                    detect = True
                    self._print_single_file(out_best_single, row_a, ptt, row)
                    if first_output["_".join([file_type, "best"])]:
                        first_output["_".join([file_type, "best"])] = False
                        self._print_title(out_best, id_file, id_folder)
                    out_best.write("\t".join([ptt, "\t".join(row_a), "\t".join(row)]) + "\n")
            f_h.close()
            if detect is False:
                os.remove(os.path.join(best_folder, ptt, "_".join([row_a[0], row_a[1] + ".csv"])))
            out_all_single.close()
            out_best_single.close()
        else:
            out_all_single = open(os.path.join(all_folder, ptt, 
                                  "_".join([row_a[0], row_a[1] + ".csv"])), "w")
            self._print_title(out_all_single, id_file, id_folder)
            self._print_single_file(out_all_single, row_a, ptt, "NA")
            if first_output["_".join([file_type, "all"])]:
                first_output["_".join([file_type, "all"])] = False
                self._print_title(out_all, id_file, id_folder)
            out_all.write("\t".join([ptt, "\t".join(row_a), "NA", "NA"]) + "\n")
            out_all_single.close()

    def _setup_folder_and_read_file(self, strain_id, pre_file, out_folder, ptts, no_specific,
                                    files, paths):
        if strain_id["file"].endswith(".ptt"):
            if strain_id["file"] != pre_file:
                self.helper.check_make_folder("_".join([self.tmp_id, strain_id["file"]]))
                paths["all"] = os.path.join(self.all_result, strain_id["file"][:-4])
                paths["best"] = os.path.join(self.best_result, strain_id["file"][:-4])
                paths["fig"] = os.path.join(self.fig, strain_id["file"][:-4])
                self.helper.check_make_folder(
                            os.path.join(self.all_result, strain_id["file"][:-4]))
                self.helper.check_make_folder(
                            os.path.join(self.best_result, strain_id["file"][:-4]))
                self.helper.check_make_folder(
                            os.path.join(self.fig, strain_id["file"][:-4]))
                self._make_subfolder(paths["all"], self.with_strain, strain_id["ptt"])
                self._make_subfolder(paths["best"], self.with_strain, strain_id["ptt"])
                self._make_subfolder(paths["fig"], self.with_strain, strain_id["ptt"])
                filename_strain = "_".join([strain_id["file"].replace(".ptt", ""), self.with_strain + ".csv"])
                files["all_specific"] = open(os.path.join(paths["all"], filename_strain), "w")
                files["best_specific"] = open(os.path.join(paths["best"], filename_strain), "w")
                if no_specific:
                    self._make_subfolder(paths["all"], self.without_strain, strain_id["ptt"])
                    self._make_subfolder(paths["best"], self.without_strain, strain_id["ptt"])
                    self._make_subfolder(paths["fig"], self.without_strain, strain_id["ptt"])
                    filename_nostrain = "_".join([strain_id["file"].replace(".ptt", ""), 
                                                  self.without_strain + ".csv"])
                    files["all_nospecific"] = open(os.path.join(paths["all"], filename_nostrain), "w")
                    files["best_nospecific"] = open(os.path.join(paths["best"], filename_nostrain), "w")
                files["id_list"] = "_".join([self.tmp_id, strain_id["file"]])
                files["id_log"] = open(os.path.join(files["id_list"], self.tmp_log), "w")
                files["action_log"] = open(os.path.join(out_folder, self.tmp_action), "w")
                files["pubmed_log"] = open(os.path.join(out_folder, self.tmp_pubmed), "w")
                pre_file = strain_id["file"]
                if strain_id["file"] in os.listdir(ptts):
                    fh = open(os.path.join(ptts, strain_id["file"]), "r")
                    genes = []
                    for row in csv.reader(fh, delimiter="\t"):
                        if (len(row) == 1) and \
                           ("-" in row[0]) and \
                           (".." in row[0]):
                            name = (row[0].split("-"))[0].strip()
                        elif (len(row) > 1) and \
                             (row[0] != "Location"):
                            genes.append({"strain": name, "locus_tag": row[5]})
                    fh.close()
            else:
                self._make_folder_no_exist(os.path.join(paths["all"], self.with_strain), strain_id["ptt"])
                self._make_folder_no_exist(os.path.join(paths["best"], self.with_strain), strain_id["ptt"])
                if no_specific:
                    self._make_folder_no_exist(os.path.join(paths["all"], self.without_strain), strain_id["ptt"])
                    self._make_folder_no_exist(os.path.join(paths["best"], self.without_strain), strain_id["ptt"])
        else:
            print("Error:wrong .ptt file!!")
            sys.exit()
        return genes

    def _wget_actions(self, files, id_file, strain_id, out_folder):
        t_h = open(os.path.join(files["id_list"], id_file), "r")
        print("Retrieving STRING actions for {0} of {1} -- {2}".format(
              id_file, strain_id["string"], strain_id["file"]))
        for row in csv.reader(t_h, delimiter="\t"):
            if row[0].startswith("stringId"): 
                continue
            else: 
                if row[1] == strain_id["string"]:
                    action_source = "http://string-db.org/api/tsv/actions?identifier={0}&species={1}".format(
                                    row[0], row[1])
                    call(["wget", action_source,
                          "-O", self.tmp_wget_action],
                          stderr=files["action_log"])
                    time.sleep(3)
                    break
        t_h.close()

    def _retrieve_actions(self, files, out_folder, strain_id, paths, 
                          score, no_specific, genes):
        for id_file in os.listdir(files["id_list"]):
            if (strain_id["protein"] == id_file) or \
               (strain_id["protein"] == "all"):
                if id_file != self.tmp_log:
                    self._wget_actions(files, id_file, strain_id, out_folder)
                    a_h = open(self.tmp_wget_action, "r") # get prefer name of gene
                    pre_row = []
                    first = True
                    detect = False
                    first_output = {"specific_all": True, "specific_best": True,
                                    "nospecific_all": True, "nospecific_best": True}
                    print("Retrieving Pubmed for {0} of {1} -- {2}".format(
                           id_file, strain_id["string"], strain_id["file"]))
                    for row_a in csv.reader(a_h, delimiter="\t"):
                        if row_a == []:
                            print("No interaction can be detected...")
                            break
                        if row_a[0].startswith("item_id_a"):
                            continue
                        else:
                            detect = True
                            if first:
                                first = False
                                mode = row_a[2]
                                actor = row_a[4]
                            else:
                                if (row_a[0] != pre_row[0]) or \
                                   (row_a[1] != pre_row[1]):
                                    self._get_pubmed(pre_row, out_folder, strain_id, 
                                                     mode, actor, score, id_file, 
                                                     first_output, no_specific, 
                                                     strain_id["ptt"], genes, files, paths)
                                    mode = row_a[2]
                                    actor = row_a[4]
                                else:
                                    mode = mode + ";" + row_a[2]
                                    actor = actor + ";" + row_a[4]
                            pre_row = row_a
                    if detect:
                        detect = False
                        self._get_pubmed(row_a, out_folder, strain_id, mode, 
                                         actor, score, id_file, first_output, 
                                         no_specific, strain_id["ptt"], genes,
                                         files, paths)

    def _plot(self, out_folder, score, size, no_specific):
        for folder in os.listdir(self.all_result): ## plot figure
            if folder in os.listdir(self.fig):
                print("plotting {0}".format(folder))
                plot_ppi(os.path.join(self.all_result, folder,
                                      "_".join([folder, self.with_strain + ".csv"])),
                         score, os.path.join(self.fig, folder, self.with_strain),
                         size)
                if no_specific:
                    plot_ppi(os.path.join(self.all_result, folder,
                                          "_".join([folder, self.without_strain + ".csv"])),
                         score, os.path.join(self.fig, folder, self.without_strain),
                         size)

    def retrieve_ppi_network(self, bin_path, ptts, strains, no_specific,
                             species, score, out_folder, size):
        strain_ids = []
        genes = []
        paths = {}
        files = {}
        for strain in strains: # import strain information
            datas = strain.split(":")
            strain_ids.append({"protein": datas[0], "file": datas[1], "ptt":datas[2],
                               "string": datas[3], "pie": datas[4]})
        strain_ids.sort(key=lambda x: x["file"])
        pre_file = ""
        file_holders = {"all_nospecific": None, "best_nospecific": None,
                        "all_specific": None, "best_specific": None,
                        "action_log": None, "pubmed_log": None, "id_log": None,
                        "all": None, "best": None, "fig": None}
        for strain_id in strain_ids:
            genes = self._setup_folder_and_read_file(strain_id, pre_file, out_folder, 
                                                     ptts, no_specific, files, paths)
            s_h = open(species, "r") ### get STRING id
            for row in csv.reader(s_h, delimiter="\t"):
                if row[0] != "##":
                    if row[0] == strain_id["string"]:
                        break
                    elif row[2] == strain_id["string"]:
                        strain_id["string"] = row[0]
                        break
                    elif row[3] == strain_id["string"]:
                        strain_id["string"] = row[0]
                        break
            self._retrieve_id(strain_id, strain_id["protein"], genes, files)
            self._retrieve_actions(files, out_folder, strain_id, paths, 
                                   score, no_specific, genes)
        self._plot(out_folder, score, size, no_specific)
        self.helper.remove_all_content(os.path.join(out_folder), "tmp", "file")
        self.helper.remove_all_content(os.path.join(out_folder), "tmp", "dir")
