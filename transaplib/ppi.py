#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
import csv
import time
from subprocess import call
import multiparser
import check_gff_attributes
class PPI_network(object):

    def _check_make_folder(self, path, folder):
        if folder in os.listdir(path):
            call(["rm", "-rf", path + folder])
        call(["mkdir", path + folder])

    def _make_folder_no_exist(self, path, folder):
        if folder not in os.listdir(path):
            call(["mkdir", path + folder])

    def _make_subfolder(self, path, strain, ptt):
        call(["mkdir", path + strain])
        call(["mkdir", path + strain + "/" + ptt])

    def _retrieve_id(self, strain_id, query_id, genes, tmp_id_list, out_id_log):
        if strain_id["protein"] == "all":
            for gene in genes:
                if gene["strain"] == strain_id["ptt"]:
                    print("Retrieving STRING ID for " + gene["locus_tag"] + " of " +
                            strain_id["string"] + "--" + strain_id["file"])
                    id_source = "http://string-db.org/api/tsv/resolve?identifier=%s&species=%s" % \
                             (gene["locus_tag"], strain_id["string"])
                    call(["wget", id_source,
                          "-O", tmp_id_list + gene["locus_tag"]],
                          stderr=out_id_log)
                    time.sleep(3)
        else:
            detect_id = False
            for gene in genes:
                if (gene["strain"] == strain_id["ptt"]) and \
                   (gene["locus_tag"] == query_id):
                    print("Retrieving STRING ID for " + gene["locus_tag"] + " of " +
                            strain_id["string"] + "--" + strain_id["file"])
                    id_source = "http://string-db.org/api/tsv/resolve?identifier=%s&species=%s" % \
                             (gene["locus_tag"], strain_id["string"])
                    call(["wget", id_source,
                          "-O", tmp_id_list + gene["locus_tag"]],
                          stderr=out_id_log)
                    time.sleep(3)
                    detect_id = True
            if detect_id is False:
                print("Error:there is no " + strain["protein"] + " in " + strain["file"])

    def _get_prefer_name(self, row_a, strain_id, id_folder, genes, out_id_log):
        prefername = ""
        filename = row_a.split(".")
        if (filename[1] not in os.listdir(id_folder)) and \
           (strain_id["protein"] != "all"):
            self._retrieve_id(strain_id, filename[1], genes, id_folder, out_id_log)
        if filename[1] in os.listdir(id_folder):
            id_h = open(id_folder + filename[1], "r")
            for row_i in csv.reader(id_h, delimiter="\t"):
                if row_a == row_i[0]:
                    prefername = row_i[3]
            id_h.close()
        return prefername

    def _print_title(self, out, id_file, id_folder):
        id_h = open(id_folder + id_file, "r")
        prefername = id_file
        for row_i in csv.reader(id_h, delimiter="\t"):
            prefername = row_i[3]
        id_h.close()
        out.write("Interaction of " + id_file + "|" + prefername + "\n")
        out.write("strain\titem_id_a\titem_id_b\tmode\taction\ta_is_acting\t"
                   "STRING_action_score\tpubmed_id\tpubmed_score\n")

    def _get_pubmed(self, row, id_folder, out_folder, strain_id,
                    mode, actor, score, id_file,
                    out_all_specific, out_best_specific,
                    out_all_nospecific, out_best_nospecific,
                    first_output, out_pubmed_log, all_path, 
                    best_path, no_specific, ptt, genes, out_id_log):
        prefer1 = self._get_prefer_name(row[0], strain_id, id_folder, genes, out_id_log)
        prefer2 = self._get_prefer_name(row[1], strain_id, id_folder, genes, out_id_log)
        ####################################
        # retrieve Pubmed from Pie         #
        ####################################
        if (len(prefer1) > 0) and (len(prefer2) > 0):
            if no_specific:
                pubmed_source = "http://www.ncbi.nlm.nih.gov/CBBresearch/Wilbur/IRET/PIE/getppi.cgi?term=%s+%s" % \
                                (prefer1, prefer2)
                call(["wget", pubmed_source,
                      "-O", out_folder + "/tmp_nospecific"],
                      stderr=out_pubmed_log)
                time.sleep(3)
            strain_id["pie"] = "+".join(strain_id["pie"].split(" "))
            pubmed_source = "http://www.ncbi.nlm.nih.gov/CBBresearch/Wilbur/IRET/PIE/getppi.cgi?term=%s+%s+%s" % \
                            (prefer1, prefer2, strain_id["pie"])
            call(["wget", pubmed_source,
                  "-O", out_folder + "/tmp_specific"],
                  stderr=out_pubmed_log)
            time.sleep(3)
            #############################################
            # merge all information to one file         #
            #############################################
            row[2] = mode
            row[4] = actor
            row[0] = prefer1
            row[1] = prefer2
            self._merge_information(first_output, out_folder + "/tmp_specific", out_all_specific,
                                    out_best_specific, row, score, id_file, id_folder, "specific",
                                    all_path + "with_strain/", best_path + "with_strain/", ptt)
            if no_specific:
                self._merge_information(first_output, out_folder + "/tmp_nospecific", out_all_nospecific,
                                        out_best_nospecific, row, score, id_file, id_folder, "nospecific",
                                        all_path + "without_strain/", best_path + "without_strain/", ptt)

    def _print_single_file(self, out_single, row_a, ptt, row):
        if row == "NA":
            out_single.write("\t".join([ptt, "\t".join(row_a), "NA", "NA"]) + "\n")
        else:
            out_single.write("\t".join([ptt, "\t".join(row_a), "\t".join(row)]) + "\n")

    def _merge_information(self, first_output, filename, out_all, out_best, 
                           row_a, score, id_file, id_folder, file_type, all_folder, best_folder, ptt):
        if os.path.getsize(filename) != 0:
            f_h = open(filename, "r")
            out_all_single = open(all_folder + ptt + "/" + row_a[0] + "_" + row_a[1] + ".csv", "w")
            out_best_single = open(best_folder + ptt + "/" + row_a[0] + "_" + row_a[1] + ".csv", "w")
            self._print_title(out_all_single, id_file, id_folder)
            self._print_title(out_best_single, id_file, id_folder)
            detect = False
            for row in csv.reader(f_h, delimiter="\t"):
                self._print_single_file(out_all_single, row_a, ptt, row)
                if first_output[file_type + "_all"]:
                    first_output[file_type + "_all"] = False
                    self._print_title(out_all, id_file, id_folder)
                out_all.write("\t".join([ptt, "\t".join(row_a), "\t".join(row)]) + "\n")
                if (float(row[1]) >= score):
                    detect = True
                    self._print_single_file(out_best_single, row_a, ptt, row)
                    if first_output[file_type + "_best"]:
                        first_output[file_type + "_best"] = False
                        self._print_title(out_best, id_file, id_folder)
                    out_best.write("\t".join([ptt, "\t".join(row_a), "\t".join(row)]) + "\n")
            f_h.close()
            if detect is False:
                call(["rm", best_folder + ptt + "/" + row_a[0] + "_" + row_a[1] + ".csv"])
            out_all_single.close()
            out_best_single.close()
        else:
            out_all_single = open(all_folder + ptt + "/" + row_a[0] + "_" + row_a[1] + ".csv", "w")
            self._print_title(out_all_single, id_file, id_folder)
            self._print_single_file(out_all_single, row_a, ptt, "NA")
            if first_output[file_type + "_all"]:
                first_output[file_type + "_all"] = False
                self._print_title(out_all, id_file, id_folder)
            out_all.write("\t".join([ptt, "\t".join(row_a), "NA", "NA"]) + "\n")
            out_all_single.close()

    def retrieve_PPI_network(self, bin_path, ptts, strains, no_specific,
                             species, score, out_folder, size):
        strain_ids = []
        ###############################
        # import strain information   #
        ###############################
        for strain in strains:
            datas = strain.split(":")
            strain_ids.append({"protein": datas[0], "file": datas[1], "ptt":datas[2],
                               "string": datas[3], "pie": datas[4]})
        strain_ids.sort(key=lambda x: x["file"])
        pre_file = ""
        out_all_nospecific = ""
        out_best_nospecific = ""
        for strain_id in strain_ids:
            if strain_id["file"].endswith(".ptt"):
                if strain_id["file"] != pre_file:
                    self._check_make_folder(out_folder + "/", "tmp_id_list_" + strain_id["file"])
                    all_path = out_folder + "/all_results/" + strain_id["file"][:-4] + "/"
                    best_path = out_folder + "/best_results/" + strain_id["file"][:-4] + "/"
                    fig_path = out_folder + "/figures/" + strain_id["file"][:-4] + "/"
                    self._check_make_folder(out_folder + "/all_results/", strain_id["file"][:-4])
                    self._check_make_folder(out_folder + "/best_results/", strain_id["file"][:-4])
                    self._check_make_folder(out_folder + "/figures/", strain_id["file"][:-4])
                    self._make_subfolder(all_path, "with_strain", strain_id["ptt"])
                    self._make_subfolder(best_path, "with_strain", strain_id["ptt"])
                    self._make_subfolder(fig_path, "with_strain", strain_id["ptt"])
                    filename_strain = strain_id["file"].replace(".ptt", "") + "_with_strain.csv"
                    out_all_specific = open(all_path + filename_strain, "w")
                    out_best_specific = open(best_path + filename_strain, "w")
                    if no_specific:
                        self._make_subfolder(all_path, "without_strain", strain_id["ptt"])
                        self._make_subfolder(best_path, "without_strain", strain_id["ptt"])
                        self._make_subfolder(fig_path, "without_strain", strain_id["ptt"])
                        filename_nostrain = strain_id["file"].replace(".ptt", "") + "_without_strain.csv"
                        out_all_nospecific = open(all_path + filename_nostrain, "w")
                        out_best_nospecific = open(best_path + filename_nostrain, "w")
                    tmp_id_list = out_folder + "/tmp_id_list_" + strain_id["file"] + "/"
                    out_id_log = open(tmp_id_list + "tmp.log", "w")
                    out_action_log = open(out_folder + "/tmp_action.log", "w")
                    out_pubmed_log = open(out_folder + "/tmp_pubmed.log", "w")
                    pre_file = strain_id["file"]
                    if strain_id["file"] in os.listdir(ptts):
                        fh = open(ptts + "/" + strain_id["file"], "r")
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
                    self._make_folder_no_exist(all_path + "with_strain/", strain_id["ptt"])
                    self._make_folder_no_exist(best_path + "with_strain/", strain_id["ptt"])
                    if no_specific:
                        self._make_folder_no_exist(all_path + "without_strain/", strain_id["ptt"])
                        self._make_folder_no_exist(best_path + "without_strain/", strain_id["ptt"])
            else:
                print("Error:wrong .ptt file!!")
                sys.exit()
            ####################################
            # get the proper STRING id number  #
            ####################################
            s_h = open(species, "r")
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
            ##################################
            # retrieve STRING ID from STRING #
            ##################################
            self._retrieve_id(strain_id, strain_id["protein"], genes, tmp_id_list, out_id_log)
            ##################################
            # retrieve actions from STRING   #
            ##################################
            for id_file in os.listdir(tmp_id_list):
                if (strain_id["protein"] == id_file) or \
                   (strain_id["protein"] == "all"):
                    if id_file != "tmp.log":
                        t_h = open(tmp_id_list + id_file, "r")
                        print("Retrieving STRING actions for " + id_file + " of " + 
                              strain_id["string"] + "--" + strain_id["file"])
                        for row in csv.reader(t_h, delimiter="\t"):
                            if row[0].startswith("stringId"):
                                continue
                            else:
                                if row[1] == strain_id["string"]:
                                    gene_ids = {"locus_tag": row[0], "prefer": row[3]}
                                    action_source = "http://string-db.org/api/tsv/actions?identifier=%s&species=%s" % \
                                                    (row[0], row[1])
                                    call(["wget", action_source,
                                          "-O", out_folder + "/tmp_action"],
                                          stderr=out_action_log)
                                    time.sleep(3)
                                    break
                        t_h.close()
                        ###########################################################
                        # get prefer name of gene in order to retrieve Pubmed     #
                        ###########################################################
                        a_h = open(out_folder + "/tmp_action", "r")
                        pre_row = []
                        first = True
                        detect = False
                        first_output = {"specific_all": True, "specific_best": True,
                                        "nospecific_all": True, "nospecific_best": True}
                        print("Retrieving Pubmed for " + id_file + " of " + 
                              strain_id["string"] + "--" + strain_id["file"])
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
                                        self._get_pubmed(pre_row, tmp_id_list, out_folder, strain_id,
                                                         mode, actor, score, id_file,
                                                         out_all_specific, out_best_specific,
                                                         out_all_nospecific, out_best_nospecific,
                                                         first_output, out_pubmed_log, 
                                                         all_path, best_path, no_specific, strain_id["ptt"],
                                                         genes, out_id_log)
                                        mode = row_a[2]
                                        actor = row_a[4]
                                    else:
                                        mode = mode + ";" + row_a[2]
                                        actor = actor + ";" + row_a[4]
                                pre_row = row_a
                        if detect:
                            detect = False
                            self._get_pubmed(row_a, tmp_id_list, out_folder, strain_id,
                                             mode, actor, score, id_file,
                                             out_all_specific, out_best_specific,
                                             out_all_nospecific, out_best_nospecific,
                                             first_output, out_pubmed_log, 
                                             all_path, best_path, no_specific, strain_id["ptt"], 
                                             genes, out_id_log)
        ########################################
        # generate figures of PPI network      #
        ########################################
        for folder in os.listdir(out_folder + "/all_results"):
            if folder in os.listdir(out_folder + "/figures"):
                print("plotting " + folder)
                call(["python", bin_path + "/plot_PPI.py",
                      "-i", out_folder + "/all_results/" + folder + "/" + folder + "_with_strain.csv",
                      "-s", str(score), "-f", out_folder + "/figures/" + folder + "/with_strain",
                      "-ns", str(size)])
                if no_specific:
                    call(["python", bin_path + "/plot_PPI.py",
                          "-i", out_folder + "/all_results/" + folder + "/" + folder + "_without_strain.csv",
                          "-s", str(score), "-f", out_folder + "/figures/" + folder + "/without_strain",
                          "-ns", str(size)])
        os.system("rm -rf " + out_folder + "/tmp*")
