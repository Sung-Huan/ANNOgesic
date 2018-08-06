import os
import sys
import csv
import time
from subprocess import call
from annogesiclib.helper import Helper
from annogesiclib.plot_PPI import plot_ppi
from annogesiclib.converter import Converter
from annogesiclib.gff3 import Gff3Parser


class PPINetwork(object):
    '''detection of PPI'''

    def __init__(self, out_folder):
        self.helper = Helper()
        self.converter = Converter()
        self.gffparser = Gff3Parser()
        self.tmp_id = os.path.join(out_folder, "tmp_id_list")
        self.all_result = os.path.join(out_folder, "all_results")
        self.best_result = os.path.join(out_folder, "best_results")
        self.fig = os.path.join(out_folder, "figures")
        self.ref_tags = {}
        self.with_strain = "with_strain"
        self.without_strain = "without_strain"
        self.tmp_files = {"log": "tmp_log", "action": "tmp_action.log",
                          "pubmed": "tmp_pubmed.log",
                          "specific": os.path.join(
                                      out_folder, "tmp_specific"),
                          "nospecific": os.path.join(
                                        out_folder, "tmp_nospecific"),
                          "wget_action": os.path.join(
                                         out_folder, "tmp_action")}

    def _make_folder_no_exist(self, path, folder):
        if folder not in os.listdir(path):
            os.mkdir(os.path.join(path, folder))

    def _make_subfolder(self, path, strain, ptt):
        os.mkdir(os.path.join(path, strain))
        os.mkdir(os.path.join(path, strain, ptt))

    def _run_wget(self, source, folder, err, log):
        log.write(" ".join(["wget", source, "-O", folder]) + "\n")
        call(["wget", source, "-O", folder], stderr=err)
        time.sleep(2)

    def _wget_id(self, strain, locus, strain_id, files, log):
        detect_id = False
        if strain == strain_id["ptt"]:
            print("Retrieving STRING ID for {0} of {1} -- {2}".format(
                   locus, strain_id["string"], strain_id["file"]))
            id_source = ("http://string-db.org/api/tsv/resolve?"
                         "identifier={0}&species={1}").format(
                         locus, strain_id["string"])
            self._run_wget(id_source, os.path.join(files["id_list"], locus),
                           files["id_log"], log)
            detect_id = True
        return detect_id

    def _retrieve_id(self, strain_id, genes, files, log):
        log.write("Retrieving STRING ID for {0}.\n".format(strain_id["ptt"]))
        for gene in genes:
            if gene["gene"] != "-":
                detect_id = self._wget_id(gene["strain"], gene["gene"],
                                          strain_id, files, log)
                self.ref_tags[gene["gene"]] = gene["locus_tag"]
            else:
                detect_id = self._wget_id(gene["strain"], gene["locus_tag"],
                                          strain_id, files, log)
                self.ref_tags[gene["locus_tag"]] = gene["locus_tag"]
            if not detect_id:
                log.write("{0} is not found in {1}.\n".format(
                    gene, strain_id["file"]))
                print("Error: There is no {0} in {1}".format(
                       gene, strain_id["file"]))
        log.write("The temporary files are generated and stored in the "
                  "following folders:\n")
        log.write("\t" + os.path.join(
            files["id_list"], gene["locus_tag"]) + "\n")

    def _get_prefer_name(self, row_a, strain_id, files, querys, log):
        prefername = ""
        filename = row_a.split(".")
        if ((filename[1] not in os.listdir(files["id_list"])) and (
                "all" not in querys)) or ("all" in querys):
            self._wget_id(strain_id["ptt"], filename[1], strain_id, files, log)
        if (filename[1] in os.listdir(files["id_list"])) or (
                "all" in querys):
            if (filename[1] in os.listdir(files["id_list"])):
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
        if prefername not in self.ref_tags.keys():
            locus = id_file
        else:
            locus = self.ref_tags[prefername]
        out.write("Interaction of {0} | {1}\n".format(
            locus, prefername))
        out.write("Genome\tItem_id_a\tItem_id_b\tMode\tAction\ta_is_acting\t"
                  "STRING_action_score\tPubmed_id\tPubmed_score\n")

    def _get_pubmed(self, row, strain_id, mode, actor, id_file, first_output,
                    ptt, files, paths, args_ppi, log):
        prefer1 = self._get_prefer_name(row[0], strain_id,
                                        files, args_ppi.querys, log)
        prefer2 = self._get_prefer_name(row[1], strain_id,
                                        files, args_ppi.querys, log)
        if (len(prefer1) > 0) and (len(prefer2) > 0):
            if args_ppi.no_specific:
                pubmed_source = (
                    "http://www.ncbi.nlm.nih.gov/CBBresearch/"
                    "Wilbur/IRET/PIE/getppi.cgi?term={0}+{1}").format(
                        prefer1, prefer2)
                self._run_wget(pubmed_source, self.tmp_files["nospecific"],
                               files["pubmed_log"], log)
            strain_id["pie"] = "+".join(strain_id["pie"].split(" "))
            pubmed_source = (
                "http://www.ncbi.nlm.nih.gov/CBBresearch/Wilbur"
                "/IRET/PIE/getppi.cgi?term={0}+{1}+{2}").format(
                    prefer1, prefer2, strain_id["pie"])
            self._run_wget(pubmed_source, self.tmp_files["specific"],
                           files["pubmed_log"], log)
            row[2] = mode
            row[4] = actor
            row[0] = prefer1
            row[1] = prefer2
            self._merge_information(
                first_output, self.tmp_files["specific"],
                files["all_specific"], files["best_specific"], row,
                args_ppi.score, id_file, files["id_list"], "specific",
                os.path.join(paths["all"], self.with_strain),
                os.path.join(paths["best"], self.with_strain), ptt)
            if args_ppi.no_specific:
                self._merge_information(
                     first_output, self.tmp_files["nospecific"],
                     files["all_nospecific"], files["best_nospecific"], row,
                     args_ppi.score, id_file, files["id_list"], "nospecific",
                     os.path.join(paths["all"], self.without_strain),
                     os.path.join(paths["best"], self.without_strain), ptt)

    def _print_single_file(self, out_single, row_a, ptt, row):
        if row == "NA":
            out_single.write("\t".join(
                             [ptt, "\t".join(row_a), "NA", "NA"]) + "\n")
        else:
            out_single.write("\t".join(
                             [ptt, "\t".join(row_a), "\t".join(row)]) + "\n")

    def _merge_information(self, first_output, filename, out_all, out_best,
                           row_a, score, id_file, id_folder, file_type,
                           all_folder, best_folder, ptt):
        if os.path.getsize(filename) != 0:
            f_h = open(filename, "r")
            out_all_single = open(os.path.join(
                all_folder, ptt, "_".join([row_a[0], row_a[1] + ".csv"])), "w")
            out_best_single = open(os.path.join(
                best_folder, ptt,
                "_".join([row_a[0], row_a[1] + ".csv"])), "w")
            self._print_title(out_all_single, id_file, id_folder)
            self._print_title(out_best_single, id_file, id_folder)
            detect = False
            for row in csv.reader(f_h, delimiter="\t"):
                self._print_single_file(out_all_single, row_a, ptt, row)
                if first_output["_".join([file_type, "all"])]:
                    first_output["_".join([file_type, "all"])] = False
                    self._print_title(out_all, id_file, id_folder)
                out_all.write("\t".join([ptt, "\t".join(row_a),
                                         "\t".join(row)]) + "\n")
                if (float(row[1]) >= score):
                    detect = True
                    self._print_single_file(out_best_single, row_a, ptt, row)
                    if first_output["_".join([file_type, "best"])]:
                        first_output["_".join([file_type, "best"])] = False
                        self._print_title(out_best, id_file, id_folder)
                    out_best.write("\t".join([ptt, "\t".join(row_a),
                                              "\t".join(row)]) + "\n")
            f_h.close()
            if not detect:
                os.remove(os.path.join(best_folder, ptt,
                          "_".join([row_a[0], row_a[1] + ".csv"])))
            out_all_single.close()
            out_best_single.close()
        else:
            out_all_single = open(os.path.join(
                all_folder, ptt, "_".join([row_a[0], row_a[1] + ".csv"])), "w")
            self._print_title(out_all_single, id_file, id_folder)
            self._print_single_file(out_all_single, row_a, ptt, "NA")
            if first_output["_".join([file_type, "all"])]:
                first_output["_".join([file_type, "all"])] = False
                self._print_title(out_all, id_file, id_folder)
            out_all.write("\t".join([ptt, "\t".join(row_a),
                                     "NA", "NA"]) + "\n")
            out_all_single.close()

    def _detect_protein(self, strain_id, args_ppi):
        fh = open(os.path.join(args_ppi.ptts, strain_id["file"]), "r")
        genes = []
        for row in csv.reader(fh, delimiter="\t"):
            if (len(row) == 1) and ("-" in row[0]) and (".." in row[0]):
                name = (row[0].split("-"))[0].strip().split(",")[0].strip()
            if ("all" in args_ppi.querys):
                if (len(row) > 1) and (row[0] != "Location"):
                    genes.append({"strain": name, "locus_tag": row[4],
                                  "gene": row[5]})
            else:
                for query in args_ppi.querys:
                    datas = query.split(":")
                    strain = datas[0]
                    start = datas[1]
                    end = datas[2]
                    strand = datas[3]
                    if (len(row) > 1) and (row[0] != "Location") and (
                            name == strain) and (
                            start == row[0].split("..")[0]) and (
                            end == row[0].split("..")[1]) and (
                            strand == row[1]):
                        genes.append({"strain": name, "locus_tag": row[4],
                                      "gene": row[5]})
        fh.close()
        return genes

    def _setup_nospecific(self, paths, strain_id, files):
        self._make_subfolder(
            paths["all"], self.without_strain, strain_id["ptt"])
        self._make_subfolder(
            paths["best"], self.without_strain, strain_id["ptt"])
        self._make_subfolder(
            paths["fig"], self.without_strain, strain_id["ptt"])
        filename_nostrain = "_".join([strain_id["file"].replace(".ptt", ""),
                                      self.without_strain + ".csv"])
        files["all_nospecific"] = open(os.path.join(paths["all"],
                                       filename_nostrain), "w")
        files["best_nospecific"] = open(os.path.join(paths["best"],
                                        filename_nostrain), "w")

    def _setup_folder_and_read_file(self, strain_id, pre_file,
                                    files, paths, args_ppi):
        if strain_id["file"].endswith(".ptt"):
            if strain_id["file"] != pre_file:
                self.helper.check_make_folder(
                     "_".join([self.tmp_id, strain_id["file"]]))
                paths["all"] = os.path.join(
                     self.all_result, strain_id["file"][:-4])
                paths["best"] = os.path.join(
                     self.best_result, strain_id["file"][:-4])
                paths["fig"] = os.path.join(
                     self.fig, strain_id["file"][:-4])
                self.helper.check_make_folder(
                     os.path.join(self.all_result, strain_id["file"][:-4]))
                self.helper.check_make_folder(
                     os.path.join(self.best_result, strain_id["file"][:-4]))
                self.helper.check_make_folder(
                     os.path.join(self.fig, strain_id["file"][:-4]))
                self._make_subfolder(
                     paths["all"], self.with_strain, strain_id["ptt"])
                self._make_subfolder(
                     paths["best"], self.with_strain, strain_id["ptt"])
                self._make_subfolder(
                     paths["fig"], self.with_strain, strain_id["ptt"])
                filename_strain = "_".join(
                     [strain_id["file"].replace(".ptt", ""),
                      self.with_strain + ".csv"])
                files["all_specific"] = open(os.path.join(
                                        paths["all"], filename_strain), "w")
                files["best_specific"] = open(os.path.join(
                                         paths["best"], filename_strain), "w")
                if args_ppi.no_specific:
                    self._setup_nospecific(paths, strain_id, files)
                files["id_list"] = "_".join([self.tmp_id, strain_id["file"]])
                files["id_log"] = open(os.path.join(files["id_list"],
                                       self.tmp_files["log"]), "w")
                files["action_log"] = open(os.path.join(args_ppi.out_folder,
                                           self.tmp_files["action"]), "w")
                files["pubmed_log"] = open(os.path.join(args_ppi.out_folder,
                                           self.tmp_files["pubmed"]), "w")
                pre_file = strain_id["file"]
                if strain_id["file"] in os.listdir(args_ppi.ptts):
                    genes = self._detect_protein(strain_id, args_ppi)
            else:
                self._make_folder_no_exist(os.path.join(paths["all"],
                                           self.with_strain), strain_id["ptt"])
                self._make_folder_no_exist(os.path.join(paths["best"],
                                           self.with_strain), strain_id["ptt"])
                if args_ppi.no_specific:
                    self._make_folder_no_exist(
                        os.path.join(paths["all"], self.without_strain),
                        strain_id["ptt"])
                    self._make_folder_no_exist(
                        os.path.join(paths["best"], self.without_strain),
                        strain_id["ptt"])
        else:
            print("Error: Wrong .ptt file!")
            sys.exit()
        return genes

    def _wget_actions(self, files, id_file, strain_id, out_folder, log):
        detect = False
        t_h = open(os.path.join(files["id_list"], id_file), "r")
        print("Retrieving STRING actions for {0} of {1} -- {2}".format(
              id_file, strain_id["string"], strain_id["file"]))
        for row in csv.reader(t_h, delimiter="\t"):
            if row[0].startswith("stringId"):
                continue
            else:
                detect = True
                if row[1] == strain_id["string"]:
                    action_source = ("http://string-db.org/api/tsv/actions?"
                                     "identifier={0}&species={1}").format(
                                     row[0], row[1])
                    self._run_wget(
                        action_source, self.tmp_files["wget_action"],
                        files["action_log"], log)
        t_h.close()
        if not detect:
            log.write(id_file + " can not be found in STRING.\n")
            print("Warning: " + id_file + " can not be found in STRING!")
        return detect

    def _retrieve_actions(self, files, strain_id, paths, args_ppi, log):
        '''get the interaction of proteins'''
        log.write("Using STRING and PIE to retrieve the interaction "
                  "information for {0}.\n".format(strain_id["ptt"]))
        for id_file in os.listdir(files["id_list"]):
            if id_file != self.tmp_files["log"]:
                detect_id = self._wget_actions(files, id_file, strain_id,
                                               args_ppi.out_folder, log)
                if detect_id:
                    a_h = open(self.tmp_files["wget_action"], "r")
                    pre_row = []
                    first = True
                    detect = False
                    first_output = {"specific_all": True,
                                    "specific_best": True,
                                    "nospecific_all": True,
                                    "nospecific_best": True}
                    print("Retrieving Pubmed IDs for {0} of {1} -- {2}".format(
                           id_file, strain_id["string"], strain_id["file"]))
                    for row_a in csv.reader(a_h, delimiter="\t"):
                        if row_a == []:
                            print("No interaction can be detected")
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
                                if (row_a[0] != pre_row[0]) or (
                                        row_a[1] != pre_row[1]):
                                    self._get_pubmed(
                                        pre_row, strain_id, mode, actor,
                                        id_file, first_output,
                                        strain_id["ptt"], files, paths,
                                        args_ppi, log)
                                    mode = row_a[2]
                                    actor = row_a[4]
                                else:
                                    mode = mode + ";" + row_a[2]
                                    actor = actor + ";" + row_a[4]
                            pre_row = row_a
                    if detect:
                        detect = False
                        self._get_pubmed(
                            row_a, strain_id, mode, actor, id_file,
                            first_output, strain_id["ptt"], files,
                            paths, args_ppi, log)
        self._list_files(args_ppi, paths, files, log)
        if detect_id:
            a_h.close()

    def _list_files(self, args_ppi, paths, files, log):
        log.write("The temporary files are generated and stored in the "
                  "following folders:\n")
        if args_ppi.no_specific:
            folders = [files["id_list"],
                       self.tmp_files["wget_action"],
                       self.tmp_files["specific"],
                       self.tmp_files["nospecific"]]
        else:
            folders = [files["id_list"],
                       self.tmp_files["wget_action"],
                       self.tmp_files["specific"]]
        for folder in folders:
            log.write("\t" + os.path.join(folder) + "\n")
        log.write("The files for storing the interaction information are "
                  "generated and stored in the following folders:\n")
        for data in (paths["all"], paths["best"]):
            for files in os.listdir(data):
                if os.path.isdir(os.path.join(data, files)):
                    for file_ in os.listdir(os.path.join(data, files)):
                        log.write("\t" + os.path.join(data, files, file_) + "\n")
        log.write("The merged tables are generated:\n")
        for data in (paths["all"], paths["best"]):
            for files in os.listdir(data):
                if os.path.isfile(os.path.join(data, files)):
                    log.write("\t" + os.path.join(data, files) + "\n")

    def _plot(self, args_ppi, files, log):
        log.write("Running plot_PPI.py to generate plots of PPI.\n")
        log.write("The figures of PPI networks are generated and stored in the "
                  "following folders:\n")
        if args_ppi.no_specific:
            files["all_nospecific"].close()
            files["best_nospecific"].close()
        files["all_specific"].close()
        files["best_specific"].close()
        for folder in os.listdir(self.all_result):
            if folder in os.listdir(self.fig):
                print("Plotting {0}".format(folder))
                out_folder_spe = os.path.join(self.fig, folder,
                                          self.with_strain)
                plot_ppi(os.path.join(self.all_result, folder,
                         "_".join([folder, self.with_strain + ".csv"])),
                         args_ppi.score, out_folder_spe, args_ppi.size)
                for file_ in os.listdir(out_folder_spe):
                    log.write("\t" + os.path.join(
                        out_folder_spe, file_) + "\n")
                if args_ppi.no_specific:
                    out_folder_nospe = os.path.join(self.fig, folder,
                                          self.without_strain)
                    plot_ppi(os.path.join(self.all_result, folder,
                             "_".join([folder, self.without_strain + ".csv"])),
                             args_ppi.score, out_folder_nospe, args_ppi.size)
                    for file_ in os.listdir(out_folder_nospe):
                        log.write("\t" + os.path.join(
                            out_folder_nospe, file_) + "\n")

    def _remove_tmps(self, args_ppi):
        self.helper.remove_all_content(os.path.join(args_ppi.out_folder),
                                       "tmp", "file")
        self.helper.remove_all_content(os.path.join(args_ppi.out_folder),
                                       "tmp", "dir")
        for file_ in os.listdir(args_ppi.ptts):
            if file_.startswith("PPI_"):
                os.remove(os.path.join(args_ppi.ptts, file_))
        self.helper.remove_all_content(os.path.join(args_ppi.out_folder),
                                       "temp", "dir")

    def check_query(self, args_ppi, log):
        if "all" not in args_ppi.querys:
            for query in args_ppi.querys:
                detect = False
                datas = query.split(":")
                for gff in os.listdir(args_ppi.ptts):
                    gff_f = open(os.path.join(args_ppi.ptts, gff), "r")
                    for entry in Gff3Parser().entries(gff_f):
                        if (entry.seq_id == datas[0]) and (
                                entry.start == int(datas[1])) and (
                                entry.end == int(datas[2])) and (
                                entry.strand == datas[3]):
                            detect = True
                            break
                    if not detect:
                        log.write(query + " is not found in gff file.\n")
                        print("Error: {0} is not found in gff file!".format(query))
                        sys.exit()

    def retrieve_ppi_network(self, args_ppi, log):
        '''retrieve PPI from STRING with PIE and draw network'''
        strain_ids = []
        paths = {}
        files = {}
        self.check_query(args_ppi, log)
        log.write("Running converter.py to generate ptt and rnt files.\n")
        log.write("The following files are generated:\n")
        for strain in args_ppi.strains:
            datas = strain.split(":")
            ptt_file = "PPI_" + datas[0].replace(".gff", ".ptt")
            rnt_file = "PPI_" + datas[0].replace(".gff", ".rnt")
            self.converter.convert_gff2rntptt(
                           os.path.join(args_ppi.ptts, datas[0]),
                           "0", os.path.join(args_ppi.ptts, ptt_file),
                           os.path.join(args_ppi.ptts, rnt_file), None, None)
            strain_ids.append({"file": ptt_file,
                               "ptt": datas[1],
                               "string": datas[2],
                               "pie": datas[3]})
            log.write("\t" + os.path.join(args_ppi.ptts, ptt_file) + "\n")
            log.write("\t" + os.path.join(args_ppi.ptts, rnt_file) + "\n")
        strain_ids.sort(key=lambda x: x["file"])
        pre_file = ""
        for strain_id in strain_ids:
            genes = self._setup_folder_and_read_file(strain_id, pre_file,
                                                     files, paths, args_ppi)
            s_h = open(args_ppi.species, "r")
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
            self._retrieve_id(strain_id, genes, files, log)
            self._retrieve_actions(files, strain_id, paths, args_ppi, log)
        self._plot(args_ppi, files, log)
        self._remove_tmps(args_ppi)
