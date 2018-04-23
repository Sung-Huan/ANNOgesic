import os
import csv
import shutil
from annogesiclib.helper import Helper
from annogesiclib.multiparser import Multiparser
from annogesiclib.gene_ontology import retrieve_uniprot, map2goslim


class GoTermFinding(object):
    '''Retrieving the GO term'''

    def __init__(self, args_go):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.out_all = os.path.join(args_go.out_folder, "all_CDSs")
        self.out_express = os.path.join(args_go.out_folder, "expressed_CDSs")
        self.result_all_path = os.path.join(self.out_all, "GO_term_results")
        self.result_express_path = os.path.join(self.out_express,
                                                "GO_term_results")
        self.gff_path = os.path.join(args_go.gffs, "tmp")
        if args_go.trans is not None:
            self.tran_path = os.path.join(args_go.trans, "tmp")
        else:
            self.tran_path = None
        self.stat_all_path = os.path.join(self.out_all, "statistics")
        self.stat_express_path = os.path.join(self.out_express,
                                              "statistics")
        self.all_strain = "all_genomes_uniprot.csv"

    def _retrieve_go(self, uniprot, out_path, type_, log):
        prefixs = []
        log.write("Running gene_ontology.py to retrieve GO terms.\n")
        for gff in os.listdir(self.gff_path):
            prefix = gff.replace(".gff", "")
            prefixs.append(prefix)
            self.helper.check_make_folder(os.path.join(out_path, prefix))
            out_file = os.path.join(out_path, prefix,
                                    "_".join([prefix, "uniprot.csv"]))
            print("Extracting GO terms of {0} from UniProt".format(prefix))
            if self.tran_path is not None:
                tran_file = os.path.join(self.tran_path,
                                         "_".join([prefix, "transcript.gff"]))
            else:
                tran_file = None
            retrieve_uniprot(uniprot, os.path.join(self.gff_path, gff),
                             out_file, tran_file, type_)
            log.write("\t" + out_file + " is generated.\n")

    def _remove_header(self, out_all):
        out = open(out_all + "_tmp", "w")
        fh = open(out_all, "r")
        out.write("\t".join(["Genome", "Strand", "Start", "End",
                             "Protein_id", "Go_term"]) + "\n")
        for row in csv.reader(fh, delimiter='\t'):
            if row[0] != "Genome":
                out.write("\t".join(row) + "\n")
        out.close()
        fh.close()
        shutil.move(out_all + "_tmp", out_all)

    def _merge_files(self, gffs, out_path, out_folder, log):
        '''merge the files according to the input genome folder'''
        folders = []
        log.write("Merging the output files based on the input genome "
                  "information.\n")
        for folder in os.listdir(gffs):
            if folder.endswith("gff_folder"):
                folder_prefix = folder.replace(".gff_folder", "")
                folder_path = os.path.join(out_folder, folder_prefix)
                self.helper.check_make_folder(folder_path)
                folders.append(folder_path)
                filenames = []
                for gff in os.listdir(os.path.join(gffs, folder)):
                    if gff.endswith(".gff"):
                        filenames.append(gff.replace(".gff", ""))
                out_all = os.path.join(folder_path, self.all_strain)
                if len(filenames) > 1:
                    if self.all_strain in os.listdir(folder_path):
                        os.remove(out_all)
                    for filename in filenames:
                        csv_file = "_".join([filename, "uniprot.csv"])
                        self.helper.merge_file(os.path.join(out_path,
                                               filename, csv_file), out_all)
                        self._remove_header(out_all)
                        shutil.copy(os.path.join(out_path, filename, csv_file),
                                    folder_path)
                else:
                    shutil.copyfile(os.path.join(out_path, filenames[0],
                                    "_".join([filenames[0], "uniprot.csv"])),
                                    out_all)
        self.helper.remove_all_content(out_path, None, "dir")
        self.helper.remove_all_content(out_path, None, "file")
        for folder in folders:
            folder_prefix = folder.split("/")[-1]
            shutil.move(folder, os.path.join(out_path, folder_prefix))
            for file_ in os.listdir(os.path.join(out_path, folder_prefix)):
                log.write("\t" + os.path.join(out_path, folder_prefix, file_) + 
                          " is generated.\n")

    def _stat(self, out_path, stat_path, go, goslim, out_folder, log):
        log.write("Running gene_ontology.py to Retrieve GOslim terms and "
                  "do statistics.\n")
        log.write("The following files are generated:\n")
        for folder in os.listdir(out_path):
            strain_stat_path = os.path.join(stat_path, folder)
            self.helper.check_make_folder(strain_stat_path)
            fig_path = os.path.join(strain_stat_path, "figs")
            if "fig" not in os.listdir(strain_stat_path):
                os.mkdir(fig_path)
            stat_file = os.path.join(strain_stat_path,
                                    "_".join(["stat", folder + ".csv"]))
            map2goslim(goslim, go,
                       os.path.join(out_path, folder, self.all_strain),
                       stat_file, out_folder)
            log.write("\t" + stat_file + "\n")
            self.helper.move_all_content(out_folder, fig_path,
                                         ["_three_roots.png"])
            self.helper.move_all_content(out_folder, fig_path,
                                         ["_molecular_function.png"])
            self.helper.move_all_content(out_folder, fig_path,
                                         ["_cellular_component.png"])
            self.helper.move_all_content(out_folder, fig_path,
                                         ["_biological_process.png"])
            for file_ in os.listdir(fig_path):
                log.write("\t" + os.path.join(fig_path, file_) + "\n")

    def run_go_term(self, args_go, log):
        for gff in os.listdir(args_go.gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(
                                                 args_go.gffs, gff))
        self.multiparser.parser_gff(args_go.gffs, None)
        if args_go.trans is not None:
            self.multiparser.parser_gff(args_go.trans, "transcript")
        print("Computing all CDSs")
        log.write("Retrieving GO terms for all CDSs.\n")
        self._retrieve_go(args_go.uniprot, self.result_all_path, "all", log)
        self._merge_files(args_go.gffs, self.result_all_path, self.out_all, log)
        self._stat(self.result_all_path, self.stat_all_path, args_go.go,
                   args_go.goslim, self.out_all, log)
        if args_go.trans is not None:
            log.write("Retrieving GO terms only for expressed CDSs.\n")
            print("Computing express CDSs")
            self._retrieve_go(args_go.uniprot, self.result_express_path,
                              "express", log)
            self._merge_files(args_go.gffs, self.result_express_path,
                              self.out_express, log)
            self._stat(self.result_express_path, self.stat_express_path,
                       args_go.go, args_go.goslim, self.out_express, log)
        self.helper.remove_tmp_dir(args_go.gffs)
        if args_go.trans is not None:
            self.helper.remove_tmp_dir(args_go.trans)
