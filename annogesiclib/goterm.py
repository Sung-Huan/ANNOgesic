#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
import time
import shutil
from subprocess import call, Popen
from annogesiclib.helper import Helper
from annogesiclib.multiparser import Multiparser
from annogesiclib.gene_ontology import retrieve_uniprot, map2goslim


class Go_term_finding(object):

    def __init__(self):
        self.multiparser = Multiparser()
        self.helper = Helper()   

    def _wait_process(self, processes):
        for p in processes:
            p.wait()
        if p.stdout:
            p.stdout.close()
        if p.stdin:
            p.stdin.close()
        if p.stderr:
            p.stderr.close()
        try:
            p.kill()
        except OSError:
            pass
        time.sleep(5)

    def _retrieve_go(self, gff_path, out_path, uniprot):
        prefixs = []
        for gff in os.listdir(gff_path):
            prefix = gff.replace(".gff", "")
            prefixs.append(prefix)
            self.helper.check_make_folder(out_path, prefix)
            out_file = os.path.join(out_path, prefix, "_".join([prefix, "uniprot.csv"]))
            print("extracting Go terms of {0} from UniProt...".format(prefix))
            retrieve_Uniprot(uniprot, os.path.join(gff_path, gff), out_file)
        return prefixs

    def _merge_files(self, gffs, out_path, out_folder):
        folders = []
        for folder in os.listdir(gffs):
            if folder.endswith("gff_folder"):
                folder_prefix = folder.replace(".gff_folder", "")
                self.helper.check_make_folder(out_folder, folder_prefix)
                folders.append(os.path.join(out_folder, folder_prefix))
                filenames = []
                for gff in os.listdir(os.path.join(gffs, folder)):
                    if gff.endswith(".gff"):
                        filenames.append(gff.replace(".gff", ""))
                out_all = os.path.join(out_folder, folder_prefix, "all_strains_uniprot.csv")
                if len(filenames) > 1:
                    if "all_strains_uniprot.csv" in os.listdir(os.path.join(out_folder, folder_prefix)):
                        os.remove(out_all)
                    for filename in filenames:
                        self.helper.merge_file(os.path.join(out_path, filename), 
                                               "_".join([filename, "uniprot.csv"]),
                                               os.path.join(out_folder, folder_prefix),
                                               "all_strains_uniprot.csv")
                        shutil.copy(os.path.join(out_path, filename, 
                                    "_".join([filename, "uniprot.csv"])),
                                    os.path.join(out_folder, folder_prefix))
                else:
                    shutil.copyfile(os.path.join(out_path, filenames[0], 
                                    "_".join([filenames[0], "uniprot.csv"])),
                                    out_all)
        self.helper.remove_all_content(out_path, None, "dir")
        for folder in folders:
            folder_prefix = folder.split("/")[-1]
            os.rename(folder, os.path.join(out_path, folder_prefix))

    def _stat(self, out_path, stat_path, go, goslim, out_folder):
        for folder in os.listdir(out_path):
            self.helper.check_make_folder(stat_path, folder)
            strain_stat_path = os.path.join(stat_path, folder)
            if "fig" not in os.listdir(strain_stat_path):
                os.mkdir(os.path.join(strain_stat_path, "figs"))
            print("Computing statistics of {0}".format(folder))
            map2goslim(goslim, go, os.path.join(out_path, folder, "all_strains_uniprot.csv"),
                       os.path.join(strain_stat_path, "_".join(["stat", folder + ".csv"])),
                       out_folder)
            self.helper.move_all_content(out_folder, 
                        os.path.join(strain_stat_path, "figs"), "_three_roots.png")
            self.helper.move_all_content(out_folder, 
                        os.path.join(strain_stat_path, "figs"), "_molecular_function.png")
            self.helper.move_all_content(out_folder, 
                        os.path.join(strain_stat_path, "figs"), "_cellular_component.png")
            self.helper.move_all_content(out_folder,
                        os.path.join(strain_stat_path, "figs"), "_biological_process.png")

    def run_go_term(self, bin_path, gffs, out_folder, uniprot, go, goslim):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(gffs, gff))
        out_path = os.path.join(out_folder, "Go_term_results")
        self.multiparser._parser_gff(gffs, None)
        gff_path = os.path.join(gffs, "tmp")
        stat_path = os.path.join(out_folder, "statistics")
        prefixs = self._retrieve_go(gff_path, out_path, uniprot)
        self._merge_files(gffs, out_path, out_folder) ## merge files based on gff file
        self._stat(out_path, stat_path, go, goslim, out_folder)
        self.helper.remove_tmp(gffs)
