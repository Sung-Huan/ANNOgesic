import os
import shutil
from annogesiclib.helper import Helper
from annogesiclib.multiparser import Multiparser
from annogesiclib.gene_ontology import retrieve_uniprot, map2goslim


class GoTermFinding(object):

    def __init__(self, out_folder, gffs):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.out_path = os.path.join(out_folder, "Go_term_results")
        self.gff_path = os.path.join(gffs, "tmp")
        self.stat_path = os.path.join(out_folder, "statistics")
        self.all_strain = "all_strains_uniprot.csv"

    def _retrieve_go(self, gff_path, out_path, uniprot):
        prefixs = []
        for gff in os.listdir(gff_path):
            prefix = gff.replace(".gff", "")
            prefixs.append(prefix)
            self.helper.check_make_folder(os.path.join(out_path, prefix))
            out_file = os.path.join(out_path, prefix,
                                    "_".join([prefix, "uniprot.csv"]))
            print("extracting Go terms of {0} from UniProt...".format(prefix))
            retrieve_uniprot(uniprot, os.path.join(gff_path, gff), out_file)

    def _merge_files(self, gffs, out_path, out_folder):
        folders = []
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
            os.rename(folder, os.path.join(out_path, folder_prefix))

    def _stat(self, out_path, stat_path, go, goslim, out_folder):
        for folder in os.listdir(out_path):
            strain_stat_path = os.path.join(stat_path, folder)
            self.helper.check_make_folder(strain_stat_path)
            fig_path = os.path.join(strain_stat_path, "figs")
            if "fig" not in os.listdir(strain_stat_path):
                os.mkdir(fig_path)
            print("Computing statistics of {0}".format(folder))
            map2goslim(goslim, go,
                       os.path.join(out_path, folder, self.all_strain),
                       os.path.join(strain_stat_path,
                                    "_".join(["stat", folder + ".csv"])),
                       out_folder)
            self.helper.move_all_content(out_folder, fig_path,
                                         ["_three_roots.png"])
            self.helper.move_all_content(out_folder, fig_path,
                                         ["_molecular_function.png"])
            self.helper.move_all_content(out_folder, fig_path,
                                         ["_cellular_component.png"])
            self.helper.move_all_content(out_folder, fig_path,
                                         ["_biological_process.png"])

    def run_go_term(self, gffs, out_folder, uniprot, go, goslim):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(gffs, gff))
        self.multiparser.parser_gff(gffs, None)
        self._retrieve_go(self.gff_path, self.out_path, uniprot)
        self._merge_files(gffs, self.out_path, out_folder)
        self._stat(self.out_path, self.stat_path, go, goslim, out_folder)
        self.helper.remove_tmp(gffs)
