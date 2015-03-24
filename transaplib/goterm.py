#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
import time
from subprocess import call, Popen
from sort_gff import Sort_GFF
import multiparser
import check_gff_attributes
class Go_term_finding(object):

    def _check_make_folder(self, path, folder):
        if folder in os.listdir(path):
            call(["rm", "-rf", path + folder])
        call(["mkdir", path + folder])

    def _remove_tmp(self, folder):
        if folder is not False:
            call(["rm", "-rf", folder + "/tmp"])
            os.system("rm -rf " + folder + "/*_folder")

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

    def run_Go_term(self, bin_path, gffs, database, out_folder, uniprot):
        blast_plus_path = os.environ["BLAST_Plus_HOME"]
        blast2go_path = os.environ["BLAST2GO_HOME"]
#        for gff in os.listdir(gffs):
#            if gff.endswith(".gff"):
#                check_gff_attributes._check_uni_attributes(gffs + "/" + gff)
        out_path = out_folder + "/Go_term_results/"
        multiparser._parser_gff(gffs, None)
        gff_path = gffs + "/tmp/"
        stat_path = out_folder + "/statistics/"
        prefixs = []
        ###############################################
        # retrieve Go terms from UniProt ID           #
        ###############################################
        for gff in os.listdir(gff_path):
            prefix = gff.replace(".gff", "")
            prefixs.append(prefix)
            self._check_make_folder(out_path, prefix)
            out = open(out_path + prefix + "/" + prefix + "_uniprot.csv", "w")
            print("extracting Go terms of " + prefix + " from UniProt...")
            if uniprot is False:
                call(["python", bin_path + "/retrievego.py",
                      "-d", database + "/idmapping_selected.tab",
                      "-g", gff_path + gff,
                      "-s", prefix], stdout=out)
            else:
                call(["python", bin_path + "/retrievego.py",
                      "-d", uniprot,
                      "-g", gff_path + gff,
                      "-s", prefix], stdout=out)
        #############################################
        # merge files based on gff files            #
        #############################################
        folders = []
        for folder in os.listdir(gffs):
            if folder.endswith("gff_folder"):
                folder_prefix = folder.replace(".gff_folder", "")
                self._check_make_folder(os.getcwd() + "/", folder_prefix)
                folders.append(folder_prefix)
                filenames = []
                for gff in os.listdir(gffs + "/" + folder):
                    if gff.endswith(".gff"):
                        filenames.append(gff.replace(".gff", ""))
                if len(filenames) > 1:
                    out_uni = open(folder_prefix + "/all_strains_uniprot.csv", "w")
                    for filename in filenames:
                        call(["cat", out_path + filename + "/" + filename + "_uniprot.csv"],
                                  stdout=out_uni)
                        call(["cp", out_path + filename + "/" + filename + "_uniprot.csv",
                              os.getcwd() + "/" + folder_prefix])
                else:
                    call(["cp", out_path + filenames[0] + "/" + filenames[0] + "_uniprot.csv",
                           folder_prefix + "/all_strains_uniprot.csv"])
        
        os.system("rm -rf " + out_path + "*")
        for folder in folders:
            call(["mv", folder, out_path])
        #######################################
        # statistics and figure               #
        #######################################
        for folder in os.listdir(out_path):
            self._check_make_folder(stat_path, folder)
            strain_stat_path = stat_path + folder + "/"
            call(["mkdir", strain_stat_path + "figs"])
            print("Computing statistics of " + folder)
            call(["python", bin_path + "/map_go2slim.py",
                  "-t", database + "/go.obo",
                  "-s", database + "/goslim.obo",
                  "-g", out_path + folder + "/all_strains_uniprot.csv",
                  "-r", strain_stat_path + "stat_" + folder + ".csv"])
            os.system("mv *_three_roots.png " + strain_stat_path + "figs")
            os.system("mv *_molecular_function.png " + strain_stat_path + "figs")
            os.system("mv *_cellular_component.png " + strain_stat_path + "figs")
            os.system("mv *_biological_process.png " + strain_stat_path + "figs")
        #######################################
        # remove temperary files and folders  #
        #######################################
        self._remove_tmp(fastas)
        self._remove_tmp(gffs)
        call(["rm", "-rf", out_folder + "/tmp"])
        call(["rm", "-rf", out_folder + "/tmp_go"])
