#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
from subprocess import call, Popen
from sort_gff import Sort_GFF
import multiparser
import gff3
import check_gff_attributes
class Ribos(object):

    def _check_make_folder(self, path, folder):
        if folder in os.listdir(path):
            call(["rm", "-rf", path + folder])
        call(["mkdir", path + folder])

    def _remove_tmp(self, folder):
        if folder is not False:
            call(["rm", "-rf", folder + "/tmp"])
            os.system("rm -rf " + folder + "/*_folder")

    def run_ribos(self, bin_path, ribos_id, gffs, fastas, Rfam, 
                  out_folder, re_scan, database, fuzzy):
        infernal_path = os.environ["INFERNAL_HOME"]
        multiparser._parser_gff(gffs, None)
        multiparser._parser_fasta(fastas)
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                check_gff_attributes._check_uni_attributes(gffs + "/" + gff)
        gff_path = gffs + "/tmp/"
        fasta_path = fastas + "/tmp/"
        stat_folder = out_folder + "/statistics/"
        gff_outfolder = out_folder + "/gffs/"
        table_folder = out_folder + "/tables/"
        scan_folder = out_folder + "/scan_Rfam/"
        #### get the data of riboswitch from Rfam. ####
        ribos_rfam = open(database + "/Rfam_riboswitch.cm", "w")
        call(["python", bin_path + "/get_Rfam_ribo.py",
              "-t", ribos_id, "-r", Rfam], stdout=ribos_rfam)
        print("compressing Rfam...")
        call([infernal_path + "/cmpress", "-F", database + "/Rfam_riboswitch.cm"])
        #### extract seq for scan. ####
        prefixs = []
        self._check_make_folder(out_folder + "/", "tmp_fasta")
        self._check_make_folder(out_folder + "/", "tmp_scan")
        self._check_make_folder(out_folder + "/", "tmp_table")
        seq_path = out_folder + "/tmp_fasta/"
        gff_parser = gff3.Gff3Parser()
        for gff in os.listdir(gff_path):
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                first_seq = open(seq_path + prefix + ".fa", "w")
                prefixs.append(prefix)
                print("extracting seq of riboswitch candidates of " + prefix)
                call(["python", bin_path + "/extract_RBS.py",
                      "-s", fasta_path + prefix + ".fa",
                      "-g", gff_path + gff], stdout=first_seq)
                first_seq.close()
                #### scan Rfam. ####
                first_scan = open(out_folder + "/tmp_scan/" + prefix + "_RBS.txt", "w")
                print("scanning Rfam for " + prefix)
                call([infernal_path + "/cmscan", "--acc", 
                      database + "/Rfam_riboswitch.cm", 
                      seq_path + prefix + ".fa"], stdout=first_scan)
                first_scan.close()
                #### generate the table of the first scan results ####
                sec_seq = open(seq_path + prefix + "_regenerate.fa", "w")
                call(["python", bin_path + "/regenerate_seq.py",
                      "-i", out_folder + "/tmp_scan/" + prefix + "_RBS.txt",
                      "-s", seq_path + prefix + ".fa",
                      "-o", out_folder + "/tmp_table/" + prefix + "_RBS.csv"], stdout=sec_seq)
                sec_seq.close()
                #### if it doesn't run re_scan, delete the regenerate seq ####
                if re_scan is False:
                    os.system("rm " + seq_path + prefix + "_regenerate.fa")
                #### re-scanning ####
                if re_scan is True:
                    print("re-scanning of " + prefix)
                    sec_scan = open(out_folder + "/tmp_scan/" + prefix + "_RBS_rescan.txt", "w")
                    call([infernal_path + "/cmscan", "--acc",
                          database + "/Rfam_riboswitch.cm",
                          seq_path + prefix + "_regenerate.fa"], stdout=sec_scan)
                    sec_scan.close()
                    call(["python", bin_path + "/reextract_CM.py",
                          "-i", out_folder + "/tmp_scan/" + prefix + "_RBS_rescan.txt",
                          "-f", out_folder + "/tmp_table/" + prefix + "_RBS.csv",
                          "-o", out_folder + "/tmp_table/" + prefix + "_RBS_rescan.csv"])
                    call(["mv", out_folder + "/tmp_table/" + prefix + "_RBS_rescan.csv", 
                                out_folder + "/tmp_table/" + prefix + "_RBS.csv"])
        #### merge the results based on the gff annotation files ####
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                print("Merge results of " + prefix)
                pre_strain = ""
                self._check_make_folder(scan_folder, prefix)
                for entry in gff_parser.entries(open(gffs + "/" + gff)):
                    if entry.seq_id != pre_strain:
                        if len(pre_strain) == 0:
                            os.system("cat " + out_folder + "/tmp_table/" + entry.seq_id + "_RBS.csv > " + \
                                      table_folder + prefix + "_RBS.csv")
                        else:
                            os.system("cat " + out_folder +"/tmp_table/" + entry.seq_id + "_RBS.csv >> " + \
                                      table_folder + prefix + "_RBS.csv")
                        call(["cp", out_folder + "/tmp_scan/" + entry.seq_id + "_RBS.txt", scan_folder + prefix])
                        if re_scan:
                            call(["cp", out_folder + "/tmp_scan/" + entry.seq_id + "_RBS_rescan.txt", scan_folder + prefix])
                        pre_strain = entry.seq_id
                #### covert to gff and do statistics ####
                out_stat = open(stat_folder + "stat_" + prefix + "_RBS.txt", "w")
                print("compute statistics of " + prefix)
                call(["python", bin_path + "/ribo_gff.py",
                      "-i", table_folder + prefix + "_RBS.csv",
                      "-r", ribos_id, "-f", str(fuzzy),
                      "-g", gff_outfolder + prefix + "_RBS.gff"], stdout=out_stat)
        #### delete temperary files ####
        self._remove_tmp(gffs)
        self._remove_tmp(fastas)
        os.system("rm -rf " + out_folder + "/tmp_scan")
        os.system("rm -rf " + out_folder + "/tmp_table")
        os.system("rm -rf " + out_folder + "/tmp_fasta")
