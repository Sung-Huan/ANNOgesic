#!/usr/bin/python
import os	
import sys
import shutil
from subprocess import call, Popen
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.get_Rfam_ribo import rbs_from_rfam
from annogesiclib.extract_RBS import extract_potential_rbs
from annogesiclib.recompute_RBS import regenerate_seq, reextract_rbs
from annogesiclib.ribo_gff import stat_and_covert2gff


class Ribos(object):

    def __init__(self, gffs, fastas, out_folder, database):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.gff_parser = Gff3Parser()
        self.gff_path = os.path.join(gffs, "tmp")
        self.fasta_path = os.path.join(fastas, "tmp")
        self.stat_folder = os.path.join(out_folder, "statistics")
        self.gff_outfolder = os.path.join(out_folder, "gffs")
        self.table_folder = os.path.join(out_folder, "tables")
        self.scan_folder = os.path.join(out_folder, "scan_Rfam")
        self.ribos_rfam = os.path.join(database, "Rfam_riboswitch.cm")
        self.tmp_fasta = os.path.join(out_folder, "tmp_fasta")
        self.tmp_scan = os.path.join(out_folder, "tmp_scan")
        self.tmp_table = os.path.join(out_folder, "tmp_table")
        self.endfix_csv = "RBS.csv"
        self.endfix_txt = "RBS.txt"
        self.endfix_rescan_txt = "RBS_rescan.txt"
        self.endfix_rescan_csv = "RBS_rescan.csv"

    def _scan_extract_Rfam(self, gff_path, seq_path, prefixs, fasta_path, out_folder,
                           infernal_path, database, re_scan):
        for gff in os.listdir(gff_path):
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                first_seq = os.path.join(seq_path, prefix + ".fa")
                prefixs.append(prefix)
                print("extracting seq of riboswitch candidates of {0}".format(prefix))
                extract_potential_rbs(os.path.join(fasta_path, prefix + ".fa"), 
                                      os.path.join(gff_path, gff), first_seq) ### extract seq for scan
                first_scan_file = os.path.join(self.tmp_scan,
                                  "_".join([prefix, self.endfix_txt]))
                first_scan = open(first_scan_file, "w")
                print("scanning Rfam for {0}".format(prefix))
                call([os.path.join(infernal_path, "cmscan"), "--acc",
                      self.ribos_rfam, first_seq], stdout=first_scan) ### scan Rfam
                first_scan.close()
                sec_seq = os.path.join(seq_path, "_".join([prefix, "regenerate.fa"]))
                first_table = os.path.join(self.tmp_table, 
                                           "_".join([prefix, self.endfix_csv]))
                #### generate the table of the first scan results and seq for scan second time
                regenerate_seq(first_scan_file, first_seq, first_table, sec_seq)
                if re_scan is False: #### if it doesn't run re_scan, delete the regenerate seq
                    os.remove(sec_seq)
                if re_scan is True: #### re-scanning
                    print("re-scanning of {0}".format(prefix))
                    sec_scan_file = os.path.join(self.tmp_scan,
                                    "_".join([prefix, self.endfix_rescan_txt]))
                    sec_scan = open(sec_scan_file, "w")
                    call([os.path.join(infernal_path, "cmscan"), "--acc",
                          self.ribos_rfam, sec_seq], stdout=sec_scan)
                    sec_scan.close()
                    sec_table = os.path.join(self.tmp_table,
                                           "_".join([prefix, self.endfix_rescan_csv]))
                    reextract_rbs(sec_scan_file, first_table, sec_table)
                    os.rename(sec_table, first_table)
        return prefixs

    def _merge_results(self, gffs, scan_folder, out_folder, table_folder, stat_folder,
                       ribos_id, fuzzy, gff_outfolder, re_scan):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                print("Merge results of {0}".format(prefix))
                pre_strain = ""
                self.helper.check_make_folder(os.path.join(scan_folder, prefix))
                for entry in self.gff_parser.entries(open(os.path.join(gffs, gff))):
                    if entry.seq_id != pre_strain:
                        if len(pre_strain) == 0:
                            shutil.copyfile(os.path.join(self.tmp_table, 
                                            "_".join([entry.seq_id, self.endfix_csv])),
                                            os.path.join(table_folder, 
                                            "_".join([prefix, self.endfix_csv])))
                        else:
                            self.helper.merge_file(os.path.join(self.tmp_table,
                                                   "_".join([entry.seq_id, self.endfix_csv])),
                                                   os.path.join(table_folder, 
                                                   "_".join([prefix, self.endfix_csv])))
                        shutil.copy(os.path.join(self.tmp_scan, 
                                    "_".join([entry.seq_id, self.endfix_txt])), 
                                    os.path.join(scan_folder, prefix))
                        if re_scan:
                            shutil.copy(os.path.join(self.tmp_scan, 
                                        "_".join([entry.seq_id, self.endfix_rescan_txt])), 
                                        os.path.join(scan_folder, prefix))
                        pre_strain = entry.seq_id
                out_stat = os.path.join(stat_folder, "_".join(["stat", prefix, self.endfix_txt]))
                print("compute statistics of {0}".format(prefix)) #### covert to gff and do statistics
                stat_and_covert2gff(os.path.join(table_folder, "_".join([prefix, self.endfix_csv])), 
                                    ribos_id, 
                                    os.path.join(gff_outfolder, "_".join([prefix, "RBS.gff"])), 
                                    fuzzy, out_stat)

    def _remove_tmp(self, gffs, fastas, out_folder):
        self.helper.remove_tmp(gffs)
        self.helper.remove_tmp(fastas)
        self.helper.remove_all_content(out_folder, "tmp", "dir")

    def run_ribos(self, infernal_path, ribos_id, gffs, fastas, Rfam, 
                  out_folder, re_scan, database, fuzzy):
        self.multiparser._parser_gff(gffs, None)
        self.multiparser._parser_fasta(fastas)
#        for gff in os.listdir(gffs):
#            if gff.endswith(".gff"):
#                self.helper.check_uni_attributes(os.path.join(gffs, gff))
        rbs_from_rfam(ribos_id, Rfam, self.ribos_rfam) #### get the data of riboswitch from Rfam.
        print("compressing Rfam...")
        call([os.path.join(infernal_path, "cmpress"), "-F", self.ribos_rfam])
        prefixs = []
        self.helper.check_make_folder(self.tmp_fasta)
        self.helper.check_make_folder(self.tmp_scan)
        self.helper.check_make_folder(self.tmp_table)
        prefixs = self._scan_extract_Rfam(self.gff_path, self.tmp_fasta, prefixs, self.fasta_path, 
                                          out_folder, infernal_path, database, re_scan)
        self._merge_results(gffs, self.scan_folder, out_folder, self.table_folder, self.stat_folder,
                            ribos_id, fuzzy, self.gff_outfolder, re_scan)  #### merge the results based on annotation files
        self._remove_tmp(gffs, fastas, out_folder)
