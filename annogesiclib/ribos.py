import os
import sys
import shutil
from subprocess import call
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.get_Rfam_ribo import rbs_from_rfam
from annogesiclib.extract_RBS import extract_potential_rbs
from annogesiclib.recompute_RBS import regenerate_seq, reextract_rbs
from annogesiclib.ribo_gff import stat_and_covert2gff
from annogesiclib.modify_rbs_table import modify_table


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
        self.tmp_files = {"fasta": os.path.join(out_folder, "tmp_fasta"),
                          "scan": os.path.join(out_folder, "tmp_scan"),
                          "table": os.path.join(out_folder, "tmp_table")}
        self.suffixs = {"csv": "riboswitch.csv", "txt": "riboswitch.txt",
                        "re_txt": "riboswitch_rescan.txt", "re_csv": "riboswitch_rescan.csv"}

    def _run_infernal(self, infernal_path, e_value, seq, type_, prefix):
        scan_file = os.path.join(self.tmp_files["scan"],
                             "_".join([prefix, self.suffixs[type_]]))
        scan = open(scan_file, "w")
        call([os.path.join(infernal_path, "cmscan"), "--incE", str(e_value), "--acc",
              self.ribos_rfam, seq], stdout=scan)
        scan.close()
        return scan_file

    def _scan_extract_rfam(self, gff_path, e_value, seq_path, prefixs, fasta_path,
                           out_folder, infernal_path, re_scan, output_all,
                           start_codons, start_rbs, end_rbs, fuzzy_rbs):
        for gff in os.listdir(gff_path):
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                first_seq = os.path.join(seq_path, prefix + ".fa")
                prefixs.append(prefix)
                print("extracting seq of riboswitch candidates of {0}".format(prefix))
                extract_potential_rbs(os.path.join(fasta_path, prefix + ".fa"),
                                      os.path.join(gff_path, gff), first_seq,
                                      start_codons, start_rbs, end_rbs, fuzzy_rbs)
                print("scanning Rfam for {0}".format(prefix))
                first_scan_file = self._run_infernal(infernal_path, e_value,
                                                     first_seq, "txt", prefix)
                sec_seq = os.path.join(seq_path,
                          "_".join([prefix, "regenerate.fa"]))
                first_table = os.path.join(self.tmp_files["table"],
                              "_".join([prefix, self.suffixs["csv"]]))
                #### the table of first scan results and seq for second time
                regenerate_seq(first_scan_file, first_seq, first_table, sec_seq)
                if not re_scan: #### if no re_scan, delete the regenerate seq
                    modify_table(first_table, output_all)
                    os.remove(sec_seq)
                if re_scan: #### re-scanning
                    print("re-scanning of {0}".format(prefix))
                    sec_scan_file = self._run_infernal(infernal_path, e_value,
                                                       sec_seq, "re_txt", prefix)
                    sec_table = os.path.join(self.tmp_files["table"],
                                "_".join([prefix, self.suffixs["re_csv"]]))
                    reextract_rbs(sec_scan_file, first_table, sec_table)
                    shutil.move(sec_table, first_table)
                    modify_table(first_table, output_all)
        return prefixs

    def _merge_results(self, gffs, scan_folder, table_folder,
                       stat_folder, ribos_id, fuzzy, gff_outfolder, re_scan):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                print("Merge results of {0}".format(prefix))
                pre_strain = ""
                self.helper.check_make_folder(os.path.join(scan_folder, prefix))
                fh = open(os.path.join(gffs, gff))
                for entry in self.gff_parser.entries(fh):
                    if entry.seq_id != pre_strain:
                        if len(pre_strain) == 0:
                            shutil.copyfile(
                                   os.path.join(self.tmp_files["table"],
                                   "_".join([entry.seq_id, self.suffixs["csv"]])),
                                   os.path.join(table_folder,
                                   "_".join([prefix, self.suffixs["csv"]])))
                        else:
                            self.helper.merge_file(
                                 os.path.join(self.tmp_files["table"],
                                 "_".join([entry.seq_id, self.suffixs["csv"]])),
                                 os.path.join(table_folder,
                                 "_".join([prefix, self.suffixs["csv"]])))
                        shutil.copy(os.path.join(self.tmp_files["scan"],
                                    "_".join([entry.seq_id, self.suffixs["txt"]])),
                                    os.path.join(scan_folder, prefix))
                        if re_scan:
                            shutil.copy(os.path.join(self.tmp_files["scan"],
                                   "_".join([entry.seq_id, self.suffixs["re_txt"]])),
                                   os.path.join(scan_folder, prefix))
                        pre_strain = entry.seq_id
                out_stat = os.path.join(stat_folder,
                           "_".join(["stat", prefix, self.suffixs["txt"]]))
                print("compute statistics of {0}".format(prefix))
                stat_and_covert2gff(os.path.join(
                     table_folder, "_".join([prefix, self.suffixs["csv"]])),
                     ribos_id,
                     os.path.join(gff_outfolder, "_".join([prefix, "riboswitch.gff"])),
                     fuzzy, out_stat)
                fh.close()

    def _remove_tmp(self, gffs, fastas, out_folder):
        self.helper.remove_tmp(gffs)
        self.helper.remove_tmp(fastas)
        self.helper.remove_all_content(out_folder, "tmp", "dir")

    def run_ribos(self, infernal_path, ribos_id, gffs, fastas, rfam,
                  out_folder, re_scan, e_value, output_all, database, fuzzy,
                  start_codons, start_rbs, end_rbs, fuzzy_rbs):
        if fuzzy_rbs > 6:
            print("Error: --fuzzy_rbs should be equal or less than 6!!")
            sys.exit()
        self.multiparser.parser_gff(gffs, None)
        self.multiparser.parser_fasta(fastas)
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(gffs, gff))
        #### get the data of riboswitch from Rfam.
        rbs_from_rfam(ribos_id, rfam, self.ribos_rfam)
        print("compressing Rfam...")
        call([os.path.join(infernal_path, "cmpress"), "-F", self.ribos_rfam])
        prefixs = []
        self.helper.check_make_folder(self.tmp_files["fasta"])
        self.helper.check_make_folder(self.tmp_files["scan"])
        self.helper.check_make_folder(self.tmp_files["table"])
        prefixs = self._scan_extract_rfam(self.gff_path, e_value,
                  self.tmp_files["fasta"], prefixs, self.fasta_path,
                  out_folder, infernal_path, re_scan, output_all,
                  start_codons, start_rbs, end_rbs, fuzzy_rbs)
        #### merge the results based on annotation files
        self._merge_results(gffs, self.scan_folder, self.table_folder,
                            self.stat_folder, ribos_id, fuzzy,
                            self.gff_outfolder, re_scan)
        self._remove_tmp(gffs, fastas, out_folder)
