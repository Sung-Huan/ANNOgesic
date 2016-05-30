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
from annogesiclib.map_ribos import mapping_ribos
from annogesiclib.rbs_overlap import rbs_overlap


class Ribos(object):

    def __init__(self, args_ribo):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.gff_parser = Gff3Parser()
        self.gff_path = os.path.join(args_ribo.gffs, "tmp")
        self.tss_path = os.path.join(args_ribo.tsss, "tmp")
        self.tran_path = os.path.join(args_ribo.trans, "tmp")
        self.fasta_path = os.path.join(args_ribo.fastas, "tmp")
        self.stat_folder = os.path.join(args_ribo.out_folder, "statistics")
        self.gff_outfolder = os.path.join(args_ribo.out_folder, "gffs")
        self.table_folder = os.path.join(args_ribo.out_folder, "tables")
        self.scan_folder = os.path.join(args_ribo.out_folder, "scan_Rfam")
        self.ribos_rfam = os.path.join(args_ribo.database,
                                       "Rfam_riboswitch.cm")
        self.tmp_files = {"fasta": os.path.join(
                                   args_ribo.out_folder, "tmp_fasta"),
                          "scan": os.path.join(
                                  args_ribo.out_folder, "tmp_scan"),
                          "table": os.path.join(
                                   args_ribo.out_folder, "tmp_table")}
        self.suffixs = {"csv": "riboswitch.csv",
                        "txt": "riboswitch_prescan.txt",
                        "re_txt": "riboswitch_scan.txt",
                        "re_csv": "riboswitch_scan.csv"}

    def _run_infernal(self, args_ribo, seq, type_, prefix):
        scan_file = os.path.join(self.tmp_files["scan"],
                                 "_".join([prefix, self.suffixs[type_]]))
        scan = open(scan_file, "w")
        call([os.path.join(args_ribo.infernal_path, "cmscan"), "--incE",
              str(args_ribo.e_value), "--acc", self.ribos_rfam, seq],
             stdout=scan)
        scan.close()
        return scan_file

    def _scan_extract_rfam(self, prefixs, args_ribo):
        for gff in os.listdir(self.gff_path):
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                first_seq = os.path.join(self.tmp_files["fasta"],
                                         prefix + ".fa")
                prefixs.append(prefix)
                print("extracting seq of riboswitch candidates of {0}".format(
                      prefix))
                extract_potential_rbs(
                      os.path.join(self.fasta_path, prefix + ".fa"),
                      os.path.join(self.gff_path, gff),
                      os.path.join(self.tss_path, prefix + "_TSS.gff"),
                      os.path.join(self.tran_path, prefix + "_transcript.gff"),
                      first_seq, args_ribo)
                print("pre-scanning of {0}".format(prefix))
                first_scan_file = self._run_infernal(args_ribo, first_seq,
                                                     "txt", prefix)
                sec_seq = os.path.join(self.tmp_files["fasta"],
                                       "_".join([prefix, "regenerate.fa"]))
                first_table = os.path.join(
                        self.tmp_files["table"],
                        "_".join([prefix, self.suffixs["csv"]]))
                regenerate_seq(first_scan_file, first_seq,
                               first_table, sec_seq)
                print("scanning of {0}".format(prefix))
                sec_scan_file = self._run_infernal(args_ribo, sec_seq,
                                                   "re_txt", prefix)
                sec_table = os.path.join(
                        self.tmp_files["table"],
                        "_".join([prefix, self.suffixs["re_csv"]]))
                reextract_rbs(sec_scan_file, first_table, sec_table)
                shutil.move(sec_table, first_table)
                modify_table(first_table, args_ribo.output_all)
        return prefixs

    def _merge_results(self, args_ribo):
        for gff in os.listdir(args_ribo.gffs):
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                print("Merge results of {0}".format(prefix))
                pre_strain = ""
                self.helper.check_make_folder(os.path.join(
                                              self.scan_folder, prefix))
                fh = open(os.path.join(args_ribo.gffs, gff))
                for entry in self.gff_parser.entries(fh):
                    if entry.seq_id != pre_strain:
                        if len(pre_strain) == 0:
                            shutil.copyfile(os.path.join(
                                self.tmp_files["table"],
                                "_".join([entry.seq_id, self.suffixs["csv"]])),
                                os.path.join(
                                    self.table_folder,
                                    "_".join([prefix, self.suffixs["csv"]])))
                        else:
                            self.helper.merge_file(os.path.join(
                                self.tmp_files["table"],
                                "_".join([entry.seq_id, self.suffixs["csv"]])),
                                os.path.join(
                                    self.table_folder,
                                    "_".join([prefix, self.suffixs["csv"]])))
                        shutil.copy(os.path.join(
                            self.tmp_files["scan"],
                            "_".join([entry.seq_id, self.suffixs["txt"]])),
                            os.path.join(self.scan_folder, prefix))
                        shutil.copy(os.path.join(
                            self.tmp_files["scan"],
                            "_".join([entry.seq_id, self.suffixs["re_txt"]])),
                            os.path.join(self.scan_folder, prefix))
                        pre_strain = entry.seq_id
                out_stat = os.path.join(
                        self.stat_folder,
                        "_".join(["stat", prefix, "riboswitch.txt"]))
                print("compute statistics of {0}".format(prefix))
                stat_and_covert2gff(os.path.join(
                    self.table_folder,
                    "_".join([prefix, self.suffixs["csv"]])),
                    args_ribo.ribos_id, os.path.join(
                        self.gff_outfolder,
                        "_".join([prefix, "riboswitch.gff"])),
                    args_ribo.fuzzy, out_stat)
                fh.close()

    def _remove_tmp(self, args_ribo):
        self.helper.remove_tmp(args_ribo.gffs)
        self.helper.remove_tmp(args_ribo.fastas)
        self.helper.remove_all_content(args_ribo.out_folder, "tmp", "dir")

    def _remove_overlap(self, gff_path):
        for gff in os.listdir(gff_path):
            if gff.endswith(".gff"):
                rbs_overlap(
                    os.path.join(os.path.join(
                        self.tmp_files["table"],
                        "_".join([gff.replace(".gff", ""),
                                  self.suffixs["csv"]]))),
                    os.path.join(gff_path, gff))

    def run_ribos(self, args_ribo):
        if args_ribo.fuzzy_rbs > 6:
            print("Error: --fuzzy_rbs should be equal or less than 6!!")
            sys.exit()
        self.multiparser.parser_gff(args_ribo.gffs, None)
        self.multiparser.parser_fasta(args_ribo.fastas)
        self.multiparser.parser_gff(args_ribo.trans, "transcript")
        self.multiparser.parser_gff(args_ribo.tsss, "TSS")
        for gff in os.listdir(args_ribo.gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(
                                                 args_ribo.gffs, gff))
        rbs_from_rfam(args_ribo.ribos_id, args_ribo.rfam, self.ribos_rfam)
        print("compressing Rfam...")
        call([os.path.join(args_ribo.infernal_path, "cmpress"),
              "-F", self.ribos_rfam])
        prefixs = []
        self.helper.check_make_folder(self.tmp_files["fasta"])
        self.helper.check_make_folder(self.tmp_files["scan"])
        self.helper.check_make_folder(self.tmp_files["table"])
        prefixs = self._scan_extract_rfam(prefixs, args_ribo)
        self._remove_overlap(self.gff_path)
        self._merge_results(args_ribo)
        mapping_ribos(self.table_folder, args_ribo.ribos_id)
        self._remove_tmp(args_ribo)
