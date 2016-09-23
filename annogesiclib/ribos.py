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
    '''detection of riboswitch and RNA thermometer'''

    def __init__(self, args_ribo):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.gff_parser = Gff3Parser()
        self.gff_path = os.path.join(args_ribo.gffs, "tmp")
        self.tss_path = os.path.join(args_ribo.tsss, "tmp")
        self.tran_path = os.path.join(args_ribo.trans, "tmp")
        self.fasta_path = os.path.join(args_ribo.fastas, "tmp")
        if (args_ribo.program == "both") or (
                args_ribo.program == "riboswitch"):
            (self.ribos_stat_folder, self.ribos_gff_outfolder,
             self.ribos_table_folder, self.ribos_scan_folder,
             self.ribos_tmp_files, self.ribos_rfam,
             self.ribos_suffixs) = self._create_out_folders(
                args_ribo.ribos_out_folder, "riboswitch",
                args_ribo.database)
        if (args_ribo.program == "both") or (
                args_ribo.program == "thermometer"):
            (self.thermo_stat_folder, self.thermo_gff_outfolder,
             self.thermo_table_folder, self.thermo_scan_folder,
             self.thermo_tmp_files, self.thermo_rfam,
             self.thermo_suffixs) = self._create_out_folders(
                args_ribo.thermo_out_folder, "RNA_thermometer",
                args_ribo.database)

    def _create_out_folders(self, out_folder, feature, database):
        stat_folder = os.path.join(out_folder, "statistics")
        gff_outfolder = os.path.join(out_folder, "gffs")
        table_folder = os.path.join(out_folder, "tables")
        scan_folder = os.path.join(out_folder, "scan_Rfam")
        tmp_files = {"fasta": os.path.join(
                              out_folder, "tmp_fasta"),
                     "scan": os.path.join(
                              out_folder, "tmp_scan"),
                     "table": os.path.join(
                              out_folder, "tmp_table")}
        rfam = os.path.join(database, "Rfam_" + feature + ".cm")
        suffixs = {"csv": feature + ".csv",
                   "txt": feature + "_prescan.txt",
                   "re_txt": feature + "_scan.txt",
                   "re_csv": feature + "_scan.csv"}
        return (stat_folder, gff_outfolder, table_folder, scan_folder,
                tmp_files, rfam, suffixs)

    def _run_infernal(self, args_ribo, seq, type_, prefix, tmp_files,
                      suffixs, rfam):
        scan_file = os.path.join(tmp_files["scan"],
                                 "_".join([prefix, suffixs[type_]]))
        scan = open(scan_file, "w")
        call([os.path.join(args_ribo.infernal_path, "cmscan"), "--incE",
              str(args_ribo.e_value), "--acc", rfam, seq], stdout=scan)
        scan.close()
        return scan_file

    def _scan_extract_rfam(self, prefixs, args_ribo, tmp_files, suffixs,
                           feature, rfam):
        '''extract the seq of candidates and scanning the candidates'''
        for gff in os.listdir(self.gff_path):
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                first_seq = os.path.join(tmp_files["fasta"],
                                         prefix + ".fa")
                prefixs.append(prefix)
                print("extracting seq of candidates of {0}".format(
                      prefix))
                extract_potential_rbs(
                      os.path.join(self.fasta_path, prefix + ".fa"),
                      os.path.join(self.gff_path, gff),
                      os.path.join(self.tss_path, prefix + "_TSS.gff"),
                      os.path.join(self.tran_path, prefix + "_transcript.gff"),
                      first_seq, args_ribo, feature)
                print("pre-scanning of {0}".format(prefix))
                first_scan_file = self._run_infernal(
                        args_ribo, first_seq, "txt", prefix, tmp_files,
                        suffixs, rfam)
                sec_seq = os.path.join(tmp_files["fasta"],
                                       "_".join([prefix, "regenerate.fa"]))
                first_table = os.path.join(
                        tmp_files["table"],
                        "_".join([prefix, suffixs["csv"]]))
                regenerate_seq(first_scan_file, first_seq,
                               first_table, sec_seq)
                print("scanning of {0}".format(prefix))
                sec_scan_file = self._run_infernal(
                        args_ribo, sec_seq, "re_txt", prefix, tmp_files,
                        suffixs, rfam)
                sec_table = os.path.join(
                        tmp_files["table"],
                        "_".join([prefix, suffixs["re_csv"]]))
                reextract_rbs(sec_scan_file, first_table, sec_table)
                shutil.move(sec_table, first_table)
                modify_table(first_table, args_ribo.output_all)
        return prefixs

    def _merge_results(self, args_ribo, scan_folder, suffixs, tmp_files,
                       table_folder, stat_folder, feature_id, gff_outfolder,
                       feature):
        '''merge the results from the results of two searching'''
        for gff in os.listdir(args_ribo.gffs):
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                print("Merge results of {0}".format(prefix))
                pre_strain = ""
                self.helper.check_make_folder(os.path.join(
                                              scan_folder, prefix))
                fh = open(os.path.join(args_ribo.gffs, gff))
                for entry in self.gff_parser.entries(fh):
                    if entry.seq_id != pre_strain:
                        if len(pre_strain) == 0:
                            shutil.copyfile(os.path.join(
                                tmp_files["table"],
                                "_".join([entry.seq_id, suffixs["csv"]])),
                                os.path.join(
                                    table_folder,
                                    "_".join([prefix, suffixs["csv"]])))
                        else:
                            self.helper.merge_file(os.path.join(
                                tmp_files["table"],
                                "_".join([entry.seq_id, suffixs["csv"]])),
                                os.path.join(
                                    table_folder,
                                    "_".join([prefix, suffixs["csv"]])))
                        shutil.copy(os.path.join(
                            tmp_files["scan"],
                            "_".join([entry.seq_id, suffixs["txt"]])),
                            os.path.join(scan_folder, prefix))
                        shutil.copy(os.path.join(
                            tmp_files["scan"],
                            "_".join([entry.seq_id, suffixs["re_txt"]])),
                            os.path.join(scan_folder, prefix))
                        pre_strain = entry.seq_id
                out_stat = os.path.join(
                        stat_folder,
                        "_".join(["stat", prefix, feature + ".txt"]))
                print("compute statistics of {0}".format(prefix))
                stat_and_covert2gff(os.path.join(
                    table_folder, "_".join([prefix, suffixs["csv"]])),
                    feature_id, os.path.join(gff_outfolder,
                        "_".join([prefix, feature + ".gff"])),
                    args_ribo.fuzzy, out_stat, feature)
                fh.close()

    def _remove_tmp(self, args_ribo):
        self.helper.remove_tmp(args_ribo.gffs)
        self.helper.remove_tmp(args_ribo.fastas)
        self.helper.remove_tmp(args_ribo.trans)
        self.helper.remove_tmp(args_ribo.tsss)

    def _remove_overlap(self, gff_path, tmp_files, suffixs):
        for gff in os.listdir(gff_path):
            if gff.endswith(".gff"):
                rbs_overlap(
                    os.path.join(os.path.join(
                        tmp_files["table"],
                        "_".join([gff.replace(".gff", ""),
                                  suffixs["csv"]]))),
                    os.path.join(gff_path, gff))

    def _core_prediction(self, args_ribo, feature_id, rfam, tmp_files,
                         table_folder, feature, scan_folder, suffixs,
                         stat_folder, gff_outfolder, out_folder):
        '''main part of detection'''
        rbs_from_rfam(feature_id, args_ribo.rfam, rfam)
        print("compressing Rfam of " + feature)
        call([os.path.join(args_ribo.infernal_path, "cmpress"),
              "-F", rfam])
        prefixs = []
        self.helper.check_make_folder(tmp_files["fasta"])
        self.helper.check_make_folder(tmp_files["scan"])
        self.helper.check_make_folder(tmp_files["table"])
        prefixs = self._scan_extract_rfam(
                prefixs, args_ribo, tmp_files, suffixs, feature, rfam)
        self._remove_overlap(self.gff_path, tmp_files, suffixs)
        self._merge_results(args_ribo, scan_folder, suffixs, tmp_files,
                            table_folder, stat_folder, feature_id,
                            gff_outfolder, feature)
        mapping_ribos(table_folder, feature_id, feature)
        self.helper.remove_all_content(out_folder, "tmp", "dir")

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
        if (args_ribo.program == "both") or (
                args_ribo.program == "riboswitch"):
            print("Detecting riboswtich now...")
            self._core_prediction(
                    args_ribo, args_ribo.ribos_id, self.ribos_rfam,
                    self.ribos_tmp_files, self.ribos_table_folder,
                    "riboswitch", self.ribos_scan_folder, self.ribos_suffixs,
                    self.ribos_stat_folder, self.ribos_gff_outfolder,
                    args_ribo.ribos_out_folder)
        if (args_ribo.program == "both") or (
                args_ribo.program == "thermometer"):
            print("Detecting RNA thermometer now...")
            self._core_prediction(
                    args_ribo, args_ribo.thermo_id, self.thermo_rfam,
                    self.thermo_tmp_files, self.thermo_table_folder,
                    "RNA_thermometer", self.thermo_scan_folder,
                    self.thermo_suffixs, self.thermo_stat_folder,
                    self.thermo_gff_outfolder, args_ribo.thermo_out_folder)
        self._remove_tmp(args_ribo)
