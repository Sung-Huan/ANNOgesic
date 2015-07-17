#!/usr/bin/python

import os
import sys
import shutil
from subprocess import call
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.format_fixer import FormatFixer
from annogesiclib.extract_psortb import extract_psortb
from annogesiclib.stat_sublocal import stat_sublocal


class SubLocal(object):

    def __init__(self, gffs, fastas, out_folder):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.fixer = FormatFixer()
        self.gff_path = os.path.join(gffs, "tmp")
        self.fasta_path = os.path.join(fastas, "tmp")
        self.tmp_path = os.path.join(out_folder, "tmp")
        self.stat_path = os.path.join(out_folder, "statistics")
        self.tmp_result = os.path.join(out_folder, "tmp_results")
        self.psortb_result = os.path.join(out_folder, "psortb_results")
        self.endfix_table = "table.csv"
        self.endfix_raw = "raw.txt"

    def _get_protein_seq(self, fasta_path, gff, tmp_path, gff_path):
        prefix = gff.replace(".gff", "")
        fasta = self.helper.get_correct_file(fasta_path, ".fa", prefix, None)
        dna_seq_file = os.path.join(tmp_path, "_".join([prefix, "dna.fa"]))
        print("Generate CDS fasta files of {0}".format(prefix))
        self.helper.get_cds_seq(os.path.join(gff_path, gff),
                                fasta, dna_seq_file)
        print("transfer DNA seq to protein seq of {0}".format(prefix))
        self.helper.translation(dna_seq_file, "tmp")
        prot_seq_file = os.path.join(tmp_path, "_".join([prefix, "protein.fa"]))
        self.fixer.fix_emboss("tmp", prot_seq_file)
        os.remove("tmp")
        return prefix

    def _run_psortb(self, prefix, out_folder, gram, psortb_path, tmp_path):
        print("Running psortb of {0}".format(prefix))
        out_err = open(os.path.join(out_folder, "tmp_log"), "w")
        out_raw = open(os.path.join(self.tmp_result,
                       "_".join([prefix, self.endfix_raw])), "w")
        prot_seq_file = os.path.join(tmp_path, "_".join([prefix, "protein.fa"]))
        if gram == "positive":
            call([psortb_path, "-p", prot_seq_file],
                  stdout=out_raw, stderr=out_err)
        elif gram == "negative":
            call([psortb_path, "-n", prot_seq_file],
                  stdout=out_raw, stderr=out_err)
        else:
            print("Error:It is not a proper bacteria type - {0}!!".format(gram))
            sys.exit()

    def _extract_result(self, merge, tmp_psortb_path, prefix, gff_file, fuzzy):
        if merge:
            print("Merge to gff...")
            extract_psortb(os.path.join(tmp_psortb_path,
                           "_".join([prefix, self.endfix_raw])),
                           os.path.join(tmp_psortb_path,
                           "_".join([prefix, self.endfix_table])),
                           gff_file, os.path.join(prefix + ".gff"), fuzzy)
            os.rename(prefix + ".gff", gff_file)
        else:
            extract_psortb(os.path.join(tmp_psortb_path,
                           "_".join([prefix, self.endfix_raw])),
                           os.path.join(tmp_psortb_path,
                           "_".join([prefix, self.endfix_table])),
                           None, None, fuzzy)

    def _merge_and_stat(self, gffs, out_folder, tmp_psortb_path, stat_path):
        for folder in os.listdir(gffs):
            if folder.endswith(".gff_folder"):
                prefix = folder.replace(".gff_folder", "")
                self.helper.check_make_folder(
                     os.path.join(self.psortb_result, prefix))
                merge_table = os.path.join(self.psortb_result, prefix,
                              "_".join([prefix, self.endfix_table]))
                for gff in os.listdir(os.path.join(gffs, folder)):
                    result = self.helper.get_correct_file(tmp_psortb_path,
                                         "_" + self.endfix_raw,
                                         gff.replace(".gff", ""), None)
                    shutil.copy(result, os.path.join(self.psortb_result, prefix))
                    result = self.helper.get_correct_file(tmp_psortb_path,
                                         "_" + self.endfix_table,
                                         gff.replace(".gff", ""), None)
                    self.helper.merge_file(result, merge_table)
                self.helper.check_make_folder(os.path.join(stat_path, prefix))
                stat_sublocal(merge_table, os.path.join(stat_path,
                              prefix, prefix), os.path.join(stat_path, prefix,
                              "_".join(["stat", prefix, "sublocal.csv"])))

    def run_sub_local(self, psortb_path, gffs, fastas, gram, fuzzy,
                      merge, out_folder):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(gffs, gff))
        self.multiparser.parser_gff(gffs, None)
        self.multiparser.parser_fasta(fastas)
        self.helper.check_make_folder(self.tmp_path)
        self.helper.check_make_folder(self.tmp_result)
        for gff in os.listdir(self.gff_path):
            prefix = self._get_protein_seq(self.fasta_path, gff, self.tmp_path,
                                           self.gff_path)
            self._run_psortb(prefix, out_folder, gram,
                             psortb_path, self.tmp_path)
            self._extract_result(merge, self.tmp_result, prefix,
                                 os.path.join(self.gff_path, gff), fuzzy)
        self._merge_and_stat(gffs, out_folder, self.tmp_result, self.stat_path)
        self.helper.remove_tmp(fastas)
        self.helper.remove_tmp(gffs)
        self.helper.remove_all_content(out_folder, "tmp", "dir")
