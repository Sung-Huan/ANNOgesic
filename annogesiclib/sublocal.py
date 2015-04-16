#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
import shutil
from subprocess import call
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.format_fixer import Format_Fixer
from annogesiclib.extract_psortb import extract_psortb
from annogesiclib.stat_sublocal import stat_sublocal


class Sub_Local(object):

    def __init__(self):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.fixer = Format_Fixer()

    def _get_protein_seq(fasta_path, gff, tmp_path, gff_path, emboss_path):
        prefix = gff.replace(".gff", "")
        fasta = self.helpe.get_correct_file(fasta_path, ".fa", prefix, None)
        dna_seq_file = os.path.join(tmp_path, "_".join([prefix, "dna.fa"]))
        print("Generate CDS fasta files of {0}".format(prefix))
        self.helper.get_cds_seq(os.path.join(gff_path, gff), 
                                fasta, dna_seq_file)
        print("transfer DNA seq to protein seq of {0}".format(prefix))
        call([emboss_path, "-sequence", dna_seq_file, "-outseq", "tmp", "-trim"])
        prot_seq_file = os.path.join(tmp_path, "_".join([prefix, "protein.fa"]))
        self.fixer.fix_emboss("tmp", prot_seq_file)
        os.remove("tmp")
        return prefix

    def _run_psortb(self, prefix, out_folder, gram, psortb_path, tmp_path):
        print("Running psortb of {0}".format(prefix))
        out_err = open(os.path.join(out_folder, "tmp_log"), "w")
        tmp_psortb_path = os.path.join(out_folder, "tmp_results")
        out_raw = open(os.path.join(tmp_psortb_path, 
                       "_".join([prefix, "raw.txt"])), "w")
        prot_seq_file = os.path.join(tmp_path, "_".join([prefix, "protein.fa"]))
        if gram == "positive":
            call([psortb_path, "-p", protein_seq_file], 
                  stdout=out_raw, stderr=out_err)
        elif gram == "negative":
            call([psortb_path, "-n", protein_seq_file], 
                  stdout=out_raw, stderr=out_err)
        else:
            print("Error:It is not a proper bacteria type - {0}!!".format(gram))
            sys.exit()

    def _extract_result(self, merge, tmp_psortb_path, prefix, gff_file):
        if merge is True:
            print("Merge to gff...")
            extract_psortb(os.path.join(tmp_psortb_path, "_".join([prefix, "raw.txt"])),
                           os.path.join(tmp_psortb_path, "_".join([prefix, "table.csv"])),
                           gff_file, os.path.join(prefix + ".gff"))
            os.rename(prefix + ".gff", gff_file)
        else:
            extract_psortb(os.path.join(tmp_psortb_path, "_".join([prefix, "raw.txt"])),
                           os.path.join(tmp_psortb_path, "_".join([prefix, "table.csv"])),
                           None, None)

    def _merge_and_stat(gffs, out_folder, tmp_psortb_path, stat_path):
        for folder in os.listdir(gffs):
            if folder.endswith(".gff_folder"):
                prefix = folder.replace(".gff_folder", "")
                self.helper.check_make_folder(
                            os.path.join(out_folder, "psortb_results"), prefix)
                merge_table = os.path.join(out_folder, "psortb_results", prefix, 
                                           "_".join([prefix, "table.csv"]))
                for gff in os.listdir(os.path.join(gffs, folder)):
                    result = self.helper.get_correct_file(tmp_psortb_path, 
                                         "_raw.txt", gff.replace(".gff", ""), None)
                    shutil.copy(result, os.path.join(out_folder, "psortb_results", prefix))
                    result = self.helper.get_correct_file(tmp_psortb_path, 
                                         "_table.csv", gff.replace(".gff", ""), None)
                    self.merge_file(tmp_psortb_path, "_raw.txt", gff.replace(".gff", ""),
                                    os.path.join(out_folder, "psortb_results", prefix),
                                    "_".join([prefix, "table.csv"]))
                self.helper.check_make_folder(stat_path, prefix) ## statistics
                stat_sublocal(merge_table, os.path.join(stat_path, prefix, prefix), 
                              os.path.join(stat_path, prefix, 
                              "_".join(["stat", prefix, "sublocal.csv"])))

    def run_sub_local(self, bin_path, gffs, fastas, gram, merge, out_folder,
                      emboss_path, psortb_path):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(gffs, gff))
        self.multiparser._parser_gff(gffs, None)
        self.multiparser._parser_fasta(fastas)
        gff_path = os.path.join(gffs, "tmp")
        fasta_path = os.path.join(fastas, "tmp")
        tmp_path = os.path.join(out_folder, "tmp")
        stat_path = os.path.join(out_folder, "statistics")
        self.helper.check_make_folder(out_folder, "tmp")
        self.helper.check_make_folder(out_folder, "tmp_results")
        for gff in os.listdir(gff_path):
            prefix = self._get_protein_seq(fasta_path, gff, tmp_path, gff_path, emboss_path)
            self._run_psortb(prefix, out_folder, gram, psortb_path, tmp_path)
            self._extract_result(merge, tmp_psortb_path, prefix, os.path.join(gff_path, gff))
        self._merge_and_stat(gffs, out_folder, tmp_psortb_path, stat_path)
        self.helper.remove_tmp(fastas)
        self.helper.remove_tmp(gffs)
        self.helper.remove_all_content(out_folder, "tmp", "dir")
