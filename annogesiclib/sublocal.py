import os
import sys
import shutil
from subprocess import call
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.format_fixer import FormatFixer
from annogesiclib.extract_psortb import extract_psortb
from annogesiclib.stat_sublocal import stat_sublocal
from annogesiclib.gff3 import Gff3Parser


class SubLocal(object):
    '''detection of subcellular localization'''

    def __init__(self, args_sub):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.fixer = FormatFixer()
        self.gff_path = os.path.join(args_sub.gffs, "tmp")
        self.fasta_path = os.path.join(args_sub.fastas, "tmp")
        if args_sub.trans is not None:
            self.tran_path = os.path.join(args_sub.trans, "tmp")
        else:
            self.tran_path = None
        self.out_all = os.path.join(args_sub.out_folder, "all_CDS")
        self.out_express = os.path.join(args_sub.out_folder, "expressed_CDS")
        self.all_tmp_path = os.path.join(self.out_all, "tmp")
        self.express_tmp_path = os.path.join(self.out_express, "tmp")
        self.all_stat_path = os.path.join(self.out_all, "statistics")
        self.express_stat_path = os.path.join(self.out_express, "statistics")
        self.all_tmp_result = os.path.join(self.out_all, "tmp_results")
        self.express_tmp_result = os.path.join(self.out_express, "tmp_results")
        self.all_result = os.path.join(self.out_all, "psortb_results")
        self.express_result = os.path.join(self.out_express, "psortb_results")
        self.endfix_table = "table.csv"
        self.endfix_raw = "raw.txt"
        self._make_folder()

    def _make_folder(self):
        self.helper.check_make_folder(self.out_all)
        self.helper.check_make_folder(self.out_express)
        self.helper.check_make_folder(self.all_stat_path)
        self.helper.check_make_folder(self.express_stat_path)
        self.helper.check_make_folder(self.all_result)
        self.helper.check_make_folder(self.express_result)

    def _compare_cds_tran(self, gff_file, tran_file):
        '''compare CDS and transcript to find the expressed CDS'''
        out = open(os.path.join(self.out_all, "tmp_cds.gff"), "w")
        cdss = []
        fh = open(gff_file)
        th = open(tran_file)
        for entry in Gff3Parser().entries(fh):
            if entry.feature == "CDS":
                cdss.append(entry)
        trans = []
        for entry in Gff3Parser().entries(th):
            trans.append(entry)
        for cds in cdss:
            for ta in trans:
                if (cds.strand == ta.strand) and (
                        cds.seq_id == ta.seq_id):
                    if ((cds.end < ta.end) and (
                             cds.end > ta.start) and (
                             cds.start <= ta.start)) or (
                            (cds.start > ta.start) and (
                             cds.start < ta.end) and (
                             cds.end >= ta.end)) or (
                            (cds.end >= ta.end) and (
                             cds.start <= ta.start)) or (
                            (cds.end <= ta.end) and (
                             cds.start >= ta.start)):
                        out.write(cds.info + "\n")
                        break
        fh.close()
        th.close()
        out.close()

    def _get_protein_seq(self, gff, tmp_path, tran_path):
        prefix = gff.replace(".gff", "")
        fasta = self.helper.get_correct_file(self.fasta_path, ".fa",
                                             prefix, None, None)
        dna_seq_file = os.path.join(tmp_path, "_".join([prefix, "dna.fa"]))
        print("Generate CDS fasta files of {0}".format(prefix))
        if tran_path is not None:
            self._compare_cds_tran(os.path.join(self.gff_path, gff),
                                   os.path.join(tran_path, "_".join([
                                       prefix, "transcript.gff"])))
            self.helper.get_cds_seq(os.path.join(self.out_all, "tmp_cds.gff"),
                                    fasta, dna_seq_file)
            os.remove(os.path.join(self.out_all, "tmp_cds.gff"))
        else:
            self.helper.get_cds_seq(os.path.join(self.gff_path, gff),
                                    fasta, dna_seq_file)
        print("transfer DNA seq to protein seq of {0}".format(prefix))
        self.helper.translation(dna_seq_file, "tmp")
        prot_seq_file = os.path.join(
                tmp_path, "_".join([prefix, "protein.fa"]))
        self.fixer.fix_emboss("tmp", prot_seq_file)
        os.remove("tmp")
        return prefix

    def _psortb(self, psortb_path, strain_type, prot_seq_file,
                out_raw, out_err):
        call([psortb_path, strain_type, prot_seq_file],
             stdout=out_raw, stderr=out_err)

    def _run_psortb(self, args_sub, prefix, out_folder, tmp_path, tmp_result):
        print("Running psortb of {0}".format(prefix))
        out_err = open(os.path.join(out_folder, "tmp_log"), "w")
        out_raw = open(os.path.join(tmp_result,
                       "_".join([prefix, self.endfix_raw])), "w")
        prot_seq_file = os.path.join(tmp_path,
                                     "_".join([prefix, "protein.fa"]))
        if args_sub.gram == "positive":
            self._psortb(args_sub.psortb_path, "-p", prot_seq_file,
                         out_raw, out_err)
        elif args_sub.gram == "negative":
            self._psortb(args_sub.psortb_path, "-n", prot_seq_file,
                         out_raw, out_err)
        else:
            print("Error:It is not a proper bacteria type - {0}!!".format(
                  args_sub.gram))
            sys.exit()
        out_err.close()
        out_raw.close()

    def _extract_result(self, args_sub, tmp_psortb_path, prefix, gff_file):
        '''extract the result of psortb'''
        if args_sub.merge:
            print("Merge to gff...")
            extract_psortb(os.path.join(
                tmp_psortb_path, "_".join([prefix, self.endfix_raw])),
                os.path.join(tmp_psortb_path, "_".join([
                    prefix, self.endfix_table])),
                gff_file, os.path.join(prefix + ".gff"),
                args_sub.fuzzy)
            shutil.move(prefix + ".gff", gff_file)
        else:
            extract_psortb(os.path.join(
                tmp_psortb_path, "_".join([prefix, self.endfix_raw])),
                os.path.join(tmp_psortb_path, "_".join([
                    prefix, self.endfix_table])),
                None, None, args_sub.fuzzy)

    def _merge_and_stat(self, gffs, tmp_psortb_path, stat_path, psortb_result):
        for folder in os.listdir(gffs):
            if folder.endswith(".gff_folder"):
                prefix = folder.replace(".gff_folder", "")
                self.helper.check_make_folder(
                     os.path.join(psortb_result, prefix))
                merge_table = os.path.join(
                        psortb_result, prefix,
                        "_".join([prefix, self.endfix_table]))
                for gff in os.listdir(os.path.join(gffs, folder)):
                    result = self.helper.get_correct_file(
                            tmp_psortb_path, "_" + self.endfix_raw,
                            gff.replace(".gff", ""), None, None)
                    shutil.copy(result, os.path.join(psortb_result, prefix))
                    result = self.helper.get_correct_file(
                            tmp_psortb_path, "_" + self.endfix_table,
                            gff.replace(".gff", ""), None, None)
                    self.helper.merge_file(result, merge_table)
                self.helper.check_make_folder(os.path.join(stat_path, prefix))
                stat_sublocal(merge_table,
                              os.path.join(
                                  stat_path, prefix, prefix),
                              os.path.join(
                                  stat_path, prefix, "_".join([
                                      "stat", prefix, "sublocal.csv"])))

    def _remove_tmps(self, args_sub):
        self.helper.remove_tmp(args_sub.fastas)
        self.helper.remove_tmp(args_sub.gffs)
        self.helper.remove_all_content(args_sub.out_folder, "tmp", "dir")
        self.helper.remove_all_content(self.out_all, "tmp", "dir")
        self.helper.remove_all_content(self.out_express, "tmp", "dir")
        os.remove(os.path.join(self.out_all, "tmp_log"))
        if args_sub.trans is not None:
            os.remove(os.path.join(self.out_express, "tmp_log"))

    def run_sub_local(self, args_sub):
        for gff in os.listdir(args_sub.gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(
                                                 args_sub.gffs, gff))
        self.multiparser.parser_gff(args_sub.gffs, None)
        self.multiparser.parser_fasta(args_sub.fastas)
        if args_sub.trans is not None:
            self.multiparser.parser_gff(args_sub.trans, "transcript")
            self.helper.check_make_folder(self.express_tmp_path)
            self.helper.check_make_folder(self.express_tmp_result)
        self.helper.check_make_folder(self.all_tmp_path)
        self.helper.check_make_folder(self.all_tmp_result)
        for gff in os.listdir(self.gff_path):
            if args_sub.trans is not None:
                print("Running expressed gene now...")
                prefix = self._get_protein_seq(gff, self.express_tmp_path,
                                               self.tran_path)
                self._run_psortb(args_sub, prefix, self.out_express,
                                 self.express_tmp_path,
                                 self.express_tmp_result)
                self._extract_result(args_sub, self.express_tmp_result, prefix,
                                     os.path.join(self.gff_path, gff))
            print("Running all gene now...")
            prefix = self._get_protein_seq(gff, self.all_tmp_path, None)
            self._run_psortb(args_sub, prefix, self.out_all,
                             self.all_tmp_path, self.all_tmp_result)
            self._extract_result(args_sub, self.all_tmp_result, prefix,
                                 os.path.join(self.gff_path, gff))
        self._merge_and_stat(args_sub.gffs, self.all_tmp_result,
                             self.all_stat_path, self.all_result)
        if args_sub.trans is not None:
            self._merge_and_stat(args_sub.gffs, self.express_tmp_result,
                                 self.express_stat_path, self.express_result)
        self._remove_tmps(args_sub)
