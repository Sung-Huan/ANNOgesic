#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
import shutil
from subprocess import call, Popen
from annogesiclib.helper import Helper
from annogesiclib.multiparser import Multiparser
from annogesiclib.converter import Converter
from annogesiclib.get_inter_seq import intergenic_seq
from annogesiclib.get_polyT import poly_t
from annogesiclib.detect_coverage_term import detect_coverage
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.stat_term import stat_term


class Terminator(object):

    def __init__(self, gffs, fastas, trans, out_folder):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.converter = Converter()
        self.gff_parser = Gff3Parser()
        self.gff_path = os.path.join(gffs, "tmp")
        self.fasta_path = os.path.join(fastas, "tmp")
        self.tran_path = os.path.join(trans, "tmp")
        self.term_outfolder = os.path.join(out_folder, "gffs")
        self.csv_outfolder = os.path.join(out_folder, "tables")
        self.term_all = os.path.join(self.term_outfolder, "all_candidates")
        self.csv_all = os.path.join(self.csv_outfolder, "all_candidates")
        self.term_express = os.path.join(self.term_outfolder, "express")
        self.csv_express = os.path.join(self.csv_outfolder, "express")
        self.term_detect = os.path.join(self.term_outfolder, "detect")
        self.csv_detect = os.path.join(self.csv_outfolder, "detect")
        self.combine_path = os.path.join(self.gff_path, "combine")
        self.tmp_transterm = os.path.join(os.getcwd(), "tmp_transterm")
        self.tmp_hp = "transtermhp"
        self.tmp_hp_gff = "transtermhp.gff"
        self.transterm_path = "tmp_transterm/tmp"
        self.tmp_term_table = os.path.join(os.getcwd(), "tmp_term_table")
        self.tmp_merge = os.path.join(os.getcwd(), "tmp_merge_gff")
        self.tmp_gff = "tmp.gff"
        self.tmp_folder = os.path.join(os.getcwd(), "tmp")
        self.endfix_gff = "term.gff"
        self.endfix_csv = "term.csv"
        self.endfix_allgff = "term_all.gff"
        self._make_gff_folder()


    def _combine_annotation(self, combine_file, files):
        with open(combine_file, 'w') as result:
            for file_ in files:
                check_start = False
                for line in open( file_, 'r' ):
                    if check_start:
                        result.write(line)
                    if "Location" in line:
                        check_start = True

    def _make_gff_folder(self):
        self.helper.check_make_folder(self.term_all)
        self.helper.check_make_folder(self.csv_all)
        self.helper.check_make_folder(self.term_detect)
        self.helper.check_make_folder(self.csv_detect)
        self.helper.check_make_folder(self.term_express)
        self.helper.check_make_folder(self.csv_express)

    def _convert_gff2rntptt(self, gff_path, prefixs, fasta_path, sRNAs, file_types):
        for gff in os.listdir(gff_path):
            if gff.endswith(".gff"):
                filename = gff.split("/")
                prefix = filename[-1][:-4]
                prefixs.append(prefix)
                gff_file = os.path.join(gff_path, gff)
                rnt_file = os.path.join(gff_path, gff.replace(".gff", ".rnt"))
                ptt_file = os.path.join(gff_path, gff.replace(".gff", ".ptt"))
                fasta = self.helper.get_correct_file(
                             fasta_path, ".fa", prefix, None)
                if fasta is None:
                    print("Error: no proper file - {0}.fa".format(prefix))
                    sys.exit()
                if sRNAs is not False:
                    self.multiparser._parser_gff(sRNAs, "sRNA")
                    srna_path = os.path.join(sRNAs, "tmp")
                    srna = self.helper.get_correct_file(
                                srna_path, "_sRNA.gff", prefix, None)
                    if (srna is not None) and (fasta is not None):
                        self.converter.convert_gff2rntptt(gff_file, fasta, 
                             ptt_file, rnt_file, srna, srna.replace(".gff", ".rnt"))
                        file_types[prefix] = "srna"
                    if (srna is None) and (fasta is not None):
                        self.converter.convert_gff2rntptt(gff_file, fasta, 
                             ptt_file, rnt_file, None, None)
                        file_types[prefix] = "normal"
                else:
                    self.converter.convert_gff2rntptt(gff_file, fasta, 
                         ptt_file, rnt_file, None, None)
                    file_types[prefix] = "normal"
        return srna_path

    def _combine_ptt_rnt(self, gff_path, file_types, srna_path):
        self.helper.check_make_folder(self.combine_path)
        for prefix, file_type in file_types.items():
            combine_file = os.path.join(self.combine_path, prefix + '.ptt')
            if file_type == "normal":
                files = [os.path.join(gff_path, prefix + ".ptt"),
                         os.path.join(gff_path, prefix + ".rnt")]
                self._combine_annotation(combine_file, files)
            elif file_type == "srna":
                files = [os.path.join(gff_path, prefix + ".ptt"),
                         os.path.join(gff_path, prefix + ".rnt"),
                         os.path.join(srna_path, "_".join([prefix, "sRNA.rnt"]))]
                self._combine_annotation(combine_file, files)

    def _run_TransTermHP(self, TransTermHP_path, combine_path, fasta_path, hp_folder):
        self.helper.check_make_folder(self.tmp_transterm)
        for file_ in os.listdir(combine_path):
            if ".ptt" in file_:
                prefix = file_.replace(".ptt", "")
                fasta = self.helper.get_correct_file(fasta_path, ".fa", prefix, None)
                if fasta is None:
                    print("Error: no proper file - {0}.fa".format(prefix))
                    sys.exit()
                out_path = os.path.join(hp_folder, prefix)
                self.helper.check_make_folder(out_path)
                out = open(os.path.join(out_path, "_".join([prefix, "terminators.txt"])), "w")
                call([os.path.join(TransTermHP_path, "transterm"), "-p", 
                      os.path.join(TransTermHP_path, "expterm.dat"),
                      fasta, os.path.join(combine_path, file_), "--t2t-perf", 
                      os.path.join(out_path, 
                      "_".join([prefix, "terminators_within_robust_tail-to-tail_regions.t2t"])),
                      "--bag-output", os.path.join(out_path, 
                      "_".join([prefix, "best_terminator_after_gene.bag"]))],
                      stdout=out)
        shutil.rmtree(combine_path)

    def _convert_to_gff(self, prefixs, hp_folder, gffs):
        for prefix in prefixs:
            for folder in os.listdir(hp_folder):
                if prefix == folder:
                    out_path = os.path.join(hp_folder, folder)
                    for file_ in os.listdir(out_path):
                        if file_.endswith(".bag"):
                            out_file = os.path.join(self.tmp_transterm, 
                                                    "_".join([prefix, self.tmp_hp_gff]))
                            self.converter.convert_transtermhp2gff(
                                                    os.path.join(out_path, file_), out_file)
        self.multiparser._combine_gff(gffs, self.tmp_transterm, None, self.tmp_hp)

    def _combine_libs_wigs(self, tlibs, flibs, tex_wigs, frag_wigs):
        if (tlibs is False) and (flibs is False):
            print("Error: please input proper libraries!!")
        if (tlibs is not False) and (flibs is not False):
            libs = tlibs + flibs
        elif (tlibs is not False):
            libs = tlibs
        elif (flibs is not False):
            libs = flibs
        if (tex_wigs is not False) and (frag_wigs is not False):
            folder = tex_wigs.split("/")
            folder = "/".join(folder[:-1])
            merge_wigs = os.path.join(folder, "merge_wigs")
            self.helper.check_make_folder(merge_wigs)
            for wig in os.listdir(tex_wigs):
                if os.path.isdir(os.path.join(tex_wigs, wig)):
                    pass
                else:
                    shutil.copy(os.path.join(tex_wigs, wig), merge_wigs)
            for wig in os.listdir(frag_wigs):
                if os.path.isdir(os.path.join(frag_wigs, wig)):
                    pass
                else:
                    shutil.copy(os.path.join(frag_wigs, wig), merge_wigs)
        elif (tex_wigs is not False):
            merge_wigs = tex_wigs
        elif (frag_wigs is not False):
            merge_wigs = frag_wigs
        else:
            print("Error: no proper wig files!!!")
            sys.exit()
        return (merge_wigs, libs)

    def _merge_sRNA(self, sRNAs, prefixs, gff_path, srna_path):
        if sRNAs is not False:
            self.multiparser._parser_gff(sRNAs, "sRNA")
            srna_path = os.path.join(sRNAs, "tmp")
            self.helper.check_make_folder(self.tmp_merge)
            for prefix in prefixs:
                if self.tmp_gff in os.listdir(self.tmp_merge):
                    os.remove(os.path.join(self.tmp_merge, self.tmp_gff))
                self.helper.merge_file(gff_path, prefix + ".gff", 
                                       self.tmp_merge, self.tmp_gff)
                self.helper.merge_file(srna_path, "_".join([prefix, "sRNA.gff"]), 
                                       self.tmp_merge, self.tmp_gff)
                self.helper.sort_gff(os.path.join(self.tmp_merge, self.tmp_gff), 
                                     os.path.join(self.tmp_merge, prefix + ".gff"))
                os.remove(os.path.join(self.tmp_merge, self.tmp_gff))
            merge_path = self.tmp_merge
        else:
            merge_path = gff_path
        return merge_path

    def _move_file(self, term_outfolder, csv_outfolder):
        for gff in os.listdir(term_outfolder):
            if gff.endswith("_term.gff"):
                self.helper.sort_gff(os.path.join(term_outfolder, gff), self.tmp_gff)
                os.rename(self.tmp_gff, os.path.join(term_outfolder, gff))
                prefix = gff.replace("_term.gff", "")
                new_gff = os.path.join(self.term_all, "_".join([prefix, self.endfix_allgff]))
                csv_file = os.path.join(os.path.join(self.csv_all, 
                                        "_".join([prefix, self.endfix_csv])))
                out = open(new_gff, "w")
                out.write("##gff-version 3\n")
                out.close()
                self.helper.merge_file(term_outfolder, gff, self.term_all, 
                                       "_".join([prefix, self.endfix_allgff]))
                os.remove(os.path.join(term_outfolder, gff))
                pre_strain = ""
                if "_".join([prefix, self.endfix_csv]) in os.listdir(self.csv_all):
                    os.remove(csv_file)
                out_csv = open(csv_file, "w")
                out_csv.write("\t".join(["strain", "name", "start", "end", "strand",
                                         "detect", "coverage_detail"]) + "\n")
                out_csv.close()
                for entry in self.gff_parser.entries(open(new_gff)):
                    if entry.seq_id != pre_strain:
                        self.helper.merge_file(self.tmp_term_table, "_".join([entry.seq_id, "term_raw.csv"]),
                                               self.csv_all, "_".join([prefix, self.endfix_csv]))
                    pre_strain = entry.seq_id

    def _compute_intersection_forward_reverse(self, RNAfold_path, prefixs, tran_path, 
                        merge_path, fasta_path, cutoff_coverage, fuzzy, wig_path, merge_wigs, 
                        libs, tex_notex, replicates, decrease, term_outfolder, csv_outfolder,
                        table_best, gffs, out_folder):
        for prefix in prefixs:
            tmp_seq = os.path.join(out_folder, "_".join(["inter_seq", prefix]))
            tmp_sec = os.path.join(out_folder, "_".join(["inter_sec", prefix]))
            ### get intergenic seq, sec
            print("Extracting seq of {0}".format(prefix))
            out_seq = open(tmp_seq, "w")
            intergenic_seq(os.path.join(fasta_path, prefix + ".fa"), 
                           os.path.join(tran_path, "_".join([prefix, "transcript.gff"])), 
                           os.path.join(merge_path, prefix + ".gff"), tmp_seq)
            print("Computing secondray structure of {0}".format(prefix))
            self.helper.check_make_folder(self.tmp_folder)
            pre_cwd = os.getcwd()
            os.chdir(self.tmp_folder)
            os.system(" ".join([RNAfold_path, "<", os.path.join("..", tmp_seq), 
                                ">", os.path.join("..", tmp_sec)]))
            os.chdir(pre_cwd)
            shutil.rmtree(self.tmp_folder)
            ### detect poly U/T tail of terminators and coverage decreasing.
            tmp_cand = os.path.join(out_folder, "_".join(["term_candidates", prefix]))
            poly_t(tmp_seq, tmp_sec, os.path.join(merge_path, prefix + ".gff"), 
                   fuzzy, tmp_cand)
            print("detection of terminator")
            detect_coverage(tmp_cand, 
                            os.path.join(merge_path, prefix + ".gff"), 
                            os.path.join(tran_path, "_".join([prefix, "transcript.gff"])),
                            os.path.join(fasta_path, prefix + ".fa"),
                            os.path.join(wig_path, "_".join([prefix, "forward.wig"])),
                            os.path.join(wig_path, "_".join([prefix, "forward.wig"])), 
                            fuzzy, cutoff_coverage, 
                            os.path.join(self.transterm_path, "_".join([prefix, self.tmp_hp_gff])), 
                            merge_wigs, libs, tex_notex, replicates, 
                            os.path.join(term_outfolder, "_".join([prefix, self.endfix_gff])), 
                            os.path.join(self.tmp_term_table, "_".join([prefix, "term_raw.csv"])), 
                            table_best, decrease)
        self.multiparser._combine_gff(gffs, term_outfolder, None, "term")
        self._move_file(term_outfolder, csv_outfolder)

    def _remove_tmp_file(self, gffs, fastas, sRNAs, tex_wigs, frag_wigs, term_outfolder, out_folder):
        self.helper.remove_tmp(gffs)
        self.helper.remove_tmp(fastas)
        if sRNAs is not False:
            self.helper.remove_tmp(sRNAs)
        self.helper.remove_tmp(tex_wigs)
        self.helper.remove_tmp(frag_wigs)
        self.helper.remove_tmp(term_outfolder)
        shutil.rmtree(self.tmp_transterm)
        shutil.rmtree(self.tmp_term_table)
        shutil.rmtree(self.tmp_merge)
        self.helper.remove_all_content(out_folder, "inter_seq_", "file")
        self.helper.remove_all_content(out_folder, "inter_sec_", "file")
        self.helper.remove_all_content(out_folder, "term_candidates_", "file")


    def _compute_stat(self, term_outfolder, csv_outfolder, stat, out_folder):
        new_prefixs = []
        for gff in os.listdir(self.term_all):
            if gff.endswith("_term_all.gff"):
                out_tmp = open(self.tmp_gff, "w")
                out_tmp.write("##gff-version 3\n")
                new_prefix = gff.replace("_term_all.gff", "")
                new_prefixs.append(gff.replace("_term_all.gff", ""))
                num = 0
                for entry in self.gff_parser.entries(
                             open(os.path.join(self.term_all, gff))):
                    name ='%0*d' % (5, num)
                    entry.attributes["ID"] = "term" + str(num)
                    entry.attributes["Name"] = "_".join(["Terminator_" + name])
                    entry.attribute_string = ";".join(
                                    ["=".join(items) for items in entry.attributes.items()])
                    out_tmp.write("\t".join([entry.info_without_attributes, entry.attribute_string]) + "\n")
                    num += 1
                os.rename(self.tmp_gff, os.path.join(self.term_all, 
                                                  "_".join([new_prefix, self.endfix_gff])))
        if stat is not False:
            stat_path = os.path.join(out_folder, "statistics")
            for prefix in new_prefixs:
                stat_term(os.path.join(self.term_all, "_".join([prefix, self.endfix_gff])), 
                          os.path.join(self.csv_all, "_".join([prefix, self.endfix_csv])), 
                          os.path.join(stat_path, "_".join(["stat", prefix + ".csv"])), 
                          os.path.join(self.term_detect, "_".join([prefix, "term"])), 
                          os.path.join(self.term_express, "_".join([prefix, "term"])))
                os.rename(os.path.join(self.term_detect, "_".join([prefix, self.endfix_csv])), 
                          os.path.join(self.csv_detect, "_".join([prefix, self.endfix_csv])))
                os.rename(os.path.join(self.term_express, "_".join([prefix, self.endfix_csv])),
                          os.path.join(self.csv_express, "_".join([prefix, self.endfix_csv])))
                os.remove(os.path.join(self.term_all, "_".join([prefix, self.endfix_allgff])))

    def _check_gff_file(self, folder):
        for file_ in os.listdir(folder):
            if file_.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(folder, file_))

    def run_terminator(self, TransTermHP_path, RNAfold_path, out_folder, fastas, gffs, 
                       trans, sRNAs, stat, tex_wigs, frag_wigs, decrease, cutoff_coverage,
                       fuzzy, hp_folder, tlibs, flibs, tex_notex, replicates, table_best):
        ### First, running TransTermHP. Before, running TransTermHP, 
        ### we need to convert annotation files to .ptt and .rnt files.
        file_types = {}
        prefixs = []
        self._check_gff_file(gffs)
        self._check_gff_file(trans)
        self.multiparser._parser_gff(gffs, None)
        self.multiparser._parser_fasta(fastas)
        if (gffs is None) or (fastas is None):
            print("Error: please assign gff annotation folder and fasta folder!!!")
            sys.exit()
        srna_path = self._convert_gff2rntptt(self.gff_path, prefixs, self.fasta_path, sRNAs, file_types)
        self._combine_ptt_rnt(self.gff_path, file_types, srna_path) ### combine .ptt and .rnt to one file.
        self._run_TransTermHP(TransTermHP_path, self.combine_path, self.fasta_path, hp_folder) ### Start to run TransTermHP.
        self._convert_to_gff(prefixs, hp_folder, gffs)
        lib_datas = self._combine_libs_wigs(tlibs, flibs, tex_wigs, frag_wigs)### combine tex_notex and frag libraries.
        merge_wigs = lib_datas[0]
        libs = lib_datas[1]
        self.multiparser._parser_gff(self.gff_path, None)
        wig_path = os.path.join(merge_wigs, "tmp")
        self.multiparser._parser_wig(merge_wigs)
        self.multiparser._combine_wig(self.gff_path, wig_path, None)
        self.helper.remove_tmp(self.gff_path)
        self.multiparser._parser_gff(trans, "transcript")
        self.helper.check_make_folder(self.tmp_term_table)
        self.multiparser._parser_gff(self.tmp_transterm, self.tmp_hp)
        merge_path = self._merge_sRNA(sRNAs, prefixs, self.gff_path, srna_path)
        ### Second, running the cross regions of forward and reverese strand. 
        self._compute_intersection_forward_reverse(RNAfold_path, prefixs, self.tran_path, merge_path, 
                    self.fasta_path, cutoff_coverage, fuzzy, wig_path, merge_wigs, libs, tex_notex, 
                    replicates, decrease, self.term_outfolder, self.csv_outfolder, table_best, gffs, out_folder)
        self._compute_stat(self.term_outfolder, self.csv_outfolder, stat, out_folder)
        self._remove_tmp_file(gffs, fastas, sRNAs, tex_wigs, frag_wigs, self.term_outfolder, out_folder)
