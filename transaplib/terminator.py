#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
from subprocess import call, Popen
from transaplib.helper import Helper
from transaplib.multiparser import Multiparser
from transaplib.converter import Converter
from transaplib.get_inter_seq import Intergenic_seq
from transaplib.get_polyT import Poly_T
from transaplib.detect_coverage_term import Detect_coverage
from transaplib.gff3 import Gff3Parser
from transaplib.stat_term import Stat


class Terminator(object):

    def __init__(self):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.converter = Converter()
        self.gff_parser = Gff3Parser()

    def _combine_annotation(self, combine_file, files):
        with open(combine_file, 'w') as result:
            for file_ in files:
                check_start = False
                for line in open( file_, 'r' ):
                    if check_start:
                        result.write(line)
                    if "Location" in line:
                        check_start = True

    def _make_gff_folder(self, term_outfolder, csv_outfolder):
        self.helper.check_make_folder(term_outfolder, "all_candidates")
        self.helper.check_make_folder(csv_outfolder, "all_candidates")
        self.helper.check_make_folder(term_outfolder, "detect")
        self.helper.check_make_folder(csv_outfolder, "detect")
        self.helper.check_make_folder(term_outfolder, "express")
        self.helper.check_make_folder(csv_outfolder, "express")

    def _convert_gff2rntptt(self, gff_path, prefixs, fasta_path, sRNAs, file_types):
        for gff in os.listdir(gff_path):
            if gff.endswith(".gff"):
                filename = gff.split("/")
                prefix = filename[-1][:-4]
                prefixs.append(prefix)
                rnt = gff[:-3] + "rnt"
                ptt = gff[:-3] + "ptt"
                fasta = self.helper.get_correct_file(
                             fasta_path, ".fa", prefix, None)
                if fasta is None:
                    print("Error: no proper file - " + prefix + ".fa")
                    sys.exit()
                if sRNAs is not False:
                    self.multiparser._parser_gff(sRNAs, "sRNA")
                    srna_path = sRNAs + "/tmp/"
                    srna = self.helper.get_correct_file(
                                srna_path, "_sRNA.gff", prefix, None)
                    if (srna is not None) and (fasta is not None):
                        self.converter.Convert_gff2rntptt(gff_path + gff, fasta, 
                             gff_path + ptt, gff_path + rnt, srna, srna[:-3] + "rnt")
                        file_types[prefix] = "srna"
                    if (srna is None) and (fasta is not None):
                        self.converter.Convert_gff2rntptt(gff_path + gff, fasta, 
                             gff_path + ptt, gff_path + rnt, None, None)
                        file_types[prefix] = "normal"
                else:
                    self.converter.Convert_gff2rntptt(gff_path + gff, fasta, 
                         gff_path + ptt, gff_path + rnt, None, None)
                    file_types[prefix] = "normal"
        return srna_path

    def _combine_ptt_rnt(self, gff_path, file_types, srna_path):
        combine_path = gff_path + "/combine/"
        self.helper.check_make_folder(gff_path, "combine")
        for prefix, file_type in file_types.items():
            combine_file = combine_path + prefix + '.ptt'
            if file_type == "normal":
                files = [gff_path + prefix + ".ptt",
                         gff_path + prefix + ".rnt"]
                self._combine_annotation(combine_file, files)
            elif file_type == "srna":
                files = [gff_path + prefix + ".ptt",
                         gff_path + prefix + ".rnt",
                         srna_path + prefix + "_sRNA.rnt"]
                self._combine_annotation(combine_file, files)
        return combine_path

    def _run_TransTermHP(self, TransTermHP_path, combine_path, fasta_path, hp_folder):
        self.helper.check_make_folder(os.getcwd() + "/", "tmp_transterm")
        for file_ in os.listdir(combine_path):
            if ".ptt" in file_:
                prefix = file_.replace(".ptt", "")
                fasta = self.helper.get_correct_file(fasta_path, ".fa", prefix, None)
                if fasta is None:
                    print("Error: no proper file - " + prefix + ".fa")
                    sys.exit()
                out_path = hp_folder + "/" + prefix + "/"
                self.helper.check_make_folder(hp_folder + "/", prefix)
                out = open(out_path + prefix + "_terminators.txt", "w")
                call([TransTermHP_path + "/transterm", "-p", 
                      TransTermHP_path + "/expterm.dat",
                      fasta, combine_path + file_, "--t2t-perf", 
                      out_path + prefix + "_terminators_within_robust_tail-to-tail_regions.t2t",
                      "--bag-output", out_path + prefix + "_best_terminator_after_gene.bag"],
                      stdout=out)
        call(["rm", "-rf", combine_path])

    def _convert_to_gff(self, prefixs, hp_folder, gffs):
        for prefix in prefixs:
            for folder in os.listdir(hp_folder):
                if prefix == folder:
                    out_path = hp_folder + "/" + folder + "/"
                    for file_ in os.listdir(out_path):
                        if file_.endswith(".bag"):
                            out_file = "tmp_transterm/" + prefix + "_transtermhp.gff"
                            self.converter.Convert_TransTermHP2gff(out_path + file_, out_file)
        self.multiparser._combine_gff(gffs, "tmp_transterm", None, "transtermhp")

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
            merge_wigs = os.getcwd() + "/merge_wigs"
            self.helper.check_make_folder(os.getcwd() + "/", "merge_wigs")
            os.system("cp " + tex_wigs + "/* merge_wigs/")
            os.system("cp " + frag_wigs + "/* merge_wigs/")
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
            srna_path = sRNAs + "/tmp/"
            self.helper.check_make_folder(os.getcwd() + "/", "tmp_merge_gff")
            for prefix in prefixs:
                os.system("cat " + gff_path + prefix + ".gff > tmp_merge_gff/tmp.gff")
                os.system("cat " + srna_path + prefix + "_sRNA.gff >> tmp_merge_gff/tmp.gff")
                self.helper.sort_gff("tmp_merge_gff/tmp.gff", "tmp_merge_gff/" + prefix + ".gff")
                os.system("rm tmp_merge_gff/tmp.gff")
            merge_path = "tmp_merge_gff/"
        else:
            merge_path = gff_path
        return merge_path

    def _move_file(self, term_outfolder, csv_outfolder):
        for gff in os.listdir(term_outfolder):
            if gff.endswith("_term.gff"):
                self.helper.sort_gff(term_outfolder + gff, "tmp.gff")
                os.system("mv tmp.gff " + term_outfolder + gff)
                prefix = gff.replace("_term.gff", "")
                os.system("echo \"##gff-version 3\" > " + term_outfolder + "all_candidates/" + prefix + "_term_all.gff")
                os.system("cat " + term_outfolder + gff + " >> " + term_outfolder + "all_candidates/" + prefix + "_term_all.gff")
                os.system("rm " + term_outfolder + gff)
                new_gff = term_outfolder + "all_candidates/" + prefix + "_term_all.gff"
                pre_strain = ""
                first = True
                for entry in self.gff_parser.entries(open(new_gff)):
                    if entry.seq_id != pre_strain:
                        if first:
                            os.system("cat tmp_term_table/" + entry.seq_id + "_term_raw.csv > " +
                                      csv_outfolder + "all_candidates/" + prefix + "_term.csv")
                            first = False
                        else:
                            os.system("cat tmp_term_table/" + entry.seq_id + "_term_raw.csv >> " +
                                      csv_outfolder + "all_candidates/" + prefix + "_term.csv")
                    pre_strain = entry.seq_id

    def _compute_intersection_forward_reverse(self, RNAfold_path, prefixs, tran_path, 
                        merge_path, fasta_path, cutoff_coverage, fuzzy, wig_path, merge_wigs, 
                        libs, tex_notex, replicates, decrease, term_outfolder, csv_outfolder,
                        table_best, gffs):
        for prefix in prefixs:
            ### get intergenic seq, sec
            print("Extracting seq of " + prefix)
            out_seq = open("inter_seq_" + prefix, "w")
            Intergenic_seq(fasta_path + prefix + ".fa", 
                           tran_path + prefix + "_transcript.gff", 
                           merge_path + prefix + ".gff", 
                           "inter_seq_" + prefix)
            print("Computing secondray structure of " + prefix)
#            self.helper.check_make_folder(os.getcwd() + "/", "tmp")
#            pre_cwd = os.getcwd()
#            os.chdir(os.getcwd() + "/tmp")
#            os.system(RNAfold_path + " < ../inter_seq_" + prefix + " > ../inter_sec_" + prefix)
#            os.chdir(pre_cwd)
#            os.system("rm -rf tmp")
#            ### detect poly U/T tail of terminators and coverage decreasing.
#            Poly_T("inter_seq_" + prefix, "inter_sec_" + prefix, 
#                   merge_path + prefix + ".gff", fuzzy, "term_candidates_" + prefix)
            print("detection of terminator")
            Detect_coverage("term_candidates_" + prefix, merge_path + prefix + ".gff", 
                            tran_path + prefix + "_transcript.gff",
                            fasta_path + prefix + ".fa",
                            wig_path + prefix + "_forward.wig",
                            wig_path + prefix + "_reverse.wig", fuzzy, cutoff_coverage, 
                            "tmp_transterm/tmp/" + prefix + "_transtermhp.gff", 
                            merge_wigs, libs, tex_notex, replicates, 
                            term_outfolder + prefix + "_term.gff", 
                            "tmp_term_table/" + prefix + "_term_raw.csv", 
                            table_best, decrease)
        self.multiparser._combine_gff(gffs, term_outfolder, None, "term")
        self._move_file(term_outfolder, csv_outfolder)

    def _remove_tmp_file(self, gffs, fastas, sRNAs, tex_wigs, frag_wigs, term_outfolder):
        self.helper.remove_tmp(gffs)
        self.helper.remove_tmp(fastas)
        self.helper.remove_tmp(sRNAs)
        self.helper.remove_tmp(tex_wigs)
        self.helper.remove_tmp(frag_wigs)
        self.helper.remove_tmp(term_outfolder)
        os.system("rm -rf merge_wigs")
        os.system("rm -rf tmp_transterm")
        os.system("rm -rf tmp_term_table")
        os.system("rm -rf tmp_merge_gff")

    def _compute_stat(self, term_outfolder, csv_outfolder, stat, out_folder):
        new_prefixs = []
        for gff in os.listdir(term_outfolder + "all_candidates"):
            if gff.endswith("_term_all.gff"):
                out_tmp = open("tmp.gff", "w")
                new_prefix = gff.replace("_term_all.gff", "")
                new_prefixs.append(gff.replace("_term_all.gff", ""))
                num = 0
                for entry in self.gff_parser.entries(open(term_outfolder + "all_candidates/" + gff)):
                    name ='%0*d' % (5, num)
                    entry.attributes["ID"] = "term" + str(num)
                    entry.attributes["Name"] = "Terminator_" + name
                    entry.attribute_string = ";".join(
                                    ["=".join(items) for items in entry.attributes.items()])
                    out_tmp.write(entry.info_without_attributes + "\t" + entry.attribute_string + "\n")
                    num += 1
                os.system("mv tmp.gff " + term_outfolder + "all_candidates/" + new_prefix + "_term.gff")
        if stat is not False:
            stat_path = out_folder + "/statistics/"
            for prefix in new_prefixs:
                Stat(term_outfolder + "all_candidates/" + prefix + "_term.gff", 
                     csv_outfolder + "all_candidates/" + prefix + "_term.csv", 
                     stat_path + "stat_" + prefix + ".csv", 
                     term_outfolder + "detect/" + prefix + "_term", 
                     term_outfolder + "express/" + prefix + "_term")
                call(["mv", term_outfolder + "detect/" + prefix + "_term.csv", csv_outfolder + "detect/"])
                call(["mv", term_outfolder + "express/" + prefix + "_term.csv", csv_outfolder + "express/"])

    def _check_gff_file(self, folder):
        for file_ in os.listdir(folder):
            if file_.endswith(".gff"):
                self.helper.check_uni_attributes(folder + "/" + file_)

    def run_Terminator(self, TransTermHP_path, RNAfold_path, out_folder, fastas, gffs, 
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
        gff_path = gffs + "/tmp/"
        fasta_path = fastas + "/tmp/"
        term_outfolder = out_folder + "/gffs/"
        csv_outfolder = out_folder + "/tables/"
        self._make_gff_folder(term_outfolder, csv_outfolder)
        if (gffs is None) or (fastas is None):
            print("Error: please assign gff annotation folder and fasta folder!!!")
            sys.exit()
        srna_path = self._convert_gff2rntptt(gff_path, prefixs, fasta_path, sRNAs, file_types)
        combine_path = self._combine_ptt_rnt(gff_path, file_types, srna_path) ### combine .ptt and .rnt to one file.
        self._run_TransTermHP(TransTermHP_path, combine_path, fasta_path, hp_folder) ### Start to run TransTermHP.
        self._convert_to_gff(prefixs, hp_folder, gffs)
        lib_datas = self._combine_libs_wigs(tlibs, flibs, tex_wigs, frag_wigs)### combine tex_notex and frag libraries.
        merge_wigs = lib_datas[0]
        libs = lib_datas[1]
        self.multiparser._parser_gff(gff_path, None)
        wig_path = merge_wigs + "/tmp/"
        self.multiparser._parser_wig(merge_wigs)
        self.multiparser._combine_wig(gff_path, wig_path, None)
        self.helper.remove_tmp(gff_path)
        self.multiparser._parser_gff(trans, "transcript")
        tran_path = trans + "/tmp/"
        self.helper.check_make_folder(os.getcwd() + "/", "tmp_term_table")
        self.multiparser._parser_gff("tmp_transterm", "transtermhp")
        transterms_path = "tmp_transterm/tmp/"
        merge_path = self._merge_sRNA(sRNAs, prefixs, gff_path, srna_path)
        ### Second, running the cross regions of forward and reverese strand. 
        self._compute_intersection_forward_reverse(RNAfold_path, prefixs, tran_path, merge_path, 
                        fasta_path, cutoff_coverage, fuzzy, wig_path, merge_wigs, libs, tex_notex, 
                        replicates, decrease, term_outfolder, csv_outfolder, table_best, gffs)
        self._compute_stat(term_outfolder, csv_outfolder, stat, out_folder)
        self._remove_tmp_file(gffs, fastas, sRNAs, tex_wigs, frag_wigs, term_outfolder)
