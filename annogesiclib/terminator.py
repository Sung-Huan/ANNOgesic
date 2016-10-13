import os
import sys
import shutil
from subprocess import call
from annogesiclib.helper import Helper
from annogesiclib.multiparser import Multiparser
from annogesiclib.converter import Converter
from annogesiclib.get_inter_seq import intergenic_seq
from annogesiclib.get_polyT import poly_t
from annogesiclib.detect_coverage_term import detect_coverage
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.stat_term import stat_term
from annogesiclib.compare_tran_term import compare_term_tran


class Terminator(object):
    '''detection of terminator'''

    def __init__(self, args_term):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.converter = Converter()
        self.gff_parser = Gff3Parser()
        self.gff_path = os.path.join(args_term.gffs, "tmp")
        self.fasta_path = os.path.join(args_term.fastas, "tmp")
        self.tran_path = os.path.join(args_term.trans, "tmp")
        self.outfolder = {"term": os.path.join(args_term.out_folder, "gffs"),
                          "csv": os.path.join(args_term.out_folder, "tables")}
        self.terms = {"all": os.path.join(self.outfolder["term"],
                                          "all_candidates"),
                      "express": os.path.join(self.outfolder["term"],
                                              "express"),
                      "best": os.path.join(self.outfolder["term"], "best"),
                      "non": os.path.join(self.outfolder["term"],
                                          "non_express")}
        self.csvs = {"all": os.path.join(self.outfolder["csv"],
                                         "all_candidates"),
                     "express": os.path.join(self.outfolder["csv"], "express"),
                     "best": os.path.join(self.outfolder["csv"], "best"),
                     "non": os.path.join(self.outfolder["csv"], "non_express")}
        self.combine_path = os.path.join(self.gff_path, "combine")
        self.tmps = {"transterm": os.path.join(os.getcwd(), "tmp_transterm"),
                     "hp": "transtermhp", "hp_gff": "transtermhp.gff",
                     "hp_path": "tmp_transterm/tmp",
                     "term_table": os.path.join(os.getcwd(), "tmp_term_table"),
                     "merge": os.path.join(os.getcwd(), "tmp_merge_gff"),
                     "gff": "tmp.gff",
                     "folder": os.path.join(os.getcwd(), "tmp")}
        self.suffixs = {"gff": "term.gff", "csv": "term.csv",
                        "allgff": "term_all.gff"}
        if args_term.srnas:
            self.srna_path = os.path.join(args_term.srnas, "tmp")
        else:
            self.srna_path = None
        self._make_gff_folder()

    def _combine_annotation(self, combine_file, files):
        with open(combine_file, 'w') as result:
            for file_ in files:
                check_start = False
                fh = open(file_, 'r')
                for line in fh:
                    if check_start:
                        result.write(line)
                    if "Location" in line:
                        check_start = True
                if "\n" not in line:
                    result.write("\n")
                fh.close()

    def _make_gff_folder(self):
        self.helper.check_make_folder(self.terms["all"])
        self.helper.check_make_folder(self.csvs["all"])
        self.helper.check_make_folder(self.terms["best"])
        self.helper.check_make_folder(self.csvs["best"])
        self.helper.check_make_folder(self.terms["express"])
        self.helper.check_make_folder(self.csvs["express"])
        self.helper.check_make_folder(self.terms["non"])
        self.helper.check_make_folder(self.csvs["non"])

    def _convert_gff2rntptt(self, gff_path, fasta_path, sRNAs):
        file_types = {}
        prefixs = []
        for gff in os.listdir(gff_path):
            if gff.endswith(".gff"):
                filename = gff.split("/")
                prefix = filename[-1][:-4]
                prefixs.append(prefix)
                gff_file = os.path.join(gff_path, gff)
                rnt_file = os.path.join(gff_path, gff.replace(".gff", ".rnt"))
                ptt_file = os.path.join(gff_path, gff.replace(".gff", ".ptt"))
                fasta = self.helper.get_correct_file(
                             fasta_path, ".fa", prefix, None, None)
                if not fasta:
                    print("Error: no proper file - {0}.fa".format(prefix))
                    sys.exit()
                if sRNAs:
                    self.multiparser.parser_gff(sRNAs, "sRNA")
                    srna = self.helper.get_correct_file(
                            self.srna_path, "_sRNA.gff", prefix, None, None)
                    if (srna) and (fasta):
                        self.converter.convert_gff2rntptt(
                            gff_file, fasta, ptt_file, rnt_file, srna,
                            srna.replace(".gff", ".rnt"))
                        file_types[prefix] = "srna"
                    if (not srna) and (fasta):
                        self.converter.convert_gff2rntptt(
                            gff_file, fasta, ptt_file, rnt_file, None, None)
                        file_types[prefix] = "normal"
                else:
                    self.converter.convert_gff2rntptt(
                        gff_file, fasta, ptt_file, rnt_file, None, None)
                    file_types[prefix] = "normal"
        return file_types, prefixs

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
                         os.path.join(srna_path,
                                      "_".join([prefix, "sRNA.rnt"]))]
                self._combine_annotation(combine_file, files)

    def _TransTermHP(self, fasta, file_, out_path, prefix, out, args_term):
        call([args_term.TransTermHP_path, "-p", args_term.expterm_path,
              fasta, os.path.join(self.combine_path, file_), "--t2t-perf",
              os.path.join(out_path, "_".join([
                  prefix,
                  "terminators_within_robust_tail-to-tail_regions.t2t"])),
              "--bag-output", os.path.join(out_path, "_".join([
                  prefix, "best_terminator_after_gene.bag"]))],
             stdout=out)

    def _run_TransTermHP(self, args_term):
        self.helper.check_make_folder(self.tmps["transterm"])
        for file_ in os.listdir(self.combine_path):
            if ".ptt" in file_:
                prefix = file_.replace(".ptt", "")
                fasta = self.helper.get_correct_file(
                             self.fasta_path, ".fa", prefix, None, None)
                if not fasta:
                    print("Error: no proper file - {0}.fa".format(prefix))
                    sys.exit()
                out_path = os.path.join(args_term.hp_folder, prefix)
                self.helper.check_make_folder(out_path)
                out = open(os.path.join(out_path,
                           "_".join([prefix, "terminators.txt"])), "w")
                self._TransTermHP(fasta, file_, out_path,
                                  prefix, out, args_term)
                out.close()
        shutil.rmtree(self.combine_path)

    def _convert_to_gff(self, prefixs, args_term):
        for prefix in prefixs:
            for folder in os.listdir(args_term.hp_folder):
                if prefix == folder:
                    out_path = os.path.join(args_term.hp_folder, folder)
                    for file_ in os.listdir(out_path):
                        if file_.endswith(".bag"):
                            out_file = os.path.join(
                                    self.tmps["transterm"],
                                    "_".join([prefix, self.tmps["hp_gff"]]))
                            self.converter.convert_transtermhp2gff(
                                 os.path.join(out_path, file_), out_file)
        self.multiparser.combine_gff(args_term.gffs, self.tmps["transterm"],
                                     None, self.tmps["hp"])

    def _combine_wigs(self, args_term):
        if (args_term.tex_wigs is not None) and (
                args_term.frag_wigs is not None):
            folder = args_term.tex_wigs.split("/")
            folder = "/".join(folder[:-1])
            merge_wigs = os.path.join(folder, "merge_wigs")
            self.helper.check_make_folder(merge_wigs)
            for wig in os.listdir(args_term.tex_wigs):
                if os.path.isdir(os.path.join(args_term.tex_wigs, wig)):
                    pass
                else:
                    shutil.copy(os.path.join(args_term.tex_wigs, wig),
                                merge_wigs)
            for wig in os.listdir(args_term.frag_wigs):
                if os.path.isdir(os.path.join(args_term.frag_wigs, wig)):
                    pass
                else:
                    shutil.copy(os.path.join(args_term.frag_wigs, wig),
                                merge_wigs)
        elif (args_term.tex_wigs is not None):
            merge_wigs = args_term.tex_wigs
        elif (args_term.frag_wigs is not None):
            merge_wigs = args_term.frag_wigs
        else:
            print("Error: no proper wig files!!!")
            sys.exit()
        return merge_wigs

    def _merge_sRNA(self, sRNAs, prefixs, gff_path):
        '''searching the terminator with sRNA information'''
        if sRNAs is not None:
            self.multiparser.parser_gff(sRNAs, "sRNA")
            self.helper.check_make_folder(self.tmps["merge"])
            for prefix in prefixs:
                tmp_gff = os.path.join(self.tmps["merge"], self.tmps["gff"])
                if self.tmps["gff"] in os.listdir(self.tmps["merge"]):
                    os.remove(tmp_gff)
                self.helper.merge_file(os.path.join(gff_path, prefix + ".gff"),
                                       tmp_gff)
                self.helper.merge_file(os.path.join(
                    self.srna_path, "_".join([prefix, "sRNA.gff"])), tmp_gff)
                self.helper.sort_gff(tmp_gff, os.path.join(
                    self.tmps["merge"], prefix + ".gff"))
                os.remove(tmp_gff)
            merge_path = self.tmps["merge"]
        else:
            merge_path = gff_path
        return merge_path

    def _move_file(self, term_outfolder, csv_outfolder):
        for gff in os.listdir(term_outfolder):
            if gff.endswith("_term.gff"):
                self.helper.sort_gff(os.path.join(term_outfolder, gff),
                                     self.tmps["gff"])
                shutil.move(self.tmps["gff"],
                            os.path.join(term_outfolder, gff))
                prefix = gff.replace("_term.gff", "")
                new_gff = os.path.join(self.terms["all"], "_".join([
                        prefix, self.suffixs["allgff"]]))
                csv_file = os.path.join(
                        os.path.join(self.csvs["all"], "_".join([
                            prefix, self.suffixs["csv"]])))
                out = open(new_gff, "w")
                out.write("##gff-version 3\n")
                out.close()
                self.helper.merge_file(
                        os.path.join(term_outfolder, gff),
                        os.path.join(
                            self.terms["all"], "_".join([
                                prefix, self.suffixs["allgff"]])))
                os.remove(os.path.join(term_outfolder, gff))
                pre_strain = ""
                if ("_".join([prefix, self.suffixs["csv"]]) in
                        os.listdir(self.csvs["all"])):
                    os.remove(csv_file)
                out_csv = open(csv_file, "w")
                out_csv.write("\t".join(["strain", "name", "start", "end",
                              "strand", "detect", "coverage_decrease",
                              "coverage_detail"]) + "\n")
                out_csv.close()
                fh = open(new_gff)
                for entry in self.gff_parser.entries(fh):
                    if entry.seq_id != pre_strain:
                        self.helper.merge_file(os.path.join(
                            self.tmps["term_table"], "_".join([
                                entry.seq_id, "term_raw.csv"])),
                            os.path.join(self.csvs["all"], "_".join([
                                prefix, self.suffixs["csv"]])))
                    pre_strain = entry.seq_id
                fh.close()

    def _run_rnafold(self, RNAfold_path, tmp_seq, tmp_sec, prefix):
        print("Computing secondray structure of {0}".format(prefix))
        self.helper.check_make_folder(self.tmps["folder"])
        pre_cwd = os.getcwd()
        os.chdir(self.tmps["folder"])
        os.system(" ".join([RNAfold_path, "<", os.path.join("..", tmp_seq),
                  ">", os.path.join("..", tmp_sec)]))
        os.chdir(pre_cwd)
        shutil.rmtree(self.tmps["folder"])

    def _compute_intersection_forward_reverse(
            self, prefixs, merge_path, wig_path, merge_wigs, args_term):
        '''the approach for searching gene converged region terminator'''
        for prefix in prefixs:
            tmp_seq = os.path.join(args_term.out_folder,
                                   "_".join(["inter_seq", prefix]))
            tmp_sec = os.path.join(args_term.out_folder,
                                   "_".join(["inter_sec", prefix]))
            tran_file = os.path.join(self.tran_path,
                                     "_".join([prefix, "transcript.gff"]))
            gff_file = os.path.join(merge_path, prefix + ".gff")
            print("Extracting seq of {0}".format(prefix))
            intergenic_seq(os.path.join(self.fasta_path, prefix + ".fa"),
                           tran_file, gff_file, tmp_seq)
            self._run_rnafold(args_term.RNAfold_path, tmp_seq, tmp_sec, prefix)
            tmp_cand = os.path.join(args_term.out_folder,
                                    "_".join(["term_candidates", prefix]))
            poly_t(tmp_seq, tmp_sec, gff_file, tran_file, tmp_cand, args_term)
            print("detection of terminator")
            detect_coverage(
                tmp_cand, os.path.join(merge_path, prefix + ".gff"),
                os.path.join(self.tran_path, "_".join([
                    prefix, "transcript.gff"])),
                os.path.join(self.fasta_path, prefix + ".fa"),
                os.path.join(wig_path, "_".join([prefix, "forward.wig"])),
                os.path.join(wig_path, "_".join([prefix, "reverse.wig"])),
                os.path.join(self.tmps["hp_path"], "_".join([
                    prefix, self.tmps["hp_gff"]])), merge_wigs,
                os.path.join(self.outfolder["term"], "_".join([
                    prefix, self.suffixs["gff"]])),
                os.path.join(self.tmps["term_table"], "_".join([
                    prefix, "term_raw.csv"])), args_term)
        self.multiparser.combine_gff(args_term.gffs, self.outfolder["term"],
                                     None, "term")
        self._move_file(self.outfolder["term"], self.outfolder["csv"])

    def _remove_tmp_file(self, merge_wigs, args_term):
        self.helper.remove_tmp(args_term.gffs)
        self.helper.remove_tmp(args_term.fastas)
        if args_term.srnas is not None:
            self.helper.remove_tmp(args_term.srnas)
            shutil.rmtree(self.tmps["merge"])
        if (args_term.tex_wigs is not None) and (
                args_term.frag_wigs is not None):
            shutil.rmtree(merge_wigs)
        self.helper.remove_tmp(args_term.trans)
        self.helper.remove_tmp(args_term.tex_wigs)
        self.helper.remove_tmp(args_term.frag_wigs)
        self.helper.remove_tmp(self.outfolder["term"])
        shutil.rmtree(self.tmps["transterm"])
        shutil.rmtree(self.tmps["term_table"])
        self.helper.remove_all_content(args_term.out_folder,
                                       "inter_seq_", "file")
        self.helper.remove_all_content(args_term.out_folder,
                                       "inter_sec_", "file")
        self.helper.remove_all_content(args_term.out_folder,
                                       "term_candidates_", "file")

    def _compute_stat(self, args_term):
        new_prefixs = []
        for gff in os.listdir(self.terms["all"]):
            if gff.endswith("_term_all.gff"):
                out_tmp = open(self.tmps["gff"], "w")
                out_tmp.write("##gff-version 3\n")
                new_prefix = gff.replace("_term_all.gff", "")
                new_prefixs.append(gff.replace("_term_all.gff", ""))
                num = 0
                fh = open(os.path.join(self.terms["all"], gff))
                for entry in self.gff_parser.entries(fh):
                    name = '%0*d' % (5, num)
                    entry.attributes["ID"] = "term" + str(num)
                    entry.attributes["Name"] = "_".join(["terminator_" + name])
                    entry.attribute_string = ";".join([
                        "=".join(items) for items in entry.attributes.items()])
                    out_tmp.write("\t".join([entry.info_without_attributes,
                                  entry.attribute_string]) + "\n")
                    num += 1
                out_tmp.close()
                fh.close()
                shutil.move(self.tmps["gff"], os.path.join(self.terms["all"],
                            "_".join([new_prefix, self.suffixs["gff"]])))
        if args_term.stat:
            stat_path = os.path.join(args_term.out_folder, "statistics")
            for prefix in new_prefixs:
                stat_term(os.path.join(self.terms["all"],
                          "_".join([prefix, self.suffixs["gff"]])),
                          os.path.join(self.csvs["all"],
                          "_".join([prefix, self.suffixs["csv"]])),
                          os.path.join(stat_path,
                          "_".join(["stat", prefix + ".csv"])),
                          os.path.join(self.terms["best"],
                          "_".join([prefix, "term"])),
                          os.path.join(self.terms["express"],
                          "_".join([prefix, "term"])),
                          os.path.join(self.terms["non"],
                          "_".join([prefix, "term"])))
                shutil.move(os.path.join(self.terms["best"],
                            "_".join([prefix, self.suffixs["csv"]])),
                            os.path.join(self.csvs["best"],
                            "_".join([prefix, self.suffixs["csv"]])))
                shutil.move(os.path.join(self.terms["express"],
                            "_".join([prefix, self.suffixs["csv"]])),
                            os.path.join(self.csvs["express"],
                            "_".join([prefix, self.suffixs["csv"]])))
                shutil.move(os.path.join(self.terms["non"],
                            "_".join([prefix, self.suffixs["csv"]])),
                            os.path.join(self.csvs["non"],
                            "_".join([prefix, self.suffixs["csv"]])))
                os.remove(os.path.join(self.terms["all"],
                          "_".join([prefix, self.suffixs["allgff"]])))

    def _check_gff_file(self, folder):
        for file_ in os.listdir(folder):
            if file_.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(folder, file_))

    def _compare_term_tran(self, args_term, prefixs):
        '''searching the associated terminator to transcript'''
        self.multiparser.combine_gff(args_term.gffs, self.tran_path,
                                     None, "transcript")
        prefixs = []
        for file_ in os.listdir(self.tran_path):
            if file_.endswith("_transcript.gff"):
                prefixs.append(file_.replace("_transcript.gff", ""))
        for type_ in ("best", "express", "all_candidates"):
            compare_term_tran(self.tran_path,
                              os.path.join(self.outfolder["term"], type_),
                              args_term.fuzzy_up_ta, args_term.fuzzy_down_ta,
                              args_term.out_folder, "terminator",
                              self.outfolder["term"], args_term.trans)
            for prefix in prefixs:
                shutil.move(
                    os.path.join(
                        args_term.out_folder, "statistics",
                        "stat_compare_transcript_terminator_" + prefix + ".csv"),
                    os.path.join(
                        args_term.out_folder, "statistics",
                        "_".join(["stat_compare_terminator_transcript", prefix,
                                  type_ + ".csv"])))

    def run_terminator(self, args_term):
        self._check_gff_file(args_term.gffs)
        self._check_gff_file(args_term.trans)
        self.multiparser.parser_fasta(args_term.fastas)
        if (not args_term.gffs) or (not args_term.fastas):
            print("Error: please assign gff annotation folder "
                  "and fasta folder!!!")
            sys.exit()
        file_types, prefixs = self._convert_gff2rntptt(
                self.gff_path, self.fasta_path, args_term.srnas)
        self._combine_ptt_rnt(self.gff_path, file_types, self.srna_path)
        self._run_TransTermHP(args_term)
        self._convert_to_gff(prefixs, args_term)
        self.helper.remove_tmp(self.gff_path)
        self.multiparser.parser_gff(args_term.trans, "transcript")
        self.helper.check_make_folder(self.tmps["term_table"])
        self.multiparser.parser_gff(self.tmps["transterm"], self.tmps["hp"])
        merge_path = self._merge_sRNA(args_term.srnas, prefixs, self.gff_path)
        self._compute_intersection_forward_reverse(
                prefixs, merge_path, args_term.wig_path,
                args_term.merge_wigs, args_term)
        self._compute_stat(args_term)
        self._compare_term_tran(args_term, prefixs)
        self._remove_tmp_file(args_term.merge_wigs, args_term)
