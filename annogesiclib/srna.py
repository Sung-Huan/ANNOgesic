import os
import sys
import shutil
from subprocess import call
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.sRNA_intergenic import intergenic_srna
from annogesiclib.sRNA_utr_derived import utr_derived_srna
from annogesiclib.merge_sRNA import merge_srna_table, merge_srna_gff
from annogesiclib.extract_sRNA_info import extract_energy, extract_blast
from annogesiclib.plot_mountain import plot_mountain_plot
from annogesiclib.sRNA_class import classify_srna
from annogesiclib.gen_srna_output import gen_srna_table, gen_best_srna
from annogesiclib.blast_class import blast_class
from annogesiclib.compare_sRNA_sORF import srna_sorf_comparison
from annogesiclib.change_db_format import change_format 
from annogesiclib.compare_srna_term import compare_srna_term
from annogesiclib.print_rank_all import print_rank_all
from annogesiclib.sRNA_filter_frag import filter_frag
from annogesiclib.sRNA_filter_min_utr import filter_utr


class sRNADetection(object):

    def __init__(self, out_folder, tsss, pros, sorf_file, fastas, trans, terms):
        self.helper = Helper()
        self.multiparser = Multiparser()
        self.gff_output = os.path.join(out_folder, "gffs")
        self.table_output = os.path.join(out_folder, "tables")
        self.stat_path = os.path.join(out_folder, "statistics")
        self.tss_path = self._check_folder_exist(tsss)
        self.pro_path = self._check_folder_exist(pros)
        self.sorf_path = self._check_folder_exist(sorf_file)
        self.fasta_path = os.path.join(fastas, "tmp")
        self.tran_path = os.path.join(trans, "tmp")
        self.term_path = self._check_folder_exist(terms)
        self.merge_wigs = os.path.join(out_folder, "merge_wigs")
        self.prefixs = {"merge": os.path.join(out_folder, "tmp_merge"),
                        "utr": os.path.join(out_folder, "tmp_utrsrna"),
                        "normal": os.path.join(out_folder, "tmp_normal"),
                        "in_cds": os.path.join(out_folder, "tmp_incds"),
                        "merge_table": os.path.join(
                                       out_folder, "tmp_merge_table"),
                        "utr_table": os.path.join(
                                     out_folder, "tmp_utrsrna_table"),
                        "normal_table": os.path.join(
                                        out_folder, "tmp_normal_table"),
                        "in_cds_table": os.path.join(
                                        out_folder, "tmp_incds_table"),
                        "basic": os.path.join(out_folder, "tmp_basic"),
                        "energy": os.path.join(out_folder, "tmp_energy")}
        self.tmps = {"nr": os.path.join(out_folder, "tmp_nr"),
                     "srna": os.path.join(out_folder, "tmp_sRNA")}
        self.best_table = os.path.join(self.table_output, "best")
        self.table_output = os.path.join(out_folder, "tables")
        self.stat_path = os.path.join(out_folder, "statistics")
        self.all_best = {"all_gff": os.path.join(
                                    self.gff_output, "all_candidates"),
                         "best_gff": os.path.join(self.gff_output, "best"),
                         "all_table": os.path.join(
                                      self.table_output, "all_candidates"),
                         "best_table": os.path.join(self.table_output, "best")}

    def _check_folder_exist(self, folder):
        if folder is not None:
            path = os.path.join(folder, "tmp")
        else:
            path = None
        return path

    def _check_gff(self, gffs):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(gffs, gff))

    def _run_format(self, blast_path, database, type_, db_file, err):
        call([os.path.join(blast_path, "makeblastdb"), "-in", database,
              "-dbtype", type_, "-out", db_file], stderr=err)

    def _formatdb(self, database, type_, out_folder, blast_path, database_type):
        err = open(os.path.join(out_folder, "log.txt"), "w")
        if (database.endswith(".fa")) or (
            database.endswith(".fna")) or (
            database.endswith(".fasta")):
            pass
        else:
            folders = database.split("/")
            filename = folders[-1]
            folder = "/".join(folders[:-1])
            for fasta in os.listdir(folder):
                if (fasta.endswith(".fa")) or (
                    fasta.endswith(".fna")) or (
                    fasta.endswith(".fasta")):
                    if ".".join(fasta.split(".")[:-1]) == filename:
                        database = os.path.join(folder, fasta)
        if database_type == "sRNA":
            change_format(database, "tmp_srna_database")
            os.remove(database)
            shutil.move("tmp_srna_database", database)
        db_file = ".".join(database.split(".")[:-1])
        self._run_format(blast_path, database, type_, db_file, err)
        err.close()

    def _merge_frag_tex_file(self, frag_wigs, tex_wigs, frag_gff, tex_gff,
                             frag_csv, tex_csv, merge_table, merge_gff):
        if (frag_wigs is not None) and (tex_wigs is not None):
            self.helper.merge_file(frag_gff, tex_gff)
            self.helper.merge_file(frag_csv, tex_csv)
            shutil.move(tex_csv, merge_table)
            self.helper.sort_gff(tex_gff, merge_gff)
            os.remove(frag_csv)
            os.remove(frag_gff)
            os.remove(tex_gff)
        elif (frag_wigs is not None):
            shutil.move(frag_csv, merge_table)
            self.helper.sort_gff(frag_gff, merge_gff)
            os.remove(frag_gff)
        elif (tex_wigs is not None):
            shutil.move(tex_csv, merge_table)
            self.helper.sort_gff(tex_gff, merge_gff)
            os.remove(tex_gff)

    def _run_normal(self, import_info, tss_path, pro_path, prefix, gff_path,
                    gff, tran, fuzzy_tss, max_len, min_len, tex_path, frag_path,
                    coverage_tex, coverage_frag, tex_notex, replicates,
                    table_best, decrease_inter, fuzzy_inter, out_folder,
                    tolerance, in_cds, hypo, tex_wigs, frag_wigs,
                    tlibs, flibs, tss_source, coverage_notex):
        if "tmp_cutoff_inter" in os.listdir(out_folder):
            os.remove(os.path.join(out_folder, "tmp_cutoff_inter"))
        frag_gff = None
        frag_csv = None
        tex_gff = None
        tex_csv = None
        if ("tss" in import_info):
            tss = self.helper.get_correct_file(tss_path, "_TSS.gff",
                                               prefix, None, None)
        else:
            tss = None
        if pro_path is not None:
            pro = self.helper.get_correct_file(pro_path,
                              "_processing.gff", prefix, None, None)
        else:
            pro = None
        if frag_wigs is not None:
            frag_gff = os.path.join(out_folder, "_".join(["tmp_frag", prefix]))
            frag_csv = os.path.join(out_folder, "_".join(["tmp_frag_table", prefix]))
            intergenic_srna(os.path.join(gff_path, gff), tran, tss, pro, fuzzy_tss,
                            max_len, min_len, os.path.join(frag_path,
                            "_".join([prefix, "forward.wig"])),
                            os.path.join(frag_path,
                            "_".join([prefix, "reverse.wig"])),
                            frag_wigs, flibs, tex_notex, replicates,
                            frag_gff, frag_csv, table_best, decrease_inter,
                            fuzzy_inter, coverage_frag, tolerance, in_cds,
                            hypo, True, prefix, out_folder, "frag", None)
        if tex_wigs is not None:
            tex_gff = os.path.join(out_folder, "_".join(["tmp_tex", prefix]))
            tex_csv = os.path.join(out_folder, "_".join(["tmp_tex_table", prefix]))
            intergenic_srna(os.path.join(gff_path, gff), tran, tss, pro, fuzzy_tss,
                            max_len, min_len, os.path.join(tex_path,
                            "_".join([prefix, "forward.wig"])),
                            os.path.join(tex_path,
                            "_".join([prefix, "reverse.wig"])),
                            tex_wigs, tlibs, tex_notex, replicates,
                            tex_gff, tex_csv, table_best, decrease_inter,
                            fuzzy_inter, coverage_tex, tolerance, in_cds, hypo,
                            tss_source, prefix, out_folder, "tex", coverage_notex)
        merge_table = "_".join([self.prefixs["normal_table"], prefix])
        merge_gff = "_".join([self.prefixs["normal"], prefix])
        self._merge_frag_tex_file(frag_wigs, tex_wigs, frag_gff, tex_gff,
                             frag_csv, tex_csv, merge_table, merge_gff)
        if "TSS_class" in os.listdir(out_folder):
            tss = os.path.join(out_folder, "TSS_class", prefix + "_TSS.gff")
        return tss

    def _run_utrsrna(self, gff_path, gff, tran, fuzzy_tss, max_len, min_len,
                     prefix, tex_wigs, frag_wigs, tss, pro, fasta_path,
                     tex_notex, flibs, tlibs, replicates, table_best, min_utr,
                     decrease_utr, fuzzy_utr, utr_tex_cover, utr_frag_cover,
                     out_folder, hypo, tex_path, frag_path, utr_notex_cover):
        if "tmp_median" in os.listdir(out_folder):
            os.remove(os.path.join(out_folder, "tmp_median"))
        frag_gff = None
        frag_csv = None
        tex_gff = None
        tex_csv = None
        if tex_wigs is not None:
            tex_gff = os.path.join(out_folder, "_".join(["tmp_utr_tex", prefix]))
            tex_csv = os.path.join(out_folder, "_".join(["tmp_utr_tex_table", prefix]))
            utr_derived_srna(os.path.join(gff_path, gff), tran, tss,
                             os.path.join(tex_path,
                             "_".join([prefix, "forward.wig"])),
                             os.path.join(tex_path,
                             "_".join([prefix, "reverse.wig"])),
                             pro, max_len, min_len, fuzzy_tss, fuzzy_utr,
                             os.path.join(fasta_path, prefix + ".fa"),
                             tex_wigs, tlibs, tex_notex, replicates,
                             tex_gff, tex_csv, table_best, decrease_utr,
                             utr_tex_cover, hypo, out_folder, utr_notex_cover)
        if frag_wigs is not None:
            frag_gff = os.path.join(out_folder, "_".join(["tmp_utr_frag", prefix]))
            frag_csv = os.path.join(out_folder, "_".join(["tmp_utr_frag_table", prefix]))
            utr_derived_srna(os.path.join(gff_path, gff), tran, tss,
                             os.path.join(frag_path,
                             "_".join([prefix, "forward.wig"])),
                             os.path.join(frag_path,
                             "_".join([prefix, "reverse.wig"])),
                             pro, max_len, min_len, fuzzy_tss, fuzzy_utr,
                             os.path.join(fasta_path, prefix + ".fa"),
                             frag_wigs, flibs, tex_notex, replicates,
                             frag_gff, frag_csv, table_best, decrease_utr,
                             utr_frag_cover, hypo, out_folder, None)
        merge_table = "_".join([self.prefixs["utr_table"], prefix])
        merge_gff = "_".join([self.prefixs["utr"], prefix])
        self._merge_frag_tex_file(frag_wigs, tex_wigs, frag_gff, tex_gff,
                                  frag_csv, tex_csv, merge_table, merge_gff)
        filter_utr(merge_gff, merge_table, min_utr)

    def _check_necessary_file(self, gffs, trans, tex_wigs, frag_wigs, utr_srna,
                              tsss, pros, sorf_file, import_info, fastas, terms):
        if (gffs is None) or (trans is None) or (
            (tex_wigs is None) and (frag_wigs is None)):
            print("Error: lack required files!!!!")
            sys.exit()
        if utr_srna:
            if (tsss is None):
                print("Error: lack required TSS files for UTR derived sRNA detection!!!!")
                sys.exit()
            if (pros is None):
                print("Warning: lack Processing site files for UTR derived sRNA detection!!!")
                print("it may effect the results!!!!")
        self._check_gff(gffs)
        self._check_gff(trans)
        if tsss is not None:
            self._check_gff(tsss)
            self.multiparser.parser_gff(tsss, "TSS")
            self.multiparser.combine_gff(gffs, self.tss_path, None, "TSS")
        if pros is not None:
            self._check_gff(pros)
            self.multiparser.parser_gff(pros, "processing")
            self.multiparser.combine_gff(gffs, self.pro_path,
                                          None, "processing")
        if sorf_file is not None:
            self._check_gff(sorf_file)
            self.multiparser.parser_gff(sorf_file, "sORF")
            self.multiparser.combine_gff(gffs, self.sorf_path, None, "sORF")
        if utr_srna or ("sec_str" in import_info) or (
           "blast_nr" in import_info) or ("blast_srna" in import_info):
            if fastas is None:
                print("Error: lack required fasta files for UTR derived sRNA detection!!!!")
                sys.exit()
            self.multiparser.parser_fasta(fastas)
            self.multiparser.combine_fasta(gffs, self.fasta_path, None)
        if terms is not None:
            self._check_gff(terms)
            self.multiparser.parser_gff(terms, "term")
            self.multiparser.combine_gff(gffs, self.term_path, None, "term")
        else:
            self.term_path = None

    def _merge_wig(self, tex_wigs, frag_wigs, tex_path, frag_path, out_folder):
        if (tex_wigs is not None) and (frag_wigs is not None):
            self.helper.check_make_folder(self.merge_wigs)
            wig_path = os.path.join(self.merge_wigs, "tmp")
            self.helper.check_make_folder(wig_path)
            for wig in os.listdir(tex_wigs):
                if os.path.isfile(os.path.join(tex_wigs, wig)):
                    shutil.copy(os.path.join(tex_wigs, wig), self.merge_wigs)
            for wig in os.listdir(frag_wigs):
                if os.path.isfile(os.path.join(frag_wigs, wig)):
                    shutil.copy(os.path.join(frag_wigs, wig), self.merge_wigs)
            for wig in os.listdir(tex_path):
                if os.path.isfile(os.path.join(tex_path, wig)):
                    shutil.copy(os.path.join(tex_path, wig), wig_path)
            for wig in os.listdir(frag_path):
                if os.path.isfile(os.path.join(frag_path, wig)):
                    self.helper.merge_file(os.path.join(frag_path, wig),
                                           os.path.join(wig_path, wig))
        elif (tex_wigs is not None):
            self.merge_wigs = tex_wigs
        elif (frag_wigs is not None):
            self.merge_wigs = frag_wigs
        return self.merge_wigs, wig_path

    def _merge_libs(self, tlibs, flibs):
        if (tlibs is None) and (flibs is None):
            print("Error: please input proper libraries!!")
        if (tlibs is not None) and (flibs is not None):
            libs = tlibs + flibs
        elif (tlibs is not None):
            libs = tlibs
        elif (flibs is not None):
            libs = flibs
        return libs

    def _run_program(self, gff_path, import_info, tran_path, tss_path, pro_path,
                     max_len, min_len, tss_source, coverage_tex, coverage_notex,
                     coverage_frag, libs, tex_notex, replicates, table_best,
                     decrease_inter, decrease_utr, fuzzy_inter, fuzzy_utr,
                     fuzzy_tsss, utr_tex_cover, utr_notex_cover, utr_frag_cover,
                     out_folder, utr_srna, fasta_path, tolerance, in_cds,
                     cutoff_overlap, hypo, tex_path, frag_path, tex_wigs, frag_wigs,
                     tlibs, flibs, merge_wigs, wig_path, min_utr):
        prefixs = []
        tss = None
        for gff in os.listdir(gff_path):
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                prefixs.append(prefix)
                print("Running sRNA detection of {0}....".format(prefix))
                tran = self.helper.get_correct_file(tran_path,
                            "_transcript.gff", prefix, None, None)
                merge_gff = "_".join([self.prefixs["merge"], prefix])
                utr_gff = "_".join([self.prefixs["utr"], prefix])
                normal_gff = "_".join([self.prefixs["normal"], prefix])
                merge_table = "_".join([self.prefixs["merge_table"], prefix])
                utr_table = "_".join([self.prefixs["utr_table"], prefix])
                normal_table = "_".join([self.prefixs["normal_table"], prefix])
                tss = self._run_normal(import_info, tss_path, pro_path, prefix,
                                 gff_path, gff, tran, fuzzy_tsss["inter"],
                                 max_len, min_len, tex_path, frag_path,
                                 coverage_tex, coverage_frag, tex_notex,
                                 replicates, table_best, decrease_inter,
                                 fuzzy_inter, out_folder, tolerance, in_cds,
                                 hypo, tex_wigs, frag_wigs, tlibs, flibs,
                                 tss_source, coverage_notex)
                if utr_srna:
                    print("Running UTR derived sRNA detection of {0}...".format(
                          prefix))
                    if tss is None:
                        tss = self.helper.get_correct_file(tss_path, "_TSS.gff",
                                                           prefix, None, None)
                    if pro_path is not None:
                        pro = self.helper.get_correct_file(pro_path,
                                          "_processing.gff", prefix, None, None)
                    else:
                        pro = None
                    if tss is not None:
                        self._run_utrsrna(gff_path, gff, tran, fuzzy_tsss,
                                 max_len, min_len, prefix, tex_wigs, frag_wigs,
                                 tss, pro, fasta_path, tex_notex, flibs, tlibs,
                                 replicates, table_best, min_utr, decrease_utr,
                                 fuzzy_utr, utr_tex_cover, utr_frag_cover,
                                 out_folder, hypo, tex_path, frag_path, utr_notex_cover)
                self._merge_srna(utr_gff, normal_gff, merge_gff, normal_table,
                                 utr_table, wig_path, prefix, merge_wigs, libs,
                                 tex_notex, replicates, table_best, merge_table,
                                 os.path.join(gff_path, gff), in_cds,
                                 cutoff_overlap, tss, coverage_tex, coverage_notex,
                                 coverage_frag, utr_tex_cover, utr_notex_cover,
                                 utr_frag_cover, fuzzy_inter, out_folder)
                filter_frag(merge_table, merge_gff)
                self.helper.sort_gff(merge_gff,
                                     "_".join([self.prefixs["basic"], prefix]))
        return prefixs

    def _merge_srna(self, utr_gff, normal_gff, merge_gff, normal_table, utr_table,
                    wig_path, prefix, merge_wigs, libs, tex_notex, replicates,
                    table_best, merge_table, gff_file, in_cds, cutoff_overlap,
                    tss, coverage_tex, coverage_notex, coverage_frag, utr_tex_cover,
                    utr_notex_cover, utr_frag_cover, fuzzy_inter, out_folder):
        print("merging data of intergenic and UTR_derived sRNA...")
        merge_srna_gff(utr_gff, normal_gff, merge_gff, in_cds,
                       cutoff_overlap, gff_file)
        merge_srna_table(merge_gff, normal_table, utr_table,
                         os.path.join(wig_path,
                         "_".join([prefix, "forward.wig"])),
                         os.path.join(wig_path,
                         "_".join([prefix, "reverse.wig"])),
                         merge_wigs, libs, tex_notex, replicates,
                         table_best, merge_table, in_cds, tss, coverage_tex,
                         coverage_notex, coverage_frag, fuzzy_inter, out_folder,
                         utr_tex_cover, utr_notex_cover, utr_frag_cover)

    def _run_RNAfold(self, seq_file, vienna_path, sec_file):
        os.system(" ".join(["cat", seq_file, "|",
                  os.path.join(vienna_path, "RNAfold"),
                  "-p", ">", sec_file]))

    def _get_seq_sec(self, fasta_path, out_folder, prefix, sec_path,
                     dot_path, vienna_path):
        detect = False
        for fasta in os.listdir(fasta_path):
            if fasta.endswith(".fa") and (
               fasta.replace(".fa", "") == prefix):
                detect = True
                break
        if detect:
            detect = False
            seq_file = os.path.join(out_folder, "_".join(["sRNA_seq", prefix]))
            sec_file = os.path.join(out_folder, "_".join(["sRNA_2d", prefix]))
            self.helper.get_seq("_".join([self.prefixs["basic"], prefix]),
                                os.path.join(fasta_path, fasta), seq_file)
        else:
            print("Error:There is not fasta file of {0}".format(prefix))
        tmp_path = os.path.join(out_folder, "tmp_srna")
        self.helper.check_make_folder(tmp_path)
        main_path = os.getcwd()
        os.chdir(tmp_path)
        sec_file = os.path.join(main_path, sec_file)
        seq_file = os.path.join(main_path, seq_file)
        ## the same as sec_path, just modify for fit os.chdir
        tmp_sec_path = os.path.join(main_path, sec_path)
        tmp_dot_path = os.path.join(main_path, dot_path)
        self._run_RNAfold(seq_file, vienna_path, sec_file)
        extract_energy(os.path.join(main_path,
                       "_".join([self.prefixs["basic"], prefix])),
                       sec_file, os.path.join(main_path,
                       "_".join([self.prefixs["energy"], prefix])))
        for ps in os.listdir(os.getcwd()):
            new_ps = ps.replace("|", "_")
            shutil.move(ps, new_ps)
        return {"sec": tmp_sec_path, "dot": tmp_dot_path, "main": main_path,
                "tmp": os.path.join(main_path, tmp_path)}

    def _run_replot(self, vienna_util, tmp_paths, file_, dot_file, rel_file):
        os.system(" ".join([os.path.join(vienna_util, "relplot.pl"),
                  os.path.join(tmp_paths["tmp"], file_),
                  os.path.join(tmp_paths["tmp"], dot_file),
                  ">", os.path.join(tmp_paths["tmp"], rel_file)]))

    def _convert_pdf(self, ps2pdf14_path, tmp_paths, file_, pdf_file):
        call([ps2pdf14_path, os.path.join(tmp_paths["tmp"], file_), pdf_file])

    def _replot_sec_to_pdf(self, vienna_util, tmp_paths, ps2pdf14_path, prefix):
        ### now we are in out_folder/tmp_srna
        for file_ in os.listdir(os.getcwd()):
            if file_.endswith("ss.ps"):
                dot_file = file_.replace("ss.ps", "dp.ps")
                rel_file = file_.replace("ss.ps", "rss.ps")
                print("replot {0}".format(file_))
                self._run_replot(vienna_util, tmp_paths, file_,
                                 dot_file, rel_file)
        for file_ in os.listdir(tmp_paths["tmp"]):
            if (file_.endswith("rss.ps")) or (file_.endswith("dp.ps")):
                pdf_file = file_.replace(".ps", ".pdf")
                print("convert {0} to pdf".format(file_))
                self._convert_pdf(ps2pdf14_path, tmp_paths,
                                  file_, pdf_file)
        os.mkdir(os.path.join(tmp_paths["sec"], prefix))
        os.mkdir(os.path.join(tmp_paths["dot"], prefix))
        self.helper.move_all_content(tmp_paths["tmp"],
                 os.path.join(tmp_paths["sec"], prefix), ["rss.pdf"])
        self.helper.move_all_content(tmp_paths["tmp"],
                 os.path.join(tmp_paths["dot"], prefix), ["dp.pdf"])

    def _run_mountain(self, vienna_util, tmp_paths, dot_file, out):
        call([os.path.join(vienna_util, "mountain.pl"),
              os.path.join(tmp_paths["tmp"], dot_file)], stdout=out)

    def _plot_mountain(self, mountain, moun_path,
                       tmp_paths, prefix, vienna_util):
        if mountain:
            tmp_moun_path = os.path.join(tmp_paths["main"], moun_path)
            os.mkdir(os.path.join(tmp_moun_path, prefix))
            txt_path = os.path.join(tmp_paths["tmp"], "tmp_txt")
            self.helper.check_make_folder(txt_path)
            print("Generating mountain plot of {0}....".format(prefix))
            for dot_file in os.listdir(tmp_paths["tmp"]):
                if dot_file.endswith("dp.ps"):
                    moun_txt = os.path.join(tmp_paths["tmp"], "mountain.txt")
                    out = open(moun_txt, "w")
                    moun_file = dot_file.replace("dp.ps", "mountain.pdf")
                    print("Generating {0}".format(moun_file))
                    self._run_mountain(vienna_util, tmp_paths, dot_file, out)
                    plot_mountain_plot(moun_txt, moun_file)
                    shutil.move(moun_file,
                                os.path.join(tmp_moun_path, prefix, moun_file))
                    out.close()
                    os.remove(moun_txt)

    def _compute_2d_and_energy(self, out_folder, prefixs, fasta_path, vienna_path,
                               vienna_util, mountain, ps2pdf14_path):
        print("Running energy calculation....")
        moun_path = os.path.join(out_folder, "mountain_plot")
        sec_path = os.path.join(out_folder, "sec_structure", "sec_plot")
        dot_path = os.path.join(out_folder, "sec_structure", "dot_plot")
        self.helper.remove_all_content(sec_path, None, "dir")
        self.helper.remove_all_content(dot_path, None, "dir")
        self.helper.remove_all_content(moun_path, None, "dir")
        for prefix in prefixs:
            tmp_paths = self._get_seq_sec(fasta_path, out_folder, prefix,
                                          sec_path, dot_path, vienna_path)
            self._replot_sec_to_pdf(vienna_util, tmp_paths,
                                    ps2pdf14_path, prefix)
            self._plot_mountain(mountain, moun_path, tmp_paths,
                                prefix, vienna_util)
            self.helper.remove_all_content(os.getcwd(), ".ps", "file")
            os.chdir(tmp_paths["main"])
            shutil.move("_".join([self.prefixs["energy"], prefix]),
                      "_".join([self.prefixs["basic"], prefix]))
            shutil.rmtree(os.path.join(out_folder, "tmp_srna"))

    def _run_blast(self, blast_path, program, database, e, seq_file, blast_file,
                   strand):
        call([os.path.join(blast_path, program), "-db", database,
              "-evalue", str(e), "-strand", strand, "-query", seq_file,
              "-out", blast_file])

    def _get_strand_fasta(self, seq_file, out_folder):
        tmp_plus = os.path.join(out_folder, "tmp_plus.fa")
        tmp_minus = os.path.join(out_folder, "tmp_minus.fa")
        out_p = open(tmp_plus, "w")
        out_m = open(tmp_minus, "w")
        strand = ""
        with open(seq_file) as sh:
            for line in sh:
                line = line.strip()
                if line.startswith(">"):
                    if line[-1] == "+":
                        out_p.write(line + "\n")
                        strand = "plus"
                    elif line[-1] == "-":
                        out_m.write(line + "\n")
                        strand = "minus"
                else:
                    if strand == "plus":
                        out_p.write(line + "\n")
                    elif strand == "minus":
                        out_m.write(line + "\n")
        out_p.close()
        out_m.close()    
        return tmp_plus, tmp_minus

    def _blast(self, database, database_format, data_type, out_folder,
               blast_path, prefixs, fasta_path, program, database_type, e):
        if (database is None):
            print("Error: No database assigned!")
        else:
            if database_format:
                self._formatdb(database, data_type, out_folder, blast_path,
                               database_type)
            for prefix in prefixs:
                blast_file = os.path.join(out_folder, "blast_result_and_misc",
                             "_".join([database_type, "blast",
                                       prefix + ".txt"]))
                srna_file = "_".join([self.prefixs["basic"], prefix])
                out_file = os.path.join(out_folder,
                           "_".join(["tmp", database_type, prefix]))
                print("Running Blast of {0}".format(prefix))
                seq_file = os.path.join(out_folder,
                           "_".join(["sRNA_seq", prefix]))
                if seq_file not in os.listdir(out_folder):
                    self.helper.get_seq(srna_file,
                         os.path.join(fasta_path, prefix + ".fa"), seq_file)
                if database_type == "nr":
                    tmp_plus, tmp_minus = self._get_strand_fasta(seq_file, out_folder)
                    tmp_blast = os.path.join("tmp_blast.txt")
                    self._run_blast(blast_path, program, database, e,
                                    tmp_plus, tmp_blast, "plus")
                    self._run_blast(blast_path, program, database, e,
                                    tmp_minus, blast_file, "minus")
                    self.helper.merge_file(tmp_blast, blast_file)
                    os.remove(tmp_blast)
                    os.remove(tmp_plus)
                    os.remove(tmp_minus)
                else:
                    self._run_blast(blast_path, program, database, e,
                                    seq_file, blast_file, "both")
                extract_blast(blast_file, srna_file, out_file,
                              out_file + ".csv", database_type)
                shutil.move(out_file, srna_file)

    def _class_srna(self, import_info, prefixs, gff_output, table_output, stat_path,
                    energy, nr_hit_num, max_len, min_len, in_cds):
        if (len(import_info) != 1) or (len(import_info) != 0):
            for prefix in prefixs:
                print("classifying sRNA of {0}".format(prefix))
                class_gff = os.path.join(gff_output, "for_class")
                class_table = os.path.join(table_output, "for_class")
                self.helper.check_make_folder(os.path.join(class_table, prefix))
                self.helper.check_make_folder(os.path.join(class_gff, prefix))
                class_gff = os.path.join(class_gff, prefix)
                class_table = os.path.join(class_table, prefix)
                self.helper.check_make_folder(class_table)
                self.helper.check_make_folder(class_gff)
                out_stat = os.path.join(stat_path,
                           "_".join(["stat_sRNA_class", prefix + ".csv"]))
                classify_srna(os.path.join(self.all_best["all_gff"],
                              "_".join([prefix, "sRNA.gff"])),
                              class_gff, energy, nr_hit_num, out_stat, in_cds)
                for srna in os.listdir(class_gff):
                    out_table = os.path.join(class_table,
                                srna.replace(".gff", ".csv"))
                    print(os.path.join(class_gff, srna))
                    gen_srna_table(os.path.join(class_gff, srna),
                        "_".join([self.prefixs["merge_table"], prefix]),
                        "_".join([self.tmps["nr"], prefix + ".csv"]),
                        "_".join([self.tmps["srna"], prefix + ".csv"]),
                        max_len, min_len, out_table)

    def _get_best_result(self, prefixs, energy, nr_hit_num, all_hit,
                         best_sorf, max_len, min_len, best_term):
        for prefix in prefixs:
            best_gff = os.path.join(self.all_best["best_gff"],
                                    "_".join([prefix, "sRNA.gff"]))
            best_table = os.path.join(self.all_best["best_table"],
                                      "_".join([prefix, "sRNA.csv"]))
            gen_best_srna(os.path.join(self.all_best["all_gff"],
                                       "_".join([prefix, "sRNA.gff"])),
                          all_hit, energy, nr_hit_num, best_sorf, best_gff,
                          best_term)
            gen_srna_table(os.path.join(self.all_best["best_gff"],
                           "_".join([prefix, "sRNA.gff"])),
                           "_".join([self.prefixs["merge_table"], prefix]),
                           "_".join([self.tmps["nr"], prefix + ".csv"]),
                           "_".join([self.tmps["srna"], prefix + ".csv"]),
                           max_len, min_len, best_table)

    def _remove_file(self, out_folder, fastas, gffs, frag_wigs, tex_wigs,
                     tsss, trans, pros, sorf_file, merge_wigs):
        self.helper.remove_all_content(out_folder, "tmp_", "dir")
        self.helper.remove_all_content(out_folder, "tmp_", "file")
        self.helper.remove_tmp(fastas)
        self.helper.remove_tmp(gffs)
        if frag_wigs is not None:
            self.helper.remove_tmp(frag_wigs)
        if tex_wigs is not None:
            self.helper.remove_tmp(tex_wigs)
        if (frag_wigs is not None) and (
            tex_wigs is not None):
            shutil.rmtree(merge_wigs)
        self.helper.remove_tmp(trans)
        if tsss is not None:
            self.helper.remove_tmp(tsss)
        if pros is not None:
            self.helper.remove_tmp(pros)
        if sorf_file is not None:
            self.helper.remove_tmp(sorf_file)
        if "tmp_median" in os.listdir(out_folder):
            os.remove(os.path.join(out_folder, "tmp_median"))

    def _filter_srna(self, import_info, out_folder, prefixs, fasta_path,
                     vienna_path, vienna_util, mountain, ps2pdf14_path,
                     nr_database, srna_format, nr_format, blast_path,
                     srna_database, stat_path, sorf_path, e_nr, e_srna):
        if "sec_str" in import_info:
            self._compute_2d_and_energy(out_folder, prefixs, fasta_path,
                                        vienna_path, vienna_util,
                                        mountain, ps2pdf14_path)
        if "blast_nr" in import_info:
            self._blast(nr_database, nr_format, "prot", out_folder,
                        blast_path, prefixs, fasta_path, "blastx", "nr", e_nr)
        if "blast_srna" in import_info:
            self._blast(srna_database, srna_format, "nucl", out_folder,
                        blast_path, prefixs, fasta_path,
                        "blastn", "sRNA", e_srna)
            for prefix in prefixs:
                out_srna_blast = os.path.join(stat_path,
                                 "_".join(["stat_sRNA_blast_class",
                                 prefix + ".csv"]))
                blast_class("_".join([self.tmps["srna"], prefix + ".csv"]),
                            out_srna_blast)
        if "sorf" in import_info:
            if "_".join([prefix, "sORF.gff"]) in os.listdir(sorf_path):
                tmp_srna = os.path.join(out_folder,
                           "".join(["tmp_srna_sorf", prefix]))
                tmp_sorf = os.path.join(out_folder,
                           "".join(["tmp_sorf_srna", prefix]))
                srna_sorf_comparison("_".join([self.prefixs["basic"], prefix]),
                                     os.path.join(sorf_path,
                                     "_".join([prefix, "sORF.gff"])),
                                     tmp_srna, tmp_sorf)
                os.remove(tmp_sorf)
                shutil.move(tmp_srna, "_".join([self.prefixs["basic"], prefix]))

    def _get_replicates(self, replicates_tex, replicates_frag):
        if (replicates_tex is not None) and (
            replicates_frag is not None):
            replicates = {"tex": int(replicates_tex),
                          "frag": int(replicates_frag)}
        elif replicates_tex is not None:
            replicates = {"tex": int(replicates_tex), "frag": -1}
        elif replicates_frag is not None:
            replicates = {"tex": -1, "frag": int(replicates_frag)}
        else:
            print("Error:No replicates number assign!!!")
            sys.exit()
        return replicates

    def _import_info_format(self, import_info):
        new_info = []
        for info in import_info:
            info = info.lower()
            new_info.append(info)
        return new_info

    def _gen_table(self, prefixs, max_len, min_len, fuzzy_b, fuzzy_a,
                   import_info):
        for prefix in prefixs:
            out_table = os.path.join(self.all_best["all_table"],
                                     "_".join([prefix, "sRNA.csv"]))
            gen_srna_table(os.path.join(self.all_best["all_gff"],
                           "_".join([prefix, "sRNA.gff"])),
                           "_".join([self.prefixs["merge_table"], prefix]),
                           "_".join([self.tmps["nr"], prefix + ".csv"]),
                           "_".join([self.tmps["srna"], prefix + ".csv"]),
                           max_len, min_len, out_table)
            if ("term" in import_info) and (self.term_path is not None):
                compare_srna_term(os.path.join(self.all_best["all_gff"],
                                  "_".join([prefix, "sRNA.gff"])),
                                  out_table, os.path.join(self.term_path,
                                  "_".join([prefix, "term.gff"])),
                                  fuzzy_b, fuzzy_a)

    def _print_rank_all(self, prefixs):
        for prefix in prefixs:
            all_table = os.path.join(self.all_best["all_table"],
                                     "_".join([prefix, "sRNA.csv"]))
            best_table = os.path.join(self.all_best["best_table"],
                                      "_".join([prefix, "sRNA.csv"]))
            print_rank_all(all_table, best_table)

    def _parser_combine_wigs(self, libs, gffs, tex_wigs, frag_wigs, out_folder):
        tex_path = None
        frag_path = None
        if (tex_wigs is not None):
            tex_path = os.path.join(tex_wigs, "tmp")
            self.multiparser.parser_wig(tex_wigs)
            self.multiparser.combine_wig(gffs, tex_path, None, libs)
            merge_wigs = tex_wigs
            wig_path = tex_path
        if frag_wigs is not None:
            frag_path = os.path.join(frag_wigs, "tmp")
            self.multiparser.parser_wig(frag_wigs)
            self.multiparser.combine_wig(gffs, frag_path, None, libs)
            merge_wigs = frag_wigs
            wig_path = frag_path
        if (tex_path is not None) and (frag_path is not None):
            merge_wigs, wig_path = self._merge_wig(tex_wigs, frag_wigs,
                                        tex_path, frag_path, out_folder)
        return tex_path, frag_path, merge_wigs, wig_path

    def _filter_min_utr(self, prefixs, min_utr):
        for prefix in prefixs:
            filter_utr(os.path.join(self.all_best["all_gff"],
                                    "_".join([prefix, "sRNA.gff"])),
                       os.path.join(self.all_best["all_table"],
                                    "_".join([prefix, "sRNA.csv"])), min_utr)

    def run_srna_detection(self, vienna_path, vienna_util, blast_path,
            ps2pdf14_path, out_folder, utr_srna, gffs, tsss, trans,
            fuzzy_inter_tss, fuzzy_5utr_tss, fuzzy_3utr_tss, fuzzy_intercds_tss,
            import_info, tex_wigs, frag_wigs, pros, fastas, mountain,
            nr_format, srna_format, srna_database, nr_database, energy, coverage_tex,
            coverage_notex, coverage_frag, tolerance, utr_tex_cover,
            utr_notex_cover, utr_frag_cover, max_len, min_len, tlibs, flibs,
            replicates_tex, replicates_frag, tex_notex, e_nr, e_srna, in_cds,
            table_best, decrease_inter, decrease_utr, fuzzy_inter, fuzzy_utr,
            nr_hits_num, sorf_file, all_hit, best_sorf, cutoff_overlap, terms,
            fuzzy_b, fuzzy_a, best_term, hypo, tss_source, min_utr):
        replicates = self._get_replicates(replicates_tex, replicates_frag)
        fuzzy_tsss = {"5utr": fuzzy_5utr_tss, "3utr": fuzzy_3utr_tss,
                      "interCDS": fuzzy_intercds_tss, "inter": fuzzy_inter_tss}
        self.multiparser.parser_gff(gffs, None)
        self._check_necessary_file(gffs, trans, tex_wigs, frag_wigs,
             utr_srna, tsss, pros, sorf_file, import_info, fastas, terms)
        libs = self._merge_libs(tlibs, flibs)
        tex_path, frag_path, merge_wigs, wig_path  = self._parser_combine_wigs(
                                                     libs, gffs, tex_wigs,
                                                     frag_wigs, out_folder)
        self.multiparser.parser_gff(trans, "transcript")
        self.multiparser.combine_gff(gffs, self.tran_path, None, "transcript")
        import_info = self._import_info_format(import_info)
        prefixs = self._run_program(gffs, import_info, self.tran_path,
                      self.tss_path, self.pro_path, max_len, min_len, tss_source,
                      coverage_tex, coverage_notex, coverage_frag, libs, tex_notex,
                      replicates, table_best, decrease_inter, decrease_utr,
                      fuzzy_inter, fuzzy_utr, fuzzy_tsss, utr_tex_cover,
                      utr_notex_cover, utr_frag_cover, out_folder, utr_srna,
                      self.fasta_path, tolerance, in_cds, cutoff_overlap, hypo,
                      tex_path, frag_path, tex_wigs, frag_wigs, tlibs, flibs,
                      merge_wigs, wig_path, min_utr)
        self._filter_srna(import_info, out_folder, prefixs, self.fasta_path,
            vienna_path, vienna_util, mountain, ps2pdf14_path, nr_database,
            srna_format, nr_format, blast_path, srna_database, self.stat_path,
            self.sorf_path, e_nr, e_srna)
        for prefix in prefixs:
            shutil.copyfile("_".join([self.prefixs["basic"], prefix]),
                   os.path.join(self.all_best["all_gff"], "_".join([prefix, "sRNA.gff"])))
        self._gen_table(prefixs, max_len, min_len, fuzzy_b, fuzzy_a, import_info)
        self._class_srna(import_info, prefixs, self.gff_output, self.table_output,
             self.stat_path, energy, nr_hits_num, max_len, min_len, in_cds)
        self._get_best_result(prefixs, energy, nr_hits_num,
                              all_hit, best_sorf, max_len, min_len, best_term)
        self._print_rank_all(prefixs)
        self._remove_file(out_folder, fastas, gffs, frag_wigs, tex_wigs,
             tsss, trans, pros, sorf_file, merge_wigs)
