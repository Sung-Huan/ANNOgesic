#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
import shutil
from subprocess import call
from transaplib.multiparser import Multiparser
from transaplib.helper import Helper
from transaplib.sRNA_intergenic import Intergenic_sRNA
from transaplib.sRNA_utr_derived import UTR_derived_sRNA
from transaplib.merge_sRNA import Merge_sRNA_table, Merge_sRNA_GFF
from transaplib.extract_sRNA_info import Extract_energy, Extract_blast
from transaplib.plot_mountain import Plot_mountain_plot
from transaplib.sRNA_class import Classify_sRNA
from transaplib.gen_srna_output import Gen_sRNA_table, Gen_best_sRNA
from transaplib.blast_class import Blast_class
from transaplib.compare_sRNA_sORF import sRNA_sORF_comparison


class sRNA_detection(object):

    def __init__(self):
        self.helper = Helper()
        self.multiparpser = Multiparser()

    def _check_gff(self, gffs):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(gffs, gff))

    def _formatdb(self, database, type_, out_folder, blast_path):
        filenames = database.split(".")
        db_file = ".".join(filenames[0:-1])
        err = open(os.path.join(out_folder, "log.txt"), "w")
        call([os.path.join(blast_path, "makeblastdb"), "-in", database,
              "-dbtype", type_, "-out", db_file], stderr=err)

    def _run_normal(self, import_info, tss_path, prefix, bin_path, gff_path,
                    gff, tran, fuzzy_tss, max_len, min_len, wig_path, coverage, 
                    merge_wigs, libs, tex_notex, replicates, table_best,
                    decrease_inter, fuzzy_inter, out_folder):
        if ("1" in import_info):
            tss = self.helper.get_correct_file(tss_path, "_TSS.gff", prefix, None)
        else:
            tss = False
        Intergenic_sRNA(os.path.join(gff_path, gff), tran, tss, fuzzy_tss, max_len, 
                        min_len, os.path.join(wig_path, "_".join([prefix, "forward.wig"])), 
                        os.path.join(wig_path, "_".join([prefix, "reverse.wig"])), 
                        merge_wigs, libs, tex_notex, replicates, 
                        os.path.join(out_folder, "_".join(["tmp_normal", prefix])), 
                        os.path.join(out_folder, "_".join(["tmp_normal_table", prefix])), 
                        table_best, decrease_inter, fuzzy_inter, coverage)

    def _run_utrsrna(self, bin_path, gff_path, gff, tran, fuzzy_tss, merge_wigs,
                     max_len, min_len, wig_path, prefix, tss, pro, fasta_path,
                     libs, tex_notex, replicates, table_best, decrease_utr, 
                     fuzzy_utr, utr5_coverage, utr3_coverage, interCDS_coverage):
        UTR_derived_sRNA(os.path.join(gff_path, gff), tran, tss, 
                         os.path.join(wig_path, "_".join([prefix, "forward.wig"])), 
                         os.path.join(wig_path, "_".join([prefix, "reverse.wig"])),
                         pro, max_len, min_len, fuzzy_tss, fuzzy_utr,
                         os.path.join(fasta_path, prefix + ".fa"),
                         merge_wigs, libs, tex_notex, replicates,
                         os.path.join(out_folder, "_".join(["tmp_utrsrna", prefix])),
                         os.path.join(out_folder, "_".join(["tmp_utrsrna_table", prefix])), 
                         table_best, decrease_utr, utr3_coverage, utr5_coverage, 
                         interCDS_coverage)

    def _check_necessary_file(self, gffs, trans, tex_wigs, frag_wigs, utr_srna, 
                              tsss, pros, sORF, import_info, fastas):
        if (gffs is None) or \
           (trans is None) or \
           ((tex_wigs is None) and (frag_wigs is None)):
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
        if tsss is not False:
            self._check_gff(tsss)
            tss_path = os.path.join(tsss, "tmp")
            self.multiparser._parser_gff(tsss, "TSS")
            self.multiparser._combine_gff(gffs, tss_path, None, "TSS")
        if pros is not False:
            self._check_gff(pros)
        if sORF is not False:
            self._check_gff(sORF)
            sorf_path = os.path.join(sORF, "tmp")
            self.multiparser._parser_gff(sORF, "TSS")
            self.multiparser._combine_gff(gffs, sorf_path, None, "sORF")
        if utr_srna or ("2" in import_info) or \
           ("3" in import_info) or ("4" in import_info):
            if fastas is False:
                print("Error: lack required fasta files for UTR derived sRNA detection!!!!")
                sys.exit()
            fasta_path = os.path.join(fastas, "tmp")
            self.multiparser._parser_fasta(fastas)
            self.multiparser._combine_fasta(gffs, fasta_path, None)
            fasta_path = self._multiparser("fasta", fastas, "fasta", gffs, None)
        return (tss_path, sorf_path, fasta_path)

    def _merge_wig(self, tex_wigs, frag_wigs, out_folder):
        if (tex_wigs is not False) and (frag_wigs is not False):
            merge_wigs = os.path.join(out_folder, "merge_wigs")
            self._check_make_folder(out_folder, "merge_wigs")
            for wig in os.listdir(tex_wigs):
                shutil.copy(os.path.join(tex_wigs, wig), os.path.join(out_folder, "merge_wigs"))
            for wig in os.listdir(frag_wigs):
                shutil.copy(os.path.join(frag_wigs, wig), os.path.join(out_folder, "merge_wigs"))
        elif (tex_wigs is not False):
            merge_wigs = tex_wigs
        elif (frag_wigs is not False):
            merge_wigs = frag_wigs
        return merge_wigs

    def _merge_libs(self, tlibs, flibs):
        if (tlibs is False) and (flibs is False):
            print("Error: please input proper libraries!!")
        if (tlibs is not False) and (flibs is not False):
            libs = tlibs + flibs
        elif (tlibs is not False):
            libs = tlibs
        elif (flibs is not False):
            libs = flibs
        return libs

    def _run_program(self, gff_path, prefixs, mport_info, tran_path, tss_path, pro_path, 
                     wig_path, max_len, min_len, coverage, merge_wigs, libs, tex_notex, 
                     replicates, table_best, decrease_inter, decrease_utr, fuzzy_inter, 
                     fuzzy_utr, fuzzy_tss, utr5_coverage, utr3_coverage, interCDS_coverage,
                     out_folder):
        prefixs = []
        for gff in os.listdir(gff_path):
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                prefixs.append(prefix)
                print("Running sRNA detection of {0}....".format(prefix))
                tran = self._get_correct_file(tran_path, "_transcript.gff", prefix)
                self._run_normal(import_info, tss_path, prefix, bin_path, gff_path,
                                 gff, tran, fuzzy_tss, max_len, min_len, wig_path, coverage,
                                 merge_wigs, libs, tex_notex, replicates, table_best,
                                 decrease_inter, fuzzy_inter)
                if utr_srna:
                    print("Running UTR derived sRNA detection of {0}....".format(prefix))
                    tss = self._get_correct_file(tss_path, "_TSS.gff", prefix)
                    pro = self._get_correct_file(pro_path, "_processing.gff", prefix)
                    self._run_utrsrna(bin_path, gff_path, gff, tran, fuzzy_tss, merge_wigs,
                                      max_len, min_len, wig_path, prefix, tss, pro, fasta_path,
                                      libs, tex_notex, replicates, table_best, decrease_utr, 
                                      fuzzy_utr, utr5_coverage, utr3_coverage, interCDS_coverage)
                    os.system("rm median")
                    # merge intergenic and UTR derived sRNA   #
                    print("merging data of intergenic and UTR_derived sRNA...")
                    merge_gff = os.path.join(out_folder, "_".join(["tmp_merge", prefix]))
                    utr_gff = os.path.join(out_folder, "_".join(["tmp_utrsrna", prefix]))
                    normal_gff = os.path.join(out_folder, "_".join(["tmp_normal", prefix]))
                    merge_table = os.path.join(out_folder, "_".join(["tmp_merge_table", prefix]))
                    utr_table = os.path.join(out_folder, "_".join(["tmp_utrsrna_table", prefix]))
                    normal_table = os.path.join(out_folder, "_".join(["tmp_normal_table", prefix]))
                    Merge_sRNA_GFF(utr_gff, normal_gff, merge_gff)
                    Merge_sRNA_table(merge_gff, normal_table, utr_table,
                                     os.path.join(wig_path, "_".join([prefix, "forward.wig"])), 
                                     os.path.join(wig_path, "_".join([prefix, "reverse.wig"])),
                                     merge_wigs, libs, tex_notex, replicates, table_best, merge_table)
                else:
                    shutil.copyfile(normal_gff, merge_gff)
                    shutil.copyfile(normal_table, merge_table)
                self.helper.sort_gff(merge_gff, os.path.join(out_folder, "_".join(["tmp_basic", prefix])))
        return prefixs

    def _get_seq_sec(self, fasta_path, detect, out_folder, prefix, sec_path, seq_path, vienna_path):
        detect = False
        for fasta in os.listdir(fasta_path):
            if fasta.endswith(".fa") and \
               (fasta.replace(".fa", "") == prefix):
                detect = True
                break
        if detect:
            detect = False
            seq_file = os.path.join(out_folder, "_".join(["sRNA_seq", prefix]))
            sec_file = os.path.join(out_folder, "_".join(["sRNA_2d", prefix]))
            self.helper.get_seq(os.path.join(out_folder, "_".join(["tmp_basic", prefix])), 
                                os.path.join(fasta_path, fasta), seq_file)
        else:
            print("Error:There is not fasta file of {0}".format(prefix))
        self._check_make_folder(out_folder, "tmp_srna")
        tmp_path = os.path.join(out_folder, "tmp_srna")
        main_path = os.getcwd()
        os.chdir(tmp_path)
        sec_file = os.path.join(main_path, sec_file)
        seq_file = os.path.join(main_path, seq_file)
        tmp_sec_path = os.path.join(main_path, sec_path) ## the same as sec_path, just modify for fit os.chdir
        tmp_dot_path = os.path.join(main_path, dot_path)
        os.system(" ".join(["cat", seq_file, "|", os.path.join(vienna_path, "Progs", "RNAfold"),
                            "-p", ">", sec_file]))
        Extract_energy(os.path.join(main_path, out_folder, "_".join(["tmp_basic", prefix])), 
                       sec_file, "_".join(["tmp_energy", prefix]))
        return {"sec": tmp_sec_path, "dot": tmp_dot_path, "main": main_path}

    def _replot_sec_to_pdf(self, vienna_path, tmp_paths, ps2pdf14_path):
        for file_ in os.listdir(os.getcwd()): ### now we are in out_folder/tmp_srna
            if file_.endswith("ss.ps"):
                dot_file = file_.replace("ss.ps", "dp.ps")
                rel_file = file_.replace("ss.ps", "rss.ps")
                out = open(rel_file, "w")
                print("replot {0}".format(file_))
                call([os.path.join(vienna_path, "Utils", "relplot.pl"), 
                      file_, dot_file], stdout=out)
        for file_ in os.listdir(os.getcwd()):
            if file_.endswith(".ps"):
                pdf_file = file_.replace(".ps", ".pdf")
                print("convert {0} to pdf".format(file_))
                call([ps2pdf14_path, file_, pdf_file])
        os.mkdir(os.path.join(tmp_paths["sec"], prefix))
        os.mkdir(os.path.join(tmp_paths["dot"], prefix))
        self.helper.move_all_content(os.getcwd(), 
                 os.path.join(tmp_paths["sec"], prefix), "rss.pdf")
        self.helper.move_all_content(os.getcwd(), 
                 os.path.join(tmp_paths["dot"], prefix), "dp.pdf")

    def _plot_mountain(self, mountain, moun_path, main_path, prefix, vienna_path):
        if mountain is True:
            tmp_moun_path = os.path.join(main_path, moun_path)
            os.mkdir(os.path.join(tmp_moun_path, prefix))
            print("Generating mountain plot of {0}....".format(prefix))
            for dot_file in os.listdir(os.getcwd()):
                if dot_file.endswith("dp.ps"):
                    out = open("mountain.txt", "w")
                    moun_file = dot_file.replace("dp.ps", "mountain.pdf")
                    print("Generating {0}".format(moun_file))
                    call([os.path.join(vienna_path, "Utils", "mountain.pl"), 
                          dot_file], stdout=out)
                    Plot_mountain_plot("mountain.txt", moun_file)
                    os.rename(moun_file, os.path.join(tmp_moun_path, prefix, moun_file))
                    out.close()
            os.remove("mountain.txt")

    def _compute_2d_and_energy(self, out_folder, prefixs, fasta_path, vienna_path, 
                               mountain, ps2pdf14_path):
        print("Running energy calculation....")
        moun_path = os.path.join(out_folder, "mountain_plot")
        sec_path = os.path.join(out_folder, "sec_structure", "sec_plot")
        dot_path = os.path.join(out_folder, "sec_structure", "dot_plot")
        self.helper.remove_all_content(sec_path, None, "dir")
        self.helper.remove_all_content(dot_path, None, "dir")
        self.helper.remove_all_content(moun_path, None, "dir")
        for prefix in prefixs:
            tmp_paths = self._get_seq_sec(fasta_path, detect, out_folder,
                                          prefix, sec_path, seq_path, vienna_path)
            self._replot_sec_to_pdf(vienna_path, tmp_paths, ps2pdf14_path)
            self._plot_mountain(mountain, moun_path, main_path, prefix, vienna_path)
            self.helper.remove_all_content(os.getcwd(), ".ps", "file")
            self.helper.remove_all_content(os.getcwd(), "ss.pdf", "file")
            os.chdir(tmp_paths("main"))
            os.rename(os.path.join(out_folder, "tmp_srna", "_".join(["tmp_energy", prefix])), 
                      os.path.join(out_folder, "_".join(["tmp_basic", prefix])))
            shutil.rmtree(os.path.join(out_folder, "tmp_srna"))

    def _blast(self, database, database_format, data_type, out_folder, blast_path,
               prefixs, fasta_path, program, database_type):
        if (database is False):
            print("Error: No database assigned!")
        else:
            blast_file = os.path.join(out_folder, "blast_result_and_misc",
                         "_".join([database_type, "blast", prefix + ".txt"]))
            srna_file = os.path.join(out_folder, "_".join(["tmp_basic", prefix]))
            out_file = os.path.join(out_folder, "_".join(["tmp", database_type, prefix]))
            if database_format is True:
                self._formatdb(database, data_type, out_folder, blast_path)
            for prefix in prefixs:
                print("Running Blast of {0}".format(prefix))
                seq_file = os.path.join(out_folder, "_".join(["sRNA_seq", prefix]))
                if seq_file not in os.listdir(out_folder):
                    seq_out = open(seq_file, "w")
                    get_seq(srna_file, os.path.join(fasta_path, prefix + ".fa"), seq_file)
                call([os.path.join(blast_path, program), "-db", database,
                      "-evalue", str(0.0001), "-query", seq_file,
                      "-out", blast_file])
                Extract_blast(blast_file, srna_file, out_file, out_file + ".csv", database_type)
                os.rename(out_file, srna_file)

    def _class_sRNA(self, import_info, prefixs, gff_output, table_output, stat_path,
                    energy, nr_hit_num, max_len, min_len):
        if (len(import_info) != 1) or ('6' not in import_info):
            for prefix in prefixs:
                print("classifying sRNA of " + prefix)
                class_gff = os.path.join(gff_output, "for_class")
                class_table = os.path.join(table_output, "for_class")
                self._check_make_folder(class_table, prefix)
                self._check_make_folder(class_gff, prefix)
                class_gff = os.path.join(class_gff, prefix)
                class_table = os.path.join(class_table, prefix)
                out_stat = os.path.join(stat_path, "_".join(["stat_sRNA_class", prefix + ".csv"]))
                Classify_sRNA(os.path.join(gff_output, "all_candidates", "_".join([prefix, "sRNA.gff"])), 
                              class_gff, energy, hit_nr_num, out_stat)
                for srna in os.listdir(class_gff):
                    out_table = os.path.join(class_table, srna.replace(".gff", ".csv"))
                    Gen_sRNA_table(os.path.join(class_gff, srna), 
                                   os.path.join(out_folder, "_".join(["tmp_merge_table", prefix])), 
                                   os.path.join(out_folder, "_".join(["tmp_nr", prefix + ".csv"])), 
                                   os.path.join(out_folder, "_".join(["tmp_sRNA", prefix + ".csv"])),
                                   max_len, min_len, out_table)

    def _get_best_result(self, prefixs, gff_output, table_output, energy, nr_hit_num,
                         all_hit, best_sORF, max_len, min_len):
        for prefix in prefixs:
            best_gff = os.path.join(gff_output, "best", "_".join([prefix, "sRNA.gff"]))
            best_table = os.path.join(table_output, "best", "_".join([prefix, "sRNA.csv"]))
            Gen_best_sRNA(os.path.join(gff_output, "all_candidates", "_".join([prefix, "sRNA.gff"])), 
                          all_hit, energy, nr_hit_num, best_sORF, best_gff)
            Gen_sRNA_table(os.path.join(gff_output, "best", "_".join([prefix, "sRNA.gff"])),
                           os.path.join(out_folder, "_".join(["tmp_merge_table", prefix])),
                           os.path.join(out_folder, "_".join(["tmp_nr", prefix + ".csv"])),
                           os.path.join(out_folder, "_".join(["tmp_sRNA", prefix + ".csv"])),
                           max_len, min_len, best_table)

    def _remove_file(self, out_folder, fastas, gffs, wigs, tsss, trans, pros, sORF):
        self.helper.remove_all_content(out_folder, "tmp_", "dir")
        shutil.rmtree(os.path.join(out_folder, "merge_wigs"))
        self.helper.remove_tmp(fastas)
        self.helper.remove_tmp(gffs)
        self.helper.remove_tmp(wigs)
        self.helper.remove_tmp(trans)
        if tsss is not False:
            self.helper.remove_tmp(tsss)
        if pros is not False:
            self.helper.remove_tmp(pros)
        if sORF is not False:
            self.helper.remove_tmp(sORF)

    def _filter_sRNA(self, import_info, out_folder, prefixs, fasta_path, vienna_path, moutain,
                     ps2pdf14path, nr_database, database_format, blast_path, srna_database,
                     stat_path, sorf_path):
        if "2" in import_info:
            self._compute_2d_and_energy(out_folder, prefixs, fasta_path, vienna_path,
                                        mountain, ps2pdf14_path)
        if "3" in import_info:
            self._blast(nr_database, database_format, "prot", out_folder, blast_path,
                        prefixs, fasta_path, "blastx", "nr")
        if "4" in import_info:
            self._blast(srna_database, database_format, "nucl", out_folder, blast_path,
                        prefixs, fasta_path, "blastn", "sRNA")
            out_srna_blast = os.path.join(stat_path, "_".join(["stat_sRNA_blast_class",
                                          prefix + ".csv"]))
            Blast_class(os.path.join(out_folder, "_".join(["tmp_sRNA_", prefix + ".csv"])),
                                     out_srna_table)
        if "5" in import_info:
            if "_".join([prefix, "sORF.gff"]) in os.listdir(sorf_path):
                tmp_srna = os.path.join(out_folder, "".join(["tmp_srna_sorf", prefix]))
                tmp_sorf = os.path.join(out_folder, "".join(["tmp_sorf_srna", prefix]))
                sRNA_sORF_comparison(os.path.join(out_folder, "_".join(["tmp_basic", prefix])),
                                     os.path.join(sorf_path, "_".join([prefix, "sORF.gff"])),
                                     tmp_srna, tmp_sorf)
                os.remove(tmp_sorf)
                os.rename(tmp_srna, os.path.join(out_folder, "_".join(["tmp_basic", prefix])))

    def run_sRNA_detection(self, bin_path, out_folder, utr_srna, gffs, tsss, trans, fuzzy_tss,
                           import_info, tex_wigs, frag_wigs, pros, fastas, mountain, database_format,
                           srna_database, nr_database, energy, coverage, utr5_coverage, utr3_coverage, 
                           interCDS_coverage, max_len, min_len, tlibs, flibs, replicates, tex_notex, 
                           table_best, decrease_inter, decrease_utr, fuzzy_inter, fuzzy_utr, nr_hits_num, 
                           sORF, all_hit, best_sORF, ps2pdf14_path, vienna_path, blast_path):
        gff_output = os.path.join(out_folder, "gffs")
        table_output = os.path.join(out_folder, "tables")
        stat_path = os.path.join(out_folder, "statistics")
        multiparser._parser_gff(gffs, None)
        paths = self._check_necessary_file(gffs, trans, tex_wigs, frag_wigs, utr_srna, tsss, pros, sORF)
        tss_path = paths[0]
        sorf_path = paths[1]
        fasta_path = paths[2]
        self._merge_wig(tex_wigs, frag_wigs)
        wig_path = os.path.join(merge_wigs, "tmp")
        multiparser._parser_wig(merge_wigs)
        multiparser._combine_wig(gffs, wig_path, None)
        tran_path = os.path.join(trans, "tmp")
        multiparser._parser_gff(trans, "transcript")
        multiparser._combine_gff(gffs, tran_path, None, "transcript")
        libs = self._merge_libs(tlibs, flibs)
        prefixs = self._run_program(gffs, prefixs, mport_info, tran_path, tss_path, pro_path,             
                                    wig_path, max_len, min_len, coverage, merge_wigs, libs, tex_notex, 
                                    replicates, table_best, decrease_inter, decrease_utr, fuzzy_inter, 
                                    fuzzy_utr, fuzzy_tss, utr5_coverage, utr3_coverage, interCDS_coverage)
        self._filter_sRNA(import_info, out_folder, prefixs, fasta_path, vienna_path, moutain,
                     ps2pdf14path, nr_database, database_format, blast_path, srna_database,
                     stat_path, sorf_path)
        for prefix in prefixs:
            shutil.copyfile(os.path.join(out_folder, "_".join(["tmp_basic", prefix])),
                            os.path.join(gff_output, "all_candidates", "_".join([prefix, "sRNA.gff"])))
        for prefix in prefixs:
            out_table = os.path.join(table_output, "all_candidates", "_".join([prefix, "sRNA.csv"]))
            Gen_sRNA_table(os.path.join(gff_output, "all_candidates", "_".join([prefix, "sRNA.gff"])),
                           os.path.join(out_folder, "_".join(["tmp_merge_table", prefix])),
                           os.path.join(out_folder, "_".join(["tmp_nr", prefix + ".csv"])),
                           os.path.join(out_folder, "_".join(["tmp_sRNA", prefix + ".csv"])),
                           max_len, min_len, out_table)
        self._class_sRNA(import_info, prefixs, gff_output, table_output, stat_path,
                         energy, nr_hit_num, max_len, min_len)
        self._get_best_result(prefixs, gff_output, table_output, energy, nr_hit_num,
                              all_hit, best_sORF, max_len, min_len)
        self._remove_file(out_folder, fastas, gffs, wigs, tsss, trans, pros, sORF)
