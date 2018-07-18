import os, gc
import sys
import shutil
import time
from subprocess import call, Popen
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.sRNA_intergenic import intergenic_srna
from annogesiclib.TSS_upstream import upstream
from annogesiclib.sRNA_utr_derived import utr_derived_srna
from annogesiclib.merge_sRNA import merge_srna_gff
from annogesiclib.merge_sRNA import merge_srna_table
from annogesiclib.extract_sRNA_info import extract_energy, extract_blast
from annogesiclib.plot_mountain import plot_mountain_plot
from annogesiclib.sRNA_class import classify_srna
from annogesiclib.gen_srna_output import gen_srna_table, gen_best_srna
from annogesiclib.blast_class import blast_class
from annogesiclib.compare_sRNA_sORF import srna_sorf_comparison
from annogesiclib.change_db_format import change_format
from annogesiclib.compare_srna_term import compare_srna_term
from annogesiclib.compare_srna_promoter import compare_srna_promoter
from annogesiclib.print_rank_all import print_rank_all
from annogesiclib.sRNA_filter_frag import filter_frag
from annogesiclib.sRNA_filter_min_utr import filter_utr
from annogesiclib.sRNA_antisense import srna_antisense
from annogesiclib.args_container import ArgsContainer
from annogesiclib.lib_reader import read_wig, read_libs
from annogesiclib.extract_sec_info import extract_info_sec, modify_header
from annogesiclib.get_srna_poly_u import get_srna_poly_u
from annogesiclib.reorganize_table import reorganize_table
from annogesiclib.check_srna_overlap import check_overlap


class sRNADetection(object):
    '''detection of sRNA'''

    def __init__(self, args_srna):
        self.args_container = ArgsContainer()
        self.helper = Helper()
        self.multiparser = Multiparser()
        self.gff_output = os.path.join(args_srna.out_folder, "gffs")
        self.table_output = os.path.join(args_srna.out_folder, "tables")
        self.stat_path = os.path.join(args_srna.out_folder, "statistics")
        self.tss_path = self._check_folder_exist(args_srna.tss_folder)
        self.pro_path = self._check_folder_exist(args_srna.pro_folder)
        self.sorf_path = self._check_folder_exist(args_srna.sorf_file)
        self.fasta_path = self._check_folder_exist(args_srna.fastas)
        self.tran_path = os.path.join(args_srna.trans, "tmp")
        self.term_path = self._check_folder_exist(args_srna.terms)
        self.merge_wigs = os.path.join(args_srna.out_folder, "merge_wigs")
        self.prefixs = {"merge": os.path.join(
                            args_srna.out_folder, "tmp_merge"),
                        "utr": os.path.join(
                            args_srna.out_folder, "tmp_utrsrna"),
                        "normal": os.path.join(
                            args_srna.out_folder, "tmp_normal"),
                        "in_cds": os.path.join(
                            args_srna.out_folder, "tmp_incds"),
                        "merge_table": os.path.join(
                            args_srna.out_folder, "tmp_merge_table"),
                        "utr_table": os.path.join(
                            args_srna.out_folder, "tmp_utrsrna_table"),
                        "normal_table": os.path.join(
                            args_srna.out_folder, "tmp_normal_table"),
                        "in_cds_table": os.path.join(
                            args_srna.out_folder, "tmp_incds_table"),
                        "basic": os.path.join(
                            args_srna.out_folder, "tmp_basic"),
                        "energy": os.path.join(
                            args_srna.out_folder, "tmp_energy")}
        self.tmps = {"nr": os.path.join(args_srna.out_folder, "tmp_nr"),
                     "srna": os.path.join(args_srna.out_folder, "tmp_sRNA")}
        self.best_table = os.path.join(self.table_output, "best_candidates")
        self.table_output = os.path.join(args_srna.out_folder, "tables")
        self.stat_path = os.path.join(args_srna.out_folder, "statistics")
        self.all_best = {"all_gff": os.path.join(
                             self.gff_output, "all_candidates"),
                         "best_gff": os.path.join(self.gff_output,
                             "best_candidates"),
                         "all_table": os.path.join(
                             self.table_output, "all_candidates"),
                         "best_table": os.path.join(self.table_output,
                             "best_candidates")}

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

    def _run_format(self, blastdb, database, type_, db_file, err, log):
        log.write("Please make sure the version of BLAST+ is at least 2.2.28+.\n")
        log.write(" ".join([blastdb, "-in", database,
                            "-dbtype", type_, "-out", db_file]) + "\n")
        call([blastdb, "-in", database,
              "-dbtype", type_, "-out", db_file], stderr=err)
        log.write("Done!\n")

    def _wait_process(self, processes):
        for p in processes:
            p.wait()
            if p.stdout:
                p.stdout.close()
            if p.stdin:
                p.stdin.close()
            if p.stderr:
                p.stderr.close()
            try:
                p.kill()
            except OSError:
                pass
        time.sleep(5)

    def _formatdb(self, database, type_, out_folder,
                  blastdb, database_type, log):
        err = open(os.path.join(out_folder, "log.txt"), "w")
        if database_type == "sRNA":
            change_format(database, "tmp_srna_database")
            os.remove(database)
            shutil.move("tmp_srna_database", database)
            log.write("Formating sRNA database.\n")
        else:
            log.write("Formating nr database.\n")
        db_file = ".".join(database.split(".")[:-1])
        self._run_format(blastdb, database, type_, db_file, err, log)
        err.close()
        if (database.endswith(".fa")) or (
                database.endswith(".fna")) or (
                database.endswith(".fasta")):
            database = ".".join(database.split(".")[:-1])
        return database

    def _merge_frag_tex_file(self, files, args_srna):
        '''merge the results of fragmented and tex treated libs'''
        if (args_srna.frag_wigs is not None) and (
                args_srna.tex_wigs is not None):
            self.helper.merge_file(files["frag_gff"], files["tex_gff"])
            self.helper.merge_file(files["frag_csv"], files["tex_csv"])
            shutil.move(files["tex_csv"], files["merge_csv"])
            self.helper.sort_gff(files["tex_gff"], files["merge_gff"])
            os.remove(files["frag_csv"])
            os.remove(files["frag_gff"])
            os.remove(files["tex_gff"])
        elif (args_srna.frag_wigs is not None):
            shutil.move(files["frag_csv"], files["merge_csv"])
            self.helper.sort_gff(files["frag_gff"], files["merge_gff"])
            os.remove(files["frag_gff"])
        elif (args_srna.tex_wigs is not None):
            shutil.move(files["tex_csv"], files["merge_csv"])
            self.helper.sort_gff(files["tex_gff"], files["merge_gff"])

    def _read_lib_wig(self, args_srna):
        libs, texs = read_libs(args_srna.input_libs, args_srna.wig_folder)
        wigs_f = read_wig(args_srna.wig_f_file, "+", libs)
        wigs_r = read_wig(args_srna.wig_r_file, "-", libs)
        return [libs, texs, wigs_f, wigs_r]

    def _run_normal(self, prefix, gff, tran, fuzzy_tss, args_srna, log):
        '''detection of intergenic and antisense sRNA'''
        tex_datas = None
        frag_datas = None
        if "tmp_cutoff_inter" in os.listdir(args_srna.out_folder):
            os.remove(os.path.join(args_srna.out_folder, "tmp_cutoff_inter"))
        files = {"frag_gff": None, "frag_csv": None,
                 "tex_gff": None, "tex_csv": None,
                 "merge_gff": None, "merge_csv": None}
        if self.tss_path is not None:
            if ("TSS_classes" in os.listdir(args_srna.out_folder)) and (
                not args_srna.source):
                tss = os.path.join(args_srna.out_folder,
                                   "TSS_classes", prefix + "_TSS.gff")
            else:
                tss = self.helper.get_correct_file(self.tss_path, "_TSS.gff",
                                                   prefix, None, None)
        else:
            tss = None
        if self.pro_path is not None:
            pro = self.helper.get_correct_file(
                    self.pro_path, "_processing.gff", prefix, None, None)
        else:
            pro = None
        if args_srna.frag_wigs is not None:
            files["frag_gff"] = os.path.join(
                    args_srna.out_folder, "_".join(["tmp_frag", prefix]))
            files["frag_csv"] = os.path.join(
                    args_srna.out_folder, "_".join(["tmp_frag_table", prefix]))
            args_srna = self.args_container.container_intersrna(
                             "frag", files, args_srna, prefix,
                             os.path.join(args_srna.gffs, gff), tran, tss,
                             pro, fuzzy_tss)
            frag_datas = self._read_lib_wig(args_srna)
            log.write("Running sRNA_intergenic.py to detecting intergenic "
                      "sRNA for {0} based on fragmented libs.\n".format(prefix))
            intergenic_srna(args_srna, frag_datas[0], frag_datas[1],
                            frag_datas[2], frag_datas[3], tss)
        if args_srna.tex_wigs is not None:
            files["tex_gff"] = os.path.join(
                    args_srna.out_folder, "_".join(["tmp_tex", prefix]))
            files["tex_csv"] = os.path.join(
                    args_srna.out_folder, "_".join(["tmp_tex_table", prefix]))
            args_srna = self.args_container.container_intersrna(
                           "tex", files, args_srna, prefix,
                           os.path.join(args_srna.gffs, gff), tran, tss,
                           pro, fuzzy_tss)
            tex_datas = self._read_lib_wig(args_srna)
            log.write("Running sRNA_intergenic.py to detecting intergenic "
                      "sRNA for {0} based on dRNA-Seq libs.\n".format(prefix))
            intergenic_srna(args_srna, tex_datas[0], tex_datas[1],
                            tex_datas[2], tex_datas[3], tss)
        files["merge_csv"] = "_".join([self.prefixs["normal_table"], prefix])
        files["merge_gff"] = "_".join([self.prefixs["normal"], prefix])
        self._merge_frag_tex_file(files, args_srna)
        return tss, frag_datas, tex_datas

    def _run_utrsrna(self, gff, tran, prefix, tss, pro, args_srna,
                     frag_datas, tex_datas, log):
        '''detection of UTR-derived sRNA'''
        if "tmp_median" in os.listdir(args_srna.out_folder):
            os.remove(os.path.join(args_srna.out_folder, "tmp_median"))
        files = {"frag_gff": None, "frag_csv": None,
                 "tex_gff": None, "tex_csv": None,
                 "merge_gff": None, "merge_csv": None}
        if args_srna.tex_wigs is not None:
            files["tex_gff"] = os.path.join(
                    args_srna.out_folder, "_".join(["tmp_utr_tex", prefix]))
            files["tex_csv"] = os.path.join(
                    args_srna.out_folder,
                    "_".join(["tmp_utr_tex_table", prefix]))
            args_srna = self.args_container.container_utrsrna(
                    os.path.join(args_srna.gffs, gff), tran, tss, files,
                    pro, os.path.join(self.fasta_path, prefix + ".fa"),
                    "tex", prefix, args_srna)
            log.write("Running sRNA_utr_derived.py to detect UTR-derived "
                      "sRNAs for {0} based on dRNA-Seq data.\n".format(prefix))
            utr_derived_srna(args_srna, tex_datas[0], tex_datas[1],
                             tex_datas[2], tex_datas[3])
        if args_srna.frag_wigs is not None:
            files["frag_gff"] = os.path.join(
                args_srna.out_folder, "_".join(["tmp_utr_frag", prefix]))
            files["frag_csv"] = os.path.join(
                args_srna.out_folder, "_".join(["tmp_utr_frag_table", prefix]))
            args_srna = self.args_container.container_utrsrna(
                    os.path.join(args_srna.gffs, gff), tran, tss, files,
                    pro, os.path.join(self.fasta_path, prefix + ".fa"),
                    "frag", prefix, args_srna)
            log.write("Running sRNA_utr_derived.py to detect UTR-derived "
                      "sRNAs for {0} based on fragmented libs.\n".format(prefix))
            utr_derived_srna(args_srna, frag_datas[0], frag_datas[1],
                             frag_datas[2], frag_datas[3])
        files["merge_csv"] = "_".join([self.prefixs["utr_table"], prefix])
        files["merge_gff"] = "_".join([self.prefixs["utr"], prefix])
        self._merge_frag_tex_file(files, args_srna)
        log.write("Running sRNA_filter_min_utr.py to filter out the "
                  "UTR-derived sRNAs which length is too short.\n")
        filter_utr(files["merge_gff"], files["merge_csv"], args_srna.min_utr)

    def _check_database(self, formatdb, database):
        if formatdb:
            if (database.endswith(".fa")) or (
                    database.endswith(".fna")) or (
                    database.endswith(".fasta")):
                return database
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
                            return database
        else:
            return database
        print("Error: The nr database or sRNA database is not in fasta "
              "format or the file name does not end with "
              ".fa or .fna or .fasta!")
        sys.exit()

    def _check_necessary_file(self, args_srna, log):
        if (args_srna.gffs is None) or (args_srna.trans is None) or (
                (args_srna.tex_wigs is None) and (
                args_srna.frag_wigs is None)):
            print("Error: Lack required files!")
            log.write("The annotation gff files, transcipt files, or wiggle "
                      "files do not be assigned.\n")
            sys.exit()
        if args_srna.utr_srna:
            if (args_srna.tss_folder is None):
                print("Error: Lack required TSS files for UTR "
                      "derived sRNA detection!")
                log.write("TSS files are required for detecting UTR-derived "
                          "sRNAs.\n")
                sys.exit()
            if (args_srna.pro_folder is None):
                print("Warning: Lack Processing site files for UTR "
                      "derived sRNA detection!")
                print("It may affect the results!")
        self._check_gff(args_srna.gffs)
        self._check_gff(args_srna.trans)
        args_srna.nr_database = self._check_database(args_srna.nr_format,
                                                args_srna.nr_database)
        args_srna.srna_database = self._check_database(args_srna.srna_format,
                                                  args_srna.srna_database)
        if args_srna.tss_folder is not None:
            self._check_gff(args_srna.tss_folder)
            self.multiparser.parser_gff(args_srna.tss_folder, "TSS")
            self.multiparser.combine_gff(args_srna.gffs, self.tss_path,
                                         None, "TSS")
        if args_srna.pro_folder is not None:
            self._check_gff(args_srna.pro_folder)
            self.multiparser.parser_gff(args_srna.pro_folder, "processing")
            self.multiparser.combine_gff(args_srna.gffs, self.pro_path,
                                         None, "processing")
        if args_srna.sorf_file is not None:
            self._check_gff(args_srna.sorf_file)
            self.multiparser.parser_gff(args_srna.sorf_file, "sORF")
            self.multiparser.combine_gff(args_srna.gffs, self.sorf_path,
                                         None, "sORF")
        if args_srna.import_info is not None:
            if args_srna.utr_srna or ("sec_str" in args_srna.import_info) or (
                   args_srna.nr_database is not None) or (
                   args_srna.srna_database is not None):
                if args_srna.fastas is None:
                    print("Error: Fasta file is not assinged!")
                    log.write("Fasta file is not assinged.\n")
                    sys.exit()
                self.multiparser.parser_fasta(args_srna.fastas)
                self.multiparser.combine_fasta(args_srna.gffs,
                                               self.fasta_path, None)
        if args_srna.terms is not None:
            self._check_gff(args_srna.terms)
            self.multiparser.parser_gff(args_srna.terms, "term")
            self.multiparser.combine_gff(args_srna.gffs, self.term_path,
                                         None, "term")
        else:
            self.term_path = None

    def _merge_tex_frag_datas(self, tex_datas, frag_datas):
        if (tex_datas is not None) and (frag_datas is not None):
            for index in [2, 3]:
                for strain, conds in frag_datas[index].items():
                    if strain not in tex_datas[index].keys():
                        tex_datas[index][strain] = conds
                    else:
                        for cond, tracks in conds.items():
                            tex_datas[index][strain][cond] = tracks
        elif (tex_datas is None) and (frag_datas is not None):
            tex_datas = frag_datas
        return tex_datas

    def _run_program(self, args_srna, log):
        prefixs = []
        tss = None
        for gff in os.listdir(args_srna.gffs):
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                prefixs.append(prefix)
                tran = self.helper.get_correct_file(
                        self.tran_path, "_transcript.gff", prefix, None, None)
                gffs = {"merge": "_".join([self.prefixs["merge"], prefix]),
                        "utr": "_".join([self.prefixs["utr"], prefix]),
                        "normal": "_".join([self.prefixs["normal"], prefix])}
                csvs = {"merge": "_".join([
                            self.prefixs["merge_table"], prefix]),
                        "utr": "_".join([self.prefixs["utr_table"], prefix]),
                        "normal": "_".join([
                            self.prefixs["normal_table"], prefix])}
                if not args_srna.source:
                    if "TSS_classes" not in os.listdir(args_srna.out_folder):
                        os.mkdir(os.path.join(args_srna.out_folder, "TSS_classes"))
                    print("Classifying TSSs of {0}".format(prefix))
                    upstream(os.path.join(self.tss_path, prefix + "_TSS.gff"),
                             None,
                             os.path.join(args_srna.gffs, prefix + ".gff"),
                             os.path.join(args_srna.out_folder, "TSS_classes",
                             "_".join([prefix, "TSS.gff"])), args_srna, prefix)
                print("Running sRNA detection of {0}".format(prefix))
                tss, frag_datas, tex_datas = self._run_normal(
                        prefix, gff, tran, args_srna.fuzzy_tsss["inter"],
                        args_srna, log)
                if args_srna.utr_srna:
                    print("Running UTR derived sRNA detection of {0}".format(
                          prefix))
                    if tss is None:
                        tss = self.helper.get_correct_file(
                                self.tss_path, "_TSS.gff", prefix, None, None)
                    if self.pro_path is not None:
                        pro = self.helper.get_correct_file(
                                self.pro_path, "_processing.gff",
                                prefix, None, None)
                    else:
                        pro = None
                    if tss is not None:
                        self._run_utrsrna(gff, tran, prefix, tss, pro,
                                          args_srna, frag_datas, tex_datas, log)
                tex_datas = self._merge_tex_frag_datas(tex_datas, frag_datas)
                del frag_datas
                gc.collect()
                self._merge_srna(args_srna, gffs, csvs, prefix,
                                 os.path.join(args_srna.gffs, gff), tss, tex_datas)
                del tex_datas
                filter_frag(csvs["merge"], gffs["merge"])
                self.helper.sort_gff(gffs["merge"],
                                     "_".join([self.prefixs["basic"], prefix]))
                log.write("\t" + "_".join([self.prefixs["basic"], prefix]) + 
                          " is generated to temporary store sRNA candidates.\n")
                log.write("\t" + csvs["merge"] + " is generated to temporary store "
                          "the detail information of sRNA candidates.\n")
        return prefixs

    def _merge_srna(self, args_srna, gffs, csvs, prefix,
                    gff_file, tss, tex_datas):
        print("Merging data of sRNA")
        merge_srna_gff(gffs, args_srna.in_cds,
                       args_srna.cutoff_overlap, gff_file, args_srna.ex_srna)
        merge_srna_table(gffs["merge"], csvs, tex_datas[2], tex_datas[3],
                         tss, args_srna)

    def _run_RNAfold(self, seq_file, rnafold, sec_file, log):
        log.write("Running RNAfold to predict secondary structure.\n")
        log.write("Please make sure the version of Vienna Packge is at "
                  "least 2.3.2.\n")
        log.write(" ".join(["cat", seq_file, "|",
                  rnafold, "-p", ">", sec_file]) + "\n")
        os.system(" ".join(["cat", seq_file, "|",
                  rnafold, "-p", ">", sec_file]))
        log.write("Done!\n")
        log.write("\t" + sec_file + " is generated.\n")

    def _get_seq_sec(self, fasta_path, out_folder, prefix, sec_path,
                     dot_path, rnafold, log):
        '''extract the sec str energy'''
        detect = False
        for fasta in os.listdir(fasta_path):
            if fasta.endswith(".fa") and (
               ".".join(fasta.split(".")[:-1]) == prefix):
                detect = True
                break
        if detect:
            detect = False
            seq_file = os.path.join(out_folder, "_".join(["sRNA_seq", prefix]))
            sec_file = os.path.join(out_folder, "_".join(["sRNA_2d", prefix]))
            index_file = os.path.join(out_folder, "_".join(
                ["sRNA_index", prefix]))
            log.write("Running helper.py to get the sequences of sRNA for "
                      "{0}.\n".format(prefix))
            self.helper.get_seq("_".join([self.prefixs["basic"], prefix]),
                                os.path.join(fasta_path, fasta), index_file)
            modify_header(seq_file, index_file)
            log.write("\t" + seq_file + " is generated.\n")
        else:
            print("Error: There is not fasta file of {0}!".format(prefix))
            print("Please check your imported information.")
            log.write("No fasta file of {0}.\n".format(prefix))
            sys.exit()
        tmp_path = os.path.join(out_folder, "tmp_srna")
        self.helper.check_make_folder(tmp_path)
        main_path = os.getcwd()
        os.chdir(tmp_path)
        sec_file = os.path.join(main_path, sec_file)
        seq_file = os.path.join(main_path, seq_file)
        index_file = os.path.join(main_path, index_file)
        tmp_sec_path = os.path.join(main_path, sec_path)
        tmp_dot_path = os.path.join(main_path, dot_path)
        self._run_RNAfold(seq_file, rnafold, sec_file, log)
        extract_info_sec(sec_file, seq_file, index_file)
        os.remove(index_file)
        log.write("Running extract_sRNA_info.py to extract the energy "
                  "information for {0}.\n".format(prefix))
        extract_energy(os.path.join(main_path,
                       "_".join([self.prefixs["basic"], prefix])),
                       sec_file, os.path.join(main_path,
                       "_".join([self.prefixs["energy"], prefix])))
        log.write("\t" + os.path.join(main_path, "_".join([
            self.prefixs["energy"], prefix])) + " is generated to temporary "
            "store energy information.\n")
        for ps in os.listdir(os.getcwd()):
            new_ps = ps.replace("|", "_")
            shutil.move(ps, new_ps)
        return {"sec": tmp_sec_path, "dot": tmp_dot_path, "main": main_path,
                "tmp": os.path.join(main_path, tmp_path)}

    def _run_replot(self, relplot_pl, tmp_paths, file_, dot_file, rel_file, log):
        log.write(" ".join([relplot_pl,
                  os.path.join(tmp_paths["tmp"], file_),
                  os.path.join(tmp_paths["tmp"], dot_file),
                  ">", os.path.join(tmp_paths["tmp"], rel_file)]) + "\n")
        os.system(" ".join([relplot_pl,
                  os.path.join(tmp_paths["tmp"], file_),
                  os.path.join(tmp_paths["tmp"], dot_file),
                  ">", os.path.join(tmp_paths["tmp"], rel_file)]))

    def _replot_sec(self, relplot_pl, tmp_paths, prefix, log):
        log.write("Running relplot.pl for {0}.\n".format(prefix))
        for file_ in os.listdir(os.getcwd()):
            if file_.endswith("ss.ps"):
                dot_file = file_.replace("ss.ps", "dp.ps")
                rel_file = file_.replace("ss.ps", "rss.ps")
                print("Relplotting {0}".format(file_))
                self._run_replot(relplot_pl, tmp_paths, file_,
                                 dot_file, rel_file, log)
        log.write("Done!\n")
        os.mkdir(os.path.join(tmp_paths["sec"], prefix))
        os.mkdir(os.path.join(tmp_paths["dot"], prefix))
        self.helper.move_all_content(
                tmp_paths["tmp"], os.path.join(tmp_paths["sec"], prefix),
                ["rss.ps"])
        self.helper.move_all_content(
                tmp_paths["tmp"], os.path.join(tmp_paths["dot"], prefix),
                ["dp.ps"])
        log.write("All plots are stored in {0} and {1}.\n".format(
                  os.path.join(tmp_paths["sec"], prefix),
                  os.path.join(tmp_paths["dot"], prefix)) + "\n")

    def _run_mountain(self, mountain_pl, dot_path, dot_file, out, log):
        log.write(" ".join([mountain_pl,
                  os.path.join(dot_path, dot_file)]) + "\n")
        call([mountain_pl,
              os.path.join(dot_path, dot_file)], stdout=out)

    def _plot_mountain(self, mountain, moun_path,
                       tmp_paths, prefix, mountain_pl, log):
        if mountain:
            tmp_moun_path = os.path.join(tmp_paths["main"], moun_path)
            os.mkdir(os.path.join(tmp_moun_path, prefix))
            txt_path = os.path.join(tmp_paths["tmp"], "tmp_txt")
            self.helper.check_make_folder(txt_path)
            print("Generating mountain plots of {0}".format(prefix))
            dot_path = os.path.join(tmp_paths["dot"], prefix)
            log.write("Running mountain.pl for {0}.\n".format(prefix))
            for dot_file in os.listdir(dot_path):
                if dot_file.endswith("dp.ps"):
                    moun_txt = os.path.join(tmp_paths["tmp"], "mountain.txt")
                    out = open(moun_txt, "w")
                    moun_file = dot_file.replace("dp.ps", "mountain.pdf")
                    print("Generating {0}".format(moun_file))
                    self._run_mountain(mountain_pl, dot_path, dot_file, out, log)
                    plot_mountain_plot(moun_txt, moun_file)
                    shutil.move(moun_file,
                                os.path.join(tmp_moun_path, prefix, moun_file))
                    out.close()
                    os.remove(moun_txt)
            log.write("Done!\n")
            log.write("All plots are stored in {0}.".format(
                      os.path.join(tmp_moun_path, prefix)))

    def _compute_2d_and_energy(self, args_srna, prefixs, log):
        print("Running energy calculation")
        moun_path = os.path.join(args_srna.out_folder, "figs",
                                 "mountain_plots")
        sec_path = os.path.join(args_srna.out_folder, "figs",
                                "sec_plots")
        dot_path = os.path.join(args_srna.out_folder, "figs",
                                "dot_plots")
        self.helper.remove_all_content(sec_path, None, "dir")
        self.helper.remove_all_content(dot_path, None, "dir")
        self.helper.remove_all_content(moun_path, None, "dir")
        for prefix in prefixs:
            tmp_paths = self._get_seq_sec(
                    self.fasta_path, args_srna.out_folder, prefix, sec_path,
                    dot_path, args_srna.rnafold, log)
            self._replot_sec(args_srna.relplot_pl, tmp_paths, prefix, log)
            self._plot_mountain(args_srna.mountain, moun_path, tmp_paths,
                                prefix, args_srna.mountain_pl, log)
            os.chdir(tmp_paths["main"])
            shutil.move("_".join([self.prefixs["energy"], prefix]),
                        "_".join([self.prefixs["basic"], prefix]))
            log.write("_".join([self.prefixs["basic"], prefix]) + " is updated, "
                      "and " + "_".join([self.prefixs["energy"], prefix]) + 
                      " is deleted.\n")
            shutil.rmtree(os.path.join(args_srna.out_folder, "tmp_srna"))

    def _run_blast(self, program, database, e, seq_file,
                   blast_file, strand, para_num, processes, log):
        if para_num == 1:
            log.write(" ".join([program, "-db", database,
                  "-evalue", str(e), "-strand", strand, "-query", seq_file,
                  "-out", blast_file]) + "\n")
            call([program, "-db", database,
                  "-evalue", str(e), "-strand", strand, "-query", seq_file,
                  "-out", blast_file])
        else:
            log.write(" ".join([program, "-db", database,
                       "-evalue", str(e), "-strand", strand, "-query", seq_file,
                       "-out", blast_file]) + "\n")
            p = Popen([program, "-db", database,
                       "-evalue", str(e), "-strand", strand, "-query", seq_file,
                       "-out", blast_file])
            processes.append(p)

    def _run_para_blast(self, program, database, e, seq_file,
                        blast_file, strand, paras, log):
        srnas = {}
        with open(seq_file) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    name = line
                    srnas[name] = ""
                else:
                    srnas[name] = line
        file_num = int(len(srnas) / paras)
        processes = []
        if (file_num == 0) or (paras == 1):
            self._run_blast(program, database, e, seq_file,
                            blast_file, strand, 1, processes, log)
        else:
            cur_para = 0
            line_num = 0
            first = True
            seq_files = []
            log.write("{0} is splited to {1} subset files.\n".format(
                seq_file, paras))
            for name, seq in srnas.items():
                if (line_num >= file_num) or first:
                    if (not first) and (cur_para < paras):
                        out.close()
                    first = False
                    if cur_para < paras:
                        out = open("_".join([seq_file, str(cur_para)]),
                                   "w")
                        seq_files.append("_".join([seq_file, str(cur_para)]))
                        line_num = 0
                        cur_para += 1
                if line_num < file_num:
                    out.write(name + "\n")
                    out.write(seq + "\n")
                if (cur_para == paras) and (line_num >= file_num):
                    out.write(name + "\n")
                    out.write(seq + "\n")
                line_num += 1
            out.close()
            for para in range(paras):
                self._run_blast(
                        program, database, e, "_".join([seq_file, str(para)]),
                        "_".join([blast_file, strand, str(para)]), strand, paras,
                        processes, log)
            self._wait_process(processes)
            for para in range(paras):
                cur_blast_file = "_".join([blast_file, strand, str(para)])
                self.helper.merge_file(cur_blast_file, blast_file)
                os.remove(cur_blast_file)
            for file_ in seq_files:
                os.remove(file_)
        log.write("Done!\n")
        if (os.path.exists(blast_file)):
            log.write("\t" + blast_file + " is generated.\n")

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

    def _blast(self, database, database_format, data_type, args_srna,
               prefixs, program, database_type, e, filters, log):
        if (database is None):
            log.write(" No database was assigned!\n")
            print("Error: No database was assigned!")
        else:
            if database_format:
                database = self._formatdb(database, data_type,
                                          args_srna.out_folder,
                                          args_srna.blastdb, database_type, log)
            for prefix in prefixs:
                blast_file = os.path.join(
                        args_srna.out_folder, "blast_results_and_misc",
                        "_".join([database_type, "blast", prefix + ".txt"]))
                if os.path.exists(blast_file):
                    os.remove(blast_file)
                srna_file = "_".join([self.prefixs["basic"], prefix])
                out_file = os.path.join(
                        args_srna.out_folder,
                        "_".join(["tmp", database_type, prefix]))
                print("Running Blast of {0} in {1}".format(prefix, database))
                seq_file = os.path.join(
                        args_srna.out_folder, "_".join(["sRNA_seq", prefix]))
                if (seq_file not in os.listdir(args_srna.out_folder)) or ((
                        database_type == "nr") and ("sec_str" in filters)):
                    log.write("Running helper.py to extract the sequences "
                              "of sRNAs.\n")
                    self.helper.get_seq(
                            srna_file,
                            os.path.join(self.fasta_path, prefix + ".fa"),
                            seq_file)
                    log.write("\t" + seq_file + " is generated.\n")
                if database_type == "nr":
                    log.write("Running BLAST+ for nr database for {0}.".format(
                        prefix))
                    log.write("Make sure the version of BLAST+ is at least 2.2.28+.\n")
                    tmp_plus, tmp_minus = self._get_strand_fasta(
                            seq_file, args_srna.out_folder)
                    tmp_blast = os.path.join(args_srna.out_folder,
                                             "blast_results_and_misc",
                                             "tmp_blast.txt")
                    if os.path.exists(tmp_blast):
                        os.remove(tmp_blast)
                    self._run_para_blast(program, database, e,
                                         tmp_plus, tmp_blast, "plus",
                                         args_srna.para_blast, log)
                    self._run_para_blast(program, database, e,
                                         tmp_minus, blast_file, "minus",
                                         args_srna.para_blast, log)
                    self.helper.merge_file(tmp_blast, blast_file)
                    os.remove(tmp_plus)
                    os.remove(tmp_minus)
                else:
                    log.write("Running BLAST+ for sRNA database for {0}.".format(
                        prefix))
                    log.write("Make sure the version of BLAST+ is at least 2.2.28+.\n")
                    self._run_para_blast(program, database, e,
                                         seq_file, blast_file, "both",
                                         args_srna.para_blast, log)
                log.write("Running extract_sRNA_info.py to extract BLAST "
                          "information.\n")
                extract_blast(blast_file, srna_file, out_file,
                              out_file + ".csv", database_type,
                              args_srna.blast_score_s,
                              args_srna.blast_score_n)
                log.write(srna_file + " is updated.\n")
                shutil.move(out_file, srna_file)

    def _class_srna(self, prefixs, args_srna, log):
        '''classify the sRNA based on the filters'''
        if (args_srna.import_info is not None) or (
                args_srna.srna_database is not None) or (
                args_srna.nr_database is not None) or (
                self.sorf_path is not None) or (
                self.tss_path is not None) or (
                self.term_path is not None) or (
                args_srna.promoter_table is not None):
            log.write("Running sRNA_class.py to classify sRNAs based on "
                      "input files and --filter_info.\n")
            log.write("The following files are generated:\n")
            for prefix in prefixs:
                print("Classifying sRNA of {0}".format(prefix))
                class_gff = os.path.join(self.gff_output, "for_classes")
                class_table = os.path.join(self.table_output, "for_classes")
                self.helper.check_make_folder(os.path.join(class_table,
                                                           prefix))
                self.helper.check_make_folder(os.path.join(class_gff, prefix))
                class_gff = os.path.join(class_gff, prefix)
                class_table = os.path.join(class_table, prefix)
                self.helper.check_make_folder(class_table)
                self.helper.check_make_folder(class_gff)
                out_stat = os.path.join(
                        self.stat_path, "_".join([
                            "stat_sRNA_class", prefix + ".csv"]))
                classify_srna(os.path.join(self.all_best["all_gff"],
                              "_".join([prefix, "sRNA.gff"])), class_gff,
                              out_stat, args_srna)
                log.write("\t" + out_stat + "\n")
                for srna in os.listdir(class_gff):
                    out_table = os.path.join(
                            class_table, srna.replace(".gff", ".csv"))
                    gen_srna_table(
                        os.path.join(class_gff, srna),
                        "_".join([self.prefixs["merge_table"], prefix]),
                        "_".join([self.tmps["nr"], prefix + ".csv"]),
                        "_".join([self.tmps["srna"], prefix + ".csv"]),
                        args_srna, out_table, self.term_path)
                for folder in (class_gff, class_table):
                    for file_ in os.listdir(folder):
                        log.write("\t" + os.path.join(folder, file_) + "\n")

    def _get_best_result(self, prefixs, args_srna, log):
        '''get the best results based on the filters'''
        log.write("Running gen_srna_output to select the best candidates.\n")
        log.write("The following files are generated:\n")
        for prefix in prefixs:
            best_gff = os.path.join(self.all_best["best_gff"],
                                    "_".join([prefix, "sRNA.gff"]))
            best_table = os.path.join(self.all_best["best_table"],
                                      "_".join([prefix, "sRNA.csv"]))
            gen_best_srna(os.path.join(self.all_best["all_gff"],
                                       "_".join([prefix, "sRNA.gff"])),
                          best_gff, args_srna)
            gen_srna_table(os.path.join(self.all_best["best_gff"],
                           "_".join([prefix, "sRNA.gff"])),
                           "_".join([self.prefixs["merge_table"], prefix]),
                           "_".join([self.tmps["nr"], prefix + ".csv"]),
                           "_".join([self.tmps["srna"], prefix + ".csv"]),
                           args_srna, best_table, self.term_path)
            log.write("\t" + best_gff + "\n")
            log.write("\t" + best_table + "\n")

    def _remove_file(self, args_srna):
        self.helper.remove_all_content(args_srna.out_folder, "tmp_", "dir")
        self.helper.remove_all_content(args_srna.out_folder, "tmp_", "file")
        self.helper.remove_tmp_dir(args_srna.fastas)
        self.helper.remove_tmp_dir(args_srna.gffs)
        self.helper.remove_tmp(self.gff_output)
        if "temp_wig" in os.listdir(args_srna.out_folder):
            shutil.rmtree(os.path.join(args_srna.out_folder, "temp_wig"))
        if (args_srna.frag_wigs is not None) and (
                args_srna.tex_wigs is not None):
            shutil.rmtree(args_srna.merge_wigs)
        self.helper.remove_tmp_dir(args_srna.trans)
        if args_srna.tss_folder is not None:
            self.helper.remove_tmp_dir(args_srna.tss_folder)
        if args_srna.pro_folder is not None:
            self.helper.remove_tmp_dir(args_srna.pro_folder)
        if args_srna.sorf_file is not None:
            self.helper.remove_tmp_dir(args_srna.sorf_file)
        if "tmp_median" in os.listdir(args_srna.out_folder):
            os.remove(os.path.join(args_srna.out_folder, "tmp_median"))
        if self.term_path is not None:
            self.helper.remove_tmp_dir(args_srna.terms)
        tmp_blast = os.path.join(args_srna.out_folder,
                                 "blast_results_and_misc",
                                 "tmp_blast.txt")
        if os.path.exists(tmp_blast):
            os.remove(tmp_blast)

    def _filter_srna(self, args_srna, prefixs, log):
        '''set the filter of sRNA'''
        if args_srna.compute_sec_str:
            self._compute_2d_and_energy(args_srna, prefixs, log)
        if args_srna.nr_database is not None:
            self._blast(args_srna.nr_database, args_srna.nr_format, "prot",
                        args_srna, prefixs, args_srna.blastx, "nr",
                        args_srna.e_nr, args_srna.import_info, log)
        if self.sorf_path is not None:
            for prefix in prefixs:
                if ("_".join([prefix, "sORF.gff"]) in
                        os.listdir(self.sorf_path)):
                    tmp_srna = os.path.join(args_srna.out_folder,
                                            "".join(["tmp_srna_sorf", prefix]))
                    tmp_sorf = os.path.join(args_srna.out_folder,
                                            "".join(["tmp_sorf_srna", prefix]))
                    log.write("Running compare_sRNA_sORF.py to compare sRNAs "
                              "and sORFs.\n")
                    srna_sorf_comparison(
                            "_".join([self.prefixs["basic"], prefix]),
                            os.path.join(self.sorf_path,
                                         "_".join([prefix, "sORF.gff"])),
                            tmp_srna, tmp_sorf)
                    os.remove(tmp_sorf)
                    shutil.move(tmp_srna,
                                "_".join([self.prefixs["basic"], prefix]))
                    log.write("_".join([self.prefixs["basic"], prefix]) + 
                              " is updated.\n")
        if args_srna.srna_database is not None:
            self._blast(args_srna.srna_database, args_srna.srna_format, "nucl",
                        args_srna, prefixs, args_srna.blastn, "sRNA",
                        args_srna.e_srna, args_srna.import_info, log)

    def _import_info_format(self, import_info):
        new_info = []
        for info in import_info:
            info = info.lower()
            new_info.append(info)
        return new_info

    def _gen_table(self, prefixs, args_srna, log):
        log.write("Running gen_srna_output.py to generate sRNA table.\n")
        log.write("The following files are generated.\n")
        for prefix in prefixs:
            print("Generating table for " + prefix)
            out_table = os.path.join(self.all_best["all_table"],
                                     "_".join([prefix, "sRNA.csv"]))
            gen_srna_table(os.path.join(self.all_best["all_gff"],
                           "_".join([prefix, "sRNA.gff"])),
                           "_".join([self.prefixs["merge_table"], prefix]),
                           "_".join([self.tmps["nr"], prefix + ".csv"]),
                           "_".join([self.tmps["srna"], prefix + ".csv"]),
                           args_srna, out_table, self.term_path)
            log.write("\t" + out_table + "\n")

    def _print_rank_all(self, prefixs, log):
        log.write("Running print_rank_all.py for ranking the sRNA candidates.\n")
        log.write("The following files are updated:\n")
        for prefix in prefixs:
            all_table = os.path.join(self.all_best["all_table"],
                                     "_".join([prefix, "sRNA.csv"]))
            best_table = os.path.join(self.all_best["best_table"],
                                      "_".join([prefix, "sRNA.csv"]))
            print_rank_all(all_table, best_table)
            log.write("\t" + all_table + "\n")
            log.write("\t" + best_table + "\n")

    def _filter_min_utr(self, prefixs, min_utr):
        '''filter out the low expressed UTR-derived sRNA'''
        for prefix in prefixs:
            filter_utr(os.path.join(self.all_best["all_gff"],
                                    "_".join([prefix, "sRNA.gff"])),
                       os.path.join(self.all_best["all_table"],
                                    "_".join([prefix, "sRNA.csv"])), min_utr)

    def _antisense(self, gffs, prefixs):
        '''detection of antisense'''
        for prefix in prefixs:
            all_table = os.path.join(self.all_best["all_table"],
                                     "_".join([prefix, "sRNA.csv"]))
            best_table = os.path.join(self.all_best["best_table"],
                                      "_".join([prefix, "sRNA.csv"]))
            all_gff = os.path.join(self.all_best["all_gff"],
                                   "_".join([prefix, "sRNA.gff"]))
            best_gff = os.path.join(self.all_best["best_gff"],
                                    "_".join([prefix, "sRNA.gff"]))
            srna_antisense(all_gff, all_table,
                           os.path.join(gffs, prefix + ".gff"))
            srna_antisense(best_gff, best_table,
                           os.path.join(gffs, prefix + ".gff"))

    def _blast_stat(self, stat_path, srna_tables, log):
        '''do statistics for blast result'''
        log.write("Running blast_class.py to do statistics for BLAST results.\n")
        for srna_table in os.listdir(os.path.join(srna_tables,
                                                  "best_candidates")):
            out_srna_blast = os.path.join(
                    stat_path, "stat_" +
                    srna_table.replace(".csv", "_blast.csv"))
            blast_class(os.path.join(srna_tables, "best_candidates", srna_table),
                        out_srna_blast)
            log.write("\t" + out_srna_blast + " is generated.\n")

    def _compare_term_promoter(self, out_table, prefix, args_srna, log):
        '''compare sRNA with terminator and promoter'''
        if self.term_path is not None:
            log.write("Running compare_srna_term.py to compare sRNAs with "
                      "terminators.\n")
            compare_srna_term(os.path.join(self.all_best["all_gff"],
                              "_".join([prefix, "sRNA.gff"])),
                              out_table, os.path.join(self.term_path,
                              "_".join([prefix, "term.gff"])),
                              args_srna.fuzzy_b, args_srna.fuzzy_a)
            log.write(os.path.join(self.all_best["all_gff"],
                      "_".join([prefix, "sRNA.gff"])) + " is updated.\n")
            log.write(out_table + " is updated.\n")
        if (args_srna.promoter_table is not None):
            log.write("Running compare_srna_term.py to compare sRNAs with "
                      "promoters.\n")
            compare_srna_promoter(os.path.join(self.all_best["all_gff"],
                                  "_".join([prefix, "sRNA.gff"])),
                                  out_table, args_srna)
            log.write(os.path.join(self.all_best["all_gff"],
                       "_".join([prefix, "sRNA.gff"])) + " is updated.\n")
            log.write(out_table + " is updated.\n")

    def _get_poly_u(self, prefixs, args_srna, log):
        print("Searching poly U tail ...")
        log.write("Running get_srna_poly_u.py to seach the poly U "
                  "tails of sRNAs.\n")
        for prefix in prefixs:
            get_srna_poly_u("_".join([self.prefixs["basic"], prefix]),
                            os.path.join(self.fasta_path, prefix + ".fa"),
                            "_".join([self.prefixs["merge_table"], prefix]),
                            args_srna)

    def _re_table(self, args_srna, prefixs, log):
        log.write("Running re_table.py to generate coverage information.\n")
        log.write("The following files are updated:\n")
        for type_ in ["all_candidates", "best_candidates"]:
            for prefix in prefixs:
                srna_table = os.path.join(args_srna.out_folder, "tables",
                                          type_, "_".join([
                                               prefix, "sRNA.csv"]))
                reorganize_table(args_srna.libs, args_srna.merge_wigs,
                                 "Track/Coverage", srna_table)
                log.write("\t" + srna_table + "\n")
        for c_table in os.listdir(os.path.join(args_srna.out_folder, "tables",
                                               "for_classes", prefix)):
            for prefix in prefixs:
                srna_table_c = os.path.join(args_srna.out_folder, "tables",
                                            "for_classes", prefix, c_table)
                reorganize_table(args_srna.libs, args_srna.merge_wigs,
                                 "Track/Coverage", srna_table_c)
                log.write("\t" + srna_table_c + "\n")

    def _check_overlap_cds(self, args_srna, prefixs, log):
        log.write("Running check_srna_overlap.py to compare sRNAs with "
                  "genome annotations.\n")
        log.write("The following files are updated:\n")
        for type_ in ["all_candidates", "best_candidates"]:
            for prefix in prefixs:
                srna_table = os.path.join(args_srna.out_folder, "tables",
                                          type_, "_".join([
                                                prefix, "sRNA.csv"]))
                gff_file = os.path.join(args_srna.gffs, prefix + ".gff")
                check_overlap(srna_table, gff_file)
                log.write("\t" + srna_table + "\n")
        for c_table in os.listdir(os.path.join(args_srna.out_folder, "tables",
                                               "for_classes", prefix)):
            for prefix in prefixs:
                gff_file = os.path.join(args_srna.gffs, prefix + ".gff")
                srna_table_c = os.path.join(args_srna.out_folder, "tables",
                                           "for_classes", prefix, c_table)
                check_overlap(srna_table_c, gff_file)
                log.write("\t" + srna_table_c + "\n")

    def run_srna_detection(self, args_srna, log):
        self._check_necessary_file(args_srna, log)
        self.multiparser.parser_gff(args_srna.trans, "transcript")
        self.multiparser.combine_gff(args_srna.gffs, self.tran_path,
                                     None, "transcript")
        if args_srna.import_info is not None:
            args_srna.import_info = self._import_info_format(args_srna.import_info)
        prefixs = self._run_program(args_srna, log)
        self._get_poly_u(prefixs, args_srna, log)
        self._filter_srna(args_srna, prefixs, log)
        for prefix in prefixs:
            shutil.copyfile("_".join([self.prefixs["basic"], prefix]),
                            os.path.join(self.all_best["all_gff"],
                            "_".join([prefix, "sRNA.gff"])))
            log.write("\t" + os.path.join(self.all_best["all_gff"],
                      "_".join([prefix, "sRNA.gff"])) + " is generated, and "
                      "_".join([self.prefixs["basic"], prefix]) + " is deleted.\n")
            self._compare_term_promoter("_".join([self.prefixs["merge_table"],
                                        prefix]), prefix, args_srna, log)
        self._gen_table(prefixs, args_srna, log)
        self._class_srna(prefixs, args_srna, log)
        self._get_best_result(prefixs, args_srna, log)
        self._print_rank_all(prefixs, log)
        if args_srna.srna_database is not None:
            if "blast_srna" in args_srna.import_info:
                self._blast_stat(self.stat_path, self.table_output, log)
        self._check_overlap_cds(args_srna, prefixs, log)
        self._re_table(args_srna, prefixs, log)
        self._remove_file(args_srna)
