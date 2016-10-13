import os, gc
import sys
import shutil
from subprocess import call
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.sRNA_intergenic import intergenic_srna
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
        self.fasta_path = os.path.join(args_srna.fastas, "tmp")
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
        self.best_table = os.path.join(self.table_output, "best")
        self.table_output = os.path.join(args_srna.out_folder, "tables")
        self.stat_path = os.path.join(args_srna.out_folder, "statistics")
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

    def _formatdb(self, database, type_, out_folder,
                  blast_path, database_type):
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

    def _run_normal(self, prefix, gff, tran, fuzzy_tss, args_srna):
        '''detection of intergenic and antisense sRNA'''
        tex_datas = None
        frag_datas = None
        if "tmp_cutoff_inter" in os.listdir(args_srna.out_folder):
            os.remove(os.path.join(args_srna.out_folder, "tmp_cutoff_inter"))
        files = {"frag_gff": None, "frag_csv": None,
                 "tex_gff": None, "tex_csv": None,
                 "merge_gff": None, "merge_csv": None}
        if self.tss_path is not None:
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
            intergenic_srna(args_srna, frag_datas[0], frag_datas[1],
                            frag_datas[2], frag_datas[3])
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
            intergenic_srna(args_srna, tex_datas[0], tex_datas[1],
                            tex_datas[2], tex_datas[3])
        files["merge_csv"] = "_".join([self.prefixs["normal_table"], prefix])
        files["merge_gff"] = "_".join([self.prefixs["normal"], prefix])
        self._merge_frag_tex_file(files, args_srna)
        if ("TSS_class" in os.listdir(args_srna.out_folder)) and (
                not args_srna.tss_source):
            tss = os.path.join(args_srna.out_folder,
                               "TSS_class", prefix + "_TSS.gff")
        return tss, frag_datas, tex_datas

    def _run_utrsrna(self, gff, tran, prefix, tss, pro, args_srna,
                     frag_datas, tex_datas):
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
            utr_derived_srna(args_srna, frag_datas[0], frag_datas[1],
                             frag_datas[2], frag_datas[3])
        files["merge_csv"] = "_".join([self.prefixs["utr_table"], prefix])
        files["merge_gff"] = "_".join([self.prefixs["utr"], prefix])
        self._merge_frag_tex_file(files, args_srna)
        filter_utr(files["merge_gff"], files["merge_csv"], args_srna.min_utr)

    def _check_necessary_file(self, args_srna):
        if (args_srna.gffs is None) or (args_srna.trans is None) or (
                (args_srna.tex_wigs is None) and (
                args_srna.frag_wigs is None)):
            print("Error: lack required files!!!!")
            sys.exit()
        if args_srna.utr_srna:
            if (args_srna.tss_folder is None):
                print("Error: lack required TSS files for UTR "
                      "derived sRNA detection!!!!")
                sys.exit()
            if (args_srna.pro_folder is None):
                print("Warning: lack Processing site files for UTR "
                      "derived sRNA detection!!!")
                print("it may effect the results!!!!")
        self._check_gff(args_srna.gffs)
        self._check_gff(args_srna.trans)
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
                    print("Error: lack required fasta files for UTR "
                          "derived sRNA detection!!!!")
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

    def _run_program(self, args_srna):
        prefixs = []
        tss = None
        for gff in os.listdir(args_srna.gffs):
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                prefixs.append(prefix)
                print("Running sRNA detection of {0}....".format(prefix))
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
                tss, frag_datas, tex_datas = self._run_normal(
                        prefix, gff, tran, args_srna.fuzzy_tsss["inter"],
                        args_srna)
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
                                          args_srna, frag_datas, tex_datas)
                tex_datas = self._merge_tex_frag_datas(tex_datas, frag_datas)
                del frag_datas
                gc.collect()
                self._merge_srna(args_srna, gffs, csvs, prefix,
                                 os.path.join(args_srna.gffs, gff), tss, tex_datas)
                del tex_datas
                filter_frag(csvs["merge"], gffs["merge"])
                self.helper.sort_gff(gffs["merge"],
                                     "_".join([self.prefixs["basic"], prefix]))
        return prefixs

    def _merge_srna(self, args_srna, gffs, csvs, prefix,
                    gff_file, tss, tex_datas):
        print("merging data of sRNA...")
        merge_srna_gff(gffs, args_srna.in_cds,
                       args_srna.cutoff_overlap, gff_file)
        merge_srna_table(gffs["merge"], csvs, tex_datas[2], tex_datas[3],
                         tss, args_srna)

    def _run_RNAfold(self, seq_file, vienna_path, sec_file):
        os.system(" ".join(["cat", seq_file, "|",
                  os.path.join(vienna_path, "RNAfold"),
                  "-p", ">", sec_file]))

    def _get_seq_sec(self, fasta_path, out_folder, prefix, sec_path,
                     dot_path, vienna_path):
        '''extract the sec str energy'''
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
            print("please check your imported information")
            sys.exit()
        tmp_path = os.path.join(out_folder, "tmp_srna")
        self.helper.check_make_folder(tmp_path)
        main_path = os.getcwd()
        os.chdir(tmp_path)
        sec_file = os.path.join(main_path, sec_file)
        seq_file = os.path.join(main_path, seq_file)
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

    def _replot_sec_to_pdf(self, vienna_util, tmp_paths,
                           ps2pdf14_path, prefix):
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
        self.helper.move_all_content(
                tmp_paths["tmp"], os.path.join(tmp_paths["sec"], prefix),
                ["rss.pdf"])
        self.helper.move_all_content(
                tmp_paths["tmp"], os.path.join(tmp_paths["dot"], prefix),
                ["dp.pdf"])

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

    def _compute_2d_and_energy(self, args_srna, prefixs):
        print("Running energy calculation....")
        moun_path = os.path.join(args_srna.out_folder, "mountain_plot")
        sec_path = os.path.join(args_srna.out_folder, "sec_structure",
                                "sec_plot")
        dot_path = os.path.join(args_srna.out_folder, "sec_structure",
                                "dot_plot")
        self.helper.remove_all_content(sec_path, None, "dir")
        self.helper.remove_all_content(dot_path, None, "dir")
        self.helper.remove_all_content(moun_path, None, "dir")
        for prefix in prefixs:
            tmp_paths = self._get_seq_sec(
                    self.fasta_path, args_srna.out_folder, prefix, sec_path,
                    dot_path, args_srna.vienna_path)
            self._replot_sec_to_pdf(args_srna.vienna_util, tmp_paths,
                                    args_srna.ps2pdf14_path, prefix)
            self._plot_mountain(args_srna.mountain, moun_path, tmp_paths,
                                prefix, args_srna.vienna_util)
            self.helper.remove_all_content(os.getcwd(), ".ps", "file")
            os.chdir(tmp_paths["main"])
            shutil.move("_".join([self.prefixs["energy"], prefix]),
                        "_".join([self.prefixs["basic"], prefix]))
            shutil.rmtree(os.path.join(args_srna.out_folder, "tmp_srna"))

    def _run_blast(self, blast_path, program, database, e, seq_file,
                   blast_file, strand):
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

    def _blast(self, database, database_format, data_type, args_srna,
               prefixs, program, database_type, e):
        if (database is None):
            print("Error: No database assigned!")
        else:
            if database_format:
                self._formatdb(database, data_type, args_srna.out_folder,
                               args_srna.blast_path, database_type)
            for prefix in prefixs:
                blast_file = os.path.join(
                        args_srna.out_folder, "blast_result_and_misc",
                        "_".join([database_type, "blast", prefix + ".txt"]))
                srna_file = "_".join([self.prefixs["basic"], prefix])
                out_file = os.path.join(
                        args_srna.out_folder,
                        "_".join(["tmp", database_type, prefix]))
                print("Running Blast of {0} in {1}".format(prefix, database))
                seq_file = os.path.join(
                        args_srna.out_folder, "_".join(["sRNA_seq", prefix]))
                if seq_file not in os.listdir(args_srna.out_folder):
                    self.helper.get_seq(
                            srna_file,
                            os.path.join(self.fasta_path, prefix + ".fa"),
                            seq_file)
                if database_type == "nr":
                    tmp_plus, tmp_minus = self._get_strand_fasta(
                            seq_file, args_srna.out_folder)
                    tmp_blast = os.path.join("tmp_blast.txt")
                    self._run_blast(args_srna.blast_path, program, database, e,
                                    tmp_plus, tmp_blast, "plus")
                    self._run_blast(args_srna.blast_path, program, database, e,
                                    tmp_minus, blast_file, "minus")
                    self.helper.merge_file(tmp_blast, blast_file)
                    os.remove(tmp_blast)
                    os.remove(tmp_plus)
                    os.remove(tmp_minus)
                else:
                    self._run_blast(args_srna.blast_path, program, database, e,
                                    seq_file, blast_file, "both")
                extract_blast(blast_file, srna_file, out_file,
                              out_file + ".csv", database_type)
                shutil.move(out_file, srna_file)

    def _class_srna(self, prefixs, args_srna):
        '''classify the sRNA based on the filters'''
        if (args_srna.import_info is not None) or (
                args_srna.srna_database is not None) or (
                args_srna.nr_database is not None) or (
                self.sorf_path is not None) or (
                self.tss_path is not None) or (
                self.term_path is not None) or (
                args_srna.promoter_table is not None):
            for prefix in prefixs:
                print("classifying sRNA of {0}".format(prefix))
                class_gff = os.path.join(self.gff_output, "for_class")
                class_table = os.path.join(self.table_output, "for_class")
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
                for srna in os.listdir(class_gff):
                    out_table = os.path.join(
                            class_table, srna.replace(".gff", ".csv"))
                    gen_srna_table(
                        os.path.join(class_gff, srna),
                        "_".join([self.prefixs["merge_table"], prefix]),
                        "_".join([self.tmps["nr"], prefix + ".csv"]),
                        "_".join([self.tmps["srna"], prefix + ".csv"]),
                        args_srna, out_table, self.term_path)

    def _get_best_result(self, prefixs, args_srna):
        '''get the best results based on the filters'''
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

    def _remove_file(self, args_srna):
        self.helper.remove_all_content(args_srna.out_folder, "tmp_", "dir")
        self.helper.remove_all_content(args_srna.out_folder, "tmp_", "file")
        self.helper.remove_tmp(args_srna.fastas)
        self.helper.remove_tmp(args_srna.gffs)
        self.helper.remove_tmp(self.gff_output)
        if args_srna.frag_wigs is not None:
            self.helper.remove_tmp(args_srna.frag_wigs)
        if args_srna.tex_wigs is not None:
            self.helper.remove_tmp(args_srna.tex_wigs)
        if (args_srna.frag_wigs is not None) and (
                args_srna.tex_wigs is not None):
            shutil.rmtree(args_srna.merge_wigs)
        self.helper.remove_tmp(args_srna.trans)
        if args_srna.tss_folder is not None:
            self.helper.remove_tmp(args_srna.tss_folder)
        if args_srna.pro_folder is not None:
            self.helper.remove_tmp(args_srna.pro_folder)
        if args_srna.sorf_file is not None:
            self.helper.remove_tmp(args_srna.sorf_file)
        if "tmp_median" in os.listdir(args_srna.out_folder):
            os.remove(os.path.join(args_srna.out_folder, "tmp_median"))
        if self.term_path is not None:
            self.helper.remove_tmp(args_srna.terms)

    def _filter_srna(self, args_srna, prefixs):
        '''set the filter of sRNA'''
        if args_srna.import_info is not None:
            if "sec_str" in args_srna.import_info:
                self._compute_2d_and_energy(args_srna, prefixs)
        if args_srna.nr_database is not None:
            self._blast(args_srna.nr_database, args_srna.nr_format, "prot",
            args_srna, prefixs, "blastx", "nr", args_srna.e_nr)
        if self.sorf_path is not None:
            for prefix in prefixs:
                if ("_".join([prefix, "sORF.gff"]) in
                        os.listdir(self.sorf_path)):
                    tmp_srna = os.path.join(args_srna.out_folder,
                                            "".join(["tmp_srna_sorf", prefix]))
                    tmp_sorf = os.path.join(args_srna.out_folder,
                                            "".join(["tmp_sorf_srna", prefix]))
                    srna_sorf_comparison(
                            "_".join([self.prefixs["basic"], prefix]),
                            os.path.join(self.sorf_path,
                                         "_".join([prefix, "sORF.gff"])),
                            tmp_srna, tmp_sorf)
                    os.remove(tmp_sorf)
                    shutil.move(tmp_srna,
                                "_".join([self.prefixs["basic"], prefix]))
        if args_srna.srna_database is not None:
            self._blast(args_srna.srna_database, args_srna.srna_format, "nucl",
                        args_srna, prefixs, "blastn", "sRNA", args_srna.e_srna)

    def _import_info_format(self, import_info):
        new_info = []
        for info in import_info:
            info = info.lower()
            new_info.append(info)
        return new_info

    def _gen_table(self, prefixs, args_srna):
        for prefix in prefixs:
            out_table = os.path.join(self.all_best["all_table"],
                                     "_".join([prefix, "sRNA.csv"]))
            gen_srna_table(os.path.join(self.all_best["all_gff"],
                           "_".join([prefix, "sRNA.gff"])),
                           "_".join([self.prefixs["merge_table"], prefix]),
                           "_".join([self.tmps["nr"], prefix + ".csv"]),
                           "_".join([self.tmps["srna"], prefix + ".csv"]),
                           args_srna, out_table, self.term_path)

    def _print_rank_all(self, prefixs):
        for prefix in prefixs:
            all_table = os.path.join(self.all_best["all_table"],
                                     "_".join([prefix, "sRNA.csv"]))
            best_table = os.path.join(self.all_best["best_table"],
                                      "_".join([prefix, "sRNA.csv"]))
            print_rank_all(all_table, best_table)

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

    def _blast_stat(self, stat_path, srna_tables):
        '''do statistics for blast result'''
        for srna_table in os.listdir(os.path.join(srna_tables, "best")):
            out_srna_blast = os.path.join(
                    stat_path, "stat_" +
                    srna_table.replace(".csv", "_blast.csv"))
        blast_class(os.path.join(srna_tables, "best", srna_table),
                    out_srna_blast)

    def _compare_term_promoter(self, out_table, prefix, args_srna):
        '''compare sRNA with terminator and promoter'''
        if self.term_path is not None:
            compare_srna_term(os.path.join(self.all_best["all_gff"],
                              "_".join([prefix, "sRNA.gff"])),
                              out_table, os.path.join(self.term_path,
                              "_".join([prefix, "term.gff"])),
                              args_srna.fuzzy_b, args_srna.fuzzy_a)
        if (args_srna.promoter_table is not None):
            compare_srna_promoter(os.path.join(self.all_best["all_gff"],
                                  "_".join([prefix, "sRNA.gff"])),
                                  out_table, args_srna)

    def run_srna_detection(self, args_srna):
        self._check_necessary_file(args_srna)
        self.multiparser.parser_gff(args_srna.trans, "transcript")
        self.multiparser.combine_gff(args_srna.gffs, self.tran_path,
                                     None, "transcript")
        if args_srna.import_info is not None:
            args_srna.import_info = self._import_info_format(args_srna.import_info)
        prefixs = self._run_program(args_srna)
        self._filter_srna(args_srna, prefixs)
        for prefix in prefixs:
            shutil.copyfile("_".join([self.prefixs["basic"], prefix]),
                            os.path.join(self.all_best["all_gff"],
                            "_".join([prefix, "sRNA.gff"])))
            self._compare_term_promoter("_".join([self.prefixs["merge_table"],
                                        prefix]), prefix, args_srna)
        self._gen_table(prefixs, args_srna)
        self._class_srna(prefixs, args_srna)
        self._get_best_result(prefixs, args_srna)
        self._print_rank_all(prefixs)
        if args_srna.srna_database is not None:
            if "blast_srna" in args_srna.import_info:
                self._blast_stat(self.stat_path, self.table_output)
        self._remove_file(args_srna)
