#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
from subprocess import call
import transaplib.multiparser
import transaplib.helper
class sRNA_detection(object):

    def _get_correct_file(self, folder, feature, prefix):
        detect = False
        for file_ in os.listdir(folder):
            if (file_.replace(feature, "") == prefix):
                detect = True
                return folder + file_
        if detect is False:
            print("Error: there is no proper file - " + prefix + feature)
            sys.exit()

    def _check_gff(self, gffs):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                check_gff_attributes._check_uni_attributes(gffs + "/" + gff)

    def _check_make_folder(self, path, folder):
        if folder in os.listdir(path):
            call(["rm", "-rf", path + folder])
        call(["mkdir", path + folder])

    def _remove_tmp(self, folder):
        if folder is not False:
            call(["rm", "-rf", folder + "/tmp"])
            os.system("rm -rf " + folder + "/*_folder")

    def _multiparser(self, type_, folder, tar_feature, ref_path, ref_feature):
        if (folder is not False) and (folder is not None):
            if type_ == "gff":
                tar_path = folder + "/tmp/"
                multiparser._parser_gff(folder, tar_feature)
                multiparser._combine_gff(ref_path, tar_path, ref_feature, tar_feature)
            elif type_ == "fasta":
                tar_path = folder + "/tmp/"
                multiparser._parser_fasta(folder)
                multiparser._combine_fasta(ref_path, tar_path, ref_feature)
            elif type_ == "wig":
                tar_path = folder + "/tmp/"
                multiparser._parser_wig(folder)
                multiparser._combine_wig(ref_path, tar_path, ref_feature)
        else:
            print("Error: no %s file!!" % tar_feature)
            sys.exit()
        return tar_path

    def _formatdb(self, database, type_, out_folder, blast_path):
        filenames = database.split(".")
        db_file = ".".join(filenames[0:-1])
        err = open(out_folder + "/log.txt", "w")
        call([blast_path + "/makeblastdb", "-in", database,
              "-dbtype", type_, "-out", db_file], stderr=err)

    def _blast_db(self, prefixs, fasta_path, out_folder, bin_path, 
                  program, type_, blast_path, db_file):
        for prefix in prefixs:
            print("Running Blast of " + prefix)
            seq_file = out_folder + "/sRNA_seq_" + prefix
            if seq_file not in os.listdir(out_folder):
                seq_out = open(seq_file, "w")
                call(["python", bin_path + "/get_seq.py",
                      "-s", "tmp_basic_" + prefix, "-f", fasta_path + prefix + ".fa"], stdout=seq_out)
            call([blast_path + "/" + program, "-db", db_file,
                  "-evalue", str(0.0001), "-query", seq_file,
                  "-out", out_folder + "/blast_result_and_misc/" + type_ + "_blast_" + prefix + ".txt"])
            call(["python", bin_path + "/extract_blast.py",
                  "-b", out_folder + "/blast_result_and_misc/" + type_ + "_blast_" + prefix + ".txt",
                  "-s", "tmp_basic_" + prefix,
                  "-o", "tmp_" + type_ + "_" + prefix,
                  "-t", "tmp_" + type_ + "_" + prefix + ".csv",
                  "-d", type_])
            os.system("mv tmp_" + type_ + "_" + prefix + " tmp_basic_" + prefix)

    def _run_normal(self, import_info, tss_path, prefix, bin_path, gff_path,
                    gff, tran, fuzzy, max_len, min_len, wig_path, coverage, 
                    merge_wigs, libs, tex_notex, replicates, table_best,
                    decrease_inter, fuzzy_inter):
        command = " ".join(["python", bin_path + "/sRNA_detect.py",
                  "-g", gff_path + gff, "-a", tran, "-f", str(fuzzy),
                  "-o", "tmp_normal_" + prefix, "-M", str(max_len),
                  "-m", str(min_len), "-wf", wig_path + prefix + "_forward.wig",
                  "-wr", wig_path + prefix + "_reverse.wig",
                  "-b", merge_wigs, "-l", " ".join(libs),
                  "-te", str(tex_notex), "-r", str(replicates),
                  "-ot", "tmp_normal_table_" + prefix, "-cf", str(coverage),
                  "-d", str(decrease_inter), "-fe", str(fuzzy_inter)])
        if ("1" in import_info):
            tss = self._get_correct_file(tss_path, "_TSS.gff", prefix)
            if table_best is True:
                os.system(command + " -t " + tss + " -tb")
            else:
                os.system(command + " -t " + tss)
        else:
            if table_best is True:
                os.system(command + " -tb")
            else:
                os.system(command)

    def _run_utrsrna(self, bin_path, gff_path, gff, tran, fuzzy, merge_wigs,
                     max_len, min_len, wig_path, prefix, tss, pro, fasta_path,
                     libs, tex_notex, replicates, table_best, decrease_utr, 
                     fuzzy_utr, utr5_coverage, utr3_coverage, interCDS_coverage, err):
        command = " ".join(["python", bin_path + "/get_utr_srna.py",
                  "-g", gff_path + gff, "-a", tran, "-f", str(fuzzy),
                  "-b", merge_wigs, "-M", str(max_len), "-m", str(min_len),
                  "-wf", wig_path + prefix + "_forward.wig", 
                  "-wr", wig_path + prefix + "_reverse.wig",
                  "-t", tss, "-p", pro, "-s", fasta_path + prefix + ".fa",
                  "-l", " ".join(libs), "-te", str(tex_notex),
                  "-r", str(replicates), "-o", "tmp_utrsrna_" + prefix,
                  "-ot", "tmp_utrsrna_table_" + prefix, "-cf5", str(utr5_coverage),
                  "-cf3", str(utr3_coverage), "-cfc", str(interCDS_coverage), 
                  "-d", str(decrease_utr), "-fe", str(fuzzy_utr)])
        if table_best is True:
            os.system(command + " -tb > err" + prefix)
        else:
            os.system(command + " > err" + prefix)

    def run_sRNA_detection(self, bin_path, out_folder, utr_srna, gffs, tsss, trans, fuzzy,
                           import_info, tex_wigs, frag_wigs, pros, fastas, mountain, database_format,
                           srna_database, nr_database, energy, coverage, utr5_coverage, utr3_coverage, interCDS_coverage,
                           max_len, min_len, tlibs, flibs, replicates, tex_notex, table_best, decrease_inter,
                           decrease_utr, fuzzy_inter, fuzzy_utr, nr_hits_num, sORF, all_hit, best_sORF):
        vienna_path = os.environ["Vienna_HOME"]
        blast_path = os.environ["BLAST_Plus_HOME"]
        gff_output = out_folder + "/gffs/"
        table_output = out_folder + "/tables/"
        stat_path = out_folder + "/statistics/"
        gff_path = gffs + "/"
        multiparser._parser_gff(gffs, None)
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
                print("Warning: lack Processing site files for UTR derived sRNA detection, it may effect the results!!!!")
#        self._check_gff(gffs)
        if tsss is not False:
            self._check_gff(tsss)
            tss_path = self._multiparser("gff", tsss, "TSS", gffs, None)
#        self._check_gff(trans)
#        if pros is not False:
#            self._check_gff(pros)
        if sORF is not False:
#            self._check_gff(sORF)
            sorf_path = self._multiparser("gff", sORF, "sORF", gffs, None)
        if "1" not in import_info:
            tss_path = None
        if pros is not False:
            pro_path = self._multiparser("gff", pros, "processing", gffs, None)
        if utr_srna or ("2" in import_info) or \
           ("3" in import_info) or ("4" in import_info):
            if fastas is False:
                print("Error: lack required fasta files for UTR derived sRNA detection!!!!")
                sys.exit()
            fasta_path = self._multiparser("fasta", fastas, "fasta", gffs, None)
        if (tex_wigs is not False) and (frag_wigs is not False):
            merge_wigs = os.getcwd() + "/merge_wigs"
            self._check_make_folder( os.getcwd() + "/", "merge_wigs")
            os.system("cp " + tex_wigs + "/* merge_wigs/")
            os.system("cp " + frag_wigs + "/* merge_wigs/")
        elif (tex_wigs is not False):
            merge_wigs = tex_wigs
        elif (frag_wigs is not False):
            merge_wigs = frag_wigs
        wig_path = self._multiparser("wig", merge_wigs, "wig", gffs, None)
        tran_path = self._multiparser("gff", trans, "transcript", gffs, None)
        ###########################################################
        # combine tex_notex and frag libraries .                  #
        ###########################################################
        if (tlibs is False) and (flibs is False):
            print("Error: please input proper libraries!!")
        if (tlibs is not False) and (flibs is not False):
            libs = tlibs + flibs
        elif (tlibs is not False):
            libs = tlibs
        elif (flibs is not False):
            libs = flibs
        ################################################################
        # Compare transcript assembly to find sRNA candidates.         #
        # if TSS wig information offered, it also add the information. #
        ################################################################
        prefixs = []
        for gff in os.listdir(gff_path):
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                prefixs.append(prefix)
                print("Running sRNA detection of " + prefix +"....")
                tran = self._get_correct_file(tran_path, "_transcript.gff", prefix)
                self._run_normal(import_info, tss_path, prefix, bin_path, gff_path,
                                 gff, tran, fuzzy, max_len, min_len, wig_path, coverage, 
                                 merge_wigs, libs, tex_notex, replicates, table_best,
                                 decrease_inter, fuzzy_inter)
                if utr_srna:
                    print("Running UTR derived sRNA detection of " + prefix +"....")
                    tss = self._get_correct_file(tss_path, "_TSS.gff", prefix)
                    pro = self._get_correct_file(pro_path, "_processing.gff", prefix)
                    self._run_utrsrna(bin_path, gff_path, gff, tran, fuzzy, merge_wigs,
                                      max_len, min_len, wig_path, prefix, tss, pro, fasta_path,
                                      libs, tex_notex, replicates, table_best,
                                      decrease_utr, fuzzy_utr, utr5_coverage, utr3_coverage, 
                                      interCDS_coverage, "err")
                    os.system("rm median")
                    ###########################################
                    # merge intergenic and UTR derived sRNA   #
                    ###########################################
                    print("merging data of intergenic and UTR_derived sRNA...")
                    merge_gff = open("tmp_merge_" + prefix, "w")
                    call(["python", bin_path + "/merge_sRNA.py",
                          "-u", "tmp_utrsrna_" + prefix,
                          "-i", "tmp_normal_" + prefix],
                          stdout=merge_gff)
                    merge_table = open("tmp_merge_table_" + prefix, "w")
                    command = " ".join(["python", bin_path + "/merge_table_sRNA.py",
                              "-g", "tmp_merge_" + prefix,
                              "-ti", "tmp_normal_table_" + prefix,
                              "-tu", "tmp_utrsrna_table_" + prefix,
                              "-wf", wig_path + prefix + "_forward.wig",
                              "-wr", wig_path + prefix + "_reverse.wig",
                              "-b", merge_wigs, "-l", " ".join(libs),
                              "-te", str(tex_notex), "-r", str(replicates)])
                    if table_best is True:
                        os.system(command + " -tb > tmp_merge_table_" + prefix)
                    else:
                        os.system(command + " > tmp_merge_table_" + prefix)
                else:
                    os.system("cp tmp_normal_" + prefix + " tmp_merge_" + prefix)
                    os.system("cp tmp_normal_" + prefix + " tmp_merge_table_" + prefix)
                Sort_GFF.sort(self, "tmp_merge_" + prefix, "tmp_basic_" + prefix)
        if "2" in import_info:
            detect = False
            print("Running energy calculation....")
            moun_path = out_folder + "/mountain_plot/"
            sec_path = out_folder + "/sec_structure/sec_plot/"
            dot_path = out_folder + "/sec_structure/dot_plot/"
            os.system("rm -rf " + sec_path + "*")
            os.system("rm -rf " + dot_path + "*")
            os.system("rm -rf " + moun_path + "*")
            for prefix in prefixs:
                for fasta in os.listdir(fasta_path):
                    if fasta.endswith(".fa") and \
                       (fasta.replace(".fa", "") == prefix):
                        detect = True
                        break
                if detect:
                    detect = False
                    seq_file = out_folder + "/sRNA_seq_" + prefix
                    sec_file = out_folder + "/sRNA_2d_" + prefix
                    seq = open(seq_file, "w")
                    sec = open(sec_file, "w")
                    call(["python", bin_path + "/get_seq.py",
                          "-s", "tmp_basic_" + prefix, "-f", fasta_path + fasta], stdout=seq)
                else:
                    print("Error:There is not fasta file of " + prefix)
                self._check_make_folder(os.getcwd() + "/", "tmp_srna")
                tmp_path = os.getcwd() + "/tmp_srna/"
                os.chdir(tmp_path)
                sec_file = "../" + sec_file
                seq_file = "../" + seq_file
                out = open("tmp_energy_" + prefix, "w")
                tmp_sec_path = "../" + sec_path
                tmp_dot_path = "../" + dot_path
                os.system("cat " + seq_file + " | " + vienna_path + "/Progs/RNAfold -p > " + sec_file)
                call(["python", bin_path + "/merge_sRNA_energy.py",
                      "-s", "../tmp_basic_" + prefix, "-d", sec_file], stdout=out)
                for file_ in os.listdir(os.getcwd()):
                    if file_[-5:] == "ss.ps":
                        dot_file = file_[0:-5] + "dp.ps"
                        rel_file = file_[0:-5] + "rss.ps"
                        out = open(rel_file, "w")
                        print("replot " + file_)
                        call([vienna_path + "/Utils/relplot.pl", file_, dot_file], stdout=out)
                for file_ in os.listdir(os.getcwd()):
                    if file_[-3:] == ".ps":
                        pdf_file = file_[0:-2] + "pdf"
                        print("convert " + file_ + " to pdf")
                        call(["ps2pdf14", file_, pdf_file])
                call(["mkdir", tmp_sec_path + prefix])
                call(["mkdir", tmp_dot_path + prefix])
                os.system("mv *rss.pdf " + tmp_sec_path + prefix)
                os.system("mv *dp.pdf " + tmp_dot_path + prefix)
                if mountain is True:
                    tmp_moun_path = "../" + moun_path
                    call(["mkdir", tmp_moun_path + prefix])
                    print("Generating mountain plot of " + prefix + "....")
                    for file_ in os.listdir(os.getcwd()):
                        if "dp.ps" in file_:
                            out = open("mountain.txt", "w")
                            moun_file = file_[0:-5] + "mountain.pdf"
                            print("Generating " + moun_file)
                            call([vienna_path + "/Utils/mountain.pl", file_],
                                  stdout=out)
                            call(["python", bin_path + "/plot_mountain.py",
                                  "-i", "mountain.txt", "-o", moun_file])
                            call(["mv", moun_file, tmp_moun_path + prefix])
                            out.close()
                    os.system("rm mountain.txt")
                os.system("rm *.ps")
                os.system("rm *ss.pdf")
                os.chdir("../")
                call(["mv", "tmp_srna/tmp_energy_" + prefix, "tmp_basic_" + prefix])
                call(["rm", "-rf", "tmp_srna"])
        if "3" in import_info:
            db_file = nr_database
            if (nr_database is False):
                print("Error: No database assigned!")
            else:
                if database_format is True:
                    self._formatdb(nr_database, "prot", out_folder, blast_path)
                self._blast_db(prefixs, fasta_path, out_folder, bin_path, 
                               "blastx", "nr", blast_path, db_file)
        if "4" in import_info:
            db_file = srna_database
            if (srna_database is False):
                print("Error: No database assigned!")
            else:
                if database_format is True:
                    self._formatdb(srna_database, "nucl")
                self._blast_db(prefixs, fasta_path, out_folder, bin_path, 
                               "blastn", "sRNA", blast_path, db_file)
            ################################
            # statistics of sRNA blast hit #
            ################################
            out_srna_blast = open(stat_path + "stat_sRNA_blast_class_" + prefix + ".csv", "w")
            call(["python", bin_path + "/blast_class.py",
                  "-i", "tmp_sRNA_" + prefix + ".csv"], stdout=out_srna_blast)
        if "5" in import_info:
            if prefix + "_sORF.gff" in os.listdir(sorf_path):
                call(["python", bin_path + "/compare_sRNA_sORF.py",
                      "-r", "tmp_basic_" + prefix,
                      "-o", sorf_path + prefix + "_sORF.gff",
                      "-or", "tmp_srna_sorf" + prefix,
                      "-oo", "tmp_sorf_srna" + prefix])
                os.system("rm tmp_sorf_srna" + prefix)
                os.system("mv tmp_srna_sorf" + prefix + " tmp_basic_" + prefix)
        for prefix in prefixs:
            os.system("cp tmp_basic_" + prefix + " " + gff_output + "all_candidates/" + prefix + "_sRNA.gff")
        #######################
        # generate table      #
        #######################
        for prefix in prefixs:
            out_table = open(table_output + "all_candidates/" + prefix + "_sRNA.csv", "w")
            call(["python", bin_path + "/gen_srna_table.py",
                  "-g", gff_output + "all_candidates/" + prefix + "_sRNA.gff",
                  "-t", "tmp_merge_table_" + prefix,
                  "-n", "tmp_nr_" + prefix + ".csv",
                  "-s", "tmp_sRNA_" + prefix + ".csv",
                  "-M", str(max_len), "-m", str(min_len)], stdout=out_table)

        #######################
        #  classification     #
        #######################
        if (len(import_info) != 1) or ('6' not in import_info):
            for prefix in prefixs:
                print("classifying sRNA of " + prefix)
                class_gff = gff_output + "for_class/"
                class_table = table_output + "for_class/"
                self._check_make_folder(class_table, prefix)
                self._check_make_folder(class_gff, prefix)
                class_gff = class_gff + prefix
                class_table = class_table + prefix
                out = open(stat_path + "/stat_sRNA_class_" + prefix + ".csv", "w")
                call(["python", bin_path + "/sRNA_class.py",
                      "-i", gff_output + "all_candidates/" + prefix + "_sRNA.gff",
                      "-f", class_gff,
                      "-e", str(energy),
                      "-cn", str(nr_hits_num)], stdout=out)
                for file_ in os.listdir(class_gff):
                    out_table = open(class_table + "/" + file_.replace(".gff", ".csv"), "w")
                    call(["python", bin_path + "/gen_srna_table.py",
                          "-g", class_gff + "/" + file_,
                          "-t", "tmp_merge_table_" + prefix,
                          "-n", "tmp_nr_" + prefix + ".csv",
                          "-s", "tmp_sRNA_" + prefix + ".csv",
                          "-M", str(max_len), "-m", str(min_len)], stdout=out_table)
        ###################################
        # select best sRNA                #
        ###################################
        for prefix in prefixs:
            best_gff = open(gff_output + "best/" + prefix + "_sRNA.gff", "w")
            best_table = open(table_output + "best/" + prefix + "_sRNA.csv", "w")
            command =["python", bin_path + "/gen_best_sRNA.py",
                      "-g", gff_output + "all_candidates/" + prefix + "_sRNA.gff",
                      "-e", str(energy), "-cn", str(nr_hits_num)]
            if (all_hit is True) and (best_sORF is True):
                call(command + ["-a", "-so"], stdout=best_gff)
            elif (all_hit is True):
                call(command + ["-a"], stdout=best_gff)
            elif (best_sORF is True):
                call(command + ["-so"], stdout=best_gff)
            else:
                call(command, stdout=best_gff)
            call(["python", bin_path + "/gen_srna_table.py",
                  "-g", gff_output + "best/" + prefix + "_sRNA.gff",
                  "-t", "tmp_merge_table_" + prefix,
                  "-n", "tmp_nr_" + prefix + ".csv",
                  "-s", "tmp_sRNA_" + prefix + ".csv",
                  "-M", str(max_len), "-m", str(min_len)], stdout=best_table)
        ######################################
        # remove temperary folders and files.#
        ######################################
        os.system("rm tmp_*")
        os.system("rm merge_wigs")
        self._remove_tmp(fastas)
        self._remove_tmp(gffs)
        self._remove_tmp(wigs)
        self._remove_tmp(tsss)
        self._remove_tmp(trans)
        self._remove_tmp(pros)
        self._remove_tmp(sORF)
        os.system("rm median")
