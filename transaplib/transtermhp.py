#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
from subprocess import call, Popen
sys.path.append(os.environ["Transap_BIN"])
import sort_gff
import multiparser
class TransTermHP(object):

    def _get_correct_file(self, datas, prefix, feature, for_wig_type):
        detect = False
        for data in os.listdir(datas):
            if for_wig_type == "":
                if feature in data:
                    file_ = data[:-1 * len(feature)]
                    if prefix == file_:
                        detect = True
                        return datas + data
            else:
                filename = data.split("_STRAIN_")
                if ("reverse" in data) and ("forward" in data):
                    print("Error: Unclear wig file. It is reverse or forward!!!")
                    sys.exit()
                elif (prefix == filename[-1][:-1 * len(feature)]) and \
                     (for_wig_type in data):
                    return datas + data
        if detect:
            detect = False
        else:
            print("Warning: no proper file - " + prefix + feature)

    def _combine_annotation(self, combine_file, files):
        with open(combine_file, 'w') as result:
            for file_ in files:
                check_start = False
                for line in open( file_, 'r' ):
                    if check_start:
                        result.write(line)
                    if "Location" in line:
                        check_start = True

    def _remove_tmp(self, folder):
        if folder is not False:
            call(["rm", "-rf", folder + "/tmp"])
            os.system("rm -rf " + folder + "/*_folder")

    def run_TransTermHP(self, bin_path, out_folder, fastas, gffs, sRNAs,
                        combine, stat, wigs, decrease, fuzzy, output_fuzzy):
        #################################################################
        # Before, running TransTermHP,                                  #
        # we need to convert annotation files to .ptt and .rnt files.   #
        #################################################################
        file_types = {}
        prefixs = []
        tran_path = os.environ["TransTermHP_HOME"]
        multiparser._parser_gff(gffs, None)
        multiparser._parser_fasta(fastas)
        gff_path = gffs + "/tmp/"
        fasta_path = fastas + "/tmp/"
        term_outfolder = out_folder + "/Term_gff/"
        if gff_path is None:
            print("Error: please assign the annotation gff folder!!")
            sys.exit()
        for gff in os.listdir(gff_path):
            if gff.endswith(".gff"):
                filename = gff.split("/")
                prefix = filename[-1][:-4]
                prefixs.append(prefix)
                rnt = gff[:-3] + "rnt"
                ptt = gff[:-3] + "ptt"
                fasta = self._get_correct_file(fasta_path, prefix, ".fa", "")
                if fasta is None:
                    print("Error: no proper file - " + prefix + ".fa")
                    sys.exit()
                if sRNAs is not False:
                    multiparser._parser_gff(sRNAs, "sRNA")
                    srna_path = sRNAs + "/tmp/"
                    srna = self._get_correct_file(srna_path, prefix, "_sRNA.gff", "")
                    if (srna is not None) and (fasta is not None):
                        call(["python", bin_path + "/gff2rntptt.py",
                              "-g", gff_path + gff, 
                              "-f", fasta,
                              "-p", gff_path + ptt, 
                              "-r", gff_path + rnt,
                              "-si", srna, "-so", srna[:-3] + "rnt"])
                        file_types[prefix] = "srna"
                    if (srna is None) and (fasta is not None):
                        call(["python", bin_path + "/gff2rntptt.py",
                              "-g", gff_path + gff,
                              "-f", fasta,
                              "-p", gff_path + ptt,
                              "-r", gff_path + rnt])
                        file_types[prefix] = "normal"
                else:
                    if fasta is not None:
                        call(["python", bin_path + "/gff2rntptt.py",
                              "-g", gff_path + gff, 
                              "-f", fasta,
                              "-p", gff_path + ptt,
                              "-r", gff_path + rnt])
                        file_types[prefix] = "normal"

        ###################################################################
        # Then, we still need to combine all .ptt and .rnt to be one file.#
        ###################################################################
        combine_path = gff_path + "/combine/"
        if "combine" in os.listdir(gff_path):
            call(["rm", "-rf", combine_path])
        call(["mkdir", combine_path])
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

        ##################################################################
        # Start to run TransTermHP.                                      #
        ##################################################################
        output_folders = []
        for file_ in os.listdir(combine_path):
            if ".ptt" in file_:
                prefix = file_.replace(".ptt", "")
                fasta = self._get_correct_file(fasta_path, prefix, ".fa", "")
                if fasta is None:
                    print("Error: no proper file - " + prefix + ".fa")
                    sys.exit()
                out_path = out_folder + "/" + prefix + "/"
                if prefix in os.listdir(out_folder):
                    call(["rm", "-rf", out_path])
                call(["mkdir", out_path])
                output_folders.append(out_path)
                out = open(out_path + prefix + "_terminators.txt", "w")
                call([tran_path + "/transterm",
                      "-p", tran_path + "/expterm.dat",
                      fasta, combine_path + file_, 
                      "--t2t-perf", out_path + prefix + "_terminators_within_robust_tail-to-tail_regions.t2t",
                      "--bag-output", out_path + prefix + "_best_terminator_after_gene.bag"],
                      stdout=out)
        call(["rm", "-rf", combine_path])
        ################################################################
        # Transfer to gff.                                             #
        ################################################################
        for prefix in prefixs:
            for folder in os.listdir(out_folder):
                if prefix == folder:
                    out_path = out_folder + "/" + folder + "/"
                    for file_ in os.listdir(out_path):
                        if file_.endswith(".bag"):
                            out_file = term_outfolder + prefix + "_term.gff"
                            out_f = open(out_file, "w")
                            call(["python", 
                                  bin_path + "/convert_transtermhp_bag_file_to_gff.py",
                                  out_path + file_], stdout=out_f)
                            out_f.close()
        multiparser._combine_gff(gffs, term_outfolder, None, "term")
        #################################################################
        # Statistics of TransTermHP                                     #
        #################################################################
        if stat is not False:
            wig_path = wigs + "/tmp/"
            stat_path = out_folder + "/statistics/"
            multiparser._parser_gff(term_outfolder, "term")
            multiparser._parser_wig(wigs)
            multiparser._combine_wig(term_outfolder, wig_path, "term")
            print("Running statistics.......\n")
            for term in os.listdir(term_outfolder):
                if term.endswith("_term.gff"):
                    prefix = term.replace("_term.gff", "")
                    wig_f = self._get_correct_file(wig_path, prefix, "_forward.wig", "forward")
                    wig_r = self._get_correct_file(wig_path, prefix, "_reverse.wig", "reverse")
                    if (wig_f) and (wig_r):
                        if output_fuzzy is not False:
                            call(["python", bin_path + "/stat_term.py",
                                  "-t", term_outfolder + term,
                                  "-wf", wig_f,
                                  "-wr", wig_r,
                                  "-o", out_folder + "/fuzzy_gff/" + prefix + "_term_no_fuzzy.gff",
                                  "-of", out_folder + "/fuzzy_gff/" + prefix + "_term_fuzzy.gff",
                                  "-s", stat_path + "stat_terminator_" + prefix + ".csv",
                                  "-d", str(decrease), "-f", str(fuzzy)])
                        else:
                            call(["python", bin_path + "/stat_term.py",
                                  "-t", term_outfolder + term,
                                  "-wf", wig_f,
                                  "-wr", wig_r,
                                  "-s", stat_path + "stat_terminator_" + prefix + ".csv",
                                  "-d", str(decrease), "-f", str(fuzzy)])
        ######################################################
        # Remove the temperary folders and files             #
        ######################################################
        self._remove_tmp(gffs)
#        self._remove_tmp(fastas)
        self._remove_tmp(sRNAs)
        self._remove_tmp(out_folder + "/Term_gff")
        self._remove_tmp(wigs)
