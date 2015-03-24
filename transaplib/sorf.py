#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
from subprocess import call
import transaplib.multiparser
import transaplib.helper
class sORF_detection(object):

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

    def run_sORF_detection(self, bin_path, out_folder, utr_detect, trans, gffs, tsss, utr_length, 
                           min_len, max_len, tex_wigs, frag_wigs, cutoff_inter, cutoff_5utr, 
                           cutoff_3utr, cutoff_interCDS, fastas, tlibs, flibs, tex_notex, 
                           replicate, table_best, sRNAs, start_coden, stop_coden, condition):
        if (gffs is None) or \
           (trans is None) or \
           ((tex_wigs is None) and (frag_wigs is None)):
            print("Error: lack required files!!!!")
            sys.exit()
        if utr_detect:
            if (tsss is None):
                print("Error: lack required files for UTR derived sRNA detection!!!!")
                sys.exit()
#        self._check_gff(gffs)
        gff_path = gffs + "/"
        multiparser._parser_gff(gffs, None)
        if tsss is not False:
            self._check_gff(tsss)
            tss_path = self._multiparser("gff", tsss, "TSS", gffs, None)
#        self._check_gff(trans)
        if sRNAs is not False:
            self._check_gff(sRNAs)
            srna_path = self._multiparser("gff", sRNAs, "sRNA", gffs, None)
        gff_output = out_folder + "/gffs/"
        table_output = out_folder + "/tables/"
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
        fasta_path = self._multiparser("fasta", fastas, "fasta", gffs, None)
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
        #################################################
        # compare transcript and annotation             #
        #################################################
        prefixs = []
        for gff in os.listdir(gff_path):
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                prefixs.append(prefix)
                print("comparing transcript and CDS of " + prefix)
                command = ["python", bin_path + "/compare_tran_gff.py",
                           "-g", gff_path + gff, 
                           "-t", tran_path + prefix + "_transcript.gff",
                           "-o", out_folder + "/" + prefix + "_inter.gff"]
                if utr_detect:
                    call(command + ["-u"])
                else:
                    call(command)
        ##########################################################
        # detect start and stop coden to get the sORF candidates #
        ##########################################################
        for prefix in prefixs:
            print("detect sORF of " + prefix)
            command = ["python", bin_path + "/start_stop_codon.py",
                       "-f", fasta_path + prefix + ".fa",
                       "-i", out_folder + "/" + prefix + "_inter.gff",
                       "-uf", str(utr_length), "-l", " ".join(libs),
                       "-te", str(tex_notex), "-r", str(replicate),
                       "-ic", str(cutoff_inter), "-u3", cutoff_3utr,
                       "-u5", cutoff_5utr, "-uc", cutoff_interCDS,
                       "-wf", wig_path + prefix + "_forward.wig",
                       "-wr", wig_path + prefix + "_reverse.wig",
                       "-b", merge_wigs, "-ac", " ".join(start_coden),
                       "-oc", " ".join(stop_coden), "-M", str(max_len),
                       "-m", str(min_len), "-c", condition, 
                       "-op", gff_output + "all_candidates/" + prefix + "_sORF"]
            if tsss is not False:
                command = command + ["-t", tss_path + prefix + "_TSS.gff"]
            if sRNAs is not False:
                command = command + ["-s", srna_path + prefix + "_sRNA.gff"]
            if table_best:
                command = command + ["-tb"]
            if utr_detect:
                command = command + ["-ud"]
            os.system(" ".join(command))
            if prefix + "_sORF_all.gff" in os.listdir(gff_output + "all_candidates/"):
                os.system(" ".join(["mv", "".join([gff_output, "all_candidates/", prefix, "_sORF_all.gff"]), \
                                    "".join([gff_output, "all_candidates/", prefix, "_sORF.gff"])]))
                os.system(" ".join(["mv", "".join([gff_output, "all_candidates/", prefix, "_sORF_best.gff"]), \
                                    "".join([gff_output, "best/", prefix, "_sORF.gff"])]))
                os.system(" ".join(["mv", "".join([gff_output, "all_candidates/", prefix, "_sORF_all.csv"]), \
                                    "".join([table_output, "all_candidates/", prefix, "_sORF.csv"])]))
                os.system(" ".join(["mv", "".join([gff_output, "all_candidates/", prefix, "_sORF_best.csv"]), \
                                    "".join([table_output, "best/", prefix, "_sORF.csv"])]))
        for sorf in os.listdir(gff_output + "all_candidates"):
            print("statistics of " + sorf)
            if sorf.endswith("_sORF.gff"):
                command = ["python", bin_path + "/stat_sorf.py",
                           "-g", gff_output + "all_candidates/" + sorf,
                           "-s", out_folder + "/statistics/stat_" + sorf.replace(".gff", ".csv")]
                if utr_detect:
                    call(command + ["-u"])
                else:
                    call(command)
        ##########################
        # remove temperary files #
        ##########################
        os.system("rm merge_wigs")
        os.system("rm " + out_folder + "/*.gff")
        self._remove_tmp(fastas)
        self._remove_tmp(gffs)
        self._remove_tmp(tsss)
        self._remove_tmp(trans)
        self._remove_tmp(sRNAs)
