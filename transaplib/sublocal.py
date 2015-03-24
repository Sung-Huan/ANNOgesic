#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
from subprocess import call
import multiparser
import check_gff_attributes
class Sub_Local(object):

    def _check_make_folder(self, path, folder):
        if folder in os.listdir(path):
            call(["rm", "-rf", path + folder])
        call(["mkdir", path + folder])

    def _get_correct_file(self, folder, feature, prefix):
        detect = False
        for file_ in os.listdir(folder):
            if (file_.replace(feature, "") == prefix):
                detect = True
                return folder + file_
        if detect is False:
            print("Error: there is no proper file - " + prefix + feature)
            sys.exit()

    def run_sub_local(self, bin_path, gffs, fastas, 
                      gram, merge, out_folder):
        emboss_path = os.environ["EMBOSS_HOME"]
        psortb_path = os.environ["PSORTB_HOME"]
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                check_gff_attributes._check_uni_attributes(gffs + "/" + gff)
        multiparser._parser_gff(gffs, None)
        multiparser._parser_fasta(fastas)
        gff_path = gffs + "/tmp/"
        fasta_path = fastas + "/tmp/"
        tmp_path = out_folder + "/" + "tmp/"
        stat_path = out_folder + "/statistics/"
        self._check_make_folder(out_folder + "/", "tmp")
        self._check_make_folder(out_folder + "/", "tmp_results")
        for gff in os.listdir(gff_path):
            prefix = gff.replace(".gff", "")
            fasta = self._get_correct_file(fasta_path, ".fa", prefix)
            dna_seq_file = tmp_path + prefix + "_dna.fa"
            out_dna_seq = open(dna_seq_file, "w")
            #### get CDS seq ####
            print("Generate CDS fasta files of " + prefix)
            call(["python", bin_path + "/get_cds_seq.py",
                  "-g", gff_path + gff,
                  "-f", fasta], stdout= out_dna_seq)
            #### transfer to protein seq ####
            print("transfer DNA seq to protein seq of " + prefix)
            call([emboss_path + "/transeq",
                  "-sequence", dna_seq_file, "-outseq", "tmp", "-trim"])
            protein_seq_file = tmp_path + prefix + "_protein.fa"
            out_protein_seq = open(protein_seq_file, "w")
            call(["python", bin_path + "/fix_emboss_format.py",
                  "-i", "tmp"], stdout=out_protein_seq)
            #### running psortb ####
            print("Running psortb of " + prefix)
            out_err = open("tmp_log", "w")
            tmp_psortb_path = out_folder + "/tmp_results/"
            out_raw = open(tmp_psortb_path + prefix + "_raw.txt", "w")
            if gram == "positive":
                call([psortb_path + "/psort",
                      "-p", protein_seq_file], stdout=out_raw, stderr=out_err)
#                print(psortb_path + "/psort" + "-p", protein_seq_file)
            elif gram == "negative":
                call([psortb_path + "/psort",
                      "-n", protein_seq_file], stdout=out_raw, stderr=out_err)
            else:
                print("Error:It is not a proper bacteria type - " + gram + "!!")
                sys.exit()
            #### extract psortb results ####
            if merge is True:
                print("Merge to gff...")
                call(["python", bin_path + "/extract_psortb.py",
                      "-p", tmp_psortb_path + prefix + "_raw.txt",
                      "-m", gff_path + gff,
                      "-op", tmp_psortb_path + prefix + "_table.csv",
                      "-om", prefix + ".gff"])
                call(["cp", prefix + ".gff", gff_path + gff])
            else:
                call(["python", bin_path + "/extract_psortb.py",
                      "-p", tmp_psortb_path + prefix + "_raw.txt",
                      "-op", tmp_psortb_path + prefix + "_table.csv"])
        #### merge the psortb results based on annotation file ####
        for folder in os.listdir(gffs):
            if folder.endswith(".gff_folder"):
                prefix = folder.replace(".gff_folder", "")
                self._check_make_folder(out_folder + "/psortb_results/", prefix)
                merge_table = out_folder + "/psortb_results/" + prefix + "/" + prefix + "_table.csv"
                for gff in os.listdir(gffs + "/" + folder):
                    result = self._get_correct_file(tmp_psortb_path, "_raw.txt", gff.replace(".gff", ""))
                    call(["cp", result, out_folder + "/psortb_results/" + prefix])
                    result = self._get_correct_file(tmp_psortb_path, "_table.csv", gff.replace(".gff", ""))
                    os.system("cat " + result + " >> " + merge_table)
                #### statistics of psortb data ####
                self._check_make_folder(stat_path, prefix)
                call(["python", bin_path + "/stat_sublocal.py",
                      "-p", merge_table,
                      "-o", stat_path + prefix + "/" + prefix,
                      "-s", stat_path + prefix + "/stat_" + prefix + "_sublocal.csv"])
#        call(["rm", "-rf", fastas + "/tmp"])
#        os.system("rm -rf " + fastas + "/*_folder")
        os.system("rm -rf " + out_folder + "/tmp*")
