#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
from subprocess import call
sys.path.append(os.environ["Transap_BIN"])
import multiparser
import check_gff_attributes
class MEME(object):

    def _run_normal_motif(self, meme_path, input_path, out_path,
                          filename, width, parallel, num_motif, fasta):
        folder = "promoter_motifs_" + filename + "_" + str(width) + "_nt"
        if folder not in os.listdir(out_path):
            call([meme_path + "/meme", "-maxsize", "1000000",
                  "-dna", "-nmotifs", str(num_motif),
                  "-w", str(width), "-maxiter", "100",
                  "-p", str(parallel), "-oc", out_path + folder,
                  input_path + fasta])
    def _run_small_motif(self, meme_path, input_path, out_path,
                         filename, width, parallel, num_motif, fasta):
        data = width.split("-")
        min_width = data[0]
        max_width = data[1]
        folder = "promoter_motifs_" + filename + "_" + str(min_width) + "-" + str(max_width) + "_nt"
        if folder not in os.listdir(out_path):
             call([meme_path + "/meme", "-maxsize", "1000000",
                   "-dna", "-nmotifs", str(num_motif),
                   "-minsites", "0", "-maxsites", "2",
                   "-minw", str(min_width), "-maxw", str(max_width),
                   "-maxiter", "100", "-p", str(parallel),
                   "-oc", out_path + folder,
                   input_path + fasta])

    def _remove_tmp(self, folder):
        if folder is not False:
            call(["rm", "-rf", folder + "/tmp"])
            os.system("rm -rf " + folder + "/*_folder")
        
    def run_MEME(self, bin_path, input_folder, output_folder, libs, tsss,
                 fastas, num_motif, widths, parallel, source, wigs, gffs):
        meme_path = os.environ["MEME_HOME"]
        multiparser._parser_fasta(fastas)
        multiparser._parser_gff(tsss, "TSS")
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                check_gff_attributes._check_uni_attributes(gffs + "/" + gff)
        tss_path = tsss + "/tmp/"
        fasta_path = fastas + "/"
        input_path = input_folder + "/"
        multiparser._combine_gff(fastas, tss_path, "fasta", "TSS")    
        os.system("rm -rf " + input_path + "*")
        ###################################################################
        # Before running MEME,                                            #
        # we have to generate fasta file which is based on the TSS types. #
        ###################################################################
        detect = False
        prefixs = []
        for tss in os.listdir(tss_path):
            prefix = tss.replace("_TSS.gff", "")
            prefixs.append(prefix)
            out_path = output_folder + "/" + prefix + "/"
            if prefix in os.listdir(output_folder):
                call(["rm", "-rf", out_path])
            call(["mkdir", out_path])
            input_path = input_folder + "/" + prefix + "/"
            if prefix in os.listdir(input_folder):
                call(["rm", "-rf", input_path])
            call(["mkdir", input_path])
            for fasta in os.listdir(fasta_path):
                if (fasta.endswith(".fa")) and \
                   (prefix == fasta.replace(".fa", "")):
                    detect = True
                    break
                elif (fasta.endswith(".fna")) and \
                     (prefix == fasta.replace(".fna", "")):
                    detect = True
                    break
                elif (fasta.endswith(".fasta")) and \
                     (prefix == fasta.replace(".fasta", "")):
                    detect = True
                    break
            if source is True:
                print("generating fasta file of " + prefix)
                call(["python", bin_path + "/get_upstream.py",
                      "-t", tss_path + tss, "-f", fasta_path + fasta])
            else:
                multiparser._parser_gff(gffs, None)
                gff_path = gffs + "/tmp/"
                multiparser._combine_gff(fastas, gff_path, "fasta", None)
                if "TSS_class" not in os.listdir(output_folder):
                    os.system("mkdir " + output_folder + "/TSS_class")
                print("classifying TSS and extracting fasta of " + prefix)
                command = " ".join(["python", bin_path + "/get_upstream.py",
                                    "-t", tss_path + tss, "-f", fasta_path + fasta,
                                    "-g", gff_path + prefix + ".gff", "-s",
                                    "-b", wigs, "-l", " ".join(libs)])
                os.system(command + " > " + output_folder + "/TSS_class/" + prefix + "_TSS.gff")
            out_all = open("all_type.fa", "w")
            out_no_orph = open("without_orphan.fa", "w")
            os.system("cat primary.fa > tmp.fa")
            os.system("cat secondary.fa >> tmp.fa")
            os.system("cat internal.fa >> tmp.fa")
            os.system("cat antisense.fa >> tmp.fa")
            os.system("cat tmp.fa > tmp_all.fa")
            os.system("cat orphan.fa >> tmp_all.fa")
            call(["python", bin_path + "/del_repeat_fasta.py", 
                  "-f", "tmp.fa"], stdout=out_no_orph)
            call(["python", bin_path + "/del_repeat_fasta.py", 
                  "-f", "tmp_all.fa"], stdout=out_all)
            call(["rm", "tmp.fa"])
            call(["rm", "tmp_all.fa"])
            os.system("mv primary.fa " + input_path + prefix + "_allstrain_primary.fa")
            os.system("mv secondary.fa " + input_path + prefix + "_allstrain_secondary.fa")
            os.system("mv internal.fa " + input_path + prefix + "_allstrain_internal.fa")
            os.system("mv antisense.fa " + input_path + prefix + "_allstrain_antisense.fa")
            os.system("mv orphan.fa " + input_path + prefix + "_allstrain_orphan.fa")
            os.system("mv all_type.fa " + input_path + prefix + "_allstrain_all_types.fa")
            os.system("mv without_orphan.fa " + input_path + prefix + "_allstrain_without_orphan.fa")
            for fasta in os.listdir(input_path):
                if fasta.endswith(".fa"):
                    pre_strain = ""
                    num_strain = 0
                    with open(input_path + fasta, "r") as f_h:
                        for line in f_h:
                            line = line.strip()
                            if line.startswith(">"):
                                datas = line.split("_")
                                strain = "_".join(datas[2:])
                                if pre_strain != strain:
                                    num_strain += 1
                                    filename = fasta.split("allstrain")
                                    out = open(input_path + filename[0] + strain + filename[-1], "w")
                                    pre_strain = strain
                                out.write(line + "\n")
                            else:
                                out.write(line + "\n")
                    if num_strain <= 1:
                        call(["rm", input_path + filename[0] + strain + filename[-1]])
            out.close()
        for prefix in prefixs:
            input_path = input_folder + "/" + prefix + "/"
            out_path = output_folder + "/" + prefix + "/"
            for fasta in os.listdir(input_path):
                filename = fasta.replace(".fa", "")
                for width in widths:
                    print("Computing promoters of " + fasta + " - " + width)
                    if "-" in width:
                        self._run_small_motif(meme_path, input_path, out_path, 
                                        filename, width, parallel, num_motif, fasta)
                    else:
                        self._run_normal_motif(meme_path, input_path, out_path,
                                        filename, width, parallel, num_motif, fasta)
        self._remove_tmp(fastas)
        self._remove_tmp(tsss)
        self._remove_tmp(gffs)
        self._remove_tmp(wigs)
