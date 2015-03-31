#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
import shutil
from subprocess import call
from transaplib.multiparser import Multiparser
from transaplib.helper import Helper
from transaplib.TSS_upstream import Upstream, Del_repeat_fasta


class MEME(object):

    def __init__(self):
        self.multiparser = Multiparser()
        self.helper = Helper()

    def _run_normal_motif(self, meme_path, input_path, out_path,
                          filename, width, parallel, num_motif, fasta):
        folder = "_".join(["promoter_motifs", filename, str(width), "nt"])
        if folder not in os.listdir(out_path):
            call([meme_path, "-maxsize", "1000000",
                  "-dna", "-nmotifs", str(num_motif),
                  "-w", str(width), "-maxiter", "100",
                  "-p", str(parallel), 
                  "-oc", os.path.join(out_path, folder),
                  os.path.join(input_path, fasta)])
    def _run_small_motif(self, meme_path, input_path, out_path,
                         filename, width, parallel, num_motif, fasta):
        data = width.split("-")
        min_width = data[0]
        max_width = data[1]
        folder = "_".join(["promoter_motifs", filename, 
                 "-".join([str(min_width), str(max_width)]), "nt"])
        if folder not in os.listdir(out_path):
             call([meme_path, "-maxsize", "1000000",
                   "-dna", "-nmotifs", str(num_motif),
                   "-minsites", "0", "-maxsites", "2",
                   "-minw", str(min_width), "-maxw", str(max_width),
                   "-maxiter", "100", "-p", str(parallel),
                   "-oc", os.path.join(out_path, folder),
                   os.path.join(input_path, fasta)])

    def _get_fasta_file(self, fasta_path, prefix):
        for fasta in os.listdir(fasta_path):
            if (fasta.endswith(".fa")) and \
               (prefix == fasta.replace(".fa", "")):
                break
            elif (fasta.endswith(".fna")) and \
                 (prefix == fasta.replace(".fna", "")):
                break
            elif (fasta.endswith(".fasta")) and \
                 (prefix == fasta.replace(".fasta", "")):
                break
        return fasta
        
    def _check_gff(self, gffs):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(gffs, gff))

    def _move_and_merge_fasta(self, input_path, prefix):
        if "tmp/all_type.fa" in os.listdir("tmp"):
            os.remove("tmp/all_type.fa")
        if "tmp/without_orphan.fa" in os.listdir("tmp"):
            os.remove("tmp/without_orphan.fa")
        shutil.copyfile("tmp/primary.fa", "tmp/tmp.fa")
        self.helper.merge_file("tmp", "secondary.fa", "tmp", "tmp.fa")
        self.helper.merge_file("tmp", "internal.fa", "tmp", "tmp.fa")
        self.helper.merge_file("tmp", "antisense.fa", "tmp", "tmp.fa")
        shutil.copyfile("tmp/tmp.fa", "tmp/tmp_all.fa")
        self.helper.merge_file("tmp", "orphan.fa", "tmp", "tmp_all.fa")
        Del_repeat_fasta("tmp/tmp.fa", "tmp/without_orphan.fa")
        Del_repeat_fasta("tmp/tmp_all.fa", "tmp/all_type.fa")
        os.remove("tmp/tmp.fa")
        os.remove("tmp/tmp_all.fa")
        out_prefix = os.path.join(input_path, prefix)
        os.rename("tmp/primary.fa", "_".join([out_prefix, "allstrain_primary.fa"]))
        os.rename("tmp/secondary.fa", "_".join([out_prefix, "allstrain_secondary.fa"]))
        os.rename("tmp/internal.fa", "_".join([out_prefix, "allstrain_internal.fa"]))
        os.rename("tmp/antisense.fa", "_".join([out_prefix, "allstrain_antisense.fa"]))
        os.rename("tmp/orphan.fa", "_".join([out_prefix, "allstrain_orphan.fa"]))
        os.rename("tmp/all_type.fa", "_".join([out_prefix, "allstrain_all_types.fa"]))
        os.rename("tmp/without_orphan.fa", "_".join([out_prefix, "allstrain_without_orphan.fa"]))

    def _split_fasta_by_strain(self, input_path):
        for fasta in os.listdir(input_path):
            if fasta.endswith(".fa"):
                pre_strain = ""
                num_strain = 0
                with open(os.path.join(input_path, fasta), "r") as f_h:
                    for line in f_h:
                        line = line.strip()
                        if line.startswith(">"):
                            datas = line.split("_")
                            strain = "_".join(datas[2:])
                            if pre_strain != strain:
                                num_strain += 1
                                filename = fasta.split("allstrain")
                                out = open(os.path.join(input_path, 
                                      "".join([filename[0], strain, filename[-1]])), "w")
                                pre_strain = strain
                            out.write(line + "\n")
                        else:
                            out.write(line + "\n")
                if num_strain <= 1:
                    os.remove(os.path.join(input_path, 
                              "".join([filename[0], strain, filename[-1]])))
        out.close()

    def run_MEME(self, meme_path, input_folder, output_folder, input_libs, tsss,
                 fastas, num_motif, widths, parallel, source, wigs, gffs):
        self.multiparser._parser_fasta(fastas)
        self.multiparser._parser_gff(tsss, "TSS")
#        self._check_gff(gffs)
        self._check_gff(tsss)
        tss_path = os.path.join(tsss, "tmp")
        self.multiparser._combine_gff(fastas, tss_path, "fasta", "TSS")    
        self.helper.remove_all_content(input_folder, None, "dir")
        self.helper.check_make_folder(os.getcwd(), "tmp")
        prefixs = []
        for tss in os.listdir(tss_path): ### generate fasta file which is based on the TSS types.
            prefix = tss.replace("_TSS.gff", "")
            prefixs.append(prefix)
            self.helper.check_make_folder(output_folder, prefix)
            out_path = os.path.join(output_folder, prefix)
            self.helper.check_make_folder(input_folder, prefix)
            input_path = os.path.join(input_folder, prefix)
            fasta = self._get_fasta_file(fastas, prefix)
            if source is True:
                print("generating fasta file of {0}".format(prefix))
                Upstream(os.path.join(tss_path, tss), os.path.join(fastas, fasta), 
                         gffs, source, wigs, input_libs, None)
            else:
                if (gffs is False) or (wigs is False) or (input_libs is False):
                    print("Error:please assign proper annotation, tex +/- wig folder and tex treated libs!!!")
                    sys.exit()
                self.multiparser._parser_gff(gffs, None)
                gff_path = os.path.join(gffs, "tmp")
                self.multiparser._combine_gff(fastas, gff_path, "fasta", None)
                if "TSS_class" not in os.listdir(output_folder):
                    os.mkdir(os.path.join(output_folder, "TSS_class"))
                print("classifying TSS and extracting fasta of {0}".format(prefix))
                Upstream(os.path.join(tss_path, tss), os.path.join(fastas, fasta),
                         os.path.join(gff_path, prefix + ".gff"), source, wigs, 
                         input_libs, os.path.join(output_folder, "TSS_class", 
                                                  "_".join([prefix, "TSS.gff"])))
            self._move_and_merge_fasta(input_path, prefix)
            self._split_fasta_by_strain(input_path)
        for prefix in prefixs: ### run MEME
            input_path = os.path.join(input_folder, prefix)
            out_path = os.path.join(output_folder, prefix)
            for fasta in os.listdir(input_path):
                filename = fasta.replace(".fa", "")
                for width in widths:
                    print("Computing promoters of {0} - {1}".format(fasta, width))
                    if "-" in width:
                        self._run_small_motif(meme_path, input_path, out_path, 
                                        filename, width, parallel, num_motif, fasta)
                    else:
                        self._run_normal_motif(meme_path, input_path, out_path,
                                        filename, width, parallel, num_motif, fasta)
        self.helper.remove_tmp(fastas)
        self.helper.remove_tmp(tsss)
        self.helper.remove_tmp(gffs)
        self.helper.remove_tmp(wigs)
        shutil.rmtree("tmp")
