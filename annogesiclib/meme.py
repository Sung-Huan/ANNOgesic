#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
import shutil
from subprocess import call
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.TSS_upstream import upstream, del_repeat_fasta


class MEME(object):

    def __init__(self, tsss, gffs):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.tss_path = os.path.join(tsss, "tmp")
        self.gff_path = os.path.join(gffs, "tmp")
        self.all_type = "all_type.fa"
        self.tmp_all = "tmp_all.py"
        self.pri = "primary.fa"
        self.sec = "secondary.fa"
        self.inter = "internal.fa"
        self.anti = "antisense.fa"
        self.orph = "orphan.fa"
        self.all_no_orph = "without_orphan.fa"
        self.tmp_fa = "tmp.fa"
        self.tmp_folder = os.path.join(os.getcwd(), "tmp")

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
        if self.all_type in os.listdir(self.tmp_folder):
            os.remove(os.path.join(self.tmp_folder, self.all_type))
        if self.all_no_orph in os.listdir(self.tmp_folder):
            os.remove(os.path.join(self.tmp_folder, self.all_no_orph))
        shutil.copyfile(os.path.join(self.tmp_folder, self.pri), 
                        os.path.join(self.tmp_fa))
        self.helper.merge_file(self.tmp_folder, self.sec, 
                               self.tmp_folder, self.tmp_fa)
        self.helper.merge_file(self.tmp_folder, self.inter, 
                               self.tmp_folder, self.tmp_fa)
        self.helper.merge_file(self.tmp_folder, self.anti, 
                               self.tmp_folder, self.tmp_fa)
        shutil.copyfile(os.path.join(self.tmp_folder, self.tmp_fa), 
                        os.path.join(self.tmp_folder, self.tmp_all))
        self.helper.merge_file(self.tmp_folder, self.orph, 
                               self.tmp_folder, self.tmp_all)
        del_repeat_fasta(os.path.join(self.tmp_folder, self.tmp_fa), 
                         os.path.join(self.tmp_folder, self.all_no_orph))
        del_repeat_fasta(os.path.join(self.tmp_folder, self.tmp_all), 
                         os.path.join(self.tmp_folder, self.all_type))
        os.remove(os.path.join(self.tmp_folder, self.tmp_fa))
        os.remove(os.path.join(self.tmp_folder, self.tmp_all))
        out_prefix = os.path.join(input_path, prefix)
        os.rename(os.path.join(self.tmp_folder, self.pri), 
                  "_".join([out_prefix, "allstrain_primary.fa"]))
        os.rename(os.path.join(self.tmp_folder, self.sec), 
                  "_".join([out_prefix, "allstrain_secondary.fa"]))
        os.rename(os.path.join(self.tmp_folder, self.inter), 
                  "_".join([out_prefix, "allstrain_internal.fa"]))
        os.rename(os.path.join(self.tmp_folder, self.anti), 
                  "_".join([out_prefix, "allstrain_antisense.fa"]))
        os.rename(os.path.join(self.tmp_folder, self.orph), 
                  "_".join([out_prefix, "allstrain_orphan.fa"]))
        os.rename(os.path.join(self.tmp_folder, self.all_type), 
                  "_".join([out_prefix, "allstrain_all_types.fa"]))
        os.rename(os.path.join(self.tmp_folder, self.all_no_orph), 
                  "_".join([out_prefix, "allstrain_without_orphan.fa"]))

    def _split_fasta_by_strain(self, input_path):
        for fasta in os.listdir(input_path):
            if "allstrain" not in fasta:
                os.remove(os.path.join(input_path, fasta))
        for fasta in os.listdir(input_path):
            if fasta.endswith(".fa"):
                pre_strain = ""
                num_strain = 0
                with open(os.path.join(input_path, fasta), "r") as f_h:
                    for line in f_h:
                        line = line.strip()
                        if line.startswith(">"):
                            datas = line.split("_")
                            strain = "_".join(datas[:-2])[1:]
                            if pre_strain != strain:
                                num_strain += 1
                                filename = fasta.split("allstrain")
                                out = open(os.path.join(input_path, 
                                      "".join([filename[0], strain, filename[-1]])), "a")
                                pre_strain = strain
                            out.write(line + "\n")
                        else:
                            out.write(line + "\n")
                if num_strain <= 1:
                    os.remove(os.path.join(input_path, 
                              "".join([filename[0], strain, filename[-1]])))
        out.close()

    def _run_program(self, prefixs, input_folder, output_folder, widths, 
                     meme_path, parallel, num_motif):
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

    def run_meme(self, meme_path, input_folder, output_folder, input_libs, tsss,
                 fastas, num_motif, widths, parallel, source, wigs, gffs):
        self.multiparser._parser_fasta(fastas)
        self.multiparser._parser_gff(tsss, "TSS")
#        self._check_gff(gffs)
#        self._check_gff(tsss)
        self.multiparser._combine_gff(fastas, self.tss_path, "fasta", "TSS")    
        self.helper.remove_all_content(input_folder, None, "dir")
        self.helper.check_make_folder(self.tmp_folder)
        prefixs = []
        for tss in os.listdir(self.tss_path): ### generate fasta file which is based on the TSS types.
            prefix = tss.replace("_TSS.gff", "")
            prefixs.append(prefix)
            self.helper.check_make_folder(os.path.join(output_folder, prefix))
            out_path = os.path.join(output_folder, prefix)
            self.helper.check_make_folder(os.path.join(input_folder, prefix))
            input_path = os.path.join(input_folder, prefix)
            fasta = self._get_fasta_file(fastas, prefix)
            if source is True:
                print("generating fasta file of {0}".format(prefix))
                upstream(os.path.join(self.tss_path, tss), os.path.join(fastas, fasta), 
                         gffs, source, wigs, input_libs, None)
            else:
                if (gffs is False) or (wigs is False) or (input_libs is False):
                    print("Error:please assign proper annotation, tex +/- wig folder and tex treated libs!!!")
                    sys.exit()
                self.multiparser._parser_gff(gffs, None)
                self.multiparser._combine_gff(fastas, self.gff_path, "fasta", None)
                if "TSS_class" not in os.listdir(output_folder):
                    os.mkdir(os.path.join(output_folder, "TSS_class"))
                print("classifying TSS and extracting fasta of {0}".format(prefix))
                print(os.path.join(self.tss_path, tss))
                upstream(os.path.join(self.tss_path, tss), os.path.join(fastas, fasta),
                         os.path.join(self.gff_path, prefix + ".gff"), source, wigs, 
                         input_libs, os.path.join(output_folder, "TSS_class", 
                                                  "_".join([prefix, "TSS.gff"])))
            self._move_and_merge_fasta(input_path, prefix)
            self._split_fasta_by_strain(input_path)
        self._run_program(prefixs, input_folder, output_folder, widths, 
                     meme_path, parallel, num_motif)
        self.helper.remove_tmp(fastas)
        self.helper.remove_tmp(tsss)
        self.helper.remove_tmp(gffs)
        self.helper.remove_tmp(wigs)
        shutil.rmtree("tmp")
