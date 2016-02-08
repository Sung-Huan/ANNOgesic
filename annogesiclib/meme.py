import os
import sys
import shutil
from subprocess import call
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.TSS_upstream import upstream, del_repeat_fasta


class MEME(object):

    def __init__(self, tsss, gffs, fastas):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.tss_path = os.path.join(tsss, "tmp")
        if gffs is not None:
            self.gff_path = os.path.join(gffs, "tmp")
        else:
            self.gff_path = None
        self.tmp_folder = os.path.join(os.getcwd(), "tmp")
        self.fastas = {"pri": os.path.join(self.tmp_folder, "primary.fa"),
                       "sec": os.path.join(self.tmp_folder, "secondary.fa"),
                       "inter": os.path.join(self.tmp_folder, "internal.fa"),
                       "anti": os.path.join(self.tmp_folder, "antisense.fa"),
                       "orph": os.path.join(self.tmp_folder, "orphan.fa"),
                       "all_no_orph": "without_orphan.fa", "all": "all_type.fa",
                       "tmp_fa": os.path.join(self.tmp_folder, "tmp.fa"),
                       "tmp_all": os.path.join(self.tmp_folder, "tmp_all.fa")}
        self.all_fasta = os.path.join(fastas, "allfasta.fa")
        self.all_tss = os.path.join(self.tss_path, "allfasta_TSS.gff")

    def _run_normal_motif(self, meme_path, input_path, out_path,
                          filename, width, num_motif, fasta):
        folder = "_".join(["promoter_motifs", filename, str(width), "nt"])
        if folder not in os.listdir(out_path):
            call([meme_path, "-maxsize", "1000000",
                  "-dna", "-nmotifs", str(num_motif),
                  "-w", str(width), "-maxiter", "100",
                  "-oc", os.path.join(out_path, folder),
                  os.path.join(input_path, fasta)])

    def _run_small_motif(self, meme_path, input_path, out_path,
                         filename, width, num_motif, fasta):
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
                  "-maxiter", "100",
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
        all_type = os.path.join(self.tmp_folder, self.fastas["all"])
        all_no_orph = os.path.join(self.tmp_folder, self.fastas["all_no_orph"])
        if self.fastas["all"] in os.listdir(self.tmp_folder):
            os.remove(all_type)
        if self.fastas["all_no_orph"] in os.listdir(self.tmp_folder):
            os.remove(all_no_orph)
        shutil.copyfile(self.fastas["pri"], self.fastas["tmp_fa"])
        self.helper.merge_file(self.fastas["sec"], self.fastas["tmp_fa"])
        self.helper.merge_file(self.fastas["inter"], self.fastas["tmp_fa"])
        self.helper.merge_file(self.fastas["anti"], self.fastas["tmp_fa"])
        shutil.copyfile(self.fastas["tmp_fa"], self.fastas["tmp_all"])
        self.helper.merge_file(self.fastas["orph"], self.fastas["tmp_all"])
        del_repeat_fasta(self.fastas["tmp_fa"], all_no_orph)
        del_repeat_fasta(self.fastas["tmp_all"], all_type)
        os.remove(self.fastas["tmp_fa"])
        os.remove(self.fastas["tmp_all"])
        out_prefix = os.path.join(input_path, prefix)
        shutil.move(self.fastas["pri"], "_".join([out_prefix,
                                      "allstrain_primary.fa"]))
        shutil.move(self.fastas["sec"], "_".join([out_prefix,
                                      "allstrain_secondary.fa"]))
        shutil.move(self.fastas["inter"], "_".join([out_prefix,
                                        "allstrain_internal.fa"]))
        shutil.move(self.fastas["anti"], "_".join([out_prefix,
                                       "allstrain_antisense.fa"]))
        shutil.move(self.fastas["orph"], "_".join([out_prefix,
                                       "allstrain_orphan.fa"]))
        shutil.move(all_type, "_".join([out_prefix,
                            "allstrain_all_types.fa"]))
        shutil.move(all_no_orph, "_".join([out_prefix,
                               "allstrain_without_orphan.fa"]))

    def _split_fasta_by_strain(self, input_path):
        for fasta in os.listdir(input_path):
            if "allstrain" not in fasta:
                os.remove(os.path.join(input_path, fasta))
        out = None
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
                                if out is not None:
                                    out.close()
                                out = open(os.path.join(input_path,
                                      "".join([filename[0], strain,
                                               filename[-1]])), "a")
                                pre_strain = strain
                            out.write(line + "\n")
                        else:
                            out.write(line + "\n")
                if num_strain <= 1:
                    os.remove(os.path.join(input_path,
                              "".join([filename[0], strain, filename[-1]])))
        out.close()

    def _run_program(self, prefixs, input_folder, output_folder, widths,
                     meme_path, num_motif):
        for prefix in prefixs: ### run MEME
            input_path = os.path.join(input_folder, prefix)
            out_path = os.path.join(output_folder, prefix)
            for fasta in os.listdir(input_path):
                filename = fasta.replace(".fa", "")
                for width in widths:
                    print("Computing promoters of {0} - {1}".format(
                          fasta, width))
                    if "-" in width:
                        self._run_small_motif(meme_path, input_path, out_path,
                                    filename, width, num_motif, fasta)
                    else:
                        self._run_normal_motif(meme_path, input_path, out_path,
                                    filename, width, num_motif, fasta)

    def _combine_file(self, source, gffs, fastas, wigs, input_libs,
                      input_folder, output_folder, prefixs, nt_before):
        if source:
            for tss in os.listdir(self.tss_path):
                if tss.endswith("_TSS.gff"):
                    self.helper.merge_file(os.path.join(
                                           self.tss_path, tss), self.all_tss)
            for fasta in os.listdir(fastas):
                if (fasta.endswith(".fa")) or (
                    fasta.endswith(".fna")) or (
                    fasta.endswith(".fasta")):
                    self.helper.merge_file(os.path.join(
                                           fastas, fasta), self.all_fasta)
        else:
            for tss in os.listdir(os.path.join(output_folder, "TSS_class")):
                if tss.endswith("_TSS.gff"):
                    self.helper.merge_file(os.path.join(
                                           self.tss_path, tss), self.all_tss)
            for fasta in os.listdir(fastas):
                if (fasta.endswith(".fa")) or \
                   (fasta.endswith(".fna")) or \
                   (fasta.endswith(".fasta")):
                    self.helper.merge_file(os.path.join(
                                           fastas, fasta), self.all_fasta)
        print("generating fasta file of all fasta files")
        prefixs.append("allfasta")
        input_path = os.path.join(input_folder, "allfasta")
        self.helper.check_make_folder(os.path.join(output_folder, "allfasta"))
        self.helper.check_make_folder(os.path.join(input_folder, "allfasta"))
        upstream(self.all_tss, self.all_fasta, gffs,
                 True, wigs, input_libs, None, nt_before)
        self._move_and_merge_fasta(input_path, "allfasta")

    def _remove_files(self, fastas, tsss, gffs, wigs):
        self.helper.remove_tmp(fastas)
        self.helper.remove_tmp(tsss)
        self.helper.remove_tmp(gffs)
        self.helper.remove_tmp(wigs)
        if "allfasta.fa" in os.listdir(fastas):
            os.remove(self.all_fasta)
        if "allfasta" in os.listdir(os.getcwd()):
            shutil.rmtree("allfasta")
        shutil.rmtree("tmp")

    def run_meme(self, meme_path, input_folder, output_folder, input_libs,
                 tsss, fastas, num_motif, nt_before, widths, source, wigs,
                 gffs, combine):
        if "allfasta.fa" in os.listdir(fastas):
            os.remove(self.all_fasta)
            if "allfasta.fa_folder" in os.listdir(fastas):
                shutil.rmtree(os.path.join(fastas, "allfasta.fa_folder"))
        self.multiparser.parser_fasta(fastas)
        self.multiparser.parser_gff(tsss, "TSS")
        if "allfasta_TSS.gff" in os.listdir(self.tss_path):
            os.remove(self.all_tss)
        if gffs is not None:
            self._check_gff(gffs)
        self._check_gff(tsss)
        self.multiparser.combine_gff(fastas, self.tss_path, "fasta", "TSS")
        self.helper.remove_all_content(input_folder, None, "dir")
        self.helper.check_make_folder(self.tmp_folder)
        prefixs = []
        for tss in os.listdir(self.tss_path):
            print(tss)
            prefix = tss.replace("_TSS.gff", "")
            prefixs.append(prefix)
            self.helper.check_make_folder(os.path.join(output_folder, prefix))
            self.helper.check_make_folder(os.path.join(input_folder, prefix))
            input_path = os.path.join(input_folder, prefix)
            fasta = self._get_fasta_file(fastas, prefix)
            if source:
                print("generating fasta file of {0}".format(prefix))
                upstream(os.path.join(self.tss_path, tss),
                         os.path.join(fastas, fasta), gffs, source,
                         wigs, input_libs, None, nt_before)
            else:
                if (gffs is None) or (wigs is None) or (input_libs is None):
                    print("Error:please assign proper annotation, tex +/- wig folder and tex treated libs!!!")
                    sys.exit()
                self.multiparser.parser_gff(gffs, None)
                self.multiparser.combine_gff(fastas, self.gff_path,
                                              "fasta", None)
                if "TSS_class" not in os.listdir(output_folder):
                    os.mkdir(os.path.join(output_folder, "TSS_class"))
                print("classifying TSS and extracting fasta {0}".format(prefix))
                upstream(os.path.join(self.tss_path, tss),
                         os.path.join(fastas, fasta),
                         os.path.join(self.gff_path, prefix + ".gff"),
                         source, wigs, input_libs,
                         os.path.join(output_folder, "TSS_class",
                         "_".join([prefix, "TSS.gff"])), nt_before)
            self._move_and_merge_fasta(input_path, prefix)
            self._split_fasta_by_strain(input_path)
        if combine:
            self._combine_file(source, gffs, fastas, wigs, input_libs,
                               input_folder, output_folder, prefixs, nt_before)
        self._run_program(prefixs, input_folder, output_folder, widths,
                          meme_path, num_motif)
        self._remove_files(fastas, tsss, gffs, wigs)
