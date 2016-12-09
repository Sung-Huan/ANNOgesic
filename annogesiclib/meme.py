import os
import sys
import shutil
from subprocess import call
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.TSS_upstream import upstream, del_repeat_fasta
from annogesiclib.gen_promoter_table import gen_promoter_table


class MEME(object):
    '''detection of promoter'''

    def __init__(self, args_pro):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.tss_path = os.path.join(args_pro.tsss, "tmp")
        if args_pro.gffs is not None:
            self.gff_path = os.path.join(args_pro.gffs, "tmp")
        else:
            self.gff_path = None
        self.out_fasta = os.path.join(args_pro.output_folder, "fasta_class")
        self.tmp_folder = os.path.join(os.getcwd(), "tmp")
        self.fastas = {"pri": os.path.join(self.tmp_folder, "primary.fa"),
                       "sec": os.path.join(self.tmp_folder, "secondary.fa"),
                       "inter": os.path.join(self.tmp_folder, "internal.fa"),
                       "anti": os.path.join(self.tmp_folder, "antisense.fa"),
                       "orph": os.path.join(self.tmp_folder, "orphan.fa"),
                       "all_no_orph": "without_orphan.fa",
                       "all": "all_type.fa",
                       "tmp_fa": os.path.join(self.tmp_folder, "tmp.fa"),
                       "tmp_all": os.path.join(self.tmp_folder, "tmp_all.fa")}
        self.all_fasta = os.path.join(args_pro.fastas, "allfasta.fa")
        self.all_tss = os.path.join(self.tss_path, "allfasta_TSS.gff")

    def _gen_and_check_folder(self, out_path, folder, type_):
        sub_out_folder = os.path.join(out_path, type_)
        if folder in os.listdir(sub_out_folder):
            shutil.rmtree(os.path.join(sub_out_folder, folder))
        return sub_out_folder

    def _run_normal_motif(self, input_path, out_path, filename,
                          fasta, width, args_pro):
        '''run MEME with specific width'''
        folder = "_".join(["promoter_motifs", filename,
                           str(width), "nt"])
        if (args_pro.program.lower() == "meme") or (
                args_pro.program.lower() == "both"):
            meme_folder = self._gen_and_check_folder(
                            out_path, folder, "MEME")
            command = [args_pro.meme_path, "-maxsize", "1000000",
                       "-dna", "-nmotifs", str(args_pro.num_motif),
                       "-w", str(width), "-maxiter", "100",
                       "-evt", str(args_pro.e_value)]
            if args_pro.para is not None:
                command = command + ["-p", args_pro.para]
            call(command + ["-oc", os.path.join(meme_folder, folder),
                            os.path.join(input_path, fasta)])
        if (args_pro.program.lower() == "glam2") or (
                args_pro.program.lower() == "both"):
            glam_folder = self._gen_and_check_folder(
                            out_path, folder, "GLAM2")
            call([args_pro.glam2_path,
                  "-O", os.path.join(glam_folder, folder), "-w",
                  str(width), "-b", str(width), "-r",
                  str(args_pro.num_motif), "-n", str(args_pro.end_run),
                  "n", os.path.join(input_path, fasta)])

    def _run_small_motif(self, input_path, out_path, filename,
                         fasta, width, args_pro):
        '''run MEME with range of width'''
        data = width.split("-")
        min_width = data[0]
        max_width = data[1]
        folder = "_".join(["promoter_motifs", filename,
                           "-".join([str(min_width), str(max_width)]), "nt"])
        if (args_pro.program.lower() == "meme") or (
                args_pro.program.lower() == "both"):
            meme_folder = self._gen_and_check_folder(
                            out_path, folder, "MEME")
            command = [args_pro.meme_path, "-maxsize", "1000000",
                       "-dna", "-nmotifs", str(args_pro.num_motif),
                       "-minsites", "0", "-maxsites", "2",
                       "-minw", str(min_width), "-maxw", str(max_width),
                       "-maxiter", "100",
                       "-evt", str(args_pro.e_value)]
            if args_pro.para is not None:
                command = command + ["-p", args_pro.para]
            call(command + ["-oc", os.path.join(meme_folder, folder),
                            os.path.join(input_path, fasta)])
        elif (args_pro.program.lower() == "glam2") or (
                args_pro.program.lower() == "both"):
            glam_folder = self._gen_and_check_folder(
                            out_path, folder, "GLAM2")
            call([args_pro.glam2_path,
                  "-O", os.path.join(glam_folder, folder), "-a",
                  str(min_width), "-b", str(max_width), "-r",
                  str(args_pro.num_motif), "-n", str(args_pro.end_run),
                  "n", os.path.join(input_path, fasta)])

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
        shutil.move(self.fastas["pri"], "_".join([
            out_prefix, "allstrain_primary.fa"]))
        shutil.move(self.fastas["sec"], "_".join([
            out_prefix, "allstrain_secondary.fa"]))
        shutil.move(self.fastas["inter"], "_".join([
            out_prefix, "allstrain_internal.fa"]))
        shutil.move(self.fastas["anti"], "_".join([
            out_prefix, "allstrain_antisense.fa"]))
        shutil.move(self.fastas["orph"], "_".join([
            out_prefix, "allstrain_orphan.fa"]))
        shutil.move(all_type, "_".join([
            out_prefix, "allstrain_all_types.fa"]))
        shutil.move(all_no_orph, "_".join([
            out_prefix, "allstrain_without_orphan.fa"]))

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
                            strain = "_".join(datas[2:])
                            if pre_strain != strain:
                                num_strain += 1
                                filename = fasta.split("allstrain")
                                if out is not None:
                                    out.close()
                                out = open(os.path.join(
                                           input_path, "".join([
                                               filename[0], strain,
                                               filename[-1]])), "a")
                                pre_strain = strain
                            out.write(line + "\n")
                        else:
                            out.write(line + "\n")
                if num_strain <= 1:
                    os.remove(os.path.join(input_path,
                              "".join([filename[0], strain, filename[-1]])))
        out.close()

    def _run_program(self, prefixs, args_pro):
        for prefix in prefixs:
            input_path = os.path.join(self.out_fasta, prefix)
            out_path = os.path.join(args_pro.output_folder, prefix)
            if args_pro.program.lower() == "both":
                self.helper.check_make_folder(os.path.join(out_path, "MEME"))
                self.helper.check_make_folder(os.path.join(out_path, "GLAM2"))
            elif args_pro.program.lower() == "meme":
                self.helper.check_make_folder(os.path.join(out_path, "MEME"))
            elif args_pro.program.lower() == "glam2":
                self.helper.check_make_folder(os.path.join(out_path, "GLAM2"))
            for fasta in os.listdir(input_path):
                filename = fasta.replace(".fa", "")
                for width in args_pro.widths:
                    print("Computing promoters of {0} - {1}".format(
                          fasta, width))
                    if "-" in width:
                        self._run_small_motif(input_path, out_path, filename,
                                              fasta, width, args_pro)
                    else:
                        self._run_normal_motif(input_path, out_path, filename,
                                               fasta, width, args_pro)

    def _combine_file(self, prefixs, args_pro):
        '''combine all TSS file in the input folder to generate the 
        global TSS for detecting the global promoter'''
        if args_pro.source:
            for tss in os.listdir(self.tss_path):
                if tss.endswith("_TSS.gff"):
                    self.helper.merge_file(os.path.join(
                         self.tss_path, tss), self.all_tss)
            for fasta in os.listdir(args_pro.fastas):
                if (fasta.endswith(".fa")) or (
                        fasta.endswith(".fna")) or (
                        fasta.endswith(".fasta")):
                    self.helper.merge_file(os.path.join(
                         args_pro.fastas, fasta), self.all_fasta)
        else:
            for tss in os.listdir(os.path.join(
                                  args_pro.output_folder, "TSS_class")):
                if tss.endswith("_TSS.gff"):
                    self.helper.merge_file(os.path.join(
                         self.tss_path, tss), self.all_tss)
            for fasta in os.listdir(args_pro.fastas):
                if (fasta.endswith(".fa")) or (
                        fasta.endswith(".fna")) or (
                        fasta.endswith(".fasta")):
                    self.helper.merge_file(os.path.join(
                         args_pro.fastas, fasta), self.all_fasta)
        print("Generating fasta file of all fasta files")
        prefixs.append("allfasta")
        input_path = os.path.join(self.out_fasta, "allfasta")
        self.helper.check_make_folder(os.path.join(
                                      args_pro.output_folder, "allfasta"))
        self.helper.check_make_folder(os.path.join(
                                      self.out_fasta, "allfasta"))
        args_pro.source = True
        upstream(self.all_tss, self.all_fasta, None,
                 None, args_pro, None)
        self._move_and_merge_fasta(input_path, "allfasta")

    def _remove_files(self, args_pro):
        self.helper.remove_tmp_dir(args_pro.fastas)
        self.helper.remove_tmp_dir(args_pro.tsss)
        self.helper.remove_tmp_dir(args_pro.gffs)
        if "tmp_wig" in os.listdir(args_pro.output_folder):
            shutil.rmtree(os.path.join(args_pro.output_folder, "tmp_wig"))
        if "allfasta" in os.listdir(os.getcwd()):
            shutil.rmtree("allfasta")
        shutil.rmtree("tmp")

    def _gen_table(self, output_folder, prefixs, combine, program):
        '''generate the promoter table'''
        if combine:
            strains = prefixs + ["allfasta"]
        else:
            strains = prefixs
        for strain in strains:
            tss_file = os.path.join(self.tss_path, strain + "_TSS.gff")
            if (program.lower() == "both") or (
                    program.lower() == "meme"):
                for folder in os.listdir(os.path.join(output_folder,
                                                      strain, "MEME")):
                    gen_promoter_table(os.path.join(output_folder, strain,
                                       "MEME", folder, "meme.txt"),
                                       os.path.join(output_folder, strain,
                                       "MEME", folder, "meme.csv"), tss_file,
                                       "meme")
            if (program.lower() == "both") or (
                    program.lower() == "glam2"):
                for folder in os.listdir(os.path.join(output_folder,
                                                      strain, "GLAM2")):
                    gen_promoter_table(os.path.join(output_folder, strain,
                                        "GLAM2", folder, "glam2.txt"),
                                        os.path.join(output_folder, strain,
                                        "GLAM2", folder, "glam2.csv"),
                                        tss_file, "glam2")

    def _get_upstream(self, args_pro, prefix, tss, fasta):
        '''get upstream sequence of TSS'''
        if args_pro.source:
            print("Generating fasta file of {0}".format(prefix))
            upstream(os.path.join(self.tss_path, tss),
                     os.path.join(args_pro.fastas, fasta),
                     None, None, args_pro, prefix)
        else:
            if (args_pro.gffs is None) or (
                    args_pro.tex_wigs is None) or (
                    args_pro.input_libs is None):
                print("Error: Please assign proper annotation, tex +/- "
                      "wig folder and tex treated libs!!!")
                sys.exit()
            if "TSS_class" not in os.listdir(args_pro.output_folder):
                os.mkdir(os.path.join(args_pro.output_folder, "TSS_class"))
            
            print("Classifying TSS and extracting fasta {0}".format(prefix))
            upstream(os.path.join(self.tss_path, tss),
                     os.path.join(args_pro.fastas, fasta),
                     os.path.join(self.gff_path, prefix + ".gff"),
                     os.path.join(args_pro.output_folder, "TSS_class",
                     "_".join([prefix, "TSS.gff"])), args_pro, prefix)

    def run_meme(self, args_pro):
        if "allfasta.fa" in os.listdir(args_pro.fastas):
            os.remove(self.all_fasta)
            if "allfasta.fa_folder" in os.listdir(args_pro.fastas):
                shutil.rmtree(os.path.join(args_pro.fastas,
                              "allfasta.fa_folder"))
        self.multiparser.parser_fasta(args_pro.fastas)
        self.multiparser.parser_gff(args_pro.tsss, "TSS")
        if "allfasta_TSS.gff" in os.listdir(self.tss_path):
            os.remove(self.all_tss)
        if args_pro.gffs is not None:
            self._check_gff(args_pro.gffs)
            self.multiparser.parser_gff(args_pro.gffs, None)
            self.multiparser.combine_gff(args_pro.fastas, self.gff_path,
                                         "fasta", None)
        self._check_gff(args_pro.tsss)
        self.multiparser.combine_gff(args_pro.fastas, self.tss_path,
                                     "fasta", "TSS")
        self.helper.check_make_folder(self.out_fasta)
        self.helper.check_make_folder(self.tmp_folder)
        prefixs = []
        for tss in os.listdir(self.tss_path):
            prefix = tss.replace("_TSS.gff", "")
            prefixs.append(prefix)
            self.helper.check_make_folder(os.path.join(args_pro.output_folder,
                                                       prefix))
            self.helper.check_make_folder(os.path.join(self.out_fasta,
                                                       prefix))
            input_path = os.path.join(self.out_fasta, prefix)
            fasta = self._get_fasta_file(args_pro.fastas, prefix)
            self._get_upstream(args_pro, prefix, tss, fasta)
            self._move_and_merge_fasta(input_path, prefix)
            self._split_fasta_by_strain(input_path)
        if args_pro.combine:
            self._combine_file(prefixs, args_pro)
        self._run_program(prefixs, args_pro)
        print("Generating the table")
        self._gen_table(args_pro.output_folder, prefixs,
                        args_pro.combine, args_pro.program)
        self._remove_files(args_pro)
