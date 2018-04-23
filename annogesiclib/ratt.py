import os
import sys
import shutil
from subprocess import call, DEVNULL
from annogesiclib.multiparser import Multiparser
from annogesiclib.converter import Converter
from annogesiclib.format_fixer import FormatFixer
from annogesiclib.helper import Helper


class RATT(object):
    '''annotation transfer'''

    def __init__(self, args_ratt):
        self.multiparser = Multiparser()
        self.converter = Converter()
        self.format_fixer = FormatFixer()
        self.helper = Helper()
        if args_ratt.ref_gbk:
            self.gbk = os.path.join(args_ratt.ref_gbk, "gbk_tmp")
            self.gbk_tmp = os.path.join(self.gbk, "tmp")
            self.embl = os.path.join(args_ratt.ref_gbk, "embls")
        if args_ratt.ref_embls:
            self.embl = args_ratt.ref_embls
        self.ratt_log = os.path.join(args_ratt.output_path, "ratt_log.txt")
        self.tmp_files = {"tar": os.path.join(args_ratt.tar_fastas, "tmp"),
                          "ref": os.path.join(args_ratt.ref_fastas, "tmp"),
                          "out_gff": os.path.join(args_ratt.gff_outfolder,
                                                  "tmp"),
                          "gff": os.path.join(args_ratt.gff_outfolder,
                                              "tmp.gff"),
                          "ptt": os.path.join(args_ratt.gff_outfolder,
                                              "tmp.ptt"),
                          "rnt": os.path.join(args_ratt.gff_outfolder,
                                              "tmp.rnt")}

    def _convert_to_pttrnt(self, gffs, files, log):
        for gff in files:
            if gff.endswith(".gff"):
                gff = os.path.join(gffs, gff)
                filename = gff.split("/")
                prefix = filename[-1][:-4]
                rnt = gff[:-3] + "rnt"
                ptt = gff[:-3] + "ptt"
                fasta = self.helper.get_correct_file(self.tmp_files["tar"],
                                                     ".fa", prefix, None, None)
                if fasta:
                    self.converter.convert_gff2rntptt(gff, fasta, ptt, rnt,
                                                      None, None)
                    log.write("\t" + ptt + " is generated.\n")
                    log.write("\t" + rnt + " is generated.\n")

    def _remove_files(self, args_ratt, out_gbk, log):
        self.helper.remove_all_content(args_ratt.gff_outfolder, ".gff", "file")
        self.helper.remove_all_content(args_ratt.gff_outfolder, ".ptt", "file")
        self.helper.remove_all_content(args_ratt.gff_outfolder, ".rnt", "file")
        log.write("Moving the final output files to {0}.\n".format(args_ratt.gff_outfolder))
        self.helper.move_all_content(self.tmp_files["out_gff"],
                                     args_ratt.gff_outfolder, None)
        log.write("Remove the temperary files.\n")
        shutil.rmtree(self.tmp_files["out_gff"])
        shutil.rmtree(self.tmp_files["tar"])
        shutil.rmtree(self.tmp_files["ref"])
        self.helper.remove_tmp_dir(args_ratt.tar_fastas)
        self.helper.remove_tmp_dir(args_ratt.ref_fastas)
        self.helper.remove_tmp_dir(args_ratt.ref_embls)
        self.helper.remove_tmp_dir(args_ratt.ref_gbk)

    def _convert_to_gff(self, ratt_result, args_ratt, files, log):
        name = ratt_result.split(".")
        filename = ".".join(name[1:-2]) + ".gff"
        output_file = os.path.join(args_ratt.output_path, filename)
        self.converter.convert_embl2gff(
             os.path.join(args_ratt.output_path, ratt_result), output_file)
        self.format_fixer.fix_ratt(output_file, ".".join(name[1:-2]),
                                   "tmp_gff")
        shutil.move("tmp_gff", output_file)
        shutil.copy(output_file, os.path.join(args_ratt.gff_outfolder,
                                              filename))
        log.write("\t" + os.path.join(args_ratt.gff_outfolder, filename) + 
                  " is generated.\n")
        files.append(filename)

    def _parser_embl_gbk(self, files):
        self.helper.check_make_folder(self.gbk)
        for file_ in files:
            close = False
            with open(file_, "r") as f_h:
                for line in f_h:
                    if (line.startswith("LOCUS")):
                        out = open(self.gbk_tmp, "w")
                        datas = line.split(" ")
                        for data in datas:
                            if (len(data) != 0) and (data != "LOCUS"):
                                filename = ".".join([data.strip(), "gbk"])
                                break
                    elif (line.startswith("VERSION")):
                        datas = line.split(" ")
                        for data in datas:
                            if (len(data) != 0) and (data != "VERSION"):
                                new_filename = ".".join([data.strip(), "gbk"])
                                break
                        if new_filename.find(filename):
                            filename = new_filename
                    if out:
                        out.write(line)
                    if line.startswith("//"):
                        out.close()
                        close = True
                        shutil.move(self.gbk_tmp,
                                    os.path.join(self.gbk, filename))
            if not close:
                out.close()
        return self.gbk

    def _convert_embl(self, ref_embls, log):
        '''convert gbk to embl'''
        detect_gbk = False
        gbks = []
        out_gbk = None
        for embl in os.listdir(ref_embls):
            if (embl.endswith(".gbk")) or (
                    embl.endswith(".gbff")) or (
                    embl.endswith(".gb")):
                detect_gbk = True
                gbks.append(os.path.join(ref_embls, embl))
        if not detect_gbk:
            log.write("--related_gbk_files is assigned, but not gbk files are detected.\n"
                      "The gbk file names need to be ended at .gbk, .gb, or .gbff. \n")
            print("Error: Please assign proper Genebank files!")
            sys.exit()
        elif detect_gbk:
            out_gbk = self._parser_embl_gbk(gbks)
            log.write("Running converter.py to convert gbk file to embl format.\n")
            self.converter.convert_gbk2embl(out_gbk)
            self.helper.check_make_folder(self.embl)
            self.helper.move_all_content(out_gbk, self.embl, [".embl"])
            log.write("\t" + self.embl + " is generated and the embl files are stored in it.\n")
        return out_gbk

    def _run_ratt(self, args_ratt, tar, ref, out, log):
        if (not os.path.exists(self.embl)) or (
                not os.path.exists(os.path.join(
                    self.tmp_files["tar"], tar + ".fa"))) or (
                not os.path.exists(os.path.join(
                    self.tmp_files["ref"], ref + ".fa"))):
            print("Error: Please check --compare_pair, the strain names "
                  "should be the same as the strain names in fasta, "
                  "genbank or embl files!")
            log.write("The strain names in --compare_pair should be the same "
                      "as the strain names in fasta, genbank, or embl files.\n")
            sys.exit()
        log.write("Make sure your RATT version is at least 1.64.\n")
        log.write("If the RATT can not run properly, please check the "
                  "RATT_HOME and PAGIT_HOME is assigned correctly.\n")
        log.write(" ".join([args_ratt.ratt_path, self.embl,
              os.path.join(self.tmp_files["tar"], tar + ".fa"),
              args_ratt.element, args_ratt.transfer_type,
              os.path.join(self.tmp_files["ref"], ref + ".fa")]) + "\n")
        call([args_ratt.ratt_path, self.embl,
              os.path.join(self.tmp_files["tar"], tar + ".fa"),
              args_ratt.element, args_ratt.transfer_type,
              os.path.join(self.tmp_files["ref"], ref + ".fa")],
             stdout=out, stderr=DEVNULL)
        log.write("Done!\n")

    def _format_and_run(self, args_ratt, log):
        print("Running RATT")
        for pair in args_ratt.pairs:
            ref = pair.split(":")[0]
            tar = pair.split(":")[1]
            out = open(self.ratt_log, "w+")
            self._run_ratt(args_ratt, tar, ref, out, log)
            log.write("The following files are generatd:\n")
            for filename in os.listdir():
                if ("final" in filename):
                    log.write("\t" + filename + "\n")
                    shutil.move(filename, os.path.join(args_ratt.output_path,
                                                       filename))
                elif (args_ratt.element in filename) or (
                      "query" in filename) or (
                      "Reference" in filename) or (
                      "Query" in filename) or (
                      "Sequences" in filename):
                    log.write("\t" + filename + "\n")
                    if os.path.isfile(filename):
                        os.remove(filename)
                    if os.path.isdir(filename):
                        shutil.rmtree(filename)
        out.close()

    def annotation_transfer(self, args_ratt, log):
        self.multiparser.parser_fasta(args_ratt.tar_fastas)
        self.multiparser.parser_fasta(args_ratt.ref_fastas)
        out_gbk = None
        if args_ratt.ref_embls is None:
            out_gbk = self._convert_embl(args_ratt.ref_gbki, log)
        self._format_and_run(args_ratt, log)
        files = []
        for data in os.listdir(args_ratt.output_path):
            if "final.embl" in data:
                log.write("Running converter.py to convert embl "
                          "files in {0} to gff, ptt, and rnt format.\n".format(data))
                self._convert_to_gff(data, args_ratt, files, log)
                self._convert_to_pttrnt(args_ratt.gff_outfolder, files, log)
        self.helper.check_make_folder(self.tmp_files["out_gff"])
        log.write("Merging the output of {0}.\n".format(data))
        for folder in os.listdir(args_ratt.tar_fastas):
            files = []
            if "_folder" in folder:
                datas = folder.split("_folder")
                prefix = ".".join(datas[0].split(".")[:-1])
                for file_ in os.listdir(os.path.join(args_ratt.tar_fastas,
                                                     folder)):
                    files.append(file_[:-3])
                for gff in os.listdir(args_ratt.gff_outfolder):
                    for file_ in files:
                        if (".gff" in gff) and (file_ == gff[:-4]):
                            self.helper.merge_file(os.path.join(
                                 args_ratt.gff_outfolder, gff),
                                 self.tmp_files["gff"])
                        if (".ptt" in gff) and (file_ == gff[:-4]):
                            self.helper.merge_file(os.path.join(
                                 args_ratt.gff_outfolder, gff),
                                 self.tmp_files["ptt"])
                        if (".rnt" in gff) and (file_ == gff[:-4]):
                            self.helper.merge_file(os.path.join(
                                 args_ratt.gff_outfolder, gff),
                                 self.tmp_files["rnt"])
                if os.path.exists(self.tmp_files["gff"]):
                    shutil.move(self.tmp_files["gff"], os.path.join(
                                self.tmp_files["out_gff"], prefix + ".gff"))
                    shutil.move(self.tmp_files["ptt"], os.path.join(
                                self.tmp_files["out_gff"], prefix + ".ptt"))
                    shutil.move(self.tmp_files["rnt"], os.path.join(
                                self.tmp_files["out_gff"], prefix + ".rnt"))
                else:
                    print("Error: Please check your fasta or "
                          "annotation files, they should only contain "
                          "the query genome. And make sure your RATT can "
                          "work properly (check $ANNOgesic/output/"
                          "annotation_transfer/ratt_log.txt).")
                    log.write("Please check your fasta or "
                              "annotation files, they should only contain "
                              "the query genome. And make sure your RATT can "
                              "work properly (check $ANNOgesic/output/"
                              "annotation_transfer/ratt_log.txt).\n")
        self._remove_files(args_ratt, out_gbk, log)
