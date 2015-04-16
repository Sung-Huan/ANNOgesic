#!/usr/bin/python
import os	
import sys
import shutil
from subprocess import call, DEVNULL
import csv
from annogesiclib.multiparser import Multiparser
from annogesiclib.converter import Converter
from annogesiclib.format_fixer import Format_Fixer
from annogesiclib.helper import Helper

class RATT(object):

    def __init__(self, ref_embls, output_path, tar_fastas, gff_outfolder):
        self.multiparser = Multiparser()
        self.converter = Converter()   
        self.format_fixer = Format_Fixer()
        self.helper = Helper()
        self.gbk = os.path.join(ref_embls, "gbk_tmp")
        self.gbk_tmp = os.path.join(self.gbk, "tmp")
        self.embl = os.path.join(ref_embls, "embls")
        self.ratt_log = os.path.join(output_path, "ratt_log.txt")
        self.tar_tmp = os.path.join(tar_fastas, "tmp")
        self.out_gff_tmp = os.path.join(gff_outfolder, "tmp")

    def _convert_to_pttrnt(self, gffs, fastas):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                gff = os.path.join(gffs, gff)
                filename = gff.split("/")
                prefix = filename[-1][:-4]
                rnt = gff[:-3] + "rnt"
                ptt = gff[:-3] + "ptt"
                fasta = self.helper.get_correct_file(self.tar_tmp, ".fa", prefix, None)
                if fasta:
                    self.converter.convert_gff2rntptt(gff, fasta, ptt, rnt, None, None)

    def _remove_files(self, gff_outfolder, tar_fastas, out_gbk):
        self.helper.remove_all_content(gff_outfolder, ".gff", "file")
        self.helper.remove_all_content(gff_outfolder, ".ptt", "file")
        self.helper.remove_all_content(gff_outfolder, ".rnt", "file")
        self.helper.move_all_content(self.out_gff_tmp, gff_outfolder, None)
        shutil.rmtree(self.out_gff_tmp)
        shutil.rmtree(self.tar_tmp)
        shutil.rmtree(self.embl)
        self.helper.remove_all_content(tar_fastas, "_folder", "dir")
        if out_gbk:
            shutil.rmtree(out_gbk)

    def _convert_to_gff(self, ratt_result, output_path, gff_outfolder, combine):
        name = ratt_result.split(".")
        filename = name[1] + ".gff"
        output_file = os.path.join(output_path, filename)
        self.converter.convert_embl2gff(os.path.join(output_path, ratt_result), 
                                        output_file)
        self.format_fixer.fix_ratt(output_file, name[1], "tmp_gff")            
        os.rename("tmp_gff", output_file)
        shutil.copy(output_file, os.path.join(gff_outfolder, filename))

    def _move_annotation(self, ref_folder, ref, tar_folder, tar):
        os.rename(os.path.join(ref_folder, ref), 
                  os.path.join(tar_folder, tar))

    def _parser_embl_gbk(self, files, ref_embls):
        self.helper.check_make_folder(self.gbk)
        for file_ in files:
            with open(file_, "r") as f_h:
                for line in f_h:
                    if (line.startswith("LOCUS")):
                        out = open(self.gbk_tmp, "w")
                        datas = line.split(" ")
                        for data in datas:
                            if (len(data) != 0) and (data != "LOCUS"):
                                filename = ".".join([data, "gbk"])
                                break
                    elif (line.startswith("VERSION")):
                        datas = line.split(" ")
                        for data in datas:
                            if (len(data) != 0) and (data != "VERSION"):
                                filename = ".".join([data, "gbk"])
                                break
                    if out:
                        out.write(line)
                    if line.startswith("//"):
                       out.close()
                       os.rename(self.gbk_tmp, 
                                 os.path.join(self.gbk, filename))
        return self.gbk

    def _convert_embl(self, ref_embls):
        detect_gbk = False
        embls = []
        gbks = []
        out_gbk = None
        for embl in os.listdir(ref_embls):
            if embl.endswith(".gbk"):
                detect_gbk = True
                gbks.append(os.path.join(ref_embls, embl))
        if (detect_gbk is False):
            print("Error: please assign proper folder for Genebank file!!!")
            sys.exit()
        elif detect_gbk:
            out_gbk = self._parser_embl_gbk(gbks, ref_embls)
            self.converter.convert_gbk2embl(out_gbk)
            self.helper.check_make_folder(self.embl)
            self.helper.move_all_content(out_gbk, self.embl, ".embl")
        return out_gbk

    def _run_ratt(self, output_path, tar_fastas, ratt_path,
                  element, transfer_type, ref_fastas):
        print("Running RATT...")
        out = open(self.ratt_log, "w+")
        for target in os.listdir(self.tar_tmp):
            print(target)
            call([os.path.join(ratt_path, "start.ratt.sh"),
                  self.embl,
                  os.path.join(self.tar_tmp, target),
                  element, transfer_type, ref_fastas],
                  stdout=out, stderr=DEVNULL)
            for filename in os.listdir():
                if ("final" in filename):
                    shutil.move(filename, os.path.join(output_path, filename))
                elif (element in filename) or \
                   ("query" in filename) or \
                   ("Reference" in filename) or \
                   ("Query" in filename) or \
                   ("Sequences" in filename):
                    if os.path.isfile(filename):
                        os.remove(filename)
                    if os.path.isdir(filename):
                        shutil.rmtree(filename)

    def annotation_transfer(self, ratt_path, pagit_folder, element, transfer_type, ref_embls, 
                            tar_fastas, ref_fastas, output_path, convert, combine, gff_outfolder):
        self.multiparser._parser_fasta(tar_fastas)
        first = True
        out_gbk = self._convert_embl(ref_embls)
        self._run_ratt(output_path, tar_fastas, ratt_path,
                       element, transfer_type, ref_fastas)
        if convert is True:
            for data in os.listdir(output_path):
                if "final.embl" in data:
                    self._convert_to_gff(data, output_path, gff_outfolder, combine)
                    self._convert_to_pttrnt(gff_outfolder, tar_fastas)
        self.helper.check_make_folder(self.out_gff_tmp)
        for folder in os.listdir(tar_fastas):
            files = []
            if "_folder" in folder:
                datas = folder.split("_folder")
                prefix = datas[0][:-3]
                for file_ in os.listdir(os.path.join(tar_fastas, folder)):
                    files.append(file_[:-3])
                for gff in os.listdir(gff_outfolder):
                    for file_ in files:
                        if (".gff" in gff) and (file_ == gff[:-4]):
                            self.helper.merge_file(gff_outfolder, gff, gff_outfolder, "tmp.gff")
                        if (".ptt" in gff) and (file_ == gff[:-4]):
                            self.helper.merge_file(gff_outfolder, gff, gff_outfolder, "tmp.ptt")
                        if (".rnt" in gff) and (file_ == gff[:-4]):
                            self.helper.merge_file(gff_outfolder, gff, gff_outfolder, "tmp.rnt")
                self._move_annotation(gff_outfolder, "tmp.gff", self.out_gff_tmp, prefix + ".gff")
                self._move_annotation(gff_outfolder, "tmp.ptt", self.out_gff_tmp, prefix + ".ptt")
                self._move_annotation(gff_outfolder, "tmp.rnt", self.out_gff_tmp, prefix + ".rnt")
        self._remove_files(gff_outfolder, tar_fastas, out_gbk) 
