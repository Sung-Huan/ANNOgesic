#!/usr/bin/python
import os	
import sys
import shutil
from subprocess import call, DEVNULL
import csv
from transaplib.multiparser import Multiparser
from transaplib.converter import Converter
from transaplib.format_fixer import Format_Fixer
from transaplib.helper import Helper

class RATT(object):

    def __init__(self):
        self.multiparser = Multiparser()
        self.converter = Converter()   
        self.format_fixer = Format_Fixer()
        self.helper = Helper()

    def _convert_to_pttrnt(self, gffs, fastas):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                gff = os.path.join(gffs, gff)
                filename = gff.split("/")
                prefix = filename[-1][:-4]
                rnt = gff[:-3] + "rnt"
                ptt = gff[:-3] + "ptt"
                fasta = self.helper.get_correct_file(os.path.join(fastas, "tmp"), ".fa", prefix, None)
                if fasta:
                    self.converter.Convert_gff2rntptt(gff, fasta, ptt, rnt, None, None)

    def _remove_files(self, gff_outfolder, ref_fastas, tar_fastas):
        self.helper.remove_all_content(gff_outfolder, ".gff", "file")
        self.helper.remove_all_content(gff_outfolder, ".ptt", "file")
        self.helper.remove_all_content(gff_outfolder, ".rnt", "file")
        self.helper.move_all_content(os.path.join(gff_outfolder, "tmp"), gff_outfolder, None)
        shutil.rmtree(os.path.join(gff_outfolder, "tmp"))
        shutil.rmtree(os.path.join(ref_fastas, "tmp"))
        self.helper.remove_all_content(ref_fastas, "_folder", "dir")
        shutil.rmtree(os.path.join(tar_fastas, "tmp"))
        self.helper.remove_all_content(tar_fastas, "_folder", "dir")

    def _convert_to_gff(self, ratt_result, output_path, gff_outfolder, combine):
        name = ratt_result.split(".")
        filename = name[1] + ".gff"
        output_file = os.path.join(output_path, filename)
        self.converter.Convert_embl2gff(os.path.join(output_path, ratt_result), 
                                        output_file)
        self.format_fixer.Fix_ratt(output_file, name[1], "tmp_gff")            
        os.rename("tmp_gff", output_file)
        shutil.copy (output_file, os.path.join(gff_outfolder, filename))

    def _move_annotation(self, ref_folder, ref, tar_folder, tar):
        os.rename(os.path.join(ref_folder, ref), 
                  os.path.join(tar_folder, tar))

    def run_RATT(self, ratt_path, pagit_folder,  compares, element, transfer_type, ref_embls, 
                 tar_fastas, ref_fastas, output_path, convert, combine, gff_outfolder):
        self.multiparser._parser_fasta(ref_fastas)
        self.multiparser._parser_fasta(tar_fastas)
        first = True
        datas = []
        for compare in compares:
            files = compare.split(":")
            datas.append({"target_id": files[0], "ref_id": files[1]})
        print("Running RATT...")
        for data in datas:
            out = open(os.path.join(output_path, "ratt_log.txt"), "w+")
            call([os.path.join(ratt_path, "start.ratt.sh"),
                  os.path.join(ref_embls, data["ref_id"] + "_embl"),
                  os.path.join(tar_fastas, "tmp", data["target_id"] + ".fa"),
                  element, transfer_type,
                  os.path.join(ref_fastas, "tmp", data["ref_id"] + ".fa")],
                  stdout=out, stderr=DEVNULL)
            out.close()
        for filename in os.listdir():
            if (element in filename) or \
               ("query" in filename) or \
               ("Reference" in filename) or \
               ("Query" in filename) or \
               ("Sequences" in filename):
                shutil.move(filename, os.path.join(output_path, filename))
        if convert is True:
            for data in datas:
                ratt_result = ".".join([element, data["target_id"], "final.embl"])
                self._convert_to_gff(ratt_result, output_path, gff_outfolder, combine)
                self._convert_to_pttrnt(gff_outfolder, tar_fastas)
        self.helper.check_make_folder(gff_outfolder, "tmp")
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
                tar_folder = os.path.join(gff_outfolder, "tmp")
                self._move_annotation(gff_outfolder, "tmp.gff", tar_folder, prefix + ".gff")
                self._move_annotation(gff_outfolder, "tmp.ptt", tar_folder, prefix + ".ptt")
                self._move_annotation(gff_outfolder, "tmp.rnt", tar_folder, prefix + ".rnt")
        self._remove_files(gff_outfolder, ref_fastas, tar_fastas) 
