#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
from subprocess import call
import multiparser
import shutil
class Screen(object):

    def _import_libs(self, texs, strand, wig_path, lib_dict):
        if strand == "+":
            tex = "ft"
            notex = "fn"
        else:
            tex = "rt"
            notex = "rn"
        for flib in texs:
            if (flib[1] == "tex"):
                detect = False
                lib_dict[tex].append(wig_path + "/" + flib[0])
                for nlib in texs:
                    if (nlib[1] == "notex") and \
                       (flib[2] == nlib[2]) and \
                       (flib[3] == nlib[3]):
                        detect = True
                        lib_dict[notex].append(wig_path + "/" + nlib[0])

    def gen_screenshot(self, bin_path, main_gff, side_gff, fasta, frag_wigs, 
                       tex_wigs, height, tlibs, flibs, present, output_folder):
        if os.path.exists(output_folder):
            print("Error: The " + output_folder + " already exist!!!")
            sys.exit()
        else:
            os.system("mkdir " + output_folder)
            os.system("mkdir " + output_folder + "/forward")
            os.system("mkdir " + output_folder + "/reverse")
        lib_dict = {"ft": [], "fn": [], "rt": [], "rn": [], "ff": [], "rf": []}
        f_texs = []
        r_texs = []
        if tlibs is not False:
            for lib in tlibs:
                lib_datas = lib.split(":")
                if lib_datas[0].endswith(".wig") is not True:
                    print("Error:Exist a not proper wig files!!")
                    sys.exit()
                else:
                    if lib_datas[-1] == "+":
                        f_texs.append(lib_datas)
                    else:
                        r_texs.append(lib_datas)
            f_texs = sorted(f_texs, key = lambda x: (x[1], x[2], x[3]))
            r_texs = sorted(r_texs, key = lambda x: (x[1], x[2], x[3]))
            self._import_libs(f_texs, "+", tex_wigs, lib_dict)
            self._import_libs(r_texs, "-", tex_wigs, lib_dict)
        if flibs is not False:
            for lib in flibs:
                lib_datas = lib.split(":")
                if lib_datas[0].endswith(".wig") is not True:
                    print("Error:Exist a not proper wig files!!")
                    sys.exit()
                else:
                    if lib_datas[-1] == "+":
                        lib_dict["ff"].append(frag_wigs + "/" + lib_datas[0])
                    else:
                        lib_dict["rf"].append(frag_wigs + "/" + lib_datas[0])
        command = ["python", bin_path + "/gen_screenshots.py",
                   "-mg", main_gff, "-sg", " ".join(side_gff),
                   "-os", output_folder, "-he", str(height),
                   "-f", fasta, "-p", present, "-of", output_folder + "/forward.txt",
                   "-or", output_folder + "/reverse.txt"]
        if (tlibs is not False) and (flibs is not False):
            command = command + ["-ft", " ".join(lib_dict["ft"]),
                                 "-fn", " ".join(lib_dict["fn"]),
                                 "-rt", " ".join(lib_dict["rt"]),
                                 "-rn", " ".join(lib_dict["rn"]),
                                 "-ff", " ".join(lib_dict["ff"]),
                                 "-rf", " ".join(lib_dict["rf"])]
        elif (tlibs is not False):
            command = command + ["-ft", " ".join(lib_dict["ft"]),
                                 "-fn", " ".join(lib_dict["fn"]),
                                 "-rt", " ".join(lib_dict["rt"]),
                                 "-rn", " ".join(lib_dict["rn"])]
        elif (flibs is not False):
            command = command + ["-ff", " ".join(lib_dict["ff"]),
                                 "-rf", " ".join(lib_dict["rf"])]
        else:
            print("Error: There are no wig file assigned!!!")
            sys.exit()
        os.system(" ".join(command))
