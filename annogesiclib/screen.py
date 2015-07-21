import os
import sys
from subprocess import call
import shutil
from annogesiclib.multiparser import Multiparser
from annogesiclib.gen_screenshots import gen_screenshot
from annogesiclib.helper import Helper


class Screen(object):

    def __init__(self, output_folder, fasta):
        self.multiparser = Multiparser()
        self.helper = Helper()
        if os.path.exists(output_folder):
            print("Error: The {0} already exist!!!".format(output_folder))
            sys.exit()
        else:
            os.mkdir(output_folder)
        filename = fasta.split("/")[-1]
        self.strain = ".".join(filename.split(".")[0:-1])
        self.helper.check_make_folder(os.path.join(output_folder, self.strain))
        self.forward_file = os.path.join(output_folder, self.strain, "forward")
        self.reverse_file = os.path.join(output_folder, self.strain, "reverse")
        os.mkdir(self.forward_file)
        os.mkdir(self.reverse_file)

    def _import_libs(self, texs, strand, wig_path, lib_dict):
        if strand == "+":
            tex = "ft"
            notex = "fn"
        else:
            tex = "rt"
            notex = "rn"
        for flib in texs:
            if (flib[1] == "tex"):
                lib_dict[tex].append(os.path.join(wig_path, flib[0]))
                for nlib in texs:
                    if (nlib[1] == "notex") and \
                       (flib[2] == nlib[2]) and \
                       (flib[3] == nlib[3]):
                        lib_dict[notex].append(os.path.join(wig_path, nlib[0]))

    def screenshot(self, main_gff, side_gffs, fasta, frag_wigs,
                   tex_wigs, height, tlibs, flibs, present, output_folder):
        lib_dict = {"ft": [], "fn": [], "rt": [], "rn": [], "ff": [], "rf": []}
        f_texs = []
        r_texs = []
        if tlibs is not None:
            for lib in tlibs:
                lib_datas = lib.split(":")
                if not lib_datas[0].endswith(".wig"):
                    print("Error:Exist a not proper wig files!!")
                    sys.exit()
                else:
                    if lib_datas[-1] == "+":
                        f_texs.append(lib_datas)
                    else:
                        r_texs.append(lib_datas)
            f_texs = sorted(f_texs, key=lambda x: (x[1], x[2], x[3]))
            r_texs = sorted(r_texs, key=lambda x: (x[1], x[2], x[3]))
            self._import_libs(f_texs, "+", tex_wigs, lib_dict)
            self._import_libs(r_texs, "-", tex_wigs, lib_dict)
        if flibs is not None:
            for lib in flibs:
                lib_datas = lib.split(":")
                if not lib_datas[0].endswith(".wig"):
                    print("Error:Exist a not proper wig files!!")
                    sys.exit()
                else:
                    if lib_datas[-1] == "+":
                        lib_dict["ff"].append(os.path.join(
                                       frag_wigs, lib_datas[0]))
                    else:
                        lib_dict["rf"].append(os.path.join(
                                       frag_wigs, lib_datas[0]))
        gen_screenshot(main_gff, self.forward_file + ".txt",
                       self.reverse_file + ".txt", output_folder, height,
                       lib_dict["ft"], lib_dict["fn"],
                       lib_dict["rt"], lib_dict["rn"],
                       lib_dict["ff"], lib_dict["rf"],
                       fasta, side_gffs, present, self.strain)
        if (tlibs is None) and (flibs is None):
            print("Error: There are no wig file assigned!!!")
            sys.exit()
