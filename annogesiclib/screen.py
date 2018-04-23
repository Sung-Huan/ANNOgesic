import os
import sys
from annogesiclib.gen_screenshots import gen_screenshot
from annogesiclib.helper import Helper


class Screen(object):
    '''generation of screenshot'''

    def __init__(self, args_sc, out_folder):
        self.helper = Helper()
        args_sc.output_folder = out_folder
        filename = args_sc.fasta.split("/")[-1]
        self.strain = ".".join(filename.split(".")[0:-1])
        self.helper.check_make_folder(os.path.join(args_sc.output_folder,
                                                   self.strain))
        self.forward_file = os.path.join(args_sc.output_folder,
                                         self.strain, "forward")
        self.reverse_file = os.path.join(args_sc.output_folder,
                                         self.strain, "reverse")
        os.mkdir(self.forward_file)
        os.mkdir(self.reverse_file)

    def _import_libs(self, texs, strand, lib_dict):
        if strand == "+":
            tex = "ft"
            notex = "fn"
        else:
            tex = "rt"
            notex = "rn"
        for flib in texs:
            if (flib[1] == "tex"):
                lib_dict[tex].append(flib[0])
                for nlib in texs:
                    if (nlib[1] == "notex") and \
                       (flib[2] == nlib[2]) and \
                       (flib[3] == nlib[3]):
                        lib_dict[notex].append(nlib[0])

    def screenshot(self, args_sc, log):
        lib_dict = {"ft": [], "fn": [], "rt": [], "rn": [], "ff": [], "rf": []}
        f_texs = []
        r_texs = []
        if args_sc.tlibs is not None:
            for lib in args_sc.tlibs:
                lib_datas = lib.split(":")
                if not lib_datas[0].endswith(".wig"):
                    log.write("Wiggle files should end with .wig.\n")
                    print("Error: Wiggle files should end with .wig!")
                    sys.exit()
                else:
                    if lib_datas[-1] == "+":
                        f_texs.append(lib_datas)
                    else:
                        r_texs.append(lib_datas)
            f_texs = sorted(f_texs, key=lambda x: (x[1], x[2], x[3]))
            r_texs = sorted(r_texs, key=lambda x: (x[1], x[2], x[3]))
            self._import_libs(f_texs, "+", lib_dict)
            self._import_libs(r_texs, "-", lib_dict)
        if args_sc.flibs is not None:
            for lib in args_sc.flibs:
                lib_datas = lib.split(":")
                if not lib_datas[0].endswith(".wig"):
                    log.write("Wiggle files should end with .wig.\n")
                    print("Error: Wiggle files should end with .wig!")
                    sys.exit()
                else:
                    if lib_datas[-1] == "+":
                        lib_dict["ff"].append(lib_datas[0])
                    else:
                        lib_dict["rf"].append(lib_datas[0])
        log.write("Running gen_screenshots.py to generate IGV batch script.\n")
        gen_screenshot(args_sc, lib_dict, self.forward_file + ".txt",
                       self.reverse_file + ".txt", self.strain)
        log.write("\t" + self.forward_file + ".txt is generated.\n")
        log.write("\t" + self.reverse_file + ".txt is generated.\n")
        if (args_sc.tlibs is None) and (args_sc.flibs is None):
            log.write("No wig files can be found.\n")
            print("Error: There is no wig file assigned!")
            sys.exit()
