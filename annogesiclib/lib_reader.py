import sys
import os
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.parser_wig import WigParser

def read_libs(input_libs, wig_folder):
    libs = []
    texs = {}
    for lib in input_libs:
        datas = lib.split(":")
        for wig in os.listdir(wig_folder):
            if wig == datas[0]:
                with open(os.path.join(wig_folder, wig), "r") as w_h:
                    for line in w_h:
                        line = line.strip()
                        if line.startswith("track"):
                            name = line.split("=")[-1][1:-1]
                            break
        if (datas[1] == "tex") or (datas[1] == "notex"):
            cond = "texnotex"
        else:
            cond = datas[1]
        libs.append({"name": name, "type": datas[1],
                     "cond": "_".join([datas[2], cond]),
                     "rep": datas[3], "strand": datas[4]})
    for lib1 in libs:
        if lib1["type"] == "frag":
            type_ = "frag"
        elif (lib1["type"] == "tex") or (lib1["type"] == "notex"):
            prefix1 = lib1["cond"].split("_")[0]
            for lib2 in libs:
                prefix2 = lib2["cond"].split("_")[0]
                if (prefix1 == prefix2) and \
                   (lib1["rep"] == lib2["rep"]) and \
                   (lib1["type"] == "tex") and \
                   (lib2["type"] == "notex") and \
                   (lib1["strand"] == lib2["strand"]):
                    texs[lib1["name"] + "@AND@" + lib2["name"]] = 0
            type_ = "tex"
        else:
            print("Error:library type (fragmented or Tex treated??)")
            sys.exit()
    return libs, texs

def read_wig(filename, strand, libs):
    wig_parser = WigParser()
    wigs = {}
    if filename is not False:
        wig_fh = open(filename)
        for entry in wig_parser.parser(wig_fh, strand):
            if entry.strain not in wigs.keys():
                strain = entry.strain
                wigs[strain] = {}
                for lib in libs:
                    if lib["cond"] not in wigs[strain]:
                        wigs[strain][lib["cond"]] = {}
            for lib in libs:
                if (lib["name"] == entry.track) and (
                    lib["strand"] == entry.strand):
                    if entry.track not in wigs[strain][lib["cond"]].keys():
                        wigs[strain][lib["cond"]][entry.track] = []
                    wigs[strain][lib["cond"]][entry.track].append({
                         "pos": entry.pos, "coverage": entry.coverage,
                         "strand": entry.strand, "type": lib["type"]})
        wig_fh.close()
    return wigs
