import sys
import os
import csv
#from annogesiclib.gff3 import Gff3Parser
#import annogesiclib.parser_wig as par_wig
from gff3 import Gff3Parser
import parser_wig as par_wig

def read_libs(libs, texs, input_libs, wig_folder):
    for lib in input_libs:
        datas = lib.split(":")
        for wig in os.listdir(wig_folder):
            if wig == datas[0]:
                with open(wig_folder + "/" + wig, "r") as w_h:
                    for line in w_h:
                        line = line.strip()
                        if line.startswith("track"):
                            name = line.split("=")[-1][1:-1]
                            break
        if (datas[1] == "tex") or (datas[1] == "notex"):
            cond = "texnotex"
        else:
            cond = datas[1]
        libs.append({"name": name, "type": datas[1], "cond": "_".join([datas[2], cond]),
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
                    texs[lib1["name"] + "_" + lib2["name"]] = 0
            type_ = "tex"
        else:
            print("Error:not correct library type (fragmented or Tex treated??)")
            sys.exit()

def read_wig(wigs, filename, strand, libs):
    wig_parser = par_wig.parser_wig()
    if filename is not False:
        for entry in wig_parser.parser(filename, strand):
            if entry.strain not in wigs.keys():
                strain = entry.strain
                wigs[strain] = {}
                for lib in libs:
                    if lib["cond"] not in wigs[strain]:
                        wigs[strain][lib["cond"]] = {}
            for lib in libs:
                if lib["name"] == entry.track:
                    if entry.track not in wigs[strain][lib["cond"]].keys():
                        wigs[strain][lib["cond"]][entry.track] = []
                    wigs[strain][lib["cond"]][entry.track].append({
                         "pos": entry.pos, "coverage": entry.coverage,
                         "strand": entry.strand, "type": lib["type"]})
