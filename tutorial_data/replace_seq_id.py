#!/usr/bin/python

import os
import shutil
import argparse

__author__ = "Sung-Huan Yu <shyu@biochem.mpg.de>"
__email__ = "shyu@biochem.mpg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_wig_folder",help="input wig file")
parser.add_argument("-n","--strain_name",help="strain_name")
args = parser.parse_args()

def main():
    for wig in os.listdir(args.input_wig_folder):
        out = open("tmp", "w")
        with open(os.path.join(args.input_wig_folder, wig)) as fh:
            for line in fh:
                if line.startswith("variableStep"):
                    data = line.split(" ")
                    choms = data[1].split("=")
                    choms[-1] = args.strain_name
                    data[1] = "=".join(choms)
                    out.write(" ".join(data))
                else:
                    out.write(line)
        out.close()
        os.remove(os.path.join(args.input_wig_folder, wig))
        shutil.move("tmp", os.path.join(args.input_wig_folder, wig))

if __name__ == "__main__":
    main()
