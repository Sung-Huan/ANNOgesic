#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
from subprocess import call

class Color_PNG(object):

    def Generate_color_png(self, bin_path, track_num, file_type, out_folder):
        imagemagick_path = os.environ["IMAGEMAGICK_HOME"]
        out_folder = out_folder + "/" + file_type + "/screenshots/"
        if (file_type != "TSS") and \
           (file_type != "processing_site"):
            print("Error: " + file_type + " is not a proper type.")
            print("Please assign TSS or processing_site.")
            sys.exit()
        for strain in os.listdir(out_folder):
            if os.path.isdir(out_folder + strain):
                for strand in ["forward", "reverse"]:
                    print("Running for " + strain + "_" + strand)
                    out_path = out_folder + strain + "/" + strand + "/"
                    #########################################
                    # convert original png file to svg file #
                    # and give color on it.                 #
                    #########################################
                    for screenshot in os.listdir(out_path):
                        if screenshot.endswith(".png"):
                            print("convert " + screenshot + " to svg and painting tracks now...")
                            svg_file = screenshot.replace(".png", ".svg")
                            call([imagemagick_path + "/convert", out_path + screenshot,
                                  out_path + svg_file])
                            with open(out_path + svg_file, "r") as f_h:
                                for line in f_h:
                                    line = line.strip()
                                    if line.startswith("<svg"):
                                        line = line.split(" ")
                                        height = line[-1].split("=")[-1][1:-2]
                                        width = line[1].split("=")[-1][1:-1]
                                        break
                            call(["python", bin_path + "/gen_svg.py",
                                  "-i", out_path + screenshot,
                                  "-t", str(track_num),
                                  "-he", str(height),
                                  "-w", str(width)])
                    ########################################
                    # convert to png file again            #
                    ########################################
                    for screenshot in os.listdir(out_path):
                        if screenshot.endswith(".svg"):
                            print("convert " + screenshot + " to png now...")
                            png_file = screenshot.replace(".svg", ".png")
                            print(out_path + png_file)
                            call([imagemagick_path + "/convert",
                                  "-background", "none", out_path + screenshot,
                                  out_path + png_file])
                    os.system("rm " + out_path + "*.svg")
