#!/usr/bin/python
import os
from subprocess import call
from annogesiclib.gen_svg import gen_svg
from annogesiclib.helper import Helper


class ColorPNG(object):

    def generate_color_png(self, track_num, out_folder, imagemagick_path):
        out_folder = os.path.join(out_folder, "screenshots")
        for strain in os.listdir(out_folder):
            if os.path.isdir(os.path.join(out_folder, strain)):
                for strand in ["forward", "reverse"]:
                    print("Running for {0}_{1}".format(strain, strand))
                    out_path = os.path.join(out_folder, strain, strand)
                    # convert original png file to svg file and give color on it.
                    for screenshot in os.listdir(out_path):
                        if screenshot.endswith(".png"):
                            print("convert {0} to svg and painting tracks now...".format(
                                  screenshot))
                            svg_file = screenshot.replace(".png", ".svg")
                            call([imagemagick_path,
                                  os.path.join(out_path, screenshot),
                                  os.path.join(out_path, svg_file)])
                            with open(os.path.join(out_path, svg_file), "r") as f_h:
                                for line in f_h:
                                    line = line.strip()
                                    if line.startswith("<svg"):
                                        line = line.split(" ")
                                        height = line[-1].split("=")[-1][1:-2]
                                        width = line[1].split("=")[-1][1:-1]
                                        break
                            gen_svg(os.path.join(out_path, screenshot),
                                    track_num, height, width)
                    # convert to png file again
                    for screenshot in os.listdir(out_path):
                        if screenshot.endswith(".svg"):
                            print("convert {0} to png now...".format(
                                  screenshot))
                            png_file = screenshot.replace(".svg", ".png")
                            print(os.path.join(out_path, png_file))
                            call([imagemagick_path, "-background", "none",
                                  os.path.join(out_path, screenshot),
                                  os.path.join(out_path, png_file)])
                    Helper().remove_all_content(out_path, ".svg", "file")
