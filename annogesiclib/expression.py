import os
import sys
import shutil
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.gene_express_analysis import gene_expression


class Expression(object):

    def __init__(self, gffs):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.out_folder = os.path.join(gffs, "for_libs")
        if os.path.exists(self.out_folder):
            shutil.rmtree(self.out_folder)
        os.mkdir(self.out_folder)
        self.stat = os.path.join(self.out_folder, "statistics")
        os.mkdir(self.stat)
        self.gff_folder = os.path.join(self.out_folder, "gffs")
        os.mkdir(self.gff_folder)
        self.merge_wigs = os.path.join(gffs, "merge_wigs")
        if os.path.exists(self.merge_wigs):
            shutil.rmtree(self.merge_wigs)

    def _get_replicates(self, replicates_tex, replicates_frag):
        if (replicates_tex is not None) and (
                replicates_frag is not None):
            replicates = {"tex": int(replicates_tex),
                          "frag": int(replicates_frag)}
        elif replicates_tex is not None:
            replicates = {"tex": int(replicates_tex), "frag": -1}
        elif replicates_frag is not None:
            replicates = {"tex": -1, "frag": int(replicates_frag)}
        else:
            print("Error:No replicates number assign!!!")
            sys.exit()
        return replicates

    def expression(self, tex_libs, frag_libs, tex_notex, replicates_tex,
                   replicates_frag, tex_wigs, frag_wigs, percent_tex,
                   percent_frag, cutoff_coverage, gffs, features,
                   cover_type, max_color, min_color):
        replicates = self._get_replicates(replicates_tex, replicates_frag)
        if (tex_libs is not None) and (frag_libs is not None):
            input_libs = tex_libs + frag_libs
        elif tex_libs is not None:
            input_libs = tex_libs
        elif frag_libs is not None:
            input_libs = frag_libs
        else:
            print("Error: plese assign the libraries!!\n")
            sys.exit()
        if (tex_wigs is not None) and (frag_wigs is not None):
            merge_wigs = self.merge_wigs
            os.mkdir(merge_wigs)
            for wig in os.listdir(tex_wigs):
                if os.path.isfile(os.path.join(tex_wigs, wig)):
                    shutil.copy(os.path.join(tex_wigs, wig), merge_wigs)
            for wig in os.listdir(frag_wigs):
                if os.path.isfile(os.path.join(frag_wigs, wig)):
                    shutil.copy(os.path.join(frag_wigs, wig), merge_wigs)
        elif tex_wigs is not None:
            merge_wigs = tex_wigs
        elif frag_wigs is not None:
            merge_wigs = frag_wigs
        else:
            print("Error: plese assign the wiggle files!!\n")
            sys.exit()
        wig_f_file = os.path.join(merge_wigs, "whole_forward.wig")
        wig_r_file = os.path.join(merge_wigs, "whole_reverse.wig")
        for wig in os.listdir(merge_wigs):
            for lib in input_libs:
                if (wig in lib) and (lib[-1] == "+"):
                    self.helper.merge_file(os.path.join(merge_wigs, wig),
                                           wig_f_file)
                elif (wig in lib) and (lib[-1] == "-"):
                    self.helper.merge_file(os.path.join(merge_wigs, wig),
                                           wig_r_file)
        print("Computing expression analysis...")
        gene_expression(input_libs, gffs, percent_tex, percent_frag,
                        wig_f_file, wig_r_file, features, merge_wigs,
                        cutoff_coverage, tex_notex, replicates, self.stat,
                        self.gff_folder, cover_type, max_color, min_color)
        os.remove(wig_f_file)
        os.remove(wig_r_file)
        if os.path.exists(self.merge_wigs):
            shutil.rmtree(self.merge_wigs)
