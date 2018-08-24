import os
import math
import matplotlib as mpl
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.helper import Helper
import numpy as np
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')


def plot(lens, out_figure):
    ticks = max(lens) / 50
    bin_num = np.arange(0, max(lens), ticks)
    n, bins, hist1 = plt.hist(lens, bin_num,
                              color="#FF9999", label='Transcript',
                              edgecolor='black', linewidth=1)
    plt.xlabel("Transcript_length (nt)")
    plt.ylabel("Amount")
    plt.savefig(out_figure)
    plt.clf()


def plot_tran(tran_folder, stat_folder, max_dist):
    lens = []
    less = []
    for tran in os.listdir(tran_folder):
        if tran.endswith(".gff"):
            prefix = tran.replace("_transcript.gff", "")
            gff_f = open(os.path.join(tran_folder, tran), "r")
            for entry in Gff3Parser().entries(gff_f):
                if entry.feature == "transcript":
                    lens.append(entry.end - entry.start)
                    if entry.end - entry.start <= max_dist:
                        less.append(entry.end - entry.start)
            plot(lens, os.path.join(stat_folder, prefix + "_length_all.png"))
            plot(less, os.path.join(stat_folder, prefix + "_length_less_" + 
                                    str(max_dist) + ".png"))

