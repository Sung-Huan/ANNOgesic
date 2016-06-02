from copy import deepcopy
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')


def fig(rowlabels, collabels, cells, filename, max_color, min_color):
    row_num = len(rowlabels) / 100
    if row_num == 0:
        row_num = 1
    col_num = len(collabels) / 8
    if col_num == 0:
        col_num = 1
    plt.figure(figsize=(18*col_num, 10*row_num), edgecolor=None)
    img = plt.imshow(cells, interpolation='none', aspect='auto', cmap="RdBu_r")
    plt.xticks(range(len(collabels)), collabels, fontsize=6)
    plt.yticks(range(len(rowlabels)), rowlabels, fontsize=6)
    plt.colorbar(fraction=0.046, pad=0.04)
    img.set_clim(vmin=min_color, vmax=max_color)
    plt.savefig(filename)


def plot_table(plots, max_color, min_color, filename):
    rowlabels = []
    collabels = []
    cells = []
    first = True
    t_num = 0
    for plot in plots:
        for key, value in plot.items():
            rowlabels.append(key)
        cell = []
        for cond, tracks in value.items():
            for track, cover in tracks.items():
                if first:
                    name = track
                    if len(track) > 16:
                        diff = int(len(name) / 16)
                        for i in range(diff):
                            name = (name[:(16)*(i+1)+i] + "\n" +
                                    name[(16)*(i+1)+i:])
                    collabels.append(name)
                cell.append(round(cover, 1))
        cells.append(deepcopy(cell))
        first = False
        if len(rowlabels) >= 500:
            plotname = (filename[:-4] + "_" + str(t_num) + "-" +
                        str(t_num + 500) + ".png")
            fig(rowlabels, collabels, cells, plotname, max_color, min_color)
            t_num = t_num + 500
            rowlabels = []
            cells = []
    if t_num == 0:
        fig(rowlabels, collabels, cells, filename, max_color, min_color)
