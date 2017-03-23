import os
import networkx as nx
import matplotlib
matplotlib.use('Agg')
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
plt.rcParams['image.cmap'] = 'RdBu_r'


def node(item, nodes, center, colors, labels1, labels2):
    if item not in nodes:
        nodes.append(item)
        if len(item) > 5:
            labels2[item] = item
            labels1[item] = ""
        else:
            labels1[item] = item
            labels2[item] = ""
        if (center["locus_tag"] == item) or \
           (center["gene_name"] == item):
            colors[item] = '#FFFF66'
        else:
            colors[item] = '#CCFFCC'


def get_largest_compare(tick, score):
    same = True
    if (score >= 20) and tick >= 20:
        pass
    else:
        if score != tick:
            same = False
    if score > tick:
        if score >= 20:
            tick = 20
        else:
            tick = score
    return (tick, same)


def add_edge(G, ppi, style, weight, colorppi):
    G.add_edge(ppi["item_a"], ppi["item_b"],
               color=float(colorppi), style=style, weight=weight)


def add_node(G, nodes):
    G.add_nodes_from(nodes)


def best_assign_attributes(check_na, G, ppi, pre_ppi, first, style):
    check_na["best"] = True
    if ppi["score"] == 0:
        if ppi["below"] >= 20:
            weight = 22
        else:
            weight = ppi["below"] + 1
    else:
        if ppi["score"] >= 20:
            weight = 22
        else:
            weight = ppi["score"] + ppi["below"] + 1
    add_edge(G, ppi, style, weight, ppi["best"])
    if not first:
        if pre_ppi["best"] != ppi["best"]:
            check_na["same_best"] = True


def create_node(ppis, scores, nodes, center, colors, labels1, labels2, edges,
                G, cutoff_score, check_na, pre_ppi):
    first = True
    for ppi in ppis:
        scores.append(ppi["score"])
        node(ppi["item_a"], nodes, center, colors, labels1, labels2)
        node(ppi["item_b"], nodes, center, colors, labels1, labels2)
        if ((ppi["item_a"], ppi["item_b"]) not in edges) or \
           ((ppi["item_b"], ppi["item_a"]) not in edges):
            edges.append((ppi["item_a"], ppi["item_b"]))
            if ppi["best"] == "NA":
                add_edge(G, ppi, 'dashed', 1, -1)
            elif float(ppi["best"]) <= cutoff_score:
                best_assign_attributes(check_na, G, ppi,
                                       pre_ppi, first, "dashdot")
            else:
                best_assign_attributes(check_na, G, ppi,
                                       pre_ppi, first, "solid")
            pre_ppi = ppi
            first = False
    add_node(G, nodes)
    return pre_ppi


def modify_label(labels2, new_labels):
    for key, value in labels2.items():
        if "_" in value:
            new_labels[key] = value.replace("_", "\n")
        else:
            new_labels[key] = value


def plot_text(check_na, plt, ppis, ppi, color_edge):
    na = False
    if check_na["na"]:
        na = True
    elif check_na["best"]:
        if len(ppis) < 2:
            na = True
        else:
            cbar = plt.colorbar(color_edge)
            cbar.ax.tick_params(labelsize=16)
    return na


def nx_node(G, pos, node_size, colors, color_list):
    '''draw the node'''
    nx.draw_networkx_nodes(G, pos, node_size=node_size, node_shape='o',
                           nodelist=colors.keys(), node_color=color_list,
                           linewidths=1)


def nx_edge(G, pos, edges, colors, styles, weights):
    '''draw the edge'''
    color_edge = (nx.draw_networkx_edges(G, pos, edges=edges,
                  edge_color=colors, style=styles, width=weights,
                  edge_vmin=-1, edge_vmax=1, vmin=-1, vmax=1))
    return color_edge


def nx_label(G, pos, labels, size):
    '''setup the label of network'''
    nx.draw_networkx_labels(G, pos, labels, font_size=size, font_weight='bold')


def nx_color_style(G, edges):
    '''setup the color of network'''
    colors = []
    styles = []
    check_na = True
    for u, v in edges:
        colors.append(G[u][v]['color'])
        styles.append(G[u][v]['style'])
        if (G[u][v]['style'] == "solid") or (
                G[u][v]['style'] == "dashdot"):
            check_na = False
    return colors, styles, check_na


def print_title(plt, na, center):
    if not na:
        plt.title("|".join([center["locus_tag"],
                           " ".join([center["gene_name"],
                                     "(based on the score of best literature)"])]),
                  fontsize="16")
    else:
        plt.title("|".join([center["locus_tag"],
                           " ".join([center["gene_name"],
                                     "(based on the score of best literature)"])]) + \
                  "\n the numbers of supported literatures in all interactions are 0",
                  fontsize="16")

def plot(ppis, center, strain, cutoff_score, node_size, out_folder):
    nodes = []
    edges = []
    labels1 = {}
    labels2 = {}
    colors = {}
    check_na = {"number": False, "best": False, "na": False,
                "same_number": False, "same_best": False}
    pre_ppi = ""
    scores = []
    weights = []
    plt.figure(figsize=(15, 15))
    G = nx.Graph()
    pre_ppi = create_node(ppis, scores, nodes, center, colors,
                          labels1, labels2, edges, G,
                          cutoff_score, check_na, pre_ppi)
    pos = nx.spring_layout(G, k=2, scale=3, iterations=20)
    color_list = []
    for color in colors.values():
        color_list.append(color)
    nx_node(G, pos, node_size, colors, color_list)
    connects = G.edges()
    for weight in G.edges(data=True):
        if weight[2]["weight"] <= 30:
            weights.append(weight[2]["weight"])
        else:
            weights.append(30)
    colors, styles, check_na["na"] = nx_color_style(G, connects)
    color_edge = nx_edge(G, pos, connects, colors, styles, weights)
    nx_label(G, pos, labels1, 12)
    new_labels = {}
    modify_label(labels2, new_labels)
    nx_label(G, pos, new_labels, 10)
    na = plot_text(check_na, plt, ppis, pre_ppi, color_edge)
    print_title(plt, na, center)
    plt.axis('off')
    if strain not in os.listdir(out_folder):
        os.mkdir(os.path.join(out_folder, strain))
    plt.savefig(os.path.join(out_folder, strain,
                "_".join([center["locus_tag"], center["gene_name"] + ".png"])),
                bbox_inches="tight")
    plt.clf()
    plt.close('all')
    return check_na

def score_compare(score, scores, cutoff_score, ppi):
    '''check the number of literatures which are pass the cutoff'''
    if score == "NA":
        ppi["score"] = 0
        ppi["below"] = 0
    elif float(score) >= cutoff_score:
        scores["score"] += 1
    else:
        scores["below"] += 1


def assign_score_below(pre_ppi, scores, ppis):
    if "score" not in pre_ppi.keys():
        pre_ppi["score"] = scores["score"]
    if "below" not in pre_ppi.keys():
        pre_ppi["below"] = scores["below"]
    ppis.append(pre_ppi)


def get_best(pre_ppi, ppi, row):
    '''get the best score of PPI'''
    if "best" not in pre_ppi.keys():
        ppi["best"] = row[8]
    else:
        if pre_ppi["best"] == "NA":
            ppi["best"] = row[8]
        else:
            if float(pre_ppi["best"]) < float(row[8]):
                ppi["best"] = row[8]
            else:
                ppi["best"] = pre_ppi["best"]


def interaction(first, pre_ppi, scores, ppis, match, center, cutoff_score,
                node_size, out_folder):
    '''check the interaction of two proteins'''
    if first:
        pass
    else:
        assign_score_below(pre_ppi, scores, ppis)
        if match:
            plot(ppis, center, pre_ppi["strain"], cutoff_score,
                 node_size, out_folder)
            match = False
        else:
            print("No interacted partner with {0} | {1}".format(
                  center["locus_tag"], center["gene_name"]))
        scores = {"score": 0, "below": 0}
        ppis = []
        first = True
    return first, scores, match, ppis


def plot_ppi(PPI_file, cutoff_score, out_folder, node_size):
    '''plot the network of PPI'''
    ppis = []
    first = True
    pre_ppi = None
    scores = {"score": 0, "below": 0}
    center = {}
    start = False
    match = False
    with open(PPI_file) as fh:
        for line in fh:
            line = line.strip()
            row = line.split("\t")
            start = True
            if row[0].startswith("Interaction"):
                first, scores, match, ppis = interaction(
                        first, pre_ppi, scores, ppis, match, center,
                        cutoff_score, node_size, out_folder)
                datas = row[0].split(" | ")
                center["locus_tag"] = datas[0].split(" ")[-1]
                center["gene_name"] = datas[-1]
                print("Plotting {0}".format(center["gene_name"]))
            elif row[0] == "Genome":
                pass
            else:
                ppi = {"strain": row[0], "item_a": row[1], "item_b": row[2],
                       "mode": row[3]}
                if (ppi["item_a"] == center["locus_tag"]) or (
                        ppi["item_a"] == center["gene_name"]) or (
                        ppi["item_b"] == center["locus_tag"]) or (
                        ppi["item_b"] == center["gene_name"]):
                    match = True
                if first:
                    first = False
                    score_compare(row[8], scores, cutoff_score, ppi)
                    ppi["best"] = row[8]
                else:
                    if (ppi["strain"] == pre_ppi["strain"]) and (
                            ppi["item_a"] == pre_ppi["item_a"]) and (
                            ppi["item_b"] == pre_ppi["item_b"]):
                        get_best(pre_ppi, ppi, row)
                        score_compare(row[8], scores, cutoff_score, ppi)
                    else:
                        assign_score_below(pre_ppi, scores, ppis)
                        scores = {"score": 0, "below": 0}
                        score_compare(row[8], scores, cutoff_score, ppi)
                        ppi["best"] = row[8]
                pre_ppi = ppi
    if start and match:
        assign_score_below(pre_ppi, scores, ppis)
        plot(ppis, center, pre_ppi["strain"],
             cutoff_score, node_size, out_folder)
    elif not start:
        print("No proper result can be retrieved in " + PPI_file)
    elif not match:
        print("No interacted partner with {0} | {1}".format(
               center["locus_tag"], center["gene_name"]))
