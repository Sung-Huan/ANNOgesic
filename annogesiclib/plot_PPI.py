#!/usr/bin/python

import os	
import sys
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import networkx as nx


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

def best_assign_attributes(check_NA, G, ppi, pre_ppi, first, style):
    check_NA["best"] = True
    if ppi["score"] == 0:
        if ppi["below"] >= 20:
            weight = 22
        else:
            weight = ppi["below"] + 2
    else:
        if ppi["score"] >= 20:
            weight = 22
        else:
            weight = ppi["score"] + ppi["below"] + 2
    G.add_edge(ppi["item_a"], ppi["item_b"], color=float(ppi["best"]), style=style, weight=weight)
    if first is False:
        if pre_ppi["best"] != ppi["best"]:
            check_NA["same_best"] = True

def create_node(ppis, scores, nodes, center, colors, labels1, labels2, edges, 
                G, cutoff_score, check_NA, pre_ppi):
    first = True
    for ppi in ppis:
        scores.append(ppi["score"])
        node(ppi["item_a"], nodes, center, colors, labels1, labels2)
        node(ppi["item_b"], nodes, center, colors, labels1, labels2)
        if ((ppi["item_a"], ppi["item_b"]) not in edges) or \
           ((ppi["item_b"], ppi["item_a"]) not in edges):
            edges.append((ppi["item_a"], ppi["item_b"]))
            if ppi["best"] == "NA":
                G.add_edge(ppi["item_a"], ppi["item_b"], color=-1, style='dotted', weight=2)
            elif float(ppi["best"]) <= cutoff_score:
                best_assign_attributes(check_NA, G, ppi, pre_ppi, first, "dashdot")
            else:
                best_assign_attributes(check_NA, G, ppi, pre_ppi, first, "solid")
            pre_ppi = ppi
            first = False
    G.add_nodes_from(nodes)
    return pre_ppi

def modify_label(labels2, new_labels):
    for key, value in labels2.items():
        if "_" in value:
            new_labels[key] = value.replace("_", "\n")
        else:
            new_labels[key] = value

def plot_text(check_NA, plt, ppis, ppi, color_edge):
    if (check_NA["best"] is True):
        if len(ppis) < 2:
            plt.text(0.0, -0.1, 
                     " ".join(["the number of literatures supported is", str(ppi["best"])]), 
                     color='blue', verticalalignment='bottom', horizontalalignment='center', 
                     fontsize=12)
        elif check_NA["same_best"] is False:
            plt.text(0.0, -0.1, 
                     " ".join(["all numbers of literature supported are", str(ppi["best"])]), 
                     color='blue', verticalalignment='bottom', horizontalalignment='center', 
                     fontsize=12)
        else:
            plt.colorbar(color_edge)

def plot(ppis, center, strain, cutoff_score, node_size, out_folder):
    fig = plt.figure(figsize=(12, 12))
    G = nx.Graph()
    nodes = []
    edges = []
    labels1 = {}
    labels2 = {}
    sizes = []
    colors = {}
    tick = 0
    check_NA = {"number": False, "best": False, "same_number": False, "same_best": False}
    pre_ppi = ""
    scores = []
    weights = []
    pre_ppi = create_node(ppis, scores, nodes, center, colors, labels1, labels2, edges, 
                          G, cutoff_score, check_NA, pre_ppi)
    pos=nx.spring_layout(G)
    color_list = []
    for color in colors.values():
        color_list.append(color)
    nx.draw_networkx_nodes(G, pos, node_size=node_size, node_shape='o', 
                    nodelist = colors.keys(), node_color = color_list, linewidths=1)
    connects = G.edges()
    for weight in G.edges(data=True):
         weights.append(weight[2]["weight"])
    colors = [G[u][v]['color'] for u,v in edges]
    styles = [G[u][v]['style'] for u,v in edges]
    color_edge = nx.draw_networkx_edges(G, pos, edges=connects, edge_color=colors, 
                                        style=styles, len=10, width=weights)
    nx.draw_networkx_labels(G,pos,labels1,font_size=12, font_weight='bold')
    new_labels = {}
    modify_label(labels2, new_labels)
    nx.draw_networkx_labels(G,pos,new_labels,font_size=10, font_weight='bold')
    plt.title("|".join([center["locus_tag"], 
              " ".join([center["gene_name"], "(based on the score of best literature)"])]), 
              fontsize="16")
    plot_text(check_NA, plt, ppis, pre_ppi, color_edge)
    plt.axis('off')
    if strain not in os.listdir(out_folder):
        os.mkdir(os.path.join(out_folder, strain))
    plt.savefig(os.path.join(out_folder, strain, 
                "_".join([center["locus_tag"], center["gene_name"] + ".png"])), 
                bbox_inches="tight")
    plt.clf()
    plt.close('all')

def score_compare(score, scores, cutoff_score, ppi):
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

def plot_ppi(PPI_file, cutoff_score, out_folder, node_size):
    ppis = []
    first = True
    scores = {"score": 0, "below": 0}
    center = {}
    fh = open(PPI_file, "r")
    print(PPI_file)
    for row in csv.reader(fh, delimiter="\t"):
        if row[0].startswith("Interaction"):
            if first:
                pass
            else:
                assign_score_below(pre_ppi, scores, ppis)
                plot(ppis, center, pre_ppi["strain"], cutoff_score, 
                     node_size, out_folder)
                scores = {"score": 0, "below": 0}
                ppis = []
                first = True
            datas = row[0].split(" | ")
            center["locus_tag"] = datas[0].split(" ")[-1]
            center["gene_name"] = datas[-1]
            print("plotting {0}".format(center["gene_name"]))
        elif row[0] == "strain":
            pass
        else:
            ppi = {"strain": row[0], "item_a": row[1], "item_b": row[2],
                   "mode": row[3]}
            if first:
                first = False
                score_compare(row[8], scores, cutoff_score, ppi)
                ppi["best"] = row[8]
            else:
                if (ppi["strain"] == pre_ppi["strain"]) and \
                   (ppi["item_a"] == pre_ppi["item_a"]) and \
                   (ppi["item_b"] == pre_ppi["item_b"]):
                    get_best(pre_ppi, ppi, row)
                    score_compare(row[8], scores, cutoff_score, ppi)
                else:
                    assign_score_below(pre_ppi, scores, ppis)
                    scores = {"score": 0, "below": 0}
                    score_compare(row[8], scores, cutoff_score, ppi)
                    ppi["best"] = row[8]
            pre_ppi = ppi
    assign_score_below(pre_ppi, scores, ppis)
    plot(ppis, center, pre_ppi["strain"], cutoff_score, node_size, out_folder)
