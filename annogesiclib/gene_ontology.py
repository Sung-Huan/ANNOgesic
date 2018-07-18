import os
import csv
from annogesiclib.gff3 import Gff3Parser
import numpy as np
import copy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')


def import_uniprot_data(entry, name_list, feature):
    '''import uniprot to dict'''
    ref_name = entry.attributes[feature]
    if ref_name not in name_list:
        name_list.add(ref_name)


def compare_cds_tran(gffs, trans):
    '''compare CDS with transcript to get the expressed CDS'''
    new_gffs = []
    for gff in gffs:
        for ta in trans:
            if (gff.seq_id == ta.seq_id) and (
                    gff.strand == ta.strand):
                if ((gff.end < ta.end) and (
                         gff.end > ta.start) and (
                         gff.start <= ta.start)) or (
                        (gff.start > ta.start) and (
                         gff.start < ta.end) and (
                         gff.end >= ta.end)) or (
                        (gff.end >= ta.end) and (
                         gff.start <= ta.start)) or (
                        (gff.end <= ta.end) and (
                         gff.start >= ta.start)):
                    new_gffs.append(gff)
                    break
    return new_gffs


def get_go_id(gffs, id_, uni_lines, gos):
    '''get the GO id of CDS'''
    detect = False
    for gff in gffs:
        if ("Name" in gff.attributes.keys()):
            if (id_ == gff.attributes["Name"]):
                detect = True
        if ("protein_id" in gff.attributes.keys()):
            if (id_ == gff.attributes["protein_id"]):
                detect = True
        if detect:
            detect = False
            gos.append({"strain": gff.seq_id, "strand": gff.strand,
                        "start": gff.start, "end": gff.end,
                        "protein_id": id_, "go": uni_lines[6]})
            gff.attributes["print"] = True


def print_go(gos, out):
    gos = sorted(gos, key=lambda x: (x["strain"], x["start"],
                                     x["end"], x["strand"]))
    pre_go = None
    for go in gos:
        if (go != pre_go) and (pre_go is not None):
            if (go["strain"] == pre_go["strain"]) and (
                    go["strand"] == pre_go["strand"]) and (
                    go["start"] == pre_go["start"]) and (
                    go["end"] == pre_go["end"]) and (
                    go["protein_id"] == pre_go["protein_id"]):
                go_ids = []
                for go_id in go["go"].split("; "):
                    if go_id not in go_ids:
                        go_ids.append(go_id)
                for go_id in pre_go["go"].split("; "):
                    if go_id not in go_ids:
                        go_ids.append(go_id)
                pre_go["go"] = "; ".join(go_ids)
                out.write("\t".join([pre_go["strain"], pre_go["strand"],
                          str(pre_go["start"]), str(pre_go["end"]),
                          pre_go["protein_id"], pre_go["go"]]) + "\n")
                pre_go["print"] = True
                go["print"] = True
            else:
                if "print" not in pre_go.keys():
                    out.write("\t".join([pre_go["strain"], pre_go["strand"],
                              str(pre_go["start"]), str(pre_go["end"]),
                              pre_go["protein_id"], pre_go["go"]]) + "\n")
                    pre_go["print"] = True
        pre_go = copy.deepcopy(go)
    if "print" not in pre_go.keys():
        out.write("\t".join([pre_go["strain"], pre_go["strand"],
                  str(pre_go["start"]), str(pre_go["end"]),
                  pre_go["protein_id"], pre_go["go"]]) + "\n")

def retrieve_uniprot(database_file, gff_file, out_file, tran_file, type_):
    '''Retrieve the GO term from Uniprot'''
    name_list = set()
    gffs = []
    out = open(out_file, "w")
    out.write("\t".join(["Genome", "Strand", "Start", "End",
                         "Protein_id", "Go_term"]) + "\n")
    for entry in Gff3Parser().entries(open(gff_file)):
        if entry.feature == "CDS":
            if ("Name" in entry.attributes.keys()) and (
                    "protein_id" in entry.attributes.keys()):
                if entry.attributes["Name"] == entry.attributes["protein_id"]:
                    import_uniprot_data(entry, name_list, "Name")
                else:
                    import_uniprot_data(entry, name_list, "Name")
                    import_uniprot_data(entry, name_list, "protein_id")
            elif ("Name" in entry.attributes.keys()):
                import_uniprot_data(entry, name_list, "Name")
            elif ("protein_id" in entry.attributes.keys()):
                import_uniprot_data(entry, name_list, "protein_id")
            gffs.append(entry)
    if (type_ == "express") and (tran_file is not None):
        trans = []
        for entry in Gff3Parser().entries(open(tran_file)):
            trans.append(entry)
        new_gffs = compare_cds_tran(gffs, trans)
        gffs = new_gffs
    idmapping = open(database_file, "r")
    gos = []
    for uni_id in idmapping:
        uni_line = uni_id.rstrip("\n")
        uni_lines = uni_line.split("\t")
        uni_ids = uni_lines[3].split(";")
        for id_ in uni_ids:
            id_ = id_.strip()
            if id_ in name_list:
                get_go_id(gffs, id_, uni_lines, gos)
    for gff in gffs:
        if "print" not in gff.attributes.keys():
            gos.append({"strain": gff.seq_id, "strand": gff.strand,
                        "start": gff.start, "end": gff.end,
                        "protein_id": id_, "go": ""})
    print_go(gos, out)
    out.close()
    idmapping.close()


def plot(total_nums, strain, filename, total, out_folder):
    '''plot the distribution of GO term by GOslim'''
    sort_total_nums = sorted(total_nums.items(),
                             key=lambda x: (x[1]), reverse=True)
    classes = []
    nums = []
    width = 0.4
    plt.figure(figsize=(16, 12))
    for total_num in sort_total_nums:
        class_ = total_num[0]
        num = total_num[1]
        if class_ != "total":
            percent = (float(num) / float(total)) * 100
            classes.append(class_.replace("_", " "))
            nums.append(num)
    ind = np.arange(len(nums))
    plt.bar(ind, nums, width, color='#FF9999')
    if filename == "three_roots":
        title = "Distribution of GO term hits in the three root classes"
    elif (filename == "molecular_function") or (
            filename == "cellular_component") or (
            filename == "biological_process"):
        tag = filename.replace("_", " ")
        title = " -- ".join(["Distribution of GO term of the class", tag])
    plt.title(title, fontsize=22)
    plt.ylabel('Amount', fontsize=16)
    plt.xlim([0, len(nums) + 1])
    plt.yticks(fontsize=16)
    plt.xticks(ind+width, classes, rotation=45, fontsize=16, ha='right')
    plt.tight_layout(3, None, None, None)
    plt.savefig(os.path.join(out_folder,
                             "_".join([strain, filename + ".png"])))


def import_obo(filename):
    '''import the information of obo file to dict'''
    obos = []
    start = False
    with open(filename, "r") as o_h:
        for line in o_h:
            line = line.strip()
            if line == "[Term]":
                obo = {}
                start = True
            elif start:
                if len(line) == 0:
                    obos.append(obo.copy())
                    start = False
                else:
                    datas = line.split(": ")
                    if datas[0] == "is_a":
                        if "is_a" not in obo.keys():
                            obo["is_a"] = []
                        obo["is_a"].append(datas[1].strip())
                    else:
                        obo[datas[0]] = datas[1].strip()
    return obos


def import_class(slim_obo, classes, strain):
    if slim_obo["name"] not in classes[strain][slim_obo["namespace"]]:
        classes[strain][slim_obo["namespace"]][slim_obo["name"]] = 0
    classes[strain][slim_obo["namespace"]][slim_obo["name"]] += 1


def import_total(slim_obo, total_nums, strain):
    total_nums[strain][slim_obo["namespace"]] += 1
    total_nums[strain]["total"] += 1


def print_file(classes, total_nums, out_folder, stat):
    out_stat = open(stat, "w")
    printed = True
    for strain, datas in classes.items():
        if (strain == "All_genome") and len(classes) <= 2:
            printed = False
        if (printed) and (total_nums[strain]["total"] != 0):
            plot(total_nums[strain], strain, "three_roots",
                 total_nums[strain]["total"], out_folder)
            out_stat.write("{0}:\n".format(strain))
            for origin, types in datas.items():
                plot(types, strain, origin,
                     total_nums[strain][origin], out_folder)
                out_stat.write("\t{0}: {1}(percentage in total: {2})\n".format(
                               origin, total_nums[strain][origin],
                               float(total_nums[strain][origin]) /
                               float(total_nums[strain]["total"])))
                for type_, num in types.items():
                    out_stat.write("\t\t{0}: {1}(percentage "
                                   "in {2}: {3})\n".format(
                                       type_, num, origin, float(num) /
                                       float(total_nums[strain][origin])))
        else:
            printed = True
    out_stat.close()


def initiate_dict(classes, total_nums, index):
    classes[index] = {"biological_process": {},
                      "cellular_component": {},
                      "molecular_function": {}}
    total_nums[index] = {"biological_process": 0,
                         "cellular_component": 0,
                         "molecular_function": 0,
                         "total": 0}


def compare_go_slim(gos, term_obos, slim_obos, classes, total_nums):
    '''Compare GO term and GOslim database'''
    detect = False
    for strain, pros in gos.items():
        for pro, go_ids in pros.items():
            pro_list = []
            for go_id in go_ids:
                target_terms = [go_id]
                for target_term in target_terms:
                    for term_obo in term_obos:
                        if target_term == term_obo["id"]:
                            if "is_a" in term_obo.keys():
                                for is_a in term_obo["is_a"]:
                                    go_a = is_a.split(" ! ")
                                    if (go_a[1] != "biological_process") and (
                                        go_a[1] != "cellular_component") and (
                                            go_a[1] != "molecular_function"):
                                        target_terms.append(go_a[0])
                            elif ("is_obsolete" in term_obo.keys()):
                                if term_obo["is_obsolete"] == "true":
                                    break
                            for slim_obo in slim_obos:
                                for target_term in target_terms:
                                    if (target_term == slim_obo["id"]) and (
                                            target_term not in pro_list):
                                        detect = True
                                        import_class(slim_obo, classes, strain)
                                        import_class(slim_obo, classes,
                                                     "All_genome")
                                        import_total(slim_obo, total_nums,
                                                     strain)
                                        import_total(slim_obo, total_nums,
                                                     "All_genome")
                                        pro_list.append(target_term)
                            break
                    if detect:
                        detect = False
                        break


def map2goslim(slim_file, term_file, go_table, stat, out_folder):
    '''For mapping the GO to GOslim'''
    gos = {}
    classes = {}
    total_nums = {}
    initiate_dict(classes, total_nums, "All_genome")
    pre_strain = ""
    g_h = open(go_table, "r")
    print("Loading go table")
    for row in csv.reader(g_h, delimiter="\t"):
        if row[0] != "Genome":
            if row[0] != pre_strain:
                gos[row[0]] = {}
                initiate_dict(classes, total_nums, row[0])
            go_terms = row[-1].split("; ")
            gos[row[0]]["\t".join(row[1:4])] = go_terms
            pre_strain = row[0]
    print("Loading obo file")
    term_obos = import_obo(term_file)
    slim_obos = import_obo(slim_file)
    print("Starting mapping")
    compare_go_slim(gos, term_obos, slim_obos, classes, total_nums)
    print("Doing statistics and ploting")
    print_file(classes, total_nums, out_folder, stat)
    g_h.close()
