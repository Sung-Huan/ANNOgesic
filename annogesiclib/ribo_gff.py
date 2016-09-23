import csv


def import_data(row):
    return{"strain": row[1], "strand": row[2],
           "associate": row[3], "start_seq": int(row[4]),
           "end_seq": int(row[5]), "rfam": row[6], "e": row[7],
           "start_align": int(row[8]), "end_align": int(row[9]),
           "info": "|".join(row[0:6]), "ID": row[0]}


def read_file(ribo_table, rfam_table):
    ribos = []
    rfams = []
    f_h = open(ribo_table, "r")
    for row in csv.reader(f_h, delimiter="\t"):
        if not row[0].startswith("#"):
            ribos.append(import_data(row))
    r_h = open(rfam_table, "r")
    for row in csv.reader(r_h, delimiter="\t"):
        rfams.append({"ID": row[0].strip(), "class": row[1].strip()})
    ribos = sorted(ribos, key=lambda x: (x["strain"], x["start_seq"]))
    f_h.close()
    r_h.close()
    return ribos, rfams


def get_overlap(pre_ribo, ribo, overlap, overlaps):
    if (pre_ribo["strain"] == ribo["strain"]) and \
       (pre_ribo["strand"] == ribo["strand"]) and \
       (pre_ribo["ID"] == ribo["ID"]):
        overlap = True
    if overlap:
        detect = False
        for over in overlaps[ribo["strain"]]:
            if pre_ribo["info"] in over:
                over = over + ";" + ribo["info"]
                detect = True
        if not detect:
            overlaps[ribo["strain"]].append(
                     pre_ribo["info"] + ";" + ribo["info"])


def print_gff(num, ribo, out, stats, strain, feature):
    attribute = ";".join(["=".join(items) for items in [
                          ("ID", "_".join([feature.lower(), str(num)])),
                          ("Name", ribo["rfam_name"]),
                          ("rfam_id", ribo["rfam"]),
                          ("e_value", ribo["e"]),
                          ("method", "infernal_to_Rfam")]])
    out.write("\t".join([str(field) for field in [
              ribo["strain"], "ANNOgesic", feature,
              str(ribo["start_seq"]), str(ribo["end_seq"]),
              ".", ribo["strand"], ".", attribute]]) + "\n")
    stats["total"]["total"] += 1
    stats[strain]["total"] += 1


def import_stat(rfams, ribo, stats, strain):
    for rfam in rfams:
        if ribo["rfam"] == rfam["ID"]:
            ribo["rfam_name"] = rfam["class"]
            if rfam["class"] not in stats["total"].keys():
                stats["total"][rfam["class"]] = 1
            else:
                stats["total"][rfam["class"]] += 1
            if rfam["class"] not in stats[strain].keys():
                stats[strain][rfam["class"]] = 1
            else:
                stats[strain][rfam["class"]] += 1


def print_number(stats, repeat, out, strain, feature):
    out.write("Total number of potential {0} are {1}\n".format(
               feature.replace("_", " "), stats[strain]["total"]))
    out.write("The number of potential {0} which "
              "have overlap region with others are {1}\n".format(
               feature.replace("_", " "), repeat,))
    out.write(feature + "_name\tnumbers\n")
    for type_, num in stats[strain].items():
        if type_ != "total":
            out.write("{0}\t{1}\n".format(type_, num))


def print_stat(stats, out_stat, overlaps, feature):
    out = open(out_stat, "w")
    print_file = False
    repeat = 0
    if len(stats) > 2:
        out.write("All strains:\n")
        print_file = True
        for strain, overs in overlaps.items():
            for over in overs:
                datas = over.split(";")
                repeat = repeat + len(datas)
        print_number(stats, repeat, out, "total", feature)
    for strain, datas in stats.items():
        repeat = 0
        if strain != "total":
            print_file = True
            out.write("{0}:\n".format(strain))
            for over in overlaps[strain]:
                datas = over.split(";")
                repeat = repeat + len(datas)
            print_number(stats, repeat, out, strain, feature)
            print_strain = strain
    if print_file:
        count = 1
        if len(stats) > 2:
            for strain, overs in overlaps.items():
                for over in overs:
                    datas = over.split(";")
                    out.write("\noverlap candidates set {0}:\n".format(count))
                    count += 1
                    for data in datas:
                        out.write("\t{0}\n".format(data))
        else:
            for over in overlaps[print_strain]:
                datas = over.split(";")
                out.write("\noverlap candidates set {0}:\n".format(count))
                count += 1
                for data in datas:
                    out.write("\t{0}\n".format(data))
    out.close()


def stat_and_covert2gff(ribo_table, rfam_table, gff_file, fuzzy, out_stat,
                        feature):
    '''do statistics and print gff file of riboswitch'''
    stats = {}
    overlaps = {}
    pre_strain = ""
    stats["total"] = {"total": 0}
    num = 0
    ribos, rfams = read_file(ribo_table, rfam_table)
    out = open(gff_file, "w")
    out.write("##gff-version 3\n")
    pre_gff = None
    for ribo in ribos:
        overlap = False
        if ribo["strain"] != pre_strain:
            overlaps[ribo["strain"]] = []
            first = True
            strain = ribo["strain"]
            pre_strain = ribo["strain"]
            stats[strain] = {"total": 0}
        if first:
            first = False
            pre_ribo = ribo
        else:
            get_overlap(pre_ribo, ribo, overlap, overlaps)
            pre_ribo = ribo
        if ribo["start_align"] > fuzzy:
            ribo["start_seq"] = ribo["start_seq"] + ribo["start_align"] - fuzzy
        if (ribo["end_seq"] - (ribo["start_seq"] + ribo["end_align"])) > fuzzy:
            ribo["end_seq"] = ribo["start_seq"] + ribo["end_align"] + fuzzy
        import_stat(rfams, ribo, stats, strain)
        if pre_gff is not None:
            if (pre_gff["strain"] == ribo["strain"]) and (
                    pre_gff["strand"] == ribo["strand"]) and (
                    pre_gff["start_seq"] == ribo["start_seq"]) and (
                    pre_gff["end_seq"] == ribo["end_seq"]):
                pre_gff["rfam_name"] = "/".join(
                        [pre_gff["rfam_name"], ribo["rfam_name"]])
                pre_gff["rfam"] = ",".join([pre_gff["rfam"], ribo["rfam"]])
                pre_gff["e"] = ",".join([pre_gff["e"], ribo["e"]])
            else:
                print_gff(num, pre_gff, out, stats, strain, feature)
                num += 1
                pre_gff = ribo
        else:
            pre_gff = ribo
    print_gff(num, pre_gff, out, stats, strain, feature)
    print_stat(stats, out_stat, overlaps, feature)
    out.close()
