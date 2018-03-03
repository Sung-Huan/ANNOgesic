import csv
import itertools


def _boolean(data):
    if data == "False":
        result = False
    else:
        result = True
    return result


def row_to_location(row):
    if row[4] == "0":
        sub = False
        nosub = True
    else:
        sub = True
        nosub = False
    tss = _boolean(row[6])
    term = _boolean(row[8])
    return {"have no sub-operons": nosub, "have sub-operons": sub,
            "start with tss": tss, "stop with terminator": term}


def plus_num(num_total, strain, type_):
    num_total["total"][type_] += 1
    num_total[strain][type_] += 1
    num_total["total"]["total"] += 1
    num_total[strain]["total"] += 1


def print_stat(operons, total_num, class_operon, out):
    num_features = {}
    out.write("Total number of operons is {0}\n".format(total_num))
    out.write("The sub operon and features:\n")
    for operon in operons:
        for it in range(1, 5):
            for features in itertools.combinations(operon.keys(), it):
                check_key = 0
                for key in features:
                    if operon[key]:
                        if it == 1:
                            if key in num_features.keys():
                                num_features[key] += 1
                            else:
                                num_features[key] = 1
                        check_key += 1
                if (check_key == it) and (it != 1):
                    key = " and ".join(features)
                    if key in num_features.keys():
                        num_features[key] += 1
                    else:
                        num_features[key] = 1
    for key, value in num_features.items():
        out.write("\tthe number of operons which {0} = {1} ({2})\n".format(
                  key, value, float(value) / float(total_num)))
    out.write("mono/polycistronic:\n")
    out.write("\tmonocistronic: {0} ({1})\n".format(
              class_operon["mono"],
              float(class_operon["mono"]) / float(class_operon["total"])))
    out.write("\tpolycistronic: {0} ({1})\n".format(
              class_operon["poly"],
              float(class_operon["poly"]) / float(class_operon["total"])))


def stat(input_file, out_file):
    out = open(out_file, "w")
    operons = {}
    operons_all = []
    tmp_id = ""
    f_h = open(input_file, "r")
    pre_seq_id = ""
    total_num = {}
    total_num_all = 0
    class_operon = {}
    class_operon["total"] = {"na": 0, "mono": 0, "poly": 0, "total": 0}
    for row in csv.reader(f_h, delimiter="\t"):
        if row[0] != "Operon_ID":
            if row[0] != tmp_id:
                if pre_seq_id != row[1]:
                    pre_seq_id = row[1]
                    operons[row[1]] = []
                    total_num[row[1]] = 0
                    class_operon[row[1]] = {"na": 0, "mono": 0,
                                            "poly": 0, "total": 0}
                operons[row[1]].append(row_to_location(row))
                operons_all.append(row_to_location(row))
                total_num[row[1]] += 1
                total_num_all += 1
                if row[-1] == "NA":
                    plus_num(class_operon, row[1], "na")
                elif len(row[-1].split(",")) == 1:
                    plus_num(class_operon, row[1], "mono")
                elif len(row[-1].split(",")) > 1:
                    plus_num(class_operon, row[1], "poly")
                tmp_id = row[0]
    if len(operons) > 1:
        out.write("All genomes:\n")
        print_stat(operons_all, total_num_all, class_operon["total"], out)
    for strain in operons.keys():
        out.write("\n" + strain + ":\n")
        print_stat(operons[strain], total_num[strain],
                   class_operon[strain], out)
    out.close()
    f_h.close()
