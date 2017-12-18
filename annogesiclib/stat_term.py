import csv
from annogesiclib.gff3 import Gff3Parser


def plus_num(nums, strain, feature):
    nums[strain][feature] += 1
    nums["total"][feature] += 1


def print_percent(out, total, fract, name):
    if total != 0:
        out.write("\t\t(percentage of total {0}terminators = {1})\n".format(
                  name, (fract / total)))
    else:
        out.write("\t\t(percentage of total {0}terminators = 0.0)\n".format(
                  name))


def print_express(out, total, fract, name):
    if total != 0:
        out.write("\t\t(percentage of total {0}terminators which have "
                  "gene expression = {1})\n".format(name, (fract / total)))
    else:
        out.write("\t\t(percentage of total {0}terminators which have "
                  "gene expression = 0.0)\n".format(name))


def print_decrease(out, total, fract, name):
    if total != 0:
        out.write("\t\t(percentage of total {0}terminators "
                  "which have dramatic coverage decreasing = {1})\n".format(
                      name, (fract / total)))
    else:
        out.write("\t\t(percentage of total {0}terminators "
                  "which have dramatic coverage decreasing = 0.0)\n".format(
                      name))


def print_method(nums, method_name, method, express, detect, only, out):
    out.write(method_name + "\n")
    out.write("\tTotal {0} terminators = {1}\n".format(
              method_name, nums[method] + nums["frhp"]))
    print_percent(out, float(nums["total"]),
                  float(nums[method] + nums["frhp"]), "")
    out.write("\tTotal terminators which only can be "
              "detected in {0} = {1}\n".format(method_name, nums[method]))
    print_percent(out, float(nums["total"]), float(nums[method]), "")
    print_percent(out, float(nums[method] + nums["frhp"]),
                  float(nums[method]), method_name + " ")
    out.write("\tTotal {0} terminators which located in gene "
              "expression region = {1}\n".format(method_name, nums[express]))
    print_percent(out, float(nums["total"]), float(nums[express]), "")
    print_percent(out, float(nums[method] + nums["frhp"]),
                  float(nums[express]), method_name + " ")
    out.write("\tTotal {0} terminators which have dramatic coverage "
              "decreasing = {1}\n".format(method_name, nums[detect]))
    print_percent(out, float(nums["total"]), float(nums[detect]), "")
    print_percent(out, float(nums[method] + nums["frhp"]),
                  float(nums[detect]), method_name + " ")
    print_express(out, float(nums["total_ex"]), float(nums[detect]), "")
    print_express(out, float(nums[express]), float(nums[detect]),
                  method_name + " ")
    out.write("\tTotal terminators which have dramatic coverage decreasing"
              "(unique in {0}) = {1}\n".format(method_name, nums[only]))
    print_decrease(out, float(nums["total_de"]), float(nums[only]), "")
    print_decrease(out, float(nums[detect]), float(nums[only]),
                   method_name + " ")
    out.write("\n")


def print_intersection_number(out, nums, type_):
    print_percent(out, float(nums["total"]), float(nums[type_]), "")
    print_percent(out, float(nums["fr"]), float(nums[type_]), "method_1 ")
    print_percent(out, float(nums["hp"]), float(nums[type_]), "method_2 ")


def print_intersection_express(out, nums, type_):
    print_express(out, float(nums["total_ex"]), float(nums[type_]), "")
    print_express(out, float(nums["ex_fr"]), float(nums[type_]), "method_1 ")
    print_express(out, float(nums["ex_hp"]), float(nums[type_]), "method_2 ")


def print_file(nums, out, strain):
    out.write(strain + ":\n")
    out.write("Combine two methods:\n")
    out.write("\tTotal terminators = {0}\n".format(nums["total"]))
    out.write("\tTotal terminators which located in gene expression "
              "region = {0}\n".format(nums["total_ex"]))
    print_percent(out, float(nums["total"]), float(nums["total_ex"]), "")
    out.write("\tTotal terminators which have dramatic coverage "
              "decreasing = {0}\n".format(nums["total_de"]))
    print_percent(out, float(nums["total"]), float(nums["total_de"]), "")
    print_express(out, float(nums["total_ex"]), float(nums["total_de"]), "")
    out.write("\n")
    print_method(nums, "method_1", "fr", "ex_fr", "de_fr", "only_de_fr", out)
    print_method(nums, "method_2", "hp", "ex_hp", "de_hp", "only_de_hp", out)
    out.write("intersection two methods:\n")
    out.write("\tTotal terminators which overlap with "
              "two methods = {0}\n".format(nums["frhp"]))
    print_intersection_number(out, nums, "frhp")
    out.write("\tTotal overlaped terminators which located in "
              "gene expression region = {0}\n".format(nums["ex_frhp"]))
    print_intersection_number(out, nums, "ex_frhp")
    print_intersection_express(out, nums, "ex_frhp")
    out.write("\tTotal overlaped terminators which have dramatic "
              "coverage decreasing = {0}\n".format(nums["de_frhp"]))
    print_intersection_number(out, nums, "de_frhp")
    print_intersection_express(out, nums, "de_frhp")
    print_express(out, float(nums["total_de"]), float(nums["de_frhp"]), "")
    print_express(out, float(nums["de_fr"]),
                  float(nums["de_frhp"]), "method_1 ")
    print_express(out, float(nums["de_hp"]),
                  float(nums["de_frhp"]), "method_2 ")


def classify_terms(terms, nums, out_d, out_e, out_n, pre_strain):
    for term in terms:
        if term.seq_id != pre_strain:
            pre_strain = term.seq_id
            strain = term.seq_id
            nums[strain] = {
                    "fr": 0, "hp": 0, "frhp": 0, "ex_fr": 0, "ex_hp": 0,
                    "ex_frhp": 0, "de_fr": 0, "de_hp": 0, "de_frhp": 0,
                    "total": 0, "total_de": 0, "total_ex": 0, "only_de_fr": 0,
                    "only_de_hp": 0, "only_ex_fr": 0, "only_ex_hp": 0,
                    "de_frhp": 0, "ex_frhp": 0}
        if term.attributes["coverage_decrease"] == "True":
            out_d.write(term.info + "\n")
        if term.attributes["express"] == "True":
            out_e.write(term.info + "\n")
        if term.attributes["express"] != "True":
            out_n.write(term.info + "\n")
        if term.attributes["method"] == "gene_converged":
            plus_num(nums, strain, "total")
            plus_num(nums, strain, "fr")
            if term.attributes["coverage_decrease"] == "True":
                plus_num(nums, strain, "de_fr")
                plus_num(nums, strain, "only_de_fr")
                plus_num(nums, strain, "total_de")
            if term.attributes["express"] == "True":
                plus_num(nums, strain, "ex_fr")
                plus_num(nums, strain, "only_ex_fr")
                plus_num(nums, strain, "total_ex")
        elif term.attributes["method"] == "TransTermHP":
            plus_num(nums, strain,  "total")
            plus_num(nums, strain,  "hp")
            if term.attributes["coverage_decrease"] == "True":
                plus_num(nums, strain, "de_hp")
                plus_num(nums, strain, "only_de_hp")
                plus_num(nums, strain, "total_de")
            if term.attributes["express"] == "True":
                plus_num(nums, strain, "ex_hp")
                plus_num(nums, strain, "only_ex_hp")
                plus_num(nums, strain, "total_ex")
        elif term.attributes["method"] == "gene_converged,TransTermHP":
            plus_num(nums, strain, "total")
            plus_num(nums, strain, "frhp")
            if term.attributes["coverage_decrease"] == "True":
                plus_num(nums, strain, "de_frhp")
                plus_num(nums, strain, "de_fr")
                plus_num(nums, strain, "de_hp")
                plus_num(nums, strain, "total_de")
            if term.attributes["express"] == "True":
                plus_num(nums, strain, "ex_frhp")
                plus_num(nums, strain, "ex_fr")
                plus_num(nums, strain, "ex_hp")
                plus_num(nums, strain, "total_ex")


def check_repeat(checks, strain, start, end, strand):
    detect = False
    try:
        term = {"strain": strain, "start": int(start),
                "strand": strand, "end": int(end)}
        if len(checks) == 0:
            checks.append(term)
            detect = True
        else:
            if term not in checks:
                detect = True
                checks.append(term)
        return detect
    except ValueError:
        return detect


def stat_term(term_gff, term_table, stat, output_decrease,
              output_expression, output_non):
    terms = []
    nums = {}
    nums["total"] = {
            "fr": 0, "hp": 0, "frhp": 0, "ex_fr": 0, "ex_hp": 0, "ex_frhp": 0,
            "de_fr": 0, "de_hp": 0, "de_frhp": 0, "total": 0, "total_de": 0,
            "total_ex": 0, "only_de_fr": 0, "only_de_hp": 0, "only_ex_fr": 0,
            "only_ex_hp": 0, "de_frhp": 0, "ex_frhp": 0}
    pre_strain = ""
    out_te = open(output_expression + ".csv", "w")
    out_td = open(output_decrease + ".csv", "w")
    out_tn = open(output_non + ".csv", "w")
    fh = open(term_table, "r")
    out_tn.write("\t".join(["Genome", "Name", "Start", "End", "Strand",
                            "Method", "Associated_gene", "Associated_transcript",
                            "Coverage_decrease", "Coverage_detail"]) + "\n")
    gh = open(term_gff)
    checks = []
    for entry in Gff3Parser().entries(gh):
        detect = check_repeat(checks, entry.seq_id, entry.start,
                              entry.end, entry.strand)
        if detect:
            terms.append(entry)
    checks = []
    for row in csv.reader(fh, delimiter="\t"):
        detect = check_repeat(checks, row[0], row[2], row[3], row[4])
        if detect:
            if (row[-1] != "NA") and (row[-1] != "No_coverage_decreasing"):
                out_td.write("\t".join(row) + "\n")
                out_te.write("\t".join(row) + "\n")
            if (row[-1] == "No_coverage_decreasing"):
                out_te.write("\t".join(row) + "\n")
            if (row[-1] == "NA"):
                out_tn.write("\t".join(row) + "\n")
    out = open(stat, "w")
    out_e = open(output_expression + ".gff", "w")
    out_d = open(output_decrease + ".gff", "w")
    out_n = open(output_non + ".gff", "w")
    out_e.write("##gff-version 3\n")
    out_d.write("##gff-version 3\n")
    out_n.write("##gff-version 3\n")
    classify_terms(terms, nums, out_d, out_e, out_n, pre_strain)
    out.write("method_1 is searching the gene converged region.\n")
    out.write("method_2 is TransTermHP.\n")
    if len(nums) > 2:
        print_file(nums["total"], out, "All genome")
    else:
        for strain, datas in nums.items():
            if strain != "total":
                print_file(datas, out, strain)
    out_te.close()
    out_td.close()
    out_tn.close()
    out.close()
    out_e.close()
    out_d.close()
    out_n.close()
    fh.close()
    gh.close()
