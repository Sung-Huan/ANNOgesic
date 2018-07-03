from annogesiclib.gff3 import Gff3Parser


def create_dict(nums, strain, utr_detect):
    nums[strain] = {}
    if utr_detect:
        types = ["all", "5'UTR_derived", "3'UTR_derived",
                 "interCDS", "intergenic", "antisense"]
    else:
        types = ["all", "intergenic", "antisense"]
    for type_ in types:
        nums[strain][type_] = {}
        for feature in ["TSS", "sRNA", "all", "RBS", "TSS_RBS",
                        "TSS_sRNA_RBS", "TSS_sRNA", "RBS_sRNA"]:
            nums[strain][type_][feature] = 0
    return nums


def plus_data(nums, strain, sorf_types, features, utr_detect):
    for sorf_type in sorf_types:
        if ((not utr_detect) and (
                (sorf_type == "intergenic") or (
                sorf_type == "antisense") or (
                sorf_type == "all"))) or (
                utr_detect):
            for feature in features:
                nums[strain][sorf_type][feature] += 1


def print_num(out, num, nums, strain, type_):
    out.write("(for genome {0}; ".format(
              float(num) / float(nums[strain]["all"]["all"])))
    if nums[strain][type_]["all"] == 0:
        out.write("for {0} - {1})\n".format(
                  type_, 0))
    else:
        out.write("for {0} - {1})\n".format(
                  type_, float(num) / float(nums[strain][type_]["all"])))


def print_stat(nums, nums_best, strain, out, utr_detect):
    out.write(strain + ":\n")
    if utr_detect:
        out.write("\ttotal sORF in this genome are {0}\n".format(
                  nums[strain]["all"]["all"]))
    for type_, features in nums[strain].items():
        out.write("\ttotal sORF of {0} sORF candidates are {1}".format(
                  type_, nums[strain][type_]["all"]))
        out.write("(for this genome - {0})\n".format(
                  float(nums[strain][type_]["all"]) /
                  float(nums[strain]["all"]["all"])))
        for feature, num in features.items():
            if feature == "TSS":
                out.write("\t\ttotal sORF which start "
                          "from TSS are {0}".format(num))
                print_num(out, num, nums, strain, type_)
            elif feature == "sRNA":
                out.write("\t\ttotal sORF without overlap with "
                          "sRNA candidates are {0}".format(num))
                print_num(out, num, nums, strain, type_)
            elif feature == "RBS":
                out.write("\t\ttotal sORF which related with "
                          "ribosomal binding site are {0}".format(num))
                print_num(out, num, nums, strain, type_)
            elif feature == "TSS_RBS":
                out.write("\t\ttotal sORF which start from TSS and related "
                          "with ribosomal binding site are {0}".format(num))
                print_num(out, num, nums, strain, type_)
            elif feature == "TSS_sRNA":
                out.write("\t\ttotal sORF which start from TSS and without "
                          "overlap with sRNA candidates are {0}".format(num))
                print_num(out, num, nums, strain, type_)
            elif feature == "RBS_sRNA":
                out.write("\t\ttotal sORF which related with "
                          "ribosomal binding site and ")
                out.write("without overlap with "
                          "sRNA candidates are {0}".format(num))
                print_num(out, num, nums, strain, type_)
            elif feature == "TSS_RBS_sRNA":
                out.write("\t\ttotal sORF which start from TSS and "
                          "related with ribosomal binding site and ")
                out.write("without overlap with "
                          "sRNA candidates are {0}".format(num))
                print_num(out, num, nums, strain, type_)
        if strain in nums_best.keys():
            out.write("\t\tThe best sORF are {0}\n".format(
                nums_best[strain][type_]["all"]))
            out.write("\t\tThe best sORF which without overlap with "
                      "sRNA are {0}".format(nums_best[strain][type_]["sRNA"]))
            print_num(out, nums_best[strain][type_]["sRNA"],
                      nums_best, strain, type_)
        else:
            out.write("\t\tThe best sORF are 0\n")
            out.write("\t\tThe best sORF which without overlap with "
                      "sRNA are 0\n")


def read_file(sorf_gff):
    sorfs = []
    fh = open(sorf_gff)
    for entry in Gff3Parser().entries(fh):
        sorfs.append(entry)
    sorfs = sorted(sorfs, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    fh.close()
    return sorfs


def get_stat_num(sorfs_all, utr_detect):
    strain = ""
    nums = {}
    create_dict(nums, "total", utr_detect)
    for sorf in sorfs_all:
        if strain != sorf.seq_id:
            create_dict(nums, sorf.seq_id, utr_detect)
            strain = sorf.seq_id
        if sorf.attributes["sORF_type"] == "intergenic":
            sorf_type = "intergenic"
        elif sorf.attributes["sORF_type"] == "antisense":
            sorf_type = "antisense"
        else:
            if "5utr" in sorf.attributes["sORF_type"]:
                sorf_type = "5'UTR_derived"
            elif "3utr" in sorf.attributes["sORF_type"]:
                sorf_type = "3'UTR_derived"
            elif "interCDS" in sorf.attributes["sORF_type"]:
                sorf_type = "interCDS"
        check_class(sorf, nums, sorf_type, utr_detect, strain)
    return nums


def check_class(sorf, nums, sorf_type, utr_detect, strain):
    if (sorf.attributes["with_TSS"] != "NA") and \
       (sorf.attributes["sRNA"] == "NA") and \
       (sorf.attributes["rbs"] != "NA"):
        plus_data(nums, "total", [sorf_type, "all"],
                  ["all", "TSS", "sRNA", "RBS", "TSS_RBS",
                   "TSS_sRNA_RBS", "TSS_sRNA", "RBS_sRNA"], utr_detect)
        plus_data(nums, strain, [sorf_type, "all"],
                  ["all", "TSS", "sRNA", "RBS", "TSS_RBS",
                   "TSS_sRNA_RBS", "TSS_sRNA", "RBS_sRNA"], utr_detect)
    elif (sorf.attributes["with_TSS"] != "NA") and \
         (sorf.attributes["sRNA"] == "NA"):
        plus_data(nums, "total", [sorf_type, "all"],
                  ["all", "TSS", "sRNA", "TSS_sRNA"], utr_detect)
        plus_data(nums, strain, [sorf_type, "all"],
                  ["all", "TSS", "sRNA", "TSS_sRNA"], utr_detect)
    elif (sorf.attributes["with_TSS"] != "NA") and \
         (sorf.attributes["rbs"] != "NA"):
        plus_data(nums, "total", [sorf_type, "all"],
                  ["all", "TSS", "RBS", "TSS_RBS"], utr_detect)
        plus_data(nums, strain, [sorf_type, "all"],
                  ["all", "TSS", "RBS", "TSS_RBS"], utr_detect)
    elif (sorf.attributes["rbs"] != "NA") and \
         (sorf.attributes["sRNA"] == "NA"):
        plus_data(nums, "total", [sorf_type, "all"],
                  ["all", "RBS", "sRNA", "RBS_sRNA"], utr_detect)
        plus_data(nums, strain, [sorf_type, "all"],
                  ["all", "RBS", "sRNA", "RBS_sRNA"], utr_detect)
    elif sorf.attributes["with_TSS"] != "NA":
        plus_data(nums, "total", [sorf_type, "all"],
                  ["all", "TSS"], utr_detect)
        plus_data(nums, strain, [sorf_type, "all"],
                  ["all", "TSS"], utr_detect)
    elif sorf.attributes["sRNA"] == "NA":
        plus_data(nums, "total", [sorf_type, "all"],
                  ["all", "sRNA"], utr_detect)
        plus_data(nums, strain, [sorf_type, "all"],
                  ["all", "sRNA"], utr_detect)
    elif sorf.attributes["rbs"] != "NA":
        plus_data(nums, "total", [sorf_type, "all"],
                  ["all", "RBS"], utr_detect)
        plus_data(nums, strain, [sorf_type, "all"],
                  ["all", "RBS"], utr_detect)
    else:
        plus_data(nums, "total", [sorf_type, "all"], ["all"], utr_detect)
        plus_data(nums, strain, [sorf_type, "all"], ["all"], utr_detect)


def stat(sorf_all, sorf_best, stat_file, utr_detect):
    sorfs_all = read_file(sorf_all)
    sorfs_best = read_file(sorf_best)
    nums = get_stat_num(sorfs_all, utr_detect)
    nums_best = get_stat_num(sorfs_best, utr_detect)
    out = open(stat_file, "w")
    out.write("The filtering condition for the best sORF: \n")
    out.write("1. If TSS file exists, it will select the "
              "sORF which start from TSS.\n")
    out.write("2. If TSS file exists, it will select the "
              "sORF which have a ribosomal binding site ")
    out.write("and the ribosomal binding site shoule after a TSS.\n")
    out.write("3. If sRNA file exists and you want to "
              "exclude sORF which overlap with sRNA, ")
    out.write("it will select sORF which have non-overlap with sRNA.\n\n")
    if len(nums) <= 2:
        for strain in nums.keys():
            if strain != "total":
                print_stat(nums, nums_best, strain, out, utr_detect)
    else:
        for strain in nums.keys():
            print_stat(nums, nums_best, strain, out, utr_detect)
    out.close()
