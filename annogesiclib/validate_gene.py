from annogesiclib.gff3 import Gff3Parser


def print_stat(feature, num, out):
    if num["_".join(["all", feature])] != 0:
        out.write("The number of {0} which is start "
                  "from TSS: {1} ({2})\n".format(
                      feature, num[feature],
                      float(num[feature]) / float(
                          num["_".join(["all", feature])])))
    else:
        out.write("The number of {0} which is start "
                  "from TSS: {1} ({2})\n".format(
                      feature, num[feature], "NA"))


def read_gff(gff_file, tss_file):
    tsss = []
    gffs = []
    gff_parser = Gff3Parser()
    fh = open(gff_file)
    for gff in gff_parser.entries(fh):
        gffs.append(gff)
    gffs = sorted(gffs, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    fh.close()
    tss_f = open(tss_file, "r")
    for tss in gff_parser.entries(tss_f):
        tsss.append(tss)
    tsss = sorted(tsss, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    tss_f.close()
    return gffs, tsss


def compare_tss(tsss, gff, utr_length, num_all, num_strain, program):
    detect = False
    for tss in tsss:
        length = utr_length
        if (gff.feature == "CDS") or (
                gff.feature == "rRNA") or (
                gff.feature == "tRNA"):
            if (gff.seq_id == tss.seq_id) and (
                    gff.start < tss.start) and (
                    gff.strand == "+") and (tss.strand == "+"):
                break
            elif (gff.seq_id == tss.seq_id) and (
                    gff.end < tss.start - utr_length) and (
                    gff.strand == "-") and (tss.strand == "-"):
                break
            if (gff.seq_id == tss.seq_id) and (
                    gff.strand == "+") and (tss.strand == "+") and (
                    gff.start - tss.start <= utr_length) and (
                    gff.start - tss.start >= 0):
                detect = True
                if (gff.start - tss.start) <= length:
                    start = tss
                    length = (gff.start - tss.start)
            elif (gff.seq_id == tss.seq_id) and (
                    gff.strand == "-") and (tss.strand == "-") and (
                    tss.start - gff.end <= utr_length) and (
                    tss.start - gff.end >= 0):
                detect = True
                if (tss.start - gff.end) <= length:
                    start = tss
                    length = (tss.start - gff.end)
    if program == "tss":
        type_ = "TSS"
    elif program == "processing":
        type_ = "Cleavage"
    if detect:
        detect = False
        gff.attributes["start_" + type_] = (
                type_ + "_" + str(start.start) + start.strand)
        if gff.feature == "CDS":
            num_all["cds"] += 1
            num_strain[gff.seq_id]["cds"] += 1
        elif gff.feature == "tRNA":
            num_all["tRNA"] += 1
            num_strain[gff.seq_id]["tRNA"] += 1
        elif gff.feature == "rRNA":
            num_all["rRNA"] += 1
            num_strain[gff.seq_id]["rRNA"] += 1


def print_file(gffs, out_cds_file, stat_file, num_all, num_strain):
    out_cds = open(out_cds_file, "w")
    out_cds.write("##gff-version 3\n")
    for gff in gffs:
        attribute_string = ";".join(
            ["=".join(items) for items in gff.attributes.items()])
        out_cds.write("\t".join([str(field) for field in [
                        gff.seq_id, gff.source, gff.feature, gff.start,
                        gff.end, gff.score, gff.strand, gff.phase,
                        attribute_string]]) + "\n")
    out = open(stat_file, "w")
    out.write("All genomes:\n")
    print_stat("cds", num_all, out)
    print_stat("tRNA", num_all, out)
    print_stat("rRNA", num_all, out)
    if len(num_strain) > 1:
        for strain in num_strain.keys():
            out.write(strain + ":\n")
            print_stat("cds", num_strain[strain], out)
            print_stat("tRNA", num_strain[strain], out)
            print_stat("rRNA", num_strain[strain], out)
    out_cds.close()
    out.close()


def validate_gff(tss_file, gff_file, stat_file, out_cds_file, utr_length,
                 program):
    num_all = {"all_cds": 0, "all_tRNA": 0, "all_rRNA": 0,
               "cds": 0, "tRNA": 0, "rRNA": 0}
    num_strain = {}
    pre_seq_id = ""
    gffs, tsss = read_gff(gff_file, tss_file)
    for gff in gffs:
        if gff.seq_id != pre_seq_id:
            num_strain[gff.seq_id] = {"all_cds": 0, "all_tRNA": 0,
                                      "all_rRNA": 0, "cds": 0,
                                      "tRNA": 0, "rRNA": 0}
            pre_seq_id = gff.seq_id
        if gff.feature == "CDS":
            num_all["all_cds"] += 1
            num_strain[gff.seq_id]["all_cds"] += 1
        elif gff.feature == "tRNA":
            num_all["all_tRNA"] += 1
            num_strain[gff.seq_id]["all_tRNA"] += 1
        elif gff.feature == "rRNA":
            num_all["all_rRNA"] += 1
            num_strain[gff.seq_id]["all_rRNA"] += 1
        compare_tss(tsss, gff, utr_length, num_all, num_strain, program)
    print_file(gffs, out_cds_file, stat_file, num_all, num_strain)
