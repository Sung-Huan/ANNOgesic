from annogesiclib.splice_parser import SpliceParser
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.helper import Helper


def get_feature(cds):
    '''get proper feature name'''
    if "locus_tag" in cds.attributes.keys():
        feature = cds.attributes["locus_tag"]
    elif "protein_id" in cds.attributes.keys():
        feature = cds.attributes["protein_id"]
    elif "ID" in cds.attributes.keys():
        strand = Helper().get_strand_name(cds.strand)
        feature = "".join([cds.attributes["ID"], ":",
                           str(cds.start), "-", str(cds.end),
                           "_", strand])
    else:
        strand = Helper().get_strand_name(cds.strand)
        feature = "".join([cds.feature, ":",
                           str(cds.start), "-", str(cds.end),
                           "_", strand])
    return feature


def detect_conflict(gffs, circ, num, out, out_best, args_circ):
    '''remove the false positive which overlap with known annotation'''
    detect = False
    gff = None
    for gff in gffs:
        if (gff.seq_id == circ.strain) and (
                gff.strand == circ.strand):
            if ((gff.start <= circ.start) and (
                 gff.end >= circ.start) and (
                 gff.end <= circ.end)) or (
                (gff.start >= circ.start) and (
                 gff.end <= circ.end)) or (
                (gff.start >= circ.start) and (
                 gff.start <= circ.end) and (
                 gff.end >= circ.end)) and (
                (gff.start <= circ.start) and (
                 gff.end >= circ.end)):
                detect = True
                break
    if detect:
        feature = get_feature(gff)
        out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                  circ.strain, circ.strand,
                  circ.start, circ.end, feature, circ.supported_reads,
                  float(circ.supported_reads) / float(circ.start_site_reads),
                  float(circ.supported_reads) / float(circ.end_site_reads)))
    else:
        start_read = float(circ.supported_reads) / float(circ.start_site_reads)
        end_read = float(circ.supported_reads) / float(circ.end_site_reads)
        if (circ.supported_reads >= args_circ.support) and (
                start_read >= args_circ.start_ratio) and (
                end_read >= args_circ.end_ratio):
            out_best.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}"
                           "\t{6}\t{7}\n".format(
                            circ.strain,
                            circ.strand, circ.start, circ.end, "NA",
                            circ.supported_reads, start_read, end_read))
        out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                  circ.strain, circ.strand,
                  circ.start, circ.end, "NA", circ.supported_reads,
                  start_read, end_read))
    return detect


def import_num(support, nums, strain):
    if support not in nums[strain].keys():
        nums[strain][support] = 0
    if support not in nums["all"].keys():
        nums["all"][support] = 0
    nums[strain][support] += 1
    nums["all"][support] += 1


def print_file(nums, stat, strain):
    for key in sorted(nums[strain].keys()):
        stat.write("\tthe number of potential circular RNAs, ")
        stat.write("more than {0} supported it = {1}\n".format(
                   key, nums[strain][key]))


def read_file(input_file, gff_file, hypo):
    circs = []
    gffs = []
    ps = SpliceParser()
    high = 0
    splice_fh = open(input_file)
    for entry in ps.parser(splice_fh):
        if entry.supported_reads > high:
            high = entry.supported_reads
        circs.append(entry)
    gff_parser = Gff3Parser()
    for entry in gff_parser.entries(open(gff_file)):
        if ("product" in entry.attributes.keys()) and (hypo):
            if "hypothetical protein" not in entry.attributes["product"]:
                gffs.append(entry)
        else:
            gffs.append(entry)
    gffs = sorted(gffs, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    circs = sorted(circs, key=lambda x: (x.strain, x.supported_reads),
                   reverse=True)
    splice_fh.close()
    return circs, gffs, high


def get_circrna(circs, gffs, high, out, out_best, args_circ):
    '''search the splice data to find the potential circRNA'''
    num_circular = {}
    num_circular["all"] = 0
    num_support = {}
    num_support["all"] = {}
    num_conflict = {}
    num_conflict["all"] = {}
    pre_seq_id = ""
    num = 0
    for circ in circs:
        if pre_seq_id != circ.strain:
            num_support[circ.strain] = {}
            num_conflict[circ.strain] = {}
            num_circular[circ.strain] = 0
        if (circ.situation != "F") and (circ.splice_type == "C"):
            num_circular[circ.strain] += 1
            num_circular["all"] += 1
            detect = detect_conflict(gffs, circ, num, out, out_best, args_circ)
            for support in range(0, high + 5, 5):
                if circ.supported_reads >= int(support):
                    import_num(support, num_support, circ.strain)
                    if detect is False:
                        if (float(circ.supported_reads) /
                            float(circ.start_site_reads) >=
                            args_circ.start_ratio) and (
                            float(circ.supported_reads) /
                            float(circ.end_site_reads) >=
                                args_circ.end_ratio):
                            import_num(support, num_conflict, circ.strain)
            num += 1
        pre_seq_id = circ.strain
    return {"circular": num_circular, "support": num_support,
            "conflict": num_conflict}


def detect_circrna(input_file, gff_file, output_file, args_circ, statistics):
    circs, gffs, high = read_file(input_file, gff_file, args_circ.hypo)
    out = open(output_file, "w")
    out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
               "Genome", "Strand", "Start", "End", "Annotation_overlap",
               "Supported_reads", "Supported_reads/Reads_at_start",
               "Supported_reads/Reads_at_end"))
    out_best = open(output_file.replace("_all.csv", "_best.csv"), "w")
    out_best.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
               "Genome", "Strand", "Start", "End", "Annotation_overlap",
               "Supported_reads", "Supported_reads/Reads_at_start",
               "Supported_reads/Reads_at_end"))
    nums = get_circrna(circs, gffs, high, out, out_best, args_circ)
    stat = open(statistics, "w")
    stat.write("All genomes:\n")
    stat.write("\tBefore filtering:\n")
    stat.write("\tthe number of all circular RNAs = {0}\n".format(
               nums["circular"]["all"]))
    print_file(nums["support"], stat, "all")
    stat.write("\n\tAfter filtering:\n")
    stat.write("\t\twithout conflict with annotation\n")
    stat.write("\t\tsupport read ratio of starting "
               "point is larger than {0}\n".format(
                   args_circ.start_ratio))
    stat.write("\t\tsupport read ratio of end point "
               "is larger than {0}\n".format(
                   args_circ.end_ratio))
    print_file(nums["conflict"], stat, "all")
    if len(nums["circular"]) > 2:
        for strain in nums["circular"].keys():
            if strain != "all":
                stat.write("\n{0}:\n".format(strain))
                stat.write("\tBefore filtering:\n")
                stat.write("\tthe number of all circular RNAs = {0}\n".format(
                               nums["circular"][strain]))
                print_file(nums["support"], stat, strain)
                stat.write("\n\tAfter filtering:\n")
                stat.write("\t\twithout conflict with annotation\n")
                stat.write("\t\tsupport read ratio of starting point "
                           "is larger than {0}\n".format(
                               args_circ.start_ratio))
                stat.write("\t\tsupport read ratio of end point "
                           "is larger than {0}\n".format(
                               args_circ.end_ratio))
                print_file(nums["conflict"], stat, strain)
    out.close()
    out_best.close()
    stat.close()
