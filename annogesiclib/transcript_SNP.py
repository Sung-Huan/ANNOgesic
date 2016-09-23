import csv
import copy
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')


def plot_bar(cutoffs, strain, out_snp, type_):
    name = []
    for index in range(0, len(cutoffs) + 1):
        name.append(index * 10)
    ind = np.arange(len(cutoffs))
    width = 0.5
    plt.figure(figsize=(20, 15))
    plt.bar(ind, cutoffs, width, color='#FF9999')
    plt.ylabel('the number of SNPs', fontsize=20)
    plt.xlabel('QUAL of SNP in transcripts', fontsize=20)
    plt.xlim([0, len(cutoffs) + 1])
    plt.xticks(ind+width-0.75, name, fontsize=18, rotation=40)
    plt.yticks(fontsize=18)
    plt.savefig("_".join([out_snp, strain, "SNP_QUAL", type_]))
    plt.clf()


def row_in_list(row):
    info = ""
    filt = "."
    snps = {"strain": None, "pos": -1, "id": "", "ref": "",
            "alt": "", "qual": -1, "filter": "", "info": "",
            "depth": -1, "all_info": "", "indel": -1, "frac": -1,
            "dp4_sum": -1, "dp4_frac": -1}
    if len(row) >= 8:
        snps = {"strain": row[0], "pos": int(row[1]), "id": row[2],
                "ref": row[3], "alt": row[4], "qual": float(row[5]),
                "filter": filt, "info": "", "depth": -1,
                "all_info": "\t".join(row), "indel": -1, "frac": -1}
        infos = row[7].split(";")
        snps["info"] = infos
        for info in infos:
            datas = info.split("=")
            if len(datas) > 1:
                if datas[0] == "DP":
                    snps["depth"] = int(datas[1])
                if datas[0] == "IDV":
                    snps["indel"] = int(datas[1])
                if datas[0] == "IMF":
                    snps["frac"] = float(datas[1])
                if datas[0] == "DP4":
                    dp4s = datas[1].split(",")
                    snps["dp4_sum"] = int(dp4s[2]) + int(dp4s[3])
                    total = (int(dp4s[0]) + int(dp4s[1]) + 
                             int(dp4s[2]) + int(dp4s[3]))
                    snps["dp4_frac"] = float(snps["dp4_sum"] / total)
        return snps
    else:
        return snps


def gen_ref(snps, pos, refs, num):
    if num == 1:
        for snp in snps:
            refs.append(":".join([str(pos), snp["alt"]]))
    else:
        new_refs = []
        for snp in snps:
            for ref in refs:
                new_refs.append(ref + "_" + str(pos) + ":" + snp["alt"])
        refs = copy.deepcopy(new_refs)
    return refs


def change(snp, seq):
    start_point = snp["pos"] - 1 + seq["num_mod"]
    end_point = snp["pos"] - 1 + len(snp["ref"]) + seq["num_mod"]
    if len(snp["ref"]) == len(snp["alt"]):
        if seq["seq"][start_point: end_point].upper() == snp["ref"].upper():
            seq["seq"] = (seq["seq"][:start_point] +
                          snp["alt"].lower() + seq["seq"][end_point:])
    if len(snp["ref"]) > len(snp["alt"]):
        if seq["seq"][start_point: end_point].upper() == snp["ref"].upper():
            seq["seq"] = (seq["seq"][:start_point] +
                          snp["alt"].lower() + seq["seq"][end_point:])
            seq["num_mod"] = seq["num_mod"] - (
                    len(snp["ref"]) - len(snp["alt"]))
    if len(snp["ref"]) < len(snp["alt"]):
        if seq["seq"][start_point: end_point].upper() == snp["ref"].upper():
            seq["seq"] = (seq["seq"][:start_point] +
                          snp["alt"].lower() + seq["seq"][end_point:])
            seq["num_mod"] = seq["num_mod"] - (
                    len(snp["ref"]) - len(snp["alt"]))


def get_n_a_value(para, depth_file, min_sample):
    '''get the corresponding number of cutoff'''
    tag = para.split("_")[0]
    value = float(para.split("_")[-1])
    if tag == "a":
        fh = open(depth_file, "r")
        total = 0
        num = 0
        for row in csv.reader(fh, delimiter='\t'):
            total = total + int(row[-1])
            num += 1
        avg_dep = float(total / num)
        fh.close()
        value = avg_dep * value
    elif tag == "n":
        value = (min_sample * value)
    elif tag == "r":
        pass
    return value

def apply_filter(filters, snp, snps):
    '''apply the filter to remove the false positives'''
    exclude = False
    for filt in filters:
        tag = filt.split("_")[0]
        cutoff = filt.split("_")[-1]
        for info in snp["info"]:
            if tag in info:
                value = float(info.split("=")[-1])
                if cutoff[0] == "b":
                    if value < float(cutoff[1:]):
                        exclude = True
                if cutoff[0] == "s":
                    if value > float(cutoff[1:]):
                        exclude = True
    if not exclude:
        snps.append(snp)

def import_data(snp_file, args_snp, bam_number, depth_file):
    snps = []
    raw_snps = []
    dess = []
    max_quals = {}
    max_quals["All_strain"] = 0
    pre_strain = ""
    cutoff_sum = get_n_a_value(args_snp.dp4_sum, depth_file,
                               args_snp.min_sample)
    cutoff_frac = float(args_snp.dp4_frac)
    depth_s = get_n_a_value(args_snp.depth_s, depth_file,
                            args_snp.min_sample)
    depth_b = get_n_a_value(args_snp.depth_b, depth_file,
                            args_snp.min_sample)
    idv = get_n_a_value(args_snp.idv, depth_file,
                        args_snp.min_sample)
    imf = float(args_snp.imf)
    fh = open(snp_file, "r")
    for row in csv.reader(fh, delimiter="\t"):
        if row[0].startswith("##"):
            dess.append(row[0])
        if row[0].startswith("#"):
            continue
        else:
            snp = row_in_list(row)
            raw_snps.append(snp)
            if snp["strain"]:
                if snp["strain"] != pre_strain:
                    pre_strain = snp["strain"]
                    max_quals[snp["strain"]] = 0
                if snp["qual"] > max_quals[snp["strain"]]:
                    max_quals[snp["strain"]] = snp["qual"]
                if snp["qual"] > max_quals["All_strain"]:
                    max_quals["All_strain"] = snp["qual"]
                if (snp["depth"] >= depth_s) and (
                        snp["depth"] <= depth_b) and (
                        snp["dp4_sum"] >= cutoff_sum) and (
                        snp["dp4_frac"] >= cutoff_frac):
                    if snp["indel"] == -1:
                        apply_filter(args_snp.filters, snp, snps)
                    else:
                        if (snp["frac"] >= imf) and (
                                snp["indel"] >= idv):
                            apply_filter(args_snp.filters, snp, snps)
    fh.close()
    return max_quals, snps, dess, raw_snps


def check_overlap(new_snps, overlaps):
    count = 0
    element = 0
    first_overlap = True
    printeds = []
    count_overlap = len(overlaps)
    for key, value in new_snps.items():
        if first_overlap:
            for overlap in overlaps:
                if "print" in overlap.keys():
                    element += 1
                    printeds.append(overlap)
                else:
                    break
            first_overlap = False
        if ("print" not in overlaps[element].keys()):
            if len(printeds) == 0:
                value.append(overlaps[element])
                count += 1
            else:
                for printed in printeds:
                    if printed not in value:
                        value.append(overlaps[element])
                        count += 1
        if count_overlap != 0:
            if count == (len(new_snps.keys()) / count_overlap):
                overlaps[element]["print"] = True
                element += 1
                if element >= len(overlaps):
                    break
                count = 0


def overlap_position(qual_snps):
    '''deal with the conflict position of SNPs'''
    first = True
    qual_nooverlap_snps = {}
    num_overlap = 1
    qual_nooverlap_snps[num_overlap] = []
    conflicts = []
    for snp1 in qual_snps:
        overlaps = []
        overlaps.append(snp1)
        for snp2 in qual_snps:
            if (snp1 != snp2) and (snp1["strain"] == snp2["strain"]) and (
                    (snp2["pos"] - snp1["pos"] < len(snp1["ref"]))) and (
                    "print" not in snp2.keys()):
                overlaps.append(snp2)
        if len(overlaps) != 1:
            conflicts.append(overlaps)
            if first:
                for overlap in overlaps:
                    qual_nooverlap_snps[num_overlap] = []
                    qual_nooverlap_snps[num_overlap].append(overlap)
                    num_overlap += 1
                    overlap["print"] = True
                num_overlap = 1
                first = False
            else:
                new_snps = qual_nooverlap_snps.copy()
                index = len(qual_nooverlap_snps.keys())
                repeat = 0
                for overlap in overlaps:
                    if "print" in overlap.keys():
                        repeat += 1
                for times in range(1, len(overlaps) - repeat):
                    for key, value in qual_nooverlap_snps.items():
                        new_snps[key + index * times] = list(value)
                check_overlap(new_snps, overlaps)
                qual_nooverlap_snps = new_snps.copy()
        else:
            if "print" not in snp1.keys():
                if first:
                    qual_nooverlap_snps[num_overlap].append(snp1)
                    first = False
                else:
                    for key, value in qual_nooverlap_snps.items():
                        value.append(snp1)
                snp1["print"] = True
        if first:
            first = False
    return conflicts, qual_nooverlap_snps


def print_file(refs, out_ref, conflicts, key, values, mod_seq_init,
               mod_seqs, out_seq):
    num_seq = 1
    num_nt = 0
    paths = []
    if len(conflicts) == 0:
        paths.append("All")
    else:
        for conflict in conflicts:
            for path in conflict:
                for value in values:
                    if path == value:
                        paths.append(str(path["pos"]))
    if len(refs) != 0:
        for seq_name, ref_datas in refs.items():
            num_ref = 1
            for ref in ref_datas:
                out_ref.write("\t".join([str(key), "_".join(paths),
                              str(num_ref), ref, seq_name]) + "\n")
                num_ref += 1
    else:
        out_ref.write("\t".join([str(key), "_".join(paths), "1", "All",
                      mod_seq_init["genome"]]) + "\n")
    if len(mod_seqs) == 0:
        out_fasta = open("_".join([out_seq, mod_seq_init["genome"],
                                  str(key), "1.fa"]), "w")
        out_fasta.write(">{0}\n".format(mod_seq_init["genome"]))
        for nt in mod_seq_init["seq"]:
            num_nt += 1
            out_fasta.write("{0}".format(nt))
            if num_nt % 60 == 0:
                out_fasta.write("\n")
        out_fasta.close()
    else:
        for seq in mod_seqs:
            num_nt = 0
            out_fasta = open("_".join([out_seq, seq["genome"], str(key),
                                       str(num_seq)]) + ".fa", "w")
            out_fasta.write(">{0}\n".format(seq["genome"]))
            for nt in seq["seq"]:
                num_nt += 1
                out_fasta.write("{0}".format(nt))
                if num_nt % 60 == 0:
                    out_fasta.write("\n")
            num_seq += 1
            out_fasta.close()


def stat(max_quals, trans_snps, bam_number, stat_prefix,
         out_snp, args_snp, type_):
    out_stat = open("_".join([stat_prefix, type_]), "w")
    printed = False
    for strain, max_qual in max_quals.items():
        max_qual = int(((max_qual / 10) + 1) * 10)
        cutoffs = []
        if (strain == "All_strain") and (len(max_quals) > 2):
            printed = True
        elif (strain != "All_strain"):
            printed = True
        if printed:
            for cutoff in range(0, max_qual, 10):
                cutoffs.append(0)
            for snp in trans_snps:
                if (snp["strain"] == strain) or (strain == "All_strain"):
                    index = int(snp["qual"] / 10)
                    cutoffs[index] += 1
            num_cutoff = 10
            num_quality = 0
            out_stat.write(strain + ":\n")
            for cutoff in cutoffs:
                if args_snp.quality <= (num_cutoff - 10):
                    num_quality = num_quality + cutoff
                out_stat.write("the number of QUAL which is between "
                               "{0} and {1} = {2}\n".format(
                                   num_cutoff - 10, num_cutoff, cutoff))
                num_cutoff = num_cutoff + 10
            out_stat.write("the total numbers of QUAL which is "
                           "higher than {0} = {1}\n".format(
                               args_snp.quality, num_quality))
            plot_bar(cutoffs, strain, out_snp,
                     type_.replace(".csv", ".png"))
            printed = False
    out_stat.close()


def read_fasta(fasta_file):
    seqs = []
    first = True
    num_index = 0
    seq_name = ""
    with open(fasta_file, "r") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if first:
                    seqs.append({line[1:]: ""})
                    first = False
                else:
                    seqs.append({line[1:]: ""})
                    num_index += 1
                seq_name = line[1:]
            else:
                seqs[num_index][seq_name] = seqs[num_index][seq_name] + line
    return seqs


def gen_new_fasta(qual_nooverlap_snps, seqs, out_ref, conflicts, out_seq):
    refs = {}
    for key, values in qual_nooverlap_snps.items():
        for seq in seqs:
            for name in seq.keys():
                break
            refs[name] = []
            num_var = 0
            for strain, fasta in seq.items():
                mod_seq_init = {"genome": strain, "seq": fasta, "num_mod": 0}
                mod_seqs = []
            for snp in values:
                if snp["strain"] in seq.keys():
                    if "," in snp["alt"]:
                        num_var += 1
                        tmps = []
                        tmp_snps = []
                        alts = snp["alt"].split(",")
                        for alt in alts:
                            tmp_snp = snp.copy()
                            tmp_snp["alt"] = alt
                            tmp_snps.append(tmp_snp)
                            if len(mod_seqs) == 0:
                                num_mod_seqs = len(mod_seqs)
                                tmps.append(mod_seq_init.copy())
                            else:
                                num_mod_seqs = len(mod_seqs)
                                for mod_seq in mod_seqs:
                                    tmps.append(mod_seq.copy())
                        mod_seqs = list(tmps)
                        num_mod = 0
                        num = 1
                        refs[name] = gen_ref(tmp_snps, snp["pos"],
                                             refs[name], num_var)
                        for mod_seq in mod_seqs:
                            change(tmp_snps[num_mod], mod_seq)
                            if num >= num_mod_seqs:
                                num_mod += 1
                                num = 0
                            num += 1
                    else:
                        if len(mod_seqs) == 0:
                            change(snp, mod_seq_init)
                        else:
                            for mod_seq in mod_seqs:
                                change(snp, mod_seq)
        print_file(refs, out_ref, conflicts, key, values,
                   mod_seq_init, mod_seqs, out_seq)


def snp_detect(fasta_file, snp_file, depth_file, out_snp, out_seq,
               bam_number, stat_prefix, args_snp):
    max_quals, snps, dess, raw_snps = import_data(
            snp_file, args_snp, bam_number, depth_file)
    out_best = open(out_snp + "_best.vcf", "w")
    out_ref = open(out_snp + "_seq_reference.csv", "w")
    out_best.write("\n".join(dess) + "\n")
    out_best.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL"
                    "\tFILTER\tINFO\tFORMAT\tBAM\n")
    best_snps = []
    for snp in snps:
        if snp["qual"] >= args_snp.quality:
            out_best.write(snp["all_info"] + "\n")
            best_snps.append(snp)
    conflicts, qual_nooverlap_snps = overlap_position(best_snps)
    stat(max_quals, raw_snps, bam_number, stat_prefix, out_snp, args_snp, "raw.csv")
    stat(max_quals, snps, bam_number, stat_prefix, out_snp, args_snp, "best.csv")
    seqs = read_fasta(fasta_file)
    gen_new_fasta(qual_nooverlap_snps, seqs, out_ref, conflicts, out_seq)
    out_best.close()
    out_ref.close()
