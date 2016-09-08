import csv
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')


def plot(subs, total, unknown, strain, prefix_name):
    nums = []
    nums_no_unknown = []
    classes = []
    classes_no_unknown = []
    width = 0.4
    tmp_unknown = ["Unknown", 0]
    sort_subs = sorted(subs.items(),
                       key=lambda x: (x[1]), reverse=True)
    for datas in sort_subs:
        if datas[0] == "Unknown":
            tmp_unknown = datas
        else:
            nums.append(datas[1])
            nums_no_unknown.append(datas[1])
            classes.append(datas[0])
            classes_no_unknown.append(datas[0])
    nums.append(tmp_unknown[1])
    classes.append(tmp_unknown[0])
    plt.figure(figsize=(12, 16))
    plt.subplot(211)
    ind = np.arange(len(nums))
    plt.bar(ind, nums, width, color='#FF9999')
    plt.title('Subcellular localization including Unknown\n', fontsize=24)
    plt.ylabel('Amount', fontsize=20)
    plt.yticks(fontsize=16)
    plt.xlim([0, len(nums) + 1])
    plt.xticks(ind+width, classes, rotation=40, fontsize=20, ha='right')
    plt.tight_layout(2, None, None, None)
    plt.subplot(212)
    ind = np.arange(len(nums_no_unknown))
    plt.bar(ind, nums_no_unknown, width, color='#FF9999')
    plt.title('Subcellular localization excluding Unknown\n', fontsize=24)
    plt.ylabel('Amount', fontsize=20)
    plt.xlim([0, len(nums_no_unknown) + 1])
    plt.xticks(ind+width, classes_no_unknown, rotation=40,
               fontsize=20, ha='right')
    plt.yticks(fontsize=16)
    plt.tight_layout(2, None, None, None)
    plt.savefig("_".join([prefix_name, strain, "sublocal.png"]))


def read_table(psortb_file):
    subs = {}
    subs["all_strain"] = {}
    total_nums = {}
    total_nums["all_strain"] = 0
    unknown_nums = {}
    unknown_nums["all_strain"] = 0
    pre_strain = ""
    f_h = open(psortb_file, "r")
    for row in csv.reader(f_h, delimiter="\t"):
        if not row[0].startswith("#"):
            if pre_strain != row[0]:
                subs[row[0]] = {}
                pre_strain = row[0]
                total_nums[row[0]] = 0
                unknown_nums[row[0]] = 0
            if row[5] not in subs[row[0]].keys():
                if row[5] == "Unknown":
                    unknown_nums[row[0]] += 1
                subs[row[0]][row[5]] = 1
                total_nums[row[0]] += 1
            else:
                if row[5] == "Unknown":
                    unknown_nums[row[0]] += 1
                subs[row[0]][row[5]] += 1
                total_nums[row[0]] += 1
            if row[5] not in subs["all_strain"].keys():
                if row[5] == "Unknown":
                    unknown_nums["all_strain"] += 1
                subs["all_strain"][row[5]] = 1
                total_nums["all_strain"] += 1
            else:
                if row[5] == "Unknown":
                    unknown_nums["all_strain"] += 1
                subs["all_strain"][row[5]] += 1
                total_nums["all_strain"] += 1
    f_h.close()
    return subs, total_nums, unknown_nums


def print_file_and_plot(sub, total_nums, unknown_nums,
                        strain, out_stat, prefix_name):
    plot(sub, total_nums[strain], unknown_nums[strain], strain, prefix_name)
    out_stat.write(strain + ":\n")
    out_stat.write("Total including Unknown is {0}; "
                   "Total excluding Unknown is {1}\n".format(
                       total_nums[strain],
                       total_nums[strain] - unknown_nums[strain]))
    for local, num in sub.items():
        if local != "Unknown":
            out_stat.write(
                "\t{0}\t{1}(including Unknown {2}; "
                "excluding Unknonwn {3})\n".format(
                    local, num, float(num) / float(total_nums[strain]),
                    float(num) / (float(total_nums[strain]) - float(
                        unknown_nums[strain]))))
        else:
            out_stat.write("\t{0}\t{1}(including Unknown {2})\n".format(
                local, num, float(num) / float(total_nums[strain])))


def stat_sublocal(psortb_file, prefix_name, stat_file):
    subs, total_nums, unknown_nums = read_table(psortb_file)
    out_stat = open(stat_file, "w")
    if len(subs) > 2:
        print_file_and_plot(subs["all_strain"], total_nums, unknown_nums,
                            "all_strain", out_stat, prefix_name)
    for strain, sub in subs.items():
        if strain != "all_strain":
            print_file_and_plot(sub, total_nums, unknown_nums, strain,
                                out_stat, prefix_name)
