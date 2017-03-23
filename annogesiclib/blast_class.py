import csv


def read_file(srna_file, nums):
    srna_f = open(srna_file, "r")
    for row in csv.reader(srna_f, delimiter="\t"):
        if (row[-6] != "NA") and (row[0] != "Rank"):
            if row[1] not in nums.keys():
                nums[row[1]] = {}
            if row[2] not in nums[row[1]].keys():
                nums[row[1]][row[2]] = 1
            else:
                nums[row[1]][row[2]] += 1
            if row[2] not in nums["total"].keys():
                nums["total"][row[2]] = 1
            else:
                nums["total"][row[2]] += 1
    srna_f.close()


def blast_class(srna_file, out_file):
    '''statistics of the results of blast sRNA database'''
    nums = {}
    nums["total"] = {}
    read_file(srna_file, nums)
    out = open(out_file, "w")
    if len(nums) > 1:
        if len(nums) > 2:
            out.write("All genomes:\n")
            out.write("sRNA_name\tamount\n")
            for blast, num in nums["total"].items():
                out.write("{0}\t{1}\n".format(blast, num))
        for strain, srna_name in nums.items():
            if strain != "total":
                out.write(strain + ":\n")
                out.write("sRNA_name\tamount\n")
                for blast, num in srna_name.items():
                    out.write("{0}\t{1}\n".format(blast, num))
    else:
        out.write("No known sRNA!!\n")
    out.close()
