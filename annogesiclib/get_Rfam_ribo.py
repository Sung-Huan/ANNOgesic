import csv


def rbs_from_rfam(ribo_table, rfam_file, out_file):
    ribos = []
    out = open(out_file, "w")
    f_h = open(ribo_table, "r")
    for row in csv.reader(f_h, delimiter="\t"):
        if not row[0].startswith("#"):
            ribos.append(row[0].strip())
    detect = False
    with open(rfam_file, "r") as r_h:
        for line in r_h:
            line = line.rstrip("\n")
            datas = line.split(" ")
            if ("INFERNAL" in datas[0]) or (
                    "HMMER" in datas[0]):
                header = line
                detect = False
            elif "NAME" in datas[0]:
                name = line
            elif ("ACC" in datas[0]):
                for ribo in ribos:
                    if datas[-1] == ribo:
                        out.write("{0}\n{1}\n{2}\n".format(header, name, line))
                        detect = True
            else:
                if (detect):
                    out.write(line + "\n")
    out.close()
    f_h.close()
