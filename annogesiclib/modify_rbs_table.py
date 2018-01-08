import shutil
import csv


def import_data(row):
    return{"strain": row[1], "strand": row[2],
           "associate": row[3], "start_seq": int(row[4]),
           "end_seq": int(row[5]), "rfam": row[6], "e": row[7],
           "score": row[8],
           "start_align": int(row[9]), "end_align": int(row[10]),
           "info": row[0:6], "ID": row[0]}


def modify_table(table, output_all):
    first = True
    rbss = []
    out = open("tmp.csv", "w")
    out.write("#ID\tGenome\tStrand\tAssociated_CDS\tStart_genome\t"
              "End_genome\tRfam\tE_value\tScore\tStart_align\tEnd_align\n")
    if output_all:
        with open(table) as fh:
            for line in fh:
                line = line.strip()
                if first:
                    first = False
                    rbss.append(line)
                    out.write(line + "\n")
                else:
                    if line not in rbss:
                        rbss.append(line)
                        out.write(line + "\n")
    else:
        fh = open(table, "r")
        for row in csv.reader(fh, delimiter='\t'):
            rbss.append(import_data(row))
        for rbs1 in rbss:
            repeat = False
            if "print" not in rbs1.keys():
                rbs1["print"] = True
                for rbs2 in rbss:
                    if (rbs1["strain"] == rbs2["strain"]) and \
                       (rbs1["strand"] == rbs2["strand"]) and \
                       (rbs1["ID"] == rbs2["ID"]):
                        if "print" not in rbs2.keys():
                            rbs2["print"] = True
                            repeat = True
                if not repeat:
                    out.write("\t".join(rbs1["info"] + [rbs1["rfam"],
                                        rbs1["e"], rbs1["score"],
                                        str(rbs1["start_align"]),
                                        str(rbs1["end_align"])]) + "\n")
        fh.close()
    out.close()
    shutil.move("tmp.csv", table)
