import os
import csv
import shutil


def mapping_ribos(table_folder, id_file, feature):
    ids = []
    ih = open(id_file, "r")
    for row in csv.reader(ih, delimiter='\t'):
        if not row[0].startswith("#"):
            ids.append({"id": row[0].strip(),
                        "name": row[1].strip(),
                        "info": row[2].strip()})
    for table_file in os.listdir(table_folder):
        if table_file.endswith("_" + feature + ".csv"):
            tmp_table = os.path.join(table_folder, "tmp" + table_file)
            table_file = os.path.join(table_folder, table_file)
            out = open(tmp_table, "w")
            tables = []
            fh = open(table_file, "r")
            out.write("#ID\tGenome\tStrand\tAssociated_CDS\tStart_genome\t"
                      "End_genome\tRfam\tE_value\tStart_align\tEnd_align\n")
            for row in csv.reader(fh, delimiter='\t'):
                if not row[0].startswith("#"):
                    tables.append({"input": row[0:6], "Rfam": row[6],
                                   "e": row[7], "start": row[8],
                                   "end": row[9]})
            for table in tables:
                for id_ in ids:
                    if table["Rfam"] == id_["id"]:
                        name = id_["name"]
                out.write("\t".join(table["input"] + [table["Rfam"], name,
                                    table["e"], table["start"],
                                    table["end"]]) + "\n")
            out.close()
            os.remove(table_file)
            shutil.move(tmp_table, table_file)
