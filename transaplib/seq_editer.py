#!/usr/bin/python

import os
import sys
import csv
from transaplib.seqmodifier import SeqModifier
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


class Seq_Editer(object):

    def _row_to_location(self, row):
        return({"target_id": row[0], "ref_id": row[1], "datas": [{"ref_nt": row[2],
                "tar_nt": row[4], "position": row[3]}]})

    def _import_data(self, mod_table_file, datas):
        first = True
        num_index = 0
        for row in csv.reader(open(mod_table_file), delimiter="\t"):
            if row[0].startswith("target_id"):
                continue
            else:
                if first:
                    datas.append(self._row_to_location(row))
                    pre_ref_id = row[1].strip()
                    pre_tar_id = row[0].strip()
                    first = False
                else:
                    if (row[1] == pre_ref_id) and \
                       (row[0] == pre_tar_id):
                        datas[num_index]["datas"].append({"ref_nt": row[2].strip(),
                               "tar_nt": row[4].strip(), "position": row[3].strip()})
                    else:
                        datas.append(self._row_to_location(row))
                        num_index += 1
                        pre_ref_id = row[1].strip()
                        pre_tar_id = row[0].strip()
    
    def modify_seq(self, fasta_folder, mod_table_file, output_folder):
        datas = []
        self._import_data(mod_table_file, datas)
        for data in datas:
            seq = ""
            if (data["ref_id"] + ".fa") in os.listdir(fasta_folder):
                filename = fasta_folder + "/" + data["ref_id"] + ".fa"
                with open (filename, "r") as fasta:
                    for line in fasta:
                        line = line.strip()
                        if line[0] != ">":
                            seq = seq + line
                seq_modifier = SeqModifier(seq)
                for change in data["datas"]:
                    if change["ref_nt"] == "-":
                        seq_modifier.insert(int(change["position"]), change["tar_nt"])
                    elif change["tar_nt"] == "-":
                        seq_modifier.remove(int(change["position"]))
                    else:
                        seq_modifier.replace(int(change["position"]), change["tar_nt"])
                record = SeqRecord(Seq(seq_modifier.seq()))
                record.id = data["target_id"]
                record.description = ""
                SeqIO.write(record, output_folder + "/" + record.id + ".fa", "fasta")

    def modify_header(self, input_file):
        first = True
        tmp_file_path = input_file + "_TMP"
        output_fh = open(input_file + "_TMP", "w")
        with open(input_file, "r") as s_h:
            for line in s_h:
                line = line.strip()
                if first:
                    first = False
                    if (line[0] != ">"):
                        print("Error: No proper header!!")
                        sys.exit()
                if line.startswith(">"):
                    mod = line.split("|")
                    folder = input_file.split("/")
                    folder = "/".join(folder[:-1])
                    if (len(mod) == 5) and (line[0] == ">"):
                        new_header = ">%s" % (mod[3])
                        print(folder + "/" + mod[3] + ".fa")
                    elif (len(mod) != 5) and (line[0] == ">"):
                        new_header = line
                    elif (line[0] != ">"):
                        print("Error: No proper header!!")
                        sys.exit()
                    line = new_header
                output_fh.write(line + "\n")
        output_fh.close()
        os.rename(tmp_file_path, input_file)
