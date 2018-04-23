import os
import shutil
import sys
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from annogesiclib.seqmodifier import SeqModifier


class SeqEditer(object):
    '''Edit the sequence if it is needed'''

    def _row_to_location(self, row, out_name):
        return({"ref_id": row[0], "target_id": "_".join([out_name, row[0]]),
                "datas": [{"ref_nt": row[3],
                           "tar_nt": row[4], "position": row[1]}]})

    def _import_data(self, mod_table_file, out_name):
        datas = []
        first = True
        num_index = 0
        fh = open(mod_table_file)
        for row in csv.reader(fh, delimiter="\t"):
            if row[0].startswith("#"):
                continue
            else:
                if first:
                    datas.append(self._row_to_location(row, out_name))
                    pre_ref_id = row[0].strip()
                    first = False
                else:
                    if (row[0] == pre_ref_id):
                        datas[num_index]["datas"].append(
                              {"ref_nt": row[3].strip(),
                               "tar_nt": row[4].strip(),
                               "position": row[1].strip()})
                    else:
                        datas.append(self._row_to_location(row, out_name))
                        num_index += 1
                        pre_ref_id = row[0].strip()
        fh.close()
        return datas

    def modify_seq(self, fasta_folder, mod_table_file, output_folder, out_name):
        datas = self._import_data(mod_table_file, out_name)
        for data in datas:
            seq = ""
            if (data["ref_id"] + ".fa") in os.listdir(fasta_folder):
                filename = os.path.join(fasta_folder, data["ref_id"] + ".fa")
                with open(filename, "r") as fasta:
                    for line in fasta:
                        line = line.strip()
                        if len(line) != 0:
                            if line[0] != ">":
                                seq = seq + line
                seq_modifier = SeqModifier(seq)
                for change in data["datas"]:
                    if change["ref_nt"] == "-":
                        seq_modifier.insert(
                                     int(change["position"]), change["tar_nt"])
                    elif change["tar_nt"] == "-":
                        seq_modifier.remove(int(change["position"]),
                                            len(change["ref_nt"]))
                    else:
                        seq_modifier.replace(
                                     int(change["position"]), change["tar_nt"])
                record = SeqRecord(Seq(seq_modifier.seq()))
                record.id = data["target_id"]
                record.description = ""
                SeqIO.write(record, os.path.join(
                            output_folder, record.id + ".fa"), "fasta")

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
                    elif (len(mod) != 5) and (line[0] == ">"):
                        new_header = line.split(" ")[0]
                    elif (line[0] != ">"):
                        print("Error: No proper header!!")
                        sys.exit()
                    line = new_header
                output_fh.write(line + "\n")
        output_fh.close()
        shutil.move(tmp_file_path, input_file)
