import os
import shutil
import csv
from annogesiclib.multiparser import Multiparser
from annogesiclib.seq_editer import SeqEditer
from annogesiclib.helper import Helper


class TargetFasta(object):
    '''detection of sRNA target interaction'''

    def __init__(self, tar_folder, ref_folder):
        self.multiparser = Multiparser()
        self.seq_editer = SeqEditer()
        self.helper = Helper()
        self.folders = {"tmp_tar": os.path.join(tar_folder, "tmp")}

    def gen_folder(self, out_folder, ref_files):
        new_ref_folder = os.path.join(out_folder, "tmp_reference")
        self.helper.check_make_folder(new_ref_folder)
        for file_ in ref_files:
            shutil.copy(file_, new_ref_folder)
        self.folders["tmp_ref"] = os.path.join(new_ref_folder, "tmp")
        self.multiparser.parser_fasta(new_ref_folder)
        if os.path.exists(os.path.join(out_folder, "fasta_files")):
            shutil.rmtree(os.path.join(out_folder, "fasta_files"))
            os.mkdir(os.path.join(out_folder, "fasta_files"))
        if os.path.exists(self.folders["tmp_tar"]):
            shutil.rmtree(self.folders["tmp_tar"])
        os.mkdir(self.folders["tmp_tar"])
        return new_ref_folder

    def get_target_fasta(self, mut_table, tar_folder, ref_files,
                         out_name, out_folder, log):
        new_ref_folder = self.gen_folder(out_folder, ref_files)
        log.write("Running seq_editor.py for updating sequence.\n")
        self.seq_editer.modify_seq(self.folders["tmp_ref"], mut_table,
                                   self.folders["tmp_tar"], out_name)
        print("Updating the reference sequences")
        mh = open(mut_table, "r")
        pre_strain = None
        out = None
        strain_num = 0
        for row in csv.reader(mh, delimiter='\t'):
            if not row[0].startswith("#"):
                if (pre_strain != row[0]):
                    strain_num = strain_num + 1
                    tmp_tar_name = "_".join([out_name, row[0]]) + ".fa"
                    fasta = os.path.join(out_folder, "fasta_files",
                                         tmp_tar_name)
                    if out is not None:
                        out.close()
                    out = open(fasta, "w")
                    if tmp_tar_name in os.listdir(self.folders["tmp_tar"]):
                        with open(os.path.join(
                                  self.folders["tmp_tar"],
                                  tmp_tar_name)) as f_h:
                            for line in f_h:
                                out.write(line)
                    else:
                        print("Error: No updated information of {0}.fa".format(
                              row[0]))
                pre_strain = row[0]
        out.close()
        out_seq = out_name + ".fa"
        if os.path.exists(out_seq):
            os.remove(out_seq)
        if strain_num == 1:
            o_s = open(out_seq, "w")
            for seq in os.listdir(os.path.join(out_folder, "fasta_files")):
                if seq.endswith(".fa"):
                    with open(os.path.join(
                            out_folder, "fasta_files", seq)) as t_h:
                        for line in t_h:
                            if len(line) != 0:
                                if line.startswith(">"):
                                    o_s.write(">" + out_name + "\n")
                                else:
                                     o_s.write(line)
                    os.remove(os.path.join(out_folder, "fasta_files", seq))
            o_s.close()
        else:
            for seq in os.listdir(os.path.join(out_folder, "fasta_files")):
                if seq.endswith(".fa"):
                    os.system(" ".join(["cat", os.path.join(
                                            out_folder, "fasta_files", seq),
                                        ">>", out_seq]))
                    os.remove(os.path.join(out_folder, "fasta_files", seq))
        shutil.move(out_seq, os.path.join(
            out_folder, "fasta_files", out_seq))
        shutil.rmtree(self.folders["tmp_tar"])
        shutil.rmtree(self.folders["tmp_ref"])
        if "tmp_reference" in os.listdir(out_folder):
            shutil.rmtree(new_ref_folder)
        log.write("\t" + os.path.join(out_folder, "fasta_files", out_seq) + 
                  " is generated.\n")
        print("Please use the new fasta files to remapping again.")
