import os
import shutil
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
        if "tmp_tar" in os.listdir(out_folder):
            shutil.rmtree(self.folders["tmp_tar"])
        os.mkdir(self.folders["tmp_tar"])
        return new_ref_folder

    def get_target_fasta(self, mut_table, tar_folder, ref_files,
                         output, out_folder):
        new_ref_folder = self.gen_folder(out_folder, ref_files)
        self.seq_editer.modify_seq(self.folders["tmp_ref"], mut_table,
                                   self.folders["tmp_tar"])
        print("Transfering to target fasta")
        for file_ in output:
            first = True
            datas = file_.split(":")
            filename = datas[0]
            strains = datas[1].split(",")
            out = open(filename, "w")
            for strain in strains:
                if strain + ".fa" in os.listdir(self.folders["tmp_tar"]):
                    if first:
                        first = False
                    else:
                        out.write("\n")
                    with open(os.path.join(
                              self.folders["tmp_tar"],
                              strain + ".fa")) as f_h:
                        for line in f_h:
                            out.write(line)
                else:
                    print("Error: No fasta information of {0}.fa".format(
                          strain))
            out.close()
        shutil.rmtree(self.folders["tmp_tar"])
        shutil.rmtree(self.folders["tmp_ref"])
        if "tmp_reference" in os.listdir(out_folder):
            shutil.rmtree(new_ref_folder)
        print("Please use the new fasta file to remapping again.")
