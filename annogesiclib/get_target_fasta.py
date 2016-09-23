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
        self.folders = {"tmp_tar": os.path.join(tar_folder, "tmp"),
                        "tmp_ref": os.path.join(ref_folder, "tmp")}

    def get_target_fasta(self, mut_table, tar_folder, ref_folder, output):
        self.multiparser.parser_fasta(ref_folder)
        if "tmp" in os.listdir(tar_folder):
            shutil.rmtree(self.folders["tmp_tar"])
        os.mkdir(self.folders["tmp_tar"])
        self.seq_editer.modify_seq(self.folders["tmp_ref"], mut_table,
                                   self.folders["tmp_tar"])
        print("transfer to target fasta...")
        if output is not None:
            for file_ in output:
                first = True
                datas = file_.split(":")
                filename = datas[0]
                strains = datas[1].split("_and_")
                out = open(os.path.join(tar_folder, filename + ".fa"), "w")
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
                        print("Error:no fasta information of {0}.fa".format(
                              strain))
                out.close()
        else:
            self.helper.move_all_content(self.folders["tmp_tar"],
                                         tar_folder, [".fa"])
        shutil.rmtree(self.folders["tmp_tar"])
        shutil.rmtree(self.folders["tmp_ref"])
        self.helper.remove_all_content(ref_folder, "_folder", "dir")
        print("please use the new fasta file to remapping again.")
        print("Then copy BAMs and wigs back to input/align_results/BAMs "
              "and input/align_results/wigs")
