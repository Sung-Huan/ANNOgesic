#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
import time
from subprocess import call, Popen
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.potential_target import potential_target
from annogesiclib.format_fixer import Fixer
from annogesiclib.merge_rnaplex_rnaup import merge_srna_target
from annogesiclib.gff3 import Gff3Parser


class sRNA_target_prediction(object):

    def __init__(self):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.fixer = Fixer()
        self.gff_parser = Gff3Parser()
        self.target_seq_path = os.path.join(out_folder, "target_seqs")
        self.srna_seq_path = os.path.join(out_folder, "sRNA_seqs")
        self.rnaplex_path = os.path.join(out_folder, "RNAplex")
        self.rnaup_path = os.path.join(out_folder, "RNAup")
        self.merge_path = os.path.join(out_folder, "merge")
        self.srna_path = os.path.join(srnas, "tmp")
        self.fasta_path = os.path.join(fastas, "tmp")
        self.gff_path = os.path.join(gffs, "tmp")
        self.tmp = "tmp"
        self.tmp_rnaup = "tmp_rnaup"
        self.tmp_log = "tmp_log"
        self.tmp_all_fa = "tmp*.fa"
        self.tmp_all_txt = "tmp*.txt"

    def _check_gff(self, gffs):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(gffs, gff))

    def _run_rnaplfold(self, vienna_path, file_type, win_size, span,
                       unstr_region, seq_path, prefix, out_path):
        current = os.getcwd()
        os.chdir(out_path)
        command = " ".join([vienna_path + "/Progs/RNAplfold",
                            "-W", str(win_size),
                            "-L", str(span),
                            "-u", str(unstr_region),
                            "-O"])
        if file_type == "sRNA":
            os.system("<".join([command, os.path.join(current, seq_path, 
                                "_".join([self.tmp ,prefix, file_type + ".fa"]))]))
        else:
            os.system("<".join([command, os.path.join(current, seq_path, 
                                "_".join([prefix, file_type + ".fa"]))]))
        os.chdir(current)
    def _wait_process(self, processes):
        for p in processes:
            p.wait()
        if p.stdout:
            p.stdout.close()
        if p.stdin:
            p.stdin.close()
        if p.stderr:
            p.stderr.close()
        try:
            p.kill()
        except OSError:
            pass
        time.sleep(5)

    def _organize_files(self, srna_prefix, target_path, file_type, folder):
        tmp_path = folder.replace(".gff_folder", "") + file_type
        self.helper.check_make_folder(tmp_path)
        for srna in os.listdir(target_path):
            if srna == srna_prefix:
                if file_type != "_merge":
                    for txt in os.listdir(os.path.join(target_path, srna)):
                        if txt.endswith("txt"):
                            os.copy(os.path.join(target_path, srna, txt), tmp_path)
                for cvs_file in os.listdir(os.path.join(target_path, srna)):
                os.copy(os.path.join(target_path, srna, cvs_file), tmp_path)

    def _sort_sRNA_fasta(self, fasta, prefix, path):
        out = open(os.path.join(path, "_".join([self.tmp, prefix, "sRNA.fa"])), "w")
        srnas = []
        with open(fasta) as f_h:
            for line in f_h:
                line = line.strip()
                if line.startswith(">"):
                    name = line[1:]
                else:
                    srnas.append({"name": name, "seq": line, "len": len(line)})
        srnas = sorted(srnas, key = lambda x: (x["len"]))
        for srna in srnas:
            out.write(">" + srna["name"].split("|")[0] + "\n")
            out.write(srna["seq"] + "\n")

    def _read_fasta(self, fasta_file):
        seq = ""
        with open(fasta_file, "r") as seq_f:
            for line in seq_f:
                line = line.strip()
                if line.startswith(">"):
                    continue
                else:
                    seq = seq + line
        return seq

    def _get_specific_seq(srna_file, seq_file, srna_out, querys):
        for query in querys:
            srna_datas = query.split(":")
            srna = {"seq_id": srna_datas[0], "strand": srna_datas[1],
                    "start": int(srna_datas[2]), "end": int(srna_datas[3])}
            gff_f = open(srna_file, "r")
            out = open(srna_out, "w")
            seq = self._read_fasta(seq_file)
            num = 0
            for entry in self.gff_parser.entries(gff_f):
                if (entry.seq_id == srna["seq_id"]) and \
                   (entry.strand == srna["strand"]) and \
                   (entry.start == srna["start"]) and \
                   (entry.end == srna["end"]):
                    if "ID" in entry.attributes.keys():
                        id_ = entry.attributes["ID"]
                    else:
                        id_ = entry.feature + str(num)
                    gene = self.heler.extract_gene(seq, entry.start, entry.end, entry.strand)
                    out.write(">{0}|{1}|{2}|{3}|{4}\n{5}\n".format(id_, entry.seq_id, entry.start,
                                                                   entry.end, entry.strand, gene))
                    num += 1
            gff_f.close()

    def _gen_seq(self, prefixs, querys):
        print("Generating sRNA fasta files...")
        for srna in os.listdir(self.srna_path): # extract sRNA seq
            if srna.endswith("_sRNA.gff"):
                prefix = srna.replace("_sRNA.gff", "")
                prefixs.append(prefix)
                srna_out = os.path.join(self.srna_seq_path, "_".join([prefix, "sRNA.fa"]))
                if "all" in querys:
                    self.helper.get_seq(os.path.join(self.srna_path, srna),
                                        os.path.join(self.fasta_path, prefix + ".fa"),
                                        srna_out)
                else:
                    self._get_specific_seq(os.path.join(self.srna_path, srna),
                                           os.path.join(self.fasta_path, prefix + ".fa"),
                                           srna_out, querys)
                self._sort_sRNA_fasta(srna_out, prefix, self.srna_seq_path)
        print("Generating target fasta files...")
        for gff in os.listdir(self.gff_path): # extract target seq
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                self.potential_target(os.path.join(self.gff_path, gff),
                                      os.path.join(self.fasta_path, prefix + ".fa"),
                                      os.path.join(self.target_seq_path))
                file_num = 1
                num = 0
                sub_prefix = os.path.join(self.target_seq_path,     
                                  "_".join([prefix, "target"])) 
                sub_out = open("".join([sub_prefix, str(file_num), ".fa"]), "w")
                with open(sub_prefix + ".fa"])), "r") as t_f:
                    for line in t_f:
                        line = line.strip()
                        if line.startswith(">"):
                            num += 1
                        if (num == 100):
                            num = 0
                            file_num += 1
                            sub_out.close()
                            sub_out = open("".join([sub_prefix, str(file_num), ".fa"]), "w")
                        sub_out.write(line + "\n")

    def _rna_plex(self, prefixs, rnaplex_path, vienna_path, win_size_s, span_s, 
                  unstr_region_rnaplex_s, srna_seq_path, win_size_t, span_t, 
                  unstr_region_rnaplex_t, target_seq_path, inter_length, energy, 
                  duplex_dist, core_plex):
        for prefix in prefixs:
            print("Running RNAplfold of {0}".format(prefix)) # Run RNAplfold first for RNAplex
            self.helper.check_make_folder(os.path.join(self.rnaplex_path, prefix))
            rnaplfold_path = os.path.join(rnaplex_path, prefix, "RNAplfold")
            os.mkdir(rnaplfold_path)
            self._run_rnaplfold(vienna_path, "sRNA", win_size_s, span_s,
                                unstr_region_rnaplex_s, srna_seq_path,
                                prefix, rnaplfold_path)
            self._run_rnaplfold(vienna_path, "target", win_size_t, span_t,
                                unstr_region_rnaplex_t, target_seq_path,
                                prefix, rnaplfold_path)
            print("Running RNAplex of {0}".format(prefix))
            num_process = 0
            processes = []
            for seq in os.listdir(target_seq_path): # Run RNAplex
                if (prefix in seq) and ("_target_" in seq):
                    print(seq)
                    out_rnaplex = open(os.path.join(rnaplex_path, prefix, 
                                  "_".join([prefix, "RNAplex", str(num_process) + ".txt"])), "w")
                    num_process += 1
                    p = Popen([os.path.join(vienna_path, "Progs", "RNAplex"),
                               "-q", os.path.join(srna_seq_path, "_".join([self.tmp, prefix, "sRNA.fa"])),
                               "-t", os.path.join(target_seq_path, seq),
                               "-l", str(inter_length),
                               "-e", str(energy), "-z", str(duplex_dist),
                               "-a", rnaplfold_path], stdout=out_rnaplex)
                    processes.append(p)
                    if num_process % core_plex == 0:
                        self._wait_process(processes)
            self._wait_process(processes)
            rnaplex_file = os.path.join(rnaplex_path, prefix, "_".join([prefix, "RNAplex.txt"]))
            if "_".join([prefix, "RNAplex.txt"]) in os.listdir(os.path.join(rnaplex_path, prefix)):
                os.remove(rnaplex_file)
            for index in range(0, num_process):
                self.helper.merge_file(os.path.join(rnaplex_path, prefix, 
                                       "_".join([prefix, "RNAplex", str(index) + ".txt"])),
                                       rnaplex_file)
            self.helper.remove_all_content(os.path.join(rnaplex_path, prefix), "_RNAplex_", "file")
            self.fixer.fix_rnaplex(rnaplex_file, self.tmp)
            os.rename(self.tmp, rnaplex_file)

    def _run_rnaup(self, num_up, index, vienna_path, unstr_region_rnaup, processes):
        for index in range(1, num_up + 1):
            out_tmp_up = open("".join([self.tmp_rnaup, str(index), ".txt"]), "w")
            out_err = open("".join([self.tmp_log, str(index), ".txt"]), "w")
            in_up = open("".join([self.tmp, str(index), ".fa"]), "r")
            p = Popen([os.path.join(vienna_path, "Progs", "RNAup"),
                       "-u", str(unstr_region_rnaup),
                       "-o", "--interaction_first"],
                       stdin = in_up, stdout=out_tmp_up, stderr=out_err)
            processes.append(p)

    def _merge_txt(self, num_up, index, out_ranup, our_log):
        for index in range(1, num_up + 1):
            self.helper.merge_file("".join([self.tmp_rnaup, str(index), ".txt"]), out_rnaup)
            self.helper.merge_file("".join([self.tmp_log, str(index), ".txt"]), out_log)

    def _rnaup(self, prefixs, rnaup_path, srna_seq_path, target_seq_path, core_up, vienna_path,
               unstr_region_rnaup):
        for prefix in prefixs:
            print("Running RNAup of {0}".format(prefix))
            self.helper.check_make_folder(os.path.join(rnaup_path, prefix))
            num_up = 0
            processes = []
            out_rnaup = os.path.join(rnaup_path, prefix, "_".join([prefix + "_RNAup.txt"]))
            out_log = os.path.join(rnaup_path, prefix, "_".join([prefix + "_RNAup.log"]))
            if "_".join([prefix, "RNAup.txt"]) in os.listdir(os.path.join(rnaup_path, prefix)):
                os.remove(out_rnaup)
                os.remove(out_log)
            with open(os.path.join(srna_seq_path, "_".join([self.tmp, prefix, "sRNA.fa"])), "r") as s_f:
                for line in s_f:
                    line = line.strip()
                    if line.startswith(">"):
                        print(line)
                        num_up += 1
                        out_up = open("".join([self.tmp, str(num_up), ".fa"]), "w")
                        out_up.write(line + "\n")
                    else:
                        out_up.write(line + "\n")
                        out_up.close()
                        self.helper.merge_file(os.path.join(target_seq_path, "_".join([prefix, "target.fa"])),
                                               "".join([self.tmp, str(num_up), ".fa"]))
                        if num_up == core_up:
                            self._run_rnaup(num_up, index, vienna_path, unstr_region_rnaup, processes)
                            os.system("rm " + self.tmp_all_fa)
                            self._merge_txt(num_up, index, out_ranup, our_log)
                            os.system("rm " + self.tmp_all_txt)
                            processes = []
                            num_up = 0
            self._run_rnaup(num_up, index, vienna_path, unstr_region_rnaup, processes)
            if len(processes) != 0:
                self._wait_process(processes)
                os.system("rm " + self.tmp_all_fa)
                self._merge_txt(num_up, index, out_ranup, our_log)
                os.system("rm " + self.tmp_all_txt)

    def _merge_rnaplex_rnaup(self, prefixs, merge_path, program, rnaplex_path, 
                             rnaup_path, srna_path, gff_path):
        for prefix in prefixs: # merge and rank the sRNA target output
            self.helper.check_make_folder(os.path.join(merge_path, prefix))
            print("Ranking {0} now...".format(prefix))
            if (program == "both") or (program == "RNAplex"):
                rnaplex_file = os.path.join(rnaplex_path, prefix, "_".join([prefix, "RNAplex.txt"]))
                out_rnaplex = os.path.join(rnaplex_path, prefix, "_".join([prefix, "RNAplex_rank.csv"]))
            elif (program == "both") or (program == "RNAup"):
                rnaup_file = os.path.join(rnaup_path, prefix, "_".join([prefix, "RNAup.txt"]))
                out_rnaup = os.path.join(rnaup_path, prefix, "_".join([prefix, "RNAup_rank.csv"]))
            merge_srna_target(rnaplex_file, rnaup_file, top, out_rnaplex, out_rnaup,
                              os.path.join(merge_path, prefix, "_".join([prefix, "merge.csv"])),
                              os.path.join(merge_path, prefix, "_".join([prefix, "overlap.csv"])), 
                              os.path.join(srna_path, "_".join([prefix, "sRNA.gff"])), 
                              os.path.join(gff_path, prefix + ".gff"))

    def _organize_gff(self, gffs, rnaplex_path, rnaup_path, merge_path):
        srna_prefixs = []
        folders = []
        for folder in os.listdir(gffs):
            if folder.endswith(".gff_folder"):
                folders.append(folder.replace(".gff_folder", ""))
                for gff in os.listdir(os.path.join(gffs, folder)):
                    if gff.endswith(".gff"):
                       srna_prefixs.append(gff.replace(".gff", ""))
                for srna_prefix in srna_prefixs:
                    self._organize_files(srna_prefix, rnaplex_path, "_RNAplex", folder)
                    self._organize_files(srna_prefix, rnaup_path, "_RNAup", folder)
                    self._organize_files(srna_prefix, merge_path, "_merge", folder)

    def run_sRNA_target_prediction(self, vienna_path, gffs, fastas, srnas, query, program, 
                                   inter_length, win_size_t, span_t, win_size_s, span_s, 
                                   unstr_region_rnaplex_t, unstr_region_rnaplex_s, 
                                   unstr_region_rnaup, energy, duplex_dist, top, 
                                   out_folder, core_plex, core_up):
        self._check_gff(gffs)
        self._check_gff(srnas)
        multiparser._parser_gff(gffs, None)
        multiparser._parser_fasta(fastas)
        multiparser._parser_gff(srnas, "sRNA")
        prefixs = []
        self._gen_seq(prefixs, query)
        if (program == "both") or \
           (program == "RNAplex"):
            self._rna_plex(prefixs, self.rnaplex_path, vienna_path, win_size_s, span_s, 
                           unstr_region_rnaplex_s, self.srna_seq_path, win_size_t, span_t, 
                           unstr_region_rnaplex_t, self.target_seq_path, inter_length, energy, 
                           duplex_dist, core_plex)
        self.remove_all_content(self.target_seq_path, "_target_", "file")
        if (program == "both") or \
           (program == "RNAup"):
            self._rnaup(prefixs, self.rnaup_path, self.srna_seq_path, 
                        self.target_seq_path, core_up, vienna_path, unstr_region_rnaup)
        self._merge_rnaplex_rnaup(prefixs, self.merge_path, program, self.rnaplex_path, 
                                  self.rnaup_path, self.srna_path, self.gff_path)
        self._organize_gff(gffs, rnaplex_path, rnaup_path, merge_path)
        if (program == "RNAplex") or \
           (program == "both"):
            for strain in os.listdir(os.path.join(out_folder, "RNAplex")):
                shutil.rmtree(os.path.join(out_folder, "RNAplex", strain, "RNAplfold"))
        self.helper.remove_all_content(os.getcwd(), self.tmp, "dir")
        self._remove_tmp(gffs)
        self._remove_tmp(srnas)
        self._remove_tmp(fastas)
