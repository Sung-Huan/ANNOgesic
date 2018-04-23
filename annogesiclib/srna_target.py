import os
import shutil
import sys
import time
from subprocess import Popen, call
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.potential_target import potential_target
from annogesiclib.format_fixer import FormatFixer
from annogesiclib.merge_rnaplex_rnaup import merge_srna_target
from annogesiclib.gff3 import Gff3Parser


class sRNATargetPrediction(object):
    '''detection of sRNA-target interaction'''

    def __init__(self, args_tar):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.fixer = FormatFixer()
        self.gff_parser = Gff3Parser()
        self.target_seq_path = os.path.join(args_tar.out_folder, "target_seqs")
        self.srna_seq_path = os.path.join(args_tar.out_folder, "sRNA_seqs")
        self.rnaplex_path = os.path.join(args_tar.out_folder, "RNAplex_results")
        self.rnaup_path = os.path.join(args_tar.out_folder, "RNAup_results")
        self.intarna_path = os.path.join(args_tar.out_folder, "IntaRNA_results")
        self.merge_path = os.path.join(args_tar.out_folder, "merged_results")
        self.srna_path = os.path.join(args_tar.srnas, "tmp")
        self.fasta_path = os.path.join(args_tar.fastas, "tmp")
        self.gff_path = os.path.join(args_tar.gffs, "tmp")
        self.tmps = {"tmp": "tmp_srna_target", "rnaup": "tmp_rnaup",
                     "log": "tmp_log",
                     "all_fa": "tmp*.fa", "all_txt": "tmp*.txt"}

    def _check_gff(self, gffs):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(gffs, gff))

    def _run_rnaplfold(self, rnaplfold_path, file_type, win_size, span,
                       unstr_region, seq_path, prefix, out_path, log):
        current = os.getcwd()
        os.chdir(out_path)
        command = " ".join([rnaplfold_path,
                            "-W", str(win_size),
                            "-L", str(span),
                            "-u", str(unstr_region),
                            "-O"])
        if file_type == "sRNA":
            log.write("<".join([command, os.path.join(current, seq_path,
                                "_".join([self.tmps["tmp"], prefix,
                                          file_type + ".fa"]))]) + "\n")
            os.system("<".join([command, os.path.join(current, seq_path,
                                "_".join([self.tmps["tmp"], prefix,
                                          file_type + ".fa"]))]))
        else:
            log.write("<".join([command, os.path.join(current, seq_path,
                                "_".join([prefix, file_type + ".fa"]))]) + "\n")
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

    def _sort_srna_fasta(self, fasta, prefix, path):
        out = open(os.path.join(path,
                   "_".join([self.tmps["tmp"], prefix, "sRNA.fa"])), "w")
        srnas = []
        with open(fasta) as f_h:
            for line in f_h:
                line = line.strip()
                if line.startswith(">"):
                    name = line[1:]
                else:
                    srnas.append({"name": name, "seq": line, "len": len(line)})
        srnas = sorted(srnas, key=lambda x: (x["len"]))
        for srna in srnas:
            out.write(">" + srna["name"].split("|")[0] + "\n")
            out.write(srna["seq"] + "\n")
        out.close()

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

    def _get_specific_seq(self, srna_file, seq_file, srna_out, querys):
        for query in querys:
            srna_datas = query.split(":")
            srna = {"seq_id": srna_datas[0], "strand": srna_datas[3],
                    "start": int(srna_datas[1]), "end": int(srna_datas[2])}
            gff_f = open(srna_file, "r")
            out = open(srna_out, "a")
            seq = self._read_fasta(seq_file)
            num = 0
            detect = False
            for entry in self.gff_parser.entries(gff_f):
                if (entry.seq_id == srna["seq_id"]) and (
                        entry.strand == srna["strand"]) and (
                        entry.start == srna["start"]) and (
                        entry.end == srna["end"]):
                    detect = True
                    if "ID" in entry.attributes.keys():
                        id_ = entry.attributes["ID"]
                    else:
                        id_ = entry.feature + str(num)
                    gene = self.helper.extract_gene(seq, entry.start,
                                                    entry.end, entry.strand)
                    out.write(">{0}|{1}|{2}|{3}|{4}\n{5}\n".format(
                              id_, entry.seq_id, entry.start,
                              entry.end, entry.strand, gene))
                    num += 1
            if not detect:
                print("Error: Some of the query sRNAs do not exist!")
                sys.exit()
            gff_f.close()
            out.close()

    def _gen_seq(self, prefixs, args_tar):
        print("Generating sRNA fasta files")
        for srna in os.listdir(self.srna_path):
            if srna.endswith("_sRNA.gff"):
                prefix = srna.replace("_sRNA.gff", "")
                prefixs.append(prefix)
                srna_out = os.path.join(self.srna_seq_path,
                                        "_".join([prefix, "sRNA.fa"]))
                if "all" in args_tar.query:
                    self.helper.get_seq(
                            os.path.join(self.srna_path, srna),
                            os.path.join(self.fasta_path, prefix + ".fa"),
                            srna_out)
                else:
                    if "_".join([prefix, "sRNA.fa"]) in os.listdir(
                       self.srna_seq_path):
                        os.remove(srna_out)
                    self._get_specific_seq(
                            os.path.join(self.srna_path, srna),
                            os.path.join(self.fasta_path, prefix + ".fa"),
                            srna_out, args_tar.query)
                self._sort_srna_fasta(srna_out, prefix, self.srna_seq_path)
        print("Generating target fasta files")
        for gff in os.listdir(self.gff_path):
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                potential_target(os.path.join(self.gff_path, gff),
                                 os.path.join(self.fasta_path, prefix + ".fa"),
                                 os.path.join(self.target_seq_path), args_tar)
                file_num = 1
                num = 0
                sub_prefix = os.path.join(self.target_seq_path,
                                          "_".join([prefix, "target"]))
                sub_out = open("_".join([sub_prefix, str(file_num) + ".fa"]),
                               "w")
                with open((sub_prefix + ".fa"), "r") as t_f:
                    for line in t_f:
                        line = line.strip()
                        if line.startswith(">"):
#                            line = line.replace("|", "_")
                            num += 1
                        if (num == 100):
                            num = 0
                            file_num += 1
                            sub_out.close()
                            sub_out = open("_".join([sub_prefix,
                                           str(file_num) + ".fa"]), "w")
                        sub_out.write(line + "\n")
                sub_out.close()

    def _run_rnaplex(self, prefix, rnaplfold_folder, args_tar, log):
        print("Running RNAplex of {0}".format(prefix))
        num_process = 0
        processes = []
        for seq in os.listdir(self.target_seq_path):
            if (prefix in seq) and ("_target_" in seq):
                print("Running RNAplex with {0}".format(seq))
                out_rnaplex = open(os.path.join(
                    self.rnaplex_path, prefix, "_".join([
                        prefix, "RNAplex", str(num_process) + ".txt"])), "w")
                num_process += 1
                log.write(" ".join([args_tar.rnaplex_path,
                           "-q", os.path.join(
                               self.srna_seq_path, "_".join([
                                   self.tmps["tmp"], prefix, "sRNA.fa"])),
                           "-t", os.path.join(self.target_seq_path, seq),
                           "-l", str(args_tar.inter_length),
                           "-e", str(args_tar.energy),
                           "-z", str(args_tar.duplex_dist),
                           "-a", rnaplfold_folder]) + "\n")
                p = Popen([args_tar.rnaplex_path,
                           "-q", os.path.join(
                               self.srna_seq_path, "_".join([
                                   self.tmps["tmp"], prefix, "sRNA.fa"])),
                           "-t", os.path.join(self.target_seq_path, seq),
                           "-l", str(args_tar.inter_length),
                           "-e", str(args_tar.energy),
                           "-z", str(args_tar.duplex_dist),
                           "-a", rnaplfold_folder], stdout=out_rnaplex)
                processes.append(p)
                if num_process % args_tar.core_plex == 0:
                    self._wait_process(processes)
        self._wait_process(processes)
        log.write("The prediction for {0} is done.\n".format(prefix))
        log.write("The following temporary files for storing results of {0} are "
                  "generated:\n".format(prefix))
        for file_ in os.listdir(os.path.join(self.rnaplex_path, prefix)):
            log.write("\t" + os.path.join(self.rnaplex_path, prefix, file_) + "\n")
        return num_process

    def _rna_plex(self, prefixs, args_tar, log):
        log.write("Using RNAplex and RNAplfold to predict sRNA targets.\n")
        log.write("Please make sure the version of Vienna RNA package is "
                  "at least 2.3.2.\n")
        for prefix in prefixs:
            print("Running RNAplfold of {0}".format(prefix))
            self.helper.check_make_folder(
                        os.path.join(self.rnaplex_path, prefix))
            rnaplfold_folder = os.path.join(self.rnaplex_path, prefix,
                                          "RNAplfold")
            os.mkdir(rnaplfold_folder)
            self._run_rnaplfold(
                args_tar.rnaplfold_path, "sRNA", args_tar.win_size_s,
                args_tar.span_s, args_tar.unstr_region_rnaplex_s,
                self.srna_seq_path, prefix, rnaplfold_folder, log)
            self._run_rnaplfold(
                args_tar.rnaplfold_path, "target", args_tar.win_size_t,
                args_tar.span_t, args_tar.unstr_region_rnaplex_t,
                self.target_seq_path, prefix, rnaplfold_folder, log)
            num_process = self._run_rnaplex(prefix, rnaplfold_folder, args_tar, log)
            rnaplex_file = os.path.join(self.rnaplex_path, prefix,
                                        "_".join([prefix, "RNAplex.txt"]))
            if ("_".join([prefix, "RNAplex.txt"]) in
                    os.listdir(os.path.join(self.rnaplex_path, prefix))):
                os.remove(rnaplex_file)
            for index in range(0, num_process):
                log.write("Using helper.py to merge the temporary files.\n")
                self.helper.merge_file(os.path.join(
                    self.rnaplex_path, prefix, "_".join([
                        prefix, "RNAplex", str(index) + ".txt"])),
                    rnaplex_file)
            log.write("\t" + rnaplex_file + " is generated.\n")
            self.helper.remove_all_content(os.path.join(
                 self.rnaplex_path, prefix), "_RNAplex_", "file")
            self.fixer.fix_rnaplex(rnaplex_file, self.tmps["tmp"])
            shutil.move(self.tmps["tmp"], rnaplex_file)
            shutil.rmtree(rnaplfold_folder)

    def _run_rnaup(self, num_up, processes, prefix, out_rnaup, out_log,
                   args_tar, log):
        for index in range(1, num_up + 1):
            out_tmp_up = open(os.path.join(
                args_tar.out_folder, "".join([self.tmps["rnaup"],
                                              str(index), ".txt"])), "w")
            out_err = open(os.path.join(
                args_tar.out_folder, "".join([self.tmps["log"],
                                              str(index), ".txt"])), "w")
            in_up = open(os.path.join(
                args_tar.out_folder, "".join([self.tmps["tmp"],
                                              str(index), ".fa"])), "r")
            log.write(" ".join([args_tar.rnaup_path,
                       "-u", str(args_tar.unstr_region_rnaup),
                       "-o", "--interaction_first"]) + "\n")
            p = Popen([args_tar.rnaup_path,
                       "-u", str(args_tar.unstr_region_rnaup),
                       "-o", "--interaction_first"],
                      stdin=in_up, stdout=out_tmp_up, stderr=out_err)
            processes.append(p)
        if len(processes) != 0:
            time.sleep(5)
            self._wait_process(processes)
            log.write("The following temporary files for storing results of {0} are "
                      "generated:\n".format(prefix))
            for file_ in os.listdir(os.path.join(args_tar.out_folder)):
                log.write("\t" + os.path.join(args_tar.out_folder, file_) + "\n")
            os.system("rm " + os.path.join(args_tar.out_folder,
                                           self.tmps["all_fa"]))
            self._merge_txt(num_up, out_rnaup, out_log, args_tar.out_folder)
            os.system("rm " + os.path.join(args_tar.out_folder,
                                           self.tmps["all_txt"]))

    def _merge_txt(self, num_up, out_rnaup, out_log, out_folder):
        for index in range(1, num_up + 1):
            self.helper.merge_file(
                os.path.join(out_folder, "".join([self.tmps["rnaup"],
                                                  str(index), ".txt"])),
                out_rnaup)
            self.helper.merge_file(
                os.path.join(out_folder, "".join([self.tmps["log"],
                                                  str(index), ".txt"])),
                out_log)

    def _get_continue(self, out_rnaup):
        '''For RNAup, it can continue running RNAup based on previous run'''
        srnas = []
        matchs = {}
        out = open("tmp.txt", "w")
        with open(out_rnaup) as f_h:
            for line in f_h:
                line = line.strip()
                if ">srna" in line:
                    srna = line[1:]
                    srnas.append(srna)
                    matchs[srna] = []
                else:
                    matchs[srna].append(line)
        srnas = srnas[:-1]
        for srna in srnas:
            out.write(">" + srna + "\n")
            for target in matchs[srna]:
                out.write(target + "\n")
        out.close()
        os.remove(out_rnaup)
        shutil.move("tmp.txt", out_rnaup)
        return srnas

    def _rnaup(self, prefixs, args_tar, log):
        log.write("Using RNAup to predict sRNA targets.\n")
        log.write("Please make sure the version of Vienna RNA package is "
                  "at least 2.3.2.\n")
        for prefix in prefixs:
            srnas = []
            print("Running RNAup of {0}".format(prefix))
            if not os.path.exists(os.path.join(self.rnaup_path, prefix)):
                os.mkdir(os.path.join(self.rnaup_path, prefix))
            num_up = 0
            processes = []
            out_rnaup = os.path.join(self.rnaup_path, prefix,
                                     "_".join([prefix + "_RNAup.txt"]))
            out_log = os.path.join(self.rnaup_path, prefix,
                                   "_".join([prefix + "_RNAup.log"]))
            if "_".join([prefix, "RNAup.txt"]) in \
                    os.listdir(os.path.join(self.rnaup_path, prefix)):
                if not args_tar.continue_rnaup:
                    os.remove(out_rnaup)
                    os.remove(out_log)
                else:
                    log.write("The data from the previous run is found.\n")
                    srnas = self._get_continue(out_rnaup)
                    log.write("The previous data is loaded.\n")
            with open(os.path.join(self.srna_seq_path, "_".join([
                    self.tmps["tmp"], prefix, "sRNA.fa"])), "r") as s_f:
                for line in s_f:
                    line = line.strip()
                    if line.startswith(">"):
                        if line[1:] in srnas:
                            start = False
                            continue
                        start = True
                        print("Running RNAup with {0}".format(line[1:]))
                        num_up += 1
                        out_up = open(os.path.join(args_tar.out_folder,
                                      "".join([self.tmps["tmp"],
                                               str(num_up), ".fa"])), "w")
                        out_up.write(line + "\n")
                    else:
                        if start:
                            out_up.write(line + "\n")
                            out_up.close()
                            self.helper.merge_file(os.path.join(
                                self.target_seq_path,
                                "_".join([prefix, "target.fa"])),
                                os.path.join(args_tar.out_folder,
                                             "".join([self.tmps["tmp"],
                                                      str(num_up), ".fa"])))
                            if num_up == args_tar.core_up:
                                self._run_rnaup(num_up, processes, prefix,
                                                out_rnaup, out_log, args_tar, log)
                                processes = []
                                num_up = 0
            self._run_rnaup(num_up, processes, prefix, out_rnaup, out_log,
                            args_tar, log)
            log.write("The prediction for {0} is done.\n".format(prefix))
            log.write("\t" + out_rnaup + " is complete generated and updated.\n")

    def _intarna(self, prefixs, args_tar, log):
        log.write("Using IntaRNA to predict sRNA targets.\n")
        log.write("Please make sure the version of IntaRNA is at least 2.0.4.\n")
        for prefix in prefixs:
            print("Running IntaRNA of {0}".format(prefix))
            intarna_file = os.path.join(self.intarna_path, prefix,
                                        prefix + "_IntaRNA.txt")
            self.helper.check_make_folder(
                        os.path.join(self.intarna_path, prefix))
            call([args_tar.intarna_path,
                  "-q", os.path.join(
                      self.srna_seq_path, "_".join([
                          self.tmps["tmp"], prefix, "sRNA.fa"])),
                  "-t", os.path.join(self.target_seq_path,
                                     prefix + "_target.fa"),
                  "--qAccW", str(args_tar.slide_win_srna),
                  "--qAccL", str(args_tar.max_loop_srna),
                  "--tAccW", str(args_tar.slide_win_target),
                  "--tAccL", str(args_tar.max_loop_target),
                  "--outMode", "C", "-m", args_tar.mode_intarna,
                  "--threads", str(args_tar.core_inta),
                  "--out", intarna_file])
            log.write("The prediction for {0} is done.\n".format(prefix))
            log.write("\t" + intarna_file + " is generated.\n")

    def _merge_rnaplex_rnaup(self, prefixs, args_tar, log):
        '''merge the result of IntaRNA, RNAup and RNAplex'''
        log.write("Running merge_rnaplex_rnaup.py to merge the results from "
                  "RNAplex, RNAup, and IntaRNA for generating finanl output.\n")
        log.write("The following files are generated:\n")
        for prefix in prefixs:
            rnaplex_file = None
            rnaup_file = None
            out_rnaplex = None
            out_rnaup = None
            intarna_file = None
            out_intarna = None
            self.helper.check_make_folder(os.path.join(
                                          self.merge_path, prefix))
            print("Ranking {0} now".format(prefix))
            if ("RNAplex" in args_tar.program):
                rnaplex_file = os.path.join(self.rnaplex_path, prefix,
                                            "_".join([prefix, "RNAplex.txt"]))
                out_rnaplex = os.path.join(
                        self.rnaplex_path, prefix,
                        "_".join([prefix, "RNAplex_rank.csv"]))
            if ("RNAup" in args_tar.program):
                rnaup_file = os.path.join(self.rnaup_path, prefix,
                                          "_".join([prefix, "RNAup.txt"]))
                out_rnaup = os.path.join(self.rnaup_path, prefix,
                                         "_".join([prefix, "RNAup_rank.csv"]))
            if ("IntaRNA" in args_tar.program):
                intarna_file = os.path.join(self.intarna_path, prefix,
                                            "_".join([prefix, "IntaRNA.txt"]))
                out_intarna = os.path.join(self.intarna_path, prefix,
                                           "_".join([prefix, "IntaRNA_rank.csv"]))
            overlap_file = os.path.join(self.merge_path, prefix,
                                        "_".join([prefix, "overlap.csv"]))
            merge_file = os.path.join(self.merge_path, prefix,
                                      "_".join([prefix, "merge.csv"]))
            merge_srna_target(rnaplex_file, rnaup_file, intarna_file, args_tar,
                              out_rnaplex, out_rnaup, out_intarna,
                              os.path.join(self.fasta_path, prefix + ".fa"),
                              merge_file, overlap_file,
                              os.path.join(self.srna_path,
                                           "_".join([prefix, "sRNA.gff"])),
                              os.path.join(self.gff_path, prefix + ".gff"))
            if ("RNAplex" in args_tar.program):
                log.write("\t" + out_rnaplex + "\n")
            if ("RNAup" in args_tar.program):
                log.write("\t" + out_rnaup + "\n")
            if ("IntaRNA" in args_tar.program):
                log.write("\t" + out_intarna + "\n")
            if (os.path.exists(merge_file)):
                log.write("\t" + merge_file + "\n")
            if (os.path.exists(overlap_file)):
                log.write("\t" + overlap_file + "\n")

    def run_srna_target_prediction(self, args_tar, log):
        self._check_gff(args_tar.gffs)
        self._check_gff(args_tar.srnas)
        self.multiparser.parser_gff(args_tar.gffs, None)
        self.multiparser.parser_fasta(args_tar.fastas)
        self.multiparser.parser_gff(args_tar.srnas, "sRNA")
        prefixs = []
        self._gen_seq(prefixs, args_tar)
        if ("RNAplex" in args_tar.program):
            self._rna_plex(prefixs, args_tar, log)
        self.helper.remove_all_content(self.target_seq_path,
                                       "_target_", "file")
        log.write("The temporary files for running RNAplex are deleted.\n")
        if ("RNAup" in args_tar.program):
            self._rnaup(prefixs, args_tar, log)
        if ("IntaRNA" in args_tar.program):
            self._intarna(prefixs, args_tar, log)
        self._merge_rnaplex_rnaup(prefixs, args_tar, log)
        self.helper.remove_all_content(args_tar.out_folder,
                                       self.tmps["tmp"], "dir")
        self.helper.remove_all_content(args_tar.out_folder,
                                       self.tmps["tmp"], "file")
        self.helper.remove_tmp_dir(args_tar.gffs)
        self.helper.remove_tmp_dir(args_tar.srnas)
        self.helper.remove_tmp_dir(args_tar.fastas)
        self.helper.remove_all_content(self.srna_seq_path, "tmp_", "file")
