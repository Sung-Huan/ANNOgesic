import os
import sys
import time
import shutil
from subprocess import call, Popen
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.converter import Converter
from annogesiclib.circRNA import detect_circrna


class CircRNADetection(object):
    '''Detection of circRNA'''

    def __init__(self, args_circ):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.converter = Converter()
        self.alignment_path = os.path.join(args_circ.output_folder,
                                           "segemehl_align")
        self.splice_path = os.path.join(args_circ.output_folder,
                                        "segemehl_splice")
        self.candidate_path = os.path.join(args_circ.output_folder,
                                           "circRNA_tables")
        self.gff_folder = os.path.join(args_circ.output_folder, "gffs")
        self.gff_path = os.path.join(args_circ.gffs, "tmp")
        self.splices = {"all_file": "splicesites_all.bed",
                        "file": "splicesites.bed",
                        "all": "splicesites_all", "splice": "splicesites"}
        self.trans = {"all_file": "transrealigned_all.bed",
                      "file": "transrealigned.bed",
                      "all": "transrealigned_all", "trans": "transrealigned"}
        self.bams = {"whole": "whole_reads.bam", "sort": "whole_reads_sort"}
        if args_circ.align:
            if args_circ.fastas is None:
                print("Error: There is no genome fasta file!!!")
                sys.exit()
            else:
                self.fasta_path = os.path.join(args_circ.fastas, "tmp")
        else:
            self.fasta_path = os.path.join(args_circ.fastas, "tmp")

    def _wait_process(self, processes):
        '''wait for the parallels to finish the process'''
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

    def _deal_zip_file(self, read_folder):
        tmp_reads = []
        for read in os.listdir(read_folder):
            if read.endswith(".bz2"):
                mod_read = read.replace(".bz2", "")
                if (".fa" not in mod_read) and (
                        ".fasta" not in mod_read) and (
                        ".fna" not in mod_read):
                    mod_read = mod_read + ".fa"
                read_out = open(os.path.join(read_folder, mod_read), "w")
                tmp_reads.append(os.path.join(read_folder, mod_read))
                print(" ".join(["unzip", read]))
                call(["bzcat", os.path.join(read_folder, read)],
                     stdout=read_out)
                read_out.close()
            elif read.endswith(".gz"):
                mod_read = read.replace(".gz", "")
                if (".fa" not in mod_read) and (
                        ".fasta" not in mod_read) and (
                        ".fna" not in mod_read):
                    mod_read = mod_read + ".fa"
                read_out = open(os.path.join(read_folder, mod_read), "w")
                tmp_reads.append(os.path.join(read_folder, mod_read))
                print(" ".join(["unzip", read]))
                call(["zcat", os.path.join(read_folder, read)],
                     stdout=read_out)
                read_out.close()
        return tmp_reads

    def _run_segemehl_fasta_index(self, segemehl_path, fasta_path,
                                  index, fasta):
        call([os.path.join(segemehl_path, "segemehl.x"),
              "-x", os.path.join(fasta_path, index),
              "-d", os.path.join(fasta_path, fasta)])

    def _run_segemehl_align(self, args_circ, index, fasta, read,
                            sam_file, log_file, fasta_prefix):
        out = open(os.path.join(self.alignment_path,
                   fasta_prefix, sam_file), "w")
        log = open(os.path.join(self.alignment_path,
                   fasta_prefix, log_file), "w")
        p = Popen([os.path.join(args_circ.segemehl_path, "segemehl.x"),
                   "-i", os.path.join(self.fasta_path, index),
                   "-d", os.path.join(self.fasta_path, fasta),
                   "-q", os.path.join(args_circ.read_folder, read), "-S"],
                  stdout=out, stderr=log)
        return p

    def _align(self, args_circ):
        '''align the read. if the bam files are provided, it can be skipped.'''
        prefixs = []
        align_files = []
        for fasta in os.listdir(self.fasta_path):
            index = fasta.replace(".fa", ".idx")
            self._run_segemehl_fasta_index(args_circ.segemehl_path,
                                           self.fasta_path, index, fasta)
            processes = []
            num_process = 0
            fasta_prefix = fasta.replace(".fa", "")
            prefixs.append(fasta_prefix)
            self.helper.check_make_folder(os.path.join(
                                self.alignment_path, fasta_prefix))
            for read in os.listdir(args_circ.read_folder):
                num_process += 1
                if read.endswith(".fa") or \
                   read.endswith(".fna") or \
                   read.endswith("fasta"):
                    filename = read.split(".")
                    read_prefix = ".".join(filename[:-1])
                    sam_file = "_".join([read_prefix, fasta_prefix + ".sam"])
                    log_file = "_".join([read_prefix, fasta_prefix + ".log"])
                    align_files.append("_".join([read_prefix, fasta_prefix]))
                    print("mapping {0}".format(sam_file))
                    p = self._run_segemehl_align(
                            args_circ, index, fasta, read,
                            sam_file, log_file, fasta_prefix)
                    processes.append(p)
                    if num_process == args_circ.cores:
                        self._wait_process(processes)
                        num_process = 0
            self._wait_process(processes)
        return align_files, prefixs

    def _run_samtools_convert_bam(self, samtools_path, pre_sam, out_bam):
        call([samtools_path, "view", "-bS", pre_sam, "-o", out_bam])

    def _convert_sam2bam(self, sub_alignment_path, samtools_path, align_files):
        bam_files = []
        convert_ones = []
        remove_ones = []
        for sam in os.listdir(sub_alignment_path):
            pre_sam = os.path.join(sub_alignment_path, sam)
            if sam.endswith(".sam"):
                bam_file = sam.replace(".sam", ".bam")
                print("Convert {0} to {1}".format(sam, bam_file))
                out_bam = os.path.join(sub_alignment_path, bam_file)
                self._run_samtools_convert_bam(samtools_path, pre_sam, out_bam)
                bam_files.append(out_bam)
                if align_files:
                    if bam_file.replace(".bam", "") not in align_files:
                        convert_ones.append(out_bam)
                    else:
                        remove_ones.append(pre_sam)
            elif sam.endswith(".bam"):
                if (pre_sam not in convert_ones) and (
                        pre_sam not in remove_ones):
                    bam_files.append(pre_sam)
            elif sam.endswith(".log"):
                os.remove(pre_sam)
        return bam_files, convert_ones, remove_ones

    def _run_samtools_merge_sort(self, samtools_path,
                                 sub_alignment_path, bam_files):
        print("Merge all bam files....")
        whole_bam = os.path.join(sub_alignment_path, self.bams["whole"])
        if len(bam_files) <= 1:
            shutil.copyfile(bam_files[0], whole_bam)
        else:
            file_line = " ".join(bam_files)
            os.system(" ".join([samtools_path, "merge",
                                whole_bam, file_line]))
        print("Sort bam files....")
        call([samtools_path, "sort", "-o", os.path.join(sub_alignment_path,
              self.bams["sort"] + ".bam"), whole_bam])
        os.remove(os.path.join(sub_alignment_path, self.bams["whole"]))

    def _run_samtools_convert_sam(self, samtools_path, sub_alignment_path):
        print("Convert whole reads bam file to sam file....")
        call([samtools_path, "view", "-h", "-o",
              os.path.join(sub_alignment_path, self.bams["sort"] + ".sam"),
              os.path.join(sub_alignment_path, self.bams["sort"] + ".bam")])

    def _merge_sort_aligment_file(self, bam_files, samtools_path,
                                  sub_alignment_path, convert_ones,
                                  tmp_reads, remove_ones):
        self._run_samtools_merge_sort(samtools_path,
                                      sub_alignment_path, bam_files)
        self._run_samtools_convert_sam(samtools_path, sub_alignment_path)
        for bam in convert_ones:
            os.remove(bam)
        for sam in remove_ones:
            os.remove(sam)
        if len(tmp_reads) != 0:
            for read in tmp_reads:
                os.remove(read)

    def _run_testrealign(self, prefix, segemehl_path, sub_alignment_path):
        self.helper.check_make_folder(os.path.join(self.splice_path, prefix))
        sub_splice_path = os.path.join(self.splice_path, prefix)
        err_log = os.path.join(sub_splice_path, prefix + ".log")
        print("Running testrealign.x for {0}".format(prefix))
        command = " ".join([
                  os.path.join(segemehl_path, "testrealign.x"),
                  "-d", os.path.join(self.fasta_path, prefix + ".fa"),
                  "-q", os.path.join(sub_alignment_path,
                                     self.bams["sort"] + ".sam"),
                  "-n"])
        os.system(command + " 2>" + err_log)
        self.helper.move_all_content(os.getcwd(), sub_splice_path, [".bed"])
        self.helper.remove_all_content(sub_alignment_path,
                                       self.bams["sort"], "file")

    def _merge_bed(self, fastas, splice_path):
        '''Merge the bed files for analysis'''
        tmp_prefixs = []
        for fasta in os.listdir(fastas):
            headers = []
            if (fasta.endswith(".fa") or fasta.endswith(".fna") or
                    fasta.endswith(".fasta")):
                with open(os.path.join(fastas, fasta), "r") as f_h:
                    for line in f_h:
                        line = line.strip()
                        if line.startswith(">"):
                            headers.append(line[1:])
                filename = fasta.split(".")
                fasta_prefix = ".".join(filename[:-1])
                tmp_prefixs.append(fasta_prefix)
                self.helper.check_make_folder(os.path.join(
                                              os.getcwd(), fasta_prefix))
                for header in headers:
                    shutil.copyfile(os.path.join(splice_path, header,
                                    self.splices["file"]),
                                    os.path.join(fasta_prefix,
                                    "_".join([self.splices["splice"],
                                              header + ".bed"])))
                    shutil.copyfile(os.path.join(splice_path, header,
                                    self.trans["file"]),
                                    os.path.join(fasta_prefix,
                                    "_".join([self.trans["trans"],
                                              header + ".bed"])))
                out_splice = os.path.join(fasta_prefix,
                                          self.splices["all_file"])
                out_trans = os.path.join(fasta_prefix,
                                         self.trans["all_file"])
                if len(headers) > 1:
                    for file_ in os.listdir(fasta_prefix):
                        if (self.splices["splice"] in file_) and (
                                self.splices["all"] not in file_):
                            self.helper.merge_file(os.path.join(
                                    fasta_prefix, file_), out_splice)
                        elif (self.trans["trans"] in file_) and (
                                self.trans["all"] not in file_):
                            self.helper.merge_file(os.path.join(
                                    fasta_prefix, file_), out_trans)
                else:
                    shutil.move(os.path.join(
                                fasta_prefix,
                                "_".join([self.splices["splice"],
                                         headers[0] + ".bed"])),
                                out_splice)
                    shutil.move(os.path.join(
                                fasta_prefix,
                                "_".join([self.trans["trans"],
                                          headers[0] + ".bed"])),
                                out_trans)
        self.helper.remove_all_content(splice_path, None, "dir")
        return tmp_prefixs

    def _stat_and_gen_gff(self, tmp_prefixs, args_circ):
        '''do statistics and print the result to gff file'''
        for prefix in tmp_prefixs:
            self.helper.check_make_folder(os.path.join(self.gff_folder,
                                                       prefix))
            shutil.copytree(prefix, os.path.join(self.splice_path, prefix))
            self.helper.check_make_folder(os.path.join(
                                          self.candidate_path, prefix))
            print("comparing with annotation of {0}".format(prefix))
            if self.splices["all_file"] in os.listdir(os.path.join(
                                           self.splice_path, prefix)):
                detect_circrna(os.path.join(self.splice_path, prefix,
                               self.splices["all_file"]), os.path.join(
                               self.gff_path, prefix + ".gff"),
                               os.path.join(self.candidate_path, prefix,
                               "_".join(["circRNA", prefix + "_all.csv"])),
                               args_circ, os.path.join(args_circ.stat_folder,
                               "_".join(["stat_circRNA", prefix + ".csv"])))
                self.converter.convert_circ2gff(
                     os.path.join(self.candidate_path, prefix,
                                  "_".join(["circRNA",
                                            prefix + "_all.csv"])),
                     args_circ, os.path.join(
                                self.gff_folder, prefix,
                                "_".join([prefix, "circRNA_all.gff"])),
                     os.path.join(self.gff_folder, prefix,
                                  "_".join([prefix, "circRNA_best.gff"])))

    def _assign_merge_bam(self, args_circ):
        '''Check and merge bam files'''
        remove_frags = []
        bam_files = []
        if (args_circ.normal_bams is not None) and (
                args_circ.frag_bams is not None):
            for frag in os.listdir(args_circ.frag_bams):
                if frag.endswith(".bam"):
                    shutil.copyfile(os.path.join(args_circ.frag_bams, frag),
                                    os.path.join(args_circ.normal_bams, frag))
                    remove_frags.append(frag)
            merge_folder = args_circ.normal_bams
        elif (args_circ.normal_bams is not None):
            merge_folder = args_circ.normal_bams
        elif (args_circ.frag_bams is not None):
            merge_folder = args_circ.frag_bams
        else:
            print("Error: please assign bam folder or do alignment!!")
            sys.exit()
        for bam in os.listdir(merge_folder):
            if bam.endswith(".bam"):
                bam_files.append(os.path.join(merge_folder, bam))
        return merge_folder, remove_frags, bam_files

    def run_circrna(self, args_circ):
        '''detection of circRNA'''
        for gff in os.listdir(args_circ.gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(
                                                 args_circ.gffs, gff))
        if args_circ.segemehl_path is None:
            print("Error: please assign segemehl folder!!")
            sys.exit()
        self.multiparser.parser_gff(args_circ.gffs, None)
        self.multiparser.combine_gff(args_circ.fastas, self.gff_path,
                                     "fasta", None)
        tmp_reads = []
        if args_circ.align:
            self.multiparser.parser_fasta(args_circ.fastas)
            tmp_reads = self._deal_zip_file(args_circ.read_folder)
            align_files, prefixs = self._align(args_circ)
        else:
            self.multiparser.parser_fasta(args_circ.fastas)
            prefixs = []
            for fasta in os.listdir(self.fasta_path):
                fasta_prefix = fasta.replace(".fa", "")
                prefixs.append(fasta_prefix)
            merge_folder, remove_frag, bam_files = self._assign_merge_bam(
                                                   args_circ)
            align_files = None
        for prefix in prefixs:
            if args_circ.align:
                sub_alignment_path = os.path.join(self.alignment_path, prefix)
                bam_files, convert_ones, remove_ones = self._convert_sam2bam(
                    sub_alignment_path, args_circ.samtools_path, align_files)
            else:
                sub_alignment_path = merge_folder
                convert_ones = []
                remove_ones = []
            self._merge_sort_aligment_file(
                bam_files, args_circ.samtools_path, sub_alignment_path,
                convert_ones, tmp_reads, remove_ones)
            self._run_testrealign(prefix, args_circ.segemehl_path,
                                  sub_alignment_path)
        tmp_prefixs = self._merge_bed(args_circ.fastas, self.splice_path)
        self.multiparser.parser_gff(args_circ.gffs, None)
        self.multiparser.combine_gff(args_circ.fastas, self.gff_path,
                                     "fasta", None)
        self._stat_and_gen_gff(tmp_prefixs, args_circ)
        self.helper.remove_tmp(args_circ.fastas)
        self.helper.remove_tmp(args_circ.gffs)
        for tmp_prefix in tmp_prefixs:
            shutil.rmtree(tmp_prefix)
        if (not args_circ.align) and (len(remove_frag) != 0):
            for frag in remove_frag:
                os.remove(os.path.join(merge_folder, frag))
