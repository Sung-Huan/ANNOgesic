import os
import sys
import time
import shutil
import copy
from glob import glob
from subprocess import call, Popen
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.converter import Converter
from annogesiclib.circRNA_detection import detect_circrna


class CircRNADetection(object):
    '''Detection of circRNA'''

    def __init__(self, args_circ):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.converter = Converter()
        self.alignment_path = os.path.join(args_circ.output_folder,
                                           "segemehl_alignment_files")
        self.splice_path = os.path.join(args_circ.output_folder,
                                        "segemehl_splice_results")
        self.candidate_path = os.path.join(args_circ.output_folder,
                                           "circRNA_tables")
        self.gff_folder = os.path.join(args_circ.output_folder, "gffs")
        self.gff_path = os.path.join(args_circ.gffs, "tmp")
        self.splices = {"file": "splicesites.bed",
                        "splice": "splicesites"}
        self.trans = {"file": "transrealigned.bed",
                      "trans": "transrealigned"}
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

    def _deal_zip_file(self, read_files, log):
        tmp_datas = []
        tmp_reads = []
        for reads in read_files:
            zips = []
            tmp_datas = reads["files"]
            for read in reads["files"]:
                if read.endswith(".bz2"):
                    mod_read = read.replace(".bz2", "")
                    if (".fa" not in mod_read) and (
                            ".fasta" not in mod_read) and (
                            ".fna" not in mod_read) and (
                            ".fq" not in mod_read) and (
                            ".fastq" not in mod_read):
                        mod_read = mod_read + ".fa"
                    read_out = open(mod_read, "w")
                    tmp_datas.append(mod_read)
                    zips.append(mod_read)
                    print(" ".join(["Uncompressing", read]))
                    log.write(" ".join(["bzcat", read]) + "\n")
                    call(["bzcat", read], stdout=read_out)
                    log.write("\t" + mod_read + " is generated.\n")
                    read_out.close()
                elif read.endswith(".gz"):
                    mod_read = read.replace(".gz", "")
                    if (".fa" not in mod_read) and (
                            ".fasta" not in mod_read) and (
                            ".fna" not in mod_read) and (
                            ".fq" not in mod_read) and (
                            ".fastq" not in mod_read):
                        mod_read = mod_read + ".fa"
                    read_out = open(mod_read, "w")
                    tmp_datas.append(mod_read)
                    zips.append(mod_read)
                    print(" ".join(["Uncompressing", read]))
                    log.write(" ".join(["zcat", read]) + "\n")
                    call(["zcat", read], stdout=read_out)
                    read_out.close()
                    log.write("\t" + mod_read + " is generated.\n")
            tmp_reads.append({"sample": reads["sample"],
                              "files": tmp_datas, "zips": zips})   
        return tmp_reads

    def _run_segemehl_fasta_index(self, segemehl_path, fasta_path,
                                  index, fasta, log):
        log.write(" ".join([segemehl_path,
                  "-x", os.path.join(fasta_path, index),
                  "-d", os.path.join(fasta_path, fasta)]) + "\n")
        call([segemehl_path,
              "-x", os.path.join(fasta_path, index),
              "-d", os.path.join(fasta_path, fasta)])

    def _run_segemehl_align(self, args_circ, index, fasta, read,
                            sam_file, log_file, fasta_prefix, log):
        out = open(os.path.join(self.alignment_path,
                   fasta_prefix, sam_file), "w")
        log = open(os.path.join(self.alignment_path,
                   fasta_prefix, log_file), "w")
        log.write(" ".join([args_circ.segemehl_path,
                   "-i", os.path.join(self.fasta_path, index),
                   "-d", os.path.join(self.fasta_path, fasta),
                   "-q", read, "-S"]) + "\n")
        p = Popen([args_circ.segemehl_path,
                   "-i", os.path.join(self.fasta_path, index),
                   "-d", os.path.join(self.fasta_path, fasta),
                   "-q", read, "-S"],
                  stdout=out, stderr=log)
        return p

    def _align(self, args_circ, read_datas, log):
        '''align the read. if the bam files are provided, it can be skipped.'''
        prefixs = []
        align_files = []
        log.write("Using segemehl to align the read.\n")
        log.write("Please make sure the version of segemehl is at least 0.1.9.\n")
        for fasta in os.listdir(self.fasta_path):
            index = fasta.replace(".fa", ".idx")
            self._run_segemehl_fasta_index(args_circ.segemehl_path,
                                           self.fasta_path, index, fasta, log)
            processes = []
            num_process = 0
            fasta_prefix = fasta.replace(".fa", "")
            prefixs.append(fasta_prefix)
            self.helper.check_make_folder(os.path.join(
                                self.alignment_path, fasta_prefix))
            log.write("Running for {0}.\n".format(fasta_prefix))
            for reads in read_datas:
                for read in reads["files"]:
                    num_process += 1
                    read_name = read.split("/")[-1]
                    if read_name.endswith(".fa") or \
                       read_name.endswith(".fna") or \
                       read_name.endswith(".fasta") or \
                       read_name.endswith(".fq") or \
                       read_name.endswith(".fastq"):
                        filename = read_name.split(".")
                        read_prefix = ".".join(filename[:-1])
                        sam_file = "_".join([read_prefix, fasta_prefix + ".sam"])
                        log_file = "_".join([read_prefix, fasta_prefix + ".log"])
                        align_files.append("_".join([read_prefix, fasta_prefix]))
                        print("Mapping {0}".format(sam_file))
                        p = self._run_segemehl_align(
                                args_circ, index, fasta, read,
                                sam_file, log_file, fasta_prefix, log)
                        processes.append(p)
                        if num_process == args_circ.cores:
                            self._wait_process(processes)
                            num_process = 0
                self._wait_process(processes)
            log.write("Done!\n")
            log.write("The following files are generated in {0}:\n".format(
                  os.path.join(self.alignment_path, fasta_prefix)))
            for file_ in os.listdir(os.path.join(
                   self.alignment_path, fasta_prefix)):
                log.write("\t" + file_ + "\n")
        return align_files, prefixs

    def _run_samtools_convert_bam(self, samtools_path, pre_sam, out_bam, log):
        log.write(" ".join([samtools_path, "view",
                            "-bS", pre_sam, "-o", out_bam]) + "\n")
        call([samtools_path, "view", "-bS", pre_sam, "-o", out_bam])

    def _convert_sam2bam(self, sub_alignment_path, samtools_path, align_files, log):
        bam_files = []
        convert_ones = []
        remove_ones = []
        log.write("Using Samtools to convert SAM files to BAM files.\n")
        log.write("Please make sure the version of Samtools is at least 1.3.1.\n")
        for sam in os.listdir(sub_alignment_path):
            pre_sam = os.path.join(sub_alignment_path, sam)
            if sam.endswith(".sam"):
                bam_file = sam.replace(".sam", ".bam")
                print("Converting {0} to {1}".format(sam, bam_file))
                out_bam = os.path.join(sub_alignment_path, bam_file)
                self._run_samtools_convert_bam(samtools_path, pre_sam,
                                               out_bam, log)
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
        log.write("Done!\n")
        log.write("The following files are generated:\n")
        for file_ in os.listdir(sub_alignment_path):
            if file_.endswith(".bam"):
                log.write("\t" + os.path.join(sub_alignment_path, file_) + "\n")
        return bam_files, convert_ones, remove_ones

    def _run_samtools_merge_sort(self, samtools_path, prefix,
                                 out_folder, bam_datas, log):
        log.write("Using Samtools for merging, sorting and converting "
                  "the BAM files.\n")
        log.write("Make sure the version Samtools is at least 1.3.1.\n")
        for bam_data in bam_datas:
            print("Merging bam files for {0} of {1}".format(
                prefix, bam_data["sample"]))
            sample_bam = os.path.join(out_folder, "_".join([
                prefix, bam_data["sample"] + ".bam"]))
            if len(bam_data["files"]) <= 1:
                shutil.copyfile(bam_data["files"][0], sample_bam)
            else:
                file_line = " ".join(bam_data["files"])
                log.write(" ".join([samtools_path, "merge",
                                    sample_bam, file_line]) + "\n")
                os.system(" ".join([samtools_path, "merge",
                                    sample_bam, file_line]))
            print("Sorting bam files for {0} of {1}".format(
                prefix, bam_data["sample"]))
            sort_sample = os.path.join(out_folder,
                  "_".join([prefix, bam_data["sample"] + "_sort.bam"]))
            log.write(" ".join([samtools_path, "sort",
                      "-o", sort_sample, sample_bam]) + "\n")
            call([samtools_path, "sort", "-o", sort_sample, sample_bam])
            os.remove(sample_bam)
            print("Converting bam files to sam files for {0} of {1}".format(
                prefix, bam_data["sample"]))
            log.write(" ".join([samtools_path, "view", "-h", "-o",
                      sort_sample.replace(".bam", ".sam"), sort_sample]) + "\n")
            call([samtools_path, "view", "-h", "-o",
                  sort_sample.replace(".bam", ".sam"), sort_sample])
        log.write("Done!\n")
        log.write("\t" + sort_sample.replace(".bam", ".sam") + " is generated.\n")

    def _merge_sort_aligment_file(
            self, bam_datas, read_datas, samtools_path,
            out_folder, convert_ones, tmp_reads, remove_ones, prefix, log):
        if bam_datas is None:
            merge_bam_datas = []
            for read_data in read_datas:
                bam_files = []
                for read in read_data["files"]:
                    if read.endswith(".gz") or read.endswith(".bz2"):
                        read = ".".join(
                                read.split("/")[-1].split(".")[:-1])
                    read_prefix = ".".join(
                        read.split("/")[-1].split(".")[:-1])
                    bam_files.append(os.path.join(
                        self.alignment_path, prefix,
                        "_".join([read_prefix, prefix + ".bam"])))
                merge_bam_datas.append({"sample": read_data["sample"],
                                        "files": bam_files})
        elif (bam_datas is not None) and (read_datas is not None):
            merge_bam_datas = copy.deepcopy(bam_datas)
            for bam_data in merge_bam_datas:
                for read_data in read_datas:
                    if bam_data["sample"] == read_data["sample"]:
                        for read in read_data["files"]:
                            read_prefix = ".".join(
                                read.split("/")[-1].split(".")[:-1])
                            bam = os.path.join(
                                self.alignment_path, prefix,
                                "_".join([read_prefix, prefix + ".bam"]))
                            if (bam not in bam_data["files"]):
                                bam_data["files"].append(bam)
        else:
            merge_bam_datas = copy.deepcopy(bam_datas)
        self._run_samtools_merge_sort(samtools_path, prefix,
                                      out_folder, merge_bam_datas, log)
        for bam in convert_ones:
            os.remove(bam)
        for sam in remove_ones:
            os.remove(sam)

    def _run_testrealign(self, prefix, testrealign_path, out_folder, log):
        log.write("Using Segemehl to detect circular RNAs.\n")
        log.write("Please make sure the version of Segemehl is at least 0.1.9.\n")
        log.write("Please make sure your testrealign.x exists. If it does not "
                  "exists, please reinstall your Segemehl via using make all.\n")
        sub_splice_path = os.path.join(self.splice_path, prefix)
        if not os.path.exists(sub_splice_path):
            os.mkdir(sub_splice_path)
        err_log = os.path.join(sub_splice_path, prefix + ".log")
        print("Running testrealign.x for {0}".format(prefix))
        for sam_file in os.listdir(out_folder):
            if sam_file.endswith("sort.sam"):
                sample_prefix = sam_file.replace("_sort.sam", "")
                command = " ".join([
                    testrealign_path,
                    "-d", os.path.join(self.fasta_path, prefix + ".fa"),
                    "-q", os.path.join(out_folder, sam_file), "-n",
                    "-U", os.path.join(sub_splice_path,
                                       sample_prefix + "_splicesites.bed"),
                    "-T", os.path.join(sub_splice_path,
                                       sample_prefix + "_transrealigned.bed")])
                log.write(command + " 2>" + err_log + "\n")
                os.system(command + " 2>" + err_log)
        log.write("Done!\n")
        log.write("The following files are generated:\n")
        for file_ in os.listdir(sub_splice_path):
            log.write("\t" + os.path.join(sub_splice_path, file_) + "\n")
        self.helper.remove_all_content(out_folder, ".sam", "file")

    def _merge_bed(self, fastas, splice_path, output_folder):
        '''Merge the bed files for analysis'''
        fa_prefixs = []
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
                fa_prefixs.append(fasta_prefix)
                bed_folder = os.path.join(
                    output_folder, fasta_prefix)
                self.helper.check_make_folder(bed_folder)
                samples = []
                for header in headers:
                    for splice in os.listdir(os.path.join(
                            splice_path, header)):
                        if splice.endswith(".bed"):
                            if self.splices["file"] in splice:
                                sample = splice.replace(header, "")
                                sample = sample.replace(
                                    self.splices["file"], "")
                                if sample not in samples:
                                    samples.append(sample)
                            shutil.copyfile(
                                os.path.join(
                                splice_path, header, splice),
                                os.path.join(
                                bed_folder, "tmp_" + splice))
                for sample in samples:
                    out_splice = os.path.join(bed_folder, "".join([
                        fasta_prefix + sample + self.splices["file"]]))
                    out_trans = os.path.join(bed_folder, "".join([
                        fasta_prefix + sample + self.trans["file"]]))
                    if os.path.exists(out_splice):
                        os.remove(out_splice)
                    if os.path.exists(out_trans):
                        os.remove(out_trans)
                    for file_ in os.listdir(bed_folder):
                        if (self.splices["splice"] in file_) and (
                                sample in file_):
                            self.helper.merge_file(os.path.join(
                                    bed_folder, file_), out_splice)
                        elif (self.trans["trans"] in file_) and (
                                sample in file_):
                            self.helper.merge_file(os.path.join(
                                    bed_folder, file_), out_trans)
        self.helper.remove_all_content(splice_path, None, "dir")
        return samples, fa_prefixs

    def _stat_and_gen_gff(self, prefixs, samples, args_circ, log):
        '''do statistics and print the result to gff file'''
        log.write("Running circRNA.py to do statistics and generate gff files.\n")
        log.write("The following files are generated:\n")
        for prefix in prefixs:
            self.helper.check_make_folder(os.path.join(self.gff_folder,
                                                       prefix))
            self.helper.check_make_folder(os.path.join(self.splice_path,
                                                       prefix))
            for bed in os.listdir(os.path.join(
                args_circ.output_folder, prefix)):
                if (bed.split("_")[0] != "tmp") and (bed.endswith(".bed")):
                    shutil.copy(
                        os.path.join(args_circ.output_folder, prefix, bed),
                        os.path.join(self.splice_path, prefix))
            self.helper.check_make_folder(os.path.join(
                                          self.candidate_path, prefix))
            print("Comparing circular RNAs with annotations of {0}".format(
                prefix))
            for sample in samples:
                splice_file = os.path.join(
                    self.splice_path, prefix,
                    "".join([prefix, sample, self.splices["file"]]))
                stat_file = os.path.join(args_circ.stat_folder,
                               "".join(["stat_", prefix, sample,
                                        "circRNA.csv"]))
                csv_all = os.path.join(self.candidate_path, prefix,
                               "".join([prefix, sample, "circRNA_all.csv"]))
                csv_best = os.path.join(self.candidate_path, prefix,
                               "".join([prefix, sample, "circRNA_best.csv"]))
                gff_all = os.path.join(self.gff_folder, prefix,
                                "".join([prefix, sample, "circRNA_all.gff"]))
                gff_best = os.path.join(self.gff_folder, prefix,
                                  "".join([prefix, sample, "circRNA_best.gff"]))
                detect_circrna(splice_file, os.path.join(
                               self.gff_path, prefix + ".gff"), csv_all,
                               args_circ, stat_file)
                self.converter.convert_circ2gff(
                     os.path.join(self.candidate_path, prefix,
                                  "".join([prefix, sample, "circRNA_all.csv"])),
                     args_circ, gff_all, gff_best)
                log.write("\t" + stat_file + "\n")
                log.write("\t" + csv_all + "\n")
                log.write("\t" + csv_best + "\n")
                log.write("\t" + gff_all + "\n")
                log.write("\t" + gff_best + "\n")

    def _extract_input_files(self, inputs):
        input_datas = []
        for input_ in inputs:
            datas = input_.split(":")
            if len(datas) != 2:
                print("Error: the format of --bam_files or "
                      "--read_files is wrong!")
                sys.exit()
            for file_ in datas[-1].split(","):
                if not os.path.exists(file_):
                    print("Error: some files in --bam_files or "
                          "--read_files do not exist!")
                    sys.exit()
            input_datas.append({"sample": datas[0],
                                "files": datas[-1].split(",")})
        return input_datas

    def _combine_read_bam(self, bam_files, bam_datas, read_datas):
        if bam_datas is not None:
            for bam_data in bam_datas:
                for read_data in read_datas:
                    if bam_data["sample"] == read_data["sample"]:
                        for read in read_data["files"]:
                            prefix = ".".join(
                                read.split("/")[-1].split(".")[:-1])
                            bam = os.path.join(self.alignment_path,
                                               prefix + ".bam")
                            if (bam in bam_files) and (
                                    bam not in bam_data["files"]):
                                bam_data["files"].append(bam)
        else:
            bam_datas = []
            for read_data in read_datas:
                bam_files = []
                for read in read_data["files"]:
                    prefix = ".".join(
                        read.split("/")[-1].split(".")[:-1])
                    bam_files.append(os.path.join(
                        self.alignment_path, prefix + ".bam"))
                bam_datas.append({"sample": read_data["sample"],
                                  "files": bam_files})
        return bam_datas

    def _remove_tmp_files(self, args_circ, fa_prefixs):
        self.helper.remove_tmp_dir(args_circ.fastas)
        self.helper.remove_tmp_dir(args_circ.gffs)
        self.helper.remove_all_content(args_circ.output_folder,
                                       ".bam", "file")
        for prefix in fa_prefixs:
            shutil.rmtree(os.path.join(args_circ.output_folder, prefix))

    def run_circrna(self, args_circ, log):
        '''detection of circRNA'''
        bam_datas = None
        read_datas = None
        if (args_circ.bams is None) and (args_circ.read_files is None):
            log.write("--bam_files and --read_files can not be both emtpy.\n")
            print("Error: --bam_files or --read_files should be assigned.")
            sys.exit()
        if args_circ.bams is not None:
            bam_datas = self._extract_input_files(args_circ.bams)
        if args_circ.read_files is not None:
            read_datas = self._extract_input_files(args_circ.read_files)
        for gff in os.listdir(args_circ.gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(
                                                 args_circ.gffs, gff))
        if args_circ.segemehl_path is None:
            log.write("segemehl does not exists.\n")
            print("Error: please assign segemehl path!!")
            sys.exit()
        self.multiparser.parser_fasta(args_circ.fastas)
        self.multiparser.parser_gff(args_circ.gffs, None)
        self.multiparser.combine_gff(args_circ.fastas, self.gff_path,
                                     "fasta", None)
        tmp_reads = []
        if args_circ.read_files:
            log.write("Raw read files are found.\n")
            tmp_reads = self._deal_zip_file(read_datas, log)
            align_files, prefixs = self._align(args_circ, tmp_reads, log)
        else:
            align_files = None
        prefixs = []
        for fasta in os.listdir(self.fasta_path):
            if fasta.endswith(".fa"):
                fasta_prefix = fasta.replace(".fa", "")
                prefixs.append(fasta_prefix)
        for prefix in prefixs:
            if args_circ.read_files:
                sub_alignment_path = os.path.join(self.alignment_path, prefix)
                bam_files, convert_ones, remove_ones = self._convert_sam2bam(
                sub_alignment_path, args_circ.samtools_path, align_files, log)
            else:
                convert_ones = []
                remove_ones = []
            self._merge_sort_aligment_file(
                bam_datas, read_datas, args_circ.samtools_path,
                args_circ.output_folder,
                convert_ones, tmp_reads, remove_ones, prefix, log)
            self._run_testrealign(prefix, args_circ.testrealign_path,
                                  args_circ.output_folder, log)
        samples, fa_prefixs = self._merge_bed(
            args_circ.fastas, self.splice_path, args_circ.output_folder)
        self._stat_and_gen_gff(fa_prefixs, samples, args_circ, log)
        if len(tmp_reads) != 0:
            for reads in tmp_reads:
                for read in reads["zips"]:
                    os.remove(read)
        self._remove_tmp_files(args_circ, fa_prefixs)
