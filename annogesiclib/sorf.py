import os
import sys
import shutil
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.sORF_intergenic import get_intergenic
from annogesiclib.sORF_detection import sorf_detection
from annogesiclib.stat_sorf import stat
from annogesiclib.reorganize_table import reorganize_table


class sORFDetection(object):
    '''detection of sORF'''

    def __init__(self, args_sorf):
        self.multiparser = Multiparser()
        self.helper = Helper()
        if args_sorf.tsss is not None:
            self.tss_path = os.path.join(args_sorf.tsss, "tmp")
        else:
            self.tss_path = None
        if args_sorf.srnas is not None:
            self.srna_path = os.path.join(args_sorf.srnas, "tmp")
        else:
            self.srna_path = None
        self.gff_output = os.path.join(args_sorf.out_folder, "gffs")
        self.table_output = os.path.join(args_sorf.out_folder, "tables")
        self.tran_path = os.path.join(args_sorf.trans, "tmp")
        self.fasta_path = os.path.join(args_sorf.fastas, "tmp")
        self.all_cand = "all_candidates"
        self.best = "best_candidates"

    def _check_gff(self, gffs):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(gffs, gff))

    def _check_necessary_files(self, args_sorf, log):
        if (args_sorf.gffs is None) or (args_sorf.trans is None) or (
               (args_sorf.tex_wigs is None) and (args_sorf.frag_wigs is None)):
            print("Error: lack required files!")
            log.write("genome annotation, transcript file or wiggle files "
                      "are not assigned.\n")
            sys.exit()
        if args_sorf.utr_detect:
            if (args_sorf.tsss is None):
                print("Error: TSS files are required for UTR derived"
                      " sORF detection!")
                log.write("TSS files are required for UTR derived"
                          " sORF detection!\n")
                sys.exit()
        self._check_gff(args_sorf.gffs)
        self.multiparser.parser_gff(args_sorf.gffs, None)
        if args_sorf.tsss is not None:
            self._check_gff(args_sorf.tsss)
            self.multiparser.parser_gff(args_sorf.tsss, "TSS")
            self.multiparser.combine_gff(args_sorf.gffs, self.tss_path,
                                         None, "TSS")
        self._check_gff(args_sorf.trans)
        if args_sorf.srnas is not None:
            self._check_gff(args_sorf.srnas)
            self.multiparser.parser_gff(args_sorf.srnas, "sRNA")
            self.multiparser.combine_gff(args_sorf.gffs, self.srna_path,
                                         None, "sRNA")

    def _start_stop_codon(self, prefixs, args_sorf, log):
        '''detect the sORF based on start and stop codon 
        and ribosome binding site'''
        log.write("Running sORF_detection.py for detecting sORFs.\n")
        log.write("The following files are generated:\n")
        for prefix in prefixs:
            print("Searching sORFs of {0}".format(prefix))
            if self.srna_path is not None:
                srna_file = os.path.join(self.srna_path,
                                         "_".join([prefix, "sRNA.gff"]))
            else:
                srna_file = None
            if self.tss_path is not None:
                tss_file = os.path.join(self.tss_path,
                                        "_".join([prefix, "TSS.gff"]))
            else:
                tss_file = None
            sorf_detection(os.path.join(self.fasta_path, prefix + ".fa"),
                           srna_file, os.path.join(args_sorf.out_folder,
                           "_".join([prefix, "inter.gff"])), tss_file,
                           os.path.join(args_sorf.wig_path,
                           "_".join([prefix, "forward.wig"])),
                           os.path.join(args_sorf.wig_path,
                           "_".join([prefix, "reverse.wig"])),
                           os.path.join(self.gff_output, self.all_cand,
                           "_".join([prefix, "sORF"])), args_sorf)
            if "_".join([prefix, "sORF_all.gff"]) in os.listdir(
                         os.path.join(self.gff_output, self.all_cand)):
                gff_all = os.path.join(self.gff_output, self.all_cand,
                                       "_".join([prefix, "sORF.gff"]))
                gff_best = os.path.join(self.gff_output, self.best,
                                        "_".join([prefix, "sORF.gff"]))
                csv_all = os.path.join(self.table_output, self.all_cand,
                                       "_".join([prefix, "sORF.csv"]))
                csv_best =  os.path.join(self.table_output, self.best,
                                         "_".join([prefix, "sORF.csv"]))
                shutil.move(os.path.join(self.gff_output, self.all_cand,
                            "_".join([prefix, "sORF_all.gff"])), gff_all)
                shutil.move(os.path.join(self.gff_output, self.all_cand,
                            "_".join([prefix, "sORF_best.gff"])), gff_best)
                shutil.move(os.path.join(self.gff_output, self.all_cand,
                            "_".join([prefix, "sORF_all.csv"])), csv_all)
                shutil.move(os.path.join(self.gff_output, self.all_cand,
                            "_".join([prefix, "sORF_best.csv"])), csv_best)
                log.write("\t" + gff_all + "\n")
                log.write("\t" + gff_best + "\n")
                log.write("\t" + csv_all + "\n")
                log.write("\t" + csv_best + "\n")

    def _remove_tmp(self, args_sorf):
        self.helper.remove_all_content(args_sorf.out_folder, ".gff", "file")
        self.helper.remove_tmp_dir(args_sorf.fastas)
        self.helper.remove_tmp_dir(args_sorf.gffs)
        self.helper.remove_tmp_dir(args_sorf.tsss)
        self.helper.remove_tmp_dir(args_sorf.trans)
        self.helper.remove_tmp_dir(args_sorf.srnas)
        if "temp_wig" in os.listdir(args_sorf.out_folder):
            shutil.rmtree(os.path.join(args_sorf.out_folder, "temp_wig"))
        if "merge_wigs" in os.listdir(args_sorf.out_folder):
            shutil.rmtree(os.path.join(args_sorf.out_folder, "merge_wigs"))

    def _compare_tran_cds(self, args_sorf, log):
        '''compare transcript and CDS to find the intergenic region'''
        prefixs = []
        log.write("Running sORF_intergenic.py to extract the sequences of "
                  "potential sORFs\n")
        for gff in os.listdir(args_sorf.gffs):
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                prefixs.append(prefix)
                print("Comparing transcripts and CDSs of {0}".format(prefix))
                get_intergenic(os.path.join(args_sorf.gffs, gff),
                               os.path.join(self.tran_path,
                               "_".join([prefix, "transcript.gff"])),
                               os.path.join(args_sorf.out_folder,
                               "_".join([prefix, "inter.gff"])),
                               args_sorf.utr_detect, args_sorf.hypo,
                               args_sorf.extend_5, args_sorf.extend_3)
                log.write("\t" + os.path.join(args_sorf.out_folder,
                          "_".join([prefix, "inter.gff"])) + 
                          " is generated to temporary store the sequences.\n")
        return prefixs

    def _re_table(self, args_sorf, prefixs, log):
        log.write("Running re_table.py for generating coverage information.\n")
        log.write("The following files are updated:\n")
        for type_ in ["all_candidates", "best_candidates"]:
            for prefix in prefixs:
                table_file = os.path.join(args_sorf.out_folder, "tables",
                                          type_, "_".join([
                                          prefix, "sORF.csv"]))
                reorganize_table(args_sorf.libs, args_sorf.merge_wigs,
                                 "Track_detail", table_file)
                log.write("\t" + table_file + "\n")

    def run_sorf_detection(self, args_sorf, log):
        if args_sorf.fuzzy_rbs > 6:
            log.write("--fuzzy_rbs should be equal or less than 6!\n")
            print("Error: --fuzzy_rbs should be equal or less than 6!")
            sys.exit()
        self._check_necessary_files(args_sorf, log)
        self.multiparser.parser_gff(args_sorf.trans, "transcript")
        self.multiparser.combine_gff(args_sorf.gffs, self.tran_path,
                                     None, "transcript")
        self.multiparser.parser_fasta(args_sorf.fastas)
        self.multiparser.combine_fasta(args_sorf.gffs, self.fasta_path, None)
        prefixs = self._compare_tran_cds(args_sorf, log)
        self._start_stop_codon(prefixs, args_sorf, log)
        log.write("Running stat_sorf.py to do statistics.\n")
        for sorf in os.listdir(os.path.join(self.gff_output, self.all_cand)):
            print("Running statistics of {0}".format(sorf))
            if sorf.endswith("_sORF.gff"):
                stat_file = os.path.join(args_sorf.out_folder, "statistics",
                            "_".join(["stat", sorf.replace(".gff", ".csv")]))
                stat(os.path.join(self.gff_output, self.all_cand, sorf),
                     os.path.join(self.gff_output, self.best, sorf), stat_file,
                     args_sorf.utr_detect)
                log.write("\t" + stat_file + " is generated.\n")
        self._re_table(args_sorf, prefixs, log)
        self._remove_tmp(args_sorf)
