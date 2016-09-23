import os
import sys
import shutil
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.sORF_intergenic import get_intergenic
from annogesiclib.sORF_detection import sorf_detection
from annogesiclib.stat_sorf import stat


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
        self.best = "best"

    def _check_gff(self, gffs):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(gffs, gff))

    def _check_necessary_files(self, args_sorf):
        if (args_sorf.gffs is None) or \
           (args_sorf.trans is None) or \
           ((args_sorf.tex_wigs is None) and (args_sorf.frag_wigs is None)):
            print("Error: lack required files!!!!")
            sys.exit()
        if args_sorf.utr_detect:
            if (args_sorf.tsss is None):
                print("Error: lack required files for UTR derived"
                      " sORF detection!!!!")
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

    def _start_stop_codon(self, prefixs, args_sorf):
        '''detect the sORF based on start and stop codon 
        and ribosome binding site'''
        for prefix in prefixs:
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
                shutil.move(os.path.join(self.gff_output, self.all_cand,
                            "_".join([prefix, "sORF_all.gff"])),
                            os.path.join(self.gff_output, self.all_cand,
                            "_".join([prefix, "sORF.gff"])))
                shutil.move(os.path.join(self.gff_output, self.all_cand,
                            "_".join([prefix, "sORF_best.gff"])),
                            os.path.join(self.gff_output, self.best,
                            "_".join([prefix, "sORF.gff"])))
                shutil.move(os.path.join(self.gff_output, self.all_cand,
                            "_".join([prefix, "sORF_all.csv"])),
                            os.path.join(self.table_output, self.all_cand,
                            "_".join([prefix, "sORF.csv"])))
                shutil.move(os.path.join(self.gff_output, self.all_cand,
                            "_".join([prefix, "sORF_best.csv"])),
                            os.path.join(self.table_output, self.best,
                            "_".join([prefix, "sORF.csv"])))

    def _remove_tmp(self, args_sorf):
        self.helper.remove_all_content(args_sorf.out_folder, ".gff", "file")
        self.helper.remove_tmp(args_sorf.fastas)
        self.helper.remove_tmp(args_sorf.gffs)
        self.helper.remove_tmp(args_sorf.tsss)
        self.helper.remove_tmp(args_sorf.trans)
        self.helper.remove_tmp(args_sorf.srnas)
        self.helper.remove_wigs(args_sorf.tex_wigs)
        self.helper.remove_wigs(args_sorf.frag_wigs)

    def _compare_tran_cds(self, args_sorf):
        '''compare transcript and CDS to find the intergenic region'''
        prefixs = []
        for gff in os.listdir(args_sorf.gffs):
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                prefixs.append(prefix)
                print("comparing transcript and CDS of {0}".format(prefix))
                get_intergenic(os.path.join(args_sorf.gffs, gff),
                               os.path.join(self.tran_path,
                               "_".join([prefix, "transcript.gff"])),
                               os.path.join(args_sorf.out_folder,
                               "_".join([prefix, "inter.gff"])),
                               args_sorf.utr_detect, args_sorf.hypo)
        return prefixs

    def run_sorf_detection(self, args_sorf):
        if args_sorf.fuzzy_rbs > 6:
            print("Error: --fuzzy_rbs should be equal or less than 6!!")
            sys.exit()
        self._check_necessary_files(args_sorf)
        self.multiparser.parser_gff(args_sorf.trans, "transcript")
        self.multiparser.combine_gff(args_sorf.gffs, self.tran_path,
                                     None, "transcript")
        self.multiparser.parser_fasta(args_sorf.fastas)
        self.multiparser.combine_fasta(args_sorf.gffs, self.fasta_path, None)
        prefixs = self._compare_tran_cds(args_sorf)
        self._start_stop_codon(prefixs, args_sorf)
        for sorf in os.listdir(os.path.join(self.gff_output, self.all_cand)):
            print("statistics of {0}".format(sorf))
            if sorf.endswith("_sORF.gff"):
                stat(os.path.join(self.gff_output, self.all_cand, sorf),
                     os.path.join(self.gff_output, self.best, sorf),
                     os.path.join(args_sorf.out_folder, "statistics",
                     "_".join(["stat", sorf.replace(".gff", ".csv")])),
                     args_sorf.utr_detect)
        self._remove_tmp(args_sorf)
