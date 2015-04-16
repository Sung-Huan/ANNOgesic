#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
from subprocess import call
from annogesiclib.multiparser import Multiparser
from annogesiclib.helper import Helper
from annogesiclib.sORF_intergenic import get_intergenic
from annogesiclib.sORF_detection import sorf_detection
from annogesiclib.stat_sorf import stat


class sORF_detection(object):

    def __init__(self):
        self.multiparser = Multiparser()
        self.helper = Helper()    

    def _check_gff(self, gffs):
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(gffs, gff))

    def _check_necessary_files(self, gffs, trans, tex_wigs, frag_wigs, utr_detect, tsss, sRNAs):
        if (gffs is None) or \
           (trans is None) or \
           ((tex_wigs is None) and (frag_wigs is None)):
            print("Error: lack required files!!!!")
            sys.exit()
        if utr_detect:
            if (tsss is None):
                print("Error: lack required files for UTR derived sRNA detection!!!!")
                sys.exit()
#        self._check_gff(gffs)
        self.multiparser._parser_gff(gffs, None)
        if tsss is not False:
            self._check_gff(tsss)
            tss_path = os.path.join(tsss, "tmp")
            self.multiparser._parser_gff(tsss, "TSS")
            self.multiparser._combine_gff(gffs, tss_path, None, "TSS")
        else:
            tss_path = None
#        self._check_gff(trans)
        if sRNAs is not False:
            self._check_gff(sRNAs)
            srna_path = os.path.join(sRNAs, "tmp")
            self.multiparser._parser_gff(sRNAs, "sRNA")
            self.multiparser._combine_gff(gffs, srna_path, None, "sRNA")
        else:
            srna_path = None
        return (tss_path, srna_path)

    def _merge_wigs(self, tex_wigs, frag_wigs):
        if (tex_wigs is not False) and (frag_wigs is not False):
            folder = tex_wigs.split("/")
            folder = "/".join(folder[:-1])
            merge_wigs = os.path.join(folder, "merge_wigs")
            self.helper.check_make_folder(folder, "merge_wigs")
            for wig in os.listdir(tex_wigs):
                if os.path.isdir(os.path.join(tex_wigs, wig)):
                    pass
                else:
                    shutil.copy(os.path.join(tex_wigs, wig), merge_wigs)
            for wig in os.listdir(frag_wigs):
                if os.path.isdir(os.path.join(frag_wigs, wig)):
                    pass
                else:
                    shutil.copy(os.path.join(frag_wigs, wig), merge_wigs)
        elif (tex_wigs is not False):
            merge_wigs = tex_wigs
        elif (frag_wigs is not False):
            merge_wigs = frag_wigs
        return merge_wigs

    def _combine_libs(self, tlibs, flibs):
        if (tlibs is False) and (flibs is False):
            print("Error: please input proper libraries!!")
        if (tlibs is not False) and (flibs is not False):
            libs = tlibs + flibs
        elif (tlibs is not False):
            libs = tlibs
        elif (flibs is not False):
            libs = flibs
        return libs

    def _start_stop_codon(self, prefixs, fasta_path, out_folder, utr_length, libs,
                          tex_notex, replicate, cutoff_inter, cutoff_3utr, cutoff_5utr,
                          cutoff_interCDS, wig_path, merge_wigs, start_codon, stop_codon,
                          max_len, min_len, condition, gff_output, tss_path, srna_path,
                          table_best, utr_detect):
        for prefix in prefixs:
            if srna_path is not None:
                srna_file = os.path.join(srna_path, "_".join([prefix, "sRNA.gff"]))
            if tss_path is not None:
                tss_file = os.path.join(tss_path, "_".join([prefix, "TSS.gff"]))
            sorf_detection(os.path.join(fasta_path, prefix + ".fa"), srna_file,
                           os.path.join(out_folder, "_".join([prefix, "inter.gff"])), 
                           tss_file, utr_length, utr_detect, input_libs, tex_notex, 
                           replicate, cutoff_inter, cutoff_3utr, cutoff_5utr, cutoff_interCDS, 
                           os.path.join(wig_path, "_".join([prefix, "forward.wig"])), 
                           os.path.join(wig_path, "_".join([prefix, "reverse.wig"])), 
                           merge_wigs, start_codon, stop_codon, table_best, max_len, min_len, condition, 
                           os.path.join(gff_output, "all_candidates", "_".join([prefix, "sORF"])))
            if "_".join([prefix, "sORF_all.gff"]) in os.listdir(os.path.join(gff_output, "all_candidates")):
                os.rename(os.path.join(gff_output, "all_candidates", "_".join([prefix, "sORF_all.gff"])), \
                          os.path.join(gff_output, "all_candidates", "_".join([prefix, "sORF.gff"])))
                os.rename(os.path.join(gff_output, "all_candidates", "_".join([prefix, "sORF_best.gff"])), \
                          os.path.join(gff_output, "best", "_".join([prefix, "sORF.gff"])))
                os.rename(os.path.join(gff_output, "all_candidates", "_".join([prefix, "sORF_all.csv"])), \
                          os.path.join(table_output, "all_candidates", "_".join([prefix, "sORF.csv"])))
                os.rename(os.path.join(gff_output, "all_candidates", "_".join([prefix, "sORF_best.csv"])), \
                          os.path.join(table_output, "best", "_".join([prefix, "sORF.csv"])))

    def run_sorf_detection(self, bin_path, out_folder, utr_detect, trans, gffs, tsss, utr_length, 
                           min_len, max_len, tex_wigs, frag_wigs, cutoff_inter, cutoff_5utr, 
                           cutoff_3utr, cutoff_interCDS, fastas, tlibs, flibs, tex_notex, 
                           replicate, table_best, sRNAs, start_coden, stop_coden, condition):
        paths = self._check_necessary_files(gffs, trans, tex_wigs, frag_wigs, utr_detect, tsss, sRNAs)
        tss_path = paths[0]
        srna_path = paths[1]
        gff_output = os.path.join(out_folder, "gffs")
        table_output = os.path.join(out_folder, "tables")
        merge_wigs = self._merge_wigs(tex_wigs, frag_wigs)
        wig_path = os.path.join(merge_wigs, "tmp")
        self.multiparser._parser_wig(merge_wigs)
        self.multiparser._combine_wig(gffs, wig_path, None)
        tran_path = os.path.join(trans, "tmp")
        self.multiparser._parser_gff(trans, "transcript")
        self.multiparser._combine_gff(gffs, tran_path, None, "transcript")
        fasta_path = os.path.join(fastas, "tmp")
        self.multiparser._parser_fasta(fastas)
        self.multiparser._combine_fasta(gffs, fasta_path, None)
        libs = self._combine_libs(tlibs, flibs)
        prefixs = []
        for gff in os.listdir(gff_path): ## compare transcript and CDS
            if gff.endswith(".gff"):
                prefix = gff.replace(".gff", "")
                prefixs.append(prefix)
                print("comparing transcript and CDS of " + prefix)
                get_intergenic(os.path.join(gff_path, gff), 
                               os.path.join(tran_path, "_".join([prefix, "transcript.gff"])), 
                               os.path.join(out_folder, "_".join([prefix, "inter.gff"])), 
                               utr_detect)
        self._start_stop_codon(prefixs, fasta_path, out_folder, utr_length, libs, tex_notex, 
                               replicate, cutoff_inter, cutoff_3utr, cutoff_5utr, cutoff_interCDS, 
                               wig_path, merge_wigs, start_codon, stop_codon, max_len, min_len, 
                               condition, gff_output, tss_path, srna_path, table_best, utr_detect)
        for sorf in os.listdir(gff_output + "all_candidates"): # stat
            print("statistics of {0}".format(sorf))
            if sorf.endswith("_sORF.gff"):
                stat(os.path.join(gff_output, "all_candidates", sorf),
                     os.path.join(out_folder, "statistics", 
                                  "_".join(["stat", sorf.replace(".gff", ".csv")])), 
                     utr_detect)
        # remove temperary files #
        os.system("rm " + out_folder + "/*.gff")
        self.helper.remove_tmp(fastas)
        self.helper.remove_tmp(gffs)
        self.helper.remove_tmp(tsss)
        self.helper.remove_tmp(trans)
        self.helper.remove_tmp(sRNAs)
        self.helper.remove_wig(tex_wigs)
        self.helper.remove_wig(frag_wigs)
