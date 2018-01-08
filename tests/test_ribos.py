import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
import annogesiclib.meme as me
from mock_helper import gen_file, import_data
import annogesiclib.ribos as rb
from mock_args_container import MockClass
from annogesiclib.ribos import Ribos


class Mock_func(object):

    def __init__(self):
        self.example = Example()

    def mock_run_cmscan(self, e_value,
                        seq, type_, prefix, test1, test2, test3):
        suffixs = {"csv": "test.csv",
                   "txt": "test_prescan.txt",
                   "re_txt": "test_scan.txt",
                   "re_csv": "test_scan.csv"}
        if type_ == "txt":
            gen_file("test_folder/output/test_RBS.txt", self.example.scan_file)
            return "test_folder/output/test_RBS.txt"
        else:
            gen_file("test_folder/output/test_RBS_rescan.txt",
                     self.example.rescan_file)
            return "test_folder/output/test_RBS_rescan.txt"

    def mock_modify_table(self, first_table, output_all):
        pass

    def mock_regenerate_seq(first_scan_file, first_seq,
                            first_table, sec_seq, test):
        if not os.path.exists('test_folder/output/tmp_table/test_test.csv'):
            gen_file('test_folder/output/tmp_table/test_test.csv', "test")

    def mock_reextract_rbs(sec_scan_file, first_table, sec_table, test, cutoff):
        gen_file('test_folder/output/tmp_table/test_test_scan.csv', "test")
        gen_file(os.path.join("test_folder/output", "tmp_fasta",
                              "test_regenerate.fa"), "test")
        pass

    def mock_stat_and_covert2gff(
            csv, gff, fuzzy, out_stat, feature, test2, test3):
        pass


class TestRibos(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.mock_args = MockClass()
        self.mock = Mock_func()
        self.test_folder = "test_folder"
        self.gffs = os.path.join(self.test_folder, "gffs")
        self.fastas = os.path.join(self.test_folder, "fastas")
        self.out_folder = os.path.join(self.test_folder, "output")
        self.database = os.path.join(self.test_folder, "database")
        self.seq_path = os.path.join(self.test_folder, "seqs")
        self.tables = os.path.join(self.out_folder, "tables")
        self.stat = os.path.join(self.out_folder, "statistics")
        self.scan = os.path.join(self.test_folder, "scan")
        self.tsss = os.path.join(self.test_folder, "tsss")
        self.trans = os.path.join(self.test_folder, "trans")
        self.out_gff = os.path.join(self.out_folder, "gffs")
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
            os.mkdir(self.tsss)
            os.mkdir(os.path.join(self.tsss, "tmp"))
            os.mkdir(self.trans)
            os.mkdir(os.path.join(self.trans, "tmp"))
            os.mkdir(self.gffs)
            os.mkdir(os.path.join(self.gffs, "tmp"))
            os.mkdir(self.fastas)
            os.mkdir(os.path.join(self.fastas, "tmp"))
            os.mkdir(self.out_folder)
            os.mkdir(self.database)
            os.mkdir(self.seq_path)
            os.mkdir(os.path.join(self.out_folder, "tmp_table"))
            os.mkdir(os.path.join(self.out_folder, "tmp_scan"))
            os.mkdir(os.path.join(self.out_folder, "tmp_fasta"))
            os.mkdir(os.path.join(self.out_folder, "scan_Rfam"))
            os.mkdir(self.tables)
            os.mkdir(self.scan)
            os.mkdir(self.stat)
            os.mkdir(self.out_gff)
        args = self.mock_args.mock()
        args.gffs = self.gffs
        args.fastas = self.fastas
        args.ribos_out_folder = self.out_folder
        args.database = self.database
        args.tsss = self.tsss
        args.trans = self.trans
        args.program = 'riboswtich'
        self.ribo = Ribos(args)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_scan_extract_rfam(self):
        self.ribo._run_cmscan = self.mock.mock_run_cmscan
        rb.modify_table = self.mock.mock_modify_table
        rb.regenerate_seq = self.mock.mock_regenerate_seq
        rb.reextract_rbs = self.mock.mock_reextract_rbs
        prefixs = []
        gen_file(os.path.join(self.gffs, "tmp/test.gff"),
                 self.example.gff_file)
        gen_file(os.path.join(self.fastas, "tmp/test.fa"),
                 self.example.fasta_file)
        gen_file(os.path.join(self.seq_path, "test.fa"),
                 self.example.fasta_file)
        gen_file(os.path.join(self.tsss, "tmp/test_TSS.gff"),
                 self.example.tss_file)
        gen_file(os.path.join(self.trans, "tmp/test_transcript.gff"),
                 self.example.tran_file)
        gen_file(os.path.join(self.out_folder, "tmp_fasta", "test.fa"),
                 self.example.fasta_file)
        args = self.mock_args.mock()
        args.start_codons = ["ATG"]
        args.fastas = self.fastas
        args.out_folder = self.out_folder
        args.gffs = self.gffs
        args.fuzzy = 5
        args.fuzzy_rbs = 2
        args.utr = True
        args.output_all = "test"
        args.cutoff = "e_0.01"
        tmp_files = {"fasta": os.path.join(self.out_folder, "tmp_fasta"),
                     "scan": "tmp_scan",
                     "table": os.path.join(self.out_folder, "tmp_table")}
        rfam = "Rfam_.cm"
        suffixs = {"csv": "test.csv",
                   "txt": "test_prescan.txt",
                   "re_txt": "test_scan.txt",
                   "re_csv": "test_scan.csv"}
        self.ribo._scan_extract_rfam(prefixs, args, tmp_files,
                                     suffixs, "test", rfam)
        self.assertListEqual(prefixs, ["test"])
        self.assertTrue(os.path.exists(os.path.join(
            self.out_folder, "tmp_fasta", "test_regenerate.fa")))

    def test_merge_results(self):
        rb.stat_and_covert2gff = self.mock.mock_stat_and_covert2gff
        gen_file(os.path.join(self.gffs, "test.gff"), self.example.gff_file) 
        gen_file(os.path.join(
            self.out_folder, "tmp_table/test_riboswitch.csv"),
                 self.example.table)
        gen_file(os.path.join(
            self.out_folder, "tmp_scan/test_riboswitch_prescan.txt"),
                 self.example.rescan_file)
        gen_file(os.path.join(
            self.out_folder, "tmp_scan/test_riboswitch_scan.txt"),
                 self.example.rescan_file)
        gen_file(os.path.join(
            self.test_folder, "ids"), self.example.ids)
        gen_file(os.path.join(
            self.tables, "test_riboswitch.csv"), self.example.table)
        gen_file('test_folder/output/tmp_table/test_test_scan.csv', "test")
        gen_file(os.path.join("test_folder/output", "tmp_fasta",
                              "test_regenerate.fa"), "test")
        gen_file('test_folder/output/tmp_scan/test_test_prescan.txt', "test")
        gen_file('test_folder/output/tmp_scan/test_test_scan.txt', "test")
        if not os.path.exists('test_folder/output/tmp_table/test_test.csv'):
            gen_file('test_folder/output/tmp_table/test_test.csv', "test")
        args = self.mock_args.mock()
        args.start_codons = ["ATG"]
        args.fastas = self.fastas
        args.out_folder = self.out_folder
        args.gffs = self.gffs
        args.ribos_id = os.path.join(self.test_folder, "ids")
        args.fuzzy = 3
        suffixs = {"csv": "test.csv",
                   "txt": "test_prescan.txt",
                   "re_txt": "test_scan.txt",
                   "re_csv": "test_scan.csv"}
        tmp_files = {"fasta": os.path.join(self.out_folder, "tmp_fasta"),
                     "scan": os.path.join(self.out_folder, "tmp_scan"),
                     "table": os.path.join(self.out_folder, "tmp_table")}
        rfam = "Rfam_.cm"
        self.ribo._merge_results(
            args, os.path.join(self.out_folder, "tmp_scan"), suffixs,
            tmp_files, os.path.join(self.out_folder, "tmp_scan"),
            os.path.join(self.out_folder, "scan_Rfam"),
            os.path.join(self.out_folder, "scan_Rfam"),
            os.path.join(self.out_folder, "gffs"), "riboswitch")

class Example(object):

    ids = """RF00162	SAM	SAM	riboswitch box leader
RF00174	Cobalamin	Cobalamin	riboswitch
RF00634	SAM-IV	S adenosyl methionine SAM	riboswitch"""
    table = """riboswitch_5\ttest\t+\tSAOUHSC_00013\t15948\t16046	RF00162	1.6e-18	1	99
riboswitch_11\ttest\t-\tSAOUHSC_00007\t27955\t28053	RF00162	1.6e-18	1	99"""
    fasta_file = """>test
AAAAGATAGCCGCTCGCTAAATAGGAGATCCCTCGAAGAATGATATGTGTG"""
    gff_file = """##gff-version 3
test	Refseq	gene	10	16	.	+	.	ID=gene0;Name=gene_0;locus_tag=AAA_00001
test	Refseq	CDS	10	16	.	+	.	ID=cds0;Name=CDS_0;locus_tag=AAA_00001
test	Refseq	gene	30	40	.	-	.	ID=gene1;Name=gene_1;locus_tag=AAA_00002
test	Refseq	CDS	30	40	.	-	.	ID=cds1;Name=CDS_1;locus_tag=AAA_00002"""
    tss_file = """##gff-version 3
test	Refseq	TSS	10	10	.	+	.	ID=tss0;Name=TSS_0
test	Refseq	TSS	30	30	.	-	.	ID=tss1;Name=TSS_1"""
    tran_file = """##gff-version 3
test	Refseq	Transcript	1	150	.	+	.	ID=ta0;Name=TA_0
test	Refseq	Transcript	30	230	.	-	.	ID=ta1;Name=TA_1"""
    rescan_file = """# cmscan :: search sequence(s) against a CM database
# INFERNAL 1.1.1 (July 2014)
# Copyright (C) 2014 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query sequence file:                   ANNOgesic/output/riboswitch/tmp_fasta/Staphylococcus_aureus_HG003.fa
# target CM database:                    ANNOgesic/input/database/Rfam_riboswitch.cm
# prefer accessions over names:          yes
# sequence inclusion threshold:          E-value <= 0.001
# number of worker threads:              24
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       riboswitch_0|Staphylococcus_aureus_HG003|+|SAOUHSC_00002|1878|2166  [L=289]
Hit scores:
 rank     E-value  score  bias  modelname  start    end   mdl trunc   gc  description
 ----   --------- ------ -----  --------- ------ ------   --- ----- ----  -----------

   [No hits detected that satisfy reporting thresholds]


Hit alignments:

   [No hits detected that satisfy reporting thresholds]


Internal CM pipeline statistics summary:
----------------------------------------
Query sequence(s):                                               1  (578 residues searched)
Query sequences re-searched for truncated hits:                  1  (944.5 residues re-searched, avg per model)
Target model(s):                                                36  (4015 consensus positions)
Windows   passing  local HMM SSV           filter:              34  (0.08743); expected (0.35)
Windows   passing  local HMM Viterbi       filter:                  (off)
Windows   passing  local HMM Viterbi  bias filter:                  (off)
Windows   passing  local HMM Forward       filter:              11  (0.02469); expected (0.02)
Windows   passing  local HMM Forward  bias filter:               9  (0.02051); expected (0.02)
Windows   passing glocal HMM Forward       filter:               4  (0.01173); expected (0.02)
Windows   passing glocal HMM Forward  bias filter:               3  (0.006459); expected (0.02)
Envelopes passing glocal HMM envelope defn filter:               3  (0.004762); expected (0.02)
Envelopes passing  local CM  CYK           filter:               0  (0); expected (0.0001)
Total CM hits reported:                                          0  (0); includes 0 truncated hit(s)

# CPU time: 0.13u 0.02s 00:00:00.15 Elapsed: 00:00:00.14
//
Query:       riboswitch_5|Staphylococcus_aureus_HG003|+|SAOUHSC_00013|15743|16116  [L=374]
Hit scores:
 rank     E-value  score  bias  modelname  start    end   mdl trunc   gc  description
 ----   --------- ------ -----  --------- ------ ------   --- ----- ----  -----------
  (1) !   6.2e-18   74.3   0.0  RF00162      206    304 +  cm    no 0.42  -
 ------ inclusion threshold ------
  (2) ?      0.24   11.1   0.0  RF01689      278    205 -  cm    no 0.43  -


Hit alignments:
>> RF00162
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
  (1) !   6.2e-18   74.3   0.0  cm        1      108 []         206         304 + .. 0.93    no 0.42

                                                                            v                                        NC
                                                                           ((((((((,,,,,,<<<---<<<_.____>>>------>>> CS
                                                               RF00162   1 cucUuAUcaAGAGaGGcgGAGGGA.cuGGCCCuaUGAAgCC 40
                                                                           : CUUAU   GAGA:G: GAGGGA +UGGCCCU UGA :C:
  riboswitch_5|Staphylococcus_aureus_HG003|+|SAOUHSC_00013|15743|16116 206 UUCUUAUUGUGAGAAGUUGAGGGAcUUGGCCCUGUGAUACU 246
                                                                           ***************************************** PP

                                                                                     vv         vv                   NC
                                                                           ,,<<<-<<<<<<_________>>>-->>>>>>,,,,,,,,< CS
                                                               RF00162  41 uCgGCAACCcccuauaauauagggAaGGUGCcAAuUCCugC 81
                                                                           UC:GCAACC:  U+UA    +  :A GGUGC:AA+ CC+ C
  riboswitch_5|Staphylococcus_aureus_HG003|+|SAOUHSC_00013|15743|16116 247 UCAGCAACCGACUUUA----UAGCACGGUGCUAAAACCAAC 283
                                                                           *********5444443....367899*************** PP

                                                                                                    v  NC
                                                                           <<<_________>>>>,,,)))))))) CS
                                                               RF00162  82 cggauauuuuuuccgGaaAgAUaAgag 108
                                                                           ::G+      U+C::GAA  AUAAG :
  riboswitch_5|Staphylococcus_aureus_HG003|+|SAOUHSC_00013|15743|16116 284 GAGU------UACUCGAAUGAUAAGUA 304
                                                                           ***8......56*************** PP
Internal CM pipeline statistics summary:
----------------------------------------
Query sequence(s):                                               1  (748 residues searched)
Query sequences re-searched for truncated hits:                  1  (943.0 residues re-searched, avg per model)
Target model(s):                                                36  (4015 consensus positions)
Windows   passing  local HMM SSV           filter:              56  (0.1655); expected (0.35)
Windows   passing  local HMM Viterbi       filter:                  (off)
Windows   passing  local HMM Viterbi  bias filter:                  (off)
Windows   passing  local HMM Forward       filter:              13  (0.04085); expected (0.02)
Windows   passing  local HMM Forward  bias filter:              12  (0.03619); expected (0.02)
Windows   passing glocal HMM Forward       filter:               6  (0.02035); expected (0.02)
Windows   passing glocal HMM Forward  bias filter:               6  (0.02035); expected (0.02)
Envelopes passing glocal HMM envelope defn filter:               6  (0.01398); expected (0.02)
Envelopes passing  local CM  CYK           filter:               5  (0.01176); expected (0.0001)
Total CM hits reported:                                          2  (0.002842); includes 0 truncated hit(s)

# CPU time: 1.55u 0.02s 00:00:01.57 Elapsed: 00:00:01.41
//"""

    scan_file = """Query:       riboswitch_11|Staphylococcus_aureus_HG003|-|SAOUHSC_00007|10446|33555  [L=23110]
Hit scores:
 rank     E-value  score  bias  modelname  start    end   mdl trunc   gc  description
 ----   --------- ------ -----  --------- ------ ------   --- ----- ----  -----------
  (1) !   3.8e-16   74.3   0.0  RF00162    17608  17510 -  cm    no 0.42  -
 ------ inclusion threshold ------
  (2) ?       5.7   14.5   0.0  RF01826    13085  13058 -  cm    no 0.25  -
  (3) ?       6.3   10.1   0.0  RF00521     7480   7418 -  cm    no 0.33  -


Hit alignments:
>> RF00162
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
  (1) !   3.8e-16   74.3   0.0  cm        1      108 []       17608       17510 - .. 0.93    no 0.42

                                                                               v                                       NC
                                                                              ((((((((,,,,,,<<<---<<<_.____>>>------>> CS
                                                                RF00162     1 cucUuAUcaAGAGaGGcgGAGGGA.cuGGCCCuaUGAAgC 39
                                                                              : CUUAU   GAGA:G: GAGGGA +UGGCCCU UGA :C
  riboswitch_11|Staphylococcus_aureus_HG003|-|SAOUHSC_00007|10446|33555 17608 UUCUUAUUGUGAGAAGUUGAGGGAcUUGGCCCUGUGAUAC 17569
                                                                              **************************************** PP

                                                                                         vv         vv                 NC
                                                                              >,,<<<-<<<<<<_________>>>-->>>>>>,,,,,,, CS
                                                                RF00162    40 CuCgGCAACCcccuauaauauagggAaGGUGCcAAuUCCu 79
                                                                              :UC:GCAACC:  U+UA    +  :A GGUGC:AA+ CC+
  riboswitch_11|Staphylococcus_aureus_HG003|-|SAOUHSC_00007|10446|33555 17568 UUCAGCAACCGACUUUA----UAGCACGGUGCUAAAACCA 17533
                                                                              **********5444443....367899************* PP

                                                                                                         v  NC
                                                                              ,<<<<_________>>>>,,,)))))))) CS
                                                                RF00162    80 gCcggauauuuuuuccgGaaAgAUaAgag 108
                                                                               C::G+      U+C::GAA  AUAAG :
  riboswitch_11|Staphylococcus_aureus_HG003|-|SAOUHSC_00007|10446|33555 17532 ACGAGU------UACUCGAAUGAUAAGUA 17510
                                                                              *****8......56*************** PP
Internal CM pipeline statistics summary:
----------------------------------------
Query sequence(s):                                               1  (46220 residues searched)
Query sequences re-searched for truncated hits:                  1  (890.1 residues re-searched, avg per model)
Target model(s):                                                36  (4015 consensus positions)
Windows   passing  local HMM SSV           filter:            2693  (0.3165); expected (0.35)
Windows   passing  local HMM Viterbi       filter:                  (off)
Windows   passing  local HMM Viterbi  bias filter:                  (off)
Windows   passing  local HMM Forward       filter:             276  (0.03347); expected (0.02)
Windows   passing  local HMM Forward  bias filter:             169  (0.02081); expected (0.02)
Windows   passing glocal HMM Forward       filter:              76  (0.01135); expected (0.02)
Windows   passing glocal HMM Forward  bias filter:              71  (0.01031); expected (0.02)
Envelopes passing glocal HMM envelope defn filter:              71  (0.004847); expected (0.02)
Envelopes passing  local CM  CYK           filter:               9  (0.00048); expected (0.0001)
Total CM hits reported:                                          3  (0.000112); includes 0 truncated hit(s)

# CPU time: 2.34u 0.12s 00:00:02.46 Elapsed: 00:00:01.41
//"""


if __name__ == "__main__":
    unittest.main()
