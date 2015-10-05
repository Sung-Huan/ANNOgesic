import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
import annogesiclib.meme as me
from mock_helper import gen_file, import_data
import annogesiclib.ribos as rb
from annogesiclib.ribos import Ribos


class Mock_func(object):

    def __init__(self):
        self.example = Example()

    def mock_run_infernal(self, infernal_path, e_value,
                          seq, type_, prefix):
        if type_ == "txt":
            gen_file("test_folder/output/test_RBS.txt", self.example.scan_file)
            return "test_folder/output/test_RBS.txt"
        else:
            gen_file("test_folder/output/test_RBS_rescan.txt", self.example.rescan_file)
            return "test_folder/output/test_RBS_rescan.txt"

    def mock_modify_table(self, first_table, output_all):
        pass

class TestRibos(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.mock = Mock_func()
        self.test_folder = "test_folder"
        self.gffs = os.path.join(self.test_folder, "gffs")
        self.fastas = os.path.join(self.test_folder, "fastas")
        self.out_folder = os.path.join(self.test_folder, "output")
        self.database = os.path.join(self.test_folder, "database")
        self.seq_path = os.path.join(self.test_folder, "seqs")
        self.tables = os.path.join(self.test_folder, "table")
        self.stat = os.path.join(self.test_folder, "stat")
        self.scan = os.path.join(self.test_folder, "scan")
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
            os.mkdir(self.gffs)
            os.mkdir(self.fastas)
            os.mkdir(self.out_folder)
            os.mkdir(self.database)
            os.mkdir(self.seq_path)
            os.mkdir(os.path.join(self.out_folder, "tmp_table"))
            os.mkdir(os.path.join(self.out_folder, "tmp_scan"))
            os.mkdir(self.tables)
            os.mkdir(self.scan)
            os.mkdir(self.stat)
        self.ribo = Ribos(self.gffs, self.fastas,
                          self.out_folder, self.database)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_scan_extract_rfam(self):
        self.ribo._run_infernal = self.mock.mock_run_infernal
        rb.modify_table = self.mock.mock_modify_table
        prefixs = []
        gen_file(os.path.join(self.gffs, "test.gff"), self.example.gff_file)
        gen_file(os.path.join(self.fastas, "test.fa"), self.example.fasta_file)
        gen_file(os.path.join(self.seq_path, "test.fa"), self.example.fasta_file)
        self.ribo._scan_extract_rfam(self.gffs, "0.001", self.seq_path, prefixs,
                           self.fastas, self.out_folder, "test", True, "test",
                           ["ATG"], 5, 9, 2)
        self.assertListEqual(prefixs, ["test"])
        self.assertTrue(os.path.exists(os.path.join(self.seq_path, "test_regenerate.fa")))

    def test_merge_results(self):
        gen_file(os.path.join(self.gffs, "test.gff"), self.example.gff_file) 
        gen_file(os.path.join(self.out_folder, "tmp_table/test_RBS.csv"), self.example.table)
        gen_file(os.path.join(self.out_folder, "tmp_scan/test_RBS.txt"), self.example.rescan_file)
        gen_file(os.path.join(self.test_folder, "ids"), self.example.ids)
        self.ribo._merge_results(self.gffs, self.scan, self.tables, self.out_folder,
                                 self.stat, os.path.join(self.test_folder, "ids"),
                                 3, self.out_folder, False)


class Example(object):

    ids = """RF00162	SAM	SAM	riboswitch box leader
RF00174	Cobalamin	Cobalamin	riboswitch
RF00634	SAM-IV	S adenosyl methionine SAM	riboswitch"""
    table = """riboswitch_5|test|+|SAOUHSC_00013|15948|16046	RF00162	1.6e-18	1	99
riboswitch_11|test|-|SAOUHSC_00007|27955|28053	RF00162	1.6e-18	1	99"""
    fasta_file = """>test
AAAAGATAGCCGCTCGCTAAATAGGAGATCCCTCGAAGAATGATATGTGTG"""
    gff_file = """##gff-version 3
test	Refseq	gene	10	16	.	+	.	ID=gene0;Name=gene_0;locus_tag=AAA_00001
test	Refseq	CDS	10	16	.	+	.	ID=cds0;Name=CDS_0;locus_tag=AAA_00001
test	Refseq	gene	30	40	.	-	.	ID=gene1;Name=gene_1;locus_tag=AAA_00002
test	Refseq	CDS	30	40	.	-	.	ID=cds1;Name=CDS_1;locus_tag=AAA_00002"""
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
