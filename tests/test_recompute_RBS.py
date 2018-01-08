import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import gen_file, import_data
import annogesiclib.recompute_RBS as rr


class TestRecomputeRBS(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_import_ribo(self):
        line = "(1) !   6.2e-18   74.3   0.0  RF00162      206    304 +  cm    no 0.42  -"
        ribos = []    
        seq_name = "test"
        rr.import_ribo(line, ribos, seq_name)
        self.assertListEqual(ribos, [
            {'e': 6.2e-18, 'score': 74.3, 'start': 206, 'name': 'RF00162',
             'detect': '!', 'seq_name': 'test', 'end': 304}])

    def test_print_file(self):
        out_t = StringIO()
        out_s = StringIO()
        ribos = [{'e': 6.2e-18, 'score': 74.3, 'start': 206, 'name': 'RF00162',
                  'detect': '!',
                  'seq_name': 'test_1|test_strain|+|AAA_0001|206|304',
                  'end': 304}]
        seqs = [{"name": "test_1|test_strain|+|AAA_0001|206|304",
                 "seq": "ATATAGCGTAGCCGACGTCGCAT"}]
        seq_name = "test_1|test_strain|+|AAA_0001|206|304"
        rr.print_file(ribos, out_t, out_s, seq_name, seqs)
        self.assertEqual(
            out_t.getvalue(),
            "test_1\ttest_strain\t+\tAAA_0001\t206\t304\tRF00162\t6.2e-18\t74.3\t206\t304\n")
        self.assertEqual(
            out_s.getvalue(),
            ">test_1|test_strain|+|AAA_0001|411|509\n\n")

    def test_regenerate_seq(self):
        out_table = os.path.join(self.test_folder, "table")
        out_seq = os.path.join(self.test_folder, "seq")
        align_file = os.path.join(self.test_folder, "align")
        seq_file = os.path.join(self.test_folder, "ribo_seq")
        gen_file(align_file, self.example.scan_file)
        gen_file(seq_file, self.example.seq_file)
        rr.regenerate_seq(align_file, seq_file, out_table, out_seq)
        data = import_data(out_table)
        self.assertEqual(
            "\n".join(data),
            "riboswitch_5\tStaphylococcus_aureus_HG003\t+\tSAOUHSC_00013\t10\t16\tRF00162\t6.2e-18\t74.3\t5\t12")
        data = import_data(out_seq)
        self.assertEqual(
            "\n".join(data), (
            ">riboswitch_5|Staphylococcus_aureus_HG003|+|SAOUHSC_00013|14|21\nATTATTAC\n"
            ">riboswitch_5|Staphylococcus_aureus_HG003|+|SAOUHSC_00013|14|21\nATTATTAC"))

    def test_compare_first_result(self):
        out = StringIO()
        cutoff = "e_0.01"
        ribos = [{
            'e': 6.2e-18, 'start': 206, 'name': 'RF00162', 'detect': '!',
            'seq_name': 'test_1|test_strain|+|AAA_0001|206|304', 'end': 304}]
        firsts = [{
            'e': 6.2e-10, 'start': 100, 'name': 'RF00160', 'detect': '!',
            'seq_name': 'test_2|test_strain|+|AAA_0002|100|150', 'end': 150}]
        seq_name = "test_strain"
        extras = [{
            'e': 6.2e-10, 'start': 100, 'name': 'RF00160', 'detect': '!',
            'seq_name': 'test_2|test_strain|+|AAA_0002|100|150', 'end': 150}]
        rr.compare_first_result(ribos, firsts, seq_name, out, extras, cutoff)
        self.assertListEqual(extras, [
            {'start': 100, 'name': 'RF00160', 'detect': '!', 'end': 150,
             'e': 6.2e-10,
             'seq_name': 'test_2|test_strain|+|AAA_0002|100|150'},
            {'start': 206, 'name': 'RF00162', 'detect': '!', 'end': 304,
             'e': 6.2e-18,
             'seq_name': 'test_1|test_strain|+|AAA_0001|206|304'}])

    def test_reextract_rbs(self):
        align_file = os.path.join(self.test_folder, "align")
        first_file = os.path.join(self.test_folder, "first")
        output_file = os.path.join(self.test_folder, "output")
        first_content = """riboswitch_5\tStaphylococcus_aureus_HG003\t+\tSAOUHSC_00013\t10\t16	RF00162	6.2e-18	74.3	5	12"""
        gen_file(align_file, self.example.scan_file)
        gen_file(first_file, first_content)
        cutoff = "e_0.01"
        rr.reextract_rbs(align_file, first_file, output_file, cutoff)
        data = import_data(output_file)
        self.assertEqual("\n".join(data), first_content)
        first_content = """riboswitch_5\tStaphylococcus_aureus_HG003\t+\tSAOUHSC_00013\t10\t16	RF00178	6.2e-20	74.3	13	17"""
        gen_file(first_file, first_content)
        rr.reextract_rbs(align_file, first_file, output_file, cutoff)
        data = import_data(output_file)
        self.assertEqual("\n".join(data), 
"""riboswitch_5\tStaphylococcus_aureus_HG003\t+\tSAOUHSC_00013\t10\t16	RF00162	6.2e-18	74.3	5	12
riboswitch_5\tStaphylococcus_aureus_HG003\t+\tSAOUHSC_00013\t10\t16	RF00178	6.2e-20	74.3	13	17""")

class Example(object):

    seq_file = """>riboswitch_5|Staphylococcus_aureus_HG003|+|SAOUHSC_00013|10|16
ATAGATTATTACCACCCCCAGATACGATGACCGTATGA"""
        
    scan_file = """Query:       riboswitch_5|Staphylococcus_aureus_HG003|+|SAOUHSC_00013|10|16  [L=7]
Hit scores:
 rank     E-value  score  bias  modelname  start    end   mdl trunc   gc  description
 ----   --------- ------ -----  --------- ------ ------   --- ----- ----  -----------
  (1) !   6.2e-18   74.3   0.0  RF00162      5    12 +  cm    no 0.42  -
 ------ inclusion threshold ------
  (2) ?      0.24   11.1   0.0  RF01689      5    12 -  cm    no 0.43  -


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

# CPU time: 1.51u 0.04s 00:00:01.55 Elapsed: 00:00:01.38
//"""

if __name__ == "__main__":
    unittest.main()

