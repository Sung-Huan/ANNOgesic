#!/usr/bin/python

import unittest
import os
import sys
import shutil
sys.path.append(".")
from io import StringIO
from mock_gff3 import Create_generator
from mock_helper import convert_dict, gen_file
from mock_args_container import MockClass
import annogesiclib.merge_rnaplex_rnaup as mrr


class TestMergeRNAplexRNAup(unittest.TestCase):

    def setUp(self):
        self.test_folder = "test_project"
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)
        os.mkdir(self.test_folder)
        self.example = Example()
        self.mock_args = MockClass()

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_detect_energy(self):
        srna = {"energy": -2}
        mrr.detect_energy(self.example.out_rna_txt, srna)
        self.assertDictEqual(srna, {'energy': -5.3})
        srna = {"energy": -8}
        mrr.detect_energy(self.example.out_rna_txt, srna)
        self.assertDictEqual(srna, {'energy': -8.0})

    def test_print_rank_one(self):
        out = StringIO()
        args_tar = self.mock_args.mock()
        args_tar.top = 2
        args_tar.tar_start = 20
        args_tar.tar_end = 15
        mrr.print_rank_one(self.example.srnas, out, "RNAplex", self.example.gffs, self.example.srna_gffs, args_tar)
        datas = convert_dict(out.getvalue().split("\n"))
        refs = convert_dict(self.example.out_print.split("\n"))
        self.assertDictEqual(datas, refs)

    def test_read_table(self):
        rnaplex = os.path.join(self.test_folder, "rnaplex")
        rnaup = os.path.join(self.test_folder, "rnaup")
        gen_file(rnaplex, self.example.rnaplex)
        gen_file(rnaup, self.example.rnaup)
        srnas = mrr.read_table(self.example.gffs, rnaplex, rnaup)
        self.assertDictEqual(srnas, {'RNAup': {'srna352': [{'target': 'srna1023', 'energy': 0},
                                    {'tar_pos': '571,576', 'target': 'SAOUHSC_00001|dnaA',
                                     'energy': -4.87, 'srna_pos': '20,25'},
                                    {'tar_pos': '14,30', 'target': 'SAOUHSC_00002',
                                     'energy': -5.91, 'srna_pos': '11,26'}]},
                                     'RNAplex': {'srna1023': [{'tar_pos': '571,576',
                                     'target': 'SAOUHSC_00001|dnaA', 'energy': -5.3, 'srna_pos': '20,25'}],
                                     'srna352': [{'tar_pos': '163,170', 'target': 'SAOUHSC_00001|dnaA',
                                     'energy': -1.91, 'srna_pos': '24,31'}]}})

    def test_get_srna_name(self):
        output = mrr.get_srna_name(self.example.srna_gffs, "srna0")
        self.assertEqual(output[0], 'sRNA_0')
        self.assertEqual(output[1].start, 6)

    def test_get_target_info(self):
        output = mrr.get_target_info(self.example.gffs, "AAA_00001")
        self.assertEqual(output.start, 100)

    def test_merge_base_rnaplex(self):
        args_tar = self.mock_args.mock()
        args_tar.top = 2
        args_tar.tar_start = 20
        args_tar.tar_end = 15
        merges = []
        overlap = mrr.merge_base_rnaplex(self.example.srnas, self.example.srna_gffs, args_tar, self.example.gffs, merges)
        output = [['sRNA_0', 'aaa', '6-15', '7-15', '7-15', '+', 'AAA_00001',
                   '100-150', '89-94', '89-94', '+', '-6.5', '1', '-6.5', '1'],
                  ['sRNA_0', 'aaa', '6-15', '7-12', '7-15', '+', 'AAA_00002|dnaA',
                   '2348-2934', '2330-2342', '2337-2342', '+', '-3.5', '2', '-3.5', '2'],
                  ['sRNA_1', 'aaa', '1258-2234', '1259-1267', '1259-1267', '+', 'AAA_00003',
                   '5544-5597', '5550-5545', '5550-5545', '-', '-10.5', '1', '-10.5', '1'],
                  ['sRNA_2', 'aaa', '3544-6517', '6508-6516', '6508-6516', '-', 'AAA_00001',
                   '100-150', '89-94', '89-94', '+', '-23.5', '1', '-23.5', '1']]
        count = 0
        for out in output:
            for data in overlap:
                if out == data:
                    count += 1
        self.assertEqual(count, 4)
        count = 0
        for out in output:
            for data in merges:
                if out == data:
                    count += 1
        self.assertEqual(count, 4)

    def test_merge_base_rnaup(self):
        args_tar = self.mock_args.mock()
        args_tar.top = 2
        args_tar.tar_start = 20
        args_tar.tar_end = 15
        srnas = {"RNAplex": {"srna0": [{"target": "AAA_00001", "energy": -6.5, "rank": 1, "srna_pos": "2,10", "tar_pos": "10,15"},
                                       {"target": "AAA_00002|dnaA", "energy": -3.5, "rank": 2, "srna_pos": "2,10", "tar_pos": "10,15"}],
                             "srna1": [{"target": "AAA_00003", "energy": -10.5, "rank": 1, "srna_pos": "2,10", "tar_pos": "10,15"}],
                             "srna2": [{"target": "AAA_00001", "energy": -23.5, "rank": 1, "srna_pos": "2,10", "tar_pos": "10,15"},
                                       {"target": "AAA_00002|dnaA", "energy": -3.43, "rank": 3, "srna_pos": "2,10", "tar_pos": "10,15"},
                                       {"target": "AAA_00003", "energy": -6.5, "rank": 2, "srna_pos": "2,10", "tar_pos": "10,15"}]},
                 "RNAup": {"srna0": [{"target": "AAA_00001", "energy": -6.5, "rank": 1, "srna_pos": "2,10", "tar_pos": "10,15"},
                                     {"target": "AAA_00002|dnaA", "energy": -3.5, "rank": 2, "srna_pos": "2,10", "tar_pos": "10,15"}],
                           "srna1": [{"target": "AAA_00003", "energy": -10.5, "rank": 1, "srna_pos": "2,10", "tar_pos": "10,15"}],
                           "srna2": [{"target": "AAA_00001", "energy": -23.5, "rank": 1, "srna_pos": "2,10", "tar_pos": "10,15"}]}}
        merges = []
        mrr.merge_base_rnaup(srnas, self.example.srna_gffs, args_tar, self.example.gffs, merges)
        output = [['sRNA_1', 'aaa', '1258-2234', '1259-1267', '1259-1267', '+', 'AAA_00003',
                   '5544-5597', '5550-5545', '5550-5545', '-', '-10.5', '1', '-10.5', '1'],
                  ['sRNA_2', 'aaa', '3544-6517', '6508-6516', '6508-6516', '-', 'AAA_00001',
                   '100-150', '89-94', '89-94', '+', '-23.5', '1', '-23.5', '1'],
                  ['sRNA_0', 'aaa', '6-15', '7-15', '7-15', '+', 'AAA_00001',
                   '100-150', '89-94', '89-94', '+', '-6.5', '1', '-6.5', '1'],
                  ['sRNA_0', 'aaa', '6-15', '7-15', '7-15', '+', 'AAA_00002|dnaA',
                   '2348-2934', '2337-2342', '2337-2342', '+', '-3.5', '2', '-3.5', '2']]
        count = 0
        for out in output:
            for data in merges:
                if out == data:
                    count += 1
        self.assertEqual(count, 4)

class Example(object):
    srnas = {"RNAplex": {"srna0": [{"target": "AAA_00001", "energy": -6.5, "rank": 1, "srna_pos": "2,10", "tar_pos": "10,15"},
                                   {"target": "AAA_00002|dnaA", "energy": -3.5, "rank": 2, "srna_pos": "2,7", "tar_pos": "3,15"}],
                         "srna1": [{"target": "AAA_00003", "energy": -10.5, "rank": 1, "srna_pos": "2,10", "tar_pos": "10,15"}],
                         "srna2": [{"target": "AAA_00001", "energy": -23.5, "rank": 1, "srna_pos": "2,10", "tar_pos": "10,15"},
                                   {"target": "AAA_00002|dnaA", "energy": -3.43, "rank": 3, "srna_pos": "2,10", "tar_pos": "10,15"},
                                   {"target": "AAA_00003", "energy": -6.5, "rank": 2, "srna_pos": "2,10", "tar_pos": "10,15"}]},
             "RNAup": {"srna0": [{"target": "AAA_00001", "energy": -6.5, "rank": 1, "srna_pos": "2,10", "tar_pos": "10,15"},
                                 {"target": "AAA_00002|dnaA", "energy": -3.5, "rank": 2, "srna_pos": "2,10", "tar_pos": "10,15"}],
                       "srna1": [{"target": "AAA_00003", "energy": -10.5, "rank": 1, "srna_pos": "2,10", "tar_pos": "10,15"}],
                       "srna2": [{"target": "AAA_00001", "energy": -23.5, "rank": 1, "srna_pos": "2,10", "tar_pos": "10,15"}]}}
    srna_dict = [{"start": 6, "end": 15, "phase": ".",
                  "strand": "+", "seq_id": "aaa", "score": ".",
                  "source": "Refseq", "feature": "sRNA"},
                 {"start": 1258, "end": 2234, "phase": ".",
                  "strand": "+", "seq_id": "aaa", "score": ".",
                  "source": "Refseq", "feature": "sRNA"},
                 {"start": 3544, "end": 6517, "phase": ".",
                  "strand": "-", "seq_id": "aaa", "score": ".",
                  "source": "Refseq", "feature": "sRNA"}]
    attributes_srna = [{"ID": "srna0", "Name": "sRNA_0"},
                       {"ID": "srna1", "Name": "sRNA_1"},
                       {"ID": "srna2", "Name": "sRNA_2"}]
    gff_dict = [{"start": 100, "end": 150, "phase": ".",
                 "strand": "+", "seq_id": "aaa", "score": ".",
                 "source": "Refseq", "feature": "CDS"},
                {"start": 2348, "end": 2934, "phase": ".",
                 "strand": "+", "seq_id": "aaa", "score": ".",
                 "source": "Refseq", "feature": "CDS"},
                {"start": 5544, "end": 5597, "phase": ".",
                 "strand": "-", "seq_id": "aaa", "score": ".",
                 "source": "Refseq", "feature": "CDS"}]
    attributes_gff = [{"ID": "cds0", "Name": "CDS_0", "locus_tag": "AAA_00001"},
                      {"ID": "cds0", "Name": "CDS_1", "locus_tag": "AAA_00002"},
                      {"ID": "cds1", "Name": "CDS_2", "locus_tag": "AAA_00003"}]
    srna_gffs = []
    gffs = []
    for index in range(0, 3):
        srna_gffs.append(Create_generator(srna_dict[index], attributes_srna[index], "gff"))
        gffs.append(Create_generator(gff_dict[index], attributes_gff[index], "gff"))    
    out_rna_txt = """>SAOUHSC_00001|dnaA
>srna1023
((((((&)))))) 571,576 :  20,25  (-5.30 = -7.89 +  0.18 +  2.41)"""
    out_print = """sRNA	strain	sRNA_position	sRNA_interacted_position_RNAplex	sRNA_strand	target	target_position	target_interacted_position_RNAplex	target_strand	energy_RNAplex	rank_RNAplex
sRNA_1	aaa	1258-2234	1259-1267	+	AAA_00003	5544-5597	5550-5545	-	-10.5	1
sRNA_2	aaa	3544-6517	6508-6516	-	AAA_00001	100-150	89-94	+	-23.5	1
sRNA_2	aaa	3544-6517	6508-6516	-	AAA_00003	5544-5597	5550-5545	-	-6.5	2
sRNA_0	aaa	6-15	7-15	+	AAA_00001	100-150	89-94	+	-6.5	1
sRNA_0	aaa	6-15	7-12	+	AAA_00002|dnaA	2348-2934	2330-2342	+	-3.5	2
"""
    rnaup = """>srna1023
>SAOUHSC_00001|dnaA
.(((((&))))). 571,576 :  20,25  (-4.87 = -8.00 + 0.31 + 2.81)
AACCUC&GGGGUU
>SAOUHSC_00002
(((..((((((((((((&)))))))))))).)))  14,30  :  11,26  (-5.91 = -13.15 + 4.20 + 3.05)
GAAGAUCCUAUUUUUAA&UUAAAAAUGGGGGUUC
"""
    rnaplex = """>SAOUHSC_00001|dnaA
>srna1023
((((((&)))))) 571,576 :  20,25  (-5.30 = -7.89 +  0.18 +  2.41)
>SAOUHSC_00001|dnaA
>srna352
((((((((&)))))))) 163,170 :  24,31  (-1.91 = -8.31 +  0.60 +  5.80)
"""    

if __name__ == "__main__":
    unittest.main()
