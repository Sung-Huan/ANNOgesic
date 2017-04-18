import sys
import math
import csv
import os
import shutil
import unittest
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data, gen_file
import annogesiclib.transcript_SNP as ts
from mock_args_container import MockClass


class TestTranscripSNP(unittest.TestCase):

    def setUp(self):
        self.mock_args = MockClass()
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_import_data(self):
        snp_file = os.path.join(self.test_folder, "snp")
        gen_file(snp_file, self.example.snp_file)
        depth_file = os.path.join(self.test_folder, "depth")
        gen_file(depth_file, self.example.depth_file)
        args = self.mock_args.mock()
        args.depth_s = "n_10"
        args.depth_b = "a_2"
        args.dp4_sum = "n_10"
        args.dp4_frac = 0.5
        args.idv = "n_10"
        args.imf = 0.5
        args.filters = ["VDB_s0.1"]
        args.min_sample = 2
        max_quals, snps, dess, raw_snps = ts.import_data(
            snp_file, args, 2, depth_file, 2)
        self.assertDictEqual(max_quals, {
            'NC_007795.1': 98.0, 'All_genome': 98.0})
        self.assertListEqual(snps, [
            {'dp4_frac': 1.0, 'strain': 'NC_007795.1', 'filter': '.',
             'indel': -1, 'pos': 1, 'id': '.',
             'all_info': ("NC_007795.1\t1\t.\tC\tA\t98\t.\tDP=89;DP4=0,0,"
                          "60,9;VDB=8.46526e-15\tGT:PL:DP\t1/1:125,184,0:87"),
             'qual': 98.0,
             'info': ['DP=89', 'DP4=0,0,60,9', 'VDB=8.46526e-15'],
             'alt': 'A', 'ref': 'C', 'frac': -1, 'depth': 89, 'dp4_sum': 69},
            {'dp4_frac': 1.0, 'strain': 'NC_007795.1', 'filter': '.',
             'indel': 22, 'pos': 6, 'id': '.',
             'all_info': ("NC_007795.1\t6\t.\tA\tAA\t26.9515\t.\tINDEL;IDV=22;"
                          "IMF=0.536585;DP=41;VDB=9.36323e-14;DP4=0,0,40,0\t"
                          "GT:PL:DP\t0/1:60,0,55:40"),
             'qual': 26.9515,
             'info': ['INDEL', 'IDV=22', 'IMF=0.536585', 'DP=41',
                      'VDB=9.36323e-14', 'DP4=0,0,40,0'],
             'alt': 'AA', 'ref': 'A', 'frac': 0.536585, 'depth': 41,
             'dp4_sum': 40}])

    def test_check_overlap(self):
        snps = {"test": []}
        overlaps = [{"test": []}]
        ts.check_overlap(snps, overlaps)
        self.assertListEqual(overlaps, [{'test': [], 'print': True}])
        self.assertDictEqual(snps, {'test': [{'test': [], 'print': True}]})

    def test_overlap_position(self):
        qual_snps = [{'filter': '.', 'pos': 22181, 'alt': 'A',
                      'frac': -1, 'depth': 89, 'indel': -1,
                      'info': 'MQ=20', 'id': '.', 'qual': 98.0,
                      'ref': 'CA', 'strain': 'NC_007795.1',
                      'all_info': ("NC_007795.1\t22181\t.\tC\tA\t98\t.\t"
                                   "DP=89;VDB=8.46526e-15;SGB=-0.693147\t"
                                   "GT:PL:DP\t1/1:125,184,0:87")},
                     {'filter': '.', 'pos': 22182, 'alt': 'C',
                      'frac': -1, 'depth': 89, 'indel': -1,
                      'info': 'MQ=20', 'id': '.', 'qual': 98.0,
                      'ref': 'A', 'strain': 'NC_007795.1', 
                      'all_info': ("NC_007795.1\t22182\t.\tC\tA\t98\t.\t"
                                   "DP=89;VDB=8.46526e-15;SGB=-0.693147\t"
                                   "GT:PL:DP\t1/1:125,184,0:87")},
                     {'filter': '.', 'pos': 30000, 'alt': 'A',
                      'frac': -1, 'depth': 89, 'indel': -1,
                      'info': 'MQ=20', 'id': '.', 'qual': 98.0,
                      'ref': 'C', 'strain': 'NC_007795.1', 
                      'all_info': ("NC_007795.1\t30000\t.\tC\tA\t98\t.\t"
                                   "DP=89;VDB=8.46526e-15;SGB=-0.693147\t"
                                   "GT:PL:DP\t1/1:125,184,0:87")}]
        conflicts, nooverlap = ts.overlap_position(qual_snps)
        self.assertListEqual(conflicts, [[
            {'strain': 'NC_007795.1', 'info': 'MQ=20',
             'indel': -1, 'qual': 98.0, 'ref': 'CA', 'frac': -1,
             'alt': 'A', 'depth': 89, 'print': True, 'pos': 22181,
             'filter': '.', 'id': '.',
             'all_info': ("NC_007795.1\t22181\t.\tC\tA\t98\t.\tDP=89;"
                          "VDB=8.46526e-15;SGB=-0.693147\tGT:PL:DP\t"
                          "1/1:125,184,0:87")},
            {'strain': 'NC_007795.1', 'info': 'MQ=20', 'indel': -1,
             'qual': 98.0, 'ref': 'A', 'frac': -1, 'alt': 'C',
             'depth': 89, 'print': True, 'pos': 22182, 'filter': '.',
             'id': '.',
             'all_info': ("NC_007795.1\t22182\t.\tC\tA\t98\t.\tDP=89;"
                          "VDB=8.46526e-15;SGB=-0.693147\tGT:PL:DP\t"
                          "1/1:125,184,0:87")}]])
        self.assertDictEqual(nooverlap, {1: [
            {'strain': 'NC_007795.1', 'info': 'MQ=20', 'indel': -1,
             'qual': 98.0, 'ref': 'CA', 'frac': -1, 'alt': 'A',
             'depth': 89, 'print': True, 'pos': 22181, 'filter': '.',
             'id': '.',
             'all_info': ("NC_007795.1\t22181\t.\tC\tA\t98\t.\tDP=89;"
                          "VDB=8.46526e-15;SGB=-0.693147\tGT:PL:DP\t"
                          "1/1:125,184,0:87")},
            {'strain': 'NC_007795.1', 'info': 'MQ=20', 'indel': -1,
             'qual': 98.0, 'ref': 'C', 'frac': -1, 'alt': 'A',
             'depth': 89, 'print': True, 'pos': 30000, 'filter': '.',
             'id': '.',
             'all_info': ("NC_007795.1\t30000\t.\tC\tA\t98\t.\tDP=89;"
                          "VDB=8.46526e-15;SGB=-0.693147\tGT:PL:DP\t"
                          "1/1:125,184,0:87")}],
                                         2: [
            {'strain': 'NC_007795.1', 'info': 'MQ=20', 'indel': -1,
             'qual': 98.0, 'ref': 'A', 'frac': -1, 'alt': 'C',
             'depth': 89, 'print': True,
             'pos': 22182, 'filter': '.', 'id': '.',
             'all_info': ("NC_007795.1\t22182\t.\tC\tA\t98\t.\tDP=89;"
                          "VDB=8.46526e-15;SGB=-0.693147\tGT:PL:DP\t"
                          "1/1:125,184,0:87")},
            {'strain': 'NC_007795.1', 'info': 'MQ=20', 'indel': -1,
             'qual': 98.0, 'ref': 'C', 'frac': -1, 'alt': 'A',
             'depth': 89, 'print': True, 'pos': 30000, 'filter': '.',
             'id': '.',
             'all_info': ("NC_007795.1\t30000\t.\tC\tA\t98\t.\tDP=89;"
                          "VDB=8.46526e-15;SGB=-0.693147\tGT:PL:DP\t"
                          "1/1:125,184,0:87")}]})

    def test_stat(self):
        stat_file = os.path.join(self.test_folder, "stat")
        max_quals = {'NC_007795.1': 98.0, 'All_genome': 98.0}
        trans_snps = [{'filter': '.', 'pos': 22181, 'alt': 'A',
                       'frac': -1, 'depth': 89, 'indel': -1,
                       'info': 'MQ=20', 'id': '.', 'qual': 98.0,
                       'ref': 'C', 'strain': 'NC_007795.1',
                       'all_info': ("NC_007795.1\t22181\t.\tC\tA\t98\t.\t"
                                    "DP=89;VDB=8.46526e-15\tGT:PL:DP\t"
                                    "1/1:125,184,0:87")}]
        args = self.mock_args.mock()
        args.depth = 50
        args.fraction = 0.3
        args.quality = 20
        ts.stat(max_quals, trans_snps, 2, stat_file,
                self.test_folder + "/test", args, "best.csv")
        datas = import_data(stat_file + "_best.csv")
        self.assertEqual("\n".join(datas), self.example.stat)

    def test_plot_bar(self):
        ts.plot_bar([3, 10, 30, 45, 50], "NC_007795.1",
                    self.test_folder + "/test", "best.png")
        self.assertTrue(os.path.exists(os.path.join(
            self.test_folder, "test_NC_007795.1_SNP_QUAL_best.png")))

    def test_read_fasta(self):
        fasta_file = os.path.join(self.test_folder, "NC_007795.1.fa")
        gen_file(fasta_file, self.example.fasta_file)
        seqs = ts.read_fasta(fasta_file)
        self.assertListEqual(seqs, [{
            'NC_007795.1': 'AAATATATCAGCACCGTAGACGATAGAGTAGTAC'}])

    def test_gen_ref(self):
        refs = []
        snps = [{'filter': '.', 'pos': 22181, 'alt': 'A',
                 'frac': -1, 'depth': 89, 'indel': -1,
                 'info': 'MQ=20', 'id': '.', 'qual': 98.0, 
                 'ref': 'C', 'strain': 'NC_007795.1',
                 'all_info': ("NC_007795.1\t22181\t.\tC\tA\t98\t.\tDP=89;"
                              "VDB=8.46526e-15\tGT:PL:DP\t1/1:125,184,0:87")},
                {'filter': '.', 'pos': 22500, 'alt': 'A',
                 'frac': -1, 'depth': 89, 'indel': -1,
                 'info': 'MQ=20', 'id': '.', 'qual': 98.0, 
                 'ref': 'C', 'strain': 'NC_007795.1',
                 'all_info': ("NC_007795.1\t22500\t.\tC\tA\t98\t.\tDP=89;"
                              "VDB=8.46526e-15\tGT:PL:DP\t1/1:125,184,0:87")}]
        refs = ts.gen_ref(snps, 1, refs, 1)
        self.assertListEqual(refs, ['1:A', '1:A'])
        snps = [{'filter': '.', 'pos': 22181, 'alt': 'A',
                 'frac': -1, 'depth': 89, 'indel': -1,
                 'info': 'MQ=20', 'id': '.', 'qual': 98.0,
                 'ref': 'C', 'strain': 'NC_007795.1',
                 'all_info': ("NC_007795.1\t22181\t.\tC\tA\t98\t.\tDP=89;"
                              "VDB=8.46526e-15\tGT:PL:DP\t1/1:125,184,0:87")},
                {'filter': '.', 'pos': 22500, 'alt': 'A',
                 'frac': -1, 'depth': 89, 'indel': -1,
                 'info': 'MQ=20', 'id': '.', 'qual': 98.0,
                 'ref': 'C', 'strain': 'NC_007795.1',
                 'all_info': ("NC_007795.1\t22500\t.\tC\tA\t98\t.\tDP=89;"
                              "VDB=8.46526e-15\tGT:PL:DP\t1/1:125,184,0:87")}]
        refs = ts.gen_ref(snps, 1, refs, 2)
        self.assertListEqual(refs, [
            '1:A_1:A', '1:A_1:A', '1:A_1:A', '1:A_1:A'])

    def test_change(self):
        snp = {'filter': '.', 'pos': 1, 'alt': 'A',
                'frac': -1, 'depth': 89, 'indel': -1,
                'info': 'MQ=20', 'id': '.', 'qual': 98.0,
                'ref': 'C', 'strain': 'NC_007795.1',
                'all_info': ("NC_007795.1\t1\t.\tC\tA\t98\t.\tDP=89;"
                             "VDB=8.46526e-15\tGT:PL:DP\t1/1:125,184,0:87")}
        seq = {"num_mod": 3, "seq": "CCCCATATCAGCACCGTAGACGATAGAGTAGTAC"}
        ts.change(snp, seq)
        self.assertDictEqual(seq, {
            'num_mod': 3, 'seq': 'CCCaATATCAGCACCGTAGACGATAGAGTAGTAC'})

    def test_print_file(self):
        refs = {'NC_007795.1': ['1:A', '1:GT']}
        conflicts = [[{'all_info': ("NC_007795.1\t1\t.\tCA\tA,GT\t98\t.\tDP=89;"
                                    "VDB=8.46526e-15\tGT:PL:DP\t1/1:125,"
                                    "184,0:87"),
                       'filter': '.', 'id': '.', 'frac': -1, 'indel': -1,
                       'alt': 'A,GT', 'info': 'VDB=8.46526e-15', 'qual': 98.0,
                       'ref': 'CA', 'strain': 'NC_007795.1', 'depth': 89,
                       'pos': 1, 'print': True},
                      {'all_info': ("NC_007795.1\t2\t.\tA\tAA\t26.9515\t.\t"
                                    "INDEL;IDV=22;IMF=0.536585;DP=41;"
                                    "VDB=9.36323e-14 GT:PL:DP\t0/1:60,0,"
                                    "55:40"),
                       'filter': '.', 'id': '.', 'frac': 0.536585, 'indel': 22,
                       'alt': 'AA', 'info': 'VDB=9.36323e-14 GT:PL:DP',
                       'qual': 26.9515, 'ref': 'A', 'strain': 'NC_007795.1',
                       'depth': 41, 'pos': 2, 'print': True}]]
        values = [{'all_info': ("NC_007795.1\t1\t.\tCA\tA,GT\t98\t.\tDP=89;"
                                "VDB=8.46526e-15\tGT:PL:DP\t1/1:125,184,0:87"),
                   'filter': '.', 'id': '.', 'frac': -1, 'indel': -1,
                   'alt': 'A,GT', 'info': 'VDB=8.46526e-15', 'qual': 98.0,
                   'ref': 'CA', 'strain': 'NC_007795.1', 'depth': 89, 'pos': 1,
                   'print': True},
                  {'all_info': ("NC_007795.1\t7\t.\tC\tA\t98\t.\tDP=89;"
                                "VDB=8.46526e-15\tGT:PL:DP\t1/1:125,184,0:87"),
                   'filter': '.', 'id': '.', 'frac': -1, 'indel': -1,
                   'alt': 'A', 'info': 'VDB=8.46526e-15', 'qual': 98.0,
                   'ref': 'C', 'strain': 'NC_007795.1', 'depth': 89, 'pos': 7,
                   'print': True}]
        mod_seq_init = {'genome': 'NC_007795.1', 'num_mod': 0,
                        'seq': 'CAGTACCCTCAGCACCGTAGACGATAGAGTAGTAC'}
        mod_seqs = [{'genome': 'NC_007795.1', 'num_mod': -1,
                     'seq': 'aGTACaCTCAGCACCGTAGACGATAGAGTAGTAC'},
                    {'genome': 'NC_007795.1', 'num_mod': 0,
                     'seq': 'gtGTACaCTCAGCACCGTAGACGATAGAGTAGTAC'}]
        out_ref = StringIO()
        out_seq = os.path.join(self.test_folder, "seq")
        ts.print_file(refs, out_ref, conflicts, 1, values, mod_seq_init,
                      mod_seqs, out_seq, "NC_007795.1")
        self.assertEqual(
            out_ref.getvalue(),
            "1\t1\t1\tNC_007795.1\tNC_007795.1\n")
        self.assertTrue(os.path.exists(os.path.join(
            self.test_folder, "seq_NC_007795.1_1_1.fa")))
        self.assertTrue(os.path.exists(os.path.join(
            self.test_folder, "seq_NC_007795.1_1_2.fa")))

    def test_gen_new_fasta(self):
        out_ref = StringIO()
        out_seq = os.path.join(self.test_folder, "seq")
        nooverlap = {1: [{'strain': 'NC_007795.1', 'print': True, 'id': '.',
                          'alt': 'A,GT', 'filter': '.', 'frac': -1,
                          'ref': 'CA', 'depth': 89, 'info': 'VDB=8.46526e-15',
                          'indel': -1, 'qual': 98.0, 'pos': 1,
                          'all_info': ("NC_007795.1\t1\t.\tCA\tA,GT\t98\t.\t"
                                       "DP=89;VDB=8.46526e-15\tGT:PL:DP\t"
                                       "1/1:125,184,0:87")},
                         {'strain': 'NC_007795.1', 'print': True, 'id': '.',
                          'alt': 'A', 'filter': '.', 'frac': -1, 'ref': 'C',
                          'depth': 89, 'info': 'VDB=8.46526e-15',
                          'indel': -1, 'qual': 98.0, 'pos': 7,
                          'all_info': ("NC_007795.1\t7\t.\tC\tA\t98\t.\tDP=89;"
                                       "VDB=8.46526e-15\tGT:PL:DP\t1/1:125,"
                                       "184,0:87")}],
                     2: [{'strain': 'NC_007795.1', 'print': True, 'id': '.',
                          'alt': 'AA', 'filter': '.', 'frac': 0.536585,
                          'ref': 'A', 'depth': 41,
                          'info': 'VDB=9.36323e-14 GT:PL:DP',
                          'indel': 22, 'qual': 26.9515, 'pos': 2,
                          'all_info': ("NC_007795.1\t2\t.\tA\tAA\t26.9515\t.\t"
                                       "INDEL;IDV=22;IMF=0.536585;DP=41;"
                                       "VDB=9.36323e-14 GT:PL:DP\t0/1:60,"
                                       "0,55:40")},
                         {'strain': 'NC_007795.1', 'print': True, 'id': '.',
                          'alt': 'A', 'filter': '.', 'frac': -1, 'ref': 'C',
                          'depth': 89, 'info': 'VDB=8.46526e-15',
                          'indel': -1, 'qual': 98.0, 'pos': 7,
                          'all_info': ("NC_007795.1\t7\t.\tC\tA\t98\t.\tDP=89;"
                                       "VDB=8.46526e-15\tGT:PL:DP\t1/1:125,"
                                       "184,0:87")}]}
        seqs = [{'NC_007795.1': 'CAGTACCCTCAGCACCGTAGACGATAGAGTAGTAC'}]
        conflicts = [[{'strain': 'NC_007795.1', 'print': True, 'id': '.',
                       'alt': 'A,GT', 'filter': '.', 'frac': -1, 'ref': 'CA',
                       'depth': 89, 'info': 'VDB=8.46526e-15', 'indel': -1,
                       'qual': 98.0, 'pos': 1,
                       'all_info': ("NC_007795.1\t1\t.\tCA\tA,GT\t98\t.\t"
                                    "DP=89;VDB=8.46526e-15\tGT:PL:DP\t1/1:125,"
                                    "184,0:87")},
                      {'strain': 'NC_007795.1', 'print': True, 'id': '.',
                       'alt': 'AA', 'filter': '.', 'frac': 0.536585, 'ref': 'A',
                       'depth': 41, 'info': 'VDB=9.36323e-14 GT:PL:DP',
                       'indel': 22, 'qual': 26.9515, 'pos': 2,
                       'all_info': ("NC_007795.1\t2\t.\tA\tAA\t26.9515\t.\t"
                                    "INDEL;IDV=22;IMF=0.536585;DP=41;"
                                    "VDB=9.36323e-14 GT:PL:DP\t0/1:60,0,"
                                    "55:40")}]]
        ts.gen_new_fasta(nooverlap, seqs, out_ref, conflicts, out_seq)
        self.assertEqual(out_ref.getvalue(),
                         ("1\t1\t1\t1:A\tNC_007795.1\n"
                          "1\t1\t2\t1:GT\tNC_007795.1\n"
                          "2\t2\t1\tAll\tNC_007795.1\n"))
        self.assertTrue(os.path.exists(os.path.join(
            self.test_folder, "seq_NC_007795.1_1_1.fa")))
        self.assertTrue(os.path.exists(os.path.join(
            self.test_folder, "seq_NC_007795.1_1_2.fa")))
        self.assertTrue(os.path.exists(os.path.join(
            self.test_folder, "seq_NC_007795.1_2_1.fa")))

    def test_snp_detect(self):
        depth_file = os.path.join(self.test_folder, "depth")
        gen_file(depth_file, self.example.depth_file)
        fasta_file = os.path.join(self.test_folder, "NC_007795.1.fa")
        gen_file(fasta_file, self.example.fasta_final)
        snp_file = os.path.join(self.test_folder, "NC_007795.1.snp")
        gen_file(snp_file, self.example.snp_final)
        out_seq = os.path.join(self.test_folder, "seq")
        out_snp = os.path.join(self.test_folder, "snp")
        stat_file = os.path.join(self.test_folder, "stat")
        args = self.mock_args.mock()
        args.depth = 5
        args.fraction = 0.3
        args.quality = 5
        args.depth_s = "n_10"
        args.depth_b = "a_2"
        args.dp4_sum = "n_10"
        args.dp4_frac = 0.5
        args.idv = "n_10"
        args.imf = 0.5
        args.filters = ["VDB_s0.1"]
        args.min_sample = 2
        ts.snp_detect(fasta_file, snp_file, depth_file, out_snp, out_seq,
                      2, stat_file, args, 2)
        self.assertTrue(os.path.exists(os.path.join(
            self.test_folder, "seq_NC_007795.1_1_1.fa")))
        self.assertTrue(os.path.exists(os.path.join(
            self.test_folder, "seq_NC_007795.1_1_2.fa")))
        self.assertTrue(os.path.exists(os.path.join(
            self.test_folder, "seq_NC_007795.1_2_1.fa")))
        self.assertTrue(os.path.exists(os.path.join(
            self.test_folder, "snp_seq_reference.csv")))
        self.assertTrue(os.path.exists(os.path.join(
            self.test_folder, "snp_best.vcf")))
        self.assertTrue(os.path.exists(os.path.join(
            self.test_folder, "snp_NC_007795.1_SNP_QUAL_best.png")))
        self.assertTrue(os.path.exists(os.path.join(
            self.test_folder, "snp_NC_007795.1_SNP_QUAL_raw.png")))


class Example(object):

    fasta_file = """>NC_007795.1
AAATATATCAGCACCGTAGACGATAGAGTAGTAC"""
    fasta_final = """>NC_007795.1
CAGTACCCTCAGCACCGTAGACGATAGAGTAGTAC"""
    snp_file = """#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	BAM
NC_007795.1	1	.	C	A	98	.	DP=89;DP4=0,0,60,9;VDB=8.46526e-15	GT:PL:DP	1/1:125,184,0:87
NC_007795.1	6	.	A	AA	26.9515	.	INDEL;IDV=22;IMF=0.536585;DP=41;VDB=9.36323e-14;DP4=0,0,40,0	GT:PL:DP	0/1:60,0,55:40"""
    snp_final = """#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	BAM
NC_007795.1	1	.	CA	A,GT	98	.	DP=89;DP4=0,0,60,9;VDB=8.46526e-15;DP4=0,0,40,0	GT:PL:DP	1/1:125,184,0:87
NC_007795.1	2	.	A	AA	26.9515	.	INDEL;IDV=22;IMF=0.536585;DP=41;VDB=9.36323e-14;DP4=0,0,40,0	GT:PL:DP	0/1:60,0,55:40
NC_007795.1	7	.	C	A	98	.	DP=89;VDB=8.46526e-15;DP4=0,0,40,0	GT:PL:DP	1/1:125,184,0:87"""
    stat = """NC_007795.1:
the number of QUAL which is between 0 and 10 = 0
the number of QUAL which is between 10 and 20 = 0
the number of QUAL which is between 20 and 30 = 0
the number of QUAL which is between 30 and 40 = 0
the number of QUAL which is between 40 and 50 = 0
the number of QUAL which is between 50 and 60 = 0
the number of QUAL which is between 60 and 70 = 0
the number of QUAL which is between 70 and 80 = 0
the number of QUAL which is between 80 and 90 = 0
the number of QUAL which is between 90 and 100 = 1
the number of QUAL which is between 100 and 110 = 0
the total numbers of QUAL are 1"""

    depth_file = """aaa	1	100
aaa	2	100
aaa	3	50
aaa	4	200"""

if __name__ == "__main__":
    unittest.main()
