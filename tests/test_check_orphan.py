import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
import annogesiclib.check_orphan as co

class Mock_Parser_wig(object):

    def __init__(self):
        self.wigs = [{"pos": 1, "coverage": 300,
                      "strand": "+", "strain": "aaa", "track": "track_1"},
                     {"pos": 2, "coverage": 200,
                      "strand": "+", "strain": "aaa", "track": "track_1"},
                     {"pos": 1, "coverage": 500,
                      "strand": "+", "strain": "bbb", "track": "track_1"}]

    def parser(self, filename, strand):
        for wig in self.wigs:
            self.pos = wig["pos"]
            self.coverage = wig["coverage"]
            self.track = wig["track"]
            self.strand = wig["strand"]
            self.strain = wig["strain"]
            yield self

class TestCheckOrphan(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        co.WigParser = Mock_Parser_wig
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_is_primary(self):
        cds_start = 100
        cds_end = 200
        tss_pos = 90
        self.assertTrue(co.is_primary(cds_start, cds_end, tss_pos, "+"))
        tss_pos = 300
        self.assertFalse(co.is_primary(cds_start, cds_end, tss_pos, "+"))
        self.assertTrue(co.is_primary(cds_start, cds_end, tss_pos, "-"))
        
    def test_is_internal(self):
        cds_start = 100
        cds_end = 200
        tss_pos = 150
        self.assertTrue(co.is_internal(cds_start, cds_end, tss_pos, "+"))
        tss_pos = 300
        self.assertFalse(co.is_internal(cds_start, cds_end, tss_pos, "+"))

    def test_is_antisense(self):
        cds_start = 100
        cds_end = 200
        tss_pos = 50
        self.assertTrue(co.is_antisense(cds_start, cds_end, tss_pos, "+"))
        tss_pos = 150
        self.assertTrue(co.is_antisense(cds_start, cds_end, tss_pos, "+"))
        tss_pos = 210
        self.assertTrue(co.is_antisense(cds_start, cds_end, tss_pos, "+"))
        tss_pos = 500
        self.assertFalse(co.is_antisense(cds_start, cds_end, tss_pos, "+"))

    def test_detect_coverage(self):
        tss = {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 2,
                "end": 2, "phase": ".", "strand": "+", "score": "."}
        ref = {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 3,
                "end": 3, "phase": ".", "strand": "+", "score": "."}
        attributes_tss = {"type": "Primary", "ID": "tss0", "Name": "TSS:2_+"}
        attributes_ref = {"type": "Primary", "ID": "tss1", "Name": "TSS:3_+"}
        tss_diff, ref_diff = co.detect_coverage(self.example.wigs_f,
                             Create_generator(tss, attributes_tss, "gff"),
                             Create_generator(ref, attributes_ref, "gff"))
        self.assertEqual(tss_diff, 100)
        self.assertEqual(ref_diff, 50)

    def test_del_repeat(self):
        tss_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 2,
                     "end": 2, "phase": ".", "strand": "+", "score": "."},
                    {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 22,
                     "end": 22, "phase": ".", "strand": "+", "score": "."},
                    {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 122,
                     "end": 122, "phase": ".", "strand": "-", "score": "."}]
        attributes = [{"type": "Primary", "ID": "tss0", "Name": "TSS:2_+", "UTR_length": "Primary_100",
                       "associated_gene": "AAA_00001"},
                      {"type": "Primary&Primary", "ID": "tss1", "Name": "TSS:22_+", "UTR_length": "Primary_20&Primary_50",
                       "associated_gene": "AAA_00004&AAA_00005"},
                       {"type": "Secondary&Internal", "ID": "tss2", "Name": "TSS:122_-",
                        "UTR_length": "Secondary_220&Internal_NA",
                        "associated_gene": "AAA_00008&AAA_00009"}]
        tsss = []
        for index in range(0, 3):
            tsss.append(Create_generator(tss_dict[index], attributes[index], "gff"))
        co.del_repeat(tsss)
        utrs = []
        for tss in tsss:
            utrs.append(tss.attributes["UTR_length"])
        self.assertEqual(set(utrs), set(["Primary_100", "Primary_20", "Internal_NA&Secondary_220"]))

    def test_fix_attributes(self):
        tss_entry = {"locus": "AAA_00003"}
        tss = Create_generator(self.example.tss_dict, self.example.attributes, "gff")
        co.fix_attributes(tss, tss_entry)
        self.assertEqual(tss.attributes["type"], "Primary&Secondary")

    def test_get_primary_locus_tag(self):
        tss = Create_generator(self.example.tss_dict, self.example.attributes, "gff")
        tsss = co.get_primary_locus_tag(tss)
        self.assertListEqual(tsss, [{'utr': 100, 'locus': 'AAA_00001', 'type': 'Primary'},
                                   {'utr': 120, 'locus': 'AAA_00003', 'type': 'Primary'}])

    def test_fix_primary_type(self):
        tss_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 2,
                     "end": 2, "phase": ".", "strand": "+", "score": "."},
                    {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 3,
                     "end": 3, "phase": ".", "strand": "+", "score": "."},
                    {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 4,
                     "end": 4, "phase": ".", "strand": "+", "score": "."}]
        attributes = [{"type": "Primary", "ID": "tss0", "Name": "TSS:2_+", "UTR_length": "Primary_10",
                       "associated_gene": "AAA_00001"},
                      {"type": "Primary&Internal", "ID": "tss1", "Name": "TSS:3_+", "UTR_length": "Primary_20&Internal_NA",
                       "associated_gene": "AAA_00001&AAA_00005"},
                      {"type": "Primary&Primary", "ID": "tss2", "Name": "TSS:4_+", "UTR_length": "Primary_40&Primary_60",
                       "associated_gene": "AAA_00001&AAA_00004"}]
        tsss = []
        for index in range(0, 3):
            tsss.append(Create_generator(tss_dict[index], attributes[index], "gff"))
        new_tsss = co.fix_primary_type(tsss, self.example.wigs_f, self.example.wigs_r)
        utrs = []
        for tss in new_tsss:
            utrs.append(tss.attributes["UTR_length"])
        self.assertEqual(set(utrs), set(["Internal_NA&Secondary_20", "Primary_60", "Primary_10"]))

    def test_read_wig(self):
        wig_file = os.path.join(self.test_folder, "test.wig")
        wig_fh = open(wig_file, "w")
        wig_fh.close()
        wigs = co.read_wig(wig_file, "+")
        self.assertDictEqual(wigs, {'aaa': {'track_1': [{'pos': 1, 'coverage': 300, 'strand': '+'},
                                                        {'pos': 2, 'coverage': 200, 'strand': '+'}]},
                                    'bbb': {'track_1': [{'pos': 1, 'coverage': 500, 'strand': '+'}]}})

    def test_compare_cds_check_orphan(self):
        tss_dict = [{"start": 517, "end": 517, "phase": ".",
                     "strand": "+", "seq_id": "aaa", "score": ".",
                     "source": "Refseq", "feature": "TSS"},
                    {"start": 1234, "end": 1234, "phase": ".",
                     "strand": "+", "seq_id": "aaa", "score": ".",
                     "source": "Refseq", "feature": "TSS"},
                    {"start": 3444, "end": 3444, "phase": ".",
                     "strand": "-", "seq_id": "aaa", "score": ".",
                     "source": "Refseq", "feature": "TSS"}]
        gff_dict = [{"start": 600, "end": 1517, "phase": ".",
                     "strand": "+", "seq_id": "aaa", "score": ".",
                     "source": "Refseq", "feature": "CDS"},
                    {"start": 1258, "end": 2234, "phase": ".",
                     "strand": "+", "seq_id": "aaa", "score": ".",
                     "source": "Refseq", "feature": "tRNA"},
                    {"start": 3544, "end": 6517, "phase": ".",
                     "strand": "-", "seq_id": "aaa", "score": ".",
                     "source": "Refseq", "feature": "CDS"}]
        attributes_tss = [{"type": "Orphan", "ID": "tss0", "Name": "TSS:517_+",
                           "UTR_length": "Orphan_NA", "associated_gene": "orphan"},
                          {"type": "Orphan", "ID": "tss1", "Name": "TSS:1234_+",
                           "UTR_length": "Orphan_NA", "associated_gene": "orphan"},
                          {"type": "Secondary", "ID": "tss2", "Name": "TSS:3444_-",
                           "UTR_length": "Secondary_220", "associated_gene": "AAA_00003"}]
        attributes_gff = [{"ID": "cds0", "Name": "CDS_0", "locus_tag": "AAA_00001"},
                          {"ID": "rna0", "Name": "tRNA_0", "locus_tag": "AAA_00002"},
                          {"ID": "cds1", "Name": "CDS_1", "locus_tag": "AAA_00003"}]
        out_tss = [{'Name': 'TSS:517_+', 'ID': 'tss0', 'associated_gene': 'AAA_00001',
                    'type': 'Primary', 'UTR_length': 'Primary_83'},
                   {'Name': 'TSS:1234_+', 'ID': 'tss1', 'associated_gene': 'AAA_00001&AAA_00002',
                    'type': 'Internal&Primary', 'UTR_length': 'Internal_NA&Primary_24'},
                   {'Name': 'TSS:3444_-', 'ID': 'tss2', 'associated_gene': 'AAA_00003',
                    'type': 'Secondary', 'UTR_length': 'Secondary_220'}]
        tsss = []
        cdss = []
        for index in range(0, 3):
            tsss.append(Create_generator(tss_dict[index], attributes_tss[index], "gff"))
            cdss.append(Create_generator(gff_dict[index], attributes_gff[index], "gff"))
        co.compare_cds_check_orphan(tsss, cdss)
        for index in range(0, 3):
            self.assertDictEqual(tsss[index].attributes, out_tss[index]) 

class Example(object):

    tss_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 2,
                "end": 2, "phase": ".", "strand": "+", "score": "."}
    attributes = {"type": "Primary&Primary", "ID": "tss0", "Name": "TSS:2_+", "UTR_length": "Primary_100&Primary_120",
                  "associated_gene": "AAA_00001&AAA_00003"}

    wigs_f = {"aaa": {"texnotex": [{"pos": 1, "coverage": 300,
                                    "strand": "+", "type": "tex"},
                                   {"pos": 2, "coverage": 400,
                                    "strand": "+", "type": "tex"},
                                   {"pos": 3, "coverage": 450,
                                    "strand": "+", "type": "tex"},
                                   {"pos": 4, "coverage": 470,
                                    "strand": "+", "type": "tex"}]}}

    wigs_r = {"aaa": {"texnotex": [{"pos": 1, "coverage": 300,
                                    "strand": "-", "type": "tex"},
                                   {"pos": 2, "coverage": 300,
                                    "strand": "-", "type": "tex"},
                                   {"pos": 3, "coverage": 330,
                                    "strand": "-", "type": "tex"},
                                   {"pos": 4, "coverage": 350,
                                    "strand": "-", "type": "tex"}]}}

if __name__ == "__main__":
    unittest.main()

