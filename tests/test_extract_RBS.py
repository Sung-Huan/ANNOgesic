import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
import annogesiclib.extract_RBS as er
from mock_args_container import MockClass


class Mock_Helper(object):

    def __init__(self):
        self.example = Example()

    def extract_gene(seq, start, end, strand):
        return seq[(int(start)-1):int(end)]

class TestExtractRBS(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        self.mock_args = MockClass()
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_detect_site(self):
        inters = [{"seq": "ATGGTGACCCAGGAGGTTGATCCCAGACGTAGGACCTGTTT"},
                  {"seq": "TTAGGACGTACTCCTCGAATGATCAACTGATACTTA"},
                  {"seq": "TTTTTTTTTAAAAAAAAAATATATATTTTTTTTTTT"}]
        args = self.mock_args.mock()
        args.start_codons = ["ATG"]
        args.end_rbs = 14
        args.start_rbs = 5
        args.fuzzy_rbs = 2
        args.without_rbs = False
        args.rbs_seq = ["AGGAGG"]
        ribos = er.detect_site(inters, args)
        self.assertListEqual(ribos, [
            {'seq': 'ATGGTGACCCAGGAGGTTGATCCCAGACGTAGGACCTGTTT'},
            {'seq': 'TTAGGACGTACTCCTCGAATGATCAACTGATACTTA'}])

    def test_extract_seq(self):
        er.helper = Mock_Helper
        inters = er.extract_seq(self.example.gffs, self.example.seq,
                                self.example.tsss, self.example.tas, 5, 300)
        self.assertDictEqual(inters[0], {
            'protein': 'AAA_00001', 'strain': 'aaa', 'start': 2,
            'seq': 'AAAATTAT', 'end': 3, 'strand': '+'})
        self.assertDictEqual(inters[1], {
            'protein': 'AAA_00001', 'strain': 'aaa', 'start': 1,
            'seq': 'AAAATTAT', 'end': 3, 'strand': '+'})


class Example(object):

    gff_dict = [
        {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 3,
         "end": 30, "phase": ".", "strand": "+", "score": "."},
        {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 14,
         "end": 35, "phase": ".", "strand": "-", "score": "."},
        {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 37,
         "end": 55, "phase": ".", "strand": "-", "score": "."},
        {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 40,
         "end": 66, "phase": ".", "strand": "+", "score": "."},
        {"seq_id": "bbb", "source": "Refseq", "feature": "CDS", "start": 4,
         "end": 25, "phase": ".", "strand": "-", "score": "."}]
    attributes_gff = [
        {"ID": "cds0", "Name": "CDS_0", "locus_tag": "AAA_00001"},
        {"ID": "cds1", "Name": "CDS_1", "locus_tag": "AAA_00002",
         "protein_id": "YP_500332.1"},
        {"ID": "cds2", "Name": "CDS_2"},
        {"ID": "cds3", "Name": "CDS_3", "locus_tag": "AAA_00003"},
        {"ID": "cds4", "Name": "CDS_4", "locus_tag": "BBB_00001"}]
    ta_dict = [{"seq_id": "aaa", "source": "Refseq",
                "feature": "Transcript", "start": 1,
                "end": 367, "phase": ".", "strand": "+", "score": "."},
               {"seq_id": "aaa", "source": "Refseq",
                "feature": "Transcript", "start": 230,
                "end": 240, "phase": ".", "strand": "+", "score": "."},
               {"seq_id": "bbb", "source": "Refseq",
                "feature": "Transcript", "start": 430,
                "end": 5167, "phase": ".", "strand": "-", "score": "."}]
    attributes_tas = [
        {"ID": "tran0", "Name": "Transcript_0", "locus_tag": "AAA_00001"},
        {"ID": "tran1", "Name": "Transcript_1", "locus_tag": "AAA_00002"},
        {"ID": "tran2", "Name": "Transcript_2", "locus_tag": "BBB_00001"}]
    tss_dict = [{"seq_id": "aaa", "source": "Refseq",
                 "feature": "TSS", "start": 2,
                 "end": 2, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq",
                 "feature": "TSS", "start": 230,
                 "end": 230, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "bbb", "source": "Refseq",
                 "feature": "TSS", "start": 5166,
                 "end": 5166, "phase": ".", "strand": "-", "score": "."}]
    attributes_tss = [{"ID": "tss0", "Name": "TSS_0", "type": "Primary",
                       "associated_gene": "AAA_00001"},
                      {"ID": "tss1", "Name": "TSS_1", "type": "Internal",
                       "associated_gene": "AAA_00002"},
                      {"ID": "tss2", "Name": "TSS_2", "type": "Orphan",
                       "associated_gene": "orphan"}]
    gffs = []
    tas = []
    tsss = []
    for index in range(0, 3):
        gffs.append(Create_generator(gff_dict[index],
                                     attributes_gff[index], "gff"))
        tas.append(Create_generator(ta_dict[index],
                                    attributes_tas[index], "gff"))
        tsss.append(Create_generator(tss_dict[index],
                                     attributes_tss[index], "gff"))
    seq = {"aaa": "AAAATTATAGGCGTAGTAACCTCTTGATAGCGATGGATATAGACCCTTATAAGGCCTCTGATTAGAAAATAGGTAGGCCCCCGGGGGTGTGTAATAGATAGAT",
           "bbb": "ATATGTACCCCGCGCCGTATAGCTATAAATTCGCTGCTTATTTTATA"}

if __name__ == "__main__":
    unittest.main()

