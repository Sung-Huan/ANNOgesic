import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
import annogesiclib.extract_RBS as er

class Mock_Helper(object):

    def __init__(self):
        self.example = Example()

    def extract_gene(seq, start, end, strand):
        return seq[(int(start)-1):int(end)]

class TestExtractRBS(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_detect_site(self):
        inters = [{"seq": "ATGGTGACCCAGGAGGTTGATCCCAGACGTAGGACCTGTTT"},
                  {"seq": "TTAGGACGTACTCCTCGAATGATCAACTGATACTTA"},
                  {"seq": "TTTTTTTTTAAAAAAAAAATATATATTTTTTTTTTT"}]
        ribos = er.detect_site(inters)
        self.assertListEqual(ribos, [{'seq': 'TTAGGACGTACTCCTCGAATGATCAACTGATACTTA'}])

    def test_extract_seq(self):
        er.helper = Mock_Helper
        inters = er.extract_seq(self.example.gffs, self.example.seq)
        refs = [{'protein': 'AAA_00001', 'strand': '+', 'start': 1,
                 'seq': 'AAAATTATAGGCG', 'end': 13, 'strain': 'aaa'},
                {'protein': 'AAA_00002', 'strand': '-', 'start': 25,
                 'seq': 'TCCATCGCTATCA', 'end': 37, 'strain': 'aaa'},
                {'protein': 'AAA_00003', 'strand': '+', 'start': 30,
                 'seq': 'GCGATGGATATAGACCCTTAT', 'end': 50, 'strain': 'aaa'},
                {'protein': 'cds2', 'strand': '-', 'start': 45,
                 'seq': 'ATCTATCTATTACACACCCCCGGGGGCCTACCTATTTTCTAATCAGAGGCCTTATAAGG', 'end': 103, 'strain': 'aaa'},
                {'protein': 'BBB_00001', 'strand': '-', 'start': 15,
                 'seq': 'TATAAAATAAGCAGCGAATTTATAGCTATACGG', 'end': 47, 'strain': 'bbb'}]
        for index in range(len(inters)):
            self.assertDictEqual(inters[index], refs[index])


class Example(object):

    gff_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 3,
                 "end": 30, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 14,
                 "end": 35, "phase": ".", "strand": "-", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 37,
                 "end": 55, "phase": ".", "strand": "-", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 40,
                 "end": 66, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "bbb", "source": "Refseq", "feature": "CDS", "start": 4,
                 "end": 25, "phase": ".", "strand": "-", "score": "."}]
    attributes_gff = [{"ID": "cds0", "Name": "CDS_0", "locus_tag": "AAA_00001"},
                      {"ID": "cds1", "Name": "CDS_1", "locus_tag": "AAA_00002", "protein_id": "YP_500332.1"},
                      {"ID": "cds2", "Name": "CDS_2"},
                      {"ID": "cds3", "Name": "CDS_3", "locus_tag": "AAA_00003"},
                      {"ID": "cds4", "Name": "CDS_4", "locus_tag": "BBB_00001"}]
    gffs = []
    for index in range(0, 5):
        gffs.append(Create_generator(gff_dict[index], attributes_gff[index], "gff"))
    seq = {"aaa": "AAAATTATAGGCGTAGTAACCTCTTGATAGCGATGGATATAGACCCTTATAAGGCCTCTGATTAGAAAATAGGTAGGCCCCCGGGGGTGTGTAATAGATAGAT",
           "bbb": "ATATGTACCCCGCGCCGTATAGCTATAAATTCGCTGCTTATTTTATA"}

if __name__ == "__main__":
    unittest.main()

