import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
import annogesiclib.get_inter_seq as gis


class TestGetInterSeq(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_get_inter(self):
        inters = gis.get_inter(self.example.cdss, self.example.seq, "cds")
        self.assertDictEqual(inters[0], {'parent_m': 'CDS_1', 'strain': 'aaa',
            'start': 10, 'end': 16, 'm_pos': '16-30', 'p_pos': '3-10',
                                         'print': False, 'parent_p': 'CDS_0'})

    def test_get_overlap_inters(self):
        inter1 = {'parent_m': 'AAA_00002', 'print': False, 'strain': 'aaa',
                  'end': 16, 'start': 10, 'parent_p': 'AAA_00001',
                  'm_pos': '16-30', 'p_pos': '3-10'}
        inter2 = {'parent_m': 'AAA_00003', 'print': False, 'strain': 'aaa',
                  'end': 28, 'start': 13, 'parent_p': 'AAA_00004',
                  'm_pos': '16-30', 'p_pos': '3-10'}
        merges = []
        gis.get_overlap_inters(inter1, inter2, merges, 0)
        self.assertDictEqual(merges[0], {
            'ID': 'inter_0', 'parent_m': 'AAA_00002',
            'strain': 'aaa', 'parent_p': 'AAA_00004',
            'print': False, 'start': 13, 'end': 16,
            'm_pos': '16-30', 'p_pos': '3-10'})

    def test_merge_inter(self):
        merges = gis.merge_inter(self.example.inters1, self.example.inters2)
        for index in range(len(merges)):
            self.assertEqual(merges[index], self.example.out_merges[index])

    def test_detect_confliction(self):
        inter = {'parent_m': 'AAA_00002', 'print': False, 'strain': 'aaa',
                 'end': 60, 'start': 30, 'parent_p': 'AAA_00001', "ID": "0",
                 'm_pos': '16-30', 'p_pos': '3-10'}
        merges = gis.detect_confliction(
            inter, self.example.cdss, self.example.seq)
        self.assertDictEqual(merges[0], {
            'strain': 'aaa', 'end': 60, 'print': False,
            'start': 1, 'parent_p': 'AAA_00001',
            'ID': 'inter_0', 'strand': '+',
            'parent_m': 'AAA_00002', 'm_pos': '16-30', 'p_pos': '3-10'})
        self.assertDictEqual(merges[1], {
            'strain': 'aaa', 'end': 88, 'print': False,
            'start': 30, 'parent_p': 'AAA_00001',
            'ID': 'inter_0', 'strand': '-',
            'parent_m': 'AAA_00002', 'm_pos': '16-30', 'p_pos': '3-10'})

class Example(object):
    seq = {"aaa": "AATTAGCCGTGTACTGTAGATAGCCCGCGCCTACTACCGATCAGATAGTCCAGTACATATAGCGATATTAGATCGGGTTTATATAAAA"}
    cds_dict = [
        {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 3,
         "end": 10, "phase": ".", "strand": "+", "score": "."},
        {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 16,
         "end": 30, "phase": ".", "strand": "-", "score": "."},
        {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 38,
         "end": 50, "phase": ".", "strand": "+", "score": "."},
        {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 54,
         "end": 60, "phase": ".", "strand": "+", "score": "."}]
    attributes_cds = [{"ID": "CDS0", "Name": "CDS_0",
                       "locus_tag": "AAA_00001"},
                      {"ID": "CDS1", "Name": "CDS_1",
                       "protein_id": "AAA_00002"},
                      {"ID": "CDS2", "Name": "CDS_2"},
                      {"ID": "CDS3", "Name": "CDS_3"}]    
    cdss = []
    for index in range(0, 4):
        cdss.append(Create_generator(
            cds_dict[index], attributes_cds[index], "gff"))
    inters1 = [{'parent_m': 'AAA_00002', 'print': False, 'strain': 'aaa',
                'end': 16, 'start': 10, 'parent_p': 'AAA_00001',
                'm_pos': '16-30', 'p_pos': '3-10'},
               {'parent_m': 'AAA_00012', 'print': False, 'strain': 'aaa',
                'end': 76, 'start': 50, 'parent_p': 'AAA_00007',
                'm_pos': '16-30', 'p_pos': '38-50'},
               {'parent_m': 'AAA_00102', 'print': False, 'strain': 'aaa',
                'end': 116, 'start': 100, 'parent_p': 'AAA_00301',
                'm_pos': 'NA', 'p_pos': '54-60'}]
    inters2 = [{'parent_m': 'AAA_00003', 'print': False, 'strain': 'aaa',
                'end': 28, 'start': 13, 'parent_p': 'AAA_00004',
                'm_pos': '16-30', 'p_pos': '3-10'},
               {'parent_m': 'AAA_00033', 'print': False, 'strain': 'aaa',
                'end': 70, 'start': 60, 'parent_p': 'AAA_00034',
                'm_pos': '16-30', 'p_pos': '38-50'},
               {'parent_m': 'AAA_02003', 'print': False, 'strain': 'aaa',
                'end': 328, 'start': 313, 'parent_p': 'AAA_04004',
                'm_pos': 'NA', 'p_pos': '54-60'}]
    out_merges = [{'strain': 'aaa', 'm_pos': '16-30', 'print': False,
                   'ID': 'inter_0', 'start': 13, 'end': 16,
                   'parent_p': 'AAA_00004', 'parent_m': 'AAA_00002',
                   'p_pos': '3-10'},
                  {'strain': 'aaa', 'm_pos': '16-30', 'print': False,
                   'ID': 'inter_1', 'start': 60, 'end': 70,
                   'parent_p': 'AAA_00034', 'parent_m': 'AAA_00033',
                   'p_pos': '38-50'},
                  {'strain': 'aaa', 'm_pos': 'NA', 'print': False,
                   'ID': 'inter_2', 'start': 100, 'end': 116,
                   'parent_p': 'AAA_00301', 'parent_m': 'AAA_00102',
                   'p_pos': '54-60'},
                  {'strain': 'aaa', 'm_pos': 'NA', 'print': False,
                   'ID': 'inter_3', 'start': 313, 'end': 328,
                   'parent_p': 'AAA_04004', 'parent_m': 'AAA_02003',
                   'p_pos': '54-60'}]
if __name__ == "__main__":
    unittest.main()
