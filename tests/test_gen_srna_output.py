import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
import annogesiclib.gen_srna_output as gso


class Mock_func(object):

    def __init__(self):
        self.example = Example()

    def mock_read_gff(self, srna_file):
        return self.example.srnas

class TestGensRNAOutput(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_merge_info(self):
        blasts = [{"strain": "aaa", "strand": "+", "start": 20, "end": 70, "hits": "111"},
                  {"strain": "aaa", "strand": "+", "start": 20, "end": 70, "hits": "222"},
                  {"strain": "aaa", "strand": "+", "start": 20, "end": 70, "hits": "333"},
                  {"strain": "aaa", "strand": "+", "start": 20, "end": 70, "hits": "444"},
                  {"strain": "bbb", "strand": "+", "start": 20, "end": 70, "hits": "555"}]
        merge = gso.merge_info(blasts)
        self.assertDictEqual(merge[0], {'hits': '111;222;333', 'start': 20, 'strand': '+', 'strain': 'aaa', 'end': 70})
        self.assertDictEqual(merge[1], {'hits': '555', 'start': 20, 'strand': '+', 'strain': 'bbb', 'end': 70})

    def test_compare_srna_table(self):
        final = {"energy": -23, "utr": "3UTR"}
        srna_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 300,
                     "end": 367, "phase": ".", "strand": "+", "score": "."}
        attributes_srna = {"ID": "srna0", "Name": "sRNA_0"}
        srna = Create_generator(srna_dict, attributes_srna, "gff")
        new_final = gso.compare_srna_table(self.example.srna_tables, srna, final, 30, 500)
        self.assertDictEqual(new_final, {'end': 367, 'start': 300, 'energy': -23, 'utr': '3UTR',
                                         'tss_pro': 'TSS:300_+', 'conds': 'tex_frag',
                                         'candidates': '300-367', 'strain': 'aaa', 'strand': '+',
                                         'type': 'TEX+/-;Fragmented'})

    def test_compare(self):
        finals = gso.compare(self.example.srnas, self.example.srna_tables,
                             self.example.nr_blasts, self.example.srna_blasts, 30, 500)
        for index in range(len(finals)):
            self.assertDictEqual(finals[index], self.example.finals[index])

    def test_gen_best_srna(self):
        gso.read_gff = Mock_func().mock_read_gff
        out_file = os.path.join(self.test_folder, "test.out")
        gso.gen_best_srna("test.srna", True, 0, 0, True, out_file)
        with open(out_file) as fh:
            for line in fh:
                if not (line.startswith("#")):
                    data = "\t".join(line.split("\t")[:-1])
        self.assertEqual(data, "aaa\tUTR_derived\tsRNA\t300\t367\t.\t+\t.")
        gso.gen_best_srna("test.srna", False, 0, 0, True, out_file)
        data = None
        with open(out_file) as fh:
            for line in fh:
                if not (line.startswith("#")):
                    data = "\t".join(line.split("\t")[:-1])
        self.assertEqual(data, None)

class Example(object):
    srna_dict = [{"seq_id": "aaa", "source": "UTR_derived", "feature": "sRNA", "start": 300,
                  "end": 367, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "UTR_derived", "feature": "sRNA", "start": 500,
                  "end": 640, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "bbb", "source": "intergenic", "feature": "sRNA", "start": 18,
                  "end": 50, "phase": ".", "strand": "-", "score": "."}]
    attributes_srna = [{"ID": "srna0", "Name": "sRNA_0", "2d_energy": "-34",
                        "nr_hit": "2", "sRNA_hit": "1", "UTR_type": "5utr", "sORF": "NA"},
                       {"ID": "srna1", "Name": "sRNA_1", "2d_energy": "-23",
                        "nr_hit": "NA", "sRNA_hit": "NA", "UTR_type": "3utr", "sORF": "sORF_1"},
                       {"ID": "srna2", "Name": "sRNA_2", "2d_energy": "-11",
                        "nr_hit": "1", "sRNA_hit": "NA", "sORF": "NA"}]
    srnas = []
    for index in range(0, 3):
        srnas.append(Create_generator(srna_dict[index], attributes_srna[index], "gff"))

    srna_tables = [{"strain": "aaa", "strand": "+", "start": 300, "end": 367, "tss_pro": "TSS:300_+", "conds": "tex_frag"},
                   {"strain": "aaa", "strand": "+", "start": 500, "end": 640, "tss_pro": "TSS:538_+;Cleavage:640_+", "conds": "tex"},
                   {"strain": "bbb", "strand": "-", "start": 18, "end": 50, "tss_pro": "Cleavage:50_-;Cleavage:18_-", "conds": "frag"}]
    nr_blasts = [{"strain": "aaa", "strand": "+", "start": 300, "end": 367, "hits": "111;222"},
                 {"strain": "aaa", "strand": "+", "start": 500, "end": 640, "hits": "NA"},
                 {"strain": "bbb", "strand": "-", "start": 18, "end": 50, "hits": "333"}]
    srna_blasts = [{"strain": "aaa", "strand": "+", "start": 300, "end": 367, "hits": "444"},
                   {"strain": "aaa", "strand": "+", "start": 500, "end": 640, "hits": "NA"},
                   {"strain": "bbb", "strand": "-", "start": 18, "end": 50, "hits": "NA"}]
    finals = [{'sORF': 'NA', 'nr_hit': '111;222', 'end': 367, 'conds': 'tex_frag', 'sRNA_hit': '444',
               'start': 300, 'utr': "5'UTR_derived", 'strand': '+', 'strain': 'aaa', 'sRNA_hit_num': '1',
               'type': 'TEX+/-;Fragmented', 'tss_pro': 'TSS:300_+', 'candidates': '300-367',
               'nr_hit_num': '2', 'energy': '-34', 'sORF': 'NA'},
              {'sORF': 'NA', 'nr_hit': 'NA', 'end': 640, 'conds': 'tex', 'sRNA_hit': 'NA',
               'start': 500, 'utr': "3'UTR_derived", 'strand': '+', 'strain': 'aaa', 'sRNA_hit_num': 'NA',
               'type': 'TEX+/-', 'tss_pro': 'TSS:538_+;Cleavage:640_+', 'candidates': '500-640;538-640',
               'nr_hit_num': 'NA', 'energy': '-23', 'sORF': 'sORF_1'},
              {'sORF': 'NA', 'nr_hit': '333', 'end': 50, 'conds': 'frag', 'sRNA_hit': 'NA',
               'start': 18, 'utr': 'Intergenic', 'strand': '-', 'strain': 'bbb', 'sRNA_hit_num': 'NA',
               'type': 'Fragmented', 'tss_pro': 'Cleavage:50_-;Cleavage:18_-', 'candidates': '18-50',
               'nr_hit_num': '1', 'energy': '-11', 'sORF': 'NA'}]


if __name__ == "__main__":
    unittest.main()

