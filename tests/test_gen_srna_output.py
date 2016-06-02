import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
import annogesiclib.gen_srna_output as gso
from mock_args_container import MockClass


class Mock_func(object):

    def __init__(self):
        self.example = Example()

    def mock_read_gff(self, srna_file):
        return self.example.srnas

class TestGensRNAOutput(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.mock_args = MockClass()
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
        args = self.mock_args.mock()
        args.min_len = 30
        args.max_len = 500
        new_final = gso.compare_srna_table(self.example.srna_tables, srna, final, args)
        self.assertDictEqual(new_final, {'end_pro': 'NA', 'strand': '+', 'strain': 'aaa',
                                         'avg': 100, 'type': 'TEX+/-;Fragmented',
                                         'conds': 'tex_frag', 'candidates': '300-367',
                                         'tss_pro': 'TSS:300_+', 'start': 300, 'utr': '3UTR',
                                         'energy': -23, 'end': 367})

    def test_compare(self):
        args = self.mock_args.mock()
        args.min_len = 30
        args.max_len = 500
        finals = gso.compare(self.example.srnas, self.example.srna_tables,
                             self.example.nr_blasts, self.example.srna_blasts, args)
        for index in range(len(finals)):
            self.assertDictEqual(finals[index], self.example.finals[index])

    def test_gen_best_srna(self):
        gso.read_gff = Mock_func().mock_read_gff
        args = self.mock_args.mock()
        args.min_len = 30
        args.max_len = 500
        args.nr_hits_num = 0
        args.energy = 0
        args.import_info = ["tss", "sec_str"]
        args.all_hit = True
        args.best_sorf = True
        args.best_promoter = True
        args.best_term = True
        out_file = os.path.join(self.test_folder, "test.out")
        gso.gen_best_srna("test.srna", out_file, args)
        with open(out_file) as fh:
            for line in fh:
                if not (line.startswith("#")):
                    data = "\t".join(line.split("\t")[:-1])
        self.assertEqual(data, "aaa\tUTR_derived\tsRNA\t300\t367\t.\t+\t.")

class Example(object):
    srna_dict = [{"seq_id": "aaa", "source": "UTR_derived", "feature": "sRNA", "start": 300,
                  "end": 367, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "UTR_derived", "feature": "sRNA", "start": 500,
                  "end": 640, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "bbb", "source": "intergenic", "feature": "sRNA", "start": 18,
                  "end": 50, "phase": ".", "strand": "-", "score": "."}]
    attributes_srna = [{"ID": "srna0", "Name": "sRNA_0", "2d_energy": "-34",
                        "nr_hit": "2", "sRNA_hit": "1", "sRNA_type": "5utr", "sORF": "NA", "with_term": "NA"},
                       {"ID": "srna1", "Name": "sRNA_1", "2d_energy": "-23",
                        "nr_hit": "NA", "sRNA_hit": "NA", "sRNA_type": "3utr", "sORF": "sORF_1",
                        "with_term": "terminator:366-378_+"},
                       {"ID": "srna2", "Name": "sRNA_2", "2d_energy": "-11", "sRNA_type": "intergenic",
                        "nr_hit": "1", "sRNA_hit": "NA", "sORF": "NA", "with_term": "terminator:5-19_-"}]
    srnas = []
    for index in range(0, 3):
        srnas.append(Create_generator(srna_dict[index], attributes_srna[index], "gff"))

    srna_tables = [{"strain": "aaa", "strand": "+", "start": 300, "end": 367, "tss_pro": "TSS:300_+", "end_pro": "NA", "conds": "tex_frag", "avg": 100},
                   {"strain": "aaa", "strand": "+", "start": 500, "end": 640, "tss_pro": "TSS:538_+;Cleavage:640_+", "end_pro": "Cleavage:1040_+", "conds": "tex", "avg": 30},
                   {"strain": "bbb", "strand": "-", "start": 18, "end": 50, "tss_pro": "Cleavage:50_-;Cleavage:18_-", "end_pro": "NA", "conds": "frag", "avg": 300}]
    nr_blasts = [{"strain": "aaa", "strand": "+", "start": 300, "end": 367, "hits": "111;222"},
                 {"strain": "aaa", "strand": "+", "start": 500, "end": 640, "hits": "NA"},
                 {"strain": "bbb", "strand": "-", "start": 18, "end": 50, "hits": "333"}]
    srna_blasts = [{"strain": "aaa", "strand": "+", "start": 300, "end": 367, "hits": "444"},
                   {"strain": "aaa", "strand": "+", "start": 500, "end": 640, "hits": "NA"},
                   {"strain": "bbb", "strand": "-", "start": 18, "end": 50, "hits": "NA"}]
    finals = [{'sRNA_hit': '444', 'end': 367, 'with_term': 'NA', 'candidates': '300-367',
               'start': 300, 'energy': '-34', 'conds': 'tex_frag', 'end_pro': 'NA', 'promoter': 'NA',
               'tss_pro': 'TSS:300_+', 'nr_hit_num': '2', 'type': 'TEX+/-;Fragmented',
               'utr': "5'UTR_derived", 'strain': 'aaa', 'nr_hit': '111;222', 'strand': '+',
               'score': 100, 'sRNA_hit_num': '1', 'sORF': 'NA', 'avg': 100},
              {'sRNA_hit': 'NA', 'end': 640, 'with_term': 'terminator:366-378_+',
               'candidates': '500-640;538-640', 'start': 500, 'energy': '-23', 'conds': 'tex',
               'end_pro': 'Cleavage:1040_+', 'promoter': 'NA', 'tss_pro': 'TSS:538_+;Cleavage:640_+',
               'nr_hit_num': 'NA', 'type': 'TEX+/-', 'utr': "3'UTR_derived", 'strain': 'aaa',
               'nr_hit': 'NA', 'strand': '+', 'score': 30, 'sRNA_hit_num': 'NA', 'sORF': 'sORF_1', 'avg': 30},
              {'sRNA_hit': 'NA', 'end': 50, 'with_term': 'terminator:5-19_-', 'candidates': '18-50',
               'start': 18, 'energy': '-11', 'conds': 'frag', 'end_pro': 'NA', 'promoter': 'NA',
               'tss_pro': 'Cleavage:50_-;Cleavage:18_-', 'nr_hit_num': '1', 'type': 'Fragmented',
               'utr': 'Intergenic', 'strain': 'bbb', 'nr_hit': '333', 'strand': '-', 'score': 300,
               'sRNA_hit_num': 'NA', 'sORF': 'NA', 'avg': 300}]

if __name__ == "__main__":
    unittest.main()

