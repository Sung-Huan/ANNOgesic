import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
import annogesiclib.merge_manual as mm


class TestGensRNAOutput(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)


    def test_get_primary_locus_tag(self):
        out = mm.get_primary_locus_tag(self.example.tsss[-1])
        self.assertDictEqual(out[0], {"type": "Primary", "utr": 25, "locus": "AAA_00004"})

    def test_detect_coverage(self):
        wigs = {"aaa": {"track_1": [{"pos": 1, "coverage": 200},
                                    {"pos": 2, "coverage": 300},
                                    {"pos": 3, "coverage": 400},
                                    {"pos": 4, "coverage": 600},
                                    {"pos": 5, "coverage": 650}]}}
        tss_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 3,
                    "end": 3, "phase": ".", "strand": "+", "score": "."}
        attributes_tss = {"ID": "CDS0", "Name": "CDS_0", "type": "Primary",
                          "associated_gene": "AAA_00001", "UTR_length": "Primary_25"}
        ref_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 4,
                    "end": 4, "phase": ".", "strand": "+", "score": "."}
        attributes_ref = {"ID": "CDS1", "Name": "CDS_1", "type": "Primary",
                           "associated_gene": "AAA_00002", "UTR_length": "Primary_40"}
        tss = Create_generator(tss_dict, attributes_tss, "gff")
        ref = Create_generator(ref_dict, attributes_ref, "gff")
        tss_cover, ref_cover = mm.detect_coverage(wigs, tss, ref)
        self.assertEqual(tss_cover, 100)
        self.assertEqual(ref_cover, 200)

    def test_fix_attributes(self):
        tss_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 3,
                    "end": 3, "phase": ".", "strand": "+", "score": "."}
        attributes_tss = {"ID": "CDS0", "Name": "CDS_0", "type": "Primary",
                          "associated_gene": "AAA_00001", "UTR_length": "Primary_25"}
        tss = Create_generator(tss_dict, attributes_tss, "gff")
        tss_entry = {"locus": "AAA_00001"}
        mm.fix_attributes(tss, tss_entry)
        self.assertEqual(tss.attributes["type"], "Secondary")

    def test_del_repeat(self):
        tss_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 3,
                     "end": 3, "phase": ".", "strand": "+", "score": "."},
                    {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 5,
                     "end": 5, "phase": ".", "strand": "+", "score": "."}]
        attributes_tss = [{"ID": "CDS0", "Name": "CDS_0", "type": "Primary&Primary",
                           "associated_gene": "AAA_00001&AAA_00002", "UTR_length": "Primary_25&Primary_200"},
                          {"ID": "CDS1", "Name": "CDS_1", "type": "Primary&Antisense",
                           "associated_gene": "AAA_00003&AAA_00004", "UTR_length": "Primary_25&Antisense_NA"}]
        tsss = []
        for index in range(0, 2):
            tsss.append(Create_generator(tss_dict[index], attributes_tss[index], "gff"))
        mm.del_repeat(tsss)
        self.assertEqual(tsss[0].attributes["type"], "Primary")
        self.assertEqual(tsss[1].attributes["type"], "Antisense&Primary")

    def test_fix_primary_type(self):
        wigs = {"aaa": {"track_1": [{"pos": 1, "coverage": 200},
                                    {"pos": 2, "coverage": 300},
                                    {"pos": 3, "coverage": 400},
                                    {"pos": 4, "coverage": 600},
                                    {"pos": 5, "coverage": 650},
                                    {"pos": 6, "coverage": 655}]}}
        tss_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 3,
                     "end": 3, "phase": ".", "strand": "+", "score": "."},
                    {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 5,
                     "end": 5, "phase": ".", "strand": "+", "score": "."}]
        attributes_tss = [{"ID": "CDS0", "Name": "CDS_0", "type": "Primary&Primary",
                           "associated_gene": "AAA_00001&AAA_00002", "UTR_length": "Primary_25&Primary_200"},
                          {"ID": "CDS1", "Name": "CDS_1", "type": "Primary&Antisense",
                           "associated_gene": "AAA_00001&AAA_00004", "UTR_length": "Primary_27&Antisense_NA"}]
        tsss = []
        for index in range(0, 2):
            tsss.append(Create_generator(tss_dict[index], attributes_tss[index], "gff"))
        mm.fix_primary_type(tsss, wigs, "test")
        self.assertEqual(tsss[0].attributes["type"], "Primary")
        self.assertEqual(tsss[1].attributes["type"], "Antisense&Secondary")

    def test_remove_primary(self):
        tss_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 3,
                    "end": 3, "phase": ".", "strand": "+", "score": "."}
        attributes_tss = {"ID": "CDS0", "Name": "CDS_0", "type": "Primary&Internal",
                          "associated_gene": "AAA_00001&AAA_00004", "UTR_length": "Primary_25&Internal_NA"}
        tss = Create_generator(tss_dict, attributes_tss, "gff")
        tss_entry = [tss.attribute_string, {"UTR_length": "Primary25&Internal_NA", "type": "Primary&Internal",
                                            "associated_gene": "AAA_00001&AAA_00004"}]
        tss_output = mm.remove_primary(tss, tss_entry)
        self.assertEqual(tss_output[0], 'UTR_length=Internal_NA;associated_gene=AAA_00004;type=Internal;Name=TSS_3+')
        self.assertDictEqual(tss_output[1], {'associated_gene': 'AAA_00004', 'type': 'Internal',
                                             'Name': 'TSS_3+', 'UTR_length': 'Internal_NA'})

class Example(object):
    tss_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 3,
                 "end": 3, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 16,
                 "end": 16, "phase": ".", "strand": "-", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 54,
                 "end": 54, "phase": ".", "strand": "+", "score": "."}]
    attributes_tss = [{"ID": "CDS0", "Name": "CDS_0", "type": "Primary", "associated_gene": "AAA_00001", "UTR_length": "Primary_25"},
                      {"ID": "CDS1", "Name": "CDS_1", "type": "Internal", "associated_gene": "AAA_00002", "UTR_length": "Internal_NA"},
                      {"ID": "CDS2", "Name": "CDS_2", "type": "Primary&Antisense",
                       "associated_gene": "AAA_00004&AAA_00006", "UTR_length": "Primary_25&Internal_NA"}]
    tsss = []
    for index in range(0, 3):
        tsss.append(Create_generator(tss_dict[index], attributes_tss[index], "gff"))

if __name__ == "__main__":
    unittest.main()

