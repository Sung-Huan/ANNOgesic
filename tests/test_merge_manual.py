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
        tss_entry = [tss.attribute_string, {"UTR_length": "Primary_25&Internal_NA", "type": "Primary&Internal",
                                            "associated_gene": "AAA_00001&AAA_00004"}]
        tss_output = mm.remove_primary(tss, tss_entry)
        self.assertEqual(tss_output[0], 'UTR_length=Internal_NA;associated_gene=AAA_00004;type=Internal;Name=TSS_3+')
        self.assertDictEqual(tss_output[1], {'associated_gene': 'AAA_00004', 'type': 'Internal',
                                             'Name': 'TSS_3+', 'UTR_length': 'Internal_NA'})

    def test_import_to_tss(self):
        tss_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 3,
                    "end": 3, "phase": ".", "strand": "+", "score": "."}
        attributes_tss = {"ID": "CDS0", "Name": "CDS_0", "type": "Primary&Internal",
                          "associated_gene": "AAA_00001&AAA_00004", "UTR_length": "Primary_25&Internal_NA"}
        tss = Create_generator(tss_dict, attributes_tss, "gff")
        tss_entry = [tss.attribute_string, {"UTR_length": "Primary_25&Internal_NA", "type": "Primary&Internal",
                                            "associated_gene": "AAA_00001&AAA_00004"}]
        tss_type = "Primary"
        cds_pos = 10
        locus_tag = "AAA_00001"
        output = mm.import_to_tss(tss_type, cds_pos, tss, locus_tag, tss_entry)
        self.assertEqual(output[0], 'UTR_length=Primary_7&Internal_NA;associated_gene=AAA_00001&AAA_00004;type=Primary&Internal;Name=TSS_3+')
        self.assertDictEqual(output[1], {'Name': 'TSS_3+', 'UTR_length': 'Primary_7&Internal_NA',
                                         'type': 'Primary&Internal', 'associated_gene': 'AAA_00001&AAA_00004'})

    def test_same_strand_tss_gene(self):
        tss_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 3,
                    "end": 3, "phase": ".", "strand": "+", "score": "."}
        attributes_tss = {"ID": "TSS0", "Name": "TSS_0", "type": "Primary&Internal",
                          "associated_gene": "AAA_00001&AAA_00004", "UTR_length": "Primary_25&Internal_NA"}
        tss = Create_generator(tss_dict, attributes_tss, "gff")
        tss_entry = [tss.attribute_string, {"UTR_length": "Primary_25&Internal_NA", "type": "Primary&Internal",
                                            "associated_gene": "AAA_00001&AAA_00004"}]
        anti_ends = {"forward": 1, "reverse": -1}
        gene_ends = {"forward": -1, "reverse": -1}
        gff_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 6,
                    "end": 12, "phase": ".", "strand": "+", "score": "."}
        attributes_gff = {"ID": "CDS0", "Name": "CDS_0", "locus_tag": "AAA_00001"}
        gene = Create_generator(gff_dict, attributes_gff, "gff")
        checks = {"orphan": False, "int_anti": False}
        output = mm.same_strand_tss_gene(gene, tss, anti_ends, gene_ends, checks, tss_entry)
        self.assertEqual(output[0], 'UTR_length=Primary_3&Internal_NA;associated_gene=AAA_00001&AAA_00004;type=Primary&Internal;Name=TSS_3+')
        self.assertDictEqual(output[1], {'Name': 'TSS_3+', 'UTR_length': 'Primary_3&Internal_NA',
                                         'type': 'Primary&Internal', 'associated_gene': 'AAA_00001&AAA_00004'})

    def test_diff_strand_tss_gene(self):
        tss_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 3,
                    "end": 3, "phase": ".", "strand": "+", "score": "."}
        attributes_tss = {"ID": "TSS0", "Name": "TSS_0", "type": "Primary&Internal",
                          "associated_gene": "AAA_00001&AAA_00004", "UTR_length": "Primary_25&Internal_NA"}
        tss = Create_generator(tss_dict, attributes_tss, "gff")
        tss_entry = [tss.attribute_string, {"UTR_length": "Primary_25", "type": "Primary",
                                            "associated_gene": "AAA_00001"}]
        anti_ends = {"forward": 1, "reverse": -1}
        gene_ends = {"forward": -1, "reverse": -1}
        gff_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 6,
                    "end": 12, "phase": ".", "strand": "-", "score": "."}
        attributes_gff = {"ID": "CDS0", "Name": "CDS_0", "locus_tag": "AAA_00005"}
        gene = Create_generator(gff_dict, attributes_gff, "gff")
        checks = {"orphan": False, "int_anti": False}
        output = mm.diff_strand_tss_gene(gene, tss, anti_ends, gene_ends, checks, tss_entry)
        self.assertEqual(output[0], 'UTR_length=Primary_25&Antisense_NA;associated_gene=AAA_00001&AAA_00005;type=Primary&Antisense;Name=TSS_3+')
        self.assertDictEqual(output[1], {'Name': 'TSS_3+', 'UTR_length': 'Primary_25&Antisense_NA',
                                         'type': 'Primary&Antisense', 'associated_gene': 'AAA_00001&AAA_00005'})

    def test_compare_tss_gene(self):
        tss_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 3,
                    "end": 3, "phase": ".", "strand": "+", "score": "."}
        attributes_tss = {"ID": "TSS0", "Name": "TSS_0", "type": "Primary&Internal",
                          "associated_gene": "AAA_00001&AAA_00004", "UTR_length": "Primary_25&Internal_NA"}
        tss = Create_generator(tss_dict, attributes_tss, "gff")
        output = mm.compare_tss_gene(tss, self.example.genes)
        self.assertEqual(output[0], 'UTR_length=Primary_3;associated_gene=AAA_00001;type=Primary;Name=TSS_3+')
        self.assertDictEqual(output[1], {'Name': 'TSS_3+', 'UTR_length': 'Primary_3',
                                         'type': 'Primary', 'associated_gene': 'AAA_00001'})

    def test_check_overlap(self):
        tss_m_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 7,
                      "end": 7, "phase": ".", "strand": "+", "score": "."}
        attributes_tss_m = {"ID": "TSS0", "Name": "TSS_0", "type": "Primary&Internal",
                            "associated_gene": "AAA_00001&AAA_00004", "UTR_length": "Primary_25&Internal_NA"}
        tss_m = Create_generator(tss_m_dict, attributes_tss_m, "gff")
        tss_p_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 7,
                      "end": 7, "phase": ".", "strand": "+", "score": "."}
        attributes_tss_p = {"ID": "TSS0", "Name": "TSS_0", "type": "Primary&Internal",
                            "associated_gene": "AAA_00001&AAA_00004", "UTR_length": "Primary_25&Internal_NA"}
        tss_p = Create_generator(tss_p_dict, attributes_tss_p, "gff")
        tss_pre_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 3,
                        "end": 3, "phase": ".", "strand": "+", "score": "."}
        attributes_tss_pre = {"ID": "TSS0", "Name": "TSS_0", "type": "Primary&Internal",
                              "associated_gene": "AAA_00001&AAA_00004", "UTR_length": "Primary_25&Internal_NA"}
        tss_pre = Create_generator(tss_pre_dict, attributes_tss_pre, "gff")
        nums = {"tss_p": 0, "tss_m": 0, "tss": 0}
        tsss = {"tsss_p":[], "tsss_m": [], "merge": []}
        num_strain = {"aaa": {"overlap": 0, "tsspredator": 0, "manual": 0}}
        overlap_num = 0
        output = mm.check_overlap(True, tss_pre, nums, False, num_strain, overlap_num,
                                  tss_m, tss_p, tsss, 1000, self.example.genes, self.example.genes)
        self.assertEqual(output, (False, 3, 1))
        output = mm.check_overlap(False, tss_pre, nums, 100, num_strain, overlap_num,
                                  tss_m, tss_p, tsss, 1000, self.example.genes, self.example.genes)
        self.assertEqual(output, (False, 1000, 0))

    def test_intersection(self):
        nums = {"tss_p": 0, "tss_m": 0, "tss": 0}
        tsss = {"tsss_m": [], "tsss_p": [], "merge": []}
        for tss in self.example.tsss:
            tss.attributes["print"] = False
        for tss in self.example.tsss2:
            tss.attributes["print"] = False
        tsss["tsss_m"] = self.example.tsss
        tsss["tsss_p"] = self.example.tsss2
        overlap_num, num_strain = mm.intersection(tsss, 3, nums, 1000,
                                  self.example.genes, self.example.genes)
        self.assertEqual(overlap_num, 2)
        self.assertDictEqual(num_strain, {'aaa': {'tsspredator': 1, 'overlap': 2, 'manual': 1}})

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
    tss2_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 3,
                  "end": 3, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 18,
                  "end": 18, "phase": ".", "strand": "-", "score": "."},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 23,
                  "end": 23, "phase": ".", "strand": "+", "score": "."}]
    attributes_tss2 = [{"ID": "CDS0", "Name": "CDS_0", "type": "Primary", "associated_gene": "AAA_00001", "UTR_length": "Primary_25"},
                       {"ID": "CDS1", "Name": "CDS_1", "type": "Internal", "associated_gene": "AAA_00002", "UTR_length": "Internal_NA"},
                       {"ID": "CDS2", "Name": "CDS_2", "type": "Primary&Antisense",
                        "associated_gene": "AAA_00004&AAA_00006", "UTR_length": "Primary_25&Internal_NA"}]
    gff_dict = [{"start": 6, "end": 15, "phase": ".",
                 "strand": "+", "seq_id": "aaa", "score": ".",
                 "source": "Refseq", "feature": "gene"},
                {"start": 1258, "end": 2234, "phase": ".",
                 "strand": "+", "seq_id": "aaa", "score": ".",
                 "source": "Refseq", "feature": "gene"},
                {"start": 3544, "end": 6517, "phase": ".",
                 "strand": "-", "seq_id": "aaa", "score": ".",
                 "source": "Refseq", "feature": "gene"}]
    attributes_gff = [{"ID": "gene0", "Name": "gene_0", "locus_tag": "AAA_00001"},
                      {"ID": "gene0", "Name": "gene_1", "locus_tag": "AAA_00002"},
                      {"ID": "gene1", "Name": "gene_2", "locus_tag": "AAA_00003"}]
    tsss = []
    tsss2 = []
    genes = []
    for index in range(0, 3):
        tsss.append(Create_generator(tss_dict[index], attributes_tss[index], "gff"))
        tsss2.append(Create_generator(tss2_dict[index], attributes_tss2[index], "gff"))
        genes.append(Create_generator(gff_dict[index], attributes_gff[index], "gff"))

if __name__ == "__main__":
    unittest.main()

