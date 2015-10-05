import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import extract_info
import annogesiclib.merge_sRNA as ms

class Mock_func(object):

    def __init__(self):
        self.example = Example()

    def mock_read_gff(self, gff_file):
        if gff_file == "UTR":
            return self.example.srnas_utr
        elif gff_file == "inter":
            return self.example.srnas_int


class TestMergesRNA(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_modify_attributes(self):
        pre_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 3,
                    "end": 33, "phase": ".", "strand": "+", "score": "."}
        attributes_pre = {"ID": "sRNA0", "Name": "srna_0", "UTR_type": "5utr"}
        tar1_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 3,
                     "end": 33, "phase": ".", "strand": "+", "score": "."}
        attributes_tar1 = {"ID": "sRNA0", "Name": "srna_0", "UTR_type": "3utr"}
        tar2_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 3,
                     "end": 33, "phase": ".", "strand": "+", "score": "."}
        attributes_tar2 = {"ID": "sRNA0", "Name": "srna_0", "UTR_type": "5utr"}
        pre = Create_generator(pre_dict, attributes_pre, "gff")
        tar1 = Create_generator(tar1_dict, attributes_tar1, "gff")
        ms.modify_attributes(pre, tar1, "UTR", "pre")
        self.assertEqual(pre.attributes["UTR_type"], "3utr&5utr")
        pre = Create_generator(pre_dict, attributes_pre, "gff")
        tar2 = Create_generator(tar2_dict, attributes_tar2, "gff")
        ms.modify_attributes(pre, tar2, "UTR", "pre")
        self.assertEqual(pre.attributes["UTR_type"], "5utr")
        pre = Create_generator(pre_dict, attributes_pre, "gff")
        tar1 = Create_generator(tar1_dict, attributes_tar1, "gff")
        ms.modify_attributes(pre, tar1, "UTR", "current")
        self.assertEqual(pre.attributes["UTR_type"], "5utr")
        self.assertEqual(tar1.attributes["UTR_type"], "3utr&5utr")

    def test_detect_overlap(self):
        pre_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 3,
                    "end": 33, "phase": ".", "strand": "+", "score": "."}
        attributes_pre = {"ID": "sRNA0", "Name": "srna_0", "UTR_type": "5utr"}
        tar1_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 3,
                     "end": 33, "phase": ".", "strand": "+", "score": "."}
        attributes_tar1 = {"ID": "sRNA0", "Name": "srna_0", "UTR_type": "3utr"}
        tar2_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 53,
                     "end": 233, "phase": ".", "strand": "+", "score": "."}
        attributes_tar2 = {"ID": "sRNA0", "Name": "srna_0", "UTR_type": "5utr"}
        pre = Create_generator(pre_dict, attributes_pre, "gff")
        tar1 = Create_generator(tar1_dict, attributes_tar1, "gff")
        tar2 = Create_generator(tar2_dict, attributes_tar2, "gff")
        overlap = False
        overlap = ms.detect_overlap(tar1, pre, "UTR", overlap)
        self.assertTrue(overlap)
        overlap = False
        overlap = ms.detect_overlap(tar2, pre, "UTR", overlap)
        self.assertFalse(overlap)

    def test_modify_overlap(self):
        pre_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 3,
                    "end": 33, "phase": ".", "strand": "+", "score": "."}
        attributes_pre = {"ID": "sRNA0", "Name": "srna_0", "UTR_type": "5utr", "with_TSS": "NA",
                          "start_cleavage": "cleavage_1&cleavage_2", "end_cleavage": "NA"}
        tar_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 5,
                    "end": 30, "phase": ".", "strand": "+", "score": "."}
        attributes_tar = {"ID": "sRNA0", "Name": "srna_0", "UTR_type": "3utr", "with_TSS": "TSS_1",
                          "start_cleavage": "cleavage3", "end_cleavage": "cleavage10"}
        pre = Create_generator(pre_dict, attributes_pre, "gff")
        tar = Create_generator(tar_dict, attributes_tar, "gff")
        pre_srna = ms.modify_overlap(pre, tar)
        self.assertEqual(pre_srna.attributes["with_TSS"], "TSS_1")
        self.assertEqual(pre_srna.attributes["start_cleavage"], "cleavage_1&cleavage_2&cleavage3")
        self.assertEqual(pre_srna.attributes["end_cleavage"], "cleavage10")
        self.assertEqual(pre_srna.start, 3)
        self.assertEqual(pre_srna.end, 33)

    def test_merge_srna(self):
        srnas = ms.merge_srna(self.example.srnas_utr, "UTR")
        self.assertEqual(len(srnas), 2)
        self.assertEqual(srnas[0].start, 3)
        self.assertEqual(srnas[1].start, 54)
        self.assertEqual(srnas[0].attributes["with_TSS"], "TSS_1")
        self.assertEqual(srnas[1].attributes["with_TSS"], "TSS_3")
        self.assertEqual(srnas[0].attributes["start_cleavage"], "cleavage_1&cleavage_2&cleavage_3")
        self.assertEqual(srnas[1].attributes["start_cleavage"], "cleavage_4")
        srnas = ms.merge_srna(self.example.srnas_int, "inter")
        self.assertEqual(srnas[0].attributes["with_TSS"], "TSS_1")
        self.assertEqual(srnas[1].attributes["with_TSS"], "NA")

    def test_merge_srna_gff(self):
        out_file = os.path.join(self.test_folder, "test_out")
        ms.read_gff = Mock_func().mock_read_gff
        ms.merge_srna_gff("UTR", "inter", out_file)
        datas, attributes = extract_info(out_file, "file")
        self.assertListEqual(datas, ['aaa\tUTR_derived\tsRNA\t3\t33\t.\t+\t.',
                                     'aaa\tintergenic\tsRNA\t3\t33\t.\t+\t.',
                                     'aaa\tUTR_derived\tsRNA\t54\t254\t.\t+\t.',
                                     'aaa\tintergenic\tsRNA\t54\t254\t.\t+\t.'])
        

    def test_compare_table(self):
        texs = {"track_tex_track_notex": 0}
        replicates = {"tex": 1, "frag": 1}
        wigs = {"aaa": {"frag_1": {"track_1": [{"pos": 1, "coverage": 100, "type": "frag"},
                                               {"pos": 2, "coverage": 30, "type": "frag"},
                                               {"pos": 3, "coverage": 23, "type": "frag"},
                                               {"pos": 4, "coverage": 21, "type": "frag"},
                                               {"pos": 5, "coverage": 21, "type": "frag"}]}}}
        tables = [{"strain": "aaa", "name": "sRNA_1",
                   "start": 3, "end": 4, "strand": "+",
                   "libs": "track_1", "detect": "True", "avg": 30,
                   "high": 100, "low": 20,
                   "detail": "detail"}]
        srna_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 3,
                     "end": 4, "phase": ".", "strand": "+", "score": "."}
        attributes_srna = {"ID": "sRNA0", "Name": "srna_0", "UTR_type": "3utr", "with_TSS": "TSS_1",
                           "start_cleavage": "cleavage3", "end_cleavage": "cleavage10"}
        srna = Create_generator(srna_dict, attributes_srna, "gff")
        out = StringIO()
        ms.compare_table(srna, tables, "utr", wigs, wigs, texs,
                         True, out, 2, replicates)
        self.assertEqual(out.getvalue(), "aaa\tsrna_0\t3\t4\t+\ttrack_1\tTrue\tTSS_1;cleavage3\tcleavage10\t30\t100\t20\tdetail\n")

    def test_get_coverage(self):
        wigs = {"aaa": {"frag_1": {"track_1": [{"strand": "+", "pos": 1, "coverage": 100, "type": "frag"},
                                               {"strand": "+", "pos": 2, "coverage": 30, "type": "frag"},
                                               {"strand": "+", "pos": 3, "coverage": 23, "type": "frag"},
                                               {"strand": "+", "pos": 4, "coverage": 21, "type": "frag"},
                                               {"strand": "+", "pos": 5, "coverage": 21, "type": "frag"}]}}}
        srna_cover = ms.get_coverage(wigs, "aaa", "+", 3, 4)
        self.assertEqual(srna_cover["frag_1"], [{'low': 21, 'track': 'track_1', 'type': 'frag', 'avg': 22.0, 'high': 23, 'pos': 0}])

    def test_get_tss_pro(self):
        srna_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 3,
                     "end": 4, "phase": ".", "strand": "+", "score": "."}
        attributes_srna = {"ID": "sRNA0", "Name": "srna_0", "UTR_type": "3utr", "with_TSS": "TSS_1",
                           "start_cleavage": "cleavage3", "end_cleavage": "cleavage10"}
        srna = Create_generator(srna_dict, attributes_srna, "gff")
        tss_pro = ms.get_tss_pro("utr", srna)
        self.assertEqual(tss_pro, ('TSS_1;cleavage3', 'cleavage10'))


class Example(object):
    srna_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 3,
                  "end": 33, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 5,
                  "end": 30, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 54,
                  "end": 254, "phase": ".", "strand": "+", "score": "."}]
    attributes_utr = [{"ID": "sRNA0", "Name": "srna_0", "UTR_type": "5utr",
                        "with_TSS": "NA", "start_cleavage": "cleavage_1&cleavage_2", "end_cleavage": "NA",},
                       {"ID": "sRNA1", "Name": "srna_1", "UTR_type": "3utr",
                        "with_TSS": "TSS_1", "start_cleavage": "cleavage_3", "end_cleavage": "cleavage_30",},
                       {"ID": "sRNA2", "Name": "srna_2", "with_TSS": "TSS_3", 
                        "start_cleavage": "cleavage_4", "end_cleavage": "cleavage_40", "UTR_type": "interCDS"}]
    attributes_int = [{"ID": "sRNA0", "Name": "srna_0", "UTR_type": "5utr",
                        "with_TSS": "NA", "end_cleavage": "cleavage_10"},
                       {"ID": "sRNA1", "Name": "srna_1", "UTR_type": "3utr",
                        "with_TSS": "TSS_1", "end_cleavage": "NA"},
                       {"ID": "sRNA2", "Name": "srna_2", "with_TSS": "NA", "end_cleavage": "NA"}]
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
    srnas_int = []
    srnas_utr = []
    genes = []
    for index in range(0, 3):
        srnas_int.append(Create_generator(srna_dict[index], attributes_int[index], "gff"))
        srnas_utr.append(Create_generator(srna_dict[index], attributes_utr[index], "gff"))
        genes.append(Create_generator(gff_dict[index], attributes_gff[index], "gff"))

if __name__ == "__main__":
    unittest.main()

