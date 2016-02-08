import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import extract_info, gen_file
import annogesiclib.merge_sRNA as ms

class Mock_func(object):

    def __init__(self):
        self.example = Example()

    def mock_read_gff(self, gff_file, type_):
        if gff_file == "UTR":
            return self.example.srnas_utr
        elif gff_file == "inter":
            return self.example.srnas_int
        elif type_ == "CDS":
            return self.example.genes

    def mock_replicate_comparison(self, srna_covers, template_texs,
                         strand, cutoff_coverage, tex_notex, replicates,
                         type_, median, coverages, utr_type, texs, notex):
        return {'low': 21, 'start': 3, 'conds': {'frag_1': '1'}, 'end': 4,
                'best': 22.0, 'high': 23, 'track': 'track_1',
                'detail': [{'low': 21, 'pos': 0, 'final_end': 4, 'final_start': 3,
                'type': 'frag', 'track': 'track_1', 'avg': 22.0, 'high': 23}]}


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
        attributes_pre = {"ID": "sRNA0", "Name": "srna_0", "sRNA_type": "5utr"}
        tar1_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 3,
                     "end": 33, "phase": ".", "strand": "+", "score": "."}
        attributes_tar1 = {"ID": "sRNA0", "Name": "srna_0", "sRNA_type": "3utr"}
        tar2_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 3,
                     "end": 33, "phase": ".", "strand": "+", "score": "."}
        attributes_tar2 = {"ID": "sRNA0", "Name": "srna_0", "sRNA_type": "5utr"}
        pre = Create_generator(pre_dict, attributes_pre, "gff")
        tar1 = Create_generator(tar1_dict, attributes_tar1, "gff")
        ms.modify_attributes(pre, tar1, "UTR", "pre")
        self.assertEqual(pre.attributes["sRNA_type"], "3utr&5utr")
        pre = Create_generator(pre_dict, attributes_pre, "gff")
        tar2 = Create_generator(tar2_dict, attributes_tar2, "gff")
        ms.modify_attributes(pre, tar2, "UTR", "pre")
        self.assertEqual(pre.attributes["sRNA_type"], "5utr")
        pre = Create_generator(pre_dict, attributes_pre, "gff")
        tar1 = Create_generator(tar1_dict, attributes_tar1, "gff")
        ms.modify_attributes(pre, tar1, "UTR", "current")
        self.assertEqual(pre.attributes["sRNA_type"], "5utr")
        self.assertEqual(tar1.attributes["sRNA_type"], "3utr&5utr")

    def test_detect_overlap(self):
        pre_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 3,
                    "end": 33, "phase": ".", "strand": "+", "score": "."}
        attributes_pre = {"ID": "sRNA0", "Name": "srna_0", "sRNA_type": "5utr"}
        tar1_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 3,
                     "end": 33, "phase": ".", "strand": "+", "score": "."}
        attributes_tar1 = {"ID": "sRNA0", "Name": "srna_0", "sRNA_type": "3utr"}
        tar2_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 53,
                     "end": 233, "phase": ".", "strand": "+", "score": "."}
        attributes_tar2 = {"ID": "sRNA0", "Name": "srna_0", "sRNA_type": "5utr"}
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
        attributes_pre = {"ID": "sRNA0", "Name": "srna_0", "sRNA_type": "5utr", "with_TSS": "NA",
                          "start_cleavage": "cleavage_1&cleavage_2", "end_cleavage": "NA"}
        tar_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 5,
                    "end": 30, "phase": ".", "strand": "+", "score": "."}
        attributes_tar = {"ID": "sRNA0", "Name": "srna_0", "sRNA_type": "3utr", "with_TSS": "TSS_1",
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
        gen_file(os.path.join(self.test_folder, "aaa.gff"), self.example.gff_file)
        ms.read_gff = Mock_func().mock_read_gff
        ms.merge_srna_gff("UTR", "inter", out_file, False, 0.5, os.path.join(self.test_folder, "aaa.gff"))
        datas, attributes = extract_info(out_file, "file")
        self.assertListEqual(datas, ['aaa\tANNOgesic\tsRNA\t54\t254\t.\t+\t.',
                                     'aaa\tANNOgesic\tsRNA\t54\t254\t.\t+\t.'])
        self.assertEqual(set(attributes[0]), set(['overlap_percent=NA', 'end_cleavage=cleavage_40',
                                                  'start_cleavage=cleavage_4', 'Name=sRNA_00000',
                                                  'with_TSS=TSS_3', 'ID=srna0', 'sRNA_type=interCDS',
                                                  'overlap_cds=NA']))
        self.assertEqual(set(attributes[1]), set(['overlap_percent=NA', 'end_cleavage=NA', 'Name=sRNA_00001',
                                                  'with_TSS=NA', 'ID=srna1', 'sRNA_type=intergenic',
                                                  'overlap_cds=NA']))

    def test_compare_table(self):
        ms.replicate_comparison = Mock_func().mock_replicate_comparison
        texs = {"track_tex_track_notex": 0}
        replicates = {"tex": 1, "frag": 1}
        wigs = {"aaa": {"frag_1": {"track_1": [{"strand": "+", "pos": 1, "coverage": 100, "type": "frag"},
                                               {"strand": "+", "pos": 2, "coverage": 30, "type": "frag"},
                                               {"strand": "+", "pos": 3, "coverage": 23, "type": "frag"},
                                               {"strand": "+", "pos": 4, "coverage": 21, "type": "frag"},
                                               {"strand": "+", "pos": 5, "coverage": 21, "type": "frag"}]}}}
        tables = [{"strain": "aaa", "name": "sRNA_1",
                   "start": 3, "end": 4, "strand": "+",
                   "libs": "track_1", "detect": "True", "avg": 30,
                   "high": 100, "low": 20,
                   "detail": "detail"}]
        srna_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 3,
                     "end": 4, "phase": ".", "strand": "+", "score": "."}
        attributes_srna = {"ID": "sRNA0", "Name": "srna_0", "sRNA_type": "3utr", "with_TSS": "TSS_1",
                           "start_cleavage": "cleavage3", "end_cleavage": "cleavage10",
                           "overlap_cds": "CDS1", "overlap_percent": "0.01415"}
        srna = Create_generator(srna_dict, attributes_srna, "gff")
        tss_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 3,
                     "end": 3, "phase": ".", "strand": "+", "score": "."}
        attributes_tss = {"ID": "tss0", "Name": "TSS_0", "type": "Orphan"}
        tss = Create_generator(tss_dict, attributes_tss, "gff")
        out = StringIO()
        cutoff_tex = [0, 0, 0, 50, 20]
        cutoff_notex = [0, 0, 0, 30, 10]
        cutoff_frag = [400, 200, 0, 50, 30]
        gen_file("tmp_median", "aaa\t3utr\ttrack_1\t10")
        ms.compare_table(srna, tables, "utr", wigs, wigs, texs,
                         True, out, 2, replicates, cutoff_tex,
                         cutoff_notex, cutoff_frag, tss, 10, os.getcwd())
        self.assertEqual(out.getvalue(), "aaa\tsrna_0\t3\t4\t+\tfrag_1\t1\tTSS_1;cleavage3\tcleavage10\t22.0\t23\t21\ttrack_1(avg=22.0;high=23;low=21)\tCDS1\t0.01415\n")
        os.remove("tmp_median")

    def test_get_coverage(self):
        wigs = {"aaa": {"frag_1": {"track_1": [{"strand": "+", "pos": 1, "coverage": 100, "type": "frag"},
                                               {"strand": "+", "pos": 2, "coverage": 30, "type": "frag"},
                                               {"strand": "+", "pos": 3, "coverage": 23, "type": "frag"},
                                               {"strand": "+", "pos": 4, "coverage": 21, "type": "frag"},
                                               {"strand": "+", "pos": 5, "coverage": 21, "type": "frag"}]}}}
        srna_cover = ms.get_coverage(wigs, "aaa", "+", 3, 4)
        self.assertEqual(srna_cover["frag_1"], [{'avg': 22.0, 'pos': 0, 'final_end': 4, 'type': 'frag', 'track': 'track_1', 'final_start': 3, 'low': 21, 'high': 23}])

    def test_get_tss_pro(self):
        srna_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 3,
                     "end": 4, "phase": ".", "strand": "+", "score": "."}
        attributes_srna = {"ID": "sRNA0", "Name": "srna_0", "sRNA_type": "3utr", "with_TSS": "TSS_1",
                           "start_cleavage": "cleavage3", "end_cleavage": "cleavage10"}
        srna = Create_generator(srna_dict, attributes_srna, "gff")
        tss_pro = ms.get_tss_pro("utr", srna)
        self.assertEqual(tss_pro, ('TSS_1;cleavage3', 'cleavage10'))


class Example(object):
    gff_file =  "aaa\tRefSeq\tCDS\t1\t100\t.\t+\t.\tID=CDS_1;Name=CDS_1"
    srna_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 3,
                  "end": 33, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 5,
                  "end": 30, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 54,
                  "end": 254, "phase": ".", "strand": "+", "score": "."}]
    attributes_utr = [{"ID": "sRNA0", "Name": "srna_0", "sRNA_type": "5utr",
                        "with_TSS": "NA", "start_cleavage": "cleavage_1&cleavage_2", "end_cleavage": "NA",},
                       {"ID": "sRNA1", "Name": "srna_1", "sRNA_type": "3utr",
                        "with_TSS": "TSS_1", "start_cleavage": "cleavage_3", "end_cleavage": "cleavage_30",},
                       {"ID": "sRNA2", "Name": "srna_2", "with_TSS": "TSS_3", 
                        "start_cleavage": "cleavage_4", "end_cleavage": "cleavage_40", "sRNA_type": "interCDS"}]
    attributes_int = [{"ID": "sRNA0", "Name": "srna_0", "sRNA_type": "5utr",
                        "with_TSS": "NA", "end_cleavage": "cleavage_10"},
                       {"ID": "sRNA1", "Name": "srna_1", "sRNA_type": "3utr",
                        "with_TSS": "TSS_1", "end_cleavage": "NA"},
                       {"ID": "sRNA2", "Name": "srna_2", "with_TSS": "NA", "end_cleavage": "NA", "sRNA_type": "intergenic"}]
    gff_dict = [{"start": 6, "end": 15, "phase": ".",
                 "strand": "+", "seq_id": "aaa", "score": ".",
                 "source": "Refseq", "feature": "CDS"},
                {"start": 1258, "end": 2234, "phase": ".",
                 "strand": "+", "seq_id": "aaa", "score": ".",
                 "source": "Refseq", "feature": "CDS"},
                {"start": 3544, "end": 6517, "phase": ".",
                 "strand": "-", "seq_id": "aaa", "score": ".",
                 "source": "Refseq", "feature": "CDS"}]
    attributes_gff = [{"ID": "CDS0", "Name": "CDS_0", "locus_tag": "AAA_00001"},
                      {"ID": "CDS0", "Name": "CDS_1", "locus_tag": "AAA_00002"},
                      {"ID": "CDS1", "Name": "CDS_2", "locus_tag": "AAA_00003"}]
    srnas_int = []
    srnas_utr = []
    genes = []
    for index in range(0, 3):
        srnas_int.append(Create_generator(srna_dict[index], attributes_int[index], "gff"))
        srnas_utr.append(Create_generator(srna_dict[index], attributes_utr[index], "gff"))
        genes.append(Create_generator(gff_dict[index], attributes_gff[index], "gff"))

if __name__ == "__main__":
    unittest.main()

