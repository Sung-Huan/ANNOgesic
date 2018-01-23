import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import extract_info, gen_file
import annogesiclib.merge_sRNA as ms
from mock_args_container import MockClass


class Mock_func(object):

    def __init__(self):
        self.example = Example()

    def mock_read_gff(self, gff_file, type_, ex_srna):
        if gff_file == "UTR":
            return self.example.srnas_utr
        elif gff_file == "inter":
            return self.example.srnas_int
        elif gff_file == "both":
            return self.example.srnas_utr
        elif type_ == "CDS":
            return self.example.genes

    def mock_replicate_comparison(
        self, srna_covers, template_texs, strand, cutoff_coverage, args,
        type_, median, coverages, utr_type, notex):
        return {'low': 21, 'start': 3, 'conds': {'frag_1': '1'}, 'end': 4,
                'best': 22.0, 'high': 23, 'track': 'track_1',
                'detail': [{'low': 21, 'pos': 0, 'final_end': 4,
                            'final_start': 3, 'type': 'frag',
                            'track': 'track_1', 'avg': 22.0, 'high': 23}]}


class TestMergesRNA(unittest.TestCase):

    def setUp(self):
        self.mock_args = MockClass()
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_modify_attributes(self):
        pre_dict = {"seq_id": "aaa", "source": "Refseq",
                    "feature": "sRNA", "start": 3,
                    "end": 33, "phase": ".", "strand": "+", "score": "."}
        attributes_pre = {"ID": "sRNA0", "Name": "srna_0",
                          "sRNA_type": "5utr"}
        tar1_dict = {"seq_id": "aaa", "source": "Refseq",
                     "feature": "sRNA", "start": 3,
                     "end": 33, "phase": ".", "strand": "+", "score": "."}
        attributes_tar1 = {"ID": "sRNA0", "Name": "srna_0",
                           "sRNA_type": "antisense"}
        tar2_dict = {"seq_id": "aaa", "source": "Refseq",
                     "feature": "sRNA", "start": 3,
                     "end": 33, "phase": ".", "strand": "+", "score": "."}
        attributes_tar2 = {"ID": "sRNA0", "Name": "srna_0",
                           "sRNA_type": "5utr"}
        pre = Create_generator(pre_dict, attributes_pre, "gff")
        tar1 = Create_generator(tar1_dict, attributes_tar1, "gff")
        ms.modify_attributes(pre, tar1, "UTR", "pre")
        self.assertEqual(pre.attributes["sRNA_type"], "5utr")
        pre = Create_generator(pre_dict, attributes_pre, "gff")
        tar2 = Create_generator(tar2_dict, attributes_tar2, "gff")
        ms.modify_attributes(pre, tar2, "UTR", "pre")
        self.assertEqual(pre.attributes["sRNA_type"], "5utr")
        pre = Create_generator(pre_dict, attributes_pre, "gff")
        tar1 = Create_generator(tar1_dict, attributes_tar1, "gff")
        ms.modify_attributes(pre, tar1, "UTR", "current")
        self.assertEqual(pre.attributes["sRNA_type"], "5utr")
        self.assertEqual(tar1.attributes["sRNA_type"], "5utr")

    def test_detect_overlap(self):
        pre_dict = {"seq_id": "aaa", "source": "Refseq",
                    "feature": "sRNA", "start": 3,
                    "end": 33, "phase": ".", "strand": "+", "score": "."}
        attributes_pre = {"ID": "sRNA0", "Name": "srna_0", "sRNA_type": "5utr"}
        tar1_dict = {"seq_id": "aaa", "source": "Refseq",
                     "feature": "sRNA", "start": 3,
                     "end": 33, "phase": ".", "strand": "+", "score": "."}
        attributes_tar1 = {"ID": "sRNA0", "Name": "srna_0", "sRNA_type": "3utr"}
        tar2_dict = {"seq_id": "aaa", "source": "Refseq",
                     "feature": "sRNA", "start": 53,
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
        pre_dict = {"seq_id": "aaa", "source": "Refseq",
                    "feature": "sRNA", "start": 3,
                    "end": 33, "phase": ".", "strand": "+", "score": "."}
        attributes_pre = {"ID": "sRNA0", "Name": "srna_0",
                          "sRNA_type": "5utr", "with_TSS": "NA",
                          "start_cleavage": "cleavage_1,cleavage_2",
                          "end_cleavage": "NA"}
        tar_dict = {"seq_id": "aaa", "source": "Refseq",
                    "feature": "sRNA", "start": 5,
                    "end": 30, "phase": ".", "strand": "+", "score": "."}
        attributes_tar = {"ID": "sRNA0", "Name": "srna_0",
                          "sRNA_type": "3utr", "with_TSS": "TSS_1",
                          "start_cleavage": "cleavage3",
                          "end_cleavage": "cleavage10"}
        pre = Create_generator(pre_dict, attributes_pre, "gff")
        tar = Create_generator(tar_dict, attributes_tar, "gff")
        pre_srna = ms.modify_overlap(pre, tar)
        self.assertEqual(pre_srna.attributes["with_TSS"], "TSS_1")
        self.assertEqual(pre_srna.attributes["start_cleavage"],
                         "cleavage_1,cleavage_2,cleavage3")
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
        self.assertEqual(srnas[0].attributes["start_cleavage"],
                         "cleavage_1,cleavage_2,cleavage_3")
        self.assertEqual(srnas[1].attributes["start_cleavage"], "cleavage_4")
        srnas = ms.merge_srna(self.example.srnas_int, "inter")
        self.assertEqual(srnas[0].attributes["with_TSS"], "TSS_1")
        self.assertEqual(srnas[1].attributes["with_TSS"], "NA")

    def test_merge_srna_gff(self):
        out_file = os.path.join(self.test_folder, "test_out")
        gen_file(os.path.join(self.test_folder, "aaa.gff"),
                 self.example.gff_file)
        ms.read_gff = Mock_func().mock_read_gff
        gffs = {"merge": out_file, "utr": "UTR", "normal": "inter"}
        ms.merge_srna_gff(gffs, False, 0.5,
                          os.path.join(self.test_folder, "aaa.gff"), False)
        datas, attributes = extract_info(out_file, "file")
        self.assertListEqual(datas,
                             ['aaa\tANNOgesic\tncRNA\t54\t254\t.\t+\t.'])
        self.assertEqual(set(attributes[0]),
                         set(['overlap_cds=NA', 'Name=sRNA_00000',
                         'ID=aaa_srna0', 'sRNA_type=intergenic',
                         'end_cleavage=cleavage_40', 'with_TSS=TSS_3',
                         'overlap_percent=NA']))

    def test_compare_table(self):
        ms.replicate_comparison = Mock_func().mock_replicate_comparison
        wigs = {"aaa": {"frag_1": {"track_1|+|frag": [100, 30, 23, 21, 21]}}}
        tables = [{"strain": "aaa", "name": "sRNA_1",
                   "start": 3, "end": 4, "strand": "+",
                   "libs": "track_1", "detect": "True", "avg": 30,
                   "high": 100, "low": 20,
                   "detail": "detail"}]
        srna_dict = {"seq_id": "aaa", "source": "Refseq",
                     "feature": "sRNA", "start": 3,
                     "end": 4, "phase": ".", "strand": "+", "score": "."}
        attributes_srna = {"ID": "sRNA0", "Name": "srna_0",
                           "sRNA_type": "3utr", "with_TSS": "TSS_1",
                           "start_cleavage": "cleavage3",
                           "end_cleavage": "cleavage10",
                           "overlap_cds": "CDS1", "overlap_percent": "0.01415"}
        srna = Create_generator(srna_dict, attributes_srna, "gff")
        tss_dict = {"seq_id": "aaa", "source": "Refseq",
                    "feature": "TSS", "start": 3,
                    "end": 3, "phase": ".", "strand": "+", "score": "."}
        attributes_tss = {"ID": "tss0", "Name": "TSS_0", "type": "Orphan"}
        tss = Create_generator(tss_dict, attributes_tss, "gff")
        out = StringIO()
        cutoff_tex = [0, 0, 0, 50, 20]
        cutoff_notex = [0, 0, 0, 30, 10]
        cutoff_frag = [400, 200, 0, 50, 30]
        gen_file("tmp_median", "aaa\t3utr\ttrack_1\t10")
        args = self.mock_args.mock()
        args.replicates = replicates = {"tex": 1, "frag": 1}
        args.texs = texs = {"track_tex_track_notex": 0}
        args.out_folder = os.getcwd()
        args.table_best = True
        args.tex_notex = 2
        ms.compare_table(srna, tables, "utr", wigs, wigs, texs,
                         out, [tss], args)
        self.assertEqual(out.getvalue(),
                         "aaa\tsrna_0\t3\t4\t+\tfrag_1\t1\tTSS_1;cleavage3\tcleavage10\t22.0\ttrack_1(22.0)\tCDS1\t0.01415\n")
        os.remove("tmp_median")

    def test_get_coverage(self):
        wigs = {"aaa": {"frag_1": {"track_1|+|frag": [100, 30, 23, 21, 21]}}}
        srna_cover = ms.get_coverage(wigs, self.example.srnas_int[0])
        self.assertEqual(srna_cover["frag_1"],
                         [{'low': 21, 'track': 'track_1',
                           'avg': 1.3548387096774193, 'final_end': 33,
                           'high': 21, 'pos': 0, 'final_start': 3,
                           'type': 'frag'}])

    def test_get_tss_pro(self):
        srna_dict = {"seq_id": "aaa", "source": "Refseq",
                     "feature": "sRNA", "start": 3,
                     "end": 4, "phase": ".", "strand": "+", "score": "."}
        attributes_srna = {"ID": "sRNA0", "Name": "srna_0",
                           "sRNA_type": "3utr", "with_TSS": "TSS_1",
                           "start_cleavage": "cleavage3",
                           "end_cleavage": "cleavage10"}
        srna = Create_generator(srna_dict, attributes_srna, "gff")
        tss_pro = ms.get_tss_pro("utr", srna)
        self.assertEqual(tss_pro, ('TSS_1;cleavage3', 'cleavage10'))


class Example(object):
    gff_file =  "aaa\tRefSeq\tCDS\t1\t100\t.\t+\t.\tID=CDS_1;Name=CDS_1"
    srna_dict = [
        {"seq_id": "aaa", "source": "Refseq", "feature": "ncRNA", "start": 3,
         "end": 33, "phase": ".", "strand": "+", "score": "."},
        {"seq_id": "aaa", "source": "Refseq", "feature": "ncRNA", "start": 5,
         "end": 30, "phase": ".", "strand": "+", "score": "."},
        {"seq_id": "aaa", "source": "Refseq", "feature": "ncRNA", "start": 54,
         "end": 254, "phase": ".", "strand": "+", "score": "."}]
    attributes_utr = [
        {"ID": "sRNA0", "Name": "srna_0", "sRNA_type": "5utr",
         "with_TSS": "NA", "start_cleavage": "cleavage_1,cleavage_2",
         "end_cleavage": "NA",},
        {"ID": "sRNA1", "Name": "srna_1", "sRNA_type": "3utr",
         "with_TSS": "TSS_1", "start_cleavage": "cleavage_3",
         "end_cleavage": "cleavage_30",},
        {"ID": "sRNA2", "Name": "srna_2", "with_TSS": "TSS_3", 
         "start_cleavage": "cleavage_4", "end_cleavage": "cleavage_40",
         "sRNA_type": "interCDS"}]
    attributes_int = [
        {"ID": "sRNA0", "Name": "srna_0", "sRNA_type": "5utr",
         "with_TSS": "NA", "end_cleavage": "cleavage_10"},
        {"ID": "sRNA1", "Name": "srna_1", "sRNA_type": "3utr",
         "with_TSS": "TSS_1", "end_cleavage": "NA"},
        {"ID": "sRNA2", "Name": "srna_2", "with_TSS": "NA",
         "end_cleavage": "NA", "sRNA_type": "intergenic"}]
    gff_dict = [{"start": 6, "end": 15, "phase": ".",
                 "strand": "+", "seq_id": "aaa", "score": ".",
                 "source": "Refseq", "feature": "CDS"},
                {"start": 1258, "end": 2234, "phase": ".",
                 "strand": "+", "seq_id": "aaa", "score": ".",
                 "source": "Refseq", "feature": "CDS"},
                {"start": 3544, "end": 6517, "phase": ".",
                 "strand": "-", "seq_id": "aaa", "score": ".",
                 "source": "Refseq", "feature": "CDS"}]
    attributes_gff = [
        {"ID": "CDS0", "Name": "CDS_0", "locus_tag": "AAA_00001"},
        {"ID": "CDS0", "Name": "CDS_1", "locus_tag": "AAA_00002"},
        {"ID": "CDS1", "Name": "CDS_2", "locus_tag": "AAA_00003"}]
    srnas_int = []
    srnas_utr = []
    genes = []
    for index in range(0, 3):
        srnas_int.append(Create_generator(
            srna_dict[index], attributes_int[index], "gff"))
        srnas_utr.append(Create_generator(
            srna_dict[index], attributes_utr[index], "gff"))
        genes.append(Create_generator(
            gff_dict[index], attributes_gff[index], "gff"))
if __name__ == "__main__":
    unittest.main()
