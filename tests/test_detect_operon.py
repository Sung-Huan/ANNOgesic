import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data
import annogesiclib.detect_operon as op

class Mock_func(object):

    def __init__(self):
        self.example = Example()

    def mock_read_gff(self, ta_file, gff_file, tss_file, terminator_file):
        return self.example.tas, self.example.gffs, \
               self.example.tsss, self.example.terms

class TestDetectOperon(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_get_term_feature(self):
        features = {"num": 0, "detect": False}
        datas = []
        jumps = []
        for ta in self.example.tas:
            for term in self.example.terms:
                if ta.strand == "+":
                    jump = op.get_term_feature(ta, term, 5, features, datas,
                                                  ta.end, term.start, term.end - 5)
                else:
                    jump = op.get_term_feature(ta, term, 5, features, datas,
                                                   ta.start, term.end + 5, term.start)
                jumps.append(jump)
        self.assertEqual(set(jumps), set([False, True, False, True, True, False, False, False, False]))
        strands = []
        starts = []
        for data in datas:
            strands.append(data.strand)
            starts.append(data.start)
        self.assertEqual(set(strands), set(["+", "-"]))
        self.assertEqual(set(starts), set([360, 420]))

    def test_get_tss_feature(self):
        features = {"num": 0, "detect": False}
        datas = []
        for ta in self.example.tas:
            for tss in self.example.tsss:
                if ta.strand == "+":
                    op.get_tss_feature(ta, tss, features, 3, datas, ta.start,
                                           tss.start + 3, tss.end)
                else:
                    op.get_tss_feature(ta, tss, features, 3, datas, ta.end,
                                           tss.start, tss.end - 3)
        strands = []
        starts = []
        for data in datas:
            strands.append(data.strand)
            starts.append(data.start)
        self.assertEqual(set(strands), set(["+", "+", "+", "-"]))
        self.assertEqual(set(starts), set([140, 230, 230, 5166]))

    def test_detect_features(self):
        datas = []
        for ta in self.example.tas:
            datas.append(op.detect_features(ta, self.example.gffs, "gene", 5, 3))
        features = []
        detects = []
        for data in datas:
            features.append(data["num_feature"])
            detects.append(data["with_feature"])
        self.assertEqual(set(features), set([2, 0, 0]))
        self.assertEqual(set(detects), set([True, False, False]))

    def test_sub_operon(self):
        tss_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 140,
                     "end": 140, "phase": ".", "strand": "+", "score": "."},
                    {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 200,
                     "end": 200, "phase": ".", "strand": "+", "score": "."}]
        attributes_tss = [{"ID": "tss0", "Name": "TSS_0", "locus_tag": "AAA_00001"},
                          {"ID": "tss1", "Name": "TSS_1", "locus_tag": "BBB_00001"}]
        gff_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 540,
                     "end": 640, "phase": ".", "strand": "+", "score": "."},
                    {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 166,
                     "end": 166, "phase": ".", "strand": "-", "score": "."}]
        attributes_gff = [{"ID": "tss0", "Name": "TSS_0", "locus_tag": "AAA_00001"},
                          {"ID": "tss1", "Name": "TSS_1", "locus_tag": "BBB_00001"}] 
        tsss = {"with_feature": True, "num_feature": 2, "data_list": []}
        genes = {"data_list": []}
        for index in range(0, 2):
            genes["data_list"].append(Create_generator(gff_dict[index], attributes_gff[index], "gff"))
            tsss["data_list"].append(Create_generator(tss_dict[index], attributes_tss[index], "gff"))
        operons = op.sub_operon("+", tsss, 141, 300, genes, 30)
        self.assertDictEqual(operons[0], {'end': 199, 'strand': '+', 'start': 141})
        self.assertDictEqual(operons[1], {'end': 300, 'strand': '+', 'start': 200})

    def test_operon(self):
        op.read_gff = Mock_func().mock_read_gff
        out_file = os.path.join(self.test_folder, "test.out")
        op.operon("test_ta", "test_tss", "test_gff", "test_term", 3,
                  5, 5, out_file)
        datas = import_data(out_file)
        self.assertEqual(set(datas), set(self.example.out_file.split("\n")))

class Example(object):
    ta_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "Transcript", "start": 140,
                "end": 367, "phase": ".", "strand": "+", "score": "."},
               {"seq_id": "aaa", "source": "Refseq", "feature": "Transcript", "start": 230,
                "end": 240, "phase": ".", "strand": "+", "score": "."},
               {"seq_id": "bbb", "source": "Refseq", "feature": "Transcript", "start": 430,
                "end": 5167, "phase": ".", "strand": "-", "score": "."}]
    attributes_tas = [{"ID": "tran0", "Name": "Transcript_0", "locus_tag": "AAA_00001"},
                      {"ID": "tran1", "Name": "Transcript_1", "locus_tag": "AAA_00002"},
                      {"ID": "tran2", "Name": "Transcript_2", "locus_tag": "BBB_00001"}]
    term_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "Terminator", "start": 360,
                  "end": 367, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "Terminator", "start": 530,
                  "end": 540, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "bbb", "source": "Refseq", "feature": "Terminator", "start": 420,
                  "end": 429, "phase": ".", "strand": "-", "score": "."}]
    attributes_term = [{"ID": "term0", "Name": "Terminator_0", "locus_tag": "AAA_00001"},
                       {"ID": "term1", "Name": "Terminator_1", "locus_tag": "AAA_00002"},
                       {"ID": "term2", "Name": "Terminator_2", "locus_tag": "BBB_00001"}]
    tss_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 140,
                 "end": 140, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 230,
                 "end": 230, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "bbb", "source": "Refseq", "feature": "TSS", "start": 5166,
                 "end": 5166, "phase": ".", "strand": "-", "score": "."}]
    attributes_tss = [{"ID": "tss0", "Name": "TSS_0", "locus_tag": "AAA_00001"},
                      {"ID": "tss1", "Name": "TSS_1", "locus_tag": "AAA_00002"},
                      {"ID": "tss2", "Name": "TSS_2", "locus_tag": "BBB_00001"}]
    gff_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "gene", "start": 148,
                 "end": 360, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 148,
                 "end": 360, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "gene", "start": 220,
                 "end": 239, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "rRNA", "start": 220,
                 "end": 239, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "bbb", "source": "Refseq", "feature": "gene", "start": 400,
                 "end": 5200, "phase": ".", "strand": "-", "score": "."},
                {"seq_id": "bbb", "source": "Refseq", "feature": "CDS", "start": 400,
                 "end": 5200, "phase": ".", "strand": "-", "score": "."},
                {"seq_id": "bbb", "source": "Refseq", "feature": "gene", "start": 5400,
                 "end": 5800, "phase": ".", "strand": "-", "score": "."},
                {"seq_id": "bbb", "source": "Refseq", "feature": "CDS", "start": 5400,
                 "end": 5800, "phase": ".", "strand": "-", "score": "."}]
    attributes_gff = [{"ID": "gene0", "Name": "Gene_0", "locus_tag": "AAA_00001"},
                      {"ID": "cds0", "Name": "CSD_0", "locus_tag": "AAA_00001", "protein_id": "YP_000001"},
                      {"ID": "gene1", "Name": "Gene_1", "locus_tag": "AAA_00002"},
                      {"ID": "rna0", "Name": "rRNA_0", "locus_tag": "AAA_00002"},
                      {"ID": "gene2", "Name": "Gene_2", "locus_tag": "BBB_00001"},
                      {"ID": "cds1", "Name": "CDS_1", "locus_tag": "BBB_00001", "protein_id": "YP_000002"},
                      {"ID": "gene3", "Name": "Gene_3", "locus_tag": "BBB_00002"},
                      {"ID": "cds2", "Name": "CDS_2", "locus_tag": "BBB_00002", "protein_id": "YP_000003"}]
    out_file = """Operon_ID	strain	Operon_position	strand	Number_of_suboperon	Position_of_suboperon	Start_with_TSS	Number_of_TSS	Terminated_with_terminator	Number_of_terminator	Number_of_gene_associated_suboperon	Number_of_gene_associated_operon	Associated_genes_with_suboperon	Associated_genes_with_whole_operon
Operon1	aaa	140-367	+	0	NA	True	2	True	1	NA	2	NA	AAA_00001, AAA_00002
Operon2	aaa	230-240	+	0	NA	True	1	False	0	NA	0	NA	NA
Operon3	bbb	430-5167	-	0	NA	True	1	True	1	NA	0	NA	NA"""

    tas = []
    tsss = []
    terms = []
    for index in range(0, 3):
        tas.append(Create_generator(ta_dict[index], attributes_tas[index], "gff"))
        terms.append(Create_generator(term_dict[index], attributes_term[index], "gff"))
        tsss.append(Create_generator(tss_dict[index], attributes_tss[index], "gff"))
    gffs = []
    for index in range(0, 8):
        gffs.append(Create_generator(gff_dict[index], attributes_gff[index], "gff"))

if __name__ == "__main__":
    unittest.main()

