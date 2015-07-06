import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data, extract_info
import annogesiclib.fill_gap as fg

class Mock_gff_parser(object):

    def __init__(self):
        self.example = Example()

    def entries(self, fh):
        fh.close()
        for element in self.example.longs:
            yield element

class TestFillGap(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_uni(self):
        out = StringIO()
        fg.uni(self.example.tas, self.example.gffs, out)
        datas, attributes = extract_info(out.getvalue(), "string")
        refs, attributes_ref = extract_info(self.example.out_uni, "string")
        self.assertEqual(set(datas), set(refs))

    def test_overlap(self):
        out = StringIO()
        print_list = []
        fg.overlap(self.example.tas, self.example.gffs, print_list, out)
        datas, attributes = extract_info(out.getvalue(), "string")
        refs, attributes_ref = extract_info(self.example.out_overlap, "string")
        self.assertEqual(set(datas), set(refs))

    def test_longer_ta(self):
        fg.Gff3Parser = Mock_gff_parser
        out_file = os.path.join(self.test_folder, "test.out")
        ta_file = os.path.join(self.test_folder, "test.ta")
        with open(ta_file, "w") as fh:
            fh.write("ta")
        fg.longer_ta(ta_file, 30, out_file)
        datas, attributes = extract_info(out_file, "file")
        refs, attributes_ref = extract_info(self.example.out_long, "string")
        self.assertEqual(set(datas), set(refs))

class Example(object):
    ta_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "Transcript", "start": 140,
                "end": 367, "phase": ".", "strand": "+", "score": "."},
               {"seq_id": "aaa", "source": "Refseq", "feature": "Transcript", "start": 230,
                "end": 240, "phase": ".", "strand": "+", "score": "."},
               {"seq_id": "bbb", "source": "Refseq", "feature": "Transcript", "start": 430,
                "end": 5167, "phase": ".", "strand": "-", "score": "."},
               {"seq_id": "bbb", "source": "Refseq", "feature": "Transcript", "start": 5269,
                "end": 6767, "phase": ".", "strand": "-", "score": "."},
               {"seq_id": "bbb", "source": "Refseq", "feature": "Transcript", "start": 7890,
                "end": 8767, "phase": ".", "strand": "-", "score": "."}]
    ta_long = [{"seq_id": "aaa", "source": "Refseq", "feature": "Transcript", "start": 140,
                "end": 367, "phase": ".", "strand": "+", "score": "."},
               {"seq_id": "aaa", "source": "Refseq", "feature": "Transcript", "start": 230,
                "end": 300, "phase": ".", "strand": "+", "score": "."},
               {"seq_id": "bbb", "source": "Refseq", "feature": "Transcript", "start": 430,
                "end": 440, "phase": ".", "strand": "-", "score": "."},
               {"seq_id": "bbb", "source": "Refseq", "feature": "Transcript", "start": 5269,
                "end": 6767, "phase": ".", "strand": "-", "score": "."},
               {"seq_id": "bbb", "source": "Refseq", "feature": "Transcript", "start": 5890,
                "end": 8767, "phase": ".", "strand": "-", "score": "."}]
    attributes_tas = [{"ID": "tran0", "Name": "Transcript_0"},
                      {"ID": "tran1", "Name": "Transcript_1"},
                      {"ID": "tran2", "Name": "Transcript_2"},
                      {"ID": "tran3", "Name": "Transcript_3"},
                      {"ID": "tran4", "Name": "Transcript_4"}]
    gff_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "gene", "start": 150,
                 "end": 200, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 1230,
                 "end": 1240, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "bbb", "source": "Refseq", "feature": "CDS", "start": 100,
                 "end": 6167, "phase": ".", "strand": "-", "score": "."},
                {"seq_id": "bbb", "source": "Refseq", "feature": "CDS", "start": 3100,
                 "end": 6167, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "bbb", "source": "Refseq", "feature": "CDS", "start": 7100,
                 "end": 9167, "phase": ".", "strand": "-", "score": "."}]
    attributes_gff = [{"ID": "gene0", "Name": "Gene_0", "locus_tag": "AAA_00001"},
                      {"ID": "cds1", "Name": "CDS_1", "locus_tag": "AAA_00002"},
                      {"ID": "cds2", "Name": "CDS_2", "locus_tag": "BBB_00001"},
                      {"ID": "cds3", "Name": "CDS_3", "locus_tag": "BBB_00002"},
                      {"ID": "cds4", "Name": "CDS_4", "locus_tag": "BBB_00003"}] 
    gffs = []
    tas = []
    longs = []
    for index in range(0, 5):
        tas.append(Create_generator(ta_dict[index], attributes_tas[index], "gff"))
        gffs.append(Create_generator(gff_dict[index], attributes_gff[index], "gff"))
        longs.append(Create_generator(ta_long[index], attributes_tas[index], "gff"))
    out_overlap = """bbb	Refseq	Transcript	430	6767	.	-	.	ID=tran3;Name=Transcript_3
bbb	Refseq	Transcript	7890	8767	.	-	.	ID=tran4;Name=Transcript_4"""
    out_uni = """aaa	Refseq	Transcript	140	367	.	+	.	Name=Transcript_0;ID=tran0
aaa	Refseq	Transcript	230	240	.	+	.	Name=Transcript_1;ID=tran1"""
    out_long = """##gff-version 3
aaa	Refseq	Transcript	140	367	.	+	.	ID=tran0;Name=Transcript_00000
bbb	Refseq	Transcript	5269	8767	.	-	.	ID=tran1;Name=Transcript_00001"""
if __name__ == "__main__":
    unittest.main()

