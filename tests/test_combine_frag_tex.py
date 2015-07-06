import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import extract_info
import annogesiclib.combine_frag_tex as cft

class Mock_Gff_parser(object):

    def __init__(self):
        self.example = Example()

    def entries(self, fh):
        for line in fh:
            if "frag" in line:
                lists = self.example.frag_dict
            elif "tex" in line:
                lists = self.example.tex_dict
        for index in range(0, 3):
            yield Create_generator(lists[index],
                                   self.example.attributes_gffs[index], "gff")

class TestCombineFragTex(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_compare(self):
        data1_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "transcript", "start": 140,
                      "end": 367, "phase": ".", "strand": "+", "score": "."}
        data2_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "transcript", "start": 180,
                      "end": 400, "phase": ".", "strand": "+", "score": "."}
        data3_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "transcript", "start": 50,
                      "end": 138, "phase": ".", "strand": "+", "score": "."}
        data4_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "transcript", "start": 5650,
                      "end": 7100, "phase": ".", "strand": "+", "score": "."}
        attributes = {"ID": "tran0", "Name": "Tran_0", "locus_tag": "AAA_00001"}
        overlap = False
        data1 = Create_generator(data1_dict, attributes, "gff")
        data2 = Create_generator(data2_dict, attributes, "gff")
        data3 = Create_generator(data3_dict, attributes, "gff")
        data4 = Create_generator(data4_dict, attributes, "gff")
        overlap12 = cft.compare(data1, data2, overlap, 5)
        self.assertEqual(data1.start, 140)
        self.assertEqual(data1.end, 400)
        overlap13 = cft.compare(data1, data3, overlap, 5)
        self.assertEqual(data1.start, 50)
        self.assertEqual(data1.end, 400)
        overlap14 = cft.compare(data1, data4, overlap, 5)
        self.assertEqual(data1.start, 50)
        self.assertEqual(data1.end, 400)
        self.assertTrue(overlap12)
        self.assertTrue(overlap13)
        self.assertFalse(overlap14)

    def test_combine(self):
        cft.Gff3Parser = Mock_Gff_parser
        output_file = os.path.join(self.test_folder, "test.out")
        frag_file = os.path.join(self.test_folder, "frag.gff")
        tex_file = os.path.join(self.test_folder, "tex.gff")
        with open(frag_file, "w") as fh:
            fh.write("frag")
        with open(tex_file, "w") as th:
            th.write("tex")
        cft.combine(frag_file, tex_file, 5, output_file)
        trans = []
        outs, attributes_out = extract_info(output_file, "file")
        refs, attributes_ref = extract_info(self.example.out_tran, "string")
        self.assertEqual(set(outs), set(refs))

class Example(object):

    frag_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "Transcript", "start": 140,
                  "end": 367, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "Transcript", "start": 30,
                  "end": 40, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "bbb", "source": "Refseq", "feature": "Transcript", "start": 430,
                  "end": 567, "phase": ".", "strand": "-", "score": "."}]
    tex_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "Transcript", "start": 120,
                 "end": 367, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "Transcript", "start": 3,
                 "end": 38, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "bbb", "source": "Refseq", "feature": "Transcript", "start": 1430,
                 "end": 1667, "phase": ".", "strand": "-", "score": "."}]
    attributes_gffs = [{"ID": "CDS0", "Name": "CDS_0", "locus_tag": "AAA_00001"},
                       {"ID": "CDS1", "Name": "CDS_1", "locus_tag": "AAA_00002"},
                       {"ID": "CDS2", "Name": "CDS_2", "locus_tag": "BBB_00001"}]

    out_tran = """##gff-version 3
aaa	fragmented_and_normal	Transcript	3	40	.	+	.	Name=Tran_00000;ID=tran0
aaa	fragmented_and_normal	Transcript	120	367	.	+	.	Name=Tran_00001;ID=tran1
bbb	fragmented	Transcript	430	567	.	-	.	Name=Tran_00002;ID=tran2
bbb	normal	Transcript	1430	1667	.	-	.	Name=Tran_00003;ID=tran3"""

if __name__ == "__main__":
    unittest.main()

