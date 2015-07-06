import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import extract_info, read_dict, gen_file
import annogesiclib.combine_gff as c_gff

class Mock_Gff_parser(object):

    def __init__(self):
        self.example = Example()

    def entries(self, fh):
        for line in fh:
            if "gff" in line:
                lists = self.example.gff_dict
                attributes = self.example.attributes_gff
                num = 3
            elif "tran" in line:
                lists = self.example.tran_dict
                attributes = self.example.attributes_tran
                num = 3
            elif "term" in line:
                lists = self.example.term_dict
                attributes = self.example.attributes_term
                num = 3
            elif "tss" in line:
                lists = self.example.tss_dict
                attributes = self.example.attributes_tss
                num = 3
            elif "utr5" in line:
                lists = self.example.utr5_dict
                attributes = self.example.attributes_utr5
                num = 2
            elif "utr3" in line:
                lists = self.example.utr3_dict
                attributes = self.example.attributes_utr3
                num = 2
        for index in range(0, num):
            yield Create_generator(lists[index], attributes[index], "gff")

class TestCombineGff(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_del_attributes(self):
        gff_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 140,
                   "end": 367, "phase": ".", "strand": "+", "score": "."}
        attributes = {"ID": "CDS0", "Name": "CDS_0", "locus_tag": "AAA_00001", "Parent_tran": "tran0"}
        gff = Create_generator(gff_dict, attributes, "gff")
        c_gff.del_attributes(gff)
        self.assertDictEqual(gff.attributes, {'ID': 'CDS0', 'Name': 'CDS_0', 'locus_tag': 'AAA_00001'})
        
    def test_compare_tran(self):
        tran_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "Transcript", "start": 100,
                     "end": 500, "phase": ".", "strand": "+", "score": "."}
        attributes_tran = {"ID": "tran0", "Name": "Tran_0"}
        out = StringIO()
        gffs = read_dict(3, self.example.gff_dict, self.example.attributes_gff)
        tran = Create_generator(tran_dict, attributes_tran, "gff")
        c_gff.compare_tran(gffs, tran, out)
        datas, attributes = extract_info(out.getvalue(), "string")
        parents = []
        for attribute in attributes:
            for element in attribute:
                if "Parent_tran" in element:
                    parents.append(element)
        self.assertEqual(set(datas), set(["aaa\tRefseq\tCDS\t160\t300\t.\t+\t."]))
        self.assertEqual(set(parents), set(["Parent_tran=tran0"]))
        out.close()

    def test_compare_tran_term(self):
        trans = read_dict(3, self.example.tran_dict, self.example.attributes_tran)
        terms = read_dict(3, self.example.term_dict, self.example.attributes_term)
        out = StringIO()
        for tran in trans:
            for term in terms:
                c_gff.compare_tran_term(term, tran, out, 3)
        datas, attributes = extract_info(out.getvalue(), "string")
        parents = []
        for attribute in attributes:
            for element in attribute:
                if "Parent_tran" in element:
                    parents.append(element)
        self.assertEqual(set(datas), set(["aaa\tRefseq\tTerminator\t350\t367\t.\t+\t.",
                                          "bbb\tRefseq\tTerminator\t420\t429\t.\t-\t."]))
        self.assertEqual(set(parents), set(["Parent_tran=tran0", "Parent_tran=tran2"]))
        out.close()

    def test_combine_gff(self):
        c_gff.Gff3Parser = Mock_Gff_parser
        gff_file = os.path.join(self.test_folder, "gff")
        ta_file = os.path.join(self.test_folder, "tran")
        tss_file = os.path.join(self.test_folder, "tss")
        term_file = os.path.join(self.test_folder, "term")
        utr3_file = os.path.join(self.test_folder, "utr3")
        utr5_file = os.path.join(self.test_folder, "utr5")
        ref_file = os.path.join(self.test_folder, "ref")
        gen_file(gff_file, "gff")
        gen_file(ta_file, "tran")
        gen_file(tss_file, "tss")
        gen_file(term_file, "term")
        gen_file(utr3_file, "utr3")
        gen_file(utr5_file, "utr5")
        gen_file(ref_file, self.example.out_file)
        out_file = os.path.join(self.test_folder, "test.out")
        c_gff.combine_gff(gff_file, ta_file, tss_file, utr5_file, utr3_file,
                          term_file, 5, 5, out_file)
        datas, attributes = extract_info(out_file, "file")
        refs, attributes = extract_info(ref_file, "file")
        self.assertEqual(set(datas), set(refs))

class Example(object):

    tran_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "Transcript", "start": 140,
                  "end": 367, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "Transcript", "start": 30,
                  "end": 40, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "bbb", "source": "Refseq", "feature": "Transcript", "start": 430,
                  "end": 567, "phase": ".", "strand": "-", "score": "."}]
    gff_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 160,
                 "end": 300, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 3,
                 "end": 38, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "bbb", "source": "Refseq", "feature": "CDS", "start": 420,
                 "end": 577, "phase": ".", "strand": "-", "score": "."}]
    term_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "Terminator", "start": 350,
                  "end": 367, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "bbb", "source": "Refseq", "feature": "Terminator", "start": 420,
                  "end": 429, "phase": ".", "strand": "-", "score": "."},
                 {"seq_id": "bbb", "source": "Refseq", "feature": "Terminator", "start": 1420,
                  "end": 2429, "phase": ".", "strand": "-", "score": "."}]
    tss_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 138,
                 "end": 138, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 330,
                 "end": 330, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "bbb", "source": "Refseq", "feature": "TSS", "start": 568,
                 "end": 568, "phase": ".", "strand": "-", "score": "."}]
    utr5_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "5UTR", "start": 140,
                  "end": 160, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "bbb", "source": "Refseq", "feature": "5UTR", "start": 500,
                  "end": 567, "phase": ".", "strand": "-", "score": "."}]
    utr3_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "3UTR", "start": 300,
                  "end": 367, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "3UTR", "start": 38,
                  "end": 40, "phase": ".", "strand": "+", "score": "."}]
    attributes_gff = [{"ID": "CDS0", "Name": "CDS_0", "locus_tag": "AAA_00001"},
                      {"ID": "CDS1", "Name": "CDS_1", "locus_tag": "AAA_00002"},
                      {"ID": "CDS2", "Name": "CDS_2", "locus_tag": "BBB_00001"}]
    attributes_tran = [{"ID": "tran0", "Name": "Tran_0"},
                       {"ID": "tran1", "Name": "Tran_1"},
                       {"ID": "tran2", "Name": "Tran_2"}]
    attributes_term = [{"ID": "term0", "Name": "Term_0"},
                       {"ID": "term1", "Name": "Term_1"},
                       {"ID": "term2", "Name": "Term_2"}]
    attributes_tss = [{"ID": "tss0", "Name": "TSS_0"},
                      {"ID": "tss1", "Name": "TSS_1"},
                      {"ID": "tss2", "Name": "TSS_2"}]
    attributes_utr5 = [{"ID": "utr5_0", "Name": "UTR5_0"},
                       {"ID": "utr5_1", "Name": "UTR5_1"}]
    attributes_utr3 = [{"ID": "utr3_0", "Name": "UTR3_0"},
                       {"ID": "utr3_1", "Name": "UTR3_1"}]
    out_file = """aaa	Refseq	Transcript	30	40	.	+	.	ID=tran1;Name=Tran_1
aaa	Refseq	3UTR	38	40	.	+	.	ID=utr3_1;Name=UTR3_1;Parent_tran=tran1
aaa	Refseq	Transcript	140	367	.	+	.	ID=tran0;Name=Tran_0
aaa	Refseq	TSS	138	138	.	+	.	ID=tss0;Name=TSS_0;Parent_tran=tran0
aaa	Refseq	TSS	330	330	.	+	.	ID=tss1;Name=TSS_1;Parent_tran=tran0
aaa	Refseq	5UTR	140	160	.	+	.	ID=utr5_0;Name=UTR5_0;Parent_tran=tran0
aaa	Refseq	CDS	160	300	.	+	.	Name=CDS_0;locus_tag=AAA_00001;ID=CDS0;Parent_tran=tran0
aaa	Refseq	3UTR	300	367	.	+	.	ID=utr3_0;Name=UTR3_0;Parent_tran=tran0
aaa	Refseq	Terminator	350	367	.	+	.	ID=term0;Name=Term_0;Parent_tran=tran0
bbb	Refseq	Transcript	430	567	.	-	.	ID=tran2;Name=Tran_2
bbb	Refseq	TSS	568	568	.	-	.	ID=tss2;Name=TSS_2;Parent_tran=tran2
bbb	Refseq	5UTR	500	567	.	-	.	ID=utr5_1;Name=UTR5_1;Parent_tran=tran2
bbb	Refseq	Terminator	420	429	.	-	.	ID=term1;Name=Term_1;Parent_tran=tran2
aaa	Refseq	CDS	3	38	.	+	.	Name=CDS_1;locus_tag=AAA_00002;ID=CDS1
bbb	Refseq	CDS	420	577	.	-	.	Name=CDS_2;locus_tag=BBB_00001;ID=CDS2
bbb	Refseq	Terminator	1420	2429	.	-	.	ID=term2;Name=Term_2"""

if __name__ == "__main__":
    unittest.main()

