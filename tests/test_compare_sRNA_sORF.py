import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import gen_file, extract_info
import annogesiclib.compare_sRNA_sORF as cro

class Mock_Gff_parser(object):

    def __init__(self):
        self.example = Example()

    def entries(self, fh):
        for line in fh:
            if "srna" in line:
                lists = self.example.srna_dict
                attributes = self.example.attributes_srna
            elif "sorf" in line:
                lists = self.example.sorf_dict
                attributes = self.example.attributes_sorf
        for index in range(0, 3):
            yield Create_generator(lists[index], attributes[index], "gff")
        fh.close()

class TestComparesRNAsORF(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_srna_sorf_comparison(self):
        cro.Gff3Parser = Mock_Gff_parser
        sRNA_file = os.path.join(self.test_folder, "sRNA.gff")
        sORF_file = os.path.join(self.test_folder, "sORF.gff")
        gen_file(sRNA_file, "srna")
        gen_file(sORF_file, "sorf")
        sRNA_out = os.path.join(self.test_folder, "sRNA.out")
        sORF_out = os.path.join(self.test_folder, "sORF.out")
        cro.srna_sorf_comparison(sRNA_file, sORF_file, sRNA_out, sORF_out)
        srnas, attribute_srnas = extract_info(sRNA_out, "file")
        refs, attribute_refs = extract_info(self.example.srna_out, "string")
        self.assertEqual(set(srnas), set(refs))
        self.assertEqual(set(attribute_srnas[2]), set(attribute_refs[2]))
        sorfs, attribute_sorfs = extract_info(sORF_out, "file")
        refs, attribute_refs = extract_info(self.example.sorf_out, "string")
        self.assertEqual(set(sorfs), set(refs))
        self.assertEqual(set(attribute_sorfs[2]), set(attribute_refs[2]))

class Example(object):

    srna_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 140,
                  "end": 367, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 30,
                  "end": 40, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "bbb", "source": "Refseq", "feature": "sRNA", "start": 430,
                  "end": 567, "phase": ".", "strand": "-", "score": "."}]
    sorf_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "sORF", "start": 160,
                  "end": 300, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "sORF", "start": 3,
                  "end": 38, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "bbb", "source": "Refseq", "feature": "sORF", "start": 420,
                  "end": 577, "phase": ".", "strand": "-", "score": "."}]
    attributes_srna = [{"ID": "srna0", "Name": "sRNA_0"},
                       {"ID": "srna1", "Name": "sRNA_1"},
                       {"ID": "srna2", "Name": "sRNA_2"}]
    attributes_sorf = [{"ID": "sorf0", "Name": "sORF_0"},
                       {"ID": "sorf1", "Name": "sORF_1"},
                       {"ID": "sorf2", "Name": "sORF_2"}]
    sorf_out = """##gff-version 3
aaa	Refseq	sORF	3	38	.	+	.	sRNA=NA;ID=sorf1;Name=sORF_1
aaa	Refseq	sORF	160	300	.	+	.	sRNA=srna0:140-367_f;ID=sorf0;Name=sORF_0
bbb	Refseq	sORF	420	577	.	-	.	sRNA=NA;ID=sorf2;Name=sORF_2"""

    srna_out = """##gff-version 3
aaa	Refseq	sRNA	30	40	.	+	.	ID=srna1;sORF=NA;Name=sRNA_1
aaa	Refseq	sRNA	140	367	.	+	.	ID=srna0;sORF=sorf0:160-300_f;Name=sRNA_0
bbb	Refseq	sRNA	430	567	.	-	.	ID=srna2;sORF=NA;Name=sRNA_2"""

if __name__ == "__main__":
    unittest.main()

