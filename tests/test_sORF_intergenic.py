import sys
import os
import unittest
import shutil
import copy
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data, gen_file
import annogesiclib.sORF_intergenic as si

class TestsORFIntergenic(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_get_type(self):
        inter = {"strain": "aaa", "strand": "+", "start": 100, "end": 149}
        si.get_type(inter, self.example.gffs)
        self.assertDictEqual(inter, {'strand': '+', 'source': '5utr', 'end': 149,
                                     'strain': 'aaa', 'start': 100})

    def test_compare_tran_cds(self):
        inters = si.compare_tran_cds(self.example.trans, self.example.gffs)
        self.assertListEqual(inters, [{'strand': '+', 'strain': 'aaa', 'start': 10, 'end': 100},
                                      {'strand': '+', 'strain': 'aaa', 'start': 1100, 'end': 1229}])

    def test_get_intergenic(self):
        out_file = os.path.join(self.test_folder, "out")
        gff_file = os.path.join(self.test_folder, "anno.gff")
        tran_file = os.path.join(self.test_folder, "tran.gff")
        gen_file(gff_file, self.example.gff_file)
        gen_file(tran_file, self.example.tran_file)
        si.get_intergenic(gff_file, tran_file, out_file, True)
        datas = import_data(out_file)
        self.assertEqual("\n".join(datas), self.example.out_file)


class Example(object):
    gff_file = """aaa	Refseq	CDS	112	623	.	+	.	ID=cds0;Name=CDS_00000
aaa	Refseq	CDS	800	1100	.	+	.	ID=trna1;Name=tRNA_00001"""
    tran_file = """aaa	Refseq	transcript	19	100	.	+	.	ID=tran0;Name=Tran_00000
aaa	Refseq	transcript	600	1800	.	+	.	ID=tran1;Name=Tran_00001"""
    out_file = """aaa	intergenic	sORF	19	100	.	+	.	ID=sorf0;Name=sORF_00000
aaa	UTR_derived	sORF	624	799	.	+	.	ID=sorf1;Name=sORF_00001;UTR_type=interCDS
aaa	UTR_derived	sORF	1101	1800	.	+	.	ID=sorf2;Name=sORF_00002;UTR_type=3utr"""
    gff_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "gene", "start": 150,
                 "end": 200, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 1230,
                 "end": 1240, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 7100,
                 "end": 9167, "phase": ".", "strand": "-", "score": "."}]
    attributes_gff = [{"ID": "gene0", "Name": "Gene_0", "locus_tag": "AAA_00001"},
                      {"ID": "cds1", "Name": "CDS_1", "locus_tag": "AAA_00002"},
                      {"ID": "cds4", "Name": "CDS_4", "locus_tag": "BBB_00003"}] 
    tran_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "transcript", "start": 10,
                  "end": 100, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "transcript", "start": 1100,
                  "end": 1240, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "transcript", "start": 8000,
                  "end": 900, "phase": ".", "strand": "-", "score": "."}]
    attributes_tran = [{"ID": "tran0", "Name": "tran_0"},
                       {"ID": "tran1", "Name": "tran_1"},
                       {"ID": "tran4", "Name": "tran_4"}]
    gffs = []
    trans = []
    for index in range(0, 3):
        gffs.append(Create_generator(gff_dict[index], attributes_gff[index], "gff"))
        trans.append(Create_generator(tran_dict[index], attributes_tran[index], "gff"))

if __name__ == "__main__":
    unittest.main()

