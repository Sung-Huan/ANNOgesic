import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data
import annogesiclib.circRNA as circ

class Mock_read_file(object):

    def __init__(self):
        self.example = Example()

    def read_file(self, gff_file, input_file):
        self.circs = []
        self.gffs = []
        for index in range(0, 5):
            self.circs.append(Create_generator(self.example.circ_dict[index],
                                          self.example.attributes_circ[index], "circ"))
        for index in range(0, 3):
            self.gffs.append(Create_generator(self.example.gffs_dict[index],
                                         self.example.attributes_gffs[index], "gff"))
        return self.circs, self.gffs, 50

class TestCircRNA(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_get_feature(self):
        attributes_cds = {"ID": "cds0", "Name": "CDS_0", "locus_tag": "AAA_00001",
                           "protein_id": "YP_918384.3"}
        attributes = circ.get_feature(Create_generator(self.example.cds_dict,
                                                       attributes_cds, "gff"))
        self.assertEqual(attributes, "AAA_00001")
        attributes_cds = {"ID": "cds0", "Name": "CDS_0", "protein_id": "YP_918384.3"}
        attributes = circ.get_feature(Create_generator(self.example.cds_dict,
                                                       attributes_cds, "gff"))
        self.assertEqual(attributes, "YP_918384.3")
        attributes_cds = {"ID": "cds0", "Name": "CDS_0"}
        attributes = circ.get_feature(Create_generator(self.example.cds_dict,
                                                       attributes_cds, "gff"))
        self.assertEqual(attributes, "cds0:122-267_f")

    def test_detect_conflict(self):
        circ_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "circRNA", "start": 100,
                     "end": 467, "phase": ".", "strand": "+", "score": ".", "support": 30,
                     "start_site": 30, "end_site": 35, "situation": "P", "splice_type": "C"}
        attributes_circ = {"ID": "circrna0", "Name": "circRNA_0"}
        circrna = Create_generator(circ_dict, attributes_circ, "circ")
        gffs = [Create_generator(self.example.cds_dict, self.example.attributes_cds, "gff")]
        out = StringIO()
        circ.detect_conflict(gffs, circrna, 0, out)
        self.assertEqual(out.getvalue(), "circRNA_0	aaa	+	100	467	AAA_00001	30	1.0	0.8571428571428571\n")
        out.close()

    def test_get_circrna(self):
        circs = []
        gffs = []
        for index in range(0, 5):
            circs.append(Create_generator(self.example.circ_dict[index],
                                          self.example.attributes_circ[index], "circ"))
        for index in range(0, 3):
            gffs.append(Create_generator(self.example.gffs_dict[index],
                                         self.example.attributes_gffs[index], "gff"))
        out = StringIO()
        nums = circ.get_circrna(circs, gffs, 50, 0.5, 0.5, out)
        self.assertDictEqual(nums["support"], {'aaa': {0: 2, 20: 1, 5: 2, 25: 1, 10: 2, 30: 1, 15: 1},
                                               'all': {0: 3, 20: 1, 5: 3, 25: 1, 10: 2, 30: 1, 15: 1},
                                               'bbb': {0: 1, 5: 1}})
        self.assertDictEqual(nums["circular"], {'bbb': 1, 'aaa': 2, 'all': 3})
        self.assertDictEqual(nums["conflict"], {'bbb': {0: 1, 5: 1},
                                                'aaa': {},
                                                'all': {0: 1, 5: 1}})

    def test_detect_circrna(self):
        out_file = os.path.join(self.test_folder, "out.csv")
        stat_file = os.path.join(self.test_folder, "stat.csv")
        circ.read_file = Mock_read_file().read_file
        circ.detect_circrna("test.circ", "test.gff", out_file,
                            0.5, 0.5, stat_file)
        circs = import_data(out_file)
        stats = import_data(stat_file)
        self.assertEqual(set(circs), set(self.example.out_file.split("\n")))
        self.assertEqual(set(stats), set(self.example.stat_file.split("\n")))

class Example(object):

    cds_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 122,
                "end": 267, "phase": ".", "strand": "+", "score": "."}
    attributes_cds = {"ID": "cds0", "Name": "CDS_0", "locus_tag": "AAA_00001"}

    circ_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "circRNA", "start": 100,
                  "end": 467, "phase": ".", "strand": "+", "score": ".", "support": 30,
                  "start_site": 30, "end_site": 35, "situation": "P", "splice_type": "C"},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "circRNA", "start": 1330,
                  "end": 1564, "phase": ".", "strand": "+", "score": ".", "support": 10,
                  "start_site": 30, "end_site": 40, "situation": "P", "splice_type": "C"},
                 {"seq_id": "bbb", "source": "Refseq", "feature": "circRNA", "start": 30,
                  "end": 167, "phase": ".", "strand": "-", "score": ".", "support": 5,
                  "start_site": 5, "end_site": 6, "situation": "P", "splice_type": "C"},
                 {"seq_id": "bbb", "source": "Refseq", "feature": "circRNA", "start": 230,
                  "end": 467, "phase": ".", "strand": "-", "score": ".", "support": 50,
                  "start_site": 60, "end_site": 60, "situation": "F", "splice_type": "C"},
                 {"seq_id": "bbb", "source": "Refseq", "feature": "circRNA", "start": 1230,
                  "end": 2167, "phase": ".", "strand": "-", "score": ".", "support": 50,
                  "start_site": 60, "end_site": 60, "situation": "P", "splice_type": "L"}]
    attributes_circ = [{"ID": "circrna0", "Name": "circRNA_0"},
                       {"ID": "circrna1", "Name": "circRNA_1"},
                       {"ID": "circrna2", "Name": "circRNA_2"},
                       {"ID": "circrna2", "Name": "circRNA_3"},
                       {"ID": "circrna2", "Name": "circRNA_4"}]
    gffs_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 140,
                  "end": 367, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 30,
                  "end": 40, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "bbb", "source": "Refseq", "feature": "circRNA", "start": 430,
                  "end": 5167, "phase": ".", "strand": "-", "score": "."}]
    attributes_gffs = [{"ID": "CDS0", "Name": "CDS_0", "locus_tag": "AAA_00001"},
                       {"ID": "CDS1", "Name": "CDS_1", "locus_tag": "AAA_00002"},
                       {"ID": "CDS2", "Name": "CDS_2", "locus_tag": "BBB_00001"}]
    out_file = """ID	strain	strand	start	end	annotation_overlap	supported_reads	supported_reads/reads_at_start	supported_reads/reads_at_end
circRNA_0	aaa	+	100	467	AAA_00001	30	1.0	0.8571428571428571
circRNA_1	aaa	+	1330	1564	NA	10	0.3333333333333333	0.25
circRNA_2	bbb	-	30	167	NA	5	1.0	0.8333333333333334"""

    stat_file = """All strains:
	the number of all circular RNAs = 3
	the number of potential circular RNAs, more than 0 supported it = 3
	the number of potential circular RNAs, more than 5 supported it = 3
	the number of potential circular RNAs, more than 10 supported it = 2
	the number of potential circular RNAs, more than 15 supported it = 1
	the number of potential circular RNAs, more than 20 supported it = 1
	the number of potential circular RNAs, more than 25 supported it = 1
	the number of potential circular RNAs, more than 30 supported it = 1

	the circular RNAs:
		without conflict with annotation
		support reat ratio of starting point is larger than 0.5
		support reat ratio of end point is larger than 0.5
	the number of potential circular RNAs, more than 0 supported it = 1
	the number of potential circular RNAs, more than 5 supported it = 1

bbb:
	the number of all circular RNAs = 1
	the number of potential circular RNAs, more than 0 supported it = 1
	the number of potential circular RNAs, more than 5 supported it = 1

	the circular RNAs:
		without conflict with annotation
		support reat ratio of starting point is larger than 0.5
		support reat ratio of end point is larger than 0.5
	the number of potential circular RNAs, more than 0 supported it = 1
	the number of potential circular RNAs, more than 5 supported it = 1

aaa:
	the number of all circular RNAs = 2
	the number of potential circular RNAs, more than 0 supported it = 2
	the number of potential circular RNAs, more than 5 supported it = 2
	the number of potential circular RNAs, more than 10 supported it = 2
	the number of potential circular RNAs, more than 15 supported it = 1
	the number of potential circular RNAs, more than 20 supported it = 1
	the number of potential circular RNAs, more than 25 supported it = 1
	the number of potential circular RNAs, more than 30 supported it = 1

	the circular RNAs:
		without conflict with annotation
		support reat ratio of starting point is larger than 0.5
		support reat ratio of end point is larger than 0.5"""

if __name__ == "__main__":
    unittest.main()

