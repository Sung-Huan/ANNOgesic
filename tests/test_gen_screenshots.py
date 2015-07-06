import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
import annogesiclib.gen_screenshots as gs

class Mock_func(object):

    def __init__(self):
        self.example = Example()

    def mock_import_wig(self, lib, wigs, strand):
        pass

class TestGenScreenshots(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_set_data_range(self):
        gff_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 3,
                    "end": 6, "phase": ".", "strand": "+", "score": "."}
        attributes_gff = {"ID": "CDS0", "Name": "CDS_0", "locus_tag": "AAA_00001"}
        gff = Create_generator(gff_dict, attributes_gff, "gff")
        out = StringIO()
        gs.set_data_range(out, gff, self.example.wigs_low, "+")
        self.assertEqual(out.getvalue(), "setDataRange 0,20\n")
        out.close()
        out = StringIO()
        gs.set_data_range(out, gff, self.example.wigs_high, "+")
        self.assertEqual(out.getvalue(), "setDataRange 0,510\n")

    def test_print_batch(self):
        out = StringIO()
        lib_t = "wig1 wig2"
        lib_n = "wig3 wig4"
        lib_f = "wig5"
        gs.print_batch(out, "+", lib_t, lib_n, lib_f, "fasta", "main_gff", "expend",
                       "side1 side2", 1000, self.test_folder, "test")
        self.assertEqual(out.getvalue(), self.example.out)

    def test_gen_batch(self):
        gs.import_wig = Mock_func().mock_import_wig
        out = StringIO()
        lib_t = "wig1 wig2"
        lib_n = "wig3 wig4"
        lib_f = "wig5"
        gff_dict = {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 3,
                    "end": 6, "phase": ".", "strand": "+", "score": "."}
        attributes_gff = {"ID": "CDS0", "Name": "CDS_0", "locus_tag": "AAA_00001"}
        gff = Create_generator(gff_dict, attributes_gff, "gff")
        gs.gen_batch(lib_t, lib_n, lib_f, "+", [gff], out)
        self.assertEqual(out.getvalue(), self.example.out_print_wig)

class Example(object):
    cover_low = [{"coverage": 1.342}, {"coverage": 2.341}, {"coverage": 2.3544},
                 {"coverage": 5.342}, {"coverage": 10.2341}, {"coverage": 6.231},
                 {"coverage": 1.432}, {"coverage": 1.342}, {"coverage": 1.342}]
    covers_low = []
    for cover in cover_low:
        covers_low.append(Create_generator(cover, {}, "wig"))
    wigs_low = {"test.wig": {"aaa": covers_low, "bbb": covers_low}}
    cover_high = [{"coverage": 100.342}, {"coverage": 20.341}, {"coverage": 20.3544},
                 {"coverage": 500.342}, {"coverage": 10.2341}, {"coverage": 60.231},
                 {"coverage": 100.432}, {"coverage": 100.342}, {"coverage": 10.342}]
    covers_high = []
    for cover in cover_high:
        covers_high.append(Create_generator(cover, {}, "wig"))
    wigs_high = {"test.wig": {"aaa": covers_high, "bbb": covers_high}}
    out = """new
genome /home/silas/ANNOgesic/fasta
load /home/silas/ANNOgesic/main_gff
expend main_gff
load /home/silas/ANNOgesic/s
expend s
load /home/silas/ANNOgesic/i
expend i
load /home/silas/ANNOgesic/d
expend d
load /home/silas/ANNOgesic/e
expend e
load /home/silas/ANNOgesic/1
expend 1
load /home/silas/ANNOgesic/ 
expend  
load /home/silas/ANNOgesic/s
expend s
load /home/silas/ANNOgesic/i
expend i
load /home/silas/ANNOgesic/d
expend d
load /home/silas/ANNOgesic/e
expend e
load /home/silas/ANNOgesic/2
expend 2
load /home/silas/ANNOgesic/w
load /home/silas/ANNOgesic/w
load /home/silas/ANNOgesic/i
load /home/silas/ANNOgesic/i
load /home/silas/ANNOgesic/g
load /home/silas/ANNOgesic/g
load /home/silas/ANNOgesic/1
load /home/silas/ANNOgesic/3
load /home/silas/ANNOgesic/ 
load /home/silas/ANNOgesic/ 
load /home/silas/ANNOgesic/w
load /home/silas/ANNOgesic/w
load /home/silas/ANNOgesic/i
load /home/silas/ANNOgesic/i
load /home/silas/ANNOgesic/g
load /home/silas/ANNOgesic/g
load /home/silas/ANNOgesic/2
load /home/silas/ANNOgesic/4
load /home/silas/ANNOgesic/w
load /home/silas/ANNOgesic/i
load /home/silas/ANNOgesic/g
load /home/silas/ANNOgesic/5
maxPanelHeight 1000
snapshotDirectory /home/silas/ANNOgesic/test_folder/test/forward
"""
    out_print_wig = """goto aaa:-197-206
setDataRange 0,10
snapshot aaa:3-6.png
"""
if __name__ == "__main__":
    unittest.main()

