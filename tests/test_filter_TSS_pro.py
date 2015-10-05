#!/usr/bin/python

import unittest
import os
import sys
import shutil
sys.path.append(".")
from io import StringIO
from mock_gff3 import Create_generator
import annogesiclib.filter_TSS_pro as ftp


class TestFilterTSSPro(unittest.TestCase):

    def setUp(self):
        self.test_folder = "test_project"
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)
        os.mkdir(self.test_folder)
        self.example = Example()

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_compare_tss_pro(self):
        out = StringIO()
        ftp.compare_tss_pro(self.example.tars, self.example.refs, out, 3)
        self.assertEqual("\t".join(out.getvalue().split("\t")[0:-1]), 
                         "aaa\tRefseq\tTSS\t24\t24\t.\t+\t.")


class Example(object):
    tar_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 3,
                 "end": 3, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 24,
                 "end": 24, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 1243,
                 "end": 1243, "phase": ".", "strand": "+", "score": "."}]
    attributes_tar = [{"coverage": "3", "ID": "tss1", "Name": "TSS:3_+"},
                      {"coverage": "340", "ID": "tss2", "Name": "TSS:24_+"},
                      {"coverage": "4440", "ID": "tss3", "Name": "TSS:1243_+"}]
    ref_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "Pro", "start": 3,
                 "end": 3, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "Pro", "start": 333,
                 "end": 333, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "Pro", "start": 1242,
                 "end": 1242, "phase": ".", "strand": "+", "score": "."}]
    attributes_ref = [{"coverage": "3", "ID": "pro1", "Name": "Pro:3_+"},
                      {"coverage": "330", "ID": "pro2", "Name": "Pro:333_+"},
                      {"coverage": "1230", "ID": "pro3", "Name": "Pro:1242_+"}]
    tars = []
    refs = []
    for index in range(0, 3):
        tars.append(Create_generator(tar_dict[index], attributes_tar[index], "gff"))
        refs.append(Create_generator(ref_dict[index], attributes_ref[index], "gff"))

if __name__ == "__main__":
    unittest.main()
