import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file, import_data
from annogesiclib.seqmodifier import SeqModifier


class TestSeqModifier(unittest.TestCase):

    def setUp(self):
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
        self.seq = SeqModifier("AATTATATAGGAAGGCCC")

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_init_pos_dict(self):
        self.seq._init_pos_dict()
        self.assertDictEqual(self.seq._org_pos_to_internal_pos,
                             {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6,
                              8: 7, 9: 8, 10: 9, 11: 10, 12: 11, 13: 12,
                              14: 13, 15: 14, 16: 15, 17: 16, 18: 17})

    def test_replace(self):
        self.seq.replace(2, "G")
        self.assertEqual(self.seq._seq, "AGTTATATAGGAAGGCCC")

    def test_remove(self):
        self.seq.remove(8)
        self.assertEqual(self.seq._seq, "AATTATAAGGAAGGCCC")

    def test_insert(self):
        self.seq.insert(5, "C")
        self.assertEqual(self.seq._seq, "AATTCATATAGGAAGGCCC")

if __name__ == "__main__":
    unittest.main()
