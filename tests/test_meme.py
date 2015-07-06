import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
import annogesiclib.meme as me
from mock_helper import gen_file
from annogesiclib.meme import MEME 


class Mock_func(object):

    def mock_del_repeat_fasta(self, tmp_fasta, all_no_orph):
        with open("tmp/all_type.fa", "w") as fh:
            fh.write("all")
        with open("tmp/without_orphan.fa", "w") as fh:
            fh.write("without_orphan")

class TestMEME(unittest.TestCase):

    def setUp(self):
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
        self.tss_folder = os.path.join(self.test_folder, "tss_folder")
        if (not os.path.exists(self.tss_folder)):
            os.mkdir(self.tss_folder)
        self.gff_folder = os.path.join(self.test_folder, "gff_folder")
        if (not os.path.exists(self.gff_folder)):
            os.mkdir(self.gff_folder)
        self.fa_folder = os.path.join(self.test_folder, "fa_folder")
        if (not os.path.exists(self.fa_folder)):
            os.mkdir(self.fa_folder)
        self.meme = MEME(self.tss_folder, self.gff_folder, self.fa_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_move_and_merge_fasta(self):
        me.del_repeat_fasta = Mock_func().mock_del_repeat_fasta
        if (not os.path.exists("tmp")):
            os.mkdir("tmp")
        gen_file("tmp/primary.fa", "primary")
        gen_file("tmp/secondary.fa", "secondary")
        gen_file("tmp/internal.fa", "internal")
        gen_file("tmp/antisense.fa", "antisense")
        gen_file("tmp/orphan.fa", "orphan")
        self.meme._move_and_merge_fasta(self.test_folder, "test")
        self.assertTrue(os.path.exists(os.path.join(self.test_folder, "test_allstrain_all_types.fa")))
        self.assertTrue(os.path.exists(os.path.join(self.test_folder, "test_allstrain_primary.fa")))
        self.assertTrue(os.path.exists(os.path.join(self.test_folder, "test_allstrain_secondary.fa")))
        self.assertTrue(os.path.exists(os.path.join(self.test_folder, "test_allstrain_internal.fa")))
        self.assertTrue(os.path.exists(os.path.join(self.test_folder, "test_allstrain_antisense.fa")))
        self.assertTrue(os.path.exists(os.path.join(self.test_folder, "test_allstrain_orphan.fa")))
        self.assertTrue(os.path.exists(os.path.join(self.test_folder, "test_allstrain_without_orphan.fa")))

    def test_split_fasta_by_strain(self):
        with open(os.path.join(self.fa_folder, "allstrain.fa"), "w") as fh:
            fh.write(""">aaa_aaa_aaa
ATTATATATA
>bbb_bbb_bbb
AATTAATTAA""")
        self.meme._split_fasta_by_strain(self.fa_folder)
        self.assertTrue(os.path.join(self.fa_folder, "aaa.fa"))
        self.assertTrue(os.path.join(self.fa_folder, "bbb.fa"))

if __name__ == "__main__":
    unittest.main()
