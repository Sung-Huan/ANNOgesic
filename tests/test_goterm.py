import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from annogesiclib.goterm import GoTermFinding


class TestGetPolyT(unittest.TestCase):

    def setUp(self):
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
        self.gffs = os.path.join(self.test_folder, "gff_folder")
        if (not os.path.exists(self.gffs)):
            os.mkdir(self.gffs)
        self.go_folder = os.path.join(self.test_folder, "go_folder")
        if (not os.path.exists(self.go_folder)):
            os.mkdir(self.go_folder)
        self.all_strain = "all_strains_uniprot.csv"
        self.go = GoTermFinding(self.test_folder, self.gffs)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_merge_files(self):
        gff_folder = os.path.join(self.gffs, "test.gff_folder")
        if (not os.path.exists(gff_folder)):
            os.mkdir(gff_folder)
        test1_folder = os.path.join(self.go_folder, "test1")
        if (not os.path.exists(test1_folder)):
            os.mkdir(test1_folder)
        test2_folder = os.path.join(self.go_folder, "test2")
        if (not os.path.exists(test2_folder)):
            os.mkdir(test2_folder)
        with open(os.path.join(gff_folder, "test1.gff"), "w") as fh:
            fh.write("test1")
        with open(os.path.join(gff_folder, "test2.gff"), "w") as fh:
            fh.write("test2")
        with open(os.path.join(test1_folder, "test1_uniprot.csv"), "w") as fh:
            fh.write("test1")
        with open(os.path.join(test2_folder, "test2_uniprot.csv"), "w") as fh:
            fh.write("test2")
        self.go._merge_files(self.gffs, self.go_folder, self.test_folder)
        out_file = os.path.join(self.go_folder, "test", self.all_strain)
        self.assertTrue(os.path.exists(out_file))
        with open(out_file) as fh:
            for line in fh:
                self.assertEqual(line, "test1test2")

if __name__ == "__main__":
    unittest.main()
