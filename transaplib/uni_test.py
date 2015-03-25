import hashlib
import unittest
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-r","--ref_folder",help="ref folder")
parser.add_argument("-t","--tar_folder",help="target folder")
args = parser.parse_args()

class MyFirstTest(unittest.TestCase):
 
   def setUp(self):
        self.default_files = []
        out = open("default_refernce.txt", "w")
        for file_ in os.listdir(args.input_folder):
            md5_ref = hashlib.md5(open(os.path.join(args.input_folder, file_)).read()).hexdigest()
            self.default_files.append(md5_ref)
            out.write("\t".join([file_, md5_ref]) + "\n") 
   def test_compare_string(self):
        compare_files = []
        for file_ in os.listdir(args.tar_folder):
            compare_files.append(hashlib.md5(open(os.path.join(tar_folder, file_)).read()).hexdigest())       
        self.assertEqual(self.default_files, compare_files)

if __name__ == '__main__':
    unittest.main()
