import hashlib
import unittest
import os

class MyFirstTest(unittest.TestCase):
 
   def setUp(self):
        self.default_files = []
        out = open("default_refernce.txt", "w")
        for file_ in os.listdir("test_folder"):
            md5_ref = hashlib.md5(open("test_folder/" + file_).read()).hexdigest()
            self.default_files.append(md5_ref)
            out.write("\t".join([file_, md5_ref]) + "\n") 
   def test_compare_string(self):
        compare_files = []
        for file_ in os.listdir("test_test"):
            compare_files.append(hashlib.md5(open("test_test/" + file_).read()).hexdigest())       
        self.assertEqual(self.default_files, compare_files)

if __name__ == '__main__':
    unittest.main()
