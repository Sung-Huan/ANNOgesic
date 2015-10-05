import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import gen_file
import annogesiclib.plot_mountain as pm


class TestPlotMountain(unittest.TestCase):

    def setUp(self):
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
        self.example = Example()

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_plot_mountain_plot(self):
        gen_file(os.path.join(self.test_folder, "test"), self.example.mountain)
        pm.plot_mountain_plot(os.path.join(self.test_folder, "test"),
                              os.path.join(self.test_folder, "out"))

        self.assertTrue(os.path.exists(os.path.join(self.test_folder, "out")))

class Example(object):

    mountain = """   1        0
   2  0.001304
   3  0.0037577
   4  0.0068858
   5  0.015473
   6  0.025351
   7  0.71432
   8   1.6366
   9    2.615
  10   3.6091
&
   1     0
   2     0
   3     0
   4     0
   5     0
   6     0
   7     1
   8     2
   9     3
  10     4
&
   1  0.018708
   2  0.035075
   3  0.043831
   4  0.093979
   5  0.10259
   6  0.96026
   7   0.4509
   8  0.18699
   9  0.062985
  10  0.0055594
&
   1     0
   2     0
   3     0
   4     0
   5     0
   6     0
   7     1
   8     2
   9     3
  10     4"""

if __name__ == "__main__":
    unittest.main()

