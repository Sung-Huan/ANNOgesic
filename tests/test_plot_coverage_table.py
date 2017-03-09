import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import gen_file
import annogesiclib.plot_coverage_table as pct


class Mock_func(object):

    def mock_fig(self, rowlabels, collabels, cells, filename,
                 max_color, min_color):
        gen_file(filename, "test")
        pass

class TestPlotCoverageTable(unittest.TestCase):

    def setUp(self):
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)


    def test_plot_table(self):
        pct.fig = Mock_func().mock_fig
        plots = [{"aaa": {"cond_1": {"track_1": 3.543, "track_2": 4.523},
                          "cond_2": {"track_1": 4.43, "track_2": 0.523}}}]
        pct.plot_table(plots, 100, 0, os.path.join(self.test_folder, "test"))
        self.assertTrue(os.path.exists(os.path.join(self.test_folder, "test")))

if __name__ == "__main__":
    unittest.main()

