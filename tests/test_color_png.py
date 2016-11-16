import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file, import_data
from annogesiclib.color_png import ColorPNG


class Mock_func(object):
    def __init__(self):
        self.color = ColorPNG()

    def mock_convert_svg(self, imagemagick_path, out_path, screenshot, svg_file):
        gen_file(os.path.join(out_path, svg_file), "<svg width=1111 hight=2222")

    def mock_convert_png(self, imagemagick_path, out_path, screenshot, png_file):
        gen_file(os.path.join(out_path, png_file), "test")
        pass

    def mock_gen_svg(out_path, track_num, height, width):
        pass

class TestColorPng(unittest.TestCase):

    def setUp(self):
        self.mock = Mock_func()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
            os.mkdir(os.path.join(self.test_folder, "screenshots"))
            os.mkdir(os.path.join(self.test_folder, "screenshots", "aaa"))
            os.mkdir(os.path.join(self.test_folder, "screenshots", "aaa", "forward"))
            os.mkdir(os.path.join(self.test_folder, "screenshots", "aaa", "reverse"))
        gen_file(os.path.join(self.test_folder, "screenshots", "aaa", "forward",
                              "test_f.png"), "None")
        gen_file(os.path.join(self.test_folder, "screenshots", "aaa", "reverse",
                              "test_r.png"), "None")
        self.color = ColorPNG()

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_generate_color_png(self):
        self.color._convert_svg = self.mock.mock_convert_svg
        self.color._convert_png = self.mock.mock_convert_png
        self.color.gen_svg = self.mock.mock_gen_svg
        self.color.generate_color_png(4, self.test_folder, "test")
        data = import_data(os.path.join(self.test_folder, "screenshots", "aaa",
                                        "forward", "test_f.png"))
        self.assertListEqual(data, ["test"])

if __name__ == "__main__":
    unittest.main()

