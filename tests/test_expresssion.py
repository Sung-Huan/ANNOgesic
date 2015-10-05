import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file, import_data
import annogesiclib.expression as express_file
from annogesiclib.expression import Expression


class Mock_func(object):

    def mock_expression(self, input_libs, gffs, percent_tex, percent_frag, wig_f_file,
                        wig_r_file, features, merge_wigs, cutoff_coverage,
                        tex_notex, replicates, stat, gff_folder,
                        cover_type, max_color, min_color):
        pass

class TestExpression(unittest.TestCase):

    def setUp(self):
        self.mock = Mock_func()
        self.test_folder = "test_folder"
        self.tex_path = os.path.join(self.test_folder, "tex")
        self.frag_path = os.path.join(self.test_folder, "frag")
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
            os.mkdir(self.tex_path)
            os.mkdir(self.frag_path)
        self.express = Expression(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_get_replicates(self):
        replicates = self.express._get_replicates(2, 1)
        self.assertDictEqual({'tex': 2, 'frag': 1}, replicates)

    def test_expression(self):
        express_file.gene_expression = self.mock.mock_expression
        tex_libs="tex_-TEX_forward.wig:notex:1:a:+ \
                  tex_-TEX_reverse.wig:notex:1:a:- \
                  tex_+TEX_forward.wig:tex:1:a:+ \
                  tex_+TEX_reverse.wig:tex:1:a:-"
        frag_libs="frag_forward.wig:frag:1:a:+ \
                   frag_reverse.wig:frag:1:a:-"
        gen_file(os.path.join(self.frag_path, "tex_-TEX_forward.wig"), "tex1")
        gen_file(os.path.join(self.frag_path, "tex_-TEX_reverse.wig"), "tex2")
        gen_file(os.path.join(self.frag_path, "tex_+TEX_forward.wig"), "tex3")
        gen_file(os.path.join(self.frag_path, "tex_+TEX_reverse.wig"), "tex4")
        gen_file(os.path.join(self.frag_path, "frag_forward.wig"), "frag1")
        gen_file(os.path.join(self.frag_path, "frag_reverse.wig"), "frag2")
        self.express.expression(tex_libs, frag_libs, 2, 2,
                                1, self.tex_path, self.frag_path, "all",
                                "all", 5, self.test_folder, "CDS",
                                "high", 100, 0)
        self.assertTrue(os.path.exists(os.path.join(self.test_folder, "for_libs")))
        self.assertTrue(os.path.exists(os.path.join(self.test_folder, "for_libs", "statistics")))
        self.assertTrue(os.path.exists(os.path.join(self.test_folder, "for_libs", "gffs")))


if __name__ == "__main__":
    unittest.main()

