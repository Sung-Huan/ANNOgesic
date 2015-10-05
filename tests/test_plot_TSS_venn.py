import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import gen_file
import annogesiclib.plot_TSS_venn as ptv


class Mock_func(object):

    def mock_plot_text(plt, xy1, xy2, tss_type, size, color_text):
        pass

class TestPlotTSSVenn(unittest.TestCase):

    def setUp(self):
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
        self.mock = Mock_func()
        self.example = Example()

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_check_tss_class(self):
        strain = "test"
        tss_type = "Internal"
        total_types = {"test": {}}
        ptv.check_tss_class(total_types, strain, self.example.tsss[0], tss_type)
        self.assertDictEqual(total_types, {'test': {'Internal': 0}})
        tss_type = "Primary"
        ptv.check_tss_class(total_types, strain, self.example.tsss[0], tss_type)
        self.assertDictEqual(total_types, {'test': {'Internal': 0, 'Primary': 1}})

    def test_import_types(self):
        tsss = {"test": self.example.tsss}
        types, total_types = ptv.import_types(tsss)
        self.assertDictEqual(types, {'test': {'Orphan': 1, 'Internal': 1,
                                     'Primary': 1}, 'all': {}})
        self.assertDictEqual(total_types, {'test': {'Orphan': 1, 'Antisense': 0,
                                           'Secondary': 0, 'Internal': 1,
                                           'Primary': 1}, 'all': {}})

    def test_read_gff(self):
        tss_file = os.path.join(self.test_folder, "test.gff")
        gen_file(tss_file, self.example.tss_file)
        tsss, tss_num = ptv.read_gff(tss_file)
        self.assertEqual(tsss["all"][0].start, 140)
        self.assertEqual(tsss["aaa"][0].start, 140)
        self.assertDictEqual(tss_num, {'all': 1, 'aaa': 1})

    def test_plot(self):
        types = {'test': {'Orphan': 1, 'Internal': 1,
                          'Primary': 1}, 'all': {}}
        total_types = {'test': {'Orphan': 1, 'Antisense': 0,
                                'Secondary': 0, 'Internal': 1,
                                'Primary': 1}, 'all': {}}
        tss_num = {'all': 0, 'test': 3}
        ptv.plot(types, "TSS", "TSS", total_types, tss_num)
        self.assertTrue(os.path.exists("TSS_venn_test.png"))
        os.remove("TSS_venn_test.png")

class Example(object):

    tss_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 140,
                 "end": 140, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "TSS", "start": 230,
                 "end": 230, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "bbb", "source": "Refseq", "feature": "TSS", "start": 5166,
                 "end": 5166, "phase": ".", "strand": "-", "score": "."}]
    attributes_tss = [{"ID": "tss0", "Name": "TSS_0", "type": "Primary", "associated_gene": "AAA_00001"},
                      {"ID": "tss1", "Name": "TSS_1", "type": "Internal", "associated_gene": "AAA_00002"},
                      {"ID": "tss2", "Name": "TSS_2", "type": "Orphan", "associated_gene": "orphan"}]
    tsss = []
    for index in range(0, 3):
        tsss.append(Create_generator(tss_dict[index], attributes_tss[index], "gff"))
    tss_file = """aaa\tRefseq\tTSS\t140\t140\t.\t+\t.\tID=TSS_0;Name=TSS_00000;associated_gene=AAA_00001;type=Primary"""

if __name__ == "__main__":
    unittest.main()

