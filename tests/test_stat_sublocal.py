import sys
import os
import unittest
import shutil
import copy
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data, gen_file, extract_info
import annogesiclib.stat_sublocal as ss


class TestStatSubLocal(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_read_table(self):
        psortb_file = os.path.join(self.test_folder, "test.csv")
        gen_file(psortb_file, self.example.table)
        subs, total_nums, unknown_nums = ss.read_table(psortb_file)
        self.assertDictEqual(subs, {'Staphylococcus_aureus_HG002': {'Unknown': 1},
                                    'Staphylococcus_aureus_HG003': {'CellWall': 1, 'Cytoplasmic': 2},
                                    'all_strain': {'Unknown': 1, 'CellWall': 1, 'Cytoplasmic': 2}})
        self.assertDictEqual(total_nums, {'Staphylococcus_aureus_HG002': 1,
                                          'Staphylococcus_aureus_HG003': 3, 'all_strain': 4})
        self.assertDictEqual(unknown_nums, {'Staphylococcus_aureus_HG002': 1,
                                            'Staphylococcus_aureus_HG003': 0, 'all_strain': 1})

    def test_print_file_and_plot(self):
        out_stat = StringIO()
        sub = {'Unknown': 1, 'CellWall': 1, 'Cytoplasmic': 2}
        total_nums = {'Staphylococcus_aureus_HG002': 1,
                      'Staphylococcus_aureus_HG003': 3, 'all_strain': 4}
        unknown_nums = {'Staphylococcus_aureus_HG002': 1,
                        'Staphylococcus_aureus_HG003': 0, 'all_strain': 1}
        ss.print_file_and_plot(sub, total_nums, unknown_nums,
                               "all_strain", out_stat, self.test_folder + "/")
        datas = out_stat.getvalue().split("\n")
        for data in datas:
            if "Total with Unknown" in data:
                self.assertEqual(data, "Total including Unknown is 4; Total excluding Unknown is 3")
            elif "CellWall" in data:
                self.assertEqual(data, "\tCellWall\t1(including Unknown 0.25; excluding Unknonwn 0.3333333333333333)")
            elif "Cytoplasmic" in data:
                self.assertEqual(data, "\tCytoplasmic\t2(including Unknown 0.5; excluding Unknonwn 0.6666666666666666)")
            else:
                if "include Unknown" in data:
                    self.assertEqual(data, "\tUnknown\t1(including Unknown 0.25)")

    def test_plot(self):
        subs = {'Unknown': 1, 'CellWall': 1, 'Cytoplasmic': 2}
        ss.plot(subs, 4, 1, "test", self.test_folder + "/")
        self.assertTrue(os.path.exists(os.path.join(self.test_folder, "_test_sublocal.png")))


class Example(object):

    table = """Staphylococcus_aureus_HG003	YP_498609.1	+	517	1878	Cytoplasmic	9.97
Staphylococcus_aureus_HG003	YP_498610.1	+	2156	3289	Cytoplasmic	9.97
Staphylococcus_aureus_HG003	YP_498611.1	+	3670	3915	CellWall	7.50
Staphylococcus_aureus_HG002	YP_498612.1	+	4676	5015	Unknown	7.50"""

if __name__ == "__main__":
    unittest.main()

