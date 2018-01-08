import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import gen_file, import_data
import annogesiclib.ribo_gff as rg


class TestRiboGff(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_read_file(self):
        ribo_table = os.path.join(self.test_folder, "table")
        rfam_table = os.path.join(self.test_folder, "rfam")
        gen_file(ribo_table, self.example.table)
        gen_file(rfam_table, self.example.rfam)
        ribos, rfams = rg.read_file(ribo_table, rfam_table)
        self.assertListEqual(ribos, self.example.ribos)
        self.assertListEqual(rfams, self.example.rfams)

    def test_get_overlap(self):
        pre_ribo = {"strain": "aaa", "strand": "+",
                    "ID": "id_1", "info": "test_1"}
        ribo = {"strain": "aaa", "strand": "+",
                "ID": "id_1", "info": "test_2"}
        overlaps = {}
        overlaps["aaa"] = []
        rg.get_overlap(pre_ribo, ribo, False, overlaps)
        self.assertDictEqual(overlaps, {'aaa': ['test_1;test_2']})

    def test_print_gff(self):
        out = StringIO()
        stats = {}
        stats["test_1"] = {"total": 0}
        stats["test_2"] = {"total": 0}
        stats["total"] = {"total": 0}
        ribo = {'associate': 'SAOUHSC_00013', 'start_seq': 15948,
                'start_align': 1, 'strain': 'test_1', 'end_align': 99,
                'ID': 'riboswitch_5',
                'info': 'riboswitch_5|Staphylococcus_aureus_HG003|+|SAOUHSC_00013|15948|16046',
                'strand': '+', 'e': '1.6e-18', 'score': '74.3', 'rfam': 'RF00162',
                'end_seq': 16046, "rfam_name": "SAM"}
        rg.print_gff(1, ribo, out, stats, "test_1", "riboswitch")
        self.assertDictEqual(stats, {
            'total': {'total': 1}, 'test_1': {'total': 1},
            'test_2': {'total': 0}})
        self.assertEqual(
            out.getvalue(),
            "test_1\tANNOgesic\triboswitch\t15948\t16046\t.\t+\t.\tID=test_1_riboswitch1;Name=SAM;rfam_id=RF00162;e_value=1.6e-18;score=74.3;method=infernal_to_Rfam\n")

    def test_import_stat(self):
        ribo = {'associate': 'SAOUHSC_00013', 'start_seq': 15948,
                'start_align': 1, 'strain': 'Staphylococcus_aureus_HG003',
                'end_align': 99, 'ID': 'riboswitch_5',
                'info': 'riboswitch_5|Staphylococcus_aureus_HG003|+|SAOUHSC_00013|15948|16046',
                'strand': '+', 'e': '1.6e-18', 'rfam': 'RF00162',
                'end_seq': 16046, "rfam_name": "RF00162"}
        stats = {}
        stats["Staphylococcus_aureus_HG003"] = {"total": 0}
        stats["test"] = {"total": 0}
        stats["total"] = {"total": 0}
        rg.import_stat(self.example.rfams, ribo, stats,
                       "Staphylococcus_aureus_HG003")
        self.assertDictEqual(stats, {
            'total': {'total': 0, 'SAM': 1},
            'Staphylococcus_aureus_HG003': {'total': 0, 'SAM': 1},
            'test': {'total': 0}})

    def test_print_number(self):
        stats = {}
        stats["Staphylococcus_aureus_HG003"] = {"total": 100}
        stats["test"] = {"total": 20}
        stats["total"] = {"total": 120}
        out = StringIO()
        rg.print_number(stats, 2, out, "Staphylococcus_aureus_HG003",
                        "riboswitch")
        ref = """Total number of potential riboswitch are 100
The number of potential riboswitch which have overlap region with others are 2
riboswitch_name\tnumbers
"""
        self.assertEqual(out.getvalue(), ref)

    def test_print_stat(self):
        out_stat = os.path.join(self.test_folder, "test")
        stats = {}
        stats["Staphylococcus_aureus_HG003"] = {"total": 100}
        stats["total"] = {"total": 120}
        overlaps = {'Staphylococcus_aureus_HG003': ['test_1;test_2']}
        rg.print_stat(stats, out_stat, overlaps, "riboswitch")
        data = import_data(out_stat)
        ref = """Staphylococcus_aureus_HG003:
Total number of potential riboswitch are 100
The number of potential riboswitch which have overlap region with others are 2
riboswitch_name	numbers

overlap candidates set 1:
	test_1
	test_2"""
        self.assertEqual("\n".join(data), ref)

    def test_stat_and_covert2gff(self):
        ribo_table = os.path.join(self.test_folder, "ribo")
        rfam_table = os.path.join(self.test_folder, "rfam")
        gen_file(ribo_table, self.example.table)
        gen_file(rfam_table, self.example.rfam)
        gff_file = os.path.join(self.test_folder, "gff")
        out_stat = os.path.join(self.test_folder, "stat")
        rg.stat_and_covert2gff(ribo_table, rfam_table, gff_file,
                               3, out_stat, "riboswitch")
        data = import_data(gff_file)
        self.assertEqual("\n".join(data), self.example.out_gff)


class Example(object):

    table = """#ID\tstrain\tstrand\tassociated_CDS\tstart\tend	Rfam	e_value	score	start_align	end_align
riboswitch_5\tStaphylococcus_aureus_HG003\t+\tSAOUHSC_00013\t15948\t16046	RF00162	1.6e-18	74.3	1	99
riboswitch_11\tStaphylococcus_aureus_HG003\t-\tSAOUHSC_00007\t27955\t28053	RF00162	1.6e-18	74.3	1	99
riboswitch_183\tStaphylococcus_aureus_HG003\t+\tSAOUHSC_00372\t377996\t378098	RF00167	2.2e-18	45.1	1	103"""

    rfam = """RF00162	SAM	SAM riboswitch box leader
RF00059	TPP	TPP riboswitch THI element
RF00167	Purine	Purine riboswitch
RF00168	Lysine	Lysine riboswitch
RF00520	ybhL	ybhL leader"""

    ribos = [{'associate': 'SAOUHSC_00013', 'start_seq': 15948,
              'start_align': 1, 'strain': 'Staphylococcus_aureus_HG003',
              'end_align': 99, 'ID': 'riboswitch_5',
              'info': 'riboswitch_5|Staphylococcus_aureus_HG003|+|SAOUHSC_00013|15948|16046',
              'strand': '+', 'e': '1.6e-18', 'score': '74.3', 'rfam': 'RF00162',
              'end_seq': 16046},
             {'associate': 'SAOUHSC_00007', 'start_seq': 27955,
              'start_align': 1, 'strain': 'Staphylococcus_aureus_HG003',
              'end_align': 99, 'ID': 'riboswitch_11',
              'info': 'riboswitch_11|Staphylococcus_aureus_HG003|-|SAOUHSC_00007|27955|28053',
              'strand': '-', 'e': '1.6e-18', 'score': '74.3', 'rfam': 'RF00162',
              'end_seq': 28053},
             {'associate': 'SAOUHSC_00372', 'start_seq': 377996,
              'start_align': 1, 'strain': 'Staphylococcus_aureus_HG003',
              'end_align': 103, 'ID': 'riboswitch_183',
              'info': 'riboswitch_183|Staphylococcus_aureus_HG003|+|SAOUHSC_00372|377996|378098',
              'strand': '+', 'e': '2.2e-18', 'score': '45.1', 'rfam': 'RF00167',
              'end_seq': 378098}]
    rfams = [{'class': 'SAM', 'ID': 'RF00162'},
             {'class': 'TPP', 'ID': 'RF00059'},
             {'class': 'Purine', 'ID': 'RF00167'},
             {'class': 'Lysine', 'ID': 'RF00168'},
             {'class': 'ybhL', 'ID': 'RF00520'}]

    out_gff = """##gff-version 3
Staphylococcus_aureus_HG003	ANNOgesic	riboswitch	15948	16046	.	+	.	ID=Staphylococcus_aureus_HG003_riboswitch0;Name=SAM;rfam_id=RF00162;e_value=1.6e-18;score=74.3;method=infernal_to_Rfam
Staphylococcus_aureus_HG003	ANNOgesic	riboswitch	27955	28053	.	-	.	ID=Staphylococcus_aureus_HG003_riboswitch1;Name=SAM;rfam_id=RF00162;e_value=1.6e-18;score=74.3;method=infernal_to_Rfam
Staphylococcus_aureus_HG003	ANNOgesic	riboswitch	377996	378098	.	+	.	ID=Staphylococcus_aureus_HG003_riboswitch2;Name=Purine;rfam_id=RF00167;e_value=2.2e-18;score=45.1;method=infernal_to_Rfam"""
    out_stat = """Staphylococcus_aureus_HG003:
Total number of potential riboswitch are 3
The number of potential riboswitch which have overlap region with others are 0
riboswitch_type	numbers
Purine	1
SAM	2"""
if __name__ == "__main__":
    unittest.main()

