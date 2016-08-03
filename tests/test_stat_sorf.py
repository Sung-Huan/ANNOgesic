import sys
import os
import unittest
import shutil
import copy
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data, gen_file, extract_info
import annogesiclib.stat_sorf as ss


class TestStatsORF(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_create_dict(self):
        nums = {}
        datas = ss.create_dict(nums, "aaa", True)
        self.assertDictEqual(datas, {'aaa': {'interCDS': {'TSS_sRNA_RBS': 0, 'sRNA': 0, 'TSS': 0,
                                     'TSS_RBS': 0, 'all': 0, 'RBS_sRNA': 0, 'RBS': 0, 'TSS_sRNA': 0},
                                     'intergenic': {'TSS_sRNA_RBS': 0, 'sRNA': 0, 'TSS': 0, 'TSS_RBS': 0,
                                     'all': 0, 'RBS_sRNA': 0, 'RBS': 0, 'TSS_sRNA': 0},
                                     'antisense': {'TSS_sRNA_RBS': 0, 'sRNA': 0, 'TSS': 0, 'TSS_RBS': 0,
                                     'all': 0, 'RBS_sRNA': 0, 'RBS': 0, 'TSS_sRNA': 0},
                                     "5'UTR_derived": {'TSS_sRNA_RBS': 0, 'sRNA': 0, 'TSS': 0,
                                     'TSS_RBS': 0, 'all': 0, 'RBS_sRNA': 0, 'RBS': 0, 'TSS_sRNA': 0},
                                     'all': {'TSS_sRNA_RBS': 0, 'sRNA': 0, 'TSS': 0, 'TSS_RBS': 0,
                                     'all': 0, 'RBS_sRNA': 0, 'RBS': 0, 'TSS_sRNA': 0},
                                     "3'UTR_derived": {'TSS_sRNA_RBS': 0, 'sRNA': 0, 'TSS': 0,
                                     'TSS_RBS': 0, 'all': 0, 'RBS_sRNA': 0, 'RBS': 0, 'TSS_sRNA': 0}}})

    def test_plus_data(self):
        nums = {"aaa": {"3utr": {"CDS": 0}, "5utr": {"CDS": 2, "tRNA": 1}}}
        ss.plus_data(nums, "aaa", ["3utr", "5utr"], ["CDS"], True)
        self.assertDictEqual(nums, {'aaa': {'5utr': {'CDS': 3, 'tRNA': 1}, '3utr': {'CDS': 1}}})

    def test_print_num(self):
        out = StringIO()
        nums = {"aaa": {"3utr": {"CDS": 0, "all": 10}, "5utr": {"CDS": 2, "tRNA": 1, "all": 10},
                        "all": {"all": 20}}}
        ss.print_num(out, 2, nums, "aaa", "5utr")
        self.assertEqual(out.getvalue(), "(for strain 0.1; for 5utr - 0.2)\n")

    def test_print_stat(self):
        out = StringIO()
        nums = {"aaa": {"test1": {"sRNA": 10, "all": 30}, "test2": {"sRNA": 2, "all": 10},
                        "all": {"all": 40, "sRNA": 12}}}
        nums_best = {"aaa": {"test1": {"sRNA": 2, "all": 20}, "test2": {"sRNA": 2, "all": 10},
                             "all": {"all": 30, "sRNA": 4}}}
        ss.print_stat(nums, nums_best, "aaa", out, True)
        datas = out.getvalue().split("\n")
        for data in datas:
            if "total sORF in this strain are " in data:
                self.assertEqual(data, "\ttotal sORF in this strain are 40")
            if "total sORF of test2 sORF candidates are " in data:
                self.assertEqual(data, "\ttotal sORF of test2 sORF candidates are 10(for this strain - 0.25)")
            if "total sORF of all sORF candidates are " in data:
                self.assertEqual(data, "\ttotal sORF of all sORF candidates are 40(for this strain - 1.0)")
            if "total sORF of test1 sORF candidates are " in data:
                self.assertEqual(data, "\ttotal sORF of test1 sORF candidates are 30(for this strain - 0.75)")

    def test_get_stat_num(self):
        nums = ss.get_stat_num(self.example.sorfs, True)
        self.assertDictEqual(nums, self.example.nums)

    def test_check_class(self):
        nums = copy.deepcopy(self.example.nums)
        ss.check_class(self.example.sorfs[0], nums, "5'UTR_derived", True, "aaa")
        self.assertDictEqual(nums, {'total': {"3'UTR_derived": {'TSS_RBS': 0,
                                    'sRNA': 0, 'TSS_sRNA_RBS': 0, 'TSS': 0,
                                    'all': 0, 'TSS_sRNA': 0, 'RBS': 0, 'RBS_sRNA': 0},
                                    'all': {'TSS_RBS': 0, 'sRNA': 2, 'TSS_sRNA_RBS': 0,
                                    'TSS': 0, 'all': 2, 'TSS_sRNA': 0, 'RBS': 2, 'RBS_sRNA': 2},
                                    'interCDS': {'TSS_RBS': 0, 'sRNA': 0, 'TSS_sRNA_RBS': 0,
                                    'TSS': 0, 'all': 0, 'TSS_sRNA': 0, 'RBS': 0, 'RBS_sRNA': 0},
                                    'antisense': {'TSS_RBS': 0, 'sRNA': 0, 'TSS_sRNA_RBS': 0,
                                    'TSS': 0, 'all': 0, 'TSS_sRNA': 0, 'RBS': 0, 'RBS_sRNA': 0},
                                    'intergenic': {'TSS_RBS': 0, 'sRNA': 0, 'TSS_sRNA_RBS': 0,
                                    'TSS': 0, 'all': 0, 'TSS_sRNA': 0, 'RBS': 0, 'RBS_sRNA': 0},
                                    "5'UTR_derived": {'TSS_RBS': 0, 'sRNA': 2, 'TSS_sRNA_RBS': 0,
                                    'TSS': 0, 'all': 2, 'TSS_sRNA': 0, 'RBS': 2, 'RBS_sRNA': 2}},
                                    'aaa': {"3'UTR_derived": {'TSS_RBS': 0, 'sRNA': 0, 'TSS_sRNA_RBS': 0,
                                    'TSS': 0, 'all': 0, 'TSS_sRNA': 0, 'RBS': 0, 'RBS_sRNA': 0},
                                    'all': {'TSS_RBS': 0, 'sRNA': 2, 'TSS_sRNA_RBS': 0, 'TSS': 0,
                                    'all': 2, 'TSS_sRNA': 0, 'RBS': 2, 'RBS_sRNA': 2},
                                    'interCDS': {'TSS_RBS': 0, 'sRNA': 0, 'TSS_sRNA_RBS': 0,
                                    'TSS': 0, 'all': 0, 'TSS_sRNA': 0, 'RBS': 0, 'RBS_sRNA': 0},
                                    'antisense': {'TSS_RBS': 0, 'sRNA': 0, 'TSS_sRNA_RBS': 0, 'TSS': 0,
                                    'all': 0, 'TSS_sRNA': 0, 'RBS': 0, 'RBS_sRNA': 0},
                                    'intergenic': {'TSS_RBS': 0, 'sRNA': 0, 'TSS_sRNA_RBS': 0,
                                    'TSS': 0, 'all': 0, 'TSS_sRNA': 0, 'RBS': 0, 'RBS_sRNA': 0},
                                    "5'UTR_derived": {'TSS_RBS': 0, 'sRNA': 2, 'TSS_sRNA_RBS': 0,
                                    'TSS': 0, 'all': 2, 'TSS_sRNA': 0, 'RBS': 2, 'RBS_sRNA': 2}}})

    def test_stat(self):
        sorf_all = os.path.join(self.test_folder, "aaa_all.gff")
        sorf_best = os.path.join(self.test_folder, "aaa_best.gff")
        gen_file(sorf_all, self.example.sorf_all)
        gen_file(sorf_best, self.example.sorf_best)
        stat_file = os.path.join(self.test_folder, "stat")
        ss.stat(sorf_all, sorf_best, stat_file, True)
        datas = import_data(stat_file)
        for data in datas:
            if "total sORF in this strain are " in data:
                self.assertEqual(data, "\ttotal sORF in this strain are 3")
            if "total sORF of intergenic sORF candidates are " in data:
                self.assertEqual(data, "\ttotal sORF of intergenic sORF candidates are 1(for this strain - 0.3333333333333333)")


class Example(object):

    sorf_best = """Staphylococcus_aureus_HG003	UTR_derived	sORF	399365	399439	.	+	.	ID=sorf0;Name=sORF_00000;start_TSS=TSS_TSS_399320+;with_TSS=TSS_399320+;sRNA=NA;rbs=RBS_399354;sORF_type=5utr"""
    sorf_all = """Staphylococcus_aureus_HG003	UTR_derived	sORF	399365	399439	.	+	.	ID=sorf0;Name=sORF_00000;start_TSS=TSS_TSS_399320+;with_TSS=TSS_399320+;sRNA=NA;rbs=RBS_399354;sORF_type=5utr
Staphylococcus_aureus_HG003	UTR_derived	sORF	330	497	.	+	.	ID=sorf0;Name=sORF_00000;start_TSS=TSS_313+;with_TSS=TSS_313+;sRNA=NA;rbs=NA;sORF_type=5utr
Staphylococcus_aureus_HG003	intergenic	sORF	4090	4158	.	-	.	ID=sorf1;Name=sORF_00001;start_TSS=TSS_4368-;with_TSS=TSS_4159-,TSS_4368-;sRNA=NA;rbs=NA;sORF_type=intergenic"""
    nums = {'aaa': {"5'UTR_derived": {'TSS_RBS': 0, 'sRNA': 1, 'TSS_sRNA': 0, 'TSS_sRNA_RBS': 0, 'all': 1, 'RBS_sRNA': 1, 'TSS': 0, 'RBS': 1}, "3'UTR_derived": {'TSS_RBS': 0, 'sRNA': 0, 'TSS_sRNA': 0, 'TSS_sRNA_RBS': 0, 'all': 0, 'RBS_sRNA': 0, 'TSS': 0, 'RBS': 0}, 'antisense': {'TSS_RBS': 0, 'sRNA': 0, 'TSS_sRNA': 0, 'TSS_sRNA_RBS': 0, 'all': 0, 'RBS_sRNA': 0, 'TSS': 0, 'RBS': 0}, 'interCDS': {'TSS_RBS': 0, 'sRNA': 0, 'TSS_sRNA': 0, 'TSS_sRNA_RBS': 0, 'all': 0, 'RBS_sRNA': 0, 'TSS': 0, 'RBS': 0}, 'all': {'TSS_RBS': 0, 'sRNA': 1, 'TSS_sRNA': 0, 'TSS_sRNA_RBS': 0, 'all': 1, 'RBS_sRNA': 1, 'TSS': 0, 'RBS': 1}, 'intergenic': {'TSS_RBS': 0, 'sRNA': 0, 'TSS_sRNA': 0, 'TSS_sRNA_RBS': 0, 'all': 0, 'RBS_sRNA': 0, 'TSS': 0, 'RBS': 0}}, 'total': {"5'UTR_derived": {'TSS_RBS': 0, 'sRNA': 1, 'TSS_sRNA': 0, 'TSS_sRNA_RBS': 0, 'all': 1, 'RBS_sRNA': 1, 'TSS': 0, 'RBS': 1}, "3'UTR_derived": {'TSS_RBS': 0, 'sRNA': 0, 'TSS_sRNA': 0, 'TSS_sRNA_RBS': 0, 'all': 0, 'RBS_sRNA': 0, 'TSS': 0, 'RBS': 0}, 'antisense': {'TSS_RBS': 0, 'sRNA': 0, 'TSS_sRNA': 0, 'TSS_sRNA_RBS': 0, 'all': 0, 'RBS_sRNA': 0, 'TSS': 0, 'RBS': 0}, 'interCDS': {'TSS_RBS': 0, 'sRNA': 0, 'TSS_sRNA': 0, 'TSS_sRNA_RBS': 0, 'all': 0, 'RBS_sRNA': 0, 'TSS': 0, 'RBS': 0}, 'all': {'TSS_RBS': 0, 'sRNA': 1, 'TSS_sRNA': 0, 'TSS_sRNA_RBS': 0, 'all': 1, 'RBS_sRNA': 1, 'TSS': 0, 'RBS': 1}, 'intergenic': {'TSS_RBS': 0, 'sRNA': 0, 'TSS_sRNA': 0, 'TSS_sRNA_RBS': 0, 'all': 0, 'RBS_sRNA': 0, 'TSS': 0, 'RBS': 0}}}
    sorf_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "sORF", "start": 5,
                "end": 8, "phase": ".", "strand": "+", "score": "."}]
    attributes_sorf = [{"ID": "sorf0", "Name": "sORF_0", "sORF_type": "5utr",
                        "with_TSS": "NA", "rbs": "RBS:1-2_+", "sRNA": "NA"}]
    sorfs = []
    sorfs.append(Create_generator(sorf_dict[0], attributes_sorf[0], "gff"))

if __name__ == "__main__":
    unittest.main()

