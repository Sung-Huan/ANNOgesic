import sys
import math
import csv
import os
import shutil
import unittest
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data, gen_file, extract_info
import annogesiclib.validate_gene as vg

#class Mock_func(object):


class TestTranscriptAssembly(unittest.TestCase):

    def setUp(self):
        self.example = Example()
#        self.mock = Mock_func()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_read_gff(self):
        gff_file = os.path.join(self.test_folder, "test.gff")
        tss_file = os.path.join(self.test_folder, "test_TSS.gff")
        gen_file(gff_file, self.example.gff_file)
        gen_file(tss_file, self.example.tss_file)
        gffs, tsss = vg.read_gff(gff_file, tss_file)
        self.assertEqual(gffs[0].start, 5)
        self.assertEqual(tsss[0].start, 3)

    def test_compare_tss(self):
        num_all = {"all_cds": 0, "all_tRNA": 0, "all_rRNA": 0,
                   "cds": 0, "tRNA": 0, "rRNA": 0}
        num_strain = {"test": {"all_cds": 0, "all_tRNA": 0,
                               "all_rRNA": 0, "cds": 0,
                               "tRNA": 0, "rRNA": 0}}
        vg.compare_tss(self.example.tsss, self.example.gffs[0], 300, num_all, num_strain)
        self.assertDictEqual(num_all, {'tRNA': 0, 'rRNA': 0, 'all_tRNA': 0,
                                       'cds': 1, 'all_cds': 0, 'all_rRNA': 0})
        self.assertDictEqual(num_strain, {'test': {'tRNA': 0, 'rRNA': 0, 'all_tRNA': 0,
                                                   'cds': 1, 'all_cds': 0, 'all_rRNA': 0}})

    def test_print_stat(self):
        num = {"all_cds": 300, "all_tRNA": 30, "all_rRNA": 20,
               "cds": 100, "tRNA": 10, "rRNA": 10}
        out = StringIO()
        vg.print_stat("cds", num, out)
        self.assertEqual(out.getvalue(), "The number of cds which is start from TSS: 100 (0.3333333333333333)\n")

    def test_print_file(self):
        num_all = {"all_cds": 600, "all_tRNA": 30, "all_rRNA": 30,
                   "cds": 250, "tRNA": 20, "rRNA": 20}
        num_strain = {"test": {"all_cds": 300, "all_tRNA": 20,
                               "all_rRNA": 20, "cds": 100,
                               "tRNA": 10, "rRNA": 10}}
        out_cds_file = os.path.join(self.test_folder, "cds_file")
        stat_file = os.path.join(self.test_folder, "stat_file")
        vg.print_file(self.example.gffs, out_cds_file, stat_file, num_all, num_strain)
        datas, attribute = extract_info(out_cds_file, "file")
        self.assertEqual("\n".join(datas), "test\tRefSeq\tCDS\t200\t270\t.\t+\t.")
        datas = import_data(stat_file)
        self.assertEqual("\n".join(datas), self.example.out_stat_test)

    def test_validate_gff(self):
        gff_file = os.path.join(self.test_folder, "test.gff")
        tss_file = os.path.join(self.test_folder, "test_TSS.gff")
        gen_file(gff_file, self.example.gff_file)
        gen_file(tss_file, self.example.tss_file)
        out_cds_file = os.path.join(self.test_folder, "cds_file")
        stat_file = os.path.join(self.test_folder, "stat_file")
        vg.validate_gff(tss_file, gff_file, stat_file, out_cds_file, 300)
        datas, attribute = extract_info(out_cds_file, "file")
        self.assertEqual("\n".join(datas), "test\tRefSeq\tCDS\t5\t10\t.\t+\t.")
        datas = import_data(stat_file)
        self.assertEqual("\n".join(datas), self.example.out_stat)



class Example(object):

    gff_file = """test	RefSeq	CDS	5	10	.	+	.	ID=cds0;Name=CDS_0"""
    tss_file = """test	RefSeq	TSS	3	3	.	+	.	ID=tss0;Name=TSS_0"""
    tss_dict = [{"seq_id": "test", "source": "intergenic", "feature": "TSS", "start": 170,
                "end": 170, "phase": ".", "strand": "+", "score": "."}]
    attributes_tsss = [{"ID": "tss0", "Name": "TSS_0"}]
    tsss = []
    tsss.append(Create_generator(tss_dict[0], attributes_tsss[0], "gff"))
    gff_dict = [{"seq_id": "test", "source": "RefSeq", "feature": "CDS", "start": 200,
                "end": 270, "phase": ".", "strand": "+", "score": "."}]
    attributes_gff = [{"ID": "cds0", "Name": "CDS_0"}]
    gffs = []
    gffs.append(Create_generator(gff_dict[0], attributes_gff[0], "gff"))
    out_stat_test = """All strains:
The number of cds which is start from TSS: 250 (0.4166666666666667)
The number of tRNA which is start from TSS: 20 (0.6666666666666666)
The number of rRNA which is start from TSS: 20 (0.6666666666666666)"""
    out_stat = """All strains:
The number of cds which is start from TSS: 1 (1.0)
The number of tRNA which is start from TSS: 0 (NA)
The number of rRNA which is start from TSS: 0 (NA)"""


if __name__ == "__main__":
    unittest.main()
