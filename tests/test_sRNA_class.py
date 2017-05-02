import sys
import os
import unittest
import shutil
import copy
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data, gen_file, extract_info
import annogesiclib.sRNA_class as sc
from mock_args_container import MockClass


class TestsRNAClass(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.mock_args = MockClass()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_initiate(self):
        out = StringIO()
        key = "test"
        key_list = ["test"]
        class_name = "test_name"
        class_num = 0
        index = {}
        content = "testtest"
        sc.initiate(key, key_list, class_name, class_num, index, out, content)
        self.assertEqual(out.getvalue(), "1testtest\n")

    def test_print_stat_title(self):
        out_stat = StringIO()
        strain = "aaa"
        checks = {"limit": False, "first": True, "utr": False, "inter": False}
        srna_datas = {"aaa": self.example.srnas, "all": self.example.srnas}
        args = self.mock_args.mock()
        args.energy = 0
        args.nr_hits_num = 0
        args.import_info = ["tss", "sec_str", "blast_nr", "blast_srna"]
        class_num, index = sc.print_stat_title(
            checks, out_stat, strain, srna_datas, 1, args)
        self.assertEqual(out_stat.getvalue(), """1 - Normalized(by length of sRNA) free energy change of the secondary structure is below to 0
2 - sRNA candidates start with TSS (3'UTR derived and interCDS sRNA also includes the sRNA candidates which start with processing site.)
3 - Running BLAST can not find the homology in nr database (the cutoff is 0).
4 - Running BLAST can not find the homology in sRNA database.
5 - Running BLAST can find the homology in sRNA database.
All genomes:
""")

        self.assertEqual(class_num, 5)
        self.assertDictEqual(index, {
            'sRNA_hit': 5, '2d_energy': 1, 'sRNA_no_hit': 4,
            'nr_no_hit': 3, 'with_TSS': 2})

    def test_import_class(self):
        index = {'sRNA_hit': 5, '2d_energy': 1, 'sRNA_no_hit': 4,
                 'nr_no_hit': 3, 'with_TSS': 2}
        num_srna = 0
        datas_srna = {}
        datas = {"aaa": self.example.srnas}
        num = sc.import_class(5, datas_srna, datas, index, num_srna, "aaa",
                              "UTR_derived", "5utr", 0, 0)
        self.assertEqual(num, 1)
        self.assertEqual(datas_srna["class_4"][0].start, 230)

    def test_import_data(self):
        datas = {"aaa": self.example.srnas, "all": self.example.srnas}
        index = {'sRNA_hit': 5, '2d_energy': 1, 'sRNA_no_hit': 4,
                 'nr_no_hit': 3, 'with_TSS': 2}
        num_srna = {"total": 0, "intergenic": 0, "5'UTR_derived": 0,
                    "3'UTR_derived": 0, "interCDS": 0, "in_CDS": 0}
        checks = {"limit": False, "first": True, "utr": True,
                  "inter": True, "in_CDS": True, "antisense": False}
        datas_rna = sc.import_data(5, datas, index, num_srna,
                                   "aaa", checks, 0, 0)
        self.assertEqual(datas_rna["5'UTR_derived"]["class_4"][0].start, 230)
        self.assertEqual(datas_rna["interCDS"]["class_1"][0].start, 140)
        self.assertEqual(datas_rna["in_CDS"]["class_1"][0].start, 6166)
        self.assertEqual(datas_rna["intergenic"]["class_5"][0].start, 5166)

    def test_print_intersection(self):
        num_srna = {"total": 3, "intergenic": 1, "5'UTR_derived": 1,
                    "3'UTR_derived": 0, "interCDS": 1}
        gff_name = os.path.join(self.test_folder, "test")
        out_stat = StringIO()
        keys = ["class_1", "class_4", "class_2", "class_3", "class_5"]
        datas = {"class_1": self.example.srnas, "class_2": self.example.srnas,
                 "class_3": self.example.srnas, "class_4": self.example.srnas,
                 "class_5": self.example.srnas}
        sc.print_intersection(datas, keys, 3, gff_name, "total", out_stat)
        self.assertEqual(out_stat.getvalue(),
                         ("\tclass_1 and class_4 and class_2 and class_3 "
                          "and class_5 = 4(1.3333333333333333)\n"))
        results, attributes = extract_info(os.path.join(
            self.test_folder, "test"), "file")
        self.assertEqual("\n".join(results), self.example.gff_info)

    def test_read_file(self):
        srna_file = os.path.join(self.test_folder, "srna.gff")
        gen_file(srna_file, self.example.gff_file)
        srna_datas, strains, checks = sc.read_file(srna_file)
        self.assertEqual(srna_datas["aaa"][0].start, 140)
        self.assertEqual(srna_datas["aaa"][1].start, 230)
        self.assertEqual(srna_datas["bbb"][0].start, 5166)
        self.assertListEqual(strains, ['all', 'aaa', 'bbb'])
        self.assertDictEqual(checks, {'inter': True, 'limit': False,
                                      'utr': True, 'antisense': False,
                                      'in_CDS': True, 'first': True})

    def test_sort_keys(self):
        keys = ["class_3", "class_1", "class_5"]
        final_keys = sc.sort_keys(keys)
        self.assertListEqual(final_keys, ['class_1', 'class_3', 'class_5'])

    def test_classify_srna(self):
        out_stat_file = os.path.join(self.test_folder, "stat")
        srna_file = os.path.join(self.test_folder, "srna.gff")
        gen_file(srna_file, self.example.gff_file)
        args = self.mock_args.mock()
        args.energy = 0
        args.nr_hits_num = 0
        args.in_cds = True
        args.import_info = ["tss", "sec_str"]
        sc.classify_srna(srna_file, self.test_folder, out_stat_file, args)


class Example(object):
    gff_file = """aaa	UTR_derived	sRNA	140	160	.	+	.	ID=srna0;Name=sRNA_0;sRNA_type=interCDS
aaa	UTR_derived	sRNA	230	280	.	+	.	ID=srna1;Name=sRNA_1;sRNA_type=5utr
bbb	intergenic	sRNA	5166	5266	.	-	.	ID=srna2;Name=sRNA_2;sRNA_type=intergenic
bbb	in_CDS	sRNA	6166	6499	.	-	.	ID=srna3;Name=sRNA_3;sRNA_type=in_CDS"""
    gff_info = """aaa	UTR_derived	sRNA	140	160	.	+	.
aaa	UTR_derived	sRNA	230	280	.	+	.
bbb	intergenic	sRNA	5166	5266	.	-	.
bbb	in_CDS	sRNA	6166	6499	.	-	."""
    srna_dict = [{"seq_id": "aaa", "source": "UTR_derived",
                  "feature": "sRNA", "start": 140,
                  "end": 160, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "aaa", "source": "UTR_derived",
                  "feature": "sRNA", "start": 230,
                  "end": 280, "phase": ".", "strand": "+", "score": "."},
                 {"seq_id": "bbb", "source": "intergenic",
                  "feature": "sRNA", "start": 5166,
                  "end": 5266, "phase": ".", "strand": "-", "score": "."},
                 {"seq_id": "bbb", "source": "in_CDS",
                  "feature": "sRNA", "start": 6166,
                  "end": 6499, "phase": ".", "strand": "-", "score": "."}]
    attributes_srna = [
        {"ID": "srna0", "Name": "sRNA_0", "with_TSS": "TSS:100_+",
         "with_cleavage": "Cleavage:145_+", "sRNA_type": "interCDS",
         "sRNA_hit": "ReaA", "2d_energy": "-0.45", "nr_hit": "NA"},
        {"ID": "srna1", "Name": "sRNA_1", "with_TSS": "TSS:200_+",
         "with_cleavage": "NA", "sRNA_type": "5utr", "sRNA_hit": "NA",
         "2d_energy": "0", "nr_hit": "3"},
        {"ID": "srna2", "Name": "sRNA_2", "with_TSS": "NA",
         "with_cleavage": "NA", "sRNA_type": "intergenic", "sRNA_hit": "sRAE",
         "2d_energy": "-0.032", "nr_hit": "NA"},
        {"ID": "srna3", "Name": "sRNA_3", "with_TSS": "TSS:6499_-",
         "with_cleavage": "NA", "sRNA_type": "in_CDS", "sRNA_hit": "NA",
         "2d_energy": "-0.042", "nr_hit": "NA"}]
    srnas = []
    for index in range(0, 4):
        srnas.append(Create_generator(
            srna_dict[index], attributes_srna[index], "gff"))

if __name__ == "__main__":
    unittest.main()

