import sys
import os
import unittest
import shutil
import copy
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data, gen_file, extract_info
import annogesiclib.stat_term as st


class TestStatTerm(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_plus_num(self):
        nums = {"aaa": {"test_1": 0}, "total": {"test_1": 2}}
        st.plus_num(nums, "aaa", "test_1")
        self.assertDictEqual(nums, {'total': {'test_1': 3}, 'aaa': {'test_1': 1}})

    def test_print_percent(self):
        out = StringIO()
        st.print_percent(out, 10, 2, "test")
        self.assertEqual(out.getvalue(), "\t\t(percentage of total testterminators = 0.2)\n")

    def test_print_express(self):
        out = StringIO()
        st.print_express(out, 10, 2, "test")
        self.assertEqual(out.getvalue(), "\t\t(percentage of total testterminators which have gene expression = 0.2)\n")

    def test_print_decrease(self):
        out = StringIO()
        st.print_decrease(out, 10, 2, "test")
        self.assertEqual(out.getvalue(), "\t\t(percentage of total testterminators which have dramatic coverage decreasing = 0.2)\n")

    def test_print_method(self):
        out = StringIO()
        nums = {"method": 2, "express": 3, "detect": 4, "only": 5,
                "total": 20, "total_ex": 15, "total_de": 10, "frhp": 1}
        st.print_method(nums, "method_1", "method", "express", "detect", "only", out)
        self.assertEqual(out.getvalue(), self.example.print_method + "\n")

    def test_print_intersection_number(self):
        out = StringIO()
        nums = {"test1": 3, "total": 10, "fr": 2, "hp":1}
        st.print_intersection_number(out, nums, "test1")
        self.assertEqual(out.getvalue(), self.example.print_inter_num)

    def test_print_intersection_express(self):
        out = StringIO()
        nums = {"test1": 2, "total_ex": 10, "ex_fr": 3, "ex_hp": 4}
        st.print_intersection_express(out, nums, "test1")
        self.assertEqual(out.getvalue(), self.example.print_inter_express)

    def test_print_file(self):
        out = StringIO()
        nums = {"fr": 2, "hp": 3, "frhp": 4, "ex_fr": 3, "ex_hp": 2, "ex_frhp": 1,
                "de_fr": 2, "de_hp": 3, "de_frhp": 4, "total": 20,
                "total_de": 10, "total_ex": 10, "only_de_fr": 6, "only_de_hp": 4,
                "only_ex_fr": 2, "only_ex_hp": 1, "de_frhp": 2, "ex_frhp": 3}
        st.print_file(nums, out, "aaa")
        self.assertEqual(out.getvalue().split("\n\n")[0], self.example.print_file)

    def test_stat_term(self):
        term_gff = os.path.join(self.test_folder, "aaa_term.gff")
        term_table = os.path.join(self.test_folder, "aaa_term.csv")
        stat = os.path.join(self.test_folder, "stat")
        output_decrease = os.path.join(self.test_folder, "decrease")
        output_expression = os.path.join(self.test_folder, "expression")
        gen_file(term_gff, self.example.gff)
        gen_file(term_table, self.example.table)
        st.stat_term(term_gff, term_table, stat, output_decrease, output_expression)
        self.assertTrue(stat)
        datas = import_data(output_decrease + ".csv")
        self.assertEqual("\n".join(datas), self.example.table)
        datas = import_data(output_expression + ".csv")
        self.assertEqual("\n".join(datas), self.example.table)


class Example(object):

    table = """Staphylococcus_aureus_HG003	Term_00000	9676	9705	+	True	ID-001873-Staph_aureus_sample_mix_forward(diff=19.0;high=33.0;low=14.0)
Staphylococcus_aureus_HG003	Term_00006	23768	23793	+	True	pMEM_OD_1_forward(diff=8.637108857702653;high=10.796386072128318;low=2.159277214425664)"""
    gff = """Staphylococcus_aureus_HG003	ANNOgesic	terminator	9676	9705	.	+	.	Name=Terminator_00000;express=True;diff_coverage=ID-001873-Staph_aureus_sample_mix_forward(high:33.0,low:14.0);ID=term0;associate=SAOUHSC_00006;coverage_decrease=True;Method=TransTermHP
Staphylococcus_aureus_HG003	ANNOgesic	terminator	23768	23793	.	+	.	Name=Terminator_00006;express=True;diff_coverage=pMEM_OD_1_forward(high:10.796386072128318,low:2.159277214425664);ID=term6;associate=SAOUHSC_00019;coverage_decrease=True;Method=TransTermHP"""
    print_method = """method_1
	Total method_1 terminators = 3
		(percentage of total terminators = 0.15)
	Total terminators which only can be detected in method_1 = 2
		(percentage of total terminators = 0.1)
		(percentage of total method_1 terminators = 0.6666666666666666)
	Total method_1 terminators which located in gene expression region = 3
		(percentage of total terminators = 0.15)
		(percentage of total method_1 terminators = 1.0)
	Total method_1 terminators which have dramatic coverage decreasing = 4
		(percentage of total terminators = 0.2)
		(percentage of total method_1 terminators = 1.3333333333333333)
		(percentage of total terminators which have gene expression = 0.26666666666666666)
		(percentage of total method_1 terminators which have gene expression = 1.3333333333333333)
	Total terminators which have dramatic coverage decreasing(unique in method_1) = 5
		(percentage of total terminators which have dramatic coverage decreasing = 0.5)
		(percentage of total method_1 terminators which have dramatic coverage decreasing = 1.25)
"""
    print_inter_num = """		(percentage of total terminators = 0.3)
		(percentage of total method_1 terminators = 1.5)
		(percentage of total method_2 terminators = 3.0)
"""
    print_inter_express = """		(percentage of total terminators which have gene expression = 0.2)
		(percentage of total method_1 terminators which have gene expression = 0.6666666666666666)
		(percentage of total method_2 terminators which have gene expression = 0.5)
"""
    print_file = """aaa:
Combine two methods:
	Total terminators = 20
	Total terminators which located in gene expression region = 10
		(percentage of total terminators = 0.5)
	Total terminators which have dramatic coverage decreasing = 10
		(percentage of total terminators = 0.5)
		(percentage of total terminators which have gene expression = 1.0)"""

if __name__ == "__main__":
    unittest.main()

