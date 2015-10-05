import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file, import_data
from annogesiclib.ppi import PPINetwork


class Mock_func(object):

    def mock_run_wget(self, source, folder, log):
        pass

class TestPPI(unittest.TestCase):

    def setUp(self):
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
            os.mkdir(os.path.join(self.test_folder, "tmp_specific"))
            os.mkdir(os.path.join(self.test_folder, "tmp_nospecific"))
            os.mkdir(os.path.join(self.test_folder, "with_strain"))
            os.mkdir(os.path.join(self.test_folder, "with_strain/test_ptt"))
            os.mkdir(os.path.join(self.test_folder, "without_strain"))
            os.mkdir(os.path.join(self.test_folder, "without_strain/test_ptt"))
            os.mkdir(os.path.join(self.test_folder, "all_results"))
            os.mkdir(os.path.join(self.test_folder, "best_results"))
            os.mkdir(os.path.join(self.test_folder, "figures"))
        self.ppi = PPINetwork(self.test_folder)
        self.mock = Mock_func()
        self.example = Example()

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_wget_id(self):
        self.ppi._run_wget = self.mock.mock_run_wget
        strain = "test_strain"
        locus = "test_locus"
        strain_id = {"ptt": "test_strain", "string": "string_test", "file": "file_test"}
        files = {"id_list": "test", "id_log": "test"}
        detect = self.ppi._wget_id(strain, locus, strain_id, files)
        self.assertTrue(detect)

    def test_retrieve_id(self):
        self.ppi._run_wget = self.mock.mock_run_wget
        strain_id = {"ptt": "test_strain", "string": "string_test", "file": "file_test"}
        files = {"id_list": "test", "id_log": "test"}
        genes = [{"strain": "test_strain", "locus_tag": "test_locus"}]
        self.ppi._retrieve_id(strain_id, genes, files)

    def test_get_prefer_name(self):
        row_a = "999.aaa"
        files = {"id_list": self.test_folder}
        gen_file(os.path.join(self.test_folder, "aaa"), "999.aaa\t222\t333\ttest_aaa")
        name = self.ppi._get_prefer_name(row_a, "test", files, "test")
        self.assertEqual(name, "test_aaa")

    def test_get_pubmed(self):
        out_all = StringIO()
        out_best = StringIO()
        out_noall = StringIO()
        out_nobest = StringIO()
        self.ppi._run_wget = self.mock.mock_run_wget
        files = {"id_list": self.test_folder, "id_log": "test", "pubmed_log": "test",
                 "all_specific": out_all, "best_specific": out_best,
                 "all_nospecific": out_noall, "best_nospecific": out_nobest}
        row = self.example.ppi_line.split("\t")
        strain_id = {"file": "test_file","ptt": "test_ptt", "string": "test_string", "pie": "test_pie"}
        mode = "interaction"
        actor = "test_A"
        score = 11241
        id_file = "SAOUHSC_01684"
        ptt = "test_ptt"
        gen_file(os.path.join(self.test_folder, "SAOUHSC_01684"), "93061.SAOUHSC_01684\t93061.SAOUHSC_01683\t333\ttest_aaa")
        gen_file(os.path.join(self.test_folder, "SAOUHSC_01683"), "93061.SAOUHSC_01683\t93061.SAOUHSC_01684\t333\ttest_bbb")
        paths = {"all": self.test_folder, "fig": self.test_folder, "best": self.test_folder}
        querys = "all"
        first_output = {"specific_all": True, "specific_best": True,
                        "nospecific_all": True, "nospecific_best": True}
        self.ppi._get_pubmed(row, self.test_folder, strain_id, mode, actor, score,
                    id_file, first_output, True, ptt, files, paths, querys)
        data = import_data("test_folder/without_strain/test_ptt/test_aaa_test_bbb.csv")
        self.assertEqual("\n".join(data), self.example.with_out)
        data = import_data("test_folder/with_strain/test_ptt/test_aaa_test_bbb.csv")
        self.assertEqual("\n".join(data), self.example.with_out)

    def test_merge_information(self):
        first_output = {"specific_all": True, "specific_best": True,
                        "nospecific_all": True, "nospecific_best": True}
        out_all = StringIO()
        out_best = StringIO()
        row_a = self.example.ppi_line.split("\t")
        score = 111
        id_file = "SAOUHSC_01684"
        id_folder = self.test_folder
        file_type = "specific"
        all_folder = os.path.join(self.test_folder, "with_strain")
        best_folder = os.path.join(self.test_folder, "without_strain")
        ptt = "test_ptt"
        filename = os.path.join(self.test_folder, "SAOUHSC_01684")
        gen_file(filename, "93061.SAOUHSC_01684\t1000\t333\ttest_aaa")
        self.ppi._merge_information(first_output, filename, out_all, out_best,
                           row_a, score, id_file, id_folder, file_type,
                           all_folder, best_folder, ptt)
        self.assertEqual(out_all.getvalue(), self.example.merge_out + "\n")
        self.assertEqual(out_best.getvalue(), self.example.merge_out + "\n")

    def test_detect_protein(self):
        gen_file(os.path.join(self.test_folder, "test"), self.example.ptt_file)
        strain_id = {"file": "test","ptt": "test_ptt", "string": "test_string", "pie": "test_pie"}
        genes = self.ppi._detect_protein(self.test_folder, strain_id, "all")
        self.assertListEqual(genes, [{'strain': 'Staphylococcus_aureus_HG003', 'locus_tag': 'SAOUHSC_00001'},
                                     {'strain': 'Staphylococcus_aureus_HG003', 'locus_tag': 'SAOUHSC_00002'},
                                     {'strain': 'Staphylococcus_aureus_HG003', 'locus_tag': 'SAOUHSC_00003'}])

    def test_setup_nospecific(self):
        out_all = StringIO()
        out_best = StringIO()
        out_noall = StringIO()
        out_nobest = StringIO()
        paths = {"all": os.path.join(self.test_folder, "all_results"),
                 "fig": os.path.join(self.test_folder, "figures"),
                 "best": os.path.join(self.test_folder, "best_results")}
        strain_id = {"file": "test","ptt": "test_ptt", "string": "test_string", "pie": "test_pie"}
        files = {"id_list": self.test_folder, "id_log": "test", "pubmed_log": "test",
                 "all_specific": out_all, "best_specific": out_best,
                 "all_nospecific": out_noall, "best_nospecific": out_nobest}
        self.ppi._setup_nospecific(paths, strain_id, files)
        files["all_nospecific"].close()
        files["best_nospecific"].close()
        self.assertTrue(os.path.exists("test_folder/all_results/without_strain/test_ptt"))
        self.assertTrue(os.path.exists("test_folder/best_results/without_strain/test_ptt"))
        self.assertTrue(os.path.exists("test_folder/figures/without_strain/test_ptt"))

    def test_setup_folder_and_read_file(self):
        paths = {"all": os.path.join(self.test_folder, "all_results"),
                 "fig": os.path.join(self.test_folder, "figures"),
                 "best": os.path.join(self.test_folder, "best_results")}
        strain_id = {"file": "test.ptt","ptt": "test_ptt", "string": "test_string", "pie": "test_pie"}
        files = {"id_list": self.test_folder, "id_log": "", "pubmed_log": "",
                 "all_specific": "", "best_specific": "",
                 "all_nospecific": "", "best_nospecific": "", "action_log": ""}
        gen_file(os.path.join(self.test_folder, "test.ptt"), self.example.ptt_file)
        genes = self.ppi._setup_folder_and_read_file(strain_id, "", self.test_folder,
                                                     self.test_folder, True, files, paths, "all")
        for index in ("all_specific", "all_nospecific", "best_specific", "best_nospecific",
                      "id_log", "action_log", "pubmed_log"):
            files[index].close()
        self.assertTrue(os.path.exists("test_folder/best_results/test"))
        self.assertTrue(os.path.exists("test_folder/all_results/test"))
        self.assertListEqual(genes, [{'locus_tag': 'SAOUHSC_00001', 'strain': 'Staphylococcus_aureus_HG003'},
                                     {'locus_tag': 'SAOUHSC_00002', 'strain': 'Staphylococcus_aureus_HG003'},
                                     {'locus_tag': 'SAOUHSC_00003', 'strain': 'Staphylococcus_aureus_HG003'}])

    def test_wget_actions(self):
        gen_file(os.path.join(self.test_folder, "test.txt"), "93061\ttest")
        self.ppi._run_wget = self.mock.mock_run_wget
        files = {"id_list": self.test_folder, "id_log": "", "pubmed_log": "",
                 "all_specific": "", "best_specific": "",
                 "all_nospecific": "", "best_nospecific": "", "action_log": ""}
        strain_id = {"file": "test.ptt","ptt": "test_ptt", "string": "test_string", "pie": "test_pie"}
        id_file = "test.txt"
        self.ppi._wget_actions(files, id_file, strain_id, self.test_folder)

    def test_retrieve_actions(self):
        self.ppi._run_wget = self.mock.mock_run_wget
        files = {"id_list": os.path.join(self.test_folder, "tmp_specific"), "id_log": "", "pubmed_log": "",
                 "all_specific": "", "best_specific": "",
                 "all_nospecific": "", "best_nospecific": "", "action_log": ""}
        strain_id = {"file": "test.ptt","ptt": "test_ptt", "string": "test_string", "pie": "test_pie"}
        paths = {"all": os.path.join(self.test_folder, "all_results"),
                 "fig": os.path.join(self.test_folder, "figures"),
                 "best": os.path.join(self.test_folder, "best_results")}
        gen_file(os.path.join(self.test_folder, "tmp_specific/test.txt"), "93061\ttest")
        gen_file(os.path.join(self.test_folder, "tmp_action"), self.example.ppi_line)
        self.ppi._retrieve_actions(files, self.test_folder, strain_id, paths,
                                   100, True, "all")

class Example(object):

    ppi_line = "93061.SAOUHSC_01684	93061.SAOUHSC_01683	catalysis		0	861"
    ppi_file = """93061.SAOUHSC_01684	93061.SAOUHSC_01683	catalysis		0	861
93061.SAOUHSC_01684	93061.SAOUHSC_01683	catalysis		1	861
93061.SAOUHSC_01684	93061.SAOUHSC_01683	binding		0	965"""
    with_out = """Interaction of SAOUHSC_01684 | test_aaa
strain	item_id_a	item_id_b	mode	action	a_is_acting	STRING_action_score	pubmed_id	pubmed_score
test_ptt	test_aaa	test_bbb	interaction		test_A	861	NA	NA"""

    merge_out = """Interaction of SAOUHSC_01684 | test_aaa
strain	item_id_a	item_id_b	mode	action	a_is_acting	STRING_action_score	pubmed_id	pubmed_score
test_ptt	93061.SAOUHSC_01684	93061.SAOUHSC_01683	catalysis		0	861	93061.SAOUHSC_01684	1000	333	test_aaa"""
    ptt_file = """Staphylococcus_aureus_HG003 - 1..2821337
2772 proteins
Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
517..1878	+	1362	YP_498609.1	dnaA	SAOUHSC_00001	-	-	chromosomal replication initiation protein
2156..3289	+	1134	YP_498610.1	-	SAOUHSC_00002	-	-	DNA polymerase III subunit beta
3670..3915	+	246	YP_498611.1	-	SAOUHSC_00003	-	-	hypothetical protein"""

if __name__ == "__main__":
    unittest.main()

