import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import gen_file
import annogesiclib.plot_PPI as pp


class Mock_func(object):

    def mock_add_node(self, G, nodes):
        pass

    def mock_add_edge(self, G, ppi, style, weight, color):
        pass

    def mock_nx_node(self, G, pos, node_size, colors, color_list):
        pass

    def mock_nx_edge(self, G, pos, connects, colors, styles, weights):
        return None

    def mock_nx_label(self, G, pos, labels, size):
        pass

    def mock_nx_color_style(self, G, edges):
        return None, None, None

    def mock_plot_text(self, check_na, plt, ppis, pre_ppi, color_edge):
        pass
        
    def mock_plot(self, ppis, center, strain, cutoff_score,
                  node_size, out_folder):
        pass

class TestPlotPPI(unittest.TestCase):

    def setUp(self):
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
        self.mock = Mock_func()
        self.example = Example()

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_node(self):
        nodes = []
        item = "test"
        labels1 = {}
        labels2= {}
        colors = {}
        center = {"locus_tag": "test"}
        pp.node(item, nodes, center, colors, labels1, labels2)
        self.assertListEqual(nodes, ["test"])
        self.assertDictEqual(colors, {'test': '#FFFF66'})
        self.assertDictEqual(labels1, {'test': 'test'})
        self.assertDictEqual(labels2, {'test': ''})

    def test_best_assign_attributes(self):
        pp.add_edge = self.mock.mock_add_edge
        pp.add_node = self.mock.mock_add_node
        check_na = {"best": False, "same_best": False}
        ppi = {"score": 30, "below": 10, "best": 22}
        pre_ppi = {"score": 20, "below": 20, "best": 21}
        first = False
        pp.best_assign_attributes(check_na, None, ppi, pre_ppi, first, "test")
        self.assertDictEqual(ppi, {'score': 30, 'best': 22, 'below': 10})
        self.assertDictEqual(check_na, {'best': True, 'same_best': True})

    def test_create_node(self):
        pp.add_edge = self.mock.mock_add_edge
        pp.add_node = self.mock.mock_add_node
        check_na = {"best": False, "same_best": False}
        pre_ppi = {"score": 20, "below": 20, "best": 0.21}
        ppis = [{"item_a": "test_1", "item_b": "test_2",
                 "score": 30, "below": 10, "best": 0.22},
                {"item_a": "test_3", "item_b": "test_4",
                 "score": 20, "below": 10, "best": -0.21}]
        labels1 = {}
        labels2= {}
        nodes = []
        colors = {}
        scores = []
        edges = []
        center = {"locus_tag": "test", "gene_name": "test_1"}
        pp.create_node(ppis, scores, nodes, center, colors, labels1, labels2,
                       edges, None, 0, check_na, pre_ppi)
        self.assertListEqual(nodes, ['test_1', 'test_2', 'test_3', 'test_4'])
        self.assertListEqual(scores, [30, 20])
        self.assertDictEqual(colors, {
            'test_1': '#FFFF66', 'test_4': '#CCFFCC',
            'test_2': '#CCFFCC', 'test_3': '#CCFFCC'})
        self.assertDictEqual(labels1, {'test_1': '', 'test_4': '',
                                       'test_2': '', 'test_3': ''})
        self.assertDictEqual(labels2, {'test_1': 'test_1', 'test_4': 'test_4',
                                       'test_2': 'test_2', 'test_3': 'test_3'})
        self.assertListEqual(edges, [('test_1', 'test_2'),
                                     ('test_3', 'test_4')])
        self.assertDictEqual(check_na, {'same_best': True, 'best': True})

    def test_modify_label(self):
        labels2 = {'test_1': 'test1', 'test_4': 'test4',
                   'test_2': 'test_2', 'test_3': 'test_3'}
        new_labels = {'test_1': "", 'test_4': "",
                      'test_2': "", 'test_3': ""}
        pp.modify_label(labels2, new_labels)
        self.assertDictEqual(new_labels, {
            'test_2': 'test\n2', 'test_4': 'test4',
            'test_1': 'test1', 'test_3': 'test\n3'})

    def test_plot(self):
        pp.plot_text = self.mock.mock_plot_text
        pp.nx_node = self.mock.mock_nx_node
        pp.nx_edge = self.mock.mock_nx_edge
        pp.nx_label = self.mock.mock_nx_label
        pp.add_edge = self.mock.mock_add_edge
        pp.add_node = self.mock.mock_add_node
        pp.nx_color_style = self.mock.mock_nx_color_style
        ppis = [{"item_a": "test_1", "item_b": "test_2", "score": 30,
                 "below": 10, "best": 0.22},
                {"item_a": "test_3", "item_b": "test_4", "score": 20,
                 "below": 10, "best": -0.21}]
        center = {"locus_tag": "test", "gene_name": "test_1"} 
        pp.plot(ppis, center, "test", 0, 1000, self.test_folder)
        self.assertTrue(os.path.exists(os.path.join(
            self.test_folder, "test", "test_test_1.png")))

    def test_score_compare(self):
        scores = {"score": 0, "below": 0}
        score = 11
        ppi = {"score": 6, "below": 6}
        pp.score_compare(score, scores, 0, ppi)
        self.assertDictEqual(scores, {'score': 1, 'below': 0})
        self.assertDictEqual(ppi, {'score': 6, 'below': 6})
        score = "NA"
        scores = {"score": 0, "below": 0}
        pp.score_compare(score, scores, 0, ppi)
        self.assertDictEqual(ppi, {'score': 0, 'below': 0})
        self.assertDictEqual(scores, {'score': 0, 'below': 0})

    def test_assign_score_below(self):
        ppis = []
        scores = {"score": 3, "below": 3}
        pre_ppi = {}
        pp.assign_score_below(pre_ppi, scores, ppis)
        self.assertListEqual(ppis, [{'score': 3, 'below': 3}])
        self.assertDictEqual(pre_ppi, {'score': 3, 'below': 3})

    def test_get_best(self):
        ppi = {}
        pre_ppi = {}
        row = ["", "", "", "", "", "", "", "", 30]
        pp.get_best(pre_ppi, ppi, row)
        self.assertDictEqual(ppi, {'best': 30})
        self.assertDictEqual(pre_ppi, {})
        pre_ppi = {"best": "NA"}
        ppi = {}
        pp.get_best(pre_ppi, ppi, row)
        self.assertDictEqual(ppi, {'best': 30})
        self.assertDictEqual(pre_ppi, {"best": "NA"})
        pre_ppi = {"best": 10}
        ppi = {}
        pp.get_best(pre_ppi, ppi, row)
        self.assertDictEqual(ppi, {'best': 30})
        self.assertDictEqual(pre_ppi, {"best": 10})
        pre_ppi = {"best": 50}
        ppi = {}
        pp.get_best(pre_ppi, ppi, row)
        self.assertDictEqual(ppi, {'best': 50})
        self.assertDictEqual(pre_ppi, {"best": 50})

    def test_plot_ppi(self):
        ppi_file = os.path.join(self.test_folder, "test_ppi")
        gen_file(ppi_file, self.example.ppi_file)
        pp.plot_ppi(ppi_file, 0, self.test_folder, 1000)
        strain = "Helicobacter pylori 26695 chromosome"
        self.assertTrue(os.path.exists(
            "test_folder/" + strain + "/HP0001_nusB.png"))

class Example(object):

    ppi_file = """Interaction of HP0001 | nusB
strain	item_id_a	item_id_b	mode	action	a_is_acting	STRING_action_score	pubmed_id	pubmed_score
Helicobacter pylori 26695 chromosome	nusG	rpoB	reaction;reaction;catalysis;catalysis;binding		1;0;0;1;0	977	NA	NA
Helicobacter pylori 26695 chromosome	ribD	ribH	binding		0	608	NA	NA
Helicobacter pylori 26695 chromosome	rpsJ	nusB	binding		0	810	NA	NA"""

if __name__ == "__main__":
    unittest.main()

