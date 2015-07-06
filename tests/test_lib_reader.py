import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
import annogesiclib.lib_reader


class Mock_parser_wig(object):

    def __init__(self):
        self.tracks = ["TSB_OD_0.2_TEX_reverse",
                       "TSB_OD_0.2_TEX_forward",
                       "TSB_OD_0.2_reverse",
                       "TSB_OD_0.2_forward",
                       "frag_forward",
                       "frag_reverse"]
        self.wigs = [{"pos": 312, "coverage":1.4041251228308191},
                     {"pos": 313, "coverage":56.867067474648174},
                     {"pos": 314, "coverage":56.867067474648174},
                     {"pos": 315, "coverage":56.867067474648174}]

    def parser(self, input_file, strand):
        for track in self.tracks:
            for wig in self.wigs:
                yield Mock_assign_value(
                      wig["pos"], wig["coverage"], "+", "aaa", track)

class Mock_assign_value(object):
    def __init__(self, pos, coverage, strand, strain, track):
        self.pos = int(pos)
        if strand == "+":
            self.coverage = float(coverage)
        else:
            self.coverage = -1 * float(coverage)
        self.strand = strand
        self.strain = strain
        self.track = track

class TestLibReader(unittest.TestCase):

    def setUp(self):
        self.lib_reader = annogesiclib.lib_reader
        self.example = Example()
        self.lib_reader.WigParser = Mock_parser_wig 
        self.libs = self.example.libs
        self.wigs = self.example.wigs
        self.lib_dict = self.example.lib_dict
        self.wig_file = self.example.wig_file
        self.libs_for_wig = self.example.libs_for_wig
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
        for track, wig in self.wigs.items():
            wig_file = os.path.join(self.test_folder, track)
            with open(wig_file, "w") as rh:
                rh.write(wig)    

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_read_libs(self):
        libs, texs = self.lib_reader.read_libs(self.libs, self.test_folder)
        self.assertDictEqual(texs,
        {'TSB_OD_0.2_TEX_reverse_TSB_OD_0.2_reverse': 0, 'TSB_OD_0.2_TEX_forward_TSB_OD_0.2_forward': 0})
        for index in range(0, 6):
            self.assertDictEqual(libs[index], self.lib_dict[index])

    def test_read_wig(self):
        merge_wig = os.path.join(self.test_folder, "merge.wig")
        for wig in os.listdir(self.test_folder):
            if wig.endswith(".wig"):
                os.system("cat " + os.path.join(self.test_folder, wig) + ">>" + merge_wig)
        wigs = self.lib_reader.read_wig(merge_wig, "+", self.lib_dict)
        for strain in wigs.keys():
            self.assertEqual(strain, "aaa")
        self.assertEqual(set(wigs["aaa"].keys()), set(['1_frag', '1_texnotex'])) 

class Example(object):
    
    libs = ["TSB_OD_0.2_TEX_reverse.wig:tex:1:a:-",
            "TSB_OD_0.2_TEX_forward.wig:tex:1:a:+",
            "TSB_OD_0.2_reverse.wig:notex:1:a:-",
            "TSB_OD_0.2_forward.wig:notex:1:a:+",
            "frag_forward.wig:frag:1:a:+",
            "frag_reverse.wig:frag:1:a:-"]
    
    wigs = {
"TSB_OD_0.2_TEX_forward.wig": """track type=wiggle_0 name="TSB_OD_0.2_TEX_forward"
variableStep chrom=aaa span=1
312 1.4041251228308191
313 56.867067474648174
314 56.867067474648174
315 56.867067474648174""",
"TSB_OD_0.2_TEX_reverse.wig": """track type=wiggle_0 name="TSB_OD_0.2_TEX_reverse"
variableStep chrom=aaa span=1
312 1.4041251228308191
313 56.867067474648174
314 56.867067474648174
315 56.867067474648174""",
"frag_forward.wig": """track type=wiggle_0 name="frag_forward"
variableStep chrom=aaa span=1
312 1.4041251228308191
313 56.867067474648174
314 56.867067474648174
315 56.867067474648174""",
"frag_reverse.wig": """track type=wiggle_0 name="frag_reverse"
variableStep chrom=aaa span=1
312 1.4041251228308191
313 56.867067474648174
314 56.867067474648174
315 56.867067474648174""",
"TSB_OD_0.2_forward.wig": """track type=wiggle_0 name="TSB_OD_0.2_forward"
variableStep chrom=aaa span=1
312 1.4041251228308191
313 56.867067474648174
314 56.867067474648174
315 56.867067474648174""",
"TSB_OD_0.2_reverse.wig": """track type=wiggle_0 name="TSB_OD_0.2_reverse"
variableStep chrom=aaa span=1
312 1.4041251228308191
313 56.867067474648174
314 56.867067474648174
315 56.867067474648174"""}

    lib_dict = [{'cond': '1_texnotex', 'strand': '-', 'type': 'tex', 'name': 'TSB_OD_0.2_TEX_reverse', 'rep': 'a'},
                {'cond': '1_texnotex', 'strand': '+', 'type': 'tex', 'name': 'TSB_OD_0.2_TEX_forward', 'rep': 'a'},
                {'cond': '1_texnotex', 'strand': '-', 'type': 'notex', 'name': 'TSB_OD_0.2_reverse', 'rep': 'a'},
                {'cond': '1_texnotex', 'strand': '+', 'type': 'notex', 'name': 'TSB_OD_0.2_forward', 'rep': 'a'},
                {'cond': '1_frag', 'strand': '+', 'type': 'frag', 'name': 'frag_forward', 'rep': 'a'},
                {'cond': '1_frag', 'strand': '-', 'type': 'frag', 'name': 'frag_reverse', 'rep': 'a'}]

    wig_file = [{"strain": "aaa", "strand": "+", "track": "TSB_OD_0.2_TEX_forward",
                "pos": "312", "coverage": "1.4041251228308191"},
                {"strain": "aaa", "strand": "+", "track": "TSB_OD_0.2_TEX_forward",
                "pos": "313", "coverage": "1.4041251228308191"},
                {"strain": "aaa", "strand": "+", "track": "TSB_OD_0.2_TEX_forward",
                "pos": "313", "coverage": "1.4041251228308191"},
                {"strain": "aaa", "strand": "-", "track": "TSB_OD_0.5_reverse",
                "pos": "312", "coverage": "1.4041251228308191"},
                {"strain": "aaa", "strand": "-", "track": "TSB_OD_0.5_reverse",
                "pos": "313", "coverage": "1.4041251228308191"},
                {"strain": "aaa", "strand": "-", "track": "TSB_OD_0.5_reverse",
                "pos": "313", "coverage": "1.4041251228308191"}]

    libs_for_wig = [{"name": "TSB_OD_0.2_TEX_forward", "type": "tex",
                     "cond": "1_texnotex", "rep": "a", "strand": "+"},
                    {"name": "TSB_OD_0.5_reverse", "type": "notex",
                     "cond": "2_texnotex", "rep": "a", "strand": "-"}]

if __name__ == "__main__":
    unittest.main()

