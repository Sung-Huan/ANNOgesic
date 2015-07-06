import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data, gen_file
import annogesiclib.gene_ontology as go


class Mock_gff_parser(object):

    def __init__(self):
        self.example = Example()

    def entries(self, fh):
        fh.close()
        for element in self.example.gffs:
            yield element

class Mock_func(object):

    def __init__(self):
        self.example = Example()

    def mock_plot(self, total_nums, strain, filename, total, out_folder):
        pass

class TestGeneOntology(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_retrieve_uniprot(self):
        go.Gff3Parser = Mock_gff_parser
        database_file = os.path.join(self.test_folder, "database")
        gen_file(database_file, self.example.idmapping)
        gff_file = os.path.join(self.test_folder, "test.gff")
        gen_file(gff_file, "test")
        out_file = os.path.join(self.test_folder, "out.gff")
        go.retrieve_uniprot(database_file, gff_file, out_file)
        datas = import_data(out_file)
        self.assertEqual(set(datas), set(self.example.out_retrieve.split("\n")))

    def test_import_obo(self):
        obo_file = os.path.join(self.test_folder, "obo.txt")
        gen_file(obo_file, self.example.obo)
        obos = go.import_obo(obo_file)
        self.assertListEqual(obos, self.example.out_obo)

    def test_compare_go_slim(self):
        classes = {'All_strain': {'cellular_component': {}, 'biological_process': {}, 'molecular_function': {}},
                   'aaa': {'cellular_component': {}, 'biological_process': {}, 'molecular_function': {}}}
        total_nums = {'All_strain': {'cellular_component': 0, 'total': 0, 'biological_process': 0, 'molecular_function': 0},
                      'aaa': {'cellular_component': 0, 'total': 0, 'biological_process': 0, 'molecular_function': 0}}
        gos = {'aaa': ['GO:0000003']}
        go.compare_go_slim(gos, self.example.out_obo, self.example.slim_obos, classes, total_nums)
        self.assertDictEqual(classes, {'All_strain': {'molecular_function': {}, 'cellular_component': {},
                                                      'biological_process': {'reproduction': 1}},
                                       'aaa': {'molecular_function': {}, 'cellular_component': {},
                                               'biological_process': {'reproduction': 1}}})
        self.assertDictEqual(total_nums, {'All_strain': {'molecular_function': 0, 'cellular_component': 0,
                                                         'biological_process': 1, 'total': 1},
                                          'aaa': {'molecular_function': 0, 'cellular_component': 0,
                                                  'biological_process': 1, 'total': 1}})

    def test_map2goslim(self):
        go.plot = Mock_func().mock_plot
        stat_file = os.path.join(self.test_folder, "stat.txt")
        term_file = os.path.join(self.test_folder, "term.txt")
        slim_file = os.path.join(self.test_folder, "slim.txt")
        go_table = os.path.join(self.test_folder, "go.csv")
        gen_file(term_file, self.example.obo)
        gen_file(slim_file, self.example.slim)
        gen_file(go_table, "aaa\t+\t150\t200\tYP_031579.1\tGO:0000003")
        go.map2goslim(slim_file, term_file, go_table, stat_file, self.test_folder)
        datas = import_data(stat_file)
        self.assertEqual(set(datas), set(self.example.out_stat.split("\n")))

class Example(object):
    idmapping = """Q6GZX4	001R_FRG3G	2947773	YP_031579.1	81941549; 47060116; 49237298		GO:0006355; GO:0046782; GO:0006351	UniRef100_Q6GZX4	UniRef90_Q6GZX4 UniRef50_Q6GZX4	UPI00003B0FD4		654924			15165820	AY548484	AAT09660.1
Q6GZX3	002L_FRG3G	2947774	YP_031580.1	49237299; 47060117; 81941548		GO:0033644; GO:0016021	UniRef100_Q6GZX3	UniRef90_Q6GZX3 UniRef50_Q6GZX3 UPI00003B0FD5		654924			15165820	AY548484	AAT09661.1
Q197F8	002R_IIV3	4156251	YP_654574.1	106073503; 109287880; 123808694			UniRef100_Q197F8	UniRef90_Q197F8 UniRef50_Q197F8	UPI0000D83464		345201			16912294	DQ643392	ABF82032.1"""
    gff_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 150,
                 "end": 200, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 1230,
                 "end": 1240, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "bbb", "source": "Refseq", "feature": "CDS", "start": 7100,
                 "end": 9167, "phase": ".", "strand": "-", "score": "."},
                {"seq_id": "bbb", "source": "Refseq", "feature": "CDS", "start": 100,
                 "end": 167, "phase": ".", "strand": "-", "score": "."}]
    attributes_gff = [{"ID": "cds0", "Name": "YP_031579.1", "locus_tag": "AAA_00001", "protein_id": "YP_031579.1"},
                      {"ID": "cds1", "Name": "CDS_1", "locus_tag": "AAA_00002", "protein_id": "YP_031580.1"},
                      {"ID": "cds2", "Name": "YP_654574.1", "locus_tag": "BBB_00001"},
                      {"ID": "cds3", "Name": "YP_031579.1", "locus_tag": "BBB_00002"}]
    gffs = []
    for index in range(0, 4):
        gffs.append(Create_generator(gff_dict[index], attributes_gff[index], "gff"))
    out_retrieve = """strain	strand	start	end	protein_id	Go_term
aaa	+	150	200	YP_031579.1	GO:0006355; GO:0046782; GO:0006351
aaa	+	1230	1240	YP_031580.1	GO:0033644; GO:0016021
bbb	-	100	167	YP_031579.1	GO:0006355; GO:0046782; GO:0006351
bbb	-	7100	9167	YP_654574.1"""

    obo = """format-version: 1.2
data-version: releases/2014-12-12
date: 11:12:2014 15:49
saved-by: tb
auto-generated-by: TermGenie 1.0
subsetdef: Cross_product_review "Involved_in"
subsetdef: goantislim_grouping "Grouping classes that can be excluded"
subsetdef: gocheck_do_not_annotate "Term not to be used for direct annotation"
subsetdef: gocheck_do_not_manually_annotate "Term not to be used for direct manual annotation"
subsetdef: goslim_aspergillus "Aspergillus GO slim"
subsetdef: goslim_candida "Candida GO slim"
subsetdef: goslim_generic "Generic GO slim"
subsetdef: goslim_goa "GOA and proteome slim"
subsetdef: goslim_metagenomics "Metagenomics GO slim"
subsetdef: goslim_pir "PIR GO slim"
subsetdef: goslim_plant "Plant GO slim"
subsetdef: goslim_pombe "Fission yeast GO slim"
subsetdef: goslim_synapse "synapse GO slim"
subsetdef: goslim_virus "Viral GO slim"
subsetdef: goslim_yeast "Yeast GO slim"
subsetdef: gosubset_prok "Prokaryotic GO subset"
subsetdef: mf_needs_review "Catalytic activity terms in need of attention"
subsetdef: termgenie_unvetted "Terms created by TermGenie that do not follow a template and require additional vetting by editors"
subsetdef: virus_checked "Viral overhaul terms"
synonymtypedef: systematic_synonym "Systematic synonym" EXACT
default-namespace: gene_ontology
remark: cvs version: $Revision: 22213 $
ontology: go

[Term]
id: GO:0000001
name: mitochondrion inheritance
namespace: biological_process
def: "The distribution of mitochondria, including the mitochondrial genome, into daughter cells after mitosis or meiosis, mediated by interactions between mitochondria and the cytoskeleton." [GOC:mcc, PMID:10873824, PMID:11389764]
synonym: "mitochondrial inheritance" EXACT []
is_a: GO:0048308 ! organelle inheritance
is_a: GO:0048311 ! mitochondrion distribution

[Term]
id: GO:0000002
name: mitochondrial genome maintenance
namespace: biological_process
def: "The maintenance of the structure and integrity of the mitochondrial genome; includes replication and segregation of the mitochondrial chromosome." [GOC:ai, GOC:vw]
is_obsolete: true
is_a: GO:0007005 ! mitochondrion organization

[Term]
id: GO:0000003
name: reproduction
namespace: biological_process
alt_id: GO:0019952
alt_id: GO:0050876
def: "The production of new individuals that contain some portion of genetic material inherited from one or more parent organisms." [GOC:go_curators, GOC:isa_complete, GOC:jl, ISBN:0198506732]
subset: goslim_generic
subset: goslim_pir
subset: goslim_plant
subset: gosubset_prok
synonym: "reproductive physiological process" EXACT []
xref: Wikipedia:Reproduction
is_a: GO:0008150 ! biological_process
disjoint_from: GO:0044848 ! biological phase

"""

    slim = """format-version: 1.2
subsetdef: Cross_product_review "Involved_in"
subsetdef: goantislim_grouping "Grouping classes that can be excluded"
subsetdef: gocheck_do_not_annotate "Term not to be used for direct annotation"
subsetdef: gocheck_do_not_manually_annotate "Term not to be used for direct manual annotation"
subsetdef: goslim_aspergillus "Aspergillus GO slim"
subsetdef: goslim_candida "Candida GO slim"
subsetdef: goslim_generic "Generic GO slim"
subsetdef: goslim_goa "GOA and proteome slim"
subsetdef: goslim_metagenomics "Metagenomics GO slim"
subsetdef: goslim_pir "PIR GO slim"
subsetdef: goslim_plant "Plant GO slim"
subsetdef: goslim_pombe "Fission yeast GO slim"
subsetdef: goslim_synapse "synapse GO slim"
subsetdef: goslim_virus "Viral GO slim"
subsetdef: goslim_yeast "Yeast GO slim"
subsetdef: gosubset_prok "Prokaryotic GO subset"
subsetdef: mf_needs_review "Catalytic activity terms in need of attention"
subsetdef: termgenie_unvetted "Terms created by TermGenie that do not follow a template and require additional vetting by editors"
subsetdef: virus_checked "Viral overhaul terms"
synonymtypedef: systematic_synonym "Systematic synonym" EXACT
ontology: go/subsets/goslim_generic

[Term]
id: GO:0000003
name: reproduction
namespace: biological_process
alt_id: GO:0019952
alt_id: GO:0050876
def: "The production of new individuals that contain some portion of genetic material inherited from one or more parent organisms." [GOC:go_curators, GOC:isa_complete, GOC:jl, ISBN:0198506732]
subset: goslim_generic
subset: goslim_pir
subset: goslim_plant
subset: gosubset_prok
synonym: "reproductive physiological process" EXACT []
xref: Wikipedia:Reproduction
is_a: GO:0008150 ! biological_process

[Term]
id: GO:0000228
name: nuclear chromosome
namespace: cellular_component
def: "A chromosome found in the nucleus of a eukaryotic cell." [GOC:mah]
subset: goslim_generic
synonym: "nuclear interphase chromosome" NARROW []
is_a: GO:0005694 ! chromosome
intersection_of: GO:0005694 ! chromosome
intersection_of: part_of GO:0005634 ! nucleus
relationship: part_of GO:0005634 ! nucleus

"""
    out_obo = [{'name': 'mitochondrion inheritance', 'id': 'GO:0000001', 'namespace': 'biological_process', 'synonym': '"mitochondrial inheritance" EXACT []', 'def': '"The distribution of mitochondria, including the mitochondrial genome, into daughter cells after mitosis or meiosis, mediated by interactions between mitochondria and the cytoskeleton." [GOC:mcc, PMID:10873824, PMID:11389764]', 'is_a': ['GO:0048308 ! organelle inheritance', 'GO:0048311 ! mitochondrion distribution']}, {'name': 'mitochondrial genome maintenance', 'id': 'GO:0000002', 'is_obsolete': 'true', 'namespace': 'biological_process', 'def': '"The maintenance of the structure and integrity of the mitochondrial genome; includes replication and segregation of the mitochondrial chromosome." [GOC:ai, GOC:vw]', 'is_a': ['GO:0007005 ! mitochondrion organization']}, {'name': 'reproduction', 'subset': 'gosubset_prok', 'disjoint_from': 'GO:0044848 ! biological phase', 'synonym': '"reproductive physiological process" EXACT []', 'id': 'GO:0000003', 'namespace': 'biological_process', 'alt_id': 'GO:0050876', 'xref': 'Wikipedia:Reproduction', 'def': '"The production of new individuals that contain some portion of genetic material inherited from one or more parent organisms." [GOC:go_curators, GOC:isa_complete, GOC:jl, ISBN:0198506732]', 'is_a': ['GO:0008150 ! biological_process']}]
    slim_obos = [{'id': 'GO:0000003', 'def': '"The production of new individuals that contain some portion of genetic material inherited from one or more parent organisms." [GOC:go_curators, GOC:isa_complete, GOC:jl, ISBN:0198506732]', 'namespace': 'biological_process', 'alt_id': 'GO:0050876', 'subset': 'gosubset_prok', 'name': 'reproduction', 'is_a': ['GO:0008150 ! biological_process'], 'xref': 'Wikipedia:Reproduction', 'synonym': '"reproductive physiological process" EXACT []'}, {'id': 'GO:0000228', 'def': '"A chromosome found in the nucleus of a eukaryotic cell." [GOC:mah]', 'namespace': 'cellular_component', 'subset': 'goslim_generic', 'name': 'nuclear chromosome', 'is_a': ['GO:0005694 ! chromosome'], 'intersection_of': 'part_of GO:0005634 ! nucleus', 'relationship': 'part_of GO:0005634 ! nucleus', 'synonym': '"nuclear interphase chromosome" NARROW []'}]

    out_stat = """aaa:
	molecular_function: 0(percentage in total: 0.0)
	cellular_component: 0(percentage in total: 0.0)
	biological_process: 1(percentage in total: 1.0)
		reproduction: 1(percentage in biological_process: 1.0)"""

if __name__ == "__main__":
    unittest.main()

