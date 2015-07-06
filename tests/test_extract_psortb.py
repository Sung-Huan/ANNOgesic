import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import extract_info
import annogesiclib.extract_psortb as ep

class TestExtractPsortb(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_get_information(self):
        psortb_table = os.path.join(self.test_folder, "psortb.txt")
        with open(psortb_table, "w") as fh:
            fh.write(self.example.psortb)
        out_p = StringIO()
        ep.get_information(psortb_table, out_p, 0.5)
        self.assertEqual(set(out_p.getvalue().split("\n")[:-1]), set(self.example.out_file.split("\n")))
        out_p.close()

    def test_get_results_unique(self):
        line = "Cytoplasmic            9.97"
        seq_name = "Staphylococcus_aureus_HG003__YP_498609.1_+_517_1878"
        scores = [{'local': 'Cytoplasmic', 'score': 9.97},
                  {'local': 'Extracellular', 'score': 0.02},
                  {'local': 'Cellwall', 'score': 0.01},
                  {'local': 'CytoplasmicMembrane', 'score': 0.0}]
        psortbs = []
        out_p = StringIO()
        local_name, local_score = ep.get_results(line, scores, psortbs, out_p, seq_name, 0.5)
        self.assertListEqual(local_name, [])
        self.assertListEqual(local_score, [])
        self.assertEqual(set(out_p.getvalue().split("\n")[:-1]), set(["Staphylococcus_aureus_HG003\tYP_498609.1\t+\t517\t1878\tCytoplasmic\t9.97"]))
        out_p.close()       
 
    def test_get_results_multi(self):
        line = "Cytoplasmic (This protein may have multiple localization sites.) 10.00"
        seq_name = "Staphylococcus_aureus_HG003__YP_500332.1_-_1732031_1733725"
        scores = [{'local': 'Cytoplasmic', 'score': 10.0},
                  {'local': 'CytoplasmicMembrane', 'score': 10.0},
                  {'local': 'Extracellular', 'score': 0.0},
                  {'local': 'Cellwall', 'score': 0.0}]
        psortbs = []
        out_p = StringIO()
        local_name, local_score = ep.get_results(line, scores, psortbs, out_p, seq_name, 0.5)
        self.assertEqual(set(local_name), set(['Cytoplasmic', 'CytoplasmicMembrane']))
        self.assertEqual(set(local_score), set(['10.0', '10.0']))
        self.assertEqual(set(out_p.getvalue().split("\n")[:-1]),
                         set(["Staphylococcus_aureus_HG003\tYP_500332.1\t-\t1732031\t1733725\tCytoplasmic/CytoplasmicMembrane\t10.0/10.0"]))
        out_p.close()

    def test_import_psortb(self):
        seq_name = "Staphylococcus_aureus_HG003__YP_500332.1_-_1732031_1733725"
        psortbs = []
        results = ['Cytoplasmic', '(This', 'protein', 'may', 'have', 'multiple', 'localization', 'sites.)', '10.00']
        local_name = ['Cytoplasmic', 'CytoplasmicMembrane']
        local_score = ['10.0', '10.0']
        ep.import_psortb(seq_name, psortbs, local_name, local_score, "multi", results)
        self.assertDictEqual(psortbs[0], {'score': '10.0/10.0', 'local': 'Cytoplasmic/CytoplasmicMembrane',
                                          'seq_id': 'Staphylococcus_aureus_HG003', 'protein_id': 'YP_500332.1',
                                          'end': 1733725, 'strand': '-', 'start': 1732031})

    def test_print_gff(self):
        out_m = StringIO()
        psortbs = [{'score': '10.0/10.0', 'local': 'Cytoplasmic/CytoplasmicMembrane',
                    'seq_id': 'aaa', 'protein_id': 'YP_500332.1',
                    'end': 140, 'strand': '+', 'start': 100}]
        ep.print_gff(self.example.gffs, psortbs, out_m)
        datas, attributes = extract_info(out_m.getvalue(), "string")
        self.assertEqual(set(datas), set(self.example.out_data.split("\n")))
        ref_attributes = self.example.out_attributes.split("\n")
        for index in range(len(attributes)):
            self.assertEqual(set(attributes[index]), set(ref_attributes[index].split(";")))
        out_m.close()

class Example(object):
    psortb = """SeqID: Staphylococcus_aureus_HG003__YP_498609.1_+_517_1878
  Analysis Report:
    CMSVM+            Unknown                       [No details]
    CWSVM+            Unknown                       [No details]
    CytoSVM+          Cytoplasmic                   [No details]
    ECSVM+            Unknown                       [No details]
    ModHMM+           Unknown                       [No internal helices found]
    Motif+            Unknown                       [No motifs found]
    Profile+          Unknown                       [Module skipped]
    SCL-BLAST+        Cytoplasmic                   [matched 118704: Chromosomal replication initiator protein dnaA]
    SCL-BLASTe+       Unknown                       [No matches against database]
    Signal+           Unknown                       [No signal peptide detected]
  Localization Scores:
    Cytoplasmic            9.97
    Extracellular          0.02
    Cellwall               0.01
    CytoplasmicMembrane    0.00
  Final Prediction:
    Cytoplasmic            9.97

-------------------------------------------------------------------------------

SeqID: Staphylococcus_aureus_HG003__YP_498615.1_-_9755_10456
  Analysis Report:
    CMSVM+            Unknown                       [No details]
    CWSVM+            Unknown                       [No details]
    CytoSVM+          Unknown                       [No details]
    ECSVM+            Unknown                       [No details]
    ModHMM+           Unknown                       [No internal helices found]
    Motif+            Unknown                       [No motifs found]
    Profile+          Unknown                       [Module skipped]
    SCL-BLAST+        Unknown                       [No matches against database]
    SCL-BLASTe+       Unknown                       [No matches against database]
    Signal+           Unknown                       [No signal peptide detected]
  Localization Scores:
    CytoplasmicMembrane    2.50
    Extracellular          2.50
    Cytoplasmic            2.50
    Cellwall               2.50
  Final Prediction:
    Unknown

-------------------------------------------------------------------------------

SeqID: Staphylococcus_aureus_HG003__YP_500332.1_-_1732031_1733725
  Analysis Report:
    CMSVM+            CytoplasmicMembrane           [No details]
    CWSVM+            Unknown                       [No details]
    CytoSVM+          Cytoplasmic                   [No details]
    ECSVM+            Unknown                       [No details]
    ModHMM+           Unknown                       [1 internal helix found]
    Motif+            Unknown                       [No motifs found]
    Profile+          Unknown                       [Module skipped]
    SCL-BLAST+        Cytoplasmic, CytoplasmicMembrane[matched 38503056: Septation ring formation regulator ezrA]
    SCL-BLASTe+       Cytoplasmic, CytoplasmicMembrane[matched 100% 38503056: Septation ring formation regulator ezrA]
    Signal+           Unknown                       [No signal peptide detected]
  Localization Scores:
    Cytoplasmic            10.00
    CytoplasmicMembrane    10.00
    Extracellular          0.00
    Cellwall               0.00
  Final Prediction:
    Cytoplasmic (This protein may have multiple localization sites.) 10.00

-------------------------------------------------------------------------------

"""
    out_file = """Staphylococcus_aureus_HG003	YP_498609.1	+	517	1878	Cytoplasmic	9.97
Staphylococcus_aureus_HG003	YP_498615.1	-	9755	10456	Unknown	Unknown
Staphylococcus_aureus_HG003	YP_500332.1	-	1732031	1733725	Cytoplasmic/CytoplasmicMembrane	10.0/10.0"""

    gff_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "gene", "start": 100,
                 "end": 140, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 100,
                 "end": 140, "phase": ".", "strand": "+", "score": "."},
                {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 5400,
                 "end": 5800, "phase": ".", "strand": "-", "score": "."}]
    attributes_gff = [{"ID": "gene0", "Name": "gene_0", "locus_tag": "AAA_00001"},
                      {"ID": "cds0", "Name": "CDS_0", "locus_tag": "AAA_00001", "protein_id": "YP_500332.1"},
                      {"ID": "cds1", "Name": "CDS_1"}]
    gffs = []
    for index in range(0, 3):
        gffs.append(Create_generator(gff_dict[index], attributes_gff[index], "gff"))
    out_data = """aaa	Refseq	gene	100	140	.	+	.
aaa	Refseq	CDS	100	140	.	+	.
aaa	Refseq	CDS	5400	5800	.	-	."""

    out_attributes = """ID=gene0;Name=gene_0;locus_tag=AAA_00001
ID=cds0;Name=CDS_0;locus_tag=AAA_00001;protein_id=YP_500332.1;subcellular_localization=Cytoplasmic/CytoplasmicMembrane
ID=cds1;Name=CDS_1"""

if __name__ == "__main__":
    unittest.main()

