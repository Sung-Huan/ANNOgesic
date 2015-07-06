import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data, extract_info, gen_file
import annogesiclib.extract_sRNA_info as esi

class Mock_func(object):

    def __init__(self):
        self.example = Example()

    def mock_read_gff(self, srna_file, data_type):
        srnas = []
        for index in range(0, 2):
            srnas.append(Create_generator(self.example.srna_dict[index],
                         self.example.attributes_srna[index], "gff"))
        return srnas

class TestExtractsRNAInfo(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_detect_nr(self):
        blasts = {"hit_num": 0, "blast": False, "name": ""}
        nr_file = os.path.join(self.test_folder, "test.nr")
        out_t = StringIO()
        with open(nr_file, "w") as fh:
            fh.write(self.example.blast_nr_protein)
        with open(nr_file) as blast_f:
            for line in blast_f:
                esi.detect_nr(line, blast_f.readlines(), out_t, blasts, "test")
        self.assertEqual(out_t.getvalue(),
                         "test	DNA replication and repair protein RecF domain protein	EHS30036.1,EHS80331.1,EID88948.1	4e-18\n")
        self.assertDictEqual(blasts, {'name': '', 'hit_num': 1, 'blast': True})
        out_t.close()
        blast_f.close()
        out_t = StringIO()
        blasts = {"hit_num": 0, "blast": False, "name": ""}
        with open(nr_file, "w") as fh:
            fh.write(self.example.blast_nr_unknown)
        with open(nr_file) as blast_f:
            for line in blast_f.readlines():
                esi.detect_nr(line, blast_f.readlines(), out_t, blasts, "test")
        self.assertEqual(out_t.getvalue(), "")
        self.assertDictEqual(blasts, {'name': '', 'hit_num': 0, 'blast': False})
        out_t.close()
        blast_f.close()

    def test_detect_srna(self):
        blasts = {"hit_num": 0, "blast": False, "name": ""}
        nr_file = os.path.join(self.test_folder, "test.srna")
        out_t = StringIO()
        with open(nr_file, "w") as fh:
            fh.write(self.example.blast_srna)
        with open(nr_file) as blast_f:
            for line in blast_f:
                esi.detect_srna(line, blast_f, out_t, blasts, "test")
        self.assertEqual(out_t.getvalue(), "test	ssau217.1|Staphylococcus aureus subsp. aureus N315|RsaK	4e-89\n")
        out_t.close()
        blast_f.close()

    def test_extract_blast(self):
        esi.read_gff = Mock_func().mock_read_gff
        nr_blast = os.path.join(self.test_folder, "nr_table")
        gen_file(nr_blast, self.example.blast_nr_all)
        srna_blast = os.path.join(self.test_folder, "srna_table")
        gen_file(srna_blast, self.example.blast_srna_all)
        output_file = os.path.join(self.test_folder, "out.gff")
        output_table = os.path.join(self.test_folder, "out.csv")
        esi.extract_blast(nr_blast, "test.srna", output_file, output_table, "nr")
        datas, attributes = extract_info(output_file, "file")
        refs, ref_attributes = extract_info(self.example.out_nr_gff, "string")
        self.assertEqual(set(datas), set(refs))
        self.assertEqual(set(attributes[0]), set(attributes[0]))
        self.assertEqual(set(attributes[1]), set(attributes[1]))
        datas = import_data(output_table)
        self.assertEqual(set(datas), set(self.example.out_nr_csv.split("\n")))
        esi.extract_blast(srna_blast, "test.srna", output_file, output_table, "sRNA")
        datas, attributes = extract_info(output_file, "file")
        refs, ref_attributes = extract_info(self.example.out_srna_gff, "string")
        self.assertEqual(set(datas), set(refs))
        self.assertEqual(set(attributes[0]), set(attributes[0]))
        self.assertEqual(set(attributes[1]), set(attributes[1]))
        datas = import_data(output_table)
        self.assertEqual(set(datas), set(self.example.out_srna_csv.split("\n")))

class Example(object):

    srna_dict = [{"start": 313, "end": 417, "phase": ".",
                  "strand": "+", "seq_id": "Staphylococcus_aureus_HG003", "score": ".",
                  "source": "Refseq", "feature": "sRNA"},
                 {"start": 4045, "end": 4159, "phase": ".",
                  "strand": "-", "seq_id": "Staphylococcus_aureus_HG003", "score": ".",
                  "source": "Refseq", "feature": "sRNA"}]
    attributes_srna = [{"ID": "srna0", "Name": "sRNA_candidate_0000"},
                       {"ID": "srna1", "Name": "sRNA_candidate_0001"}]
    blast_nr_protein = """> gi|375036980|gb|EHS30036.1| DNA replication and repair protein
RecF domain protein [Staphylococcus aureus subsp. aureus IS-111]
gi|375376817|gb|EHS80331.1| DNA replication and repair
protein RecF domain protein [Staphylococcus aureus subsp.
aureus IS-157] gi|383972925|gb|EID88948.1| DNA replication and
repair protein RecF domain protein [Staphylococcus aureus
subsp. aureus CO-23]
Length=95

 Score = 80.5 bits (197),  Expect = 4e-18, Method: Compositional matrix adjust.
 Identities = 37/37 (100%), Positives = 37/37 (100%), Gaps = 0/37 (0%)
 Frame = -3
"""

    blast_nr_unknown = """> gi|515566342|ref|WP_016999177.1| hypothetical protein [Staphylococcus
lentus]
Length=371

 Score = 80.5 bits (197),  Expect = 4e-18, Method: Compositional matrix adjust.
 Identities = 37/37 (100%), Positives = 37/37 (100%), Gaps = 0/37 (0%)
 Frame = -3
"""

    blast_srna = """> ssau217.1|Staphylococcus aureus subsp. aureus N315|RsaK
Length=209

 Score =   318 bits (172),  Expect = 4e-89
 Identities = 172/172 (100%), Gaps = 0/172 (0%)
 Strand=Plus/Minus

"""

    blast_srna_all = """BLASTN 2.2.29+


Reference: Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb
Miller (2000), "A greedy algorithm for aligning DNA sequences", J
Comput Biol 2000; 7(1-2):203-14.



Database: ANNOgesic/input/database/sRNA_database.fa
           897 sequences; 182,899 total letters



Query= srna0|Staphylococcus_aureus_HG003|313|417|+

Length=105


***** No hits found *****



Lambda      K        H
    1.33    0.621     1.12

Gapped
Lambda      K        H
    1.28    0.460    0.850

Effective search space used: 15000683


Query= srna1|Staphylococcus_aureus_HG003|4045|4159|-

Length=242
                                                                      Score     E
Sequences producing significant alignments:                          (Bits)  Value

  ssau217.1|Staphylococcus aureus subsp. aureus N315|RsaK               318   4e-89


> ssau217.1|Staphylococcus aureus subsp. aureus N315|RsaK
Length=209

 Score =   318 bits (172),  Expect = 4e-89
 Identities = 172/172 (100%), Gaps = 0/172 (0%)
 Strand=Plus/Minus

Query  1    GCAAACAAATAAGATTCAATAAGATGTGATGGTAAGCTAGATAAACTAACGTATCAGTTT  60
            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  172  GCAAACAAATAAGATTCAATAAGATGTGATGGTAAGCTAGATAAACTAACGTATCAGTTT  113

Query  61   TATATACTACTTAGTCAATTGCTTCTTTTACATTAATTACATACACCAACGTGTTACTAA  120
            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  112  TATATACTACTTAGTCAATTGCTTCTTTTACATTAATTACATACACCAACGTGTTACTAA  53

Query  121  GTAAGATTAGGCATGAGTTAAATTAAGCTGTGATGGTTACCAACACAGTCTA  172
            ||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  52   GTAAGATTAGGCATGAGTTAAATTAAGCTGTGATGGTTACCAACACAGTCTA  1



Lambda      K        H
    1.33    0.621     1.12

Gapped
Lambda      K        H
    1.28    0.460    0.850

Effective search space used: 37721250


"""
    blast_nr_all = """BLASTX 2.2.29+


Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro A.
Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J.
Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of
protein database search programs", Nucleic Acids Res. 25:3389-3402.



Database: ANNOgesic/input/database/nr.fa
           53,846,081 sequences; 19,408,592,203 total letters



Query= srna0|Staphylococcus_aureus_HG003|313|417|+

Length=105


***** No hits found *****



Lambda      K        H        a         alpha
   0.318    0.134    0.401    0.792     4.96

Gapped
Lambda      K        H        a         alpha    sigma
   0.267   0.0410    0.140     1.90     42.6     43.6

Effective search space used: 492023414324


Query= srna1|Staphylococcus_aureus_HG003|4045|4159|-

Length=115
                                                                      Score     E
Sequences producing significant alignments:                          (Bits)  Value

  gi|375036980|gb|EHS30036.1| DNA replication and repair protein ...  80.5    4e-18
  gi|334276553|gb|EGL94807.1| putative DNA replication and repair...  80.9    1e-17

> gi|375036980|gb|EHS30036.1| DNA replication and repair protein
RecF domain protein [Staphylococcus aureus subsp. aureus IS-111]
gi|375376817|gb|EHS80331.1| DNA replication and repair
protein RecF domain protein [Staphylococcus aureus subsp.
aureus IS-157] gi|383972925|gb|EID88948.1| DNA replication and
repair protein RecF domain protein [Staphylococcus aureus
subsp. aureus CO-23]
Length=95

 Score = 80.5 bits (197),  Expect = 4e-18, Method: Compositional matrix adjust.
 Identities = 37/37 (100%), Positives = 37/37 (100%), Gaps = 0/37 (0%)
 Frame = -3

Query  113  LALAKSHRTSNDKELIRFNADYAKIEGELSYRHGTMP  3
            LALAKSHRTSNDKELIRFNADYAKIEGELSYRHGTMP
Sbjct  46   LALAKSHRTSNDKELIRFNADYAKIEGELSYRHGTMP  82



> gi|334276553|gb|EGL94807.1| putative DNA replication and repair
protein RecF [Staphylococcus aureus subsp. aureus 21318]
gi|645287686|gb|KDP53072.1| AAA domain protein [Staphylococcus
aureus subsp. aureus CO-86]
Length=174

 Score = 80.9 bits (198),  Expect = 1e-17, Method: Compositional matrix adjust.
 Identities = 37/37 (100%), Positives = 37/37 (100%), Gaps = 0/37 (0%)
 Frame = -3

Query  113  LALAKSHRTSNDKELIRFNADYAKIEGELSYRHGTMP  3
            LALAKSHRTSNDKELIRFNADYAKIEGELSYRHGTMP
Sbjct  46   LALAKSHRTSNDKELIRFNADYAKIEGELSYRHGTMP  82



Lambda      K        H        a         alpha
   0.318    0.134    0.401    0.792     4.96

Gapped
Lambda      K        H        a         alpha    sigma
   0.267   0.0410    0.140     1.90     42.6     43.6

Effective search space used: 487823420006


"""
    out_nr_gff = """##gff-version 3
Staphylococcus_aureus_HG003	Refseq	sRNA	313	417	.	+	.	ID=srna0;Name=sRNA_candidate_0000;nr_hit=NA
Staphylococcus_aureus_HG003	Refseq	sRNA	4045	4159	.	-	.	ID=srna1;Name=sRNA_candidate_0001;nr_hit=2"""
    out_nr_csv = """Staphylococcus_aureus_HG003	srna0	+	313	417	NA
Staphylococcus_aureus_HG003	srna1	-	4045	4159	DNA replication and repair protein RecF domain protein	EHS30036.1,EHS80331.1,EID88948.1	4e-18
Staphylococcus_aureus_HG003	srna1	-	4045	4159	putative DNA replication and repair protein RecF	EGL94807.1	1e-17
Staphylococcus_aureus_HG003	srna1	-	4045	4159	AAA domain protein	KDP53072.1	1e-17"""
    out_srna_gff = """##gff-version 3
Staphylococcus_aureus_HG003	Refseq	sRNA	313	417	.	+	.	Name=sRNA_candidate_0000;ID=srna0;sRNA_hit=NA
Staphylococcus_aureus_HG003	Refseq	sRNA	4045	4159	.	-	.	Name=sRNA_candidate_0001;ID=srna1;sRNA_hit=1"""
    out_srna_csv = """Staphylococcus_aureus_HG003	srna0	+	313	417	NA
Staphylococcus_aureus_HG003	srna1	-	4045	4159	ssau217.1|Staphylococcus aureus subsp. aureus N315|RsaK	4e-89"""

if __name__ == "__main__":
    unittest.main()

