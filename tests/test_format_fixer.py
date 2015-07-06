#!/usr/bin/python

import os
import sys
import csv
import shutil
from io import StringIO
import unittest
sys.path.append(".")
from mock_helper import import_data
from annogesiclib.format_fixer import FormatFixer


class TestFormatFixer(unittest.TestCase):

    def setUp(self):
        self.fixer = FormatFixer()
        self.example = Example()
        self.ratt_out = self.example.ratt_out
        self.rnaplex_out = self.example.rnaplex_out
        self.emboss_out = self.example.emboss_out
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
        self.ratt_file = os.path.join(self.test_folder, "ratt.gff")
        with open(self.ratt_file, "w") as rh:
            rh.write(self.example.ratt_gff)
        self.rnaplex_file = os.path.join(self.test_folder, "rnaplex.txt")
        with open(self.rnaplex_file, "w") as rh:
            rh.write(self.example.rnaplex_file)
        self.emboss_file = os.path.join(self.test_folder, "emboss.txt")
        with open(self.emboss_file, "w") as rh:
            rh.write(self.example.emboss_file)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_fix_ratt(self):
        out = os.path.join(self.test_folder, "ratt.out")
        self.fixer.fix_ratt(self.ratt_file, "Staphylococcus_aureus_HG003", out)
        datas = import_data(out)
        self.assertEqual(set(datas), set(self.ratt_out.split("\n")))

    def test_fix_rnaplex(self):
        out_file = os.path.join(self.test_folder, "rnaplex.out")
        self.fixer.fix_rnaplex(self.rnaplex_file, out_file)
        datas = import_data(out_file)
        self.assertEqual(set(datas), set(self.rnaplex_out.split("\n")))

    def test_fix_emboss(self):
        out_file = os.path.join(self.test_folder, "emboss.out")
        self.fixer.fix_emboss(self.emboss_file, out_file)
        datas = import_data(out_file)
        self.assertEqual(set(datas), set(self.emboss_out.split("\n")))
        
class Example(object):

    ratt_gff = """##gff-version 3
chromosome.Staphylococcus_aureus_HG003.final	Refseq	source	1	2821337	.	+	.	mol_type=genomic DNA;db_xref=taxon:93061;strain=NCTC 8325;organism=Staphylococcus aureus subsp. aureus NCTC 8325;sub_species=aureus
chromosome.Staphylococcus_aureus_HG003.final	Refseq	gene	517	1878	.	+	.	gene=dnaA;db_xref=GeneID:3919798;locus_tag=SAOUHSC_00001
chromosome.Staphylococcus_aureus_HG003.final	Refseq	CDS	517	1878	.	+	.	gene=dnaA;db_xref=GI:88193824;db_xref=GeneID:3919798;transl_table=11;product=chromosomal replication initiation protein;note=binds to the dnaA-box as an ATP-bound complex at the origin of replication during the initiation of chromosomal replication, can also affect transcription of multiple genes including itself.;locus_tag=SAOUHSC_00001;protein_id=REF_uohsc:SAOUHSC00001;protein_id=YP_498609.1;codon_start=1
chromosome.Staphylococcus_aureus_HG003.final	Refseq	gene	2156	3289	.	+	.	db_xref=GeneID:3919799;locus_tag=SAOUHSC_00002
chromosome.Staphylococcus_aureus_HG003.final	Refseq	tRNA	2156	3289	.	+	.	EC_number=2.7.7.7;db_xref=GI:88193825;db_xref=GeneID:3919799;transl_table=11;product=DNA polymerase III subunit beta;note=binds the polymerase to DNA and acts as a sliding clamp;locus_tag=SAOUHSC_00002;protein_id=REF_uohsc:SAOUHSC00002;protein_id=YP_498610.1;codon_start=1"""

    ratt_out = """##gff-version 3
Staphylococcus_aureus_HG003	Refseq	source	1	2821337	.	+	.	mol_type=genomic DNA;db_xref=taxon:93061;strain=NCTC 8325;organism=Staphylococcus aureus subsp. aureus NCTC 8325;sub_species=aureus
Staphylococcus_aureus_HG003	Refseq	gene	517	1878	.	+	.	ID=gene0;Name=dnaA;gene=dnaA;db_xref=GeneID:3919798;locus_tag=SAOUHSC_00001
Staphylococcus_aureus_HG003	Refseq	CDS	517	1878	.	+	.	ID=cds0;Name=YP_498609.1;Parent=gene0;gene=dnaA;db_xref=GI:88193824;db_xref=GeneID:3919798;transl_table=11;product=chromosomal replication initiation protein;note=binds to the dnaA-box as an ATP-bound complex at the origin of replication during the initiation of chromosomal replication, can also affect transcription of multiple genes including itself.;locus_tag=SAOUHSC_00001;protein_id=REF_uohsc:SAOUHSC00001;protein_id=YP_498609.1;codon_start=1
Staphylococcus_aureus_HG003	Refseq	gene	2156	3289	.	+	.	ID=gene1;Name=SAOUHSC_00002;db_xref=GeneID:3919799;locus_tag=SAOUHSC_00002
Staphylococcus_aureus_HG003	Refseq	tRNA	2156	3289	.	+	.	ID=rna0;Name=SAOUHSC_00002;EC_number=2.7.7.7;db_xref=GI:88193825;db_xref=GeneID:3919799;transl_table=11;product=DNA polymerase III subunit beta;note=binds the polymerase to DNA and acts as a sliding clamp;locus_tag=SAOUHSC_00002;protein_id=REF_uohsc:SAOUHSC00002;protein_id=YP_498610.1;codon_start=1"""

    rnaplex_file = """>SAOUHSC_00001|dnaA
>srna1023
((((((&)))))) 571,576 :  20,25  (-5.30 = -7.89 +  0.18 +  2.41)
>SAOUHSC_00001|dnaA
>srna352
((((((((&)))))))) 163,170 :  24,31  (-1.91 = -8.31 +  0.60 +  5.80)
>SAOUHSC_00001|dnaA
>srna559
(((((((((((((&)))))))))).))) 301,313 :   4,17  (-5.43 = -9.60 +  3.14 +  1.03)
Error during initialization of the duplex in duplexfold_XS
>SAOUHSC_00002
>srna1023
((((((&)))))) 571,576 :  20,25  (-5.30 = -7.89 +  0.18 +  2.41)"""

    rnaplex_out = """>SAOUHSC_00001|dnaA
>srna1023
((((((&)))))) 571,576 :  20,25  (-5.30 = -7.89 +  0.18 +  2.41)
>SAOUHSC_00001|dnaA
>srna352
((((((((&)))))))) 163,170 :  24,31  (-1.91 = -8.31 +  0.60 +  5.80)
>SAOUHSC_00001|dnaA
>srna559
(((((((((((((&)))))))))).))) 301,313 :   4,17  (-5.43 = -9.60 +  3.14 +  1.03)
>SAOUHSC_00002
>srna1023
((((((&)))))) 571,576 :  20,25  (-5.30 = -7.89 +  0.18 +  2.41)"""

    emboss_file = """>A_1
DKSSNSFYKDLFIDFYIKILCITNKQDKVIHRLL
>B_1
NGIVPCLLSSPSILA*SALKRMSSLSLLVLLFAKAKX
>C_1
IELNHLSKQQKFGPTPYLSVVLFEESLLQYX"""

    emboss_out = """>A
DKSSNSFYKDLFIDFYIKILCITNKQDKVIHRLL
>B
NGIVPCLLSSPSILA*SALKRMSSLSLLVLLFAKAKX
>C
IELNHLSKQQKFGPTPYLSVVLFEESLLQYX"""
if __name__ == "__main__":
    unittest.main()
