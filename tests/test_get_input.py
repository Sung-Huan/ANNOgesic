#!/usr/bin/python

import unittest
import os
import sys
import shutil
sys.path.append(".")
from io import StringIO
import annogesiclib.get_input as get_input

class Mock_func(object):

    def mock_wget(self, ftp, input_folder, file_type):
        pass

    def modify_header(self):
        pass


class TestGetFile(unittest.TestCase):


    def setUp(self):
        self.test_folder = "test_project"
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)
        os.mkdir(self.test_folder)
        self.example = ExampleData()
        Seq_Editor = Mock_func
        get_input.wget = Mock_func().mock_wget

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_fasta(self):
        f_h = open(os.path.join(self.test_folder, "test.fasta"), "w")
        f_h.write(self.example.fasta_file)
        f_h.close()
        get_input.get_file(None, self.test_folder, "fna")
        self.assertTrue("test.fa" in os.listdir(self.test_folder))

    def test_gff(self):
        g_h = open(os.path.join(self.test_folder, "test.gff"), "w")
        g_h.write(self.example.gff_file)
        g_h.close()
        get_input.get_file(None, self.test_folder, "gff")
        self.assertTrue("ddd.gff" in os.listdir(self.test_folder))

    def test_gbk(self):
        g_h = open(os.path.join(self.test_folder, "test.gbk"), "w")
        g_h.write(self.example.gbk_file)
        g_h.close()
        get_input.get_file(None, self.test_folder, "gbk")
        self.assertTrue("NC_007795.1.gbk" in os.listdir(self.test_folder))

class ExampleData(object):

    gff_file = """##gff-version 3
#!gff-spec-version 1.20
#!processor NCBI annotwriter
##sequence-region ddd 1 2329769
##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=591001
ddd	Refseq	source	1	2821337	.	+	.	organism=Staphylococcus aureus subsp. aureus NCTC 8325;sub_species=aureus;db_xref=taxon:93061;mol_type=genomic DNA;strain=NCTC 8325
ddd	Refseq	gene	517	1878	.	+	.	locus_tag=SAOUHSC_00001;Name=dnaA;db_xref=GeneID:3919798;ID=gene0;gene=dnaA
ddd	Refseq	CDS	517	1878	.	+	.	gene=dnaA;transl_table=11;note=binds to the dnaA-box as an ATP-bound complex at the origin of replication during the initiation of chromosomal replication, can also affect transcription of multiple genes including itself.;ID=cds0;locus_tag=SAOUHSC_00001;db_xref=GeneID:3919798;Parent=gene0;codon_start=1;Name=YP_498609.1;start_TSS=TSS_313+;product=chromosomal replication initiation protein;Parent_tran=tran0;protein_id=YP_498609.1"""

    gbk_file = """LOCUS       NC_007795            2821361 bp    DNA     circular CON 16-MAY-2014
DEFINITION  Staphylococcus aureus subsp. aureus NCTC 8325 chromosome, complete
            genome.
ACCESSION   NC_007795
VERSION     NC_007795.1  GI:88193823
DBLINK      BioProject: PRJNA57795
KEYWORDS    RefSeq.
SOURCE      Staphylococcus aureus subsp. aureus NCTC 8325
  ORGANISM  Staphylococcus aureus subsp. aureus NCTC 8325
            Bacteria; Firmicutes; Bacilli; Bacillales; Staphylococcus.
REFERENCE   1  (bases 1 to 2821361)
  AUTHORS   Gillaspy,A.F., Worrell,V., Orvis,J., Roe,B.A., Dyer,D.W. and
            Iandolo,J.J.
  TITLE     The Staphylococcus aureus NCTC8325 Genome
  JOURNAL   (in) Fischetti,V., Novick,R., Ferretti,J., Portnoy,D. and Rood,J.
            (Eds.);
            GRAM POSITIVE PATHOGENS;
            ASM Press (2006)
REFERENCE   2  (bases 1 to 2821361)
  CONSRTM   NCBI Genome Project
  TITLE     Direct Submission
  JOURNAL   Submitted (18-FEB-2006) National Center for Biotechnology
            Information, NIH, Bethesda, MD 20894, USA
REFERENCE   3  (bases 1 to 2821361)
  AUTHORS   Gillaspy,A.F., Worrell,V., Orvis,J., Roe,B.A., Dyer,D.W. and
            Iandolo,J.J.
  TITLE     Direct Submission
  JOURNAL   Submitted (27-JAN-2006) Microbiology and Immunology, The University
            of Oklahoma Health Sciences Center, 940 Stanton L. Young Boulevard,
            Oklahoma City, OK 73104, USA
COMMENT     REVIEWED REFSEQ: This record has been curated by NCBI staff. The
            reference sequence was derived from CP000253.
            RefSeq Category: Reference Genome
                        UPR: UniProt Genome
            Staphylococcus aureus subsp. aureus NCTC 8325 is available from
            www.narsa.net.
            COMPLETENESS: full length."""

    fasta_file = """>gi|1111111111|ref|aa| test fasta
CGCAGGTTGAGTTCCTGTTCCCGATAGATCCGATAAACCCGCTTATGATTCCAGAGCTGTCCCTGCACAT
TGCGCAGATACAGGAAACACAGACCAAATCCCCATCTCCTGTGAGCCTGGGTCAGTCCCACCAGAAGAGC
GGCAATCCTGTCGTTCTCCGCTGCCAGTCGCGGACGATAGCGAAAGCAGGTCTCGGATATCCCAAAAATC
CGACAGGCCAGCGCAATGCTGACCCCATGATGCGCCACAGCTTGTGCGGCCAGTTCCCGGCGCTGGGCTG
GCCGCTTCATTTTTTTCCAAGGGCTTCCTTCAGGATATCCGTCTGCATGCTCAAATCCGCATACATGCGC"""

if __name__ == "__main__":
    unittest.main()
