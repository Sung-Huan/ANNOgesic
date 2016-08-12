import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file, import_data
from mock_args_container import MockClass
from annogesiclib.ratt import RATT


class Mock_func(object):

    def __init__(self):
        self.example = Example()

    def mock_run_ratt(self, ratt_path, tar, ref, out):
        gen_file("test_folder/gffs/tmp.gff", self.example.gff_file)
        gen_file("test_folder/gffs/tmp.ptt", self.example.ptt_file)
        gen_file("test_folder/gffs/tmp.rnt", self.example.rnt_file)
        pass


class TestRATT(unittest.TestCase):

    def setUp(self):
        self.mock_args = MockClass()
        self.test_folder = "test_folder"
        self.ref_embls = "test_folder/embls"
        self.output_path = "test_folder/output"
        self.tar_fastas = "test_folder/tar_fasta"
        self.ref_fastas = "test_folder/ref_fasta"
        self.gff_outfolder = "test_folder/gffs"
        self.ref_gbk = "test_folder/gbk"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
            os.mkdir(self.ref_embls)
            os.mkdir(self.ref_gbk)
            os.mkdir(self.output_path)
            os.mkdir(self.tar_fastas)
            os.mkdir(self.ref_fastas)
            os.mkdir(self.gff_outfolder)
        args = self.mock_args.mock()
        args.output_path = self.output_path
        args.ref_embls = self.ref_embls
        args.ref_gbk = self.ref_gbk
        args.tar_fastas = self.tar_fastas
        args.ref_fastas = self.ref_fastas
        args.gff_outfolder = self.gff_outfolder
        self.ratt = RATT(args)
        self.example = Example()

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_convert_to_pttrnt(self):
        files = ["aaa.gff"]
        gen_file(os.path.join(self.test_folder, "aaa.gff"), self.example.gff_file)
        os.mkdir(os.path.join(self.tar_fastas, "tmp"))
        gen_file(os.path.join(self.tar_fastas, "tmp/aaa.fa"), self.example.fasta_file)
        self.ratt._convert_to_pttrnt(self.test_folder, files)
        data = import_data(os.path.join(self.test_folder, "aaa.rnt"))
        self.assertEqual("\n".join(data), self.example.rnt_file)
        data = import_data(os.path.join(self.test_folder, "aaa.ptt"))
        self.assertEqual("\n".join(data), self.example.ptt_file)

    def test_convert_to_gff(self):
        files = ["aaa.gff"]
        ratt_result = "chromosome.aaa.final.embl"
        gen_file(os.path.join(self.output_path, ratt_result), self.example.embl_file)
        args = self.mock_args.mock()
        args.output_path = self.output_path
        args.gff_outfolder = self.gff_outfolder
        self.ratt._convert_to_gff(ratt_result, args, files)
#        self.ratt._convert_to_gff(ratt_result, self.output_path, self.gff_outfolder, files)
        data = import_data(os.path.join(self.output_path, "aaa.gff"))
        self.assertEqual("\n".join(data), self.example.embl_gff)
        data = import_data(os.path.join(self.gff_outfolder, "aaa.gff"))
        self.assertEqual("\n".join(data), self.example.embl_gff)

    def test_parser_embl_gbk(self):
        files = [os.path.join(self.test_folder, "aaa.gbk")]
        gen_file(os.path.join(self.test_folder, "aaa.gbk"), self.example.gbk_file)
        self.ratt._parser_embl_gbk(files)
        data = import_data(os.path.join(self.ref_gbk, "gbk_tmp/NC_007795.1.gbk"))
        self.assertEqual("\n".join(data), self.example.gbk_file.split("//")[0] + "//")
        data = import_data(os.path.join(self.ref_gbk, "gbk_tmp/NC_007799.1.gbk"))
        self.assertEqual("\n".join(data), self.example.gbk_file.split("//")[1].strip() + "\n//")

    def test_convert_embl(self):
        gen_file(os.path.join(self.test_folder, "aaa.gbk"), self.example.gbk_file.split("//")[0])
        out = self.ratt._convert_embl(self.test_folder)
        self.assertEqual(out, "test_folder/gbk/gbk_tmp")
        self.assertTrue(os.path.exists("test_folder/gbk/gbk_tmp"))

    def test_format_and_run(self):
        self.ratt._run_ratt = Mock_func().mock_run_ratt
        args = self.mock_args.mock()
        args.output_path
        args.pairs = ["NC_007795.1:Staphylococcus_aureus_HG003"]
        args.element = "chromosome"
        self.ratt._format_and_run(args)

    def test_annotation_transfer(self):    
        gen_file(os.path.join(self.ref_fastas, "aaa.fa"), self.example.fasta_file)
        gen_file(os.path.join(self.tar_fastas, "bbb.fa"), self.example.fasta_file)
        gen_file(os.path.join(self.ref_embls, "aaa.gbk"), self.example.gbk_file.split("//")[0])
        self.ratt._run_ratt = Mock_func().mock_run_ratt
        args = self.mock_args.mock()
        args.element = "element"
        args.ref_embls = self.ref_embls
        args.tar_fastas = self.tar_fastas
        args.ref_fastas = self.ref_fastas
        args.output_path = self.output_path
        args.gff_outfolder = self.gff_outfolder
        args.pairs = ["aaa:bbb"]
        args.convert = True
        self.ratt.annotation_transfer(args)
#        self.ratt.annotation_transfer("test", "element", "test_type",
#                                      self.ref_embls, self.tar_fastas,
#                                      self.ref_fastas, self.output_path,
#                                      True, self.gff_outfolder, pairs)
        self.assertTrue(os.path.exists(os.path.join(self.gff_outfolder, "bbb.gff")))
        self.assertTrue(os.path.exists(os.path.join(self.gff_outfolder, "bbb.rnt")))
        self.assertTrue(os.path.exists(os.path.join(self.gff_outfolder, "bbb.ptt")))


class Example(object):

    fasta_file = """>aaa
AAATATGCGATAGCGTTTCCCGCGGGATAGAAATAGTTCGCGTACCTATCATG"""
    gff_file = """aaa\tRefseq\tCDS\t1\t14\t.\t+\t.\tID=cds_0;Name=CDS_00000;locus_tag=AAA_00001;protein_id=YP.90090
aaa\tRefseq\ttRNA\t16\t24\t.\t+\t.\tID=trna_0;Name=tRNA_00000;locus_tag=AAA_00002
aaa\tRefseq\ttRNA\t29\t40\t.\t+\t.\tID=rrna_0;Name=rRNA_00000;locus_tag=AAA_00003"""
    rnt_file = """aaa - 1..53
2 proteins
Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
16..24	+	9	-	-	AAA_00002	-	-	-
29..40	+	12	-	-	AAA_00003	-	-	-"""
    ptt_file = """aaa - 1..53
1 proteins
Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
1..14	+	14	YP.90090	-	AAA_00001	-	-	-"""
    embl_file = """ID                   aaa ; ; ; ; ; 2821337 BP.
FH   Key             Location/Qualifiers
FH
FT   source          1..2821337
FT                   /strain="NCTC 8325"
FT                   /sub_species="aureus"
FT                   /mol_type="genomic DNA"
FT                   /organism="Staphylococcus aureus subsp. aureus NCTC 8325"
FT                   /db_xref="taxon:93061"
FT   gene            517..1878
FT                   /db_xref="GeneID:3919798"
FT                   /gene="dnaA"
FT                   /locus_tag="SAOUHSC_00001"
FT   CDS             517..1878
FT                   /product="chromosomal replication initiation protein"
FT                   /db_xref="GI:88193824"
FT                   /db_xref="GeneID:3919798"
FT                   /transl_table=11
FT                   /codon_start=1
FT                   /note="binds to the dnaA-box as an ATP-bound complex at the
FT                   origin of replication during the initiation of chromosomal
FT                   replication; can also affect transcription of multiple
FT                   genes including itself."
FT                   /locus_tag="SAOUHSC_00001"
FT                   /protein_id="REF_uohsc:SAOUHSC00001"
FT                   /protein_id="YP_498609.1"
FT                   /gene="dnaA"
FT   misc_feature    520..1872
FT                   /db_xref="CDD:234667"
FT                   /note="chromosomal replication initiation protein;
FT                   Reviewed; Region: dnaA; PRK00149"
FT                   /gene="dnaA"
FT                   /locus_tag="SAOUHSC_00001
FT   gene            3670..3915
FT                   /db_xref="GeneID:3919176"
FT                   /locus_tag="SAOUHSC_00003"
FT   CDS             3670..3915
FT                   /product="hypothetical protein"
FT                   /db_xref="GI:88193826"
FT                   /db_xref="GeneID:3919176"
FT                   /transl_table=11
FT                   /codon_start=1
FT                   /note="conserved hypothetical protein"
FT                   /locus_tag="SAOUHSC_00003"
FT                   /protein_id="REF_uohsc:SAOUHSC00003"
FT                   /protein_id="YP_498611.1"
FT   misc_feature    3688..3861
FT                   /db_xref="CDD:132033"
FT                   /note="S4 domain protein YaaA; Region: YaaA_near_RecF;
FT                   TIGR02988"
FT                   /locus_tag="SAOUHSC_00003"""""
    embl_gff = """##gff-version 3
aaa	Refseq	source	1	2821337	.	+	.	strain=NCTC 8325;sub_species=aureus;mol_type=genomic DNA;organism=Staphylococcus aureus subsp. aureus NCTC 8325;db_xref=taxon:93061
aaa	Refseq	gene	517	1878	.	+	.	ID=gene0;Name=dnaA;db_xref=GeneID:3919798;gene=dnaA;locus_tag=SAOUHSC_00001
aaa	Refseq	CDS	517	1878	.	+	.	ID=cds0;Name=YP_498609.1;Parent=gene0;product=chromosomal replication initiation protein;db_xref=GI:88193824;db_xref=GeneID:3919798;transl_table=11;codon_start=1;note=binds to the dnaA-box as an ATP-bound complex at the origin of replication during the initiation of chromosomal replication, can also affect transcription of multiple genes including itself.;locus_tag=SAOUHSC_00001;protein_id=REF_uohsc:SAOUHSC00001;protein_id=YP_498609.1;gene=dnaA
aaa	Refseq	gene	3670	3915	.	+	.	ID=gene1;Name=SAOUHSC_00003;db_xref=GeneID:3919176;locus_tag=SAOUHSC_00003
aaa	Refseq	CDS	3670	3915	.	+	.	ID=cds1;Name=YP_498611.1;Parent=gene1;product=hypothetical protein;db_xref=GI:88193826;db_xref=GeneID:3919176;transl_table=11;codon_start=1;note=conserved hypothetical protein;locus_tag=SAOUHSC_00003;protein_id=REF_uohsc:SAOUHSC00003;protein_id=YP_498611.1"""
    gbk_file = """LOCUS       NC_007795            2821361 bp    DNA     circular CON 16-MAY-2014
DEFINITION  Staphylococcus aureus subsp. aureus NCTC 8325 chromosome, complete
            genome.
ACCESSION   NC_007795
VERSION     NC_007795.1  GI:88193823
FEATURES             Location/Qualifiers
     source          1..2821361
                     /organism="Staphylococcus aureus subsp. aureus NCTC 8325"
                     /mol_type="genomic DNA"
                     /strain="NCTC 8325"
                     /sub_species="aureus"
                     /db_xref="taxon:93061"
     gene            517..1878
                     /gene="dnaA"
                     /locus_tag="SAOUHSC_00001"
                     /db_xref="GeneID:3919798"
     CDS             517..1878
                     /gene="dnaA"
                     /locus_tag="SAOUHSC_00001"
                     /note="binds to the dnaA-box as an ATP-bound complex at
                     the origin of replication during the initiation of
                     chromosomal replication; can also affect transcription of
                     multiple genes including itself."
                     /codon_start=1
                     /transl_table=11
                     /product="chromosomal replication initiation protein"
                     /protein_id="REF_uohsc:SAOUHSC00001"
                     /protein_id="YP_498609.1"
                     /db_xref="GI:88193824"
                     /db_xref="GeneID:3919798"
                     /translation="MSEKEIWEKVLEIAQEKLSAVSYSTFLKDTELYTIKDGEAIVLS
                     SIPFNANWLNQQYAEIIQAILFDVVGYEVKPHFITTEELANYSNNETATPKETTKPST
                     ETTEDNHVLGREQFNAHNTFDTFVIGPGNRFPHAASLAVAEAPAKAYNPLFIYGGVGL
                     GKTHLMHAIGHHVLDNNPDAKVIYTSSEKFTNEFIKSIRDNEGEAFRERYRNIDVLLI
                     DDIQFIQNKVQTQEEFFYTFNELHQNNKQIVISSDRPPKEIAQLEDRLRSRFEWGLIV
                     DITPPDYETRMAILQKKIEEEKLDIPPEALNYIANQIQSNIRELEGALTRLLAYSQLL
                     GKPITTELTAEALKDIIQAPKSKKITIQDIQKIVGQYYNVRIEDFSAKKRTKSIAYPR
                     QIAMYLSRELTDFSLPKIGEEFGGRDHTTVIHAHEKISKDLKEDPIFKQEVENLEKEI
                     RNV"
//
LOCUS       NC_007799            2821361 bp    DNA     circular CON 16-MAY-2014
DEFINITION  Staphylococcus aureus subsp. aureus NCTC 8325 chromosome, complete
            genome.
ACCESSION   NC_007795
VERSION     NC_007799.1  GI:88193823
     source          1..2821361
                     /organism="Staphylococcus aureus subsp. aureus NCTC 8325"
                     /mol_type="genomic DNA"
                     /strain="NCTC 8325"
                     /sub_species="aureus"
                     /db_xref="taxon:93061"
     gene            517..1878
                     /gene="dnaA"
                     /locus_tag="SAOUHSC_00001"
                     /db_xref="GeneID:3919798"
     CDS             517..1878
                     /gene="dnaA"
                     /locus_tag="SAOUHSC_00001"
                     /note="binds to the dnaA-box as an ATP-bound complex at
                     the origin of replication during the initiation of
                     chromosomal replication; can also affect transcription of
                     multiple genes including itself."
                     /codon_start=1
                     /transl_table=11
                     /product="chromosomal replication initiation protein"
                     /protein_id="REF_uohsc:SAOUHSC00001"
                     /protein_id="YP_498609.1"
                     /db_xref="GI:88193824"
                     /db_xref="GeneID:3919798"
                     /translation="MSEKEIWEKVLEIAQEKLSAVSYSTFLKDTELYTIKDGEAIVLS
                     SIPFNANWLNQQYAEIIQAILFDVVGYEVKPHFITTEELANYSNNETATPKETTKPST
                     ETTEDNHVLGREQFNAHNTFDTFVIGPGNRFPHAASLAVAEAPAKAYNPLFIYGGVGL
                     GKTHLMHAIGHHVLDNNPDAKVIYTSSEKFTNEFIKSIRDNEGEAFRERYRNIDVLLI
                     DDIQFIQNKVQTQEEFFYTFNELHQNNKQIVISSDRPPKEIAQLEDRLRSRFEWGLIV
                     DITPPDYETRMAILQKKIEEEKLDIPPEALNYIANQIQSNIRELEGALTRLLAYSQLL
                     GKPITTELTAEALKDIIQAPKSKKITIQDIQKIVGQYYNVRIEDFSAKKRTKSIAYPR
                     QIAMYLSRELTDFSLPKIGEEFGGRDHTTVIHAHEKISKDLKEDPIFKQEVENLEKEI
                     RNV"
     misc_feature    520..1872
                     /gene="dnaA"
                     /locus_tag="SAOUHSC_00001"
                     /note="chromosomal replication initiation protein;
                     Reviewed; Region: dnaA; PRK00149"
                     /db_xref="CDD:234667"
//"""

if __name__ == "__main__":
    unittest.main()

