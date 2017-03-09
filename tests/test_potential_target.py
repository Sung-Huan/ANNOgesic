import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data, gen_file
from mock_args_container import MockClass
import annogesiclib.potential_target as pt


class TestPotentialTarget(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.mock_args = MockClass()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_read_file(self):
        seq_file = os.path.join(self.test_folder, "seq")
        gff_file = os.path.join(self.test_folder, "gff")
        gen_file(seq_file, self.example.seq_file)
        gen_file(gff_file, self.example.gff_file)
        fasta, cdss_f, cdss_r, genes = pt.read_file(
            seq_file, gff_file, "test", ["CDS"])
        self.assertEqual(
            fasta,
            "AGGATAGTCCGATACGTATACTGATAAAGACCGAAAATATTAGCGCGTAGC")
        self.assertEqual(cdss_f[0].start, 1)
        self.assertEqual(cdss_f[0].feature, "CDS")
        self.assertEqual(cdss_r[0].start, 14)
        self.assertEqual(cdss_r[0].feature, "CDS")
        self.assertEqual(len(genes), 2)
        self.assertEqual(genes[0].start, 1)
        self.assertEqual(genes[1].start, 14)

    def test_deal_cds_forward(self):
        pt.deal_cds_forward(self.example.cdss_f, self.test_folder,
                            self.example.fasta, self.example.genes, 2, 10)
        data = import_data(os.path.join(self.test_folder, "aaa_target.fa"))
        self.assertTrue("\n".join(data), self.example.cdsf_result)

    def test_deal_cds_reverse(self):
        pt.deal_cds_reverse(self.example.cdss_r, self.test_folder,
                            self.example.fasta, self.example.genes, 2, 10)
        data = import_data(os.path.join(self.test_folder, "aaa_target.fa"))
        self.assertTrue("\n".join(data), self.example.cdsf_result)

    def test_potential_target(self):
        seq_file = os.path.join(self.test_folder, "seq")
        gff_file = os.path.join(self.test_folder, "gff")
        gen_file(seq_file, self.example.seq_file)
        gen_file(gff_file, self.example.gff_file)
        args = self.mock_args.mock()
        args.tar_start = 2
        args.tar_end = 10
        args.features = ["CDS"]
        pt.potential_target(gff_file, seq_file, self.test_folder, args)
        data = import_data(os.path.join(self.test_folder, "aaa_target.fa"))
        self.assertTrue("\n".join(data), self.example.all_result)

class Example(object):

    seq_file = """>aaa
AGGATAGTCCGATACGTATACTGATAAAGACCGAAAATATTAGCGCGTAGC"""
    gff_file = """aaa\tRefseq\tgene\t1\t12\t.\t+\t.\tID=gene_0;Name=GENE_0;locus_tag=AAA_00001
aaa\tRefseq\tCDS\t1\t12\t.\t+\t.\tID=cds_0;Name=CDS_0;locus_tag=AAA_00001;protein_id="YP.00001
aaa\tRefseq\tgene\t14\t34\t.\t-\t.\tID=gene_1;Name=gene_1;locus_tag=AAA_00002
aaa\tRefseq\tCDS\t14\t34\t.\t-\t.\tID=cds_1;Name=CDS_1;locus_tag=AAA_00002;protein_id="YP.00002"""
    gene_dict = [
        {"seq_id": "aaa", "source": "Refseq", "feature": "gene", "start": 1,
         "end": 10, "phase": ".", "strand": "+", "score": "."},
        {"seq_id": "aaa", "source": "Refseq", "feature": "gene", "start": 12,
         "end": 23, "phase": ".", "strand": "+", "score": "."},
        {"seq_id": "aaa", "source": "Refseq", "feature": "gene", "start": 25,
         "end": 30, "phase": ".", "strand": "-", "score": "."},
        {"seq_id": "aaa", "source": "Refseq", "feature": "gene", "start": 33,
         "end": 43, "phase": ".", "strand": "-", "score": "."}]
    cdsf_dict = [
        {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 1,
         "end": 10, "phase": ".", "strand": "+", "score": "."},
        {"seq_id": "aaa", "source": "Refseq", "feature": "rRNA", "start": 12,
         "end": 23, "phase": ".", "strand": "+", "score": "."}]
    cdsr_dict = [
        {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 25,
         "end": 30, "phase": ".", "strand": "+", "score": "."},
        {"seq_id": "aaa", "source": "Refseq", "feature": "rRNA", "start": 33,
         "end": 43, "phase": ".", "strand": "+", "score": "."}]
    attributes_gene = [
        {"ID": "gene0", "Name": "danA", "locus_tag": "AAA_00001"},
        {"ID": "gene1", "Name": "AAA_00002", "locus_tag": "AAA_00002"},
        {"ID": "gene2", "Name": "AAA_00003", "locus_tag": "AAA_00003"},
        {"ID": "gene3", "Name": "hrcA", "locus_tag": "AAA_00004"}]
    attributes_cdsf = [
        {"ID": "cds0", "Name": "CDS_0", "locus_tag": "AAA_00001",
         "protein_id": "YP_000001", "Parent": "gene0"},
        {"ID": "cds1", "Name": "CDS_1", "locus_tag": "AAA_00002"}]
    attributes_cdsr = [
        {"ID": "cds2", "Name": "CDS_2", "locus_tag": "AAA_00003",
         "protein_id": "YP_000004", "Parent": "gene2"},
        {"ID": "cds3", "Name": "CDS_3", "locus_tag": "AAA_00004"}]
    fasta = "AGGATAGTCCGATACGTATACTGATAAAGACCGAAAATATTAGCGCGTAGC"
    genes = []
    cdss_f = []
    cdss_r = []
    for index in range(0, 2):
        cdss_f.append(Create_generator(
            cdsf_dict[index], attributes_cdsf[index], "gff"))
        cdss_r.append(Create_generator(
            cdsr_dict[index], attributes_cdsr[index], "gff"))
    for index in range(0, 4):
        genes.append(Create_generator(
            gene_dict[index], attributes_gene[index], "gff"))

    cdsf_result = """>AAA_00001|danA
AGGATAGTCCGATACGTATACTGATAAAGACCGAAAATATTAGCGCGTAGC
>AAA_00002|CDS_1
AGGATAGTCCGATACGTATACTGATAAAGACCGAAAATATTAGCGCGTAGC"""
    cdsr_result = """>AAA_00003
AGGATAGTCCGATACGTATACTGATAAAGACCGAAAATATTAGCGCGTAGC
>AAA_00004|CDS_3
AGGATAGTCCGATACGTATACTGATAAAGACCGAAAATATTAGCGCGTAGC"""
    all_result = """>AAA_00001|CDS_0
AGGATAGTCCGATACGTATACTGATAAAGACCGAAAATATTAGCGCGTAGC
>AAA_00002|CDS_1
GCTACGCGCTAATATTTTCGGTCTTTATCAGTATACGTATCGGACTATCCT"""

if __name__ == "__main__":
    unittest.main()

