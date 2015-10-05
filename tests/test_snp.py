import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file, import_data
import annogesiclib.snp as sn
from annogesiclib.snp import SNPCalling

class Mock_func(object):

    def __init__(self):
        self.example = Example()

    def mock_run_tools(self, samtools_path, bcftools_path,
                       fasta_file, out_bcf, out_raw_prefix, type_):
        pass

    def mock_transcript_snp(self, fasta, snp, out_table_prefix, qual,
                        seq_path, fasta_prefix, read_depth, stat_path,
                        types, bam_number, fraction, table_path):
        pass

    def mock_run_sub(self, samtools_path, fasta_file, type_, out_raw_prefix,
                 bcftools_path, seq_path, prefix, out_table_prefix, quality,
                 depth, stat_path, bam_number, fraction, table_path):
        gen_file("test_folder/test", "test")

    def mock_run_bam(self, samtools_path, sub_command, bam_file):
        pass

    def mock_get_header(self, samtools_path):
        pass

class TestSNPCalling(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        self.fasta = os.path.join(self.test_folder, "fasta")
        self.snp_folder = os.path.join(self.test_folder, "snp")
        self.table = os.path.join(self.test_folder, "table")
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
            os.mkdir(self.fasta)
            os.mkdir(self.snp_folder)
            os.mkdir(self.table)
        self.snp = SNPCalling("reference", self.test_folder, self.fasta)
        self.mock = Mock_func()

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_import_bam(self):
        gen_file(os.path.join(self.test_folder, "test_1.bam"), "test")
        gen_file(os.path.join(self.test_folder, "test_2.bam"), "test")
        bams = []
        num_bams = self.snp._import_bam(self.test_folder, bams)
        self.assertEqual(num_bams, 2)
        self.assertListEqual(bams, ['test_folder/test_1.bam', 'test_folder/test_2.bam'])

    def test_transcript_snp(self):
        fasta = os.path.join(self.test_folder, "NC_007795.1.fa")
        gen_file(fasta, self.example.fasta)
        snp = os.path.join(self.test_folder, "NC_007795.1.csv")
        gen_file(snp, self.example.snp)
        self.snp._transcript_snp(fasta, snp, "test", 10,
                        self.test_folder, "test", 5, self.snp_folder,
                        "reference", 2, 0.3, self.table)
        datas = import_data(os.path.join(self.snp_folder, "stat_test_reference_SNP.csv"))
        self.assertEqual("\n".join(datas), self.example.out_stat)
        datas = import_data(os.path.join(self.test_folder, "test_NC_007795.1_1_1.fa"))
        self.assertEqual("\n".join(datas), ">NC_007795.1\nAaTTGaaTCCCGAACGACAGTTAT")
        os.remove("test_seq_reference.csv")
        os.remove("test_depth_only.vcf")
        os.remove("test_depth_quality.vcf")
        os.remove("test_NC_007795.1_SNP_QUAL.png")    

    def test_run_sub(self):
        self.snp._run_tools = self.mock.mock_run_tools
        self.snp._transcript_snp = self.mock.mock_transcript_snp
        os.mkdir(os.path.join(self.test_folder, "with_BAQ"))
        self.snp._run_sub("test", "fasta", "with", "no",
                          "test", self.test_folder, "test", "no", 10,
                          3, "stat", 2, 0.5, "table")
        self.assertTrue(os.path.exists(os.path.join(self.test_folder, "with_BAQ/test")))

    def test_run_program(self):
        self.snp._run_sub = self.mock.mock_run_sub
        self.snp._run_program(["5"], "test", "test", "fasta", "test", "test", 10, "seq",
                              "test", 10, "stat", 2, 0.5, "out", "table")
        self.assertFalse(os.path.exists(os.path.join(self.test_folder, "test")))
        self.snp._run_program(["1"], "test", "test", "fasta", "test", "test", 10, "seq",
                              "test", 10, "stat", 2, 0.5, "out", "table")
        self.assertTrue(os.path.exists(os.path.join(self.test_folder, "test")))

    def test_detect_fasta(self):
        datas = self.snp._detect_fasta("test.fa")
        self.assertEqual(datas, (True, 'test'))

    def test_merge_bams(self):
        self.snp._run_bam = self.mock.mock_run_bam
        normal_bams = os.path.join(self.test_folder, "tex_bams")
        frag_bams = os.path.join(self.test_folder, "frag_bams")
        os.mkdir(normal_bams)
        os.mkdir(frag_bams)
        gen_file(os.path.join(normal_bams, "tex.bam"), "test")
        gen_file(os.path.join(normal_bams, "notex.bam"), "test")
        gen_file(os.path.join(frag_bams, "farg.bam"), "test")
        num = self.snp._merge_bams(normal_bams, frag_bams, "test", self.test_folder)
        self.assertEqual(num, 3)

    def test_modify_header(self):
        gen_file(os.path.join(self.fasta, "test.fa"), ">AAA|BBB|CCC|DDD|EEE\nAATTAATTGGCC")
        self.snp._modify_header(self.fasta)
        datas = import_data(os.path.join(self.fasta, "test.fa"))
        self.assertEqual("\n".join(datas), ">DDD\nAATTAATTGGCC")

    def test_get_genome_name(self):
        self.snp._get_header = self.mock.mock_get_header
        gen_file(os.path.join(self.test_folder, "header"), self.example.bam)
        seq_names = self.snp._get_genome_name("test")

    def test_run_snp_calling(self):
        self.snp._get_header = self.mock.mock_get_header
        self.snp._run_bam = self.mock.mock_run_bam
        self.snp._run_sub = self.mock.mock_run_sub
        self.snp._run_tools = self.mock.mock_run_tools
        self.snp._transcript_snp = self.mock.mock_transcript_snp
        normal_bams = os.path.join(self.test_folder, "tex_bams")
        frag_bams = os.path.join(self.test_folder, "frag_bams")
        os.mkdir(normal_bams)
        os.mkdir(frag_bams)
        gen_file(os.path.join(normal_bams, "tex.bam"), "test")
        gen_file(os.path.join(normal_bams, "notex.bam"), "test")
        gen_file(os.path.join(frag_bams, "farg.bam"), "test")
        gen_file(os.path.join(self.fasta, "test.fa"), ">AAA|BBB|CCC|DDD|EEE\nAATTAATTGGCC")
        gen_file(os.path.join(self.test_folder, "header"), self.example.bam)
        gen_file(os.path.join(self.test_folder, "whole_reads.bam"), "test")
        gen_file(os.path.join(self.test_folder, "whole_reads_sorted.bam"), "test")
        self.snp.run_snp_calling("test", "test", "reference", ["1"],
                        self.fasta, normal_bams, frag_bams, 10, 5,
                        self.test_folder, 0.5)


class Example(object):
    snp = """NC_007795.1	2	.	C	A	98	.	DP=89;VDB=8.46526e-15;SGB=-0.693147;RPB=0.00816879;MQB=1;MQSB=1;BQB=0.760636;MQ0F=0;AC=2;AN=2;DP4=0,4,79,4;MQ=20	GT:PL:DP        1/1:125,184,0:87
NC_007795.1	6	.	A	AA	26.9515	.	INDEL;IDV=22;IMF=0.536585;DP=41;VDB=9.36323e-14;SGB=-0.692562;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=18,0,22,0;MQ=20	GT:PL:DP	0/1:60,0,55:40"""
    fasta = """>NC_007795.1
ACTTGATCCCGAACGACAGTTAT"""
    out_stat = """NC_007795.1:
Read depth should be higher than 5:
The fraction of Maximum read depth of insertion or deletion should be higher than 0.3:
the number of QUAL which is between 0 and 10 = 0
the number of QUAL which is between 10 and 20 = 0
the number of QUAL which is between 20 and 30 = 1
the number of QUAL which is between 30 and 40 = 0
the number of QUAL which is between 40 and 50 = 0
the number of QUAL which is between 50 and 60 = 0
the number of QUAL which is between 60 and 70 = 0
the number of QUAL which is between 70 and 80 = 0
the number of QUAL which is between 80 and 90 = 0
the number of QUAL which is between 90 and 100 = 1
the number of QUAL which is between 100 and 110 = 0
the total numbers of QUAL which is higher than 10 = 2"""
    bam = """@HD     VN:1.0
@SQ     SN:NC_007795.1  LN:2821361
@PG     ID:segemehl     VN:0.1.4-$Rev: 380 $ ($Date: 2012-12-17 12:43:51 +0100 (Mon, 17 Dec 2012) $)    CL:segemehl --query RAPL_analysis/output/read_alignments-processed_reads/pMEM_OD_0.2.fa_processed.fa --index RAPL_analysis/output/read_alignments-index/index.idx --database RAPL_analysis/input/reference_sequences/NC_007795.fa --outfile RAPL_analysis/output/read_alignments-alignments/pMEM_OD_0.2.fa_alignments.sam --hitstrategy 1 --accuracy 95 --evalue 5 --threads 24 --nomatchfilename RAPL_analysis/output/read_alignments-unaligned_reads/pMEM_OD_0.2.fa_unaligned.fa"""

if __name__ == "__main__":
    unittest.main()
