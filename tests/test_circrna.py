import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file
from annogesiclib.circrna import CircRNADetection
from mock_args_container import MockClass


class Mock_segemehl(object):

    def mock_fasta_index(self, segemehl_path, fasta_path,
                         index, fasta):
        pass

    def mock_align(self, args,
                   index, fasta, read,
                   sam_file, log_file,
                   fasta_prefix):
        return "test"

    def mock_wait_processes(self, processes):
        pass

class Mock_samtools(object):

    def mock_covert_bam(self, samtools_path, out_bam, pre_sam):
        pass

class TestCircRNADetection(unittest.TestCase):

    def setUp(self):
        self.segemehl = Mock_segemehl()
        self.samtools = Mock_samtools()
        self.mock_args = MockClass()
        self.example = Example()
        self.test_folder = "test_folder"
        self.fasta_folder = os.path.join(self.test_folder, "fasta")
        self.gff_folder = os.path.join(self.test_folder, "gff")
        self.out_folder = os.path.join(self.test_folder, "output")
        self.read_folder = os.path.join(self.test_folder, "read")
        self.splice_folder = os.path.join(self.test_folder, "splice")
        self.alignment_path = os.path.join(self.out_folder,
                                  "segemehl_align")
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)
        if (not os.path.exists(self.fasta_folder)):
            os.mkdir(self.fasta_folder)
            os.mkdir(os.path.join(self.fasta_folder, "tmp"))
        if (not os.path.exists(self.gff_folder)):
            os.mkdir(self.gff_folder)
        if (not os.path.exists(self.out_folder)):
            os.mkdir(self.out_folder)
        if (not os.path.exists(self.read_folder)):
            os.mkdir(self.read_folder)
        if (not os.path.exists(self.splice_folder)):
            os.mkdir(self.splice_folder)
        args = self.mock_args.mock()
        args.output_folder = self.out_folder
        args.gffs = self.gff_folder
        args.align = True
        args.fastas = self.fasta_folder
        self.circ = CircRNADetection(args)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)
        if os.path.exists("test1"):
            if os.path.isfile("test1"):
                os.remove("test1")
            if os.path.isdir("test1"):
                shutil.rmtree("test1")
        if os.path.exists("test2"):
            if os.path.isfile("test2"):
                os.remove("test2")
            if os.path.isdir("test2"):
                shutil.rmtree("test2")

    def test_deal_zip_file(self):
        out1 = os.path.join(self.test_folder, "test1.fa")
        out2 = os.path.join(self.test_folder, "test2")
        gen_file(out1, self.example.fasta_file)
        gen_file(out2, self.example.fasta_file)
        os.system("gzip " + out1)
        os.system("bzip2 -z " + out2)
        reads = self.circ._deal_zip_file([out1 + ".gz", out2 + ".bz2"])
        self.assertEqual(set(reads), set([out1, out2 + ".fa"]))
        self.assertTrue(os.path.exists(out1))
        self.assertTrue(os.path.exists(out2 + ".fa"))

    def test_align(self):
        self.circ._run_segemehl_fasta_index = self.segemehl.mock_fasta_index
        self.circ._run_segemehl_align = self.segemehl.mock_align
        self.circ._wait_process = self.segemehl.mock_wait_processes
        fasta1 = os.path.join(os.path.join(self.fasta_folder, "tmp/test1.fa"))
        fasta2 = os.path.join(os.path.join(self.fasta_folder, "tmp/test2.fa"))
        read1 = os.path.join(self.read_folder, "read1.fa")
        read2 = os.path.join(self.read_folder, "read2.fa")
        gen_file(fasta1, self.example.fasta_file)
        gen_file(fasta2, self.example.fasta_file)
        gen_file(read1, self.example.fasta_file)
        gen_file(read2, self.example.fasta_file)
        os.mkdir(os.path.join(self.out_folder, "segemehl_align"))
        args = self.mock_args.mock()
        args.output_folder = self.out_folder
        args.gffs = self.gff_folder
        args.align = True
        args.fastas = self.fasta_folder
        args.segemehl_path = None
        args.read_files = [read1, read2]
        args.cores = 2
        read_datas = [{"sample": "test", "files": [read1, read2]}]
        align_results, prefixs = self.circ._align(args, read_datas)
        self.assertEqual(set(align_results), set(['read1_test1', 'read2_test1',
                                                  'read1_test2', 'read2_test2']))
        self.assertEqual(set(prefixs), set(['test1', 'test2']))

    def test_convert_sam2bam(self):
        self.circ._run_samtools_convert_bam = self.samtools.mock_covert_bam
        sam1 = os.path.join(self.test_folder, "test1.sam")
        sam2 = os.path.join(self.test_folder, "test2.sam")
        bam = os.path.join(self.test_folder, "test3.bam")
        gen_file(sam1, self.example.align_file)
        gen_file(sam2, self.example.align_file)
        gen_file(bam, self.example.align_file)
        align_files = ["test1"]
        bam_files, convert_ones, remove_ones = self.circ._convert_sam2bam(
            self.test_folder, None, align_files)
        self.assertEqual(set(bam_files), set([bam, sam1.replace("sam", "bam"),
                                              sam2.replace("sam", "bam")]))
        self.assertEqual(set(convert_ones), set([sam2.replace("sam", "bam")]))
        self.assertEqual(set(remove_ones), set([sam1]))
        align_files = ["test3"]
        bam_files, convert_ones, remove_ones = self.circ._convert_sam2bam(
            self.test_folder, None, align_files)
        self.assertEqual(set(convert_ones), set([sam2.replace("sam", "bam"),
                                                 sam1.replace("sam", "bam")]))
        self.assertEqual(set(remove_ones), set([]))

    def test_merge_bed(self):
        fasta1 = os.path.join(self.fasta_folder, "test1.fa")
        fasta2 = os.path.join(self.fasta_folder, "test2.fa")
        header1 = os.path.join(self.splice_folder,
                               "Staphylococcus_aureus_HG003")
        header2 = os.path.join(self.splice_folder, "aaa")
        header3 = os.path.join(self.splice_folder, "bbb")
        os.mkdir(header1)
        os.mkdir(header2)
        os.mkdir(header3)
        splice1 = os.path.join(header1, "Staphylococcus_aureus_HG003_a1_splicesites.bed")
        splice2 = os.path.join(header2, "aaa_a1_splicesites.bed")
        splice3 = os.path.join(header3, "bbb_a1_splicesites.bed")
        tran1 = os.path.join(header1, "Staphylococcus_aureus_HG003_a1_transrealigned.bed")
        tran2 = os.path.join(header2, "aaa_a1_transrealigned.bed")
        tran3 = os.path.join(header3, "bbb_a1_transrealigned.bed")
        gen_file(fasta1, self.example.fasta_file)
        gen_file(fasta2, self.example.multi_fasta_file)
        gen_file(splice1, self.example.splice_file)
        gen_file(splice2, self.example.splice_file)
        gen_file(splice3, self.example.splice_file)
        gen_file(tran1, self.example.tran_file)
        gen_file(tran2, self.example.tran_file)
        gen_file(tran3, self.example.tran_file)
        prefixs = self.circ._merge_bed(
            self.fasta_folder, self.splice_folder, self.out_folder)
        self.assertEqual(set(prefixs[1]), set(["test1", "test2"]))
        self.assertEqual(prefixs[0][0], "_a1_")
        self.assertTrue(os.path.exists(os.path.join(
            self.out_folder, "test1", "test1_a1_splicesites.bed")))
        self.assertTrue(os.path.exists(os.path.join(
            self.out_folder, "test1", "test1_a1_transrealigned.bed")))
        self.assertTrue(os.path.exists(os.path.join(
            self.out_folder, "test2", "test2_a1_splicesites.bed")))
        self.assertTrue(os.path.exists(os.path.join(
            self.out_folder, "test2", "test2_a1_transrealigned.bed")))
        self.assertTrue(os.path.exists(os.path.join(
            self.out_folder, "test2", "tmp_bbb_a1_splicesites.bed")))
        self.assertTrue(os.path.exists(os.path.join(
            self.out_folder, "test2", "tmp_aaa_a1_splicesites.bed")))
        self.assertTrue(os.path.exists(os.path.join(
            self.out_folder, "test2", "tmp_aaa_a1_transrealigned.bed")))
        self.assertTrue(os.path.exists(os.path.join(
            self.out_folder, "test2", "tmp_bbb_a1_transrealigned.bed")))

    def test_combine_read_bam(self):
        bam_datas = [{"sample": "aaa", "files": [
            os.path.join(self.out_folder, "segemehl_align", "aaa1.bam"),
            "aaa2.bam"]},
                     {"sample": "bbb", "files": ["bbb1.bam", "bbb2.bam"]}]
        read_datas = [{"sample": "aaa", "files": [
            "aaa1.fa", "aaa3.fa", "aaa4.fa"]}]
        bam_files = [os.path.join(self.out_folder,
                                  "segemehl_align", "aaa1.bam"),
                     os.path.join(self.out_folder,
                                  "segemehl_align", "aaa3.bam")]
        self.circ._combine_read_bam(bam_files, bam_datas, read_datas)
        self.assertDictEqual(bam_datas[0], {'files': [
            'test_folder/output/segemehl_align/aaa1.bam', 'aaa2.bam',
            'test_folder/output/segemehl_align/aaa3.bam'], 'sample': 'aaa'})

class Example(object):

    fasta_file = """>Staphylococcus_aureus_HG003
GCAGGTTGAGTTCCTGTTCCCGATAGATCCGATAAACCCGCTTATGATTCCAGAGCTGTCCCTGCACATA
GCCGCTTCATTTTTTTCCAAGGGCTTCCTTCAGGATATCCGTCTGCATGCTCAAATCCGCATACATGCGG""" 

    multi_fasta_file = """>aaa
GCAGGTTGAGTTCCTGTTCCCGATAGATCCGATAAACCCGCTTATGATTCCAGAGCTGTCCCTGCACATA
GCCGCTTCATTTTTTTCCAAGGGCTTCCTTCAGGATATCCGTCTGCATGCTCAAATCCGCATACATGCGG
>bbb
GCAGGTTGAGTTCCTGTTCCCGATAGATCCGATAAACCCGCTTATGATTCCAGAGCTGTCCCTGCACATA
GCCGCTTCATTTTTTTCCAAGGGCTTCCTTCAGGATATCCGTCTGCATGCTCAAATCCGCATACATGCGG"""

    align_file = """HWI-ST863:185:C1K2FACXX:7:1101:2763:1989 1:N:0:TGACCAAT 0       Staphylococcus_aureus_HG003     1821329 255     88M1D6M *       0       0       GTCAGATGACATTAAATAGCATCTCCTCGTGTTGATTATTTTGGTTGGCTGACCAATATTTATTCTAGCACGTAGAGATGCATTTTTTAAAAAA       *       NM:i:4  MD:Z:88^G0T0G1C2        NH:i:1  XI:i:0  XA:Z:Q
HWI-ST863:185:C1K2FACXX:7:1101:3037:1980 1:N:0:TGACCAAT 16      Staphylococcus_aureus_HG003     1923770 255     94M     *       0       0       CGTCAGCAACGTCATATGAATTCTCAGTTCATGTTGTGGTGACACTTTAAACGGTCTGTGCCAGTAGCGACCGAGTCATTTCAAGAATGACCAT       *       NM:i:0  MD:Z:94 NH:i:1  XI:i:0  XA:Z:Q
HWI-ST863:185:C1K2FACXX:7:1101:4006:1911 1:N:0:TGACCAAT 0       Staphylococcus_aureus_HG003     498013  255     94M     *       0       0       TCTGGTGACTATAGCAAGGAGGTCACACCTGTTCCCATGCCGAACACAGAAGTTAAGCTCCTTAGCGTCGATGGTAGTCGAACTTACGTTCCGC       *       NM:i:0  MD:Z:94 NH:i:5  XI:i:0  XA:Z:Q
HWI-ST863:185:C1K2FACXX:7:1101:4006:1911 1:N:0:TGACCAAT 0       Staphylococcus_aureus_HG003     492125  255     94M     *       0       0       TCTGGTGACTATAGCAAGGAGGTCACACCTGTTCCCATGCCGAACACAGAAGTTAAGCTCCTTAGCGTCGATGGTAGTCGAACTTACGTTCCGC       *       NM:i:0  MD:Z:94 NH:i:5  XI:i:1  XA:Z:Q"""

    splice_file = """Staphylococcus_aureus_HG003     5086    5767    splits:1:1:1:N:P        0       +
Staphylococcus_aureus_HG003     7713    8599    splits:1:1:1:N:P        0       +
Staphylococcus_aureus_HG003     8561    9362    splits:3:3:3:N:P        0       +
Staphylococcus_aureus_HG003     12668   12703   splits:1:14:1:N:P       0       +
Staphylococcus_aureus_HG003     17647   17667   splits:1:1:1:N:F        0       +"""

    tran_file = """Staphylococcus_aureus_HG003     1       1       distsplice:Staphylococcus_aureus_HG003:2821337:2:2:2:L:P        0       +
Staphylococcus_aureus_HG003     2821337 2821337 distsplice:Staphylococcus_aureus_HG003:1:2:2:2:R:P      0       +
Staphylococcus_aureus_HG003     2278    2278    distsplice:Staphylococcus_aureus_HG003:450417:1:1:10:L:P        0       +"""

if __name__ == "__main__":
    unittest.main()

