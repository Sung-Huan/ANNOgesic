import sys
import os
import unittest
import shutil
import copy
from io import StringIO
sys.path.append(".")
from mock_gff3 import Create_generator
from mock_helper import import_data, gen_file, extract_info
import annogesiclib.stat_TA_comparison as stc


class TestStatTaComparison(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_read_tas_file(self):
        tss_file = os.path.join(self.test_folder, "aaa_TSS.gff")
        ta_file = os.path.join(self.test_folder, "aaa_transcript.gff")
        gen_file(tss_file, self.example.tss)
        gen_file(ta_file, self.example.ta)
        tsss_uni, tsss, tas_uni, tas = stc.read_tas_file(tss_file, ta_file)
        self.assertEqual(tsss_uni["aaa"][0].start, 2131)
        self.assertEqual(tsss[0].start, 2131)
        self.assertEqual(tas_uni["aaa"][0].start, 313)
        self.assertEqual(tas[0].start, 313)

    def test_detect_tas_region(self):
        out = StringIO()
        out_tss = StringIO()
        stat = stc.detect_tas_region(self.example.tsss, self.example.tas, out, out_tss, 5)
        self.assertDictEqual(stat, {'TSS_with_tran': 1, 'TSS_no_tran': 0, 'no_TSS': 0, 'with_TSS': 1})

    def test_del_attributes(self):
        string = stc.del_attributes(self.example.tsss[0], ["UTR_length", "libs"])
        self.assertDictEqual(string, {'type': 'Primary', 'associated_gene': 'SAOUHSC_00002', 'Name': 'TSS:2131_f', 'ID': 'tss3'})

    def test_assign_tss(self):
        tsss = copy.deepcopy(self.example.tsss)
        trans = copy.deepcopy(self.example.tas)
        stc.assign_tss(tsss[0], trans[0])
        self.assertDictEqual(tsss[0].attributes, {'associated_gene': 'SAOUHSC_00002',
                                                  'Name': 'TSS:2131_f',
                                                  'libs': 'TSB_OD_0.2&TSB_OD_0.5&TSB_t0&pMEM_OD_0.2&pMEM_OD_0.5&pMEM_t2',
                                                  'type': 'Primary', 'Parent_tran': 'tran0',
                                                  'UTR_length': 'Primary_25', 'ID': 'tss3'})
        self.assertDictEqual(trans[0].attributes, {'type': 'cover_CDS&cover_CDS',
                                                   'Name': 'Transcript_00000',
                                                   'associated_cds': 'YP_498609.1&YP_498610.1',
                                                   'associated_tss': 'TSS:313_+&TSS:1641_+&TSS:2128_+&TSS:2131_+&TSS:2131_f',
                                                   'ID': 'tran0'})

    def test_print_tas_stat(self):
        out = StringIO()
        stat = {'TSS_with_tran': 1, 'TSS_no_tran': 0, 'no_TSS': 0, 'with_TSS': 1}
        stc.print_tas_stat(stat, out)
        self.assertEqual(out.getvalue(), self.example.print_tas + "\n")

    def test_compare_tran_tss(self):
        out = StringIO()
        tas = copy.deepcopy(self.example.tas)
        tsss = copy.deepcopy(self.example.tsss)
        stat = {'TSS_with_tran': 1, 'TSS_no_tran': 0, 'no_TSS': 0, 'with_TSS': 1}
        stc.compare_tran_tss(tas, tsss, 5, stat, out)
        self.assertDictEqual(stat, {'with_TSS': 2, 'no_TSS': 0, 'TSS_no_tran': 0, 'TSS_with_tran': 1})

    def test_stat_ta_tss(self):
        tss_file = os.path.join(self.test_folder, "aaa_TSS.gff")
        ta_file = os.path.join(self.test_folder, "aaa_transcript.gff")
        gen_file(tss_file, self.example.tss)
        gen_file(ta_file, self.example.ta)
        stat_file = os.path.join(self.test_folder, "stat")
        out_ta_file = os.path.join(self.test_folder, "out_ta.gff")
        out_tss_file = os.path.join(self.test_folder, "out_tss.gff")
        stc.stat_ta_tss(ta_file, tss_file, stat_file, out_ta_file, out_tss_file, 5)
        datas = import_data(stat_file)
        self.assertEqual("\n".join(datas), "All strains:\n" + self.example.print_tas)
        datas, attributes = extract_info(out_ta_file, "file")
        self.assertListEqual(datas, ['aaa\tfragmented_and_normal\tTranscript\t313\t3344\t.\t+\t.'])
        for attribute in attributes:
            if "associated_tss" in attribute:
                self.assertEqual("associated_tss=TSS:2131_f")
        datas, attributes = extract_info(out_tss_file, "file")
        self.assertListEqual(datas, ['aaa\tTSSpredator\tTSS\t2131\t2131\t.\t+\t.'])
        for attribute in attributes:
            if "Parent_tran" in attribute:
                self.assertEqual(attribute, "Parent_tran=tran0")

    def test_read_tag_file(self):
        gff_file = os.path.join(self.test_folder, "aaa.gff")
        ta_file = os.path.join(self.test_folder, "aaa_transcript.gff")
        gen_file(gff_file, self.example.gff)
        gen_file(ta_file, self.example.ta)
        gffs, tas, stats = stc.read_tag_file(gff_file, ta_file, "gene")
        self.assertEqual(gffs[0].start, 517)
        self.assertEqual(tas[0].start, 313)
        self.assertEqual(stats, {'All': {'gene': 1, 'other': 0, 'bsbe': 0, 'asbe': 0, 'asae': 0, 'bsae': 0},
                                 'aaa': {'gene': 1, 'other': 0, 'bsbe': 0, 'asbe': 0, 'asae': 0, 'bsae': 0}})

    def test_detect_tag_region(self):
        stats = {'All': {'gene': 1, 'other': 0, 'bsbe': 0, 'asbe': 0, 'asae': 0, 'bsae': 0},
                 'aaa': {'gene': 1, 'other': 0, 'bsbe': 0, 'asbe': 0, 'asae': 0, 'bsae': 0}}
        out_t = StringIO()
        out_g = StringIO()
        stc.detect_tag_region(self.example.gffs, self.example.tas, stats, out_t, out_g, "gene")
        self.assertDictEqual(stats, {'All': {'asbe': 1, 'asae': 0, 'other': 0, 'bsbe': 0, 'gene': 1, 'bsae': 0},
                                     'aaa': {'asbe': 1, 'asae': 0, 'other': 0, 'bsbe': 0, 'gene': 1, 'bsae': 0}})

    def test_compare_ta_gff(self):
        tran_type = []
        check = [0, 0, 0, 0, 0]
        stats = {'All': {'asbe': 1, 'asae': 0, 'other': 0, 'bsbe': 0, 'gene': 1, 'bsae': 0},
                 'aaa': {'asbe': 1, 'asae': 0, 'other': 0, 'bsbe': 0, 'gene': 1, 'bsae': 0}}
        stc.compare_ta_gff(self.example.gffs, self.example.tas[0], check, tran_type, False, stats, "gene")
        self.assertListEqual(check, [0, 0, 0, 1, 0])
        self.assertDictEqual(stats, {'aaa': {'asae': 0, 'bsbe': 0, 'asbe': 2, 'bsae': 0, 'other': 0, 'gene': 1},
                                     'All': {'asae': 0, 'bsbe': 0, 'asbe': 2, 'bsae': 0, 'other': 0, 'gene': 1}})

    def test_assign_parent(self):
        gffs = copy.deepcopy(self.example.gffs)
        trans = copy.deepcopy(self.example.tas)
        stc.assign_parent(gffs[0], trans[0])
        self.assertDictEqual(gffs[0].attributes, {'Parent_tran': 'tran0',
                                                  'locus_tag': 'SAOUHSC_00001', 'protein_id': 'YP_498609.1',
                                                  'gene': 'dnaA', 'Name': 'YP_498609.1', 'ID': 'gene0'})
        self.assertDictEqual(trans[0].attributes, {'associated_gene': 'SAOUHSC_00001',
                                                   'associated_tss': 'TSS:313_+&TSS:1641_+&TSS:2128_+&TSS:2131_+',
                                                   'associated_cds': 'YP_498609.1&YP_498610.1', 'type': 'cover_CDS&cover_CDS',
                                                   'Name': 'Transcript_00000', 'ID': 'tran0'})


    def test_print_tag_stat(self):
        out = StringIO()
        stats = {'asae': 0, 'bsbe': 0, 'asbe': 1, 'bsae': 0, 'other': 0, 'gene': 1}
        stc.print_tag_stat(stats, out, 1, "gene")
        self.assertEqual(out.getvalue(), self.example.print_tag + "\n")

    def test_stat_ta_gff(self):
        gff_file = os.path.join(self.test_folder, "aaa.gff")
        ta_file = os.path.join(self.test_folder, "aaa_transcript.gff")
        gen_file(gff_file, self.example.gff)
        gen_file(ta_file, self.example.ta)
        stat_file = os.path.join(self.test_folder, "stat")
        out_ta_file = os.path.join(self.test_folder, "out_ta.gff")
        out_gff_file = os.path.join(self.test_folder, "out.gff")
        stc.stat_ta_gff(ta_file, gff_file, stat_file, out_ta_file, out_gff_file, "gene")
        datas = import_data(stat_file)
        self.assertEqual("\n".join(datas), "All strains:\nThe transcriptome assembly information compares with annotation gff file:\n" + \
                                            self.example.print_tag)
        datas, attributes = extract_info(out_ta_file, "file")
        self.assertListEqual(datas, ['aaa\tfragmented_and_normal\tTranscript\t313\t3344\t.\t+\t.'])
        for attribute in attributes:
            if "type" in attribute:
                self.assertEqual(attribute, "type=cover_CDS")
            if "associated_cds=" in attribute:
                self.assertEqual(attribute, "associated_cds=YP_498609.1")
        datas, attributes = extract_info(out_gff_file, "file")
        self.assertListEqual(datas, ['aaa\tRefseq\tgene\t517\t1878\t.\t+\t.',
                                     'aaa\tRefseq\tCDS\t517\t1878\t.\t+\t.'])
        for attribute in attributes:
            if "Parent_tran" in attribute:
                self.assertEqual(attribute, "Parent_tran=tran0")


class Example(object):

    tss = """aaa	TSSpredator	TSS	2131	2131	.	+	.	UTR_length=Primary_25;type=Primary;ID=tss3;libs=TSB_OD_0.2&TSB_OD_0.5&TSB_t0&pMEM_OD_0.2&pMEM_OD_0.5&pMEM_t2;associated_gene=SAOUHSC_00002;Name=TSS:2131_f"""
    ta = """aaa	fragmented_and_normal	Transcript	313	3344	.	+	.	associated_tss=TSS:313_+&TSS:1641_+&TSS:2128_+&TSS:2131_+;type=cover_CDS&cover_CDS;Name=Transcript_00000;ID=tran0;associated_cds=YP_498609.1&YP_498610.1"""
    gff = """aaa	Refseq	gene	517	1878	.	+	.	protein_id=YP_498609.1;gene=dnaA;ID=cds0;Parent=gene0;locus_tag=SAOUHSC_00001;Name=YP_498609.1
aaa	Refseq	CDS	517	1878	.	+	.	protein_id=YP_498609.1;gene=dnaA;ID=cds0;Parent=gene0;locus_tag=SAOUHSC_00001;Name=YP_498609.1"""
    tss_dict = [{"seq_id": "aaa", "source": "TSSpredator", "feature": "TSS", "start": 2131,
                 "end": 2131, "phase": ".", "strand": "+", "score": "."}]
    attributes_tss = [{"ID": "tss3", "Name": "TSS:2131_f", "UTR_length": "Primary_25",
                       "type": "Primary", "associated_gene": "SAOUHSC_00002",
                       "libs": "TSB_OD_0.2&TSB_OD_0.5&TSB_t0&pMEM_OD_0.2&pMEM_OD_0.5&pMEM_t2"}]
    tsss = []
    tsss.append(Create_generator(tss_dict[0], attributes_tss[0], "gff"))
    ta_dict = [{"seq_id": "aaa", "source": "fragmented_and_normal", "feature": "transcript",
                "start": 313, "end": 3344, "phase": ".", "strand": "+", "score": "."}]
    attributes_ta = [{"ID": "tran0", "Name": "Transcript_00000", "associated_tss": "TSS:313_+&TSS:1641_+&TSS:2128_+&TSS:2131_+",
                      "type": "cover_CDS&cover_CDS", "associated_cds": "YP_498609.1&YP_498610.1"}]
    tas = []
    tas.append(Create_generator(ta_dict[0], attributes_ta[0], "gff"))
    gff_dict = [{"seq_id": "aaa", "source": "RefSeq", "feature": "gene",
                "start": 517, "end": 1878, "phase": ".", "strand": "+", "score": "."}]
    attributes_gff = [{"ID": "gene0", "Name": "YP_498609.1", "locus_tag": "SAOUHSC_00001",
                       "gene": "dnaA", "protein_id": "YP_498609.1"}]
    gffs = []
    gffs.append(Create_generator(gff_dict[0], attributes_gff[0], "gff"))
    print_tas = """	Transcript starts or overlap with TSS:1 (1.0)
	Transcript has no relationship with TSS:0 (0.0)
	TSS starts or overlap with transcript:1 (1.0)
	TSS has no relationship with transcript:0 (0.0)"""
    print_tag = """	Transcript starts before and ends after gene:1 (1.0)
	Transcript starts after and ends before gene:0 (0.0)
	Transcript starts before and ends before gene:0 (0.0)
	Transcript starts after and ends after gene:0 (0.0)
	Transcript has no overlap of gene:0 (0.0)
	Total genes which have gene expression:1 (1.0)"""


if __name__ == "__main__":
    unittest.main()

