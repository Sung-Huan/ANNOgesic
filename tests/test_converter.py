import sys
import os
import csv
import shutil
import unittest
from collections import defaultdict
from io import StringIO
sys.path.append(".")
from mock_helper import import_data
from annogesiclib.converter import Converter

class Mock_gff3_parser(object):

    def __init__(entry_dict):
        entry_dict = Example().gff_dict

    def entries(self, entry_dict):  
        return Mock_gff3_entry(entry_dict)


class Mock_gff3_entry(object):

    def __init__(self, entry_dict):
        self.seq_id = entry_dict["seq_id"]
        self.source = entry_dict["source"]
        self.feature = entry_dict["feature"]
        # 1-based coordinates
        # Make sure that start <= end
        start, end = sorted([int(entry_dict["start"]), int(entry_dict["end"])])
        self.start = start
        self.end = end
        self.score = entry_dict["score"]
        self.strand = entry_dict["strand"]
        self.phase = entry_dict["phase"]
        self.attributes = entry_dict["attributes"]

class Mock_TSSPredatorReader(object):


    def entries(self, input_fh):
        for row in csv.reader(input_fh, delimiter="\t"):
            if not row[0].startswith("SuperPos"):
                yield Mock_TSSPredatorEntry(row)

class Mock_TSSPredatorEntry(object):

    def __init__(self, row):
        assert(len(row) == 28)
        self.super_pos = int(row[0])
        self.super_strand = row[1]
        self.map_count = int(row[2])
        self.det_count = int(row[3])
        self.genome = row[4]
        self.is_detected = True if row[5] == "1" else False
        self.is_enriched = True if row[6] == "1" else False
        self.step_heigth = row[7]
        self.step_factor = row[8]
        self.enrichment_factor = row[9]
        self.class_count = int(row[10])
        self.pos = int(row[11])
        self.strand = row[12]
        self.locus_tag = row[13]
        self.srna_asrna = row[14]
        self.product = row[15]
        self.utr_length = row[16]
        self.gene_length = row[17]
        self.is_primary = True if row[18] == "1" else False
        self.is_secondary = True if row[19] == "1" else False
        self.is_internal = True if row[20] == "1" else False
        self.is_antisense = True if row[21] == "1" else False
        self.is_automated = True if row[22] == "1" else False
        self.is_manual = True if row[23] == "1" else False
        self.is_putative_srna = True if row[24] == "1" else False
        self.is_putative_asrna = True if row[25] == "1" else False
        self.comment = row[26]
        self.seq = row[27]
        self.is_orphan = False
        if (self.is_primary is False and self.is_secondary is False and
            self.is_internal is False and self.is_antisense is False):
            self.is_orphan = True

class Mock_func(object):

    def __init__(self):
        self.gff_dict = Example().gff_dict

    def print_rntptt_title(self, out, num, seq_id, length):
        pass

    def mock_read_file(self, gff_file, fasta_file, rnas, cdss, genes):
        seq = "GCAGGTTGAGTTCCTGTTCCCGATAGATCCGATAAACCCGCTTATGATTCCA"
        num_cds = 0
        num_rna = 0
        for gff in self.gff_dict:
            if gff["feature"] == "gene":
                genes.append(Mock_gff3_parser().entries(gff))
            elif gff["feature"] == "CDS":
                cdss.append(Mock_gff3_parser().entries(gff))
                num_cds += 1
            elif gff["feature"] == "tRNA":
                rnas.append(Mock_gff3_parser().entries(gff))
                num_rna += 1
        return (num_cds, num_rna, seq)

class TestConverter(unittest.TestCase):

    def setUp(self):
        self.converter = Converter()
        self.example = Example()
        self.converter.gff3parser = Mock_gff3_parser
        self.converter._print_rntptt_title = Mock_func().print_rntptt_title
        self.converter.tsspredator = Mock_TSSPredatorReader()
        self.converter._read_file = Mock_func().mock_read_file
        self.gff_file = self.example.gff_file
        self.ptt_out = self.example.ptt_out
        self.rnt_out = self.example.rnt_out
        self.srna_out = self.example.srna_out
        self.embl_file = self.example.embl_file
        self.embl_out = self.example.embl_out
        self.multi_embl = self.example.multi_embl
        self.gff_out = self.example.gff_out
        self.mastertable = self.example.mastertable
        self.tss_file = self.example.tss_file
        self.fasta_file = self.example.fasta_file
        self.transterm = self.example.transterm
        self.term_file = self.example.term_file
        self.circ_file = self.example.circrna_table
        self.circ_all = self.example.circrna_all
        self.circ_best = self.example.circrna_best
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_print_rntptt_file(self):
        cdss = []
        genes = []
        rnas = []
        gff_dict = Example().gff_dict
        for gff in gff_dict:
            if gff["feature"] == "gene":
                genes.append(self.converter.gff3parser.entries(self, gff))
            elif gff["feature"] == "CDS":
                cdss.append(self.converter.gff3parser.entries(self, gff))
            elif gff["feature"] == "tRNA":
                rnas.append(self.converter.gff3parser.entries(self, gff))
        out_p = StringIO()
        out_r = StringIO()
        self.converter._print_rntptt_file(out_p, cdss, genes)
        self.converter._print_rntptt_file(out_r, rnas, genes)
        self.assertEqual(out_p.getvalue().split("\n")[:-1], self.ptt_out.split("\n"))
        self.assertEqual(out_r.getvalue().split("\n")[:-1], self.rnt_out.split("\n"))
        out_p.close()
        out_r.close()

    def test_srna2pttrnt(self):
        srna_input_file = os.path.join(self.test_folder, "srna.gff")
        srna_output_file = os.path.join(self.test_folder, "srna.out")
        with open(srna_input_file, "w") as fh:
            fh.write(self.gff_file)
        srnas = []
        self.converter._srna2rntptt(srna_input_file, srna_output_file, srnas, 1234567)
        datas = import_data(srna_output_file)
        self.assertEqual(set(datas), set(self.srna_out.split("\n")))

    def test_multi_embl_pos(self):
        embls = []
        for line in self.embl_file.split("\n"):
            datas = self.converter._multi_embl_pos(line.strip())
            if datas != "Wrong":
                embls.append(datas)
        for index in range(0, 7):
            self.assertDictEqual(embls[index], self.embl_out[index])
        for index in range(0, 2):
            self.assertDictEqual(embls[-1]["pos"][index], self.multi_embl[index])
        
    def test_parser_embl_data(self):
        embl_file = os.path.join(self.test_folder, "test.embl")
        embl_out = os.path.join(self.test_folder, "test.embl_out")
        out = StringIO()
        with open(embl_file, "w") as eh:
            for line in self.embl_file.split("\n"):
                eh.write(line + "\n")
        info = self.converter._parser_embl_data(embl_file, out)
        datas = out.getvalue().split("\n")
        self.assertEqual(set(datas[:-1]), set(self.gff_out.split("\n")))
        self.assertEqual(info[0], "NC_007795.1")
        for index in range(0, 2):
            self.assertDictEqual(info[1]["pos"][index], self.multi_embl[index])
        out.close()

    def test_multi_tss_class(self):
        nums = {"tss": 0, "tss_uni": 0, "class": 1}
        utrs = {"total": [], "pri": [], "sec": []}
        tss_features = {"tss_types": [], "locus_tags": [], "utr_lengths": []}
        tss_index = defaultdict(lambda: 0)
        master_file = os.path.join(self.test_folder, "test.tsv")
        fh = StringIO(self.mastertable)
        for tss in self.converter.tsspredator.entries(fh):
            self.converter._multi_tss_class(tss, tss_index, tss_features, nums, utrs)
        fh.close()
        self.assertDictEqual(nums, {'tss_uni': 0, 'class': 5, 'tss': 2})

    def test_convert_mastertable2gff(self):
        master_file = os.path.join(self.test_folder, "test.tsv")
        with open(master_file, "w") as th:
            th.write(self.mastertable)
        out_gff = os.path.join(self.test_folder, "test.tsv_out")
        self.converter.convert_mastertable2gff(master_file, "TSSpredator", "TSS",
                                               "aaa", out_gff)
        datas = import_data(out_gff)
        self.assertEqual(set(datas), set(self.tss_file.split("\n")))

    def test_convert_gff2rntptt(self):
        srna_input_file = os.path.join(self.test_folder, "srna.gff")
        srna_output_file = os.path.join(self.test_folder, "srna.out")
        gff_file = os.path.join(self.test_folder, "test.gff")
        rnt_file = os.path.join(self.test_folder, "test.rnt")
        ptt_file = os.path.join(self.test_folder, "test.ptt")
        fasta_file = os.path.join(self.test_folder, "test.fa")
        with open(srna_input_file, "w") as fh:
            fh.write(self.gff_file)
        with open(gff_file, "w") as fh:
            fh.write(self.gff_file)
        with open(fasta_file, "w") as fh:
            fh.write(self.fasta_file)
        self.converter.convert_gff2rntptt(gff_file, fasta_file, ptt_file, rnt_file,
                                          srna_input_file, srna_output_file)
        self.assertTrue(srna_output_file)
        self.assertTrue(rnt_file)
        self.assertTrue(ptt_file)


    def test_convert_embl2gff(self):
        embl_file = os.path.join(self.test_folder, "test.embl")
        gff_file = os.path.join(self.test_folder, "test.embl_out")
        with open(embl_file, "w") as eh:
            for line in self.embl_file.split("\n"):
                eh.write(line + "\n")
        self.converter.convert_embl2gff(embl_file, gff_file)
        datas = import_data(gff_file)
        self.assertEqual(set(datas[1:-2]), set(self.gff_out.split("\n")))

    def test_convert_transtermhp2gff(self):
        transterm_file = os.path.join(self.test_folder, "test_best_terminator_after_gene.bag")
        gff_file = os.path.join(self.test_folder, "transterm.gff")
        with open(transterm_file, "w") as th:
            th.write(self.transterm)
        self.converter.convert_transtermhp2gff(transterm_file, gff_file)
        datas = import_data(gff_file)
        self.assertEqual(set(datas), set(self.term_file.split("\n")))

    def test_convert_circ2gff(self):
        circ_file = os.path.join(self.test_folder, "circ.csv")
        out_all = os.path.join(self.test_folder, "all.gff")
        out_filter = os.path.join(self.test_folder, "best.gff")  
        with open(circ_file, "w") as ch:
            ch.write(self.circ_file)
        self.converter.convert_circ2gff(circ_file, 5, 0.5, 0.5,
                                        out_all, out_filter)
        datas = import_data(out_all)
        self.assertEqual(set(datas), set(self.circ_all.split("\n")))
        datas = import_data(out_filter)
        self.assertEqual(set(datas), set(self.circ_best.split("\n")))

class Example(object):

    gff_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "gene", "start": 517, "end": 1878,
                 "phase": ".", "strand": "+", "score": ".", "attributes": {"Name": "dnaA", "locus_tag": "AAA_00001",
                 "gene": "dnaA", "ID": "gene0", "dbxref": "GeneID:3919798"}},
                {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 517, "end": 1878,
                 "phase": ".", "strand": "+", "score": ".", "attributes": {"Name": "YP_498609.1",
                 "ID": "cds0", "Parent": "gene0", "product": "chromosomal replication initiation protein",
                 "protein_id": "YP_498609.1"}},
                {"seq_id": "aaa", "source": "Refseq", "feature": "gene", "start": 2156, "end": 3289,
                 "phase": ".", "strand": "-", "score": ".", "attributes": {"Name": "AAA_00002", "locus_tag": "AAA_00002",
                 "ID": "gene1", "dbxref": "GeneID:3919799"}},
                {"seq_id": "aaa", "source": "Refseq", "feature": "CDS", "start": 2156, "end": 3289,
                 "phase": ".", "strand": "-", "score": ".", "attributes": {"Name": "YP_498610.1",
                 "ID": "cds1", "locus_tag": "AAA_00002", "protein_id": "YP_498610.1"}},
                {"seq_id": "aaa", "source": "Refseq", "feature": "gene", "start": 4444, "end": 5444,
                 "phase": ".", "strand": "+", "score": ".", "attributes": {"Name": "AAA_T00004", "locus_tag": "AAA_T00004",
                 "ID": "gene2"}},
                {"seq_id": "aaa", "source": "Refseq", "feature": "tRNA", "start": 4444, "end": 5444,
                 "phase": ".", "strand": "+", "score": ".", "attributes": {"Name": "AAA_T00004",
                 "ID": "rna0", "locus_tag": "AAA_T00004"}}]

    srna_dict = [{"seq_id": "aaa", "source": "Refseq", "feature": "sRNA", "start": 12, "end": 331,
                 "phase": ".", "strand": "+", "score": ".", "attributes": {"Name": "sRNA_candidate_0001", "ID": "srna0"}},
                 {"seq_id": "aaa", "source": "Refseq", "feature": "tRNA", "start": 100, "end": 244,
                  "phase": ".", "strand": "+", "score": ".", "attributes": {"Name": "sRNA_candidate_0002", "ID": "rna0"}}]

    ptt_out = """517..1878	+	1362	YP_498609.1	dnaA	AAA_00001	-	-	chromosomal replication initiation protein
2156..3289	-	1134	YP_498610.1	-	AAA_00002	-	-	-"""

    rnt_out = """4444..5444	+	1001	-	-	AAA_T00004	-	-	-"""

    gff_file = """
aaa	Refseq	gene	517	1878	.	+	.	Name=dnaA;locus_tag=AAA_00001;gene=dnaA;ID=gene0;db_xref=GeneID:3919798
aaa	Refseq	CDS	517	1878	.	+	.	protein_id=YP_498609.1;ID=cds0;Name=YP_498609.1;product=chromosomal replication initiation protein;Parent=gene0
aaa	Refseq	gene	2156	3289	.	-	.	Name=AAA_00002;locus_tag=AAA_00002;ID=gene1;db_xref=GeneID:3919799
aaa	Refseq	CDS	2156	3289	.	-	.	protein_id=YP_498610.1;ID=cds1;Name=YP_498610.1;locus_tag=AAA_00002
aaa	Refseq	gene	4444	5444	.	+	.	Name=AAA_T00004;locus_tag=AAA_T00004;ID=gene2
aaa	Refseq	tRNA	4444	5444	.	+	.	Name=AAA_T00018;locus_tag=AAA_T00004;ID=rna0"""

    srna_out = """517..1878\t+\t1362\tncRNA_00002\t-\tncRNA_00002\t-\t-\tsRNA
517..1878\t+\t1362\tncRNA_00001\t-\tncRNA_00001\t-\t-\tsRNA
2156..3289\t-\t1134\tncRNA_00003\t-\tncRNA_00003\t-\t-\tsRNA
4444..5444\t+\t1001\tncRNA_00005\t-\tncRNA_00005\t-\t-\tsRNA
4444..5444\t+\t1001\tncRNA_00006\t-\tncRNA_00006\t-\t-\tsRNA
2156..3289\t-\t1134\tncRNA_00004\t-\tncRNA_00004\t-\t-\tsRNA"""

    embl_out = [{'pos': [{'end': '2821361', 'start': '1'}], 'source': 'source', 'strand': '+'},
{'pos': [{'end': '1878', 'start': '517'}], 'source': 'gene', 'strand': '+'},
{'pos': [{'end': '1878', 'start': '517'}], 'source': 'CDS', 'strand': '+'},
{'pos': [{'end': '1872', 'start': '520'}], 'source': 'misc_feature', 'strand': '+'},
{'pos': [{'end': '10456', 'start': '9755'}], 'source': 'gene', 'strand': '-'},
{'pos': [{'end': '10456', 'start': '9755'}], 'source': 'CDS', 'strand': '-'},
{'pos': [{'end': '1491304', 'start': '1491014'}], 'source': 'gene', 'strand': '-'}]

    multi_embl = [{'end': '1491245', 'start': '1491015'}, {'end': '1491304', 'start': '1491251'}]

    srna_file = """Staphylococcus_aureus_HG003	intergenic	sRNA	313	417	.	+	.	ID=srna0;Name=sRNA_candidate_00001"""

    fasta_file = """>Staphylococcus_aureus_HG003
GCAGGTTGAGTTCCTGTTCCCGATAGATCCGATAAACCCGCTTATGATTCCAGAGCTGTCCCTGCACATA
GCCGCTTCATTTTTTTCCAAGGGCTTCCTTCAGGATATCCGTCTGCATGCTCAAATCCGCATACATGCGG"""

    embl_file = """
ID   NC_007795; SV 1; ; DNA; ; UNC; 2821361 BP.
XX
AC   NC_007795;
XX
DE   Staphylococcus aureus subsp. aureus NCTC 8325 chromosome, complete genome.
XX
KW   RefSeq
XX
OS   Staphylococcus aureus subsp. aureus NCTC 8325
OC   Bacteria; Firmicutes; Bacilli; Bacillales; Staphylococcus.
XX
RN   [1]
RP   1-2821361
RA   Gillaspy,A.F., Worrell,V., Orvis,J., Roe,B.A., Dyer,D.W. and Iandolo,J.J.;
RT   "The Staphylococcus aureus NCTC8325 Genome";
RL   (in) Fischetti,V., Novick,R., Ferretti,J., Portnoy,D. and Rood,J. (Eds.);
RL   GRAM POSITIVE PATHOGENS; ASM Press (2006)
XX
RN   [2]
RP   1-2821361
RG   NCBI Genome Project
RT   "Direct Submission";
RL   Submitted (18-FEB-2006) National Center for Biotechnology Information, NIH,
RL   Bethesda, MD 20894, USA
XX
RN   [3]
RP   1-2821361
RA   Gillaspy,A.F., Worrell,V., Orvis,J., Roe,B.A., Dyer,D.W. and Iandolo,J.J.;
RT   "Direct Submission";
RL   Submitted (27-JAN-2006) Microbiology and Immunology, The University of
RL   Oklahoma Health Sciences Center, 940 Stanton L. Young Boulevard, Oklahoma
RL   City, OK 73104, USA
XX
CC   REVIEWED REFSEQ: This record has been curated by NCBI staff. The
CC   reference sequence was derived from CP000253.
CC   RefSeq Category: Reference Genome
CC   UPR: UniProt Genome
CC   Staphylococcus aureus subsp. aureus NCTC 8325 is available from
CC   www.narsa.net.
CC   COMPLETENESS: full length.
XX
FH   Key             Location/Qualifiers
FT   source          1..2821361
FT                   /db_xref="taxon:93061"
FT                   /mol_type="genomic DNA"
FT                   /strain="NCTC 8325"
FT                   /sub_species="aureus"
FT                   /organism="Staphylococcus aureus subsp. aureus NCTC 8325"
FT   gene            517..1878
FT                   /db_xref="GeneID:3919798"
FT                   /locus_tag="SAOUHSC_00001"
FT                   /gene="dnaA"
FT   CDS             517..1878
FT                   /codon_start=1
FT                   /product="chromosomal replication initiation protein"
FT                   /db_xref="GI:88193824"
FT                   /db_xref="GeneID:3919798"
FT                   /note="binds to the dnaA-box as an ATP-bound complex at the
FT                   origin of replication during the initiation of chromosomal
FT                   replication; can also affect transcription of multiple
FT                   genes including itself."
FT                   /translation="MSEKEIWEKVLEIAQEKLSAVSYSTFLKDTELYTIKDGEAIVLSS
FT                   IPFNANWLNQQYAEIIQAILFDVVGYEVKPHFITTEELANYSNNETATPKETTKPSTET
FT                   TEDNHVLGREQFNAHNTFDTFVIGPGNRFPHAASLAVAEAPAKAYNPLFIYGGVGLGKT
FT                   HLMHAIGHHVLDNNPDAKVIYTSSEKFTNEFIKSIRDNEGEAFRERYRNIDVLLIDDIQ
FT                   FIQNKVQTQEEFFYTFNELHQNNKQIVISSDRPPKEIAQLEDRLRSRFEWGLIVDITPP
FT                   DYETRMAILQKKIEEEKLDIPPEALNYIANQIQSNIRELEGALTRLLAYSQLLGKPITT
FT                   ELTAEALKDIIQAPKSKKITIQDIQKIVGQYYNVRIEDFSAKKRTKSIAYPRQIAMYLS
FT                   RELTDFSLPKIGEEFGGRDHTTVIHAHEKISKDLKEDPIFKQEVENLEKEIRNV"
FT                   /protein_id="REF_uohsc:SAOUHSC00001"
FT                   /protein_id="YP_498609.1"
FT                   /gene="dnaA"
FT                   /locus_tag="SAOUHSC_00001"
FT                   /transl_table=11
FT   misc_feature    520..1872
FT                   /db_xref="CDD:234667"
FT                   /locus_tag="SAOUHSC_00001"
FT                   /gene="dnaA"
FT                   /note="chromosomal replication initiation protein;
FT                   Reviewed; Region: dnaA; PRK00149"
FT   gene            complement(9755..10456)
FT                   /db_xref="GeneID:3919180"
FT                   /locus_tag="SAOUHSC_00007"
FT   CDS             complement(9755..10456)
FT                   /codon_start=1
FT                   /product="hypothetical protein"
FT                   /db_xref="GI:88193830"
FT                   /db_xref="GeneID:3919180"
FT                   /note="conserved hypothetical protein"
FT                   /translation="MLAARACVFSGSGLITVATHPTNHSALHSRCPEAMVIDINDTKML
FT                   TKMIEMTDSILIGPGLGVDFKGNNAITFLLQNIQPHQNLIVDGDAITIFSKLKPQLPTC
FT                   RVIFTPHLKEWERLSGIPIEEQTYERNREAVDRLGATVVLKKHGTEIFFKDEDFKLTIG
FT                   SPAMATGGMGDTLAGMITSFVGQFDNLKEAVMSATYTHSFIGENLAKDMYVVPPSRLIN
FT                   EIPYAMKQLES"
FT                   /protein_id="REF_uohsc:SAOUHSC00007"
FT                   /protein_id="YP_498615.1"
FT                   /locus_tag="SAOUHSC_00007"
FT                   /transl_table=11
FT   gene            complement(1491014..1491304)
FT                   /locus_tag="SAOUHSC_01545"
FT                   /db_xref="GeneID:3920588"
FT   CDS             complement(join(1491015..1491245,1491251..1491304))
FT                   /locus_tag="SAOUHSC_01545"
FT                   /protein_id="REF_uohsc:SAOUHSC01545"
FT                   /protein_id="YP_500061.1"
FT                   /transl_table=11
FT                   /product="hypothetical protein"
FT                   /db_xref="GI:88195258"
FT                   /db_xref="GeneID:3920588"
FT                   /note="conserved hypothetical phage protein"
FT                   /codon_start=1
"""

    gff_out = """NC_007795.1	Refseq	source	1	2821361	.	+	.	db_xref=taxon:93061;mol_type=genomic DNA;strain=NCTC 8325;sub_species=aureus;organism=Staphylococcus aureus subsp. aureus NCTC 8325
NC_007795.1	Refseq	gene	517	1878	.	+	.	db_xref=GeneID:3919798;locus_tag=SAOUHSC_00001;gene=dnaA
NC_007795.1	Refseq	CDS	517	1878	.	+	.	codon_start=1;product=chromosomal replication initiation protein;db_xref=GI:88193824;db_xref=GeneID:3919798;note=binds to the dnaA-box as an ATP-bound complex at the origin of replication during the initiation of chromosomal replication, can also affect transcription of multiple genes including itself.;translation=MSEKEIWEKVLEIAQEKLSAVSYSTFLKDTELYTIKDGEAIVLSS IPFNANWLNQQYAEIIQAILFDVVGYEVKPHFITTEELANYSNNETATPKETTKPSTET TEDNHVLGREQFNAHNTFDTFVIGPGNRFPHAASLAVAEAPAKAYNPLFIYGGVGLGKT HLMHAIGHHVLDNNPDAKVIYTSSEKFTNEFIKSIRDNEGEAFRERYRNIDVLLIDDIQ FIQNKVQTQEEFFYTFNELHQNNKQIVISSDRPPKEIAQLEDRLRSRFEWGLIVDITPP DYETRMAILQKKIEEEKLDIPPEALNYIANQIQSNIRELEGALTRLLAYSQLLGKPITT ELTAEALKDIIQAPKSKKITIQDIQKIVGQYYNVRIEDFSAKKRTKSIAYPRQIAMYLS RELTDFSLPKIGEEFGGRDHTTVIHAHEKISKDLKEDPIFKQEVENLEKEIRNV;protein_id=REF_uohsc:SAOUHSC00001;protein_id=YP_498609.1;gene=dnaA;locus_tag=SAOUHSC_00001;transl_table=11
NC_007795.1	Refseq	gene	9755	10456	.	-	.	db_xref=GeneID:3919180;locus_tag=SAOUHSC_00007
NC_007795.1	Refseq	CDS	9755	10456	.	-	.	codon_start=1;product=hypothetical protein;db_xref=GI:88193830;db_xref=GeneID:3919180;note=conserved hypothetical protein;translation=MLAARACVFSGSGLITVATHPTNHSALHSRCPEAMVIDINDTKML TKMIEMTDSILIGPGLGVDFKGNNAITFLLQNIQPHQNLIVDGDAITIFSKLKPQLPTC RVIFTPHLKEWERLSGIPIEEQTYERNREAVDRLGATVVLKKHGTEIFFKDEDFKLTIG SPAMATGGMGDTLAGMITSFVGQFDNLKEAVMSATYTHSFIGENLAKDMYVVPPSRLIN EIPYAMKQLES;protein_id=REF_uohsc:SAOUHSC00007;protein_id=YP_498615.1;locus_tag=SAOUHSC_00007;transl_table=11
NC_007795.1	Refseq	gene	1491014	1491304	.	-	.	locus_tag=SAOUHSC_01545;db_xref=GeneID:3920588"""

    mastertable = """SuperPos	SuperStrand	mapCount	detCount	Genome	detected	enriched	stepHeight	stepFactor	enrichmentFactor	classCount	Pos	Strand	Locus_tag	sRNA/asRNA	Product	UTRlength	GeneLength	Primary	Secondary	Internal	Antisense	Automated	Manual	Putative sRNA	Putative asRNA	Comment	Sequence -50 nt upstream + TSS (51nt)
313	+	2	2	TSB_OD_0.2	1	1	83.64	>100	4.34	1	313	+	SAOUHSC_00001		chromosomal replication initiation protein	204	1362	1	0	0	0	1	0	0	0		ATATTCACAGGTATTTGACATATAGAGAACTGAAAAAGTATAATTGTGTGG
313	+	2	2	TSB_OD_0.5	1	1	84.94	>100	4.81	1	313	+	SAOUHSC_00001		chromosomal replication initiation protein	204	1362	1	0	0	0	1	0	0	0		ATATTCACAGGTATTTGACATATAGAGAACTGAAAAAGTATAATTGTGTGG
2131	+	2	1	TSB_OD_0.2	1	1	27.88	4.58	0.4	2	2131	+	SAOUHSC_00002		DNA polymerase III subunit beta	25	1134	1	0	0	0	1	0	0	0		ACAGCACCTACTACTATTACTAAGAACTTAAAACCTATATAATTATATATA
2131	+	2	1	TSB_OD_0.2	1	1	27.88	4.58	0.4 	2	2131	+	SAOUHSC_00003		DNA polymerase III subunit beta	NA	1134	0	0	1	0	1	0	0	0		ACAGCACCTACTACTATTACTAAGAACTTAAAACCTATATAATTATATATA
2131	+	2	1	TSB_OD_0.5	0	1	31.77	5.45	0.48	2	2131	+	SAOUHSC_00002		DNA polymerase III subunit beta	25	1134	1	0	0	0	1	0	0	0		ACAGCACCTACTACTATTACTAAGAACTTAAAACCTATATAATTATATATA
2131	+	2	1	TSB_OD_0.5	0	1	31.77	5.45	0.48	2	2131	+	SAOUHSC_00003		DNA polymerase III subunit beta	NA	1134	0	0	1	0	1	0	0	0		ACAGCACCTACTACTATTACTAAGAACTTAAAACCTATATAATTATATATA"""

    tss_file = """##gff-version 3
aaa	TSSpredator	TSS	313	313	.	+	.	Name=TSS:313_f;ID=tss0;type=Primary;UTR_length=Primary_204;associated_gene=SAOUHSC_00001;libs=TSB_OD_0.2&TSB_OD_0.5
aaa	TSSpredator	TSS	2131	2131	.	+	.	Name=TSS:2131_f;ID=tss1;type=Primary&Internal;UTR_length=Primary_25&Internal_NA;associated_gene=SAOUHSC_00002&SAOUHSC_00003;libs=TSB_OD_0.2"""

    transterm = """        recF NONE
SAOUHSC_00005 NONE
SAOUHSC_00006    9676 ..    9705 + -14.5  -4.91109   AAGAATAATAAAAAA       TAAGACTTCCCTA TATG TAGGGGAGTCTTA        TTTTTATGCTAGAAA 100 33
SAOUHSC_00009   14156 ..   14177 +  -8.7  -5.69473   ATAATATTTTAAAAA           GTGGTGACG AAGC TGTCGCCAC            TTTTTTTGTGCTGTA 100 109
SAOUHSC_00010 NONE
SAOUHSC_00014 NONE"""

    term_file = """##gff-version 3
test	TransTermHP	terminator	9676	9705	.	+	.	associated_gene=SAOUHSC_00006;ID=term0;Name=terminator_00000
test	TransTermHP	terminator	14156	14177	.	+	.	associated_gene=SAOUHSC_00009;ID=term1;Name=terminator_00001"""

    circrna_table = """ID	strain	strand	start	end	annotation_overlap	supported_reads	supported_reads/reads_at_start	supported_reads/reads_at_end
circRNA_0	Staphylococcus_aureus_HG003	+	497897	498038	SAOUHSC_R0007	36	0.56822429906542055	0.52040133779264214
circRNA_1	Staphylococcus_aureus_HG003	+	492193	495089	NA	32	1.0	1.0
circRNA_2	Staphylococcus_aureus_HG003	+	492	4956	NA	2	0.5106796116504854	0.11940298507462686
circRNA_3	Staphylococcus_aureus_HG003	+	23442	49504	NA	52	0.2106796116504854	0.81940298507462686
"""

    circrna_all = """##gff-version 3
Staphylococcus_aureus_HG003	segemehl	circRNA	492	4956	.	+	.	ID=circrna2;name=circRNA_2;support_reads=2;read_at_start=0.5106796116504854;read_at_end=0.11940298507462686;confliction=NA
Staphylococcus_aureus_HG003	segemehl	circRNA	23442	49504	.	+	.	ID=circrna3;name=circRNA_3;support_reads=52;read_at_start=0.2106796116504854;read_at_end=0.8194029850746268;confliction=NA
Staphylococcus_aureus_HG003	segemehl	circRNA	492193	495089	.	+	.	ID=circrna1;name=circRNA_1;support_reads=32;read_at_start=1.0;read_at_end=1.0;confliction=NA
Staphylococcus_aureus_HG003	segemehl	circRNA	497897	498038	.	+	.	ID=circrna0;name=circRNA_0;support_reads=36;read_at_start=0.5682242990654206;read_at_end=0.5204013377926422;confliction=SAOUHSC_R0007"""

    circrna_best = """##gff-version 3
Staphylococcus_aureus_HG003	segemehl	circRNA	492193	495089	.	+	.	ID=circrna1;name=circRNA_1;support_reads=32;read_at_start=1.0;read_at_end=1.0;confliction=NA"""
if __name__ == "__main__":
    unittest.main()

