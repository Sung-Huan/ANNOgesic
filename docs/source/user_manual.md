# User manual of ANNOgesic
Authors: Sung-Huan Yu …

## Introduction
RNA-Seq is a potent method mainly used to measure expression
levels of organisms. However, it has also become a powerful
tool to improve genome annotations of numerous organisms and 
the detection of new non-coding transcripts which are otherwise
hard to predict computationally. Still, there is a need for tools
and platforms for genome annotation that directly integrate all 
information to produce a precise and high-resolution annotation.
In this study, we tested and created several tools that cover
different aspects of the genome annotation. In order to get the 
best results, we also optimized the parameters of some tools.
We developed a tool-kit – ANNOgesic, that should cover all required
functions for the annotation of a genome/transcriptome. TransAP
can automatically generate high-quality annotation information for 
query strains. TransAP is modular and its subcommands can be separately
used. Currently, the pipeline already can detect or integrate a) fasta, CDS,
tRNA, rRNA and genes of query genome, b) transcription starting sites
(TSSs), c) rho-independent terminators, d) transcript assembly,
e) untranslated region (UTRs), f) sRNA, g) promoters, h) processing sites,
i) circular RNAs, j) protein-protein interaction networks, k) potential
sRNA target, l) single-nucleotide polymorphism (SNP), m) operons,
n) GO terms, o) subcellular localization, p) riboswitch, q) potential sORF.

## Installation

### Download ANNOgesic

git clone https://github.com/Sung-Huan/ANNOgesic.git

### Pre-required tools, packages and database

Python <https://www.python.org/> : version higher or equal than 3.4

BioPython <http://biopython.org/wiki/Main_Page>

BioPerl <http://www.bioperl.org/wiki/Main_Page>

Blast+ <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/> : version higher or equal than 2.2.28+ (for sRNA_detection)

EMBOSS <http://emboss.sourceforge.net/> : version higher or equal than 6.6.0 (for subcellular localization detection)

ImageMagick <http://www.imagemagick.org/script/index.php> : version higher or equal than 6.9.0-0 (for generate screenshot)

Infernal <http://infernal.janelia.org/> : version higher or equal than 1.1.1 (for riboswitch detection)

Matplotlib <http://matplotlib.org/> : version higher or equal than 1.4.3 (for generate all statistics figure)

MEME <http://meme.nbcr.net/meme/> : version higher or equal than 4.10.0 (for promoter detection)

PAGIT and RATT <http://www.sanger.ac.uk/resources/software/pagit/> : version higher or equal than 1.64 (for annotation transfer)

Psortb <http://www.psort.org/psortb/> : version higher or equal than 3.0 (for subcellular localization detection)

Samtools <https://github.com/samtools> : version higher or equal than 1.0 (using htslib 1.0) (for SNP calling, CircRNA_detection)

Bcftools <https://github.com/samtools> : version higher or equal than 1.0 (using htslib 1.0) (for SNP calling)

Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/> : version higher or equal than 0.1.9 (for CircRNA_detection)

TranstermHP <http://transterm.cbcb.umd.edu/> : version higher or equal than 2.09 (for rho-independent terminator prediction)

TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190> : version higher or equal than 1.04 (for TSS and processing site prediction, TSS optimization)

ViennaRNA <http://www.tbi.univie.ac.at/RNA/> : version higher or equal than 2.1.7 (for rho-independent terminator prediction, sRNA detection, sRNA_target detection)

Ps2pdf14 <http://pages.cs.wisc.edu/~ghost/doc/AFPL/6.50/Ps2pdf.htm> (for sRNÁ detection)

BSRD <http://www.bac-srna.org/BSRD/index.jsp> (for sRNA detection)

nr database <ftp://ftp.ncbi.nih.gov/blast/db/FASTA/> (for sRNA detection)

Rfam <http://rfam.xfam.org/> (for sRNA riboswitch)

idmapping_selected.tab from Uniprot <http://www.uniprot.org/downloads> (for Go term detection)

goslim.obo <http://geneontology.org/page/go-slim-and-subset-guide> (for Go term detection) 

go.obo <http://geneontology.org/page/download-ontology> (for Go term detection)

species.vXXXX.txt from STRING <http://string-db.org/newstring_cgi/show_download_page.pl?UserId=ReWbu8uLrfAN&sessionId=_FAQBbatf7RX> (for protein-protein interaction prediction)

## Modules
ANNOgesic is modular pipeline which is constructed by many useful modules.
Users can select some some modules for their specific purpose.

### Setup the analysis folder

The command is:

*python ANNOgesic.py create $ANNOgesic_folder*

The parameter behind *create* is the name of the ANNOgesic folder which will store the results of ANNOgesic.
In this case is ANNOgesic. After create the ANNOgesic folder, user need to put the required files in this folder.
Every module has its pre-required information. We will explain it in description of each module.

### Get input files

The command is like the following:

*python ANNOgesic.py get_input_files *
*-F $FTP_SOURCE*
*-g -f -e -k *
*$ANNOgesic_folder*

*get_input_files* is the subcommand for download required files from NCBI. Therefore, user need to assign the FTP 
address of the reference genome. For example, ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Staphylococcus_aureus_NCTC_8325_uid57795 
Then, user can assign which kinds of files he/she want to download. The final parameter *$ANNOgesic_folder* is the name of 
folder which stores the results of ANNOgesic. If user make it blank, it will assign the current path.

The parameters:

| Full name      | Short flag | Default    | Description                                                        |
|----------------|------------|---------------------------------------------------------------------------------|
| --FTP_path     | -F         | NA         | Path of website where can download the required files.             |
| --ref_fasta    | -f         | False      | Download fasta files of reference.(True/False)                     |
| --ref_gff      | -g         | False      | Download gff files of reference.(True/False)                       |
| --ref_gbk      | -k         | False      | Download genbank files of reference.(True/False)                   |
| --convert_embl | -e         | False      | Convert gbk to embl files of reference.(True/False)                |

### Get target fasta file

Before getting the fasta of target genome, user need to store and assign some files:

1. Fasta file of reference genome.
2. Mutation table which indicate the information of mutations between reference and target genome.

The command is like the following:

*python ANNOgesic.py get_target_fasta -r ref_fasta_folder -o output_file_name:strain_name -m mutation_table $ANNOgesic_folder*

*get_target_fasta* is the subcommand for generating target fasta file from referenece genome files. It is based on the 
mutation table to modify the reference fasta to target fasta. Therefore, the similarity of reference genome and target genome
should be close. The example of format of mutation table is folowing:

| #target_id | reference_id | reference_nt | position | target_nt | impact of correction | locus tag     | gene | Description |
|-----------------------------------------------------------------------------------------------------------------------------|
| HG003      | NC_007795.1  | a            | 333      | c         | replacemen           | SAOUHSC_00002 | dnaA | XXXXXX      |
| HG003      | NC_007795.1  | t            | 543      | -         | deletion             | SAOUHSC_00003 |      | YYYYYY      |
| HG003      | NC_007795.1  | -            | 600      | g         | insertion            | SAOUHSC_00132 |      |             |

If user wants to put the name of column in the top, it needs to start from #. Each column is seperated by tab.
If the mutation type if deletion or insertion, user can put - to represent them. The information of target_id, reference_id,
reference_nt, position, target_nt is required. The others can leave them blank. However, please still use tab to seperate all 
blank columns.

If user has no mutation information between the reference genome and target genome, user can also use *SNP_calling* (one module of
ANNOPgesic) to compute it. Please refer to the description of SNP_calling.

The parameters:

| Full name           | Short flag | Default                                       | Description                                                  |
|-------------------------------------------------------------------------------------------------------------------------------------------------|
| --ref_fasta_folder  | -r         | ANNOgesic/input/reference/target/fasta        | The path of reference fasta folder.                          |
| --mutation_table    | -m         | ANNOgesic/input/mutation_table                | The path of mutation table file.                             |
| --output_format     | -o         | None                                          | Please assign the output filename and which strain should be |
                                                                                     included in it. For example: FILE1:strain1,strain2. FILE1 is
                                                                                     a output fasta file which include the information of strain1
                                                                                     and strain2. You can assign more than one output file. Just 
                                                                                     use space to seperated them.
