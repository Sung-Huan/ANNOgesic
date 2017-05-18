.. _ANNOgesic's subcommands:

ANNOgesic's subcommands
=======================

In general, the subcommands need at least one argument - the analysis
folder. If it is not given, ANNOgesic assumes the current
folder as the analysis folder.

.. _The format of filename:

The format of filename
----------------------
In order to recognize file types as well as relation of genome name and features, 
please use following principle of filename designation:

The genome filenames should be the same as the annotation file designation, i.e.
``NC_007795.fa, NC_007795.gff, NC_007795.ptt, NC_007795.rnt, NC_007705.gbk``.

For all input gff files for ANNOgesic please use following format:
``$GENOME_$FEATURE.gff``. For an example, ``NC_007795_TSS.gff, NC_007795_transcript.gff``.

Possible feature names are:

===============  ===========================
Feature name     meaning
---------------  --------------------------- 
TSS              transcription starting site
processing       processing site
transcript       transcript
sRNA             small RNA
sORF             small open reading frame
term             terminator
5UTR             5'UTR
3UTR             3'UTR
circRNA          circular RNA
riboswitch       riboswitch
RNA_thermometer  RNA_thermometer
CRISPR           CRISPR
===============  ===========================

Please avoid ``|`` in the filename and genome name of Gff3 files or fasta file.

.. _The input format of libraries for running ANNOgesic:

The input format of libraries for running ANNOgesic
---------------------------------------------------

Some ``ANNOgesic`` modules require certain library information. Please use following format:

``$LIBRARY_FILENAME:$LIBRARY_TYPE:$CONDITION:$REPLICATE:$STRAND``

``$LIBRARY_FILENAME`` means the ``.wig`` file.

``$LIBRARY_TYPE`` can be ``tex`` (TEX+) or ``notex`` (TEX-) or ``frag`` (fragmented/conventional).

``$CONDITION`` is the index of conditions. Please use 1, 2, 3, ... to represent different conditions.

``$REPLICATE`` is the index of replicates. Please use a, b, c, ... to represent different replicates.

``$STRAND`` is the strand of wiggle file. Please use + or -.

All libraries should be seperated by colons and listed in one line, like shown in the example above, i.e.
for TEX +/- treated libraries:

::

  TSB_OD_0.2_TEX_reverse.wig:tex:1:a:- \
  TSB_OD_0.5_TEX_reverse.wig:tex:2:a:- \
  TSB_OD_0.2_TEX_forward.wig:tex:1:a:+ \
  TSB_OD_0.5_TEX_forward.wig:tex:2:a:+ \
  TSB_OD_0.2_reverse.wig:notex:1:a:- \
  TSB_OD_0.5_reverse.wig:notex:2:a:- \
  TSB_OD_0.2_forward.wig:notex:1:a:+ \
  TSB_OD_0.5_forward.wig:notex:2:a:+

or for fragmented libraries (RNA-Seq generated after transcript fragmentation):

::

  fragmented_forward.wig:frag:1:a:+ fragmented_reverse.wig:frag:1:a:-

If only conventional RNA-seq data without fragmentation or TEX treated can be provided, 
it can still be assigned to fragmented libraries. However, it may influence the results.

.. _create:

create (create analysis folder)
-------------------------------

``create`` generates the folders for analysis. Once created, please move the required files 
into the corresponding folders.

The folders are following:

**input:** Stores all input files.

	**BAMs:** For ``.bam`` files. ``BAMs_map_related_genomes`` 
	is for the ``.bam`` files which are mapped on closely related genomes of our real query genomes.
	``BAMs_map_query_genomes`` is for the ``.bam`` files which are mapped on our query genomes.

	**databases:** For all databases.

	**manual_TSSs:** If the manual detected transcription starting sites (TSSs) can be provided,
	it can be stored here for running ``TSS_optimization`` or merging 
	the automatic predicted ones and manual detected ones. Please use gff3 format.

	**manual_processing_sites:** It is similar to ``manual_TSS``, but for 
	processing sites.

	**mutation_table:** If the mutation table between the closely ralated genomes and 
	query genomes is provided, please put the file here. Please check 
	the section of :ref:`update_genome_fasta` for the format of 
	mutation table.

	**reads:** For running ``circrna`` with mapping reads by ANNOgesic,
	please put the reads here. ``.bzip2`` and ``.gzip`` as input is accepted.
       
	**references:** For annotation files and fasta files. 
	If they can be downloaded from NCBI, the files can also be obtained via running :ref:`get_input_files`.

	**riboswitch_ID_file:** For storing the file which contains all the Rfam IDs of riboswitches.
	For format details, please check the section of 
	:ref:`riboswitch_thermometer`.

	**RNA_thermometer_ID_file:** For storing the file which contains all the Rfam IDs of RNA thermometer.
	For format details, please check the section of
	:ref:`riboswitch_thermometer`.

	**wigs:** For wiggle files. Based on the methods of RNA-Seq, wiggle files can be stored in  
	``fragment`` (fragmented/conventional libraries) or ``tex_notex`` (TEX +/- treated libraries).

**output:** Stores all output files.

- **Arguments**

::

    usage: annogesic create [-h] [--project_path PROJECT_PATH]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Name/path of the project

.. _get_input_files:

get_input_files (download required files)
-----------------------------------------

``get_input_files`` is the subcommand for downloading required files (fasta, annotation files) from NCBI. 
Therefore, the web address of the reference genome in NCBI needs to be assigned. For an example,
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000013425.1_ASM1342v1
Then, the user can assign the file type for download.


- **Reqired information**

**FTP source:** The IP of NCBI.

- **Arguments**


::

    usage: annogesic get_input_files [-h] --project_path PROJECT_PATH
                                     [--ftp_path FTP_PATH] [--ref_fasta]
                                     [--ref_gff] [--ref_gbk] [--ref_ptt]
                                     [--ref_rnt] [--convert_embl]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --ftp_path FTP_PATH, -F FTP_PATH
                            Path of folder on the NCBI FTP server where the
                            required files are located.
      --ref_fasta, -f       Download fasta files of the reference. Default is
                            False.
      --ref_gff, -g         Download gff files of the reference. Default is False.
      --ref_gbk, -k         Download genbank files of the reference. Default is
                            False.
    
    additional arguments:
      --ref_ptt, -p         Download ptt files of the reference. Default is False.
      --ref_rnt, -r         Download rnt files of the reference. Default is False.
      --convert_embl, -e    Convert gbk to embl files of the reference. Default is
                            False.

- **Output files**

Output files will be stored in ``$ANNOgesic_folder/input/reference``

Output folder names are following:

**fasta:** Fasta files.

**annotation:** Annotation files.

.. _update_genome_fasta:

update_genome_fasta (update reference genome fasta file)
--------------------------------------------------------

If fasta files of the query genomes do not exist, ``update_genome_fasta`` can 
update fasta files from the closely related genomes of the real query genomes to our real query ones 
via searching the mutations. 
Therefore, a table of mutation information is required. For the format of the table, please check 
`mutation table <https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/tutorial_data/mutation.csv>`_.
Titles of the columns are presented on the top and they need to start with ``#``. 
Each column is separated by ``tab``. If the mutation type is deletion or insertion, 
the user can type ``-`` to represent them. The information of ``Target_ids`` 
(the query genome names), ``Reference_ids``, (the names of closely related genomes) 
``Reference_nts`` (the nucleotides of the closely related genomes), ``Positions``, ``Target_nts`` 
(the nucleotides of the query genomes) are required. The other columns can be blank. 
Please use ``tab`` to separate all columns including blank ones.

If no mutation information is provided, ``snp`` can be used for detecting mutations. 
(one module of ``ANNOgesic``). Please check the section of :ref:`snp`.

- **Required files**

**Fasta files of reference genome**

**Mutation table:** Contains the information of mutations between related and query genomes.
For an example, please check `mutation table <https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/tutorial_data/mutation.csv>`_.

- **Arguments**

::

    usage: annogesic update_genome_fasta [-h] --project_path PROJECT_PATH
                                         --related_fasta_files RELATED_FASTA_FILES
                                         [RELATED_FASTA_FILES ...]
                                         --mutation_table MUTATION_TABLE
                                         [--combine_to_one_fasta]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --related_fasta_files RELATED_FASTA_FILES [RELATED_FASTA_FILES ...], -c RELATED_FASTA_FILES [RELATED_FASTA_FILES ...]
                            Path of the genome fasta files of the closely related
                            species.
      --mutation_table MUTATION_TABLE, -m MUTATION_TABLE
                            Path of the mutation table which stores the mutation
                            information between the query genome and genome of the
                            closely related species. For an example check
                            https://github.com/Sung-
                            Huan/ANNOgesic/blob/master/tutorial_data/mutation.csv
      --combine_to_one_fasta, -cm
                            Combinine all updated sequences in --mutation_table to
                            one fasta file. Default is False.

- **Output files**

**Fasta files of updated genome**: The updated fasta files are stored in ``$ANNOgesic_folder/output/updated_references/fasta_files``.

.. _annotation_transfer:

annotation_transfer (annotation transfer)
-----------------------------------------

``annotation transfer`` is the subcommand for transferring the annotation from the closely related genomes 
to the real query genomes. To achieve this, `RATT <http://www.sanger.ac.uk/resources/software/pagit/>`_ 
is integrated in ANNOgesic. The higher similarity between closely related genomes and query genomes are, 
the more precise the performance is. Before running ``annotation transfer``, 
please run ``source $PAGIT_HOME/sourceme.pagit`` first. it will modify the path for executing RATT. 
If you use Dockerfile to execute ANNOgesic, the path modification can be skipped.

- **Required tools**

`RATT <http://www.sanger.ac.uk/resources/software/pagit/>`_.

- **Required files**

**Annotation files of the closely related genomes**: Genbank/embl files of the closely related genomes.

**Fasta files of the closely related genomes**

**Fasta files of the updated genomes**

- **Arguments**

::

    usage: annogesic annotation_transfer [-h] --project_path PROJECT_PATH
                                         --compare_pair COMPARE_PAIR
                                         [COMPARE_PAIR ...]
                                         [--related_embl_files RELATED_EMBL_FILES [RELATED_EMBL_FILES ...]]
                                         [--related_gbk_files RELATED_GBK_FILES [RELATED_GBK_FILES ...]]
                                         --related_fasta_files RELATED_FASTA_FILES
                                         [RELATED_FASTA_FILES ...]
                                         --updated_fasta_files UPDATED_FASTA_FILES
                                         [UPDATED_FASTA_FILES ...]
                                         [--ratt_path RATT_PATH] --element ELEMENT
                                         [--transfer_type TRANSFER_TYPE]
                                         [--convert_to_gff_rnt_ptt]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --compare_pair COMPARE_PAIR [COMPARE_PAIR ...], -p COMPARE_PAIR [COMPARE_PAIR ...]
                            Please assign the sequence identifier of genome pairs,
                            e.g. NC_007795:NEW_NC_007795. The related genome is
                            NC_007795 and the target genome is NEW_NC_007795. The
                            assigned names are the headers of the fasta file
                            (start with ">"), not the filename of fasta file. If
                            multiple sequences need to be assigned, please use
                            spaces to separate them.
      --related_embl_files RELATED_EMBL_FILES [RELATED_EMBL_FILES ...], -ce RELATED_EMBL_FILES [RELATED_EMBL_FILES ...]
                            The paths of the embl files of the related species.
      --related_gbk_files RELATED_GBK_FILES [RELATED_GBK_FILES ...], -cg RELATED_GBK_FILES [RELATED_GBK_FILES ...]
                            The paths of the genbank files of the related species.
                            The genbank can be ended by .gbk, .gbff or .gb
      --related_fasta_files RELATED_FASTA_FILES [RELATED_FASTA_FILES ...], -cf RELATED_FASTA_FILES [RELATED_FASTA_FILES ...]
                            The paths of the fasta files of the related species.
      --updated_fasta_files UPDATED_FASTA_FILES [UPDATED_FASTA_FILES ...], -uf UPDATED_FASTA_FILES [UPDATED_FASTA_FILES ...]
                            The paths of updated fasta files.
    
    additional arguments:
      --ratt_path RATT_PATH
                            Path of the start.ratt.sh file of RATT folder. Default
                            is start.ratt.sh.
      --element ELEMENT, -e ELEMENT
                            --element will become the prefix of all output file.
      --transfer_type TRANSFER_TYPE, -t TRANSFER_TYPE
                            The transfer type for running RATT. (For the details,
                            please refer to the manual of RATT.) Default is
                            Strain.
      --convert_to_gff_rnt_ptt, -g
                            Convert the annotation to gff, rnt and ptt. Default is
                            False.

- **Output files**

Output files from `RATT <http://www.sanger.ac.uk/resources/software/pagit/>`_
will be stored in ``$ANNOgesic_folder/output/annotation_transfer``.

**Annotation files** (``.gff``, ``.ptt``, ``.rnt``) will be stored in ``$ANNOgesic_folder/output/updated_references/annotations``.

.. _snp:

snp (SNP calling)
-----------------

``snp`` can analyze the alignment files and fasta files to detect mutations by running 
`Samtools <https://github.com/samtools>`_ and `Bcftools <https://github.com/samtools>`_. 
There are multiple programs which can be applied to detect mutations 
(with BAQ, without BAQ and extend BAQ) and there are multiple flag options to set filters
(QUAL, DP, DP4, etc.). Moreover, ``snp`` can also be used for generating the fasta files of 
query genomes if it is necessary.

- **Required tools**

`Samtools <https://github.com/samtools>`_.

`Bcftools <https://github.com/samtools>`_.

- **Required files**

**BAM files:** BAM files from fragmented/conventional libraries or TEX +/- treated libraries both can be accepted.
For assigning the files, please follow the format -- ``$SET_NAME:$BAMFILE1,$BAMFILE2,...``. 
For an example, the user has four bam files of one genome. Then the input will be 
``set1:sample1.bam,sample2.bam,sample3.bam,sample4.bam``.

**Fasta files of the closely related genomes** or **Fasta files of the query genomes**

- **Arguments**

::

    usage: annogesic snp [-h] --project_path PROJECT_PATH --bam_type
                         {related_genome,query_genome} --program
                         {with_BAQ,without_BAQ,extend_BAQ}
                         [{with_BAQ,without_BAQ,extend_BAQ} ...] --fasta_files
                         FASTA_FILES [FASTA_FILES ...] --bam_files BAM_FILES
                         [BAM_FILES ...] [--samtools_path SAMTOOLS_PATH]
                         [--bcftools_path BCFTOOLS_PATH] [--quality QUALITY]
                         [--read_depth_range READ_DEPTH_RANGE]
                         [--ploidy {haploid,diploid}] [--rg_tag] [--caller {c,m}]
                         [--dp4_cutoff DP4_CUTOFF]
                         [--indel_fraction INDEL_FRACTION]
                         [--filter_tag_info FILTER_TAG_INFO [FILTER_TAG_INFO ...]]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --bam_type {related_genome,query_genome}, -t {related_genome,query_genome}
                            If the BAM files are produced by mapping to a related
                            genome, please assign "related_genome". the mutations
                            between the related genome and the query genome can be
                            detected for generating sequence of the query genome.
                            If the BAM files are produced by mapping to the query
                            genome, please assign "query_genome". The mutations of
                            query genome can be detected.
      --program {with_BAQ,without_BAQ,extend_BAQ} [{with_BAQ,without_BAQ,extend_BAQ} ...], -p {with_BAQ,without_BAQ,extend_BAQ} [{with_BAQ,without_BAQ,extend_BAQ} ...]
                            The program for detecting SNPs: "with_BAQ",
                            "without_BAQ", "extend_BAQ". Multi-programs can be
                            executed at the same time (separated by spaces). For
                            example, with_BAQ without_BAQ extend_BAQ.
      --fasta_files FASTA_FILES [FASTA_FILES ...], -f FASTA_FILES [FASTA_FILES ...]
                            Paths of the genome fasta files.
      --bam_files BAM_FILES [BAM_FILES ...], -b BAM_FILES [BAM_FILES ...]
                            Path of input BAM files. Required format:
                            $SET_NAME:$BAM1,$BAM2,... . If multiple sets need to
                            be assigned, please use spaces to separate them.
    
    additional arguments:
      --samtools_path SAMTOOLS_PATH
                            Path of samtools.
      --bcftools_path BCFTOOLS_PATH
                            Path of bcftools.
      --quality QUALITY, -q QUALITY
                            The minimum quality of a mutation. Default is 40.
      --read_depth_range READ_DEPTH_RANGE, -d READ_DEPTH_RANGE
                            The range of read depth for a mutation. The format is
                            $MIN,$MAX. It can be assigned by different types: 1.
                            real number ("r"), 2. times of the number of bam files
                            (counted from --bam_files) ("n") or 3. times of the
                            average read depth ("a"). For example, n_10,a_3 is
                            assigned. If the average read depth is 70 and 4 bam
                            files are provided, n_10 will be 40 and a_3 will be
                            140 (average read depth * 3). Based on the same
                            example, if this value is r_10,a_3, the minimum read
                            depth will become exact 10 reads. If "none" is
                            assigned, read depth will not be considered. Default
                            is n_10,none.
      --ploidy {haploid,diploid}, -pl {haploid,diploid}
                            The query bacteria is haploid or diploid. Default is
                            haploid.
      --rg_tag, -R          For one BAM file which includes multiple samples
                            (opposite of --ignore-RG in samtools). Default is
                            False.
      --caller {c,m}, -c {c,m}
                            The types of caller - consensus-caller or
                            multiallelic-caller. For details, please check
                            documentation of bcftools. "c" represents consensus-
                            caller. "m" represents multiallelic-caller. Default is
                            m.
      --dp4_cutoff DP4_CUTOFF, -D DP4_CUTOFF
                            The cutoff of DP4. The format is
                            $MIN_SNP_READS_NUMBER:$MIN_SNP_READS_RATIO.
                            $MIN_SNP_READS_NUMBER can be assigned by three types:
                            1. real number ("r"), 2. times of the number of bam
                            files (counted from --bam_files) ("n") or 3. times of
                            average read depth ("a"). For example, n_10,0.8. If
                            the BAM files are 4, it means the minimum mapped reads
                            of a SNP is 40 (4 * 10), and the minimum ratio of
                            mapped read of a SNP (mapped reads of a SNP / total
                            reads) is 0.8. Default is n_10,0.8.
      --indel_fraction INDEL_FRACTION, -if INDEL_FRACTION
                            The minimum IDV and IMF which supports for insertion
                            of deletion. The minimum IDV can be assigned by
                            different types: 1. real number ("r"), 2. times of the
                            number of bam files (assigned by --bam_files) ("n") or
                            3. times of the average read depth ("a"). The input
                            format is $MIN_IDF:$MIN_IMF. For example, The value is
                            n_10,0.8 and 4 BAM files are assigned. The minimum IDV
                            is 40, and the minimum IMF is 0.8. Default is
                            n_10,0.8.
      --filter_tag_info FILTER_TAG_INFO [FILTER_TAG_INFO ...], -ft FILTER_TAG_INFO [FILTER_TAG_INFO ...]
                            For using more filters to improve the detection.
                            Please assign 1. the name of tag, 2. bigger ("b") or
                            smaller ("s") and 3. the value of the filter. For
                            example, "RPB_b0.1,MQ0F_s0" means that RPB should be
                            bigger than 0.1 and MQ0F should be smaller than 0.
                            Default is RPB_b0.1,MQSB_b0.1,MQB_b0.1,BQB_b0.1.

- **Output files**

If ``bam_type`` is ``related_genome``, 
the results will be stored in ``$ANNOgesic/output/SNP_calling/compare_related_and_query_references``. 
If ``bam_type`` is ``query_genome``, the results will be stored in ``$ANNOgesic/output/SNP_calling/mutations_of_query_genomes``.

The output folders and results are following:

**SNP_raw_output:** Stores output tables which be only considered read depth and QUAL.

	**VCF Table (only consider read depth and QUAL):** Filename is ``$GENOME_$PROGRAM_$SET.vcf``.

**SNP_table:** Stores two types of output tables

        **VCF Table (consider all filters):** Filename is ``$GENOME_$PROGRAM_$SET_best.vcf``.

        **Index of fasta files:**: Filename is ``$GENOME_$PROGRAM_$SET_seq_reference.csv``.
        The meaning of this file is like following example:

::

  Staphylococcus_aureus_HG003     1632629 .       AaA     AA      57      .
  Staphylococcus_aureus_HG003     1632630 .       aA      a       57      .
  Staphylococcus_aureus_HG003     1499572 .       T       TT,TTTTT        43.8525 .

The example contains "position conflict" and "mutation conflict".
As a result, the conflicts will affect the other mutation's positions.
Therefore, it will generate four different fasta files. The first two lines are "position conflict", and 
the last line is "mutation conflict".
``$GENOME_$PROGRAM_$SET_seq_reference.csv`` is the index for these four fasta files.

::

   1       1632629 1       1499572:TT      Staphylococcus_aureus_HG003
   1       1632629 2       1499572:TTTTT   Staphylococcus_aureus_HG003
   2       1632630 1       1499572:TT      Staphylococcus_aureus_HG003
   2       1632630 2       1499572:TTTTT   Staphylococcus_aureus_HG003

The first column is the index of the "position conflict". 
The second column is the selected position.
The third one is the index of the "mutations conflict". 
The fourth one is the selected position and nucleotides. 
The last column is the genome name.

**Potential fasta files**: Filename is ``$FASTANAME_$SET_$STRIANNAME_$INDEXofPOSITIONCONNFLICT_$INDEXofMUTATIONCONFLICT.fa``, 
and it is stored in ``$ANNOgesic/output/SNP_calling/$BAM_TYPE/seqs``.
Based on the example in **Index of fasta files**, ``Staphylococcus_aureus_HG003_set1_Staphylococcus_aureus_HG003_1_1.fa``
will be generated based on the first line of ``$GENOME_$PROGRAM_seq_reference.csv``.
``Staphylococcus_aureus_HG003_set1_Staphylococcus_aureus_HG003_1_2.fa`` and will be generated based on the second line of 
``$GENOME_$PROGRAM_seq_reference.csv`` and so forth.

**statistics**: Stores the statistic files and figures, ex: the distribution of SNPs based on QUAL.

.. _tss_ps:

tss_ps (TSS and processing site prediction)
-------------------------------------------

``tss_ps`` can generate the TSS and processing sites via running  
`TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_. Since the parameters can affect the 
results strongly, ``optimize_tss_ps`` can obtain the optimized parameters of 
`TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_. Please check the section 
:ref:`optimize_tss_ps` for details.

- **Required tools**

`TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_.

- **Required files**

**Wiggle files of TEX +/-:** Please check the section :ref:`The input format of libraries for running ANNOgesic` for assigning correct format.

**Fasta files of the reference genomes**

**GFF files of the reference genomes**

- **Optional input files**

**Gff files of the manual detected TSSs:** If gff file of the manual detected TSSs can be provided, ``tss_ps`` can merge the manual detected TSSs
and TSSpredator predicted ones.

**Gff files of transcripts:** If comparing TSSs with transcripts is required, gff files of the transcripts need to be assigned.
For the transcripts, please check the section :ref:`transcript`.

- **Arguments**

::

    usage: annogesic tss_ps [-h] --project_path PROJECT_PATH [--program {TSS,PS}]
                            --fasta_files FASTA_FILES [FASTA_FILES ...]
                            --annotation_files ANNOTATION_FILES
                            [ANNOTATION_FILES ...] --tex_notex_libs TEX_NOTEX_LIBS
                            [TEX_NOTEX_LIBS ...]
                            [--replicate_tex REPLICATE_TEX [REPLICATE_TEX ...]]
                            --condition_names CONDITION_NAMES
                            [CONDITION_NAMES ...]
                            [--tsspredator_path TSSPREDATOR_PATH]
                            [--height HEIGHT]
                            [--height_reduction HEIGHT_REDUCTION]
                            [--factor FACTOR]
                            [--factor_reduction FACTOR_REDUCTION]
                            [--enrichment_factor ENRICHMENT_FACTOR]
                            [--processing_factor PROCESSING_FACTOR]
                            [--base_height BASE_HEIGHT] [--utr_length UTR_LENGTH]
                            [--fuzzy FUZZY] [--cluster CLUSTER]
                            [--manual_files MANUAL_FILES [MANUAL_FILES ...]]
                            [--curated_sequence_length CURATED_SEQUENCE_LENGTH [CURATED_SEQUENCE_LENGTH ...]]
                            [--validate_gene]
                            [--compare_transcript_files COMPARE_TRANSCRIPT_FILES [COMPARE_TRANSCRIPT_FILES ...]]
                            [--re_check_orphan] [--remove_overlap_feature]
                            [--compare_overlap_gff COMPARE_OVERLAP_GFF [COMPARE_OVERLAP_GFF ...]]
                            [--remove_low_expression REMOVE_LOW_EXPRESSION]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --program {TSS,PS}, -p {TSS,PS}
                            The feature to predict. Please assign "TSS" or "PS".
                            Default is "TSS".
      --fasta_files FASTA_FILES [FASTA_FILES ...], -f FASTA_FILES [FASTA_FILES ...]
                            Paths of the query genome fasta files.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the query genome gff files.
      --tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...], -tl TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]
                            TEX+/- wig files for TSSpredator. The format is:
                            wig_file_path:TEX+/-(tex or notex):condition_id(intege
                            r):replicate_id(alphabet):strand(+ or -). If multiple
                            wig files need to be assigned, please use spaces to
                            separate the wig files. For example,
                            my_lib_tex_forward.wig:tex:1:a:+
                            my_lib_tex_reverse.wig:tex:1:a:-.
      --replicate_tex REPLICATE_TEX [REPLICATE_TEX ...], -rt REPLICATE_TEX [REPLICATE_TEX ...]
                            This value is the minimal number of replicates that a
                            TSS has to be detected in. The format is
                            $NUMBERofCONDITION_$NUMBERofREPLICATE. If different
                            --replicate_tex values need to be assigned to
                            different conditions, please use spaces to separate
                            them. For example, 1_2 2_2 3_3 means that
                            --replicate_tex is 2 in number 1 and number 2
                            conditions. In number 3 condition, --replicate_tex is
                            3. For assigning the same --replicate_tex to all
                            conditions, just use like all_1 (--replicate_tex is 1
                            in all conditions). Default is all_1.
      --condition_names CONDITION_NAMES [CONDITION_NAMES ...], -cn CONDITION_NAMES [CONDITION_NAMES ...]
                            The output prefix of all conditions. If multiple
                            conditions need to be assigned, please use spaces to
                            separate them. For an example, prefix_condition1
                            prefix_condition2.
    
    additional arguments:
      --tsspredator_path TSSPREDATOR_PATH
                            Path of TSSpredator. Default is
                            /usr/local/bin/TSSpredator.jar
      --height HEIGHT, -he HEIGHT
                            This value relates to the minimal number of read
                            starts at a certain genomic position to be considered
                            as a TSS candidate. Default is 0.3.
      --height_reduction HEIGHT_REDUCTION, -rh HEIGHT_REDUCTION
                            When comparing different genomes/conditions and the
                            step height threshold is reached in at least one
                            genome/condition, the threshold is reduced for the
                            other genomes/conditions by the value set here. This
                            value must be smaller than the step height threshold.
                            Default is 0.2.
      --factor FACTOR, -fa FACTOR
                            The minimal factor by which the TSS height has to
                            exceed the local expression background. Default is
                            2.0.
      --factor_reduction FACTOR_REDUCTION, -rf FACTOR_REDUCTION
                            When comparing different genomes/conditions and the
                            step factor threshold is reached in at least one
                            genome/condition, the threshold is reduced for the
                            other genomes/conditions by the value set here. This
                            value must be smaller than the step factor threshold.
                            Default is 0.5.
      --enrichment_factor ENRICHMENT_FACTOR, -ef ENRICHMENT_FACTOR
                            The minimal enrichment factor. Default is 2.0.
      --processing_factor PROCESSING_FACTOR, -pf PROCESSING_FACTOR
                            The minimal processing factor. If the value for the
                            untreated library is higher than the treated library
                            the positionsis considered as a processing site and
                            not annotated as detected. Default is 1.5.
      --base_height BASE_HEIGHT, -bh BASE_HEIGHT
                            The minimal number of reads should be mapped on TSS.
                            Default is 0.0.
      --utr_length UTR_LENGTH, -u UTR_LENGTH
                            The length of UTRs. Default is 300.
      --fuzzy FUZZY, -fu FUZZY
                            If --compare_transcript_files is provided, please
                            assign the fuzzy for comparing TSS and transcript.
                            Default is 5.
      --cluster CLUSTER, -c CLUSTER
                            This value defines the maximal distance (nucleotides)
                            between TSS candidates have to be clustered together.
                            If the distance between these multiple TSSs is smaller
                            or equal to this value, only one of them will be
                            printed out. Default is 2.
      --manual_files MANUAL_FILES [MANUAL_FILES ...], -m MANUAL_FILES [MANUAL_FILES ...]
                            If gff files of the manual checked TSS are provided,
                            this function will merge manual checked ones and
                            TSSpredator predicted ones. Please assign the path of
                            manual-checked TSS gff files.
      --curated_sequence_length CURATED_SEQUENCE_LENGTH [CURATED_SEQUENCE_LENGTH ...], -le CURATED_SEQUENCE_LENGTH [CURATED_SEQUENCE_LENGTH ...]
                            The length of the sequence used for the manual set of
                            TSS/PS. This value is required to calculate the
                            accurracy. If the whole genome was used write "all".
                            Otherwise use the name of the reference sequence in
                            the folowing format: $GENOME:SLENGTH. Multiple entries
                            are accepted. For an example, test.gff contains two
                            sequences s1 and s2. For s1 100 kb were checked while
                            for s2 the whole sequence was curated. The value of
                            this argument would be s1:100000 s2:all. Per default
                            all the full length of all sequences will be used.
      --validate_gene, -v   Using TSS candidates to validate genes in annotation
                            file. it will be store in statistics folder. Default
                            is False.
      --compare_transcript_files COMPARE_TRANSCRIPT_FILES [COMPARE_TRANSCRIPT_FILES ...], -ta COMPARE_TRANSCRIPT_FILES [COMPARE_TRANSCRIPT_FILES ...]
                            If the paths of transcript gff files are provided,
                            this function will compare TSS and transcript to
                            obtain the overlap information. Default is False.
      --re_check_orphan, -ro
                            If there is no information of gene or locus_tag in
                            genome annotation gff file, all TSSs will be assigned
                            to orphan TSSs by TSSpredator. The function can
                            compare TSSs with CDSs to classify the TSS correctly.
                            Default is False.
      --remove_overlap_feature, -of
                            If a processing site and a TSS are overlaping, keep
                            "TSS", The predicted feature (based on --program) will
                            be removed. Default is False.
      --compare_overlap_gff COMPARE_OVERLAP_GFF [COMPARE_OVERLAP_GFF ...], -rg COMPARE_OVERLAP_GFF [COMPARE_OVERLAP_GFF ...]
                            If --overlap_feature is "TSS" or "PS",
                            --reference_gff_files need to be assigned. For TSS
                            prediction, please assign the path of processing site.
                            For processing site prediction, please assign the path
                            of TSS. Don't use this flag if --overlap_feature is
                            "both".
      --remove_low_expression REMOVE_LOW_EXPRESSION, -rl REMOVE_LOW_EXPRESSION
                            For removing the low expressed TSS by comparing the
                            manual detected TSSs and predicted ones. Please assign
                            the manual-checked TSS in gff format.

- **Output files**

The results of TSS are stored in ``$ANNOgesic/output/TSSs``, and the results of processing site 
are stored in ``$ANNOgesic/output/processing_sites``.

The output folders are following:

**MasterTables:** MasterTable from `TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_.

**statistics:** Statistic files.

	**Venn Figures of TSS types:** Filename is ``TSS_venn_$GENOME.png``.

	**TSS types with corresponding amounts:** Table is ``stat_TSS_class_$GENOME.csv``, and Figure is ``TSS_class_$GENOME.png``.

	**Conditions with corresponding amounts:** ``stat_TSS_libs_$GENOME.csv`` stores all combination of conditions with corresponding amounts.
	``TSSstatistics.tsv`` stores the number of TSS which can be detected or missing in each condition.

	**Comparing TSSs with other features:** ``stat_compare_TSS_transcript_$GENOME.csv`` is for comparing TSSs with transcripts.
	``stat_gene_vali_$GENOME.csv`` is for comparing TSS with genome annotations.

	**Comparing manual detected TSSs and predicted TSSs:** In ``stat_compare_TSSpredator_manual_$GENOME.csv``, the accuracy of TSS prediction can be found.

**configs**: Configuration files for running TSSpredator.

**gffs**: Output gff files of TSSs. Some useful information can be found in the tags of the attributes within the TSS gff files. 
Based on this information, we can know the details of the specific TSS. The tags are as following:

	**method:** Stores the information that this TSS is detected by manual detection or `TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_.
	
	**type:** TSS type of this TSS. It could be Primary, Secondary, Internal, Antisense or Orphan.
	
	**utr_length:** UTR length of this TSS.
	
	**associated_gene**: Which genes are associated with this TSS.
	
	**Parent:** Presents the parent transcripts of this TSS, if the user has compared TSS with the transcript.
	
	**libs:** Shows in which libraries the TSS can be detected.

.. _transcript:

transcript (transcript detection)
---------------------------------

``transcript`` can detect transcripts based on the coverage. Most of the transcript assembly tools are
focus on eukaryotic transcript. Due to this, we constructed a subcommand which is based on the nucleotide coverage data, 
given gene annotations and several parameters that can be set by the user.

- **Required files**

**Wiggle files of fragmented/conventional libraries or TEX+/- treated libraries:** For importing the information about libraries, please check the section 
:ref:`The input format of libraries for running ANNOgesic`.

- **Optional input files**

**TSS gff files:** If the user wants to compare transcripts with TSSs, TSS gff files are required.

**Genome anntation gff files:** If the user wants to compare transcripts with genome annotations or modify transcript by genome annotations, 
genome annotation gff files are required. There are four options for modification of transcripts:

	**merge_overlap:** If multiple transcripts overlap the same gene, they will be merged as one complete transcript.

	**extend_3end:** If the transcript starts at the upstream of the gene and ends within the gene, 
	the end point of the transcript will be extended to the end point of gene.

	**extend_5end:** If the transcript starts within the gene and ends at the downstream of gene, 
	the starting point of the transcript will be extended to the starting point of the gene.

        **within_extend_ends:** If the transcript is within the gene, the two ends of the transcript will be 
	extended to the two ends of gene.

	**none:** Transcripts will not be modified by the genome annotations

- **Arguments**

::

    usage: annogesic transcript [-h] --project_path PROJECT_PATH
                                [--annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...]]
                                [--modify_transcript {merge_overlap,extend_3end,extend_5end,within_extend_ends,none} [{merge_overlap,extend_3end,extend_5end,within_extend_ends,none} ...]]
                                [--tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]]
                                [--frag_libs FRAG_LIBS [FRAG_LIBS ...]]
                                [--replicate_tex REPLICATE_TEX [REPLICATE_TEX ...]]
                                [--replicate_frag REPLICATE_FRAG [REPLICATE_FRAG ...]]
                                [--tex_notex {1,2}] [--length LENGTH]
                                [--height HEIGHT] [--width WIDTH]
                                [--tolerance TOLERANCE]
                                [--tolerance_coverage TOLERANCE_COVERAGE]
                                [--tss_files TSS_FILES [TSS_FILES ...]]
                                [--compare_feature_genome COMPARE_FEATURE_GENOME [COMPARE_FEATURE_GENOME ...]]
                                [--tss_fuzzy TSS_FUZZY] [--table_best]
                                [--terminator_files TERMINATOR_FILES [TERMINATOR_FILES ...]]
                                [--terminator_fuzzy TERMINATOR_FUZZY]
                                [--max_length_distribution MAX_LENGTH_DISTRIBUTION]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            If paths of the genome annotation gff files.
      --modify_transcript {merge_overlap,extend_3end,extend_5end,within_extend_ends,none} [{merge_overlap,extend_3end,extend_5end,within_extend_ends,none} ...], -mt {merge_overlap,extend_3end,extend_5end,within_extend_ends,none} [{merge_overlap,extend_3end,extend_5end,within_extend_ends,none} ...]
                            If --annotation_files is provided, the post-
                            modification of transcript based on genome annotations
                            can be assigned. There are five opetions. 1.
                            "merge_overlap": if multiple transcripts overlap with
                            the same gene, they will be merged as one complete
                            transcript. 2. "extend_3end": if a transcript starts
                            at the upstream of a gene and ends within the gene,
                            the end point of the transcript will be extended to
                            the end point of the gene. 3. "extend_5end": if a
                            transcript starts within a gene and ends at the
                            downstream of gene, the starting point of the
                            transcript will be extended to the starting point of
                            the gene. 4. "within_extend_ends": if a transcript is
                            within a gene, the two ends of the transcript will be
                            extended to the two ends of gene. 5. "none": the
                            transcript will not be modified by the genome
                            annotations. For using mutliple modifications, please
                            separated them by spaces. Default is merge_overlapped.
      --tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...], -tl TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]
                            TEX+/- wig files. The format is:
                            wig_file_path:TEX+/-(tex or notex):condition_id(intege
                            r):replicate_id(alphabet):strand(+ or -). If multiple
                            wig files need to be assigned, please use spaces to
                            separate the wig files. For example,
                            my_lib_tex_forward.wig:tex:1:a:+
                            my_lib_tex_reverse.wig:tex:1:a:-.
      --frag_libs FRAG_LIBS [FRAG_LIBS ...], -fl FRAG_LIBS [FRAG_LIBS ...]
                            Wig files of RNA-Seq with transcript fragmented. The
                            format is: wig_file_path:frag:condition_id(integer):re
                            plicate_id(alphabet):strand(+ or -). If multiple wig
                            files need to be assigned, please use spaces to
                            separate the wig files. For example,
                            my_lib_frag_forward.wig:frag:1:a:+
                            my_lib_frag_reverse.wig:frag:1:a:-.
      --replicate_tex REPLICATE_TEX [REPLICATE_TEX ...], -rt REPLICATE_TEX [REPLICATE_TEX ...]
                            This value is the minimal number of replicates that a
                            TSS has to be detected in. The format is
                            $NUMBERofCONDITION_$NUMBERofREPLICATE. If different
                            --replicate_tex values need to be assigned to
                            different conditions, please use spaces to separate
                            them. For example, 1_2 2_2 3_3 means that
                            --replicate_tex is 2 in number 1 and number 2
                            conditions. In number 3 condition, --replicate_tex is
                            3. For assigning the same --replicate_tex to all
                            conditions, just use like all_1 (--replicate_tex is 1
                            in all conditions). Default is all_1.
      --replicate_frag REPLICATE_FRAG [REPLICATE_FRAG ...], -rf REPLICATE_FRAG [REPLICATE_FRAG ...]
                            It is similar to --replicates_tex. This value is for
                            fragmented (or conventional) libraries.
      --tex_notex {1,2}, -te {1,2}
                            The value is for TEX+/- libraries to decide the
                            transcript should be detected in both (TEX+ and TEX-)
                            or can be detected in only one library (TEX+ or TEX-).
                            Please assign 1 or 2. Default is 1.
    
    additional arguments:
      --length LENGTH, -l LENGTH
                            The minimum length of the transcript after modifying
                            base on genome annotation. Default is 20.
      --height HEIGHT, -he HEIGHT
                            The minimum coverage of the transcript. If --tex_notex
                            is 2, the average coverage of TEX+ and TEX- libraries
                            should be higher than this value. The default is 10.
      --width WIDTH, -w WIDTH
                            The minimum length of the transcript without modifying
                            by genome annotation. The default is 20.
      --tolerance TOLERANCE, -t TOLERANCE
                            The number of nucleotides which coveraes can drop
                            below the --height in a transcript. The default is 5.
      --tolerance_coverage TOLERANCE_COVERAGE, -tc TOLERANCE_COVERAGE
                            The minimum coverage of the nucleotides which match
                            the situation of --tolerance, Default is 0.
      --tss_files TSS_FILES [TSS_FILES ...], -ct TSS_FILES [TSS_FILES ...]
                            Paths of TSS files for comparing transcripts and TSSs.
      --compare_feature_genome COMPARE_FEATURE_GENOME [COMPARE_FEATURE_GENOME ...], -cf COMPARE_FEATURE_GENOME [COMPARE_FEATURE_GENOME ...]
                            If --compare_genome_annotation is provided, please
                            assign the feature for comparing. Multiple features
                            can be separated by spaces. Default is None.
      --tss_fuzzy TSS_FUZZY, -tf TSS_FUZZY
                            The fuzzy value for comparing TSS with transcript.
                            Default is 5.
      --table_best, -tb     The output table only includes the information of the
                            highest expressed library. Default is False.
      --terminator_files TERMINATOR_FILES [TERMINATOR_FILES ...], -e TERMINATOR_FILES [TERMINATOR_FILES ...]
                            Paths of terminator gff files for comparing
                            transcripts and terminators. Default is None.
      --terminator_fuzzy TERMINATOR_FUZZY, -ef TERMINATOR_FUZZY
                            If --terminator_files is assigned, please assign the
                            fuzzy value. Default is 30.
      --max_length_distribution MAX_LENGTH_DISTRIBUTION, -mb MAX_LENGTH_DISTRIBUTION
                            For generating the figure of distribution of
                            transcript length, please assign the maximum length.
                            Default is 2000.

- **Output files**

Output files are stored in ``$ANNOgesic/output/transcripts``.

The generated output folders are as following:

**tables:** Table of transcript with more details. The meanings of the columns in the table are following:

	**Genome:** Genome name.

	**Name:** Transcript name in the gff file.

	**Start:** Starting point of this transcript.

	**End:** End point of this transcript.

	**Strand:** Strand of this transcript.

	**Detect_lib_type:** This transcript can be detected in fragmented/conventional or TEX+/- libraries.

	**Associated_gene:** Which genes are associated with this transcript.

	**Associated_tss:** Which TSSs are located on this transcript.

	**Associated_term:** Which terminators are associated with this transcript.

	**Coverage_details:** Stores the average coverage information of all libraries about this transcript.

**statistics:** Stores statistic files.

	**Comparing transcript with other features:** ``stat_compare_transcript_genome_$GENOMENAME.csv`` is 
	for comparing transcript with genome annotation, ``stat_compare_transcript_TSS_$GENOMENAME.csv`` is for comparing 
	transcript with TSS, and ``stat_compare_transcript_terminator_$GENOMENAME.csv`` is for comparing
        transcript with terminator.

	**Figure of the distribution of transcript length:** ``$GENOME_length_all.png`` is for analyzing of all transcript length. 
	``$GENOME_length_less_$LENGTH.png`` is for the analyzing of the assigned length.

**gffs:** Stores gff files of transcripts. Some useful information can be found in the tags of the attributes within the transcript gff file.
Based on this information, we can know the details of the specific transcript. The tags are as following:

	**compare_$FEATURE:** State of overlap between transcripts and features
	(If ``--compare_feature_genome`` and ``--annotation_files`` are assigned). The value may be "cover", "right_shift", "left_shift", "within" or "no_related".

	**associated_tss:** Shows which TSSs are located on this transcript (If ``--tss_files`` is assigned).

	**associated_term:** Shows which terminators are located on this transcript (If ``--terminator_files`` is assigned).

	**associated_$FEATURE:** Shows that the features are located on this transcript
	(If ``--compare_feature_genome`` and ``--annotation_files`` are assigned). 

	**detect_lib:** This transcript is detected by Tex-treated libraries or fragmented/conventional libraries.

	**best_avg_coverage:** The average coverage of the highest expressed library within this transcript.

.. _terminator:

terminator (terminator detection)
---------------------------------

``terminator`` will predict the rho-independent terminators. ``ANNOgesic`` combines the results of 
two methods in order to get more reliable candidates. The first method is using `TranstermHP <http://transterm.cbcb.umd.edu/>`_.
The other one detects the specific secondary structure between converging pairs  
of transcripts and CDSs. ``ANNOgesic`` can check the coverages in order to generate the terminators 
which have coverage significant decrease.

- **Required tools**

`TranstermHP <http://transterm.cbcb.umd.edu/>`_

**RNAfold** of `ViennaRNA <http://www.tbi.univie.ac.at/RNA/>`_.

- **Required files**

**Gff files of the genome annotations**

**Fasta files of the genome sequences**

**Wiggle files of TEX +/- treated libraries or fragmented/conventional libraries**

**Gff files of the transcripts**

- **Arguments**

::

    usage: annogesic terminator [-h] --project_path PROJECT_PATH --fasta_files
                                FASTA_FILES [FASTA_FILES ...] --annotation_files
                                ANNOTATION_FILES [ANNOTATION_FILES ...]
                                --transcript_files TRANSCRIPT_FILES
                                [TRANSCRIPT_FILES ...]
                                [--tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]]
                                [--frag_libs FRAG_LIBS [FRAG_LIBS ...]]
                                [--tex_notex {1,2}]
                                [--replicate_tex REPLICATE_TEX [REPLICATE_TEX ...]]
                                [--replicate_frag REPLICATE_FRAG [REPLICATE_FRAG ...]]
                                [--transterm_path TRANSTERM_PATH]
                                [--expterm_path EXPTERM_PATH]
                                [--rnafold_path RNAFOLD_PATH]
                                [--srna_files SRNA_FILES [SRNA_FILES ...]]
                                [--decrease DECREASE]
                                [--fuzzy_detect_coverage FUZZY_DETECT_COVERAGE]
                                [--fuzzy_within_transcript FUZZY_WITHIN_TRANSCRIPT]
                                [--fuzzy_downstream_transcript FUZZY_DOWNSTREAM_TRANSCRIPT]
                                [--fuzzy_within_gene FUZZY_WITHIN_GENE]
                                [--fuzzy_downstream_gene FUZZY_DOWNSTREAM_GENE]
                                [--highest_coverage HIGHEST_COVERAGE]
                                [--table_best] [--window_size WINDOW_SIZE]
                                [--window_shift WINDOW_SHIFT]
                                [--min_loop_length MIN_LOOP_LENGTH]
                                [--max_loop_length MAX_LOOP_LENGTH]
                                [--min_stem_length MIN_STEM_LENGTH]
                                [--max_stem_length MAX_STEM_LENGTH]
                                [--miss_rate MISS_RATE]
                                [--min_u_tail_length MIN_U_TAIL_LENGTH]
                                [--range_u_tail RANGE_U_TAIL] [--keep_multi_term]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --fasta_files FASTA_FILES [FASTA_FILES ...], -f FASTA_FILES [FASTA_FILES ...]
                            Paths of the genome fasta files.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the genome annotation gff files.
      --transcript_files TRANSCRIPT_FILES [TRANSCRIPT_FILES ...], -a TRANSCRIPT_FILES [TRANSCRIPT_FILES ...]
                            Paths of the transcript gff files.
      --tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...], -tl TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]
                            TEX+/- wig files. The format is:
                            wig_file_path:TEX+/-(tex or notex):condition_id(intege
                            r):replicate_id(alphabet):strand(+ or -). If multiple
                            wig files need to be assigned, please use spaces to
                            separate the wig files. For example,
                            my_lib_tex_forward.wig:tex:1:a:+
                            my_lib_tex_reverse.wig:tex:1:a:-.
      --frag_libs FRAG_LIBS [FRAG_LIBS ...], -fl FRAG_LIBS [FRAG_LIBS ...]
                            Wig files of RNA-Seq with transcript fragmented. The
                            format is: wig_file_path:frag:condition_id(integer):re
                            plicate_id(alphabet):strand(+ or -). If multiple wig
                            files need to be assigned, please use spaces to
                            separate the wig files. For example,
                            my_lib_frag_forward.wig:frag:1:a:+
                            my_lib_frag_reverse.wig:frag:1:a:-.
      --tex_notex {1,2}, -te {1,2}
                            The value is for TEX+/- libraries to decide the
                            terminator should be detected in both (TEX+ and TEX-)
                            or can be detected in only one library (TEX+ or TEX-).
                            Please assign 1 or 2. Default is 1.
      --replicate_tex REPLICATE_TEX [REPLICATE_TEX ...], -rt REPLICATE_TEX [REPLICATE_TEX ...]
                            This value is the minimal number of replicates that a
                            TSS has to be detected in. The format is
                            $NUMBERofCONDITION_$NUMBERofREPLICATE. If different
                            --replicate_tex values need to be assigned to
                            different conditions, please use spaces to separate
                            them. For example, 1_2 2_2 3_3 means that
                            --replicate_tex is 2 in number 1 and number 2
                            conditions. In number 3 condition, --replicate_tex is
                            3. For assigning the same --replicate_tex to all
                            conditions, just use like all_1 (--replicate_tex is 1
                            in all conditions). Default is all_1.
      --replicate_frag REPLICATE_FRAG [REPLICATE_FRAG ...], -rf REPLICATE_FRAG [REPLICATE_FRAG ...]
                            It is similar to --replicates_tex. This value is for
                            fragmented (or conventional) libraries.
    
    additional arguments:
      --transterm_path TRANSTERM_PATH
                            Path of "transterm" in TransTermHP.
      --expterm_path EXPTERM_PATH
                            Path of expterm.dat for TransTermHP. Default is
                            /usr/local/bin/expterm.dat
      --rnafold_path RNAFOLD_PATH
                            Path of "RNAfold" of Vienna package.
      --srna_files SRNA_FILES [SRNA_FILES ...], -sr SRNA_FILES [SRNA_FILES ...]
                            Paths of sRNA gff files if sRNA information need to be
                            considered as well.
      --decrease DECREASE, -d DECREASE
                            The maximum ratio -- (lowest coverage / highest
                            coverage) within (or nearby) the terminator. If the
                            ratio is smaller than --decrease, the candidate will
                            be considered as highly-confidence terminator. Default
                            is 0.5.
      --fuzzy_detect_coverage FUZZY_DETECT_COVERAGE, -fc FUZZY_DETECT_COVERAGE
                            The extended region (nucleotides) of the terminators
                            for detecting coverage significant drop. For example,
                            the location of terminator is 300-400, and
                            --fuzzy_detect_coverage is 30. If the coverage
                            decrease is detected within 270-430, this candidate is
                            still considered as the terminator which have coverage
                            dramatic decrease. Default is 30.
      --fuzzy_within_transcript FUZZY_WITHIN_TRANSCRIPT, -fut FUZZY_WITHIN_TRANSCRIPT
                            If the candidates are within transcript and the
                            distance (nucleotides) between the end of transcript
                            and terminator is within this value, the candidate
                            will be considered as a terminator. Otherwise, it will
                            be removed. Default is 30.
      --fuzzy_downstream_transcript FUZZY_DOWNSTREAM_TRANSCRIPT, -fdt FUZZY_DOWNSTREAM_TRANSCRIPT
                            The meaning is similar to --fuzzy_within_transcript.
                            This value is for the candidates which are at the
                            downstream of transcript. Default is 30.
      --fuzzy_within_gene FUZZY_WITHIN_GENE, -fuc FUZZY_WITHIN_GENE
                            The meaning is similar to --fuzzy_within_transcript.
                            This value is for gene in stead of transcript. Default
                            is 10.
      --fuzzy_downstream_gene FUZZY_DOWNSTREAM_GENE, -fdg FUZZY_DOWNSTREAM_GENE
                            The meaning is similar to
                            --fuzzy_downstream_transcript. This value is for gene
                            in stead of transcript. Default is 310.
      --highest_coverage HIGHEST_COVERAGE, -hc HIGHEST_COVERAGE
                            The minimum value of the highest coverage of
                            terminator. The low expressed terminator are not
                            included in "best_candidates", but are still in
                            "all_candidates". Default is 10.
      --table_best, -tb     Output table only contains the information of the
                            library which has most significant coverage decrease.
                            Default is False.
      --window_size WINDOW_SIZE, -wz WINDOW_SIZE
                            Window size for searching secondary structure in
                            intergenic region. Default is 100 nts.
      --window_shift WINDOW_SHIFT, -ws WINDOW_SHIFT
                            The number of nucleotides for window shift. Default is
                            20 nts.
      --min_loop_length MIN_LOOP_LENGTH, -ml MIN_LOOP_LENGTH
                            The minimum loop length of terminator. Default is 3
                            nts.
      --max_loop_length MAX_LOOP_LENGTH, -Ml MAX_LOOP_LENGTH
                            The maximum loop length of terminator. Default is 10
                            nts.
      --min_stem_length MIN_STEM_LENGTH, -ms MIN_STEM_LENGTH
                            The minimum stem length of terminator. Default is 4
                            nts.
      --max_stem_length MAX_STEM_LENGTH, -Ms MAX_STEM_LENGTH
                            The maximum stem length of terminator. Default is 20
                            nts.
      --miss_rate MISS_RATE, -mr MISS_RATE
                            The percentage of nucleotides which can be no pair in
                            the stem. Default is 0.25.
      --min_u_tail_length MIN_U_TAIL_LENGTH, -mu MIN_U_TAIL_LENGTH
                            The minimum U-tail length of terminator. Default is 3
                            nts.
      --range_u_tail RANGE_U_TAIL, -ru RANGE_U_TAIL
                            The range (nucleotides) for detection of U-tail. For
                            example, if --range_u_tail is 6 and
                            --min_u_tail_length is 3, and there are 3 Us within 6
                            nts, This candidate will be assigned as the terminator
                            which has poly U-tail. Default is 6.
      --keep_multi_term, -kp
                            Sometimes, one gene is associated with multiple
                            terminators In default, it will only keep the highly-
                            confidence one. This flag can keep all terminators
                            which are associated with the same gene. Default is
                            False.

- **Output files**

Output files are stored in ``$ANNOgesic/output/terminators``. 

The output folders are as following:

**statistics:** Stores statistic files.

	**Terminator detection method with corresponding amounts:** Filename is ``stat_$GENOME.csv``.

	**Comparing terminators with transcripts:** Based on different types of terminators, 
	the files are ``stat_compare_terminator_transcript_$GENOME_all_candidates.csv``, 
	``stat_comparison_terminator_transcript_$GENOME_best.csv`` and ``stat_comparison_terminator_transcript_$GENOME_express.csv``

**transtermhp_results:** Store any output of `TranstermHP <http://transterm.cbcb.umd.edu/>`_.

**gffs:** Store gff files of terminators.

There are four different sub-folders for storing different gff files.

	**all_candidates:** Stores all terminators which ``ANNOgesic`` can detect.

	**expressed_candidates:** Stores the terminators revealing gene expression.

	**best_candidates:** Stores the terminators which reveal gene expression and show dramatic decrease of its coverage.

	**non_expressed_candidates:** Stores the terminators which has no gene expression.

Some useful information can be found in the tags of the attributes within the terminator gff file.
Based on this information, we can know the details of the specific terminator. The tags are as following:

	**method:** By which method the terminator is detected.

	**coverage_decrease:** The terminators coverage reveals dramatic decrease or not.

	**express:** The terminator reveals gene expression or not.

	**diff_coverage:** This value shows the library which reveals strongest coverage decreasing.

	**associated_gene:** Which genes are associated with this terminator.

	**Parent:** This tag presents the parent transcript of the terminator.

**tables:** Stores tables of terminators with more details.

There are four different sub-folders for storing different tables.

	**all_candidates:** Stores all terminators which ``ANNOgesic`` can detect.

        **express_candidates:** Stores the terminators revealing gene expression.

        **best_candidates:** Stores the terminators which reveal gene expression and show dramatic decrease of its coverage.

        **non_expressed_candidates:** Stores the terminators which has no gene expression.

The meanings of the columns are as following:

	**Genome:** Genome name.

	**Name:** Name of this terminator in the gff file.

	**Start:** Staring point of this terminator.

	**End:** End point of this terminator.

	**Strand:** Strand of this terminator.

	**Detect:** This terminator is detected by which method.

	**Associated_gene:** Which genes are associated with this terminator.

	**Associated_transcript:** The parent transcript of this terminator.

	**Coverage_decrease:** This terminator shows dramatic decrease of its coverage or not.

	**Coverage_detail:** Shows the coverage information of the libraries about this terminator. "high" means the highest coverage of the libraries, 
	"low" means the lowest coverage of the libraries, and "diff" represents the difference between "high" and "low". If "No_coverage_decreasing" is showed, 
	it means this terminator reveal gene expression but no coverage decrease. If "NA" is showed, it means that this terminator has no gene expression.

.. _utr:

utr (UTR detection)
-------------------

``utr`` can compare TSSs, CDSs/tRNAs/sRNAs, transcripts and terminators
to generate 5'UTR and 3'UTR. 5'UTRs are based on detecting the regions between TSSs and CDSs/tRNAs/sRNAs. 
3'UTRs are based on detecting the 
regions between the end of the transcripts and CDSs/tRNAs/sRNAs. If the input gff files of TSSs are not computed by 
ANNOgesic, please use ``--tss_source`` to classify TSSs for the analysis.

- **Required files**

**Gff files of the genome annotations**

**Gff files of the TSSs**

**Gff files of the transcripts**

- **Optional input files**

**Gff files of the terminators:** If the information of terminators is needed, the gff files of terminators are required.

- **Arguments**

::

    usage: annogesic utr [-h] --project_path PROJECT_PATH --annotation_files
                         ANNOTATION_FILES [ANNOTATION_FILES ...] --tss_files
                         TSS_FILES [TSS_FILES ...] --transcript_files
                         TRANSCRIPT_FILES [TRANSCRIPT_FILES ...]
                         [--terminator_files TERMINATOR_FILES [TERMINATOR_FILES ...]]
                         [--tss_source] [--base_5utr {both,transcript,TSS}]
                         [--utr_length UTR_LENGTH]
                         [--base_3utr {both,transcript,terminator}]
                         [--terminator_fuzzy TERMINATOR_FUZZY]
                         [--fuzzy_3utr FUZZY_3UTR] [--fuzzy_5utr FUZZY_5UTR]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the genome annotation gff files.
      --tss_files TSS_FILES [TSS_FILES ...], -t TSS_FILES [TSS_FILES ...]
                            Paths of the TSS files.
      --transcript_files TRANSCRIPT_FILES [TRANSCRIPT_FILES ...], -a TRANSCRIPT_FILES [TRANSCRIPT_FILES ...]
                            Paths of the transcript gff files.
      --terminator_files TERMINATOR_FILES [TERMINATOR_FILES ...], -e TERMINATOR_FILES [TERMINATOR_FILES ...]
                            If the paths of terminator files are assigned,
                            terminator will be used to detect 3'UTR.
    
    additional arguments:
      --tss_source, -s      The TSS gff file is generated by ANNOgesic or not.
                            Default is True (from ANNOgesic).
      --base_5utr {both,transcript,TSS}, -b5 {both,transcript,TSS}
                            The information for detection of 5'UTR. It can be
                            "TSS" or "transcript" or "both". Default is both.
      --utr_length UTR_LENGTH, -l UTR_LENGTH
                            The maximum UTR length. Default is 300.
      --base_3utr {both,transcript,terminator}, -b3 {both,transcript,terminator}
                            please assign the information for detection of 3'UTR.
                            It can be "transcript" or "terminator" or "both".
                            Default is transcript.
      --terminator_fuzzy TERMINATOR_FUZZY, -ef TERMINATOR_FUZZY
                            The minimum fuzzy value for comparing transcript and
                            terminator. Default is 30.
      --fuzzy_3utr FUZZY_3UTR, -f3 FUZZY_3UTR
                            The fuzzy value of transcript for 3'UTR detection.
      --fuzzy_5utr FUZZY_5UTR, -f5 FUZZY_5UTR
                            The fuzzy value of transcript for 5'UTR detection.

- **Output files**

Output files of 5'UTRs are stored in ``$ANNOgesic/output/UTRs/5UTRs``.

Output files of 3'UTRs are stored in ``$ANNOgesic/output/UTRs/3UTRs``.

The output folders are as following:

**gffs:** Stores gff files of the 5'UTR/3'UTR. 
Some useful information can be found in the tags of the attributes within the UTR gff file. 
Based on this information, we can know the details of the specific UTR. The tags are as following:

	**length:** UTR length.
	
	**associated_cds:** Which CDSs/rRNAs/tRNAs are associated with this UTR.
	
	**associated_gene:** Which genes are associated with this UTR.
	
	**Parent:** Shows the parent transcript of this UTR.
	
	**associated_tss:** Which TSSs are associated with this 5'UTR.
	
	**tss_type:** What types of TSSs are associated with this 5'UTR.
	
	**associated_term:** Which terminators are associated with this 3'UTR.

**statiatics:** ``$GFFNAME_$GENOME_$UTRTYPE_length.png`` is the distribution of the UTR length.

.. _srna:

srna (sRNA detection)
---------------------
``srna`` can predict different types of sRNAs. For intergenic and antisense sRNA, it 
is detected via comparison of the transcripts and annotation profiles, as well as coverage files. 
For UTR-derived sRNA, the detection is based on the TSSs, processing sites, 
transcripts, genome annotations and coverage files. Further filters like folding free energy change, 
BLAST to nr database and sRNA database can be set as well.

- **Required files**

**Gff files of the genome annotations**

**Gff files of the transcripts**

**Wiggle files of the fragmented/conventional or TEX+/- libraries:** Please check the section 
:ref:`The input format of libraries for running ANNOgesic`.

- **Optional input files**

**Gff files of the TSSs:** If you want to detect the UTR-derived sRNAs, it is necessary to input
TSS information. If you don't want to detect UTR-derived sRNAs, TSS information still can be provided as a filter.
We strongly recommend input this file.

**Gff files of processing sites:** For checking the sRNAs which end with processing sites. Moreover,
Some 3'UTR-derived and interCDS-derived sRNA candidates start
from processing sites not TSSs. If you don't want to detect UTR-derived sRNAs,
This information still can be provided to increase the accuracy, especially for some
long non-coding regions. We strongly recommend input this file if you want to detect UTR-derived sRNAs.

**Promoter tables:** Information of the promoter motifs can be used for prioritizing sRNA candidates via 
promoters and sRNA coverage. The format should be as following:

===========  ============  ==========  =======
Genome       TSS_position  TSS_strand  Motif
-----------  ------------  ----------  -------
NC_000915.1  237118        \-          MOTIF_1
NC_000915.1  729009        \-          MOTIF_1
===========  ============  ==========  =======

First irow is header of the table, the last column is the name of promoter motif.
If subcommand ``promoter`` was implemented before, the table will be generated automatically.
Please refer to the section :ref:`promoter`.

- **Filers with the corresponding input files and tools**

There are some filters which can improve the prediction. The user can assign the information to remove false positive. 
If the information is not assigned to be a filter, it still can input to the module. Then, the information will 
be shown in the output files, but this information is not considered as a filter. For an example, if terminator association 
is not assigned to be a filter, the user still can specify the path of terminator gff files. The associated terminators 
will be shown in output gff files and tables, but the sRNA candidates which are not associated with terminators will 
still be included. Following is the filter names with the required files and tools.

**Secondary structure:** Remove the false positive by checking the folding energy change of secondary structure.

	**Required tools:**

		`ViennaRNA <http://www.tbi.univie.ac.at/RNA/>`_

	**Required files:**

		**Fasta files of genome sequences**

**TSS:** Remove the candidates which are not associated with TSSs.

	**Required files:**

		**Gff files of TSSs**

**Searching sRNA candidate in sRNA database:** If homology of this sRNA candidate can be found in sRNA database, 
this candidate will be included to the result without considering other filters.

	**Required tools:**

		`Blast+ <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`_

	**Required files:**

		**sRNA database:** Such as `BSRD <http://www.bac-srna.org/BSRD/index.jsp>`_. 
		Format of the header should be ``$ID|$GENOME|$SRNANAME``. For an example, 
		``srn_4840|S._aureus_NCTC8325|RsaOV`` The ID is srn_4840, 
		the strain of this sRNA is S._aureus_NCTC8325 and the name of sRNA is RsaOV.
		If the format of the header is not correct, an error or non-sense results will occur.
		If you want to use BSRD with proper headers, you can download it from our
		`Git repository <https://github.com/Sung-Huan/ANNOgesic/tree/master/database>`_ easily.


**Searching sRNA candidate in nr database:** If homologs of this sRNA candidates can be found in nr database and the hit numbers are more than ``--cutoff_nr_hit``,
this candidates will be removed.

	**Required tools:**

		`Blast+ <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`_

	**Required files:**

		**nr database:** The file can be download from `nr database <ftp://ftp.ncbi.nih.gov/blast/db/FASTA/>`_.
	
**Terminator:** Remove the candidates which are not associated with terminators.

	**Required files:**

		**Gff files of the terminators**

**sORF:** Remove the candidates which overlap sORF.

	**Required files:**

		**Gff files of the sORFs**

**Promoter:** Remove the candidates which are not associated with promoter motif.

	**Required files:**

		**Tables of the promoters:** Please check the Promoter Tables of this section.

- **Arguments**

::

    usage: annogesic srna [-h] --project_path PROJECT_PATH [--utr_derived_srna]
                          [--filter_info {tss,sec_str,blast_nr,blast_srna,sorf,term,promoter,none} [{tss,sec_str,blast_nr,blast_srna,sorf,term,promoter,none} ...]]
                          --transcript_files TRANSCRIPT_FILES
                          [TRANSCRIPT_FILES ...] --annotation_files
                          ANNOTATION_FILES [ANNOTATION_FILES ...]
                          [--tss_files TSS_FILES [TSS_FILES ...]]
                          [--processing_site_files PROCESSING_SITE_FILES [PROCESSING_SITE_FILES ...]]
                          [--terminator_files TERMINATOR_FILES [TERMINATOR_FILES ...]]
                          [--fasta_files FASTA_FILES [FASTA_FILES ...]]
                          [--compute_sec_structures]
                          [--promoter_tables PROMOTER_TABLES [PROMOTER_TABLES ...]]
                          [--promoter_names PROMOTER_NAMES [PROMOTER_NAMES ...]]
                          [--sorf_files SORF_FILES [SORF_FILES ...]]
                          [--srna_database_path SRNA_DATABASE_PATH]
                          [--nr_database_path NR_DATABASE_PATH]
                          [--tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]]
                          [--frag_libs FRAG_LIBS [FRAG_LIBS ...]]
                          [--tex_notex {1,2}]
                          [--replicate_tex REPLICATE_TEX [REPLICATE_TEX ...]]
                          [--replicate_frag REPLICATE_FRAG [REPLICATE_FRAG ...]]
                          [--rnafold_path RNAFOLD_PATH]
                          [--relplot_path RELPLOT_PATH]
                          [--mountain_path MOUNTAIN_PATH]
                          [--blastn_path BLASTN_PATH] [--blastx_path BLASTX_PATH]
                          [--makeblastdb_path MAKEBLASTDB_PATH]
                          [--parallel_blast PARALLEL_BLAST] [--tss_source]
                          [--tss_intergenic_fuzzy TSS_INTERGENIC_FUZZY]
                          [--tss_5utr_fuzzy TSS_5UTR_FUZZY]
                          [--tss_3utr_fuzzy TSS_3UTR_FUZZY]
                          [--tss_intercds_fuzzy TSS_INTERCDS_FUZZY]
                          [--terminator_fuzzy_in_srna TERMINATOR_FUZZY_IN_SRNA]
                          [--terminator_fuzzy_out_srna TERMINATOR_FUZZY_OUT_SRNA]
                          [--min_length MIN_LENGTH] [--max_length MAX_LENGTH]
                          [--min_intergenic_tex_coverage MIN_INTERGENIC_TEX_COVERAGE]
                          [--min_intergenic_notex_coverage MIN_INTERGENIC_NOTEX_COVERAGE]
                          [--min_intergenic_fragmented_coverage MIN_INTERGENIC_FRAGMENTED_COVERAGE]
                          [--min_complete_5utr_transcript_coverage MIN_COMPLETE_5UTR_TRANSCRIPT_COVERAGE]
                          [--min_antisense_tex_coverage MIN_ANTISENSE_TEX_COVERAGE]
                          [--min_antisense_notex_coverage MIN_ANTISENSE_NOTEX_COVERAGE]
                          [--min_antisense_fragmented_coverage MIN_ANTISENSE_FRAGMENTED_COVERAGE]
                          [--min_utr_tex_coverage MIN_UTR_TEX_COVERAGE]
                          [--min_utr_notex_coverage MIN_UTR_NOTEX_COVERAGE]
                          [--min_utr_fragmented_coverage MIN_UTR_FRAGMENTED_COVERAGE]
                          [--min_all_utr_coverage MIN_ALL_UTR_COVERAGE]
                          [--cutoff_energy CUTOFF_ENERGY] [--mountain_plot]
                          [--nr_format] [--srna_format] [--table_best]
                          [--decrease_intergenic_antisense DECREASE_INTERGENIC_ANTISENSE]
                          [--decrease_utr DECREASE_UTR]
                          [--fuzzy_intergenic_antisense FUZZY_INTERGENIC_ANTISENSE]
                          [--fuzzy_utr FUZZY_UTR] [--cutoff_nr_hit CUTOFF_NR_HIT]
                          [--blast_e_nr BLAST_E_NR] [--blast_e_srna BLAST_E_SRNA]
                          [--detect_srna_in_cds]
                          [--overlap_percent_cds OVERLAP_PERCENT_CDS]
                          [--ignore_hypothetical_protein]
                          [--ranking_time_promoter RANKING_TIME_PROMOTER]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --utr_derived_srna, -u
                            Assign to detect UTR-derived sRNA. Default is False.
      --filter_info {tss,sec_str,blast_nr,blast_srna,sorf,term,promoter,none} [{tss,sec_str,blast_nr,blast_srna,sorf,term,promoter,none} ...], -d {tss,sec_str,blast_nr,blast_srna,sorf,term,promoter,none} [{tss,sec_str,blast_nr,blast_srna,sorf,term,promoter,none} ...]
                            The filters for improving sRNA detection: 1. tss (sRNA
                            has to start with a TSS), 2. sec_str (free energy
                            change of secondary structure (normalized by length)
                            has to be smaller than --cutoff_energy), 3. blast_nr
                            (the number of the homology can not be found more than
                            --cutoff_nr_hit in the non-redundant database), 4.
                            blast_srna (as long as the homology can be found in
                            sRNA database, the candidates will be included to best
                            candidtes without considering other filters), 5. sorf
                            (sRNA must not overlap sORFs), 6. term (sRNA has to be
                            associated with a terminator), 7. promoter (sRNA has
                            to be associated with a promoter motif). For using
                            multiple filters, please separated them by spaced. If
                            blast_srna was assigned, the headers of sequences in
                            sRNA database should be $ID|$GENOME|$SRNANAME. "tss
                            sec_str blast_nr blast_srna" are recommended to be
                            used. If "none" is assigned, no filters are applied.
                            Default is tss sec_str blast_nr blast_srna.
      --transcript_files TRANSCRIPT_FILES [TRANSCRIPT_FILES ...], -a TRANSCRIPT_FILES [TRANSCRIPT_FILES ...]
                            Paths of the transcript files.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the genome annotation gff files.
      --tss_files TSS_FILES [TSS_FILES ...], -t TSS_FILES [TSS_FILES ...]
                            Paths of TSS gff files. For detecting UTR-derived sRNA
                            or "tss" in --filter_info, TSS gff files MUST be
                            provided.
      --processing_site_files PROCESSING_SITE_FILES [PROCESSING_SITE_FILES ...], -p PROCESSING_SITE_FILES [PROCESSING_SITE_FILES ...]
                            Paths of processing site gff files. It can improve the
                            detection of UTR-derived sRNAs.
      --terminator_files TERMINATOR_FILES [TERMINATOR_FILES ...], -e TERMINATOR_FILES [TERMINATOR_FILES ...]
                            Paths of terminator gff files.
      --fasta_files FASTA_FILES [FASTA_FILES ...], -f FASTA_FILES [FASTA_FILES ...]
                            paths of fasta files of reference genome, If
                            "sec_str", "blast_nr" or "blast_srna" is assigned to
                            --filter_info, fasta files are required.
      --compute_sec_structures, -cs
                            Computing secondary structures of sRNAs. Default is
                            False.
      --promoter_tables PROMOTER_TABLES [PROMOTER_TABLES ...], -pt PROMOTER_TABLES [PROMOTER_TABLES ...]
                            Paths of promoter tables, The format of the table is
                            $GENOME $TSS_POSITION $TSS_STRAND $PROMOTER_NAME.
      --promoter_names PROMOTER_NAMES [PROMOTER_NAMES ...], -pn PROMOTER_NAMES [PROMOTER_NAMES ...]
                            If --promoter_tables is provided, please assign the
                            promoter name (the last column of promoter table). For
                            multiple promoters, please put spaces between the
                            promoters. Default is None.
      --sorf_files SORF_FILES [SORF_FILES ...], -O SORF_FILES [SORF_FILES ...]
                            Paths of sORF gff files
      --srna_database_path SRNA_DATABASE_PATH, -sd SRNA_DATABASE_PATH
                            Path of sRNA database with proper headers of
                            sequences. Format of the header should be
                            $ID|$GENOME|$NAME. Please check
                            https://github.com/Sung-Huan/ANNOgesic/blob/master/dat
                            abase/sRNA_database_BSRD.fa
      --nr_database_path NR_DATABASE_PATH, -nd NR_DATABASE_PATH
                            Path of nr database
      --tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...], -tl TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]
                            TEX+/- wig files. The format is:
                            wig_file_path:TEX+/-(tex or notex):condition_id(intege
                            r):replicate_id(alphabet):strand(+ or -). If multiple
                            wig files need to be assigned, please use spaces to
                            separate the wig files. For example,
                            my_lib_tex_forward.wig:tex:1:a:+
                            my_lib_tex_reverse.wig:tex:1:a:-.
      --frag_libs FRAG_LIBS [FRAG_LIBS ...], -fl FRAG_LIBS [FRAG_LIBS ...]
                            Wig files of RNA-Seq with transcript fragmented. The
                            format is: wig_file_path:frag:condition_id(integer):re
                            plicate_id(alphabet):strand(+ or -). If multiple wig
                            files need to be assigned, please use spaces to
                            separate the wig files. For example,
                            my_lib_frag_forward.wig:frag:1:a:+
                            my_lib_frag_reverse.wig:frag:1:a:-.
      --tex_notex {1,2}, -te {1,2}
                            If TEX+/- libraries is assigned, this value is that a
                            sRNA should be detected in both (TEX+ and TEX-) or can
                            be detected in only one library (TEX+ or TEX-). Please
                            assign 1 or 2. Default is 2.
      --replicate_tex REPLICATE_TEX [REPLICATE_TEX ...], -rt REPLICATE_TEX [REPLICATE_TEX ...]
                            This value is the minimal number of replicates that a
                            TSS has to be detected in. The format is
                            $NUMBERofCONDITION_$NUMBERofREPLICATE. If different
                            --replicate_tex values need to be assigned to
                            different conditions, please use spaces to separate
                            them. For example, 1_2 2_2 3_3 means that
                            --replicate_tex is 2 in number 1 and number 2
                            conditions. In number 3 condition, --replicate_tex is
                            3. For assigning the same --replicate_tex to all
                            conditions, just use like all_1 (--replicate_tex is 1
                            in all conditions). Default is all_1.
      --replicate_frag REPLICATE_FRAG [REPLICATE_FRAG ...], -rf REPLICATE_FRAG [REPLICATE_FRAG ...]
                            It is similar to --replicates_tex. This value is for
                            fragmented (or conventional) libraries.
    
    additional arguments:
      --rnafold_path RNAFOLD_PATH
                            Path of RNAfold.
      --relplot_path RELPLOT_PATH
                            Path of relplot.pl in Vienna package.
      --mountain_path MOUNTAIN_PATH
                            Path of mountain.pl in Vienna package.
      --blastn_path BLASTN_PATH
                            Path of blastn in BLAST+ package.
      --blastx_path BLASTX_PATH
                            Path of blastx in BLAST+ package.
      --makeblastdb_path MAKEBLASTDB_PATH
                            Path of makeblastdb in BLAST+ package.
      --parallel_blast PARALLEL_BLAST, -pb PARALLEL_BLAST
                            The number of parallel runs. Default is 10.
      --tss_source, -ts     The TSS gff files are generated from ANNOgesic or not.
                            Default is True (from ANNOgesic).
      --tss_intergenic_fuzzy TSS_INTERGENIC_FUZZY, -ft TSS_INTERGENIC_FUZZY
                            The fuzzy value of comparing TSS and transcript for
                            detecting intergenic sRNA. Default is 3.
      --tss_5utr_fuzzy TSS_5UTR_FUZZY, -f5 TSS_5UTR_FUZZY
                            The fuzzy value of comparing TSS and transcript for
                            detecting 5'UTR-derived sRNA.The input type can be
                            percentage ("p") or the real amount of reads ("n").
                            For example, p_0.05 means the fuzzy is 5 percent of
                            the length of 5'UTR. n_10 means the fuzzy is 10 nts.
                            Default is n_3.
      --tss_3utr_fuzzy TSS_3UTR_FUZZY, -f3 TSS_3UTR_FUZZY
                            It is similar to --tss_5utr_fuzzy. This value is for
                            3'UTR-derived sRNA instead of 5'UTR-derived sRNA.
                            Default is p_0.04.
      --tss_intercds_fuzzy TSS_INTERCDS_FUZZY, -fc TSS_INTERCDS_FUZZY
                            It is similar to --tss_5utr_fuzzy. This value is for
                            interCDS-derived sRNA instead of 5'UTR-derived sRNA.
                            Default is p_0.04.
      --terminator_fuzzy_in_srna TERMINATOR_FUZZY_IN_SRNA, -efi TERMINATOR_FUZZY_IN_SRNA
                            The fuzzy value for comparing terminator and
                            transcript for detecting the terminator which is
                            within sRNA. Default is 30.
      --terminator_fuzzy_out_srna TERMINATOR_FUZZY_OUT_SRNA, -efo TERMINATOR_FUZZY_OUT_SRNA
                            It is similar to --terminator_fuzzy_in_sRNA. This
                            value is for the terminator which is outside of sRNA.
                            Default is 30.
      --min_length MIN_LENGTH, -lm MIN_LENGTH
                            The minimum sRNA length. Default is 30.
      --max_length MAX_LENGTH, -lM MAX_LENGTH
                            The maximum sRNA length. Default is 500.
      --min_intergenic_tex_coverage MIN_INTERGENIC_TEX_COVERAGE, -it MIN_INTERGENIC_TEX_COVERAGE
                            The minimum average coverage of intergenic sRNAs for
                            TEX+ libraries. This value is based on different types
                            of TSSs. The order of numbers is
                            "Primary,Secondary,Internal,Antisense,Orphan". For
                            example, 0,0,0,50,10 means that antisense TSS (minimum
                            coverage is 50) and orphan TSS (minimum coverage is
                            10) are used for sRNA prediction. The other types of
                            TSSs will not be used for the detection (assign to 0).
                            If TSS information is not provided, the lowest value
                            would be the general cutoff for the prediction.
                            Default is 0,0,0,40,20.
      --min_intergenic_notex_coverage MIN_INTERGENIC_NOTEX_COVERAGE, -in MIN_INTERGENIC_NOTEX_COVERAGE
                            It is similar to --min_intergenic_tex_coverage. This
                            value is for TEX- libraries. Default is 0,0,0,30,10.
      --min_intergenic_fragmented_coverage MIN_INTERGENIC_FRAGMENTED_COVERAGE, -if MIN_INTERGENIC_FRAGMENTED_COVERAGE
                            It is similar to --min_intergenic_tex_coverage. This
                            value is for fragmented (or conventional) libraries.
                            Default is 400,200,0,50,20.
      --min_complete_5utr_transcript_coverage MIN_COMPLETE_5UTR_TRANSCRIPT_COVERAGE, -ib MIN_COMPLETE_5UTR_TRANSCRIPT_COVERAGE
                            Several primary/secondary TSSs are also associated
                            with a complete transcript containing no
                            CDSs/tRNA/rRNA in 5'UTR of the following CDS which is
                            located in another transcript. In order to detect the
                            sRNA candidates in these transcripts, please assign
                            the minimum average coverage of the sRNA candidates.
                            The format is $TEX,$NOTEX,$FRAG. For example,
                            200,100,100 means that the minimum average coverage is
                            200 for TEX+ libraries, 100 for TEX- and fragmented
                            (or conventional) libraries. Default is 30,20,30.
      --min_antisense_tex_coverage MIN_ANTISENSE_TEX_COVERAGE, -at MIN_ANTISENSE_TEX_COVERAGE
                            It is similar to --min_intergenic_tex_coverage. This
                            value is for antisense in stead of intergenic. Default
                            is 0,0,0,40,20.
      --min_antisense_notex_coverage MIN_ANTISENSE_NOTEX_COVERAGE, -an MIN_ANTISENSE_NOTEX_COVERAGE
                            It is similar to --min_intergenic_notex_coverage. This
                            value is for antisense in stead of intergenic. Default
                            is 0,0,0,30,10.
      --min_antisense_fragmented_coverage MIN_ANTISENSE_FRAGMENTED_COVERAGE, -af MIN_ANTISENSE_FRAGMENTED_COVERAGE
                            It is similar to --min_intergenic_fragmented_coverage.
                            This value is for antisense in stead of intergenic.
                            Default is 400,200,0,50,20.
      --min_utr_tex_coverage MIN_UTR_TEX_COVERAGE, -ut MIN_UTR_TEX_COVERAGE
                            The minimum average coverage of UTR-derived sRNA
                            candidates in TEX+ libraries. The input can be
                            assigned by the percentile ("p") or real number of
                            coverage ("n"). The order of numbers is
                            "5'UTR,3'UTR,interCDS". For example,
                            "p_0.7,p_0.5,p_0.5" means that 70 percentile of all
                            5'UTR coverages is used for the minimum coverage of
                            5'UTR-derived sRNA, median of all 3'UTR and interCDS
                            coverages is used for minimum coverage of 3'UTR and
                            interCDS-derived sRNA. Default is p_0.8,p_0.6,p_0.7.
      --min_utr_notex_coverage MIN_UTR_NOTEX_COVERAGE, -un MIN_UTR_NOTEX_COVERAGE
                            It is similar to --min_utr_tex_coverage. This value is
                            for TEX- libraries. Default is p_0.7,p_0.5,p_0.6.
      --min_utr_fragmented_coverage MIN_UTR_FRAGMENTED_COVERAGE, -uf MIN_UTR_FRAGMENTED_COVERAGE
                            It is similar to --min_utr_tex_coverage. This value is
                            for fragmented (or conventional) libraries. Default is
                            p_0.7,p_0.5,p_0.6.
      --min_all_utr_coverage MIN_ALL_UTR_COVERAGE, -mu MIN_ALL_UTR_COVERAGE
                            The minimum coverage of UTR-derived sRNA. The coverage
                            of UTR-derived sRNA should not only fit the
                            --min_utr_TEX_coverage, --min_utr_noTEX_coverage and
                            --min_utr_fragmented_coverage, but also this value.
                            Default is 50.
      --cutoff_energy CUTOFF_ENERGY, -ce CUTOFF_ENERGY
                            If "sec_str" is included in --filter_info, please
                            assign the maximum folding energy change (normalized
                            by length of gene). Default is -0.05.
      --mountain_plot, -m   Generating mountain plot of sRNA candidate. Default is
                            False.
      --nr_format, -nf      Format nr database. Default is False.
      --srna_format, -sf    Format sRNA database. Default is False.
      --table_best, -tb     The output table of sRNA candidates only contains the
                            information of the highest expressed library. Default
                            is False.
      --decrease_intergenic_antisense DECREASE_INTERGENIC_ANTISENSE, -di DECREASE_INTERGENIC_ANTISENSE
                            This value is for detecting the coverage decrease in
                            intergenic/antisense transcript. If the length of
                            intergenic transcript is longer than the --max_length,
                            it will searching the coverages of the transcript. If
                            the ratio -- (the lowest coverage / the highest
                            coverage) of the transcript is smaller than this value
                            and the length is within a given range, the transcript
                            will be considered as a sRNA as well.Default is 0.1.
      --decrease_utr DECREASE_UTR, -du DECREASE_UTR
                            It is similar to --decrease_intergenic_antisense. This
                            value is for UTR-derived sRNA. Default is 0.05.
      --fuzzy_intergenic_antisense FUZZY_INTERGENIC_ANTISENSE, -fi FUZZY_INTERGENIC_ANTISENSE
                            The fuzzy nucleotides for detecting the coverage
                            decrease (please check
                            --decrease_intergenic_antisense). For example, the
                            location of intergenic sRNA is 300-400, and
                            --fuzzy_intergenic_antisense is 30. The searching
                            region is 270-430. Default is 10.
      --fuzzy_utr FUZZY_UTR, -fu FUZZY_UTR
                            It is similar to --fuzzy_intergenic_antisense. This is
                            for UTR-derived sRNAs. Default is 10.
      --cutoff_nr_hit CUTOFF_NR_HIT, -cn CUTOFF_NR_HIT
                            The maximum hits number in nr database. Default is 0.
      --blast_e_nr BLAST_E_NR, -en BLAST_E_NR
                            The maximum e value for searching in nr database.
                            Default is 0.0001.
      --blast_e_srna BLAST_E_SRNA, -es BLAST_E_SRNA
                            The maximum e value for searching in sRNA database.
                            Default is 0.0001.
      --detect_srna_in_cds, -ds
                            Searching sRNA in CDS (e.g. the genome annotation is
                            not correct). More sRNA candidates which overlap with
                            CDS will be detected. Default is False.
      --overlap_percent_cds OVERLAP_PERCENT_CDS, -oc OVERLAP_PERCENT_CDS
                            The maximum ratio of overlapping between CDS and sRNA
                            candidates. It only works if --detect_srna_in_cds is
                            true. Default is 0.5
      --ignore_hypothetical_protein, -ih
                            For ignoring the hypothetical proteins in genome
                            annotation file. Default is False.
      --ranking_time_promoter RANKING_TIME_PROMOTER, -rp RANKING_TIME_PROMOTER
                            If --promoter_tables is provided, the information of
                            promoter can be use for ranking sRNA candidates. The
                            ranking score is --ranking_time_promoter * average
                            coverage. For example, a sRNA candidate which is
                            associated with a promoter and its average coverage is
                            10. If --ranking_time_promoter is 2, the ranking score
                            will be 20 (2*10). For the candidate which are not
                            associated with promoter, the --ranking_time_promoter
                            will be 1. Therefore, --ranking_time_promoter can not
                            be smaller than 1. Default is 2.

- **Output files**

Output files are stored in ``$ANNOgesic/output/sRNAs``. the output folders and files are following:

**sRNA_2d_$GENOME:** The secondary structures of all sRNA candidates.

**sRNA_seq_$GENOME:** The sequences of all sRNA candidates.

**blast_results_and_misc:** Stores the results of blast.

	**nr_blast_$GENOME.txt:** output of BLAST for the nr database.

	**sRNA_blast_$GENOME.txt:** output of BLAST for the sRNA database.

**figs:** Stores the figures about secondary structures of sRNAs.

	**mountain_plots:** Stores mountain plots of the sRNA candidates. Filename is as ``srna10_NC_009839.1_335339_335435_+_mountain.pdf``.
	"srna10", "NC_009839.1", "335339", "335435", "+" are ID of sRNA gff file, genome name, starting point, end point and strand, respectively.

	**sec_plots:** Stores the secondary structure plots of sRNA candidates. 
	Filename of is as ``srna10_NC_009839.1_335339_335435_+_rss.ps``.
	"srna10", "NC_009839.1", "335339", "335435", "+" are ID of sRNA gff file, genome name, starting point, end point and strand, respectively.

	**dot_plots:** Stores the dot plots of sRNA candidates. 
	Filename of dot plot is as ``srna10_NC_009839.1_335339_335435_+_dp.ps``.
	"srna10", "NC_009839.1", "335339", "335435", "+" are ID of sRNA gff file, genome name, starting point, end point and strand, respectively.

**statistics:** Stores statistics files. ``stat_$GENOME_sRNA_blast.csv`` is the analysis result of BLAST for sRNA databases.
``stat_sRNA_class_Staphylococcus_aureus_HG003.csv`` is the classification of sRNA candidates.

**TSS_classes:** If the TSSs are not computed by ANNOgesic, ``TSS_classes`` will be generated for classification of TSS.
TSS gff files with TSS types will be stored here.

**tables:** Stores sRNA tables with more details. There are also some sub-folders:

	**for_classes:** Stores the results based on different sRNA classes. The information of sRNA classes can be found in ``stat_sRNA_class_$GENOME.csv``.

	**best_candidates:** Stores the best results of sRNAs after filtering.

	**all_candidates:** Stores all candidates without filtering.

The meanings of the columns are as following:

	**Rank:** Ranking number of this sRNA. 

	**Genome:** Genome name.

	**Name:** sRNA Name which is shown in gff file.

	**Start:** Starting point of this sRNA.

	**End:** End point of this sRNA.

	**Strand:** Strand of this sRNA.

	**Start_with_TSS/Cleavage_site:** This sRNA starts with which TSS or cleavage site.

	**End_with_cleavage:** If the sRNA ends with a cleavage site, the information of this cleavage site will be showed here.

	**Candidates:** Position of this sRNA.

	**Lib_type:** This sRNA is detected by TEX+/- or fragmented/conventional library.

	**Best_avg_coverage:** Based on coverage of all libraries, The best average coverage of this sRNA will be showed here.

	**Best_highest_coverage:** Based on coverage of all libraries, The highest average coverage of this sRNA will be showed here.

	**Best_lowest_coverage:** Based on coverage of all libraries, The lowest average coverage of this sRNA will be showed here.

	**Track/Coverage:** Shows the coverage information of the libraries about this sRNA. "high" means the highest coverage of the libraries,
        "low" means the lowest coverage of the libraries, and "avg" represents the average coverage of this sRNA.

	**Normalized_secondary_energy_change(by_length):** Secondary folding energy change (normalized by length) of this sRNA.

	**sRNA_types:** Shows the sRNA type.

	**Conflict_sORF:** If this sRNA overlaps sORF, the overlapped sORF will be showed here.

	**nr_hit_number:** The hit numbers of this sRNA in nr database.

	**sRNA_hit_number:** The hit numbers of this sRNA in sRNA database.

	**nr_hit_top3|ID|e-value:** The top 3 hits of this sRNA in nr database will be showed here. The information includes protein name, ID and e-value.

	**sRNA_hit|e-value:** If the homology of this sRNA can be found in sRNA database, the information will be showed here.

	**Overlap_CDS:** If the sRNA overlaps CDS, the information of the overlapped CDS will be showed here.

	**Overlap_percent:** If the sRNA overlaps CDS, the percentage of overlap between sRNA and CDS will be showed here.

	**End_with_terminator:** The terminator which is associated with this sRNA.

	**Associated_promoter:** The promoter which is associated with this sRNA.

	**sRNA_length:** sRNA length.

**gffs:** Stores gff files of the sRNA. There are also some sub-folders:

	**for_classes:** Stores the results based on different sRNA classes.

	**best_candidates:** Stores the best results of sRNAs after filtering.

	**all_candidates:** Stores all candidates without filtering.

Some useful information can be found in the tags of the attributes within the sRNA gff file.
Based on this information, we can know the details of the specific sRNA. The tags are as following:

	**sRNA_type:** This sRNA is from 5'UTR, 3'UTR, interCDS, intergenic, antisense or within CDS.

	**with_TSS:** Which TSSs are related to this sRNA.

	**sORF:** Which sORFs overlap this sRNA.

	**sRNA_hit:** Blast hit numbers of this sRNA in sRNA database.

	**nr_hit:** Blast hit numbers of this sRNA in nr database.

	**2d_energy**: Normalized (by the length of sRNA) folding energy change of the sRNA secondary structure.

	**with_term:** Terminators which are associated with the sRNA candidate.

	**end_cleavage:** If this sRNA ends with a cleavage site, information of the cleavage site will be showed here.

	**overlap_cds:** This sRNA overlaps CDS or not.

	**overlap_percentage:** If this sRNA overlap CDS. The percentage of the overlap between CDS and sRNA will be showed here.

	**promoter:** Promoters which are associated with the sRNA.

.. _sorf:

sorf (sORF detection)
---------------------
``sorf`` can detect sORF based on searching ribosome binding sites (Shine-Dalgarno sequence), 
start codons and stop codons within the non-coding expressed regions.
Since these regions may be sRNAs or sORFs, Comparison between sORFs and sRNAs can be done by this subcommand as well. 
If multiple sORFs overlap with each other, this subcommand will merge them to be one sORF. Therefore, one region may contain more than one sORF. 
Position of the start codon which listed in output table is assigned by the first nucleotide. The position of stop codon is assigned by the last nucleotide. 
Moreover, one region may contain different frame shifts. 
Ex: (200, 202, 203) are the positions of three start codons and (241, 243) are two stop codons in 
a small transcript. There are three sORF candidates (200-241, 203-241 and 202-243).

- **Required files**

**Gff files of the genome annotations**
**Gff files of the transcripts**

**Wiggle files of TEX+/- or fragmented/conventional libraries:** Please refer to the section :ref:`The input format of libraries for running ANNOgesic`.

**fasta files of the genome sequences**

- **Optional input files**

**Gff files of the TSSs:** For checking the sORFs start from TSS or not. We strongly recommend to input this file. 

**Gff files of sRNAs:** For checking the overlap of sRNAs and sORFs.

- **Arguments**

::

    usage: annogesic sorf [-h] --project_path PROJECT_PATH [--utr_derived_sorf]
                          --fasta_files FASTA_FILES [FASTA_FILES ...]
                          --transcript_files TRANSCRIPT_FILES
                          [TRANSCRIPT_FILES ...] --annotation_files
                          ANNOTATION_FILES [ANNOTATION_FILES ...]
                          [--tss_files TSS_FILES [TSS_FILES ...]]
                          [--srna_files SRNA_FILES [SRNA_FILES ...]]
                          [--tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]]
                          [--frag_libs FRAG_LIBS [FRAG_LIBS ...]]
                          [--tex_notex TEX_NOTEX]
                          [--replicate_tex REPLICATE_TEX [REPLICATE_TEX ...]]
                          [--replicate_frag REPLICATE_FRAG [REPLICATE_FRAG ...]]
                          [--utr_length UTR_LENGTH] [--min_length MIN_LENGTH]
                          [--max_length MAX_LENGTH]
                          [--cutoff_intergenic_coverage CUTOFF_INTERGENIC_COVERAGE]
                          [--cutoff_antisense_coverage CUTOFF_ANTISENSE_COVERAGE]
                          [--cutoff_5utr_coverage CUTOFF_5UTR_COVERAGE]
                          [--cutoff_3utr_coverage CUTOFF_3UTR_COVERAGE]
                          [--cutoff_intercds_coverage CUTOFF_INTERCDS_COVERAGE]
                          [--cutoff_base_coverage CUTOFF_BASE_COVERAGE]
                          [--table_best]
                          [--start_codon START_CODON [START_CODON ...]]
                          [--stop_codon STOP_CODON [STOP_CODON ...]]
                          [--min_rbs_distance MIN_RBS_DISTANCE]
                          [--max_rbs_distance MAX_RBS_DISTANCE]
                          [--rbs_not_after_tss] [--fuzzy_rbs FUZZY_RBS]
                          [--print_all_combination] [--best_no_srna]
                          [--best_no_tss]
                          [--ignore_hypothetical_protein IGNORE_HYPOTHETICAL_PROTEIN]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --utr_derived_sorf, -u
                            Detecting UTR-derived sORF. Default is False.
      --fasta_files FASTA_FILES [FASTA_FILES ...], -f FASTA_FILES [FASTA_FILES ...]
                            Paths of fasta files of reference genome.
      --transcript_files TRANSCRIPT_FILES [TRANSCRIPT_FILES ...], -a TRANSCRIPT_FILES [TRANSCRIPT_FILES ...]
                            Paths of the transcript gff files.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the genome annotation gff files.
      --tss_files TSS_FILES [TSS_FILES ...], -t TSS_FILES [TSS_FILES ...]
                            Paths of TSS gff files.
      --srna_files SRNA_FILES [SRNA_FILES ...], -s SRNA_FILES [SRNA_FILES ...]
                            Paths of sRNA gff files for comparing sORF and sRNA to
                            detect the overlapping.
      --tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...], -tl TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]
                            TEX+/- wig files. The format is:
                            wig_file_path:TEX+/-(tex or notex):condition_id(intege
                            r):replicate_id(alphabet):strand(+ or -). If multiple
                            wig files need to be assigned, please use spaces to
                            separate the wig files. For example,
                            my_lib_tex_forward.wig:tex:1:a:+
                            my_lib_tex_reverse.wig:tex:1:a:-.
      --frag_libs FRAG_LIBS [FRAG_LIBS ...], -fl FRAG_LIBS [FRAG_LIBS ...]
                            Wig files of RNA-Seq with transcript fragmented. The
                            format is: wig_file_path:frag:condition_id(integer):re
                            plicate_id(alphabet):strand(+ or -). If multiple wig
                            files need to be assigned, please use spaces to
                            separate the wig files. For example,
                            my_lib_frag_forward.wig:frag:1:a:+
                            my_lib_frag_reverse.wig:frag:1:a:-.
      --tex_notex TEX_NOTEX, -te TEX_NOTEX
                            If the TEX+/- libraries is provided, this value is
                            that a sORF should be detected in both (TEX+ and TEX-)
                            or can be detected in only one library (TEX+ or TEX-).
                            Please assign 1 or 2. Default is 2.
      --replicate_tex REPLICATE_TEX [REPLICATE_TEX ...], -rt REPLICATE_TEX [REPLICATE_TEX ...]
                            This value is the minimal number of replicates that a
                            TSS has to be detected in. The format is
                            $NUMBERofCONDITION_$NUMBERofREPLICATE. If different
                            --replicate_tex values need to be assigned to
                            different conditions, please use spaces to separate
                            them. For example, 1_2 2_2 3_3 means that
                            --replicate_tex is 2 in number 1 and number 2
                            conditions. In number 3 condition, --replicate_tex is
                            3. For assigning the same --replicate_tex to all
                            conditions, just use like all_1 (--replicate_tex is 1
                            in all conditions). Default is all_1.
      --replicate_frag REPLICATE_FRAG [REPLICATE_FRAG ...], -rf REPLICATE_FRAG [REPLICATE_FRAG ...]
                            It is similar to --replicates_tex. This value is for
                            fragmented (or conventional) libraries.
    
    additional arguments:
      --utr_length UTR_LENGTH, -ul UTR_LENGTH
                            The utr length for comparing TSS with sORF. The
                            default number is 300.
      --min_length MIN_LENGTH, -lm MIN_LENGTH
                            The minimum residue length of sORF. Default is 30.
      --max_length MAX_LENGTH, -lM MAX_LENGTH
                            The maximum residue length of sORF. Default is 150.
      --cutoff_intergenic_coverage CUTOFF_INTERGENIC_COVERAGE, -ci CUTOFF_INTERGENIC_COVERAGE
                            The minimum coverage of intergenic sORF candidates.
      --cutoff_antisense_coverage CUTOFF_ANTISENSE_COVERAGE, -ai CUTOFF_ANTISENSE_COVERAGE
                            The minimum coverage of antisense sORF candidates.
      --cutoff_5utr_coverage CUTOFF_5UTR_COVERAGE, -cu5 CUTOFF_5UTR_COVERAGE
                            The minimum coverage of 5'UTR derived sORF candidates.
                            This value can be assigned by percentage ("p") or the
                            amount of reads ("n"). For example, p_0.05 means that
                            the coverage of sORF candidates should be higher than
                            5 percentile of all 5'UTR transcripts. n_10 means that
                            the coverage of sORF candidates should be higher than
                            10. Default is p_0.5.
      --cutoff_3utr_coverage CUTOFF_3UTR_COVERAGE, -cu3 CUTOFF_3UTR_COVERAGE
                            It is similar to --cutoff_5utr_coverage. This value is
                            for 3'UTR. Default is p_0.5.
      --cutoff_intercds_coverage CUTOFF_INTERCDS_COVERAGE, -cuf CUTOFF_INTERCDS_COVERAGE
                            It is similar to --cutoff_5utr_coverage. This value is
                            for interCDS. Default is p_0.5.
      --cutoff_base_coverage CUTOFF_BASE_COVERAGE, -cub CUTOFF_BASE_COVERAGE
                            The general minimum coverage of all sORF candidates.
                            All candidates should fit this condition as well.
                            Default is 10.
      --table_best, -tb     The output table of sORF candidates only includes
                            information of the highest expressed library. Default
                            is False.
      --start_codon START_CODON [START_CODON ...], -ac START_CODON [START_CODON ...]
                            The types of start coden. If multiple types of start
                            codon need to be assigned, please use spaces to
                            separate them. Default is ATG.
      --stop_codon STOP_CODON [STOP_CODON ...], -oc STOP_CODON [STOP_CODON ...]
                            The types of stop codon. If multiple types of stop
                            codon need to be assigned, please use spaces to
                            separate them. Default is TTA TAG TGA.
      --min_rbs_distance MIN_RBS_DISTANCE, -mr MIN_RBS_DISTANCE
                            The minimum distance (nucleotides) between the
                            ribosome binding site (Shine-Dalgarno sequence) and
                            the start codon. Default is 3.
      --max_rbs_distance MAX_RBS_DISTANCE, -Mr MAX_RBS_DISTANCE
                            The maximum distance (nucleotides) between the
                            ribosome binding site (Shine-Dalgarno sequence) and
                            the start codon. Default is 15.
      --rbs_not_after_tss, -at
                            Including the sORFs which is not associated with
                            ribosome binding site to highly-confidence sORF list.
                            Default is False.
      --fuzzy_rbs FUZZY_RBS, -zr FUZZY_RBS
                            The number of nucleotides of ribosome binding site
                            allow to be different with AGGAGG. Default is 2.
      --print_all_combination, -pa
                            For printing all combinations of multiple start and
                            stop codons. Default is False.
      --best_no_srna, -bs   Excluding the sORFs which overlap with sRNAs to highly
                            confidence sORF list. Default is False.
      --best_no_tss, -bt    Excluding the sORFs which do not start with TSS to
                            highly confidence sORF list. Default is False.
      --ignore_hypothetical_protein IGNORE_HYPOTHETICAL_PROTEIN, -ih IGNORE_HYPOTHETICAL_PROTEIN
                            For ignoring hypothetical protein in genome annotation
                            file. Default is False.

- **Output files**

Output files are stored in ``$ANNOgesic/output/sORFs``.

**statistics**: Stores statistic files.

**tables:** Stores tables of the sORFs with more details. There are also some sub-folders:

        **best_candidates:** Stores the best results of sORFs after filtering.

        **all_candidates:** Stores all candidates without filtering.

The meaning of each column is as following:

	**Genome:** Genome name.

	**Name:** the sORF name which is also shown in gff file.

	**Start:** Starting point of this sORF.

	**End:** End point of this sORF.

	**Strand:** Strand of this sORF.

	**Type:** sORF type.

	**TSS:** TSSs which are associated with this sORF.

	**Ribosome_binding_site:** Ribosome binding site (Shine-Dalgarno sequence) of this sORF.

	**All_start_points:** Positions of all start codons which can be found in the region of this sORF.

	**All_stop_points:** Positions of all stop codons which can be found in the region of this sORF.

	**sRNA_conflict:** If this sORF overlaps sRNA, the overlapped sRNA will be showed here.

	**Frame_shift:** If there are sORF candidates which can be found by frame shift, 
	the number of frame shift will be showed here. "1" means there 
	are some candidates can be found by frame shift once. "2" means there are some candidates can be found by frame shift twice.

	**Lib_type:** This sORF can be detected in TEX+/- or fragmented/conventional libraries.

	**Best_avg_coverage:** Based on coverage of all libraries, The best average coverage of this sORF will be showed here.

        **Best_highest_coverage:** Based on coverage of all libraries, The highest average coverage of this sORF will be showed here.

        **Best_lowest_coverage:** Based on coverage of all libraries, The lowest average coverage of this sORF will be showed here.

        **Track_detail:** Shows the coverage information of the libraries about this sORF. "high" means the highest coverage of the libraries,
        "low" means the lowest coverage of the libraries, and "avg" represents the average coverage of this sORF.

	**Seq:** Sequence of this sORF.

**gffs:** Stores gff files of the sORFs. There are also some sub-folders:

        **best_candidates:** Stores the best results of sORFs after filtering.

        **all_candidates:** Stores all candidates without filtering.

Some useful information can be found in the tags of the attributes within the sORF gff file.
Based on this information, we can know the details of the specific sORF. The tags are as following:

	**start_TSS:** Shows this sORF starts with which TSS.

	**with_TSS:** Which TSSs are associated with this sORF.

	**sORF_type:** Type of the sORF (5'UTR, 3'UTR, interCDS, intergenic, antisense or within CDS).

	**sRNA:** Which sRNAs are overlap with this sORF.

	**rbs:** Ribosome binding sites (Shine-Dalgarno sequences) of this sORF.

	**frame_shift:** The number of frame shifts in the regions.

.. _promoter:

promoter (Promoter motif detection)
-----------------------------------

``promoter`` can scan the upstream of TSSs to discover the promoter motifs.
We integrated `MEME <http://meme-suite.org/tools/meme>`_ (for ungapped motifs) and 
`GLAM2 <http://meme-suite.org/tools/glam2>`_ (for gapped motifs) to predict the promoters.
Based on the tool, HTML files can be generated for visualization. If the input TSS gff file 
is not computed by ANNOgesic, please add ``--tss_source`` to classify TSSs for predicting  
promoter motifs.

- **Required tools**

`MEME <http://meme-suite.org/tools/meme>`_.

`GLAM2 <http://meme-suite.org/tools/glam2>`_

`MPICH <https://http://www.mpich.org/>`_ (if parallel runs are required)

- **Required files**

**Fasta files of the genome sequences**

**Gff files of the TSSs:** If the input TSS gff file is not generated by ANNOgesic, 
the libraries and wiggle files are necessary.
Please refer to the :ref:`The input format of libraries for running ANNOgesic` 
in order to assign the correct format.

- **Arguments**

::

    usage: annogesic promoter [-h] [--program {meme,glam2,both}] --project_path
                              PROJECT_PATH --fasta_files FASTA_FILES
                              [FASTA_FILES ...] --tss_files TSS_FILES
                              [TSS_FILES ...]
                              [--motif_width MOTIF_WIDTH [MOTIF_WIDTH ...]]
                              [--num_motifs NUM_MOTIFS]
                              [--nt_before_tss NT_BEFORE_TSS]
                              [--meme_path MEME_PATH] [--glam2_path GLAM2_PATH]
                              [--e_value E_VALUE] [--end_run END_RUN]
                              [--parallels PARALLELS] [--tss_source]
                              [--tex_libs TEX_LIBS [TEX_LIBS ...]]
                              [--annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...]]
                              [--combine_all]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --program {meme,glam2,both}, -p {meme,glam2,both}
                            Please assign the program -- meme, glam2 or both.
                            Default is both
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --fasta_files FASTA_FILES [FASTA_FILES ...], -f FASTA_FILES [FASTA_FILES ...]
                            Paths of genome fasta files.
      --tss_files TSS_FILES [TSS_FILES ...], -t TSS_FILES [TSS_FILES ...]
                            Paths of TSS gff files.
      --motif_width MOTIF_WIDTH [MOTIF_WIDTH ...], -w MOTIF_WIDTH [MOTIF_WIDTH ...]
                            The length for motif detection. For detecting a range
                            of length, please insert "-" between two values.
                            Moreover, if multiple lengths need to be assigned,
                            please use spaces to separate them. For an example, 50
                            2-10 means that the lengths of motifs are 50 and
                            within 2 to 10. The number should be less or equal
                            than --nt_before_TSS. Default is 50.
      --num_motifs NUM_MOTIFS, -n NUM_MOTIFS
                            The number of motifs. Default is 10.
      --nt_before_tss NT_BEFORE_TSS, -b NT_BEFORE_TSS
                            The number of upstream nucleotides of TSS for promoter
                            prediction. Default is 50.
    
    additional arguments:
      --meme_path MEME_PATH
                            path of MEME.
      --glam2_path GLAM2_PATH
                            path of GLAM2.
      --e_value E_VALUE, -e E_VALUE
                            The maximum e value for running MEME. Default is 0.05.
      --end_run END_RUN, -er END_RUN
                            If the result of GLAM2 is not improved after running
                            this number of iteration, it will be ended. Default is
                            10000.
      --parallels PARALLELS, -pl PARALLELS
                            The number of parallel runs.
      --tss_source, -s      The TSS gff files are generated from ANNOgesic or not.
                            Default is True (from ANNOgesic)
      --tex_libs TEX_LIBS [TEX_LIBS ...], -tl TEX_LIBS [TEX_LIBS ...]
                            If --tss_source is False, please assign the name of
                            TEX+/- library. The format is:
                            wig_file_path:TEX+/-(tex or notex):condition_id(intege
                            r):replicate_id(alphabet):strand(+ or -). If multiple
                            wig files need to be assigned, please use spaces to
                            separate the wig files. For an example,
                            $WIG_PATH_1:tex:1:a:+ $WIG_PATH_2:tex:1:a:-.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            If --tss_source is False, please assign the paths of
                            genome annotation gff files.
      --combine_all, -c     For generating global promoter motifs across all
                            genomes in TSS files. Default is False.

- **Output files**

Output files are stored in ``$ANNOgesic/output/promoters``. The output folders are following:

**allfasta:** If ``--combine_all`` is True, it will combine all TSS files in ``--tss_files`` 
to generate promoter motifs. The results will be stored in this folder.

**fasta_classes:** The fasta files of different TSS types.

**$GENOME_NAME:** Stores output of MEME and GLAM2 based on different TSS types. Two sub-folders are under this folder.

	**MEME**: Stores the results of MEME.

	**GLAM2**: Stores the results of GLAM2.

	The format of folders which under these two folders is ``promoter_motifs_$FILENAME_$GENOME_$TSSTYPE_$PROMOTERLEGNTH``,
	ex: ``promoter_motifs_NC_000915.1_allgenome_internal_45_nt``.
	"NC_000915.1", "allgenome", "primary" and "45_nt" are gff filename, genome name, TSS type and upstream nucleotides of TSS, respectively.
	If genome name is "allgenome", this means that the result is generated by the information of all genomes of gff files. 
	If there is only one genome in the gff file, the genome name will be assigned as "allgenome" as well. Several files are stored in the sub-folder:
	
		**Figures of the promoter motifs:** Contains EPS and PNG files.
	
		**Details of the promoter motifs:** Contains HTML file, XML file and TXT file. These files include the TSS information.
	
		**Promoter tables:** ``meme.csv`` or ``glam2.csv`` is the promoter table which also includes the TSS information. 
		Moreover, it can used as an input for sRNA detection (``srna``). Please check the section ``srna``.

**TSS_classes:** If the TSSs are not computed by ANNOgesic, ``TSS_classes`` will be generated for classification of TSS.
TSS gff files with TSS types will be stored here.

.. _operon:

operon (Operon detection)
-------------------------

``operon`` will group TSSs, genes/CDSs/tRNAs/rRNAs, transcripts, terminators and UTRs to operons and 
sub-operons.

- **Required files**

**Gff files of the TSSs**

**Gff files of the genome annotations**

**Gff files of the transcripts**

**Gff files of the 5'UTRs**

**Gff files of the 3'UTRs**

- **Optional input files**

**Gff files of the terminators**

- **Arguments**

::

    usage: annogesic operon [-h] --project_path PROJECT_PATH --tss_files TSS_FILES
                            [TSS_FILES ...] --annotation_files ANNOTATION_FILES
                            [ANNOTATION_FILES ...] --transcript_files
                            TRANSCRIPT_FILES [TRANSCRIPT_FILES ...] --utr5_files
                            UTR5_FILES [UTR5_FILES ...] --utr3_files UTR3_FILES
                            [UTR3_FILES ...]
                            [--terminator_files TERMINATOR_FILES [TERMINATOR_FILES ...]]
                            [--tss_fuzzy TSS_FUZZY]
                            [--terminator_fuzzy TERMINATOR_FUZZY]
                            [--min_length MIN_LENGTH]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --tss_files TSS_FILES [TSS_FILES ...], -t TSS_FILES [TSS_FILES ...]
                            Paths of the TSS gff files.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the genome annotation gff files.
      --transcript_files TRANSCRIPT_FILES [TRANSCRIPT_FILES ...], -a TRANSCRIPT_FILES [TRANSCRIPT_FILES ...]
                            Paths of the transcript gff files.
      --utr5_files UTR5_FILES [UTR5_FILES ...], -u5 UTR5_FILES [UTR5_FILES ...]
                            Paths of the 5'UTR gff files.
      --utr3_files UTR3_FILES [UTR3_FILES ...], -u3 UTR3_FILES [UTR3_FILES ...]
                            Paths of the 3'UTR gff files.
      --terminator_files TERMINATOR_FILES [TERMINATOR_FILES ...], -e TERMINATOR_FILES [TERMINATOR_FILES ...]
                            Paths of terminator gff files.
    
    additional arguments:
      --tss_fuzzy TSS_FUZZY, -tf TSS_FUZZY
                            The fuzzy value for comparing between TSS and
                            transcript. Default is 5.
      --terminator_fuzzy TERMINATOR_FUZZY, -ef TERMINATOR_FUZZY
                            The fuzzy value for comparing between terminator and
                            transcript. Default is 30.
      --min_length MIN_LENGTH, -l MIN_LENGTH
                            The minimum length of operon. Default is 20.

- **Output files**

Output files are stored in ``$ANNOgesic/output/operons``. The output folders are as following:

**gffs:** Stores gff files which are integrated the information of TSSs, annotations, 
transcripts, 5'UTRs, and 3'UTRs and assign parent transcript to all features (presented by 
**Parent** in attributes of gff files).

**tables:** The tables of operons which store all information of operons and sub-operons.

The meaning of each column is as following:

	**Operon_ID:** Operon ID.

	**Genome:** Genome name.

	**Operon_position:** Starting point and end point of the operon.

	**Strand:** Strand of the operon.

	**Number_of_suboperon:** The amount of sub-operons in this operon region.

	**Position_of_suboperon:** Starting point and end point of the sub-operons.

	**Start_with_TSS:** This operon starts with TSS or not.

	**Number_of_TSS:** The number of the TSSs which are located on this operon.

	**Terminated_with_terminator:** This operon ends with TSS or not.

	**Number_of_terminator:** The number of the terminators which are associated with this operon.

	**Number_of_gene_associated_suboperon:** The number of the genes which are associated with the sub-operon.

	**Number_of_gene_associated_operon:** The number of the genes which are associated with the operon.

	**Associated_genes_with_suboperon:** Locus tag of the genes which are associated with the sub-operon.

	**Associated_genes_with_whole_operon:** Locus tag of the genes which are associated with the operon.

**statistics:** Stores statistic file which includes the number of sub-operons, monocistronic operon, polycistronic operon, etc.

.. _circrna:

circrna (circular RNA detection)
--------------------------------

``circrna`` can detect the potential circular RNAs via `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_. 
Moreover, the false positive can be removed by checking genome annotation files and quality of splicing site detection. 
The user can assign reads for detecting circular RNAs or assign BAM files to skip mapping.
BE CAREFUL, BAM files must be mapped by `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_ 
with ``--splits`` or ``circrna`` can't find the proper candidates.

- **Required tools**

`segemehl.x and testrealign.x in Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_.
For generating testrealign.x, please refer to :ref:`Required tools or databases`.

- **Required files**

**Fasta files of reads or BAM files:** If you want to use BAM files directly, they should be 
mapped by `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_ with ``--splits``.
The input format is ``$SET_NAME:$READ1,$READ2,...`` or ``$SET_NAME:$BAM1,$BAM2,...``. 
For an example, ``set1:read1.fa,read2.fa`` means these two fasta files need to be computed together.
If your BAM files are generated by mapping reads on multiple reference genomes, 
`testrealign.x <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_ may not be able to 
handle them.

**Fasta files of the genome annotations**

**Gff files of the genome annotations**

- **Arguments**

::

    usage: annogesic circrna [-h] --project_path PROJECT_PATH
                             [--read_files READ_FILES [READ_FILES ...]]
                             [--bam_files BAM_FILES [BAM_FILES ...]] --fasta_files
                             FASTA_FILES [FASTA_FILES ...] --annotation_files
                             ANNOTATION_FILES [ANNOTATION_FILES ...]
                             [--segemehl_path SEGEMEHL_PATH]
                             [--testrealign_path TESTREALIGN_PATH]
                             [--samtools_path SAMTOOLS_PATH]
                             [--parallels PARALLELS]
                             [--support_reads SUPPORT_READS]
                             [--start_ratio START_RATIO] [--end_ratio END_RATIO]
                             [--ignore_hypothetical_proteins IGNORE_HYPOTHETICAL_PROTEINS]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --read_files READ_FILES [READ_FILES ...], -rp READ_FILES [READ_FILES ...]
                            Paths of read fasta files. ANNOgesic will map the
                            reads via segemehl (with -S). Required format:
                            $SET_NAME:$READ1,$READ2,... If multiple data sets need
                            to be assigned, please separated them by spaces. For
                            using BAM files, please check --bam_files.
      --bam_files BAM_FILES [BAM_FILES ...], -b BAM_FILES [BAM_FILES ...]
                            Path of input BAM files. Required format:
                            $SET_NAME:$BAM1,$BAM2,... . BAM files need to be
                            generated using the mapper segemehl with the parameter
                            "-S". If multiple data sets need to be assigned,
                            please separated them by spaces. For using fasta files
                            of reads, please check --read_files.
      --fasta_files FASTA_FILES [FASTA_FILES ...], -f FASTA_FILES [FASTA_FILES ...]
                            Paths of the genome fasta files.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the genome annotation gff files.
    
    additional arguments:
      --segemehl_path SEGEMEHL_PATH
                            Path of segemehl.x in segemehl package.
      --testrealign_path TESTREALIGN_PATH
                            Path of testrealign.x in segemehl package.
      --samtools_path SAMTOOLS_PATH
                            Path of samtools.
      --parallels PARALLELS, -p PARALLELS
                            The number of parallels runs for mapping if
                            --read_files is assigned. Default is 10.
      --support_reads SUPPORT_READS, -s SUPPORT_READS
                            The minimum reads for supporting a circular RNA.
                            Default is 10.
      --start_ratio START_RATIO, -sr START_RATIO
                            The minimum ratio -- (reads of circRNA / all reads) at
                            starting point of candidate. Default is 0.5.
      --end_ratio END_RATIO, -er END_RATIO
                            I is similar to --start_ratio. This value is for end
                            point of candidate. Default is 0.5.
      --ignore_hypothetical_proteins IGNORE_HYPOTHETICAL_PROTEINS, -ih IGNORE_HYPOTHETICAL_PROTEINS
                            For ignoring hypothetical protein in genome annotation
                            file. Default is False.

- **Output files**

Output files are stored in ``$ANNOgesic/output/circRNAs``. The output folders are following:

**segemehl_alignment_files:** If read files are assigned, the folder is for results of mapping.

**segemehl_splice_results:** The results of splicing detection. For understanding the splicing tables, please 
refer to `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_.

**gffs:** Stores gff files of the circular RNAs. ``$GENOME_$SET_circRNA_best.gff`` is gff files for the best results after checking 
genome annotation and quality of splicing. ``$GENOME_$SET_circRNA_all.gff`` is for all candidates without filtering.
Some useful information can be found in the tags of the attributes within the circular RNA gff file.
Based on this information, we can know the details of the specific circular RNA. The tags are as following:

	**support_reads:** The number of reads which support the circular RNA.

	**read_at_start:** (supported reads / total reads) at starting point of the circular RNA.

	**read_at_end:** (supported reads / total reads) at end point of the circular RNA.

	**conflict:** The circular RNA overlap genome annotation or not.

	**method:** The circular RNA is detected by which method.

**circRNA_tables:** Stores tables of the circular RNAs with more details. The meaning of each column is as following:

	**ID:** ID of this circular RNA.

	**Genome:** Genome name.

	**Strand:** Strand of the circular RNA.

	**Start:** Starting point of the circular RNA.

	**End:** End point of the circular RNA.

	**Annotation_overlap:** If there is genome annotation which overlap this circular RNA, the overlapped feature will be showed here.

	**Supported_reads:** The number of reads which support the circular RNA.

	**Supported_reads/Reads_at_start:** (supported reads / total reads) at starting point of the circular RNA.

	**Supported_reads/Reads_at_end:** (supported reads / total reads) at end point of the circular RNA.	

.. _go_term:

go_term (GO term retrieving)
----------------------------

``go_term`` can retrieve the information of Gene Ontology from Uniprot.
Some analyses of GO terms can be done as well.

- **Required files**

**Uniprot mapping table:** `idmapping_selected.tab from Uniprot <http://www.uniprot.org/downloads>`_.

**GOslim file:** `goslim.obo <http://geneontology.org/page/go-slim-and-subset-guide>`_.

**GO file:** `go.obo <http://geneontology.org/page/download-ontology>`_.

**Gff files of the genome annotations**

- **Optional input files**

**Gff files of the transcripts:** For detecting the GO terms only based on expressed CDSs.

- **Arguments**

::

    usage: annogesic go_term [-h] --project_path PROJECT_PATH --annotation_files
                             ANNOTATION_FILES [ANNOTATION_FILES ...]
                             [--transcript_files TRANSCRIPT_FILES [TRANSCRIPT_FILES ...]]
                             --uniprot_id UNIPROT_ID --go_obo GO_OBO --goslim_obo
                             GOSLIM_OBO
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the genome annotation gff files.
      --transcript_files TRANSCRIPT_FILES [TRANSCRIPT_FILES ...], -a TRANSCRIPT_FILES [TRANSCRIPT_FILES ...]
                            Paths of the transcript gff files which can be used to
                            retrieve GO terms based on expressed CDS and all CDS.
      --uniprot_id UNIPROT_ID, -u UNIPROT_ID
                            Path of the UniProt ID mapping database.
      --go_obo GO_OBO, -go GO_OBO
                            Path of go.obo.
      --goslim_obo GOSLIM_OBO, -gs GOSLIM_OBO
                            Path of goslim.obo.

- **Output files**

Output files are stored in ``ANNOgesic/output/GO_terms``. If gff files of the transcript are assigned, 
two sub-folders will be generated. Results of the expressed genes are stored in ``ANNOgesic/output/GO_term/expressed_CDSs`` and 
results of all CDSs are stored in ``ANNOgesic/output/GO_term/all_CDSs``.

**GO_term_results:** Stores tables of the GO terms information. The meaning of each column is as following:

	**Genome:** Genome name.

	**Strand:** Strand of the gene.

	**Start:** Starting point of this CDS.

	**End:** End point of this CDS.

	**Protein_id** Protein ID of this CDS.

	**GO_term** GO term of This CDS.

**statistics:** Stores statistic files and figures.

	**GO term with corresponding amount:** Format of the filename is ``stat_$GENOME.csv``.

	**Figures of the three GO term classes:** All figures are stored in the sub-folder - ``figs``. 
	``$GENOME_biological_process.png``, ``$GENOME_cellular_component.png``, 
	``$GENOME_molecular_function.png`` and ``$GENOME_three_roots.png`` are the figures of 
	biological process, cellular component, molecular_function and three roots of GO term classes, respectively.

.. _srna_target:

srna_target (sRNA target prediction)
------------------------------------

``srna_target`` can search potential targets of the sRNA via different programs 
(RNAup or RNAplex or both). We recommend running with both 
programs. ``srna_target`` can also compare the results of both programs and provide the best ones.

- **Required tools**

`ViennaRNA <http://www.tbi.univie.ac.at/RNA/>`_ .

- **Required files**

**Gff files of the genome annotations**

**Gff files of the sRNAs**

**Fasta files of the genomes**

- Arguments

::

    usage: annogesic srna_target [-h] --project_path PROJECT_PATH
                                 --annotation_files ANNOTATION_FILES
                                 [ANNOTATION_FILES ...] --fasta_files FASTA_FILES
                                 [FASTA_FILES ...] --srna_files SRNA_FILES
                                 [SRNA_FILES ...]
                                 [--query_srnas QUERY_SRNAS [QUERY_SRNAS ...]]
                                 [--program {RNAplex,RNAup,both}] [--top TOP]
                                 [--rnaplfold_path RNAPLFOLD_PATH]
                                 [--rnaplex_path RNAPLEX_PATH]
                                 [--rnaup_path RNAUP_PATH]
                                 [--interaction_length INTERACTION_LENGTH]
                                 [--window_size_target WINDOW_SIZE_TARGET]
                                 [--span_target SPAN_TARGET]
                                 [--window_size_srna WINDOW_SIZE_SRNA]
                                 [--span_srna SPAN_SRNA]
                                 [--unstructured_region_rnaplex_target UNSTRUCTURED_REGION_RNAPLEX_TARGET]
                                 [--unstructured_region_rnaplex_srna UNSTRUCTURED_REGION_RNAPLEX_SRNA]
                                 [--unstructured_region_rnaup UNSTRUCTURED_REGION_RNAUP]
                                 [--energy_threshold ENERGY_THRESHOLD]
                                 [--duplex_distance DUPLEX_DISTANCE]
                                 [--parallels_rnaplex PARALLELS_RNAPLEX]
                                 [--parallels_rnaup PARALLELS_RNAUP]
                                 [--continue_rnaup]
                                 [--potential_target_start POTENTIAL_TARGET_START]
                                 [--potential_target_end POTENTIAL_TARGET_END]
                                 [--target_feature TARGET_FEATURE [TARGET_FEATURE ...]]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the genome annotation gff files.
      --fasta_files FASTA_FILES [FASTA_FILES ...], -f FASTA_FILES [FASTA_FILES ...]
                            Paths of the genome fasta files.
      --srna_files SRNA_FILES [SRNA_FILES ...], -s SRNA_FILES [SRNA_FILES ...]
                            Paths of the sRNA gff files.
      --query_srnas QUERY_SRNAS [QUERY_SRNAS ...], -q QUERY_SRNAS [QUERY_SRNAS ...]
                            The query sRNA. The input format is
                            $GENOME:$START:$END:$STRAND. If multiple sRNAs need to
                            be assigned, please use spaces to separate them. For
                            example, NC_007795.1:200:534:+
                            NC_007795.1:6767:6900:-. If all sRNAs of the reference
                            genome need to be used, please assign "all". Default
                            is all.
      --program {RNAplex,RNAup,both}, -p {RNAplex,RNAup,both}
                            The program for detecting sRNA-mRNA interaction.
                            Please assign "RNAplex" or "RNAup" or "both". Default
                            is both.
      --top TOP, -t TOP     The ranking number of targets which will be included
                            to final output. The ranking is based on the binding
                            energy. Default is 20.
    
    additional arguments:
      --rnaplfold_path RNAPLFOLD_PATH
                            Path of RNAplfold in Vienna package.
      --rnaplex_path RNAPLEX_PATH
                            Path of RNAplex in Vienna package.
      --rnaup_path RNAUP_PATH
                            Path of RNAup in Vienna package.
      --interaction_length INTERACTION_LENGTH, -i INTERACTION_LENGTH
                            Maximum length of an interaction. Default is 30.
      --window_size_target WINDOW_SIZE_TARGET, -wt WINDOW_SIZE_TARGET
                            The average of the pair probabilities over windows for
                            mRNA target. It only works for "RNAplex". Default is
                            240.
      --span_target SPAN_TARGET, -st SPAN_TARGET
                            The maximum allowed separation of a base pair to span
                            for mRNA target. It only works for "RNAplex". Default
                            is 160.
      --window_size_srna WINDOW_SIZE_SRNA, -ws WINDOW_SIZE_SRNA
                            It is similar to --window_size_target, but for sRNA.
                            Default is 30.
      --span_srna SPAN_SRNA, -ss SPAN_SRNA
                            It is similar to --span_target, but for sRNA. Default
                            is 30.
      --unstructured_region_rnaplex_target UNSTRUCTURED_REGION_RNAPLEX_TARGET, -ut UNSTRUCTURED_REGION_RNAPLEX_TARGET
                            Calculate the mean probability of the unpaired region
                            for mRNA target. It only works for "RNAplex". Default
                            is 30.
      --unstructured_region_rnaplex_srna UNSTRUCTURED_REGION_RNAPLEX_SRNA, -us UNSTRUCTURED_REGION_RNAPLEX_SRNA
                            It is similar to --unstructured_region_rnaplex_target,
                            but for sRNA. Default is 30.
      --unstructured_region_rnaup UNSTRUCTURED_REGION_RNAUP, -uu UNSTRUCTURED_REGION_RNAUP
                            Compute the mean probability of unpaired region. It
                            only works for "RNAup". Default is 30.
      --energy_threshold ENERGY_THRESHOLD, -e ENERGY_THRESHOLD
                            The minimum energy for a duplex. It only works for
                            "RNAplex". Default is -8.
      --duplex_distance DUPLEX_DISTANCE, -d DUPLEX_DISTANCE
                            Distance between target 3'ends of two consecutive
                            duplexes. It works for "RNAplex". Default is 20.
      --parallels_rnaplex PARALLELS_RNAPLEX, -pp PARALLELS_RNAPLEX
                            The number of parallel runs for running RNAplex.
                            Default is 5.
      --parallels_rnaup PARALLELS_RNAUP, -pu PARALLELS_RNAUP
                            The number of parallel runs for running RNAup. Default
                            is 20.
      --continue_rnaup, -cr
                            For running RNAup based on the previous intermediate
                            results if the previous process was crushed. Default
                            is False.
      --potential_target_start POTENTIAL_TARGET_START, -ps POTENTIAL_TARGET_START
                            Extracting the upstream nucleotides of
                            --target_feature. Default is 200.
      --potential_target_end POTENTIAL_TARGET_END, -pe POTENTIAL_TARGET_END
                            Extracting the nucleotides directly behind staring
                            point of --target_feature. Default is 150.
      --target_feature TARGET_FEATURE [TARGET_FEATURE ...], -tf TARGET_FEATURE [TARGET_FEATURE ...]
                            The feature name of potential targets. If multiple
                            features need to be assigned, please use spaces to
                            separate them. For example, CDS exon. Default is CDS.

- **Output files**

Output files are stored in ``$ANNOgesic/output/sRNA_targets``.

**RNAplex_results:** Stores all results of RNAplex. ``$GENOME_RNAplex.txt`` is raw results of RNAplex.
``$GENOME_RNAplex_rank.csv`` is the tables with details, and the targets are sorted by binding energy. 
The meaning of each column in ``$GENOME_RNAplex_rank.csv`` is as following:

	**sRNA:** sRNA name which is shown in sRNA gff file.

	**Genome:** Genome name.

	**sRNA_position:** Starting point and end point of this sRNA.

	**sRNA_interacted_position_RNAplex:** The interaction region of this sRNA.

	**sRNA_strand:** Strand of this sRNA.

	**Target:** Locus tag or gene name of the target mRNA.

	**Target_position:** Starting point and end point of this mRNA.

	**Target_interacted_position_RNAplex:** The interaction region of this mRNA.

	**Target_strand:** Strand of this target mRNA.

	**Energy_RNAplex:** Interaction energy change of this interaction.

	**Rank_RNAplex:** Ranking of the interaction (the ranking is based on the binding energy).

**RNAup_results:** Stored all results of RNAup. ``$GENOME_RNAup.txt`` is raw results of RNAup.
``$GENOME_RNAup_rank.csv`` is the tables with details, and the targets are 
sorted by binding energy. The meaning of each column is similar to the table of RNAplex.

**merged_results:** Store the results which are merged by the results of ``RNAplex_results`` and ``RNAup_results``. 
``$GENOME_merge.csv`` contains all candidates of the both programs. 
``$GENOME_overlap.csv`` contains the results which are top 20 (default) in the both methods. 
The meaning of each column is similar to the table of RNAplex.

**sRNA_seqs:** Stores fasta sequences of the sRNAs.

**target_seqs:** Stores fasta sequences of the potential targets.

.. _ppi_network:

ppi_network (protein-protein interaction network detection)
-----------------------------------------------------------

``ppi_network`` can retrieve the PPI data from `STRING <http://string-db.org/>`_. 
Then using `PIE <http://www.ncbi.nlm.nih.gov/CBBresearch/Wilbur/IRET/PIE/>`_ to search 
the supported literature of the protein-protein interaction networks. 

- **Required files**

**Species table of STRING:** `species.v${VERSION}.txt from STRING <http://string-db.org/cgi/download.pl>`_.

**Gff files of the genome annotations**

- **Arguments**

::

    usage: annogesic ppi_network [-h] --project_path PROJECT_PATH
                                 --annotation_files ANNOTATION_FILES
                                 [ANNOTATION_FILES ...] --species_string
                                 SPECIES_STRING --query_strains QUERY_STRAINS
                                 [QUERY_STRAINS ...] [--query QUERY [QUERY ...]]
                                 [--without_strain_pubmed] [--score SCORE]
                                 [--node_size NODE_SIZE]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the genome annotation gff files with proper
                            locus_tag item in the attributes.
      --species_string SPECIES_STRING, -d SPECIES_STRING
                            Path of the species table of STRING
                            (species.$VERSION.txt).
      --query_strains QUERY_STRAINS [QUERY_STRAINS ...], -s QUERY_STRAINS [QUERY_STRAINS ...]
                            The name of input file of the query genomes. Required
                            format: $GFF_FILE:$STRAIN_IN_GFF:$STRAIN_IN_STRING:$ST
                            RAIN_FOR_PUBMED. $GFF_FILE is the name of the gff
                            file, $STRAIN_IN_GFF is the name of the strain in gff
                            file, $STRAIN_IN_STRING is the strain name in species
                            table of STRING (species.$VERSION.txt), and
                            $STRAIN_FOR_PUBMED is the strain name for searching in
                            Pubmed. If the strain is not available in STRING
                            database, it can be relaced by a related strain. For
                            example, Staphylococcus_aureus_HG003.gff:Staphylococcu
                            s_aureus_HG003:"Staphylococcus aureus NCTC
                            8325":"Staphylococcus aureus" (Staphylococcus aureus
                            NCTC 8325 is a related strain of HG003 since HG003 is
                            not available in STRING). if the assigned name is with
                            spaces, please use Double quotes. For assigning
                            multiple strains, please separated them by spaces.
      --query QUERY [QUERY ...], -q QUERY [QUERY ...]
                            The query proteins. Required format:
                            $GENOME_NAME_OF_GFF:$START_POINT:$END_POINT:$STRAND.
                            For assigning multiple proteins, please use spaces to
                            separate them. For example,
                            Staphylococcus_aureus_HG003:345:456:+
                            Staphylococcus_aureus_HG003:2000:3211:-. For computing
                            all proteins in gff files, just type "all". Default is
                            all.
      --without_strain_pubmed, -n
                            For retrieving the literature from Pubmed only based
                            on protein name without assigning strains. Default is
                            False.
    
    additional arguments:
      --score SCORE, -ps SCORE
                            The minimum PIE score for searching literature. The
                            value is from -1 (worst) to 1 (best). Default is 0.
      --node_size NODE_SIZE, -ns NODE_SIZE
                            The size of nodes in figure, default is 4000.

- **Output files**

Output files are stored in ``$ANNOgesic/output/PPI_networks``. The output folders are as following:

**best_results:** Stores the results which the scores of `PIE <http://www.ncbi.nlm.nih.gov/CBBresearch/Wilbur/IRET/PIE/>`_
for supported literature are higher than ``--score``.

**all_results:** Stores the results of all protein-protein interactions
(including the low score(`PIE <http://www.ncbi.nlm.nih.gov/CBBresearch/Wilbur/IRET/PIE/>`_) literature).

Under "best_results" and "all_results", several files and folders are generated:

	**Results of searching literatures without assigning a specific strain:** ``$STRAIN_without_strain.csv``. 

	**Results of searching literatures with assigning a specific strain:** ``$STRAIN_with_strain.csv``.
 
	**without_strain:** Stores all interaction information which is searched without assigning a specific strain. 

	**with_strain:** Stores all interaction information which is searched with assigning a specific strain. 

**figures:** Stores the protein-protein networks of the query proteins. 
Thickness represents how many literature can be found for the interactions. 
Solid line means that strong supported literature can be found. Dash-dot line 
means that the supported literature are very weak. Dot line means that no supported literature can be found. 
Color is the best score of the supported literature of the interactions.

.. _subcellular_localization:

localization (subcellular localization prediction)
--------------------------------------------------------------

``localization`` can predict the subcellular localization of proteins. 
Some statistics and visualization files are provided as well.

- **Required tools**

`Psortb <http://www.psort.org/psortb/>`_.

- **Required files**

**Gff files of the genome annotations**

**Fasta files of the genome sequences**

- **Optional input files**

**Gff files of the transcripts:** For detecting subcellular localization only based on expressed CDSs.

- **Arguments**

::

    usage: annogesic localization [-h] --project_path PROJECT_PATH
                                  --annotation_files ANNOTATION_FILES
                                  [ANNOTATION_FILES ...] --fasta_files FASTA_FILES
                                  [FASTA_FILES ...]
                                  [--transcript_files TRANSCRIPT_FILES [TRANSCRIPT_FILES ...]]
                                  --bacteria_type {positive,negative}
                                  [--psortb_path PSORTB_PATH]
                                  [--difference_multi DIFFERENCE_MULTI]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of genome annotation gff files.
      --fasta_files FASTA_FILES [FASTA_FILES ...], -f FASTA_FILES [FASTA_FILES ...]
                            Paths of genome fasta files.
      --transcript_files TRANSCRIPT_FILES [TRANSCRIPT_FILES ...], -a TRANSCRIPT_FILES [TRANSCRIPT_FILES ...]
                            Paths of the transcript gff files for detecting the
                            subcellular localization based on expressed CDS and
                            all CDS.
      --bacteria_type {positive,negative}, -b {positive,negative}
                            The type of bacteria (Gram-positive or Gram-negative).
                            Please assign 'positive' or 'negative'.
    
    additional arguments:
      --psortb_path PSORTB_PATH
                            Path of Psortb.
      --difference_multi DIFFERENCE_MULTI, -d DIFFERENCE_MULTI
                            For the protein which have multiple predicted
                            locations, if the difference of psorb scores is
                            smaller than this value, all locations will be printed
                            out. Default is 0.5. The maximum value is 10.

- **Output files**

Output files are stored in ``$ANNOgesic/output/subcellular_localization``. If gff files of the transcripts are assigned,
two sub-folders will be generated. Results of the expressed CDSs are stored in 
``ANNOgesic/output/subcellular_localization/expressed_CDSs`` and
results of all CDS are stored in ``ANNOgesic/output/subcellular_localization/all_CDSs``.

**psortb_results**: Stores the results of Psortb:

	**Raw output data of Psortb:** Format of the filename is ``$GENOME_raw.txt``.

	**Table of subcellular localization:** Format of the filename is ``$GENOME_table.csv``. 
	The meaning of each column is as following:

		**Genome:** Genome name.

		**Protein:** Protein ID.

		**Strand:** Strand of this protein.

		**Start:** Starting point of this protein.

		**End:** End point of this protein.

		**Location:** Predicted subcellular localization of this protein.

		**Score:** Psortb score.

**statistics:** Stores statistic files and figures.

	**Subcellular localization with corresponding amounts:** Format of the filename is ``stat_$GENOME_sublocal.csv``.

	**Figure of Subcellular localization with corresponding amounts:** Format of the filename is ``$FILENAME_$GENOME_sublocal.png``.

.. _riboswitch_thermometer:

riboswitch_thermometer (riboswitch and RNA thermometer detection)
-----------------------------------------------------------------

``riboswitch_thermometer`` can search riboswitches and RNA thermometers between 
TSSs（the starting point of transcript was assigned if no TSS was detected) 
and its downstream CDSs, as well as associated ribosome binding sites 
(Shine-Dalgarno sequence). Then using `Infernal <http://infernal.janelia.org/>`_ to scan 
the region in `Rfam <http://rfam.xfam.org/>`_.

- **Required tools**

`Infernal <http://infernal.janelia.org/>`_.

- **Required files**

`Rfam <http://rfam.xfam.org/>`_.

**Gff files of the genome annotations**

**Gff files of the transcripts**

**Gff files of the TSSs**

**Fasta files of the genome sequences**

**Rfam ID files of the riboswitch or RNA thermometer:** 
The file should contain Rfam IDs, 
name and description of riboswitches or RNA thermometers as following. 

======== ==== ==========================
#Rfam_ID Name Description
-------- ---- --------------------------
RF00162  SAM  SAM riboswitch box leader
RF00059  TPP  TPP riboswitch THI element
======== ==== ==========================

All columns are separated by ``tab``. You can also download 
`riboswitch and RNA thermometer data <https://github.com/Sung-Huan/ANNOgesic/tree/master/database>`_
from our Git repository.

- **Arguments**

::

    usage: annogesic riboswitch_thermometer [-h] --project_path PROJECT_PATH
                                            [--program {riboswich,thermometer,both}]
                                            [--riboswitch_id_file RIBOSWITCH_ID_FILE]
                                            [--rna_thermometer_id_file RNA_THERMOMETER_ID_FILE]
                                            --annotation_files ANNOTATION_FILES
                                            [ANNOTATION_FILES ...] --tss_files
                                            TSS_FILES [TSS_FILES ...]
                                            --transcript_files TRANSCRIPT_FILES
                                            [TRANSCRIPT_FILES ...] --fasta_files
                                            FASTA_FILES [FASTA_FILES ...]
                                            --rfam_path RFAM_PATH
                                            [--cmscan_path CMSCAN_PATH]
                                            [--cmpress_path CMPRESS_PATH]
                                            [--utr_length UTR_LENGTH]
                                            [--e_value E_VALUE] [--output_all]
                                            [--fuzzy FUZZY]
                                            [--fuzzy_rbs FUZZY_RBS]
                                            [--start_codon START_CODON [START_CODON ...]]
                                            [--max_dist_rbs MAX_DIST_RBS]
                                            [--min_dist_rbs MIN_DIST_RBS]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --program {riboswich,thermometer,both}, -p {riboswich,thermometer,both}
                            Please assign the feature for detection. The options
                            can be "riboswitch", "thermometer", "both". Default is
                            both.
      --riboswitch_id_file RIBOSWITCH_ID_FILE, -ri RIBOSWITCH_ID_FILE
                            Path of the file whic contain the information of
                            riboswitch in Rfam. Required format of the file:
                            $RFAM_ID{tab}$RIBOSWITCH_NAME{tab}$DESCRIPTION. Please
                            check an exmple in https://github.com/Sung-Huan/ANNOge
                            sic/blob/master/database/Rfam_riboswitch_ID.csv
      --rna_thermometer_id_file RNA_THERMOMETER_ID_FILE, -ti RNA_THERMOMETER_ID_FILE
                            It is similar to -riboswitch_id_file, but for RNA
                            thermometer. Please check an example in
                            https://github.com/Sung-Huan/ANNOgesic/blob/master/dat
                            abase/Rfam_RNA_thermometer_ID.csv
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the annotation gff files.
      --tss_files TSS_FILES [TSS_FILES ...], -t TSS_FILES [TSS_FILES ...]
                            Paths of the TSS gff files.
      --transcript_files TRANSCRIPT_FILES [TRANSCRIPT_FILES ...], -a TRANSCRIPT_FILES [TRANSCRIPT_FILES ...]
                            Paths of the transcript gff files.
      --fasta_files FASTA_FILES [FASTA_FILES ...], -f FASTA_FILES [FASTA_FILES ...]
                            Paths of the genome fasta files.
      --rfam_path RFAM_PATH, -R RFAM_PATH
                            Path of the Rfam CM database.
    
    additional arguments:
      --cmscan_path CMSCAN_PATH, -cs CMSCAN_PATH
                            Path of cmscan in Infernal package.
      --cmpress_path CMPRESS_PATH, -cp CMPRESS_PATH
                            Path of cmpress in Infernal package.
      --utr_length UTR_LENGTH, -u UTR_LENGTH
                            The UTR length. Default is 300.
      --e_value E_VALUE, -e E_VALUE
                            The maximum e value. Default is 0.001.
      --output_all, -o      One query sequence may fit multiple riboswitches or
                            RNA thermometers. It can print multiple riboswitches
                            or RNA thermometers. Otherwise, only the highest
                            confident one will be printed. Default is False.
      --fuzzy FUZZY, -z FUZZY
                            The fuzzy nucleotides for extracting the sequences of
                            potential riboswitches or RNA thermometers. Default is
                            10.
      --fuzzy_rbs FUZZY_RBS, -zr FUZZY_RBS
                            The number of nucleotides of ribosome binding site
                            allow to be different with AGGAGG. Default is 2.
      --start_codon START_CODON [START_CODON ...], -ac START_CODON [START_CODON ...]
                            The types of start codon. If multiple types need to be
                            assigned, please use spaces to separate them. Default
                            is ATG.
      --max_dist_rbs MAX_DIST_RBS, -Mr MAX_DIST_RBS
                            The maximum distance (nucleotides) between the
                            ribosome binding site (Shine-Dalgarno sequence) and
                            the start codon. Default is 14.
      --min_dist_rbs MIN_DIST_RBS, -mr MIN_DIST_RBS
                            The minimum distance (nucleotides) between the
                            ribosome binding site (Shine-Dalgarno sequence) and
                            the start codon. Default is 5.

- **Output files**

Output files of the riboswitches are stored in ``$ANNOgesic/output/riboswitches`` and 
output files of the RNA thermometers are stored in ``$ANNOgesic/output/RNA_thermometers``.

Names of the output folders are as following:

**scan_Rfam_results:** Stores the results of searching to Rfam with ``cmscan`` (`Infernal <http://infernal.janelia.org/>`_). 

**gffs:** Stores gff files of riboswitches/RNA_thermometers. 
Some useful information can be found in the tags of the attributes within the gff file.
Based on this information, we can know the details of the specific riboswitches/RNA_thermometers. 
The tags are as following:

	**rfam_id:** Rfam ID of this riboswitch/RNA_thermometer.

	**e-value:** E-value of searching this riboswitch/RNA_thermometer to Rfam.

	**method:** This riboswitch/RNA_thermometer is detected by which method.

**tables:** Stores tables of riboswiches/RNA_thermometers with more details. The meaning of each column in this table is as following:

	**ID:** Riboswich/RNA_thermometer ID.

	**Genome:** Genome name.

	**Strand:** Strand of the riboswich/RNA_thermometer.

	**Associated_CDS:** Downstream CDS of the riboswich/RNA_thermometer.

	**Start_genome:** This riboswich/RNA_thermometer starts from which position of the genome.

	**End_genome:** This riboswich/RNA_thermometer ends to which position of the genome.

	**Rfam:** Rfam ID of this riboswich/RNA_thermometer.

	**E_value:** E-value of searching this riboswitch/RNA_thermometer to Rfam.

	**Start_align:** Position of this riboswich/RNA_thermometer can be aligned to the genome.

	**End_align:** Position this riboswich/RNA_thermometer can be aligned to the genome.

**statistics:** Stores the file which contains the riboswich/RNA_thermometer with corresponding amount.

.. _crispr:

crispr (CRISPR detection)
-------------------------
``crispr`` integrates CRISPR Recognition Tool (`CRT <http://www.room220.com/crt/>`_) which can detect the repeat 
units and spacers of CRISPR. Moreover, the false positive can be removed by comparing candidates with genome annotation.

- **Required tools**

`CRT <http://www.room220.com/crt/>`_.

- **Required files**

**Fasta files of the genome sequences**

- **Optional input files**

**Gff files of the genome annotations:** This file can be used for removing false positive.

- **Arguments**

::

    usage: annogesic crispr [-h] --project_path PROJECT_PATH --fasta_files
                            FASTA_FILES [FASTA_FILES ...]
                            [--annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...]]
                            [--crt_path CRT_PATH] [--window_size WINDOW_SIZE]
                            [--min_number_repeats MIN_NUMBER_REPEATS]
                            [--min_length_repeat MIN_LENGTH_REPEAT]
                            [--Max_length_repeat MAX_LENGTH_REPEAT]
                            [--min_length_spacer MIN_LENGTH_SPACER]
                            [--Max_length_spacer MAX_LENGTH_SPACER]
                            [--ignore_hypothetical_protein]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --fasta_files FASTA_FILES [FASTA_FILES ...], -f FASTA_FILES [FASTA_FILES ...]
                            Paths of the genome fasta files.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the genome gff files for comparing CRISPRs
                            and the genome annotation to remove the false
                            positives. Default is None.
    
    additional arguments:
      --crt_path CRT_PATH   Path of CRT.jar. Default is /usr/local/bin/CRT.jar
      --window_size WINDOW_SIZE, -w WINDOW_SIZE
                            Length of the window size for searching CRISPR (range:
                            6-9). Default is 8.
      --min_number_repeats MIN_NUMBER_REPEATS, -mn MIN_NUMBER_REPEATS
                            Minimum number of repeats that a CRISPR must contain.
                            Default is 3.
      --min_length_repeat MIN_LENGTH_REPEAT, -ml MIN_LENGTH_REPEAT
                            Minimum length of CRISPR repeats. Default is 23.
      --Max_length_repeat MAX_LENGTH_REPEAT, -Ml MAX_LENGTH_REPEAT
                            Maximum length of CRISPR repeats. Default is 47.
      --min_length_spacer MIN_LENGTH_SPACER, -ms MIN_LENGTH_SPACER
                            Minimum length of CRISPR spacers. Default is 26.
      --Max_length_spacer MAX_LENGTH_SPACER, -Ms MAX_LENGTH_SPACER
                            Maximum length of CRISPR spacers. Default is 50.
      --ignore_hypothetical_protein, -in
                            For ignoring the hypothetical proteins. Default is
                            False.

- **Output files**

Output files are stored in ``$ANNOgesic/output/crisprs``. The folders which are generated by the subcommand are as following:

**CRT_results:** Stores the output of `CRT <http://www.room220.com/crt/>`_.

**gffs:** Stores CRSIPR gff files. Two sub-folders are under this folder:

	**all_candidates:** Stores gff files which contains all CRISPRs.
 
	**best_candidates:** Stores gff files which contains the CRISPRs without overlapping genome annotation.

**statistics:** Stores statistic files.

.. _optimize_tss_ps:

optimize_tss_ps (optimization of TSS and processing site detection)
-------------------------------------------------------------------

``optimize_tss_ps`` can adapt the parameter set of `TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_. 
For running it, please manual detect TSSs around 200kb and find at least 50 TSSs (using gff format).
If there are less than 50 TSSs within 200kb, please continue checking until 50 TSSs are detected.
Then ``optimize_tss_ps`` can scan whole genome based on the manual detected set to get optimized parameters.

- **Required tools**

`TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_.

- **Required files**

**Wiggle files of TEX +/-:** Please check the section :ref:`The input format of libraries for running ANNOgesic`.

**Fasta files of the genome sequences**

**Gff files of the genome annotations**

**Gff files of the manual-detected TSSs**

- **Arguments**

::

    usage: annogesic optimize_tss_ps [-h] --project_path PROJECT_PATH
                                     [--program {TSS,PS}] --fasta_files
                                     FASTA_FILES [FASTA_FILES ...]
                                     --annotation_files ANNOTATION_FILES
                                     [ANNOTATION_FILES ...] --manual_files
                                     MANUAL_FILES [MANUAL_FILES ...]
                                     [--curated_sequence_length CURATED_SEQUENCE_LENGTH [CURATED_SEQUENCE_LENGTH ...]]
                                     --tex_notex_libs TEX_NOTEX_LIBS
                                     [TEX_NOTEX_LIBS ...]
                                     [--replicate_tex REPLICATE_TEX [REPLICATE_TEX ...]]
                                     --condition_names CONDITION_NAMES
                                     [CONDITION_NAMES ...]
                                     [--tsspredator_path TSSPREDATOR_PATH]
                                     [--max_height MAX_HEIGHT]
                                     [--max_height_reduction MAX_HEIGHT_REDUCTION]
                                     [--max_factor MAX_FACTOR]
                                     [--max_factor_reduction MAX_FACTOR_REDUCTION]
                                     [--max_base_height MAX_BASE_HEIGHT]
                                     [--max_enrichment_factor MAX_ENRICHMENT_FACTOR]
                                     [--max_processing_factor MAX_PROCESSING_FACTOR]
                                     [--utr_length UTR_LENGTH] [--cluster CLUSTER]
                                     [--parallels PARALLELS] [--steps STEPS]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --program {TSS,PS}, -p {TSS,PS}
                            The feature for optimization. Please assign "TSS" or
                            "PS". Default is TSS.
      --fasta_files FASTA_FILES [FASTA_FILES ...], -f FASTA_FILES [FASTA_FILES ...]
                            Paths of the fasta file of the reference genome.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the gff file of the reference genome.
      --manual_files MANUAL_FILES [MANUAL_FILES ...], -m MANUAL_FILES [MANUAL_FILES ...]
                            Paths of the manual-checked TSS or PS in gff format.
                            It is used for benchmarking the prediction during the
                            optimization. The set should comprise roughly 50
                            TSS/PS or more.
      --curated_sequence_length CURATED_SEQUENCE_LENGTH [CURATED_SEQUENCE_LENGTH ...], -le CURATED_SEQUENCE_LENGTH [CURATED_SEQUENCE_LENGTH ...]
                            The length of the sequence used for the manual set of
                            TSS/PS. This value is required to calculate the
                            accurracy. If the whole genome was used write "all".
                            Otherwise use the name of the reference sequence in
                            the folowing format: $GENOME:SLENGTH. Multiple entries
                            are accepted. For an example, test.gff contains two
                            sequences s1 and s2. For s1 100 kb were checked while
                            for s2 the whole sequence was curated. The value of
                            this argument would be s1:100000 s2:all. Per default
                            all the full length of all sequences will be used.
      --tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...], -tl TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]
                            TEX+/- wig files for TSSpredator. The format is:
                            wig_file_path:TEX+/-(tex or notex):condition_id(intege
                            r):replicate_id(alphabet):strand(+ or -). If multiple
                            wig files need to be assigned, please use spaces to
                            separate the wig files. For example,
                            my_lib_tex_forward.wig:tex:1:a:+
                            my_lib_tex_reverse.wig:tex:1:a:-.
      --replicate_tex REPLICATE_TEX [REPLICATE_TEX ...], -rt REPLICATE_TEX [REPLICATE_TEX ...]
                            This value is the minimal number of replicates that a
                            TSS has to be detected in. The format is
                            $NUMBERofCONDITION_$NUMBERofREPLICATE. If different
                            --replicate_tex values need to be assigned to
                            different conditions, please use spaces to separate
                            them. For example, 1_2 2_2 3_3 means that
                            --replicate_tex is 2 in number 1 and number 2
                            conditions. In number 3 condition, --replicate_tex is
                            3. For assigning the same --replicate_tex to all
                            conditions, just use like all_1 (--replicate_tex is 1
                            in all conditions). Default is all_1.
      --condition_names CONDITION_NAMES [CONDITION_NAMES ...], -cn CONDITION_NAMES [CONDITION_NAMES ...]
                            The prefixs of output filename.If multiple conditions
                            need to be assigned, please use spaces to separate
                            them. For an example, prefix_condition1
                            prefix_condition2.
    
    additional arguments:
      --tsspredator_path TSSPREDATOR_PATH
                            Path of TSSpredator. Default is
                            /usr/local/bin/TSSpredator.jar
      --max_height MAX_HEIGHT, -he MAX_HEIGHT
                            This value relates to the minimum number of read
                            starts at a certain genomic position to be considered
                            as a TSS candidate. During optimization, --max_height
                            will be never larger than this value. Default is 2.5.
      --max_height_reduction MAX_HEIGHT_REDUCTION, -rh MAX_HEIGHT_REDUCTION
                            When comparing different genomes/conditions and the
                            step height threshold is reached in at least one
                            genome/condition, the threshold is reduced for the
                            other genomes/conditions by the value set here. This
                            value must be smaller than the step height threshold.
                            During optimization, --max_height_reduction will be
                            never larger than this value. Default is 2.4.
      --max_factor MAX_FACTOR, -fa MAX_FACTOR
                            The minimum factor by which the TSS height has to
                            exceed the local expression background. During
                            optimization, --max_factor will be never larger than
                            this value. Default is 10.
      --max_factor_reduction MAX_FACTOR_REDUCTION, -rf MAX_FACTOR_REDUCTION
                            When comparing different genomes/conditions and the
                            step factor threshold is reached in at least one
                            genome/condition, the threshold is reduced for the
                            other genomes/conditions by the value set here. This
                            value must be smaller than the step factor threshold.
                            During optimization, --max_factor_reduction will be
                            never larger than this value. Default is 9.9.
      --max_base_height MAX_BASE_HEIGHT, -bh MAX_BASE_HEIGHT
                            The minimum number of reads should be mapped on TSS.
                            During optimization, --max_base_height will be never
                            larger than this value. Default is 0.06.
      --max_enrichment_factor MAX_ENRICHMENT_FACTOR, -ef MAX_ENRICHMENT_FACTOR
                            The minimum enrichment factor. During optimization,
                            --max_enrichment_factor will be never larger than this
                            value. Default is 6.0.
      --max_processing_factor MAX_PROCESSING_FACTOR, -pf MAX_PROCESSING_FACTOR
                            The minimum processing factor. If the value for the
                            untreated library is higher than the treated library
                            the positionsis considered as a processing site and
                            not annotated as detected. During optimization,
                            --max_processing_factor will be never larger than this
                            value. Default is 6.0
      --utr_length UTR_LENGTH, -u UTR_LENGTH
                            The length of UTR. Default is 300.
      --cluster CLUSTER, -cu CLUSTER
                            This value defines the maximal distance (nucleotides)
                            between TSS candidates have to be clustered together.
                            Default is 2.
      --parallels PARALLELS, -c PARALLELS
                            Parallel runs for optimization. Default is 4.
      --steps STEPS, -s STEPS
                            The total runs for optimization. Default is 4000 runs.

- **Output files**

Based on the programs (TSS/processing site), Output files are stored in 
``$ANNOgesic/output/TSSs/optimized_TSSpredator`` or ``$ANNOgesic/output/processing_sites/optimized_TSSpredator``. 
Two output files are as following:

**stat_$GENOME.csv:** Stores the information of every run. The first column is the number of run.
The second column is the parameter set. ``he`` represents height; ``rh`` represents 
height reduction; ``fa`` means factor; ``rf`` means factor reduction; ``bh`` indicates 
base height; ``ef`` indicates enrichment factor; ``pf`` means processing factor. About the details 
of parameters, please refer to `TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_.
For an example, ``he_2.0_rh_1.8_fa_4.4_rf_2.8_bh_0.08_ef_3.0_pf_2.6`` means that height is 2.0, 
height reduction is 1.8, factor is 4.4, factor reduction is 2.8, base height is 0.08, 
enrichment factor is 3.0 and processing factor is 2.6. The third column is 
the number of true positive. The fourth column is true positive rate. The fifth 
columns is the number of false positive. The sixth is false positive rate. 
The seventh column is the number of false negatives. The eighth 
column is missing rate.

**best_$GENOME.csv:** Stores the best parameter set. The meanings of all columns are the same as ``stat_$GENOME.csv``.

.. _screenshot:

screenshot (screenshot generation)
----------------------------------

``screenshot`` can generate batch files for producing screenshot of `IGV <https://www.broadinstitute.org/igv>`_. 
Generating screenshots can reduce the time for checking the results in genome browser.
When the batch files is produced, the user just needs to open `IGV <https://www.broadinstitute.org/igv>`_, then presses ``tools`` 
on the top tags and choose ``run batch script``. The program will automatically produce screenshots. 

- **Required tools**

`IGV <https://www.broadinstitute.org/igv>`_.

- **Required files**

**Gff files that the user wants to produce screenshots:** All screenshots will be produced based on the positions of ``--main_gff``. 
If comparing ``--main_gff`` with other features is required, please assign gff files of other features to ``--side_gffs``. 

**Fasta files of the genomes**

**Wiggle files of TEX+/- or fragmented/conventional libraries:** Please check the section ``The format of libraries for import to ANNOgesic``.

- **Arguments**

::

    usage: annogesic screenshot [-h] --project_path PROJECT_PATH --fasta_file
                                FASTA_FILE --main_gff MAIN_GFF
                                [--side_gffs SIDE_GFFS [SIDE_GFFS ...]]
                                [--tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]]
                                [--frag_libs FRAG_LIBS [FRAG_LIBS ...]]
                                --output_folder OUTPUT_FOLDER [--height HEIGHT]
                                [--present {expand,collapse,squish}]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --fasta_file FASTA_FILE, -f FASTA_FILE
                            Path of the genome fasta file.
      --main_gff MAIN_GFF, -mg MAIN_GFF
                            Screenshot will be generated based on the positions of
                            genomic features in this gff file.
      --side_gffs SIDE_GFFS [SIDE_GFFS ...], -sg SIDE_GFFS [SIDE_GFFS ...]
                            The gff files of other genomic feature (besides
                            --main_gff). For assigning multiple files, please use
                            space to separated them.
      --tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...], -tl TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]
                            TEX+/- wig files. The format is:
                            wig_file_path:TEX+/-(tex or notex):condition_id(intege
                            r):replicate_id(alphabet):strand(+ or -). If multiple
                            wig files need to be assigned, please use spaces to
                            separate the wig files. For example,
                            my_lib_tex_forward.wig:tex:1:a:+
                            my_lib_tex_reverse.wig:tex:1:a:-.
      --frag_libs FRAG_LIBS [FRAG_LIBS ...], -fl FRAG_LIBS [FRAG_LIBS ...]
                            Wig files of RNA-Seq with transcript fragmented. The
                            format is: wig_file_path:frag:condition_id(integer):re
                            plicate_id(alphabet):strand(+ or -). If multiple wig
                            files need to be assigned, please use spaces to
                            separate the wig files. For example,
                            my_lib_frag_forward.wig:frag:1:a:+
                            my_lib_frag_reverse.wig:frag:1:a:-.
      --output_folder OUTPUT_FOLDER, -o OUTPUT_FOLDER
                            Path of the output folder. It will create a sub-folder
                            "screenshots" under --output_folder to store the
                            results.
    
    additional arguments:
      --height HEIGHT, -he HEIGHT
                            Height of the screenshot. Default is 1500.
      --present {expand,collapse,squish}, -p {expand,collapse,squish}
                            The presentation types (expand, collapse, or squish)
                            of the features in the screenshot. Default is expand.

- **Output files**

Based on the paths of ``--main_gff``, ``screenshot`` will generate a folder - ``screenshots`` under the 
folder of ``--main_gff``. Output files will be stored in this folder. the output files and folders are 
as following:

**forward.txt:** Batch file of the forward strand for running on IGV.

**reverse.txt:** Batch file of reverse strand for running on IGV.

**forward:** Folder for storing screenshots of the forward strand.

**reverse:** Folder for storing screenshots of the reverse strand.

When batch files are executed on IGV, the screenshots will be automatically stored in the folder called ``forward`` and ``reverse``. 
Format of the filenames will be ``$GENOME:$START-$END.png``. For an example, ``NC_007795:1051529-1051696.png`` 
means the genome is NC_007795, the feature's start point is 1051529 and end point is 
1051696.

.. _color_png:

colorize_screenshot_tracks (colorize the tracks of screenshots)
---------------------------------------------------------------

``colorize_screenshot_tracks`` is a following procedure of ``screenshot``. If numerous samples are included in one figure, 
Tracks will be difficult to check. ``colorize_screenshot_tracks`` can color the tracks based on TEX +/- libraries 
for improving the checking process.

- **Required tools**

`ImageMagick <http://www.imagemagick.org/script/index.php>`_.

- **Required files**

**The screenshots folder:** Please make sure the folders of ``forward`` and ``reverse`` 
exist in the folder of ``screenshots``.

- **Arguments**

::

    usage: annogesic colorize_screenshot_tracks [-h] --project_path PROJECT_PATH
                                                --screenshot_folder
                                                SCREENSHOT_FOLDER --track_number
                                                TRACK_NUMBER
                                                [--imagemagick_covert_path IMAGEMAGICK_COVERT_PATH]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --screenshot_folder SCREENSHOT_FOLDER, -f SCREENSHOT_FOLDER
                            The folder containing "screenshots" which generated
                            from subcommand "screenshot" (storage of screenshots).
      --track_number TRACK_NUMBER, -t TRACK_NUMBER
                            The number of tracks.
    
    additional arguments:
      --imagemagick_covert_path IMAGEMAGICK_COVERT_PATH, -m IMAGEMAGICK_COVERT_PATH
                            Path of "convert" in ImageMagick package.

- **Output files**

The new screenshots will replace the previous ones automatically.

.. _merge_feature:

merge_features (merge all annotation features)
----------------------------------------------
If merging multiple features of the annotation to one gff file is needed, ``merge_features`` can achieve this purpose. 
``merge_features`` can merge all features (the user assigned) to one gff file, and search the parent transcript to each feature.

- **Required files**

**Gff files of features that the user wants to merge:** 
If transcript gff files can be provided, this module will search the parent transcripts to other input features.

- **Arguments**

::

    usage: annogesic merge_features [-h] --project_path PROJECT_PATH
                                    --output_prefix OUTPUT_PREFIX
                                    [--transcript_file TRANSCRIPT_FILE]
                                    [--other_features_files OTHER_FEATURES_FILES [OTHER_FEATURES_FILES ...]]
                                    [--terminator_fuzzy TERMINATOR_FUZZY]
                                    [--tss_fuzzy TSS_FUZZY]
    
    optional arguments:
      -h, --help            show this help message and exit
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --output_prefix OUTPUT_PREFIX, -op OUTPUT_PREFIX
                            The prefix name of output gff file. The filename will
                            be $OUTPUT_PREFIX_merge_features.gff.
      --transcript_file TRANSCRIPT_FILE, -a TRANSCRIPT_FILE
                            Path of transcript gff file. The parent transcripts
                            ("Parent" in gff attributes) of all features will be
                            generated.
      --other_features_files OTHER_FEATURES_FILES [OTHER_FEATURES_FILES ...], -of OTHER_FEATURES_FILES [OTHER_FEATURES_FILES ...]
                            Paths of the gff files (besides transcript gff file).
                            For assigning multiple gff files, please use spaces to
                            separate them.
    
    additional arguments:
      --terminator_fuzzy TERMINATOR_FUZZY, -ef TERMINATOR_FUZZY
                            The fuzzy nucleotides for comparing transcripts and
                            terminators. Default is 30.
      --tss_fuzzy TSS_FUZZY, -tf TSS_FUZZY
                            The fuzzy value for comparing TSSs and transcripts.
                            Default is 5.

- **Output files**

Output gff files are stored in ``$ANNOgesic/output/merge_all_features``. The tag - ``Parent`` in the attributes of 
gff file shows the parent transcript.
