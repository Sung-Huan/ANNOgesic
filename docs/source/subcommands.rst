.. _ANNOgesic's subcommands:

ANNOgesic's subcommands
=======================

In general, the subcommands need at least one argument - the analysis
folder. If it is not given, ANNOgesic assumes the current
folder as the analysis folder. All the default settings were tested 
on several different organisms and produced in general good results. 
A user can perform an analysis rather easily which those default values 
and step-wise adapt the parameters to optimize the results if needed.

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

log file for storing the details of each module
===============================================

If the user needs to check the details of the running of each module, ``log.txt`` can 
be found in the output folder of each module.

.. _create:

create (create analysis folder)
-------------------------------

``create`` generates the folders for analysis. Once created, please move the required files 
into the corresponding folders.

The folders are following:

**input:** Stores all input files.

	**BAMs:** For ``.bam`` files. ``BAMs_map_related_genomes`` 
	is for the ``.bam`` files which are mapped on closely related genomes of the reference genomes.
	``BAMs_map_reference_genomes`` is for the ``.bam`` files which are mapped on the reference genomes.

	**databases:** For all databases.

	**manual_TSSs:** If the manual detected transcription starting sites (TSSs) can be provided,
	it can be stored here for running ``TSS_optimization`` or merging 
	the automatic predicted ones and manual detected ones. Please use gff3 format.

	**manual_processing_sites:** It is similar to ``manual_TSS``, but for 
	processing sites.

	**mutation_table:** If the mutation table between the closely ralated genomes and 
	reference genomes is provided, please put the file here. Please check 
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

    usage: annogesic create --project_path PROJECT_PATH
    
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

- **Basic arguments**


::

    usage: annogesic get_input_files --project_path PROJECT_PATH
                                     --ftp_path FTP_PATH [--ref_fasta]
                                     [--ref_gff] [--ref_gbk]
    
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

- **Additional arguments**

::

    additional arguments:
      --ref_ptt, -p         Download ptt files of the reference. Default is False.
      --ref_rnt, -r         Download rnt files of the reference. Default is False.
      --convert_embl, -e    Convert gbk to embl files of the reference. Default is
                            False.

    optional arguments:
      -h, --help            show this help message and exit

- **Output files**

Output files will be stored in ``$ANNOgesic_folder/input/reference``

Output folder names are following:

**fasta:** Fasta files.

**annotation:** Annotation files.

.. _update_genome_fasta:

update_genome_fasta (update reference genome fasta file)
--------------------------------------------------------

If fasta files of the reference genomes do not exist, ``update_genome_fasta`` can 
update fasta files from the closely related genomes to our reference genomes 
via searching the mutations. 
Therefore, a table of mutation information is required. For the format of the table, please check 
`mutation table <https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/tutorial_data/mutation.csv>`_.
Titles of the columns are presented on the top and they need to start with ``#``. 
The format is basically the same as VCF format.

If no mutation information is provided, ``snp`` can be used for detecting mutations. 
(one module of ``ANNOgesic``). Please check the section of :ref:`snp`.

- **Required files**

**Fasta files of reference genome**

**Mutation table:** Contains the information of mutations between related and reference genomes.
For an example, please check `mutation table <https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/tutorial_data/mutation.csv>`_.

- **Arguments**

::
    
    usage: annogesic update_genome_fasta [-h] --project_path PROJECT_PATH
                                         --related_fasta_files RELATED_FASTA_FILES
                                         [RELATED_FASTA_FILES ...]
                                         --mutation_table MUTATION_TABLE
                                         --updated_seq_name UPDATED_SEQ_NAME
    
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
                            information between the reference genome and genome of
                            the closely related species. For an example check
                            https://github.com/Sung-
                            Huan/ANNOgesic/blob/master/tutorial_data/mutation.csv
      --updated_seq_name UPDATED_SEQ_NAME, -u UPDATED_SEQ_NAME
                            The file name of the updated sequence. The output
                            fasta file name will be --updated_seq_name.fa.

- **Output files**

**Fasta files of updated genome**: The updated fasta files are stored in ``$ANNOgesic_folder/output/updated_references/fasta_files``.

.. _annotation_transfer:

annotation_transfer (annotation transfer)
-----------------------------------------

``annotation transfer`` is the subcommand for transferring the annotation from the closely related genomes 
to the reference genomes. To achieve this, `RATT <http://www.sanger.ac.uk/resources/software/pagit/>`_ 
is integrated in ANNOgesic. The higher similarity between closely related genomes and reference genomes are, 
the more precise the performance is. Before running ``annotation transfer``, 
please run ``source $PAGIT_HOME/sourceme.pagit`` first. it will modify the path for executing RATT. 
If you use Dockerfile to execute ANNOgesic, the path modification can be skipped.
If the error message related to 'defined(@array)' shows, please check :ref:`FAQ`.

- **Required tools**

`RATT <http://www.sanger.ac.uk/resources/software/pagit/>`_.

- **Required files**

**Annotation files of the closely related genomes**: Genbank/embl files of the closely related genomes.

**Fasta files of the closely related genomes**

**Fasta files of the updated genomes**

- **Basic arguments**

::

    usage: annogesic annotation_transfer --project_path PROJECT_PATH
                                     --compare_pair COMPARE_PAIR
                                     [COMPARE_PAIR ...]
                                     [--related_embl_files RELATED_EMBL_FILES [RELATED_EMBL_FILES ...]]
                                     [--related_gbk_files RELATED_GBK_FILES [RELATED_GBK_FILES ...]]
                                     --related_fasta_files RELATED_FASTA_FILES
                                     [RELATED_FASTA_FILES ...]
                                     --target_fasta_files TARGET_FASTA_FILES
                                     [TARGET_FASTA_FILES ...]
                                     [--additional arguments]
    
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
                            the headers contain space or '|', only the string from
                            '>' to the first space or '|' will be considered as
                            the name for --compare_pair (normally this part is the
                            accession number). If multiple sequences need to be
                            assigned, please use spaces to separate them.
      --related_embl_files RELATED_EMBL_FILES [RELATED_EMBL_FILES ...], -ce RELATED_EMBL_FILES [RELATED_EMBL_FILES ...]
                            The paths of the embl files of the related species. If
                            --related_embl_files is assigned, --related_gbk_files
                            is not needed.
      --related_gbk_files RELATED_GBK_FILES [RELATED_GBK_FILES ...], -cg RELATED_GBK_FILES [RELATED_GBK_FILES ...]
                            The paths of the genbank files of the related species.
                            The genbank can be ended by .gbk, .gbff or .gb. If
                            --related_gbk_files is assigned, --related_embl_files
                            is not needed.
      --related_fasta_files RELATED_FASTA_FILES [RELATED_FASTA_FILES ...], -cf RELATED_FASTA_FILES [RELATED_FASTA_FILES ...]
                            The paths of the fasta files of the related species.
      --target_fasta_files TARGET_FASTA_FILES [TARGET_FASTA_FILES ...], -tf TARGET_FASTA_FILES [TARGET_FASTA_FILES ...]
                            The paths of the target fasta files.

- **Additional arguments**

::

    additional arguments:
      --ratt_path RATT_PATH
                            Path of the start.ratt.sh file of RATT folder. Default
                            is start.ratt.sh.
      --element ELEMENT, -e ELEMENT
                            --element will become the prefix of all output
                            file.Default is "chromosome".
      --transfer_type TRANSFER_TYPE, -t TRANSFER_TYPE
                            The transfer type for running RATT. (For the details,
                            please refer to the manual of RATT.) Default is
                            Strain.

    optional arguments:
      -h, --help            show this help message and exit


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
reference genomes if it is necessary.

- **Required tools**

`Samtools <https://github.com/samtools>`_.

`Bcftools <https://github.com/samtools>`_.

- **Required files**

**BAM files:** BAM files from fragmented/conventional libraries or TEX +/- treated libraries both can be accepted.
For assigning the files, please follow the format -- ``$SET_NAME:$BAMFILE1,$BAMFILE2,...``. 
For an example, the user has four bam files of one genome. Then the input will be 
``set1:sample1.bam,sample2.bam,sample3.bam,sample4.bam``.

**Fasta files of the closely related genomes** or **Fasta files of the reference genomes**

- **Basic arguments**

::

    usage: annogesic snp [-h] --project_path PROJECT_PATH --bam_type
                         {related_genome,reference_genome} --program
                         {with_BAQ,without_BAQ,extend_BAQ}
                         [{with_BAQ,without_BAQ,extend_BAQ} ...] --fasta_files
                         FASTA_FILES [FASTA_FILES ...] --bam_files BAM_FILES
                         [BAM_FILES ...] [--additional arguments]
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --bam_type {related_genome,reference_genome}, -t {related_genome,reference_genome}
                            If the BAM files are produced by mapping to a related
                            genome, please assign "related_genome". the mutations
                            between the related genome and the refernce genome can
                            be detected for generating sequence of the query
                            genome. If the BAM files are produced by mapping to
                            the reference genome, please assign
                            "reference_genome". The mutations of reference genome
                            can be detected.
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

- **Additional arguments**

::

    additional arguments:
      --samtools_path SAMTOOLS_PATH
                            Path of samtools.
      --bcftools_path BCFTOOLS_PATH
                            Path of bcftools.
      --quality QUALITY, -q QUALITY
                            The minimum quality score of a mutation. Default is
                            40.
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
                            The query genome is haploid or diploid. Default is
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
                            1. fix number ("r"), 2. times of the number of bam
                            files (counted from --bam_files) ("n") or 3. times of
                            average read depth ("a"). For example, n_10,0.8. If
                            the BAM files are 4, it means the minimum mapped reads
                            of a SNP is 40 (4 * 10), and the minimum ratio of
                            mapped read of a SNP (mapped reads of a SNP / total
                            reads) is 0.8. Default is n_10,0.8.
      --indel_fraction INDEL_FRACTION, -if INDEL_FRACTION
                            The minimum IDV and IMF which supports for insertion
                            of deletion. The minimum IDV can be assigned by
                            different types: 1. fix number ("r"), 2. times of the
                            number of bam files (assigned by --bam_files) ("n") or
                            3. times of the average read depth ("a"). The input
                            format is $MIN_IDF:$MIN_IMF. For example, The value is
                            n_10,0.8 and 4 BAM files are assigned. The minimum IDV
                            is 40, and the minimum IMF is 0.8. Default is
                            n_10,0.8.
      --filter_tag_info FILTER_TAG_INFO [FILTER_TAG_INFO ...], -ft FILTER_TAG_INFO [FILTER_TAG_INFO ...]
                            For using more filters to improve the detection.
                            Please assign 1. the name of a tag, 2. bigger ("b") or
                            smaller ("s") and 3. the value of the filter. For
                            example, "RPB_b0.1,MQ0F_s0" means that RPB should be
                            bigger than 0.1 and MQ0F should be smaller than 0.
                            Default is RPB_b0.1,MQSB_b0.1,MQB_b0.1,BQB_b0.1.

    optional arguments:
      -h, --help            show this help message and exit

- **Output files**

If ``bam_type`` is ``related_genome``, 
the results will be stored in ``$ANNOgesic/output/SNP_calling/compare_related_and_reference_genomes``. 
If ``bam_type`` is ``reference_genome``, the results will be stored in ``$ANNOgesic/output/SNP_calling/mutations_of_reference_genomes``.

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

The example contains ``position conflict`` and ``mutation conflict``.
As a result, the conflicts will affect the other mutation's positions.
Therefore, it will generate four different fasta files. The first two lines are ``position conflict``, and 
the last line is ``mutation conflict``.
``$GENOME_$PROGRAM_$SET_seq_reference.csv`` is the index for these four fasta files.

::

   1       1632629 1       1499572:TT      Staphylococcus_aureus_HG003
   1       1632629 2       1499572:TTTTT   Staphylococcus_aureus_HG003
   2       1632630 1       1499572:TT      Staphylococcus_aureus_HG003
   2       1632630 2       1499572:TTTTT   Staphylococcus_aureus_HG003

The first column is the index of the ``position conflict``. 
The second column is the selected position.
The third one is the index of the ``mutations conflict``. 
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

- **Basic arguments**

::

    usage: annogesic tss_ps --project_path PROJECT_PATH --program {TSS,PS}
                            --fasta_files FASTA_FILES [FASTA_FILES ...]
                            --annotation_files ANNOTATION_FILES
                            [ANNOTATION_FILES ...] --tex_notex_libs TEX_NOTEX_LIBS
                            [TEX_NOTEX_LIBS ...]
                            --condition_names CONDITION_NAMES
                            [CONDITION_NAMES ...]
                            [--additional arguments]
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --program {TSS,PS}, -p {TSS,PS}
                            The feature to predict. Please assign "TSS" or "PS".
                            Default is "TSS".
      --fasta_files FASTA_FILES [FASTA_FILES ...], -f FASTA_FILES [FASTA_FILES ...]
                            Paths of the genome fasta files.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the genome annotation gff files containing
                            CDSs, tRNAs, rRNAs, etc.
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

- **Additional arguments**

::

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
      --tolerance TOLERANCE, -to TOLERANCE
                            The 5'ends of transcripts will be extended or withdrew
                            by this value (nucleotides) for searching the
                            associated TSSs (--compare_transcript_files is
                            provided). Default is 5.
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

    optional arguments:
      -h, --help            show this help message and exit

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
	``stat_gene_vali_$GENOME.csv`` is for comparing TSS with genome annotations like CDSs.

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

**Gff files of genome anntations containing CDSs, tRNAs, rRNAs, etc:** If the user wants to compare transcripts with genome annotations or modify transcript based on genome annotations 
like CDSs, tRNAs, rRNAs, genome annotation gff files are required. There are four options for modification of transcripts:

	**merge_overlap:** If multiple transcripts overlap the same gene, they will be merged as one complete transcript.

	**extend_3end:** If the transcript starts at the upstream of the gene and ends within the gene, 
	the end point of the transcript will be extended to the end point of gene.

	**extend_5end:** If the transcript starts within the gene and ends at the downstream of gene, 
	the starting point of the transcript will be extended to the starting point of the gene.

        **within_extend_ends:** If the transcript is within the gene, the two ends of the transcript will be 
	extended to the two ends of gene.

	**none:** Transcripts will not be modified by the genome annotations

- **Basic arguments**

::

    usage: annogesic transcript --project_path PROJECT_PATH
                            [--annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...]]
                            [--modify_transcript {merge_overlap,extend_3end,extend_5end,within_extend_ends,none} [{merge_overlap,extend_3end,extend_5end,within_extend_ends,none} ...]]
                            [--tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]]
                            [--frag_libs FRAG_LIBS [FRAG_LIBS ...]]
                            [--replicate_tex REPLICATE_TEX [REPLICATE_TEX ...]]
                            [--replicate_frag REPLICATE_FRAG [REPLICATE_FRAG ...]]
                            [--tex_notex {1,2}] [--additional arguments]
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the genome annotation gff files containing
                            CDSs, tRNAs, rRNAs, etc. TSS gff files and terminator
                            gff files need to be separately assigned to
                            --tss_files and --terminator_files, respectively.
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
                            Similar to --replicates_tex. This value is for
                            fragmented (or conventional) libraries.
      --tex_notex {1,2}, -te {1,2}
                            The value is for TEX+/- libraries to decide the
                            transcript should be detected in both (TEX+ and TEX-)
                            or can be detected in only one library (TEX+ or TEX-).
                            Please assign 1 or 2. Default is 1.

- **Additional arguments**

::

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
      --tss_tolerance TSS_TOLERANCE, -tt TSS_TOLERANCE
                            The 5'ends of transcripts will be extended or withdrew
                            by this value (nucleotides) for searching the
                            associated TSSs (--tss_files is provided). Default is
                            5.
      --terminator_files TERMINATOR_FILES [TERMINATOR_FILES ...], -e TERMINATOR_FILES [TERMINATOR_FILES ...]
                            Paths of terminator gff files for comparing
                            transcripts and terminators. Default is None.
      --terminator_tolerance TERMINATOR_TOLERANCE, -et TERMINATOR_TOLERANCE
                            The 3'ends of transcripts will be extended or withdrew
                            by this value (nucleotides) for searching the
                            associated terminators. Default is 30.
      --max_length_distribution MAX_LENGTH_DISTRIBUTION, -mb MAX_LENGTH_DISTRIBUTION
                            For generating the figure of distribution of
                            transcript length, please assign the maximum length.
                            Default is 2000.

    optional arguments:
      -h, --help            show this help message and exit

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

	**Avg_coverage:$LIB_NAME:** Stores the average coverage information of the libraries about this transcript.

**statistics:** Stores statistic files.

	**Comparing transcript with other features:** ``stat_compare_transcript_genome_$GENOMENAME.csv`` is 
	for comparing transcript with genome annotation like CDSs, ``stat_compare_transcript_TSS_$GENOMENAME.csv`` is for comparing 
	transcript with TSS, and ``stat_compare_transcript_terminator_$GENOMENAME.csv`` is for comparing
        transcript with terminator.

	**Figure of the distribution of transcript length:** ``$GENOME_length_all.png`` is for analyzing of all transcript length. 
	``$GENOME_length_less_$LENGTH.png`` is for the analyzing of the assigned length.

**gffs:** Stores gff files of transcripts. Some useful information can be found in the tags of the attributes within the transcript gff file.
Based on this information, we can know the details of the specific transcript. The tags are as following:

	**compare_$FEATURE:** State of overlap between transcripts and features
	(If ``--compare_feature_genome`` and ``--annotation_files`` are assigned). The value may be ``cover``, ``right_shift``, ``left_shift``, ``within`` or ``no_related``.

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

**Gff files of the genome annotations containing CDSs, tRNAs, rRNAs, etc**

**Fasta files of the genome sequences**

**Wiggle files of TEX +/- treated libraries or fragmented/conventional libraries**

**Gff files of the transcripts**

- **Basic arguments**

::

    usage: annogesic terminator --project_path PROJECT_PATH --fasta_files
                                FASTA_FILES [FASTA_FILES ...] --annotation_files
                                ANNOTATION_FILES [ANNOTATION_FILES ...]
                                --transcript_files TRANSCRIPT_FILES
                                [TRANSCRIPT_FILES ...]
                                [--tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]]
                                [--frag_libs FRAG_LIBS [FRAG_LIBS ...]]
                                [--tex_notex {1,2}]
                                [--replicate_tex REPLICATE_TEX [REPLICATE_TEX ...]]
                                [--replicate_frag REPLICATE_FRAG [REPLICATE_FRAG ...]]
                                [--additional arguments]
    
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --fasta_files FASTA_FILES [FASTA_FILES ...], -f FASTA_FILES [FASTA_FILES ...]
                            Paths of the genome fasta files.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the genome annotation gff files containing
                            CDSs, tRNAs, rRNAs, etc. Transcript gff files and sRNA
                            gff files need to be separately assigned to
                            --transcript_files and --srna_files, respectively.
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
                            Similar to --replicates_tex. This value is for
                            fragmented (or conventional) libraries.

- **Additional arguments**

::

    additional arguments:
      --transterm_path TRANSTERM_PATH
                            Path of "transterm" in TransTermHP.
      --expterm_path EXPTERM_PATH
                            Path of expterm.dat for TransTermHP. Default is
                            /usr/local/bin/expterm.dat
      --rnafold_path RNAFOLD_PATH
                            Path of RNAfold of Vienna package.
      --srna_files SRNA_FILES [SRNA_FILES ...], -sr SRNA_FILES [SRNA_FILES ...]
                            Paths of sRNA gff files if sRNA information need to be
                            considered as well.
      --decrease DECREASE, -d DECREASE
                            The maximum ratio -- (lowest coverage / highest
                            coverage) within (or nearby) the terminator. If the
                            ratio is smaller than --decrease, the candidate will
                            be considered as highly-confidence terminator. Default
                            is 0.5.
      --tolerance_detect_coverage TOLERANCE_DETECT_COVERAGE, -tc TOLERANCE_DETECT_COVERAGE
                            The extended region (nucleotides) of the terminators
                            for detecting coverage significant drop. For example,
                            the location of terminator is 300-400, and
                            --tolerance_detect_coverage is 30. If the coverage
                            decrease is detected within 270-430, this candidate is
                            still considered as the terminator which have coverage
                            dramatic decrease. Default is 30.
      --tolerance_within_transcript TOLERANCE_WITHIN_TRANSCRIPT, -tut TOLERANCE_WITHIN_TRANSCRIPT
                            If the candidates are within transcript and the
                            distance (nucleotides) between the end of transcript
                            and terminator is within this value, the candidate
                            will be considered as a terminator. Otherwise, it will
                            be removed. Default is 30.
      --tolerance_downstream_transcript TOLERANCE_DOWNSTREAM_TRANSCRIPT, -tdt TOLERANCE_DOWNSTREAM_TRANSCRIPT
                            The meaning is similar to
                            --tolerance_within_transcript. This value is for the
                            candidates which are at the downstream of transcript.
                            Default is 30.
      --tolerance_within_gene TOLERANCE_WITHIN_GENE, -twg TOLERANCE_WITHIN_GENE
                            The meaning is similar to
                            --tolerance_within_transcript. This value is for gene
                            in stead of transcript. Default is 10.
      --tolerance_downstream_gene TOLERANCE_DOWNSTREAM_GENE, -tdg TOLERANCE_DOWNSTREAM_GENE
                            The meaning is similar to
                            --tolerance_downstream_transcript. This value is for
                            gene in stead of transcript. Default is 310.
      --highest_coverage HIGHEST_COVERAGE, -hc HIGHEST_COVERAGE
                            The minimum value of the highest coverage of
                            terminator. The low expressed terminator are not
                            included in "best_candidates", but are still in
                            "all_candidates". Default is 10.
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
      --min_u_tail MIN_U_TAIL, -mu MIN_U_TAIL
                            The minimum number of U in poly U-tail of terminator.
                            Default is 5.
      --mutation_u_tail MUTATION_U_TAIL, -uu MUTATION_U_TAIL
                            The number of nts which are not U can be tolerated.
                            Default is 2.
      --keep_multi_term, -kp
                            Sometimes, one gene is associated with multiple
                            terminators In default, it will only keep the highly-
                            confidence one. This flag can keep all terminators
                            which are associated with the same gene. Default is
                            False.
    optional arguments:
      -h, --help            show this help message and exit

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

	**Coverage_$LIB_NAME:** Shows the coverage information of the libraries about this terminator. ``high`` means the highest coverage of the libraries, 
	``low`` means the lowest coverage of the libraries, and ``diff`` represents the difference between ``high`` and ``low``. If ``No_coverage_decreasing`` is showed, 
	it means this terminator reveal gene expression but no coverage decrease. If ``NA`` is showed, it means that this terminator has no gene expression.

.. _utr:

utr (UTR detection)
-------------------

``utr`` can compare TSSs, CDSs/tRNAs/sRNAs, transcripts and terminators
to generate 5'UTR and 3'UTR. 5'UTRs are based on detecting the regions between TSSs and CDSs/tRNAs/sRNAs. 
3'UTRs are based on detecting the 
regions between the end of the transcripts and CDSs/tRNAs/sRNAs. If the input gff files of TSSs are not computed by 
ANNOgesic, please use ``--tss_source`` to classify TSSs for the analysis.

- **Required files**

**Gff files of the genome annotations containing CDSs, tRNAs, rRNAs, etc**

**Gff files of the TSSs**

**Gff files of the transcripts**

- **Optional input files**

**Gff files of the terminators:** If the information of terminators is needed, the gff files of terminators are required.

- **Basic Arguments**

::

    usage: annogesic utr --project_path PROJECT_PATH --annotation_files
                         ANNOTATION_FILES [ANNOTATION_FILES ...] --tss_files
                         TSS_FILES [TSS_FILES ...] --transcript_files
                         TRANSCRIPT_FILES [TRANSCRIPT_FILES ...]
                         [--terminator_files TERMINATOR_FILES [TERMINATOR_FILES ...]]
                         [--additional arguments]
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the genome annotation gff files containing
                            CDSs, tRNAs, rRNAs, etc. Gff files of TSSs,
                            terminators, and transcripts need to be separately
                            assigned to --tss_files, terminator_files, and
                            transcript_files, respectively.
      --tss_files TSS_FILES [TSS_FILES ...], -t TSS_FILES [TSS_FILES ...]
                            Paths of the TSS files.
      --transcript_files TRANSCRIPT_FILES [TRANSCRIPT_FILES ...], -a TRANSCRIPT_FILES [TRANSCRIPT_FILES ...]
                            Paths of the transcript gff files.
      --terminator_files TERMINATOR_FILES [TERMINATOR_FILES ...], -e TERMINATOR_FILES [TERMINATOR_FILES ...]
                            If the paths of terminator files are assigned,
                            terminator will be used to detect 3'UTR.

- **Additional arguments**

::

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
      --terminator_tolerance TERMINATOR_TOLERANCE, -et TERMINATOR_TOLERANCE
                            The 3'ends of transcripts will be extended or withdrew
                            by this value (nucleotides) for searching the
                            associated terminators. Default is 30.
      --tolerance_3utr TOLERANCE_3UTR, -t3 TOLERANCE_3UTR
                            The length of 3'UTR can be extended or withdrew by
                            this value (nucleotides). It only works when
                            transcript information is provided. Default is 10.
      --tolerance_5utr TOLERANCE_5UTR, -t5 TOLERANCE_5UTR
                            The length of 5'UTR can be extended or withdrew by
                            this value (nucleotides). It only works when
                            transcript information is provided. Default is 5.

    optional arguments:
      -h, --help            show this help message and exit

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

**Gff files of the genome annotations containing CDSs, tRNAs, rRNAs, etc**

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
this candidate will be included to the result without considering other filters. ``--blast_score_srna`` 
and ``--blast_e_srna`` can be used for adjustment of the prediction.

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
this candidates will be removed. ``--blast_score_nr`` and ``--blast_e_nr`` can be used for adjustment of the prediction.

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

- **Basic arguments**

::

    usage: annogesic srna --project_path PROJECT_PATH [--utr_derived_srna]
                          [--filter_info {tss,sec_str,blast_nr,blast_srna,sorf,term,promoter,none} [{tss,sec_str,blast_nr,blast_srna,sorf,term,promoter,none} ...]]
                          --transcript_files TRANSCRIPT_FILES
                          [TRANSCRIPT_FILES ...] --annotation_files
                          ANNOTATION_FILES [ANNOTATION_FILES ...]
                          [--tss_files TSS_FILES [TSS_FILES ...]]
                          [--fasta_files FASTA_FILES [FASTA_FILES ...]]
                          [--compute_sec_structures]
                          [--srna_database_path SRNA_DATABASE_PATH]
                          [--nr_database_path NR_DATABASE_PATH]
                          [--tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]]
                          [--frag_libs FRAG_LIBS [FRAG_LIBS ...]]
                          [--tex_notex {1,2}]
                          [--replicate_tex REPLICATE_TEX [REPLICATE_TEX ...]]
                          [--replicate_frag REPLICATE_FRAG [REPLICATE_FRAG ...]]
                          [--additional arguments]
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --utr_derived_srna, -u
                            Assign to detect UTR-derived sRNA. Default is False.
      --filter_info {tss,sec_str,blast_nr,blast_srna,sorf,term,promoter,none} [{tss,sec_str,blast_nr,blast_srna,sorf,term,promoter,none} ...], -d {tss,sec_str,blast_nr,blast_srna,sorf,term,promoter,none} [{tss,sec_str,blast_nr,blast_srna,sorf,term,promoter,none} ...]
                            The filters for improving the sRNA detection: 1. tss
                            (sRNA has to start with a TSS), 2. sec_str (free
                            energy change of secondary structure (normalized by
                            length) has to be smaller than --cutoff_energy), 3.
                            blast_nr (the number of the homologs in the non-
                            redundant database has to be below the --cutoff_nr_hit
                            ), 4. blast_srna (as long as the homologs can be found
                            in the sRNA database, the candidates will be included
                            to the best candidtes without considering other
                            filters), 5. sorf (sRNA must not overlap with sORFs),
                            6. term (sRNA has to be associated with a terminator),
                            7. promoter (sRNA has to be associated with a promoter
                            motif). For using multiple filters, please separated
                            them by spaces. If blast_srna was assigned, the
                            headers of sequences in sRNA database should be
                            $ID|$GENOME|$SRNANAME. "tss sec_str blast_nr
                            blast_srna" are recommended to be used. If "none" is
                            assigned, no filters are applied. Default is tss
                            sec_str blast_nr blast_srna.
      --transcript_files TRANSCRIPT_FILES [TRANSCRIPT_FILES ...], -a TRANSCRIPT_FILES [TRANSCRIPT_FILES ...]
                            Paths of the transcript files.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the genome annotation gff files containing
                            CDSs, tRNAs, rRNAs, etc. Gff files of transcripts,
                            TSSs, terminators, processing sites, and sORFs need to
                            be separately assigned to --transcript_files,
                            --tss_files, --terminator_files,
                            --processing_site_files, and --sorf_files,
                            respectively.
      --tss_files TSS_FILES [TSS_FILES ...], -t TSS_FILES [TSS_FILES ...]
                            Paths of TSS gff files. For detecting UTR-derived sRNA
                            or "tss" in --filter_info, TSS gff files MUST be
                            provided.
      --fasta_files FASTA_FILES [FASTA_FILES ...], -f FASTA_FILES [FASTA_FILES ...]
                            paths of fasta files of reference genome, If
                            "sec_str", "blast_nr" or "blast_srna" is assigned to
                            --filter_info or --search_poly_u is not 0, fasta files
                            are required.
      --compute_sec_structures, -cs
                            Computing secondary structures of sRNAs. Default is
                            False.
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
                            Wig files of RNA-Seq with fragmented transcripts. The
                            format is: wig_file_path:frag:condition_id(integer):re
                            plicate_id(alphabet):strand(+ or -). If multiple wig
                            files need to be assigned, please use spaces to
                            separate the wig files. For example,
                            my_lib_frag_forward.wig:frag:1:a:+
                            my_lib_frag_reverse.wig:frag:1:a:-.
      --tex_notex {1,2}, -te {1,2}
                            If TEX+/- libraries are assigned, a sRNA should be
                            detected in both (TEX+ and TEX-) or needs to be
                            detected in only one library (TEX+ or TEX-). Default
                            is 2.
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
                            conditions, use all_1 (--replicate_tex is 1 in all
                            conditions). Default is all_1.
      --replicate_frag REPLICATE_FRAG [REPLICATE_FRAG ...], -rf REPLICATE_FRAG [REPLICATE_FRAG ...]
                            Similar to --replicates_tex. This value is for
                            libraries with fragmented transcripts.

- **Additional argument**

::

    additional arguments:
      --rnafold_path RNAFOLD_PATH
                            Path of RNAfold of the Vienna package
      --relplot_path RELPLOT_PATH
                            Path of relplot.pl of the Vienna package.
      --mountain_path MOUNTAIN_PATH
                            Path of mountain.pl of the Vienna package.
      --blastn_path BLASTN_PATH
                            Path of blastn of the BLAST+ package.
      --blastx_path BLASTX_PATH
                            Path of blastx of the BLAST+ package.
      --makeblastdb_path MAKEBLASTDB_PATH
                            Path of makeblastdb of the BLAST+ package.
      --processing_site_files PROCESSING_SITE_FILES [PROCESSING_SITE_FILES ...], -p PROCESSING_SITE_FILES [PROCESSING_SITE_FILES ...]
                            Paths of the processing site gff files. It can improve
                            the detection of UTR-derived sRNAs.
      --terminator_files TERMINATOR_FILES [TERMINATOR_FILES ...], -e TERMINATOR_FILES [TERMINATOR_FILES ...]
                            Paths of the terminator gff files.
      --promoter_tables PROMOTER_TABLES [PROMOTER_TABLES ...], -pt PROMOTER_TABLES [PROMOTER_TABLES ...]
                            If promoter tables can be provided, please assign the
                            paths of promoter tables, The format of the table is
                            $GENOME $TSS_POSITION $TSS_STRAND $PROMOTER_NAME.
      --promoter_names PROMOTER_NAMES [PROMOTER_NAMES ...], -pn PROMOTER_NAMES [PROMOTER_NAMES ...]
                            If --promoter_tables is provided, please assign the
                            promoter name (the last column of promoter table). For
                            multiple promoters, please put spaces between the
                            promoters. Default is None.
      --sorf_files SORF_FILES [SORF_FILES ...], -O SORF_FILES [SORF_FILES ...]
                            If comparison between sRNAs and sORFs needs to be
                            done, Please assign the paths of sORF gff files
      --parallel_blast PARALLEL_BLAST, -pb PARALLEL_BLAST
                            The number of parallel jobs. Default is 10.
      --tss_source, -ts     The TSS gff files are generated from ANNOgesic or not.
                            Default is True (from ANNOgesic).
      --tss_intergenic_antisense_tolerance TSS_INTERGENIC_ANTISENSE_TOLERANCE, -tit TSS_INTERGENIC_ANTISENSE_TOLERANCE
                            The 5'ends of intergenic and antisense sRNA candidates
                            will be extended or withdrew by this value
                            (nucleotides) for searching the associated TSSs.
                            Default is 3.
      --tss_5utr_tolerance TSS_5UTR_TOLERANCE, -t5 TSS_5UTR_TOLERANCE
                            The 5'ends of 5'UTR-derived sRNAs will be extended or
                            withdrew by this value (nucleotides) for searching the
                            associated TSSs. The input type can be percentage
                            ("p") or the real amount of reads ("n"). For example,
                            p_0.05 means this value is 5 percent of the length of
                            5'UTR. n_10 means this value is 10 nts. Default is
                            n_3.
      --tss_3utr_tolerance TSS_3UTR_TOLERANCE, -t3 TSS_3UTR_TOLERANCE
                            Similar to --tss_5utr_tolerance. This value is for
                            3'UTR-derived sRNAs. Default is p_0.04.
      --tss_intercds_tolerance TSS_INTERCDS_TOLERANCE, -tc TSS_INTERCDS_TOLERANCE
                            Similar to --tss_5utr_tolerance. This value is for
                            interCDS-derived sRNAs. Default is p_0.04.
      --terminator_tolerance_in_srna TERMINATOR_TOLERANCE_IN_SRNA, -eti TERMINATOR_TOLERANCE_IN_SRNA
                            The 3'ends of sRNA candidates will be withdrew by this
                            value (nucleotides) for searching the associated
                            terminators which are within sRNAs. Default is 30.
      --terminator_tolerance_out_srna TERMINATOR_TOLERANCE_OUT_SRNA, -eto TERMINATOR_TOLERANCE_OUT_SRNA
                            The 3'ends of sRNA candidates will be extended by this
                            value (nucleotides) for searching the associated
                            terminators which are behind of sRNAs. Default is 30.
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
                            Similar to --min_intergenic_tex_coverage. This value
                            is for TEX- libraries. Default is 0,0,0,30,10.
      --min_intergenic_fragmented_coverage MIN_INTERGENIC_FRAGMENTED_COVERAGE, -if MIN_INTERGENIC_FRAGMENTED_COVERAGE
                            Similar to --min_intergenic_tex_coverage. This value
                            is for fragmented (or conventional) libraries. Default
                            is 400,200,0,50,20.
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
                            Similar to --min_intergenic_tex_coverage. This value
                            is for antisense in stead of intergenic. Default is
                            0,0,0,40,20.
      --min_antisense_notex_coverage MIN_ANTISENSE_NOTEX_COVERAGE, -an MIN_ANTISENSE_NOTEX_COVERAGE
                            Similar to --min_intergenic_notex_coverage. This value
                            is for antisense in stead of intergenic. Default is
                            0,0,0,30,10.
      --min_antisense_fragmented_coverage MIN_ANTISENSE_FRAGMENTED_COVERAGE, -af MIN_ANTISENSE_FRAGMENTED_COVERAGE
                            Similar to --min_intergenic_fragmented_coverage. This
                            value is for antisense in stead of intergenic. Default
                            is 400,200,0,50,20.
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
                            Similar to --min_utr_tex_coverage. This value is for
                            TEX- libraries. Default is p_0.7,p_0.5,p_0.6.
      --min_utr_fragmented_coverage MIN_UTR_FRAGMENTED_COVERAGE, -uf MIN_UTR_FRAGMENTED_COVERAGE
                            Similar to --min_utr_tex_coverage. This value is for
                            fragmented (or conventional) libraries. Default is
                            p_0.7,p_0.5,p_0.6.
      --min_all_utr_coverage MIN_ALL_UTR_COVERAGE, -mu MIN_ALL_UTR_COVERAGE
                            The minimum coverage of UTR-derived sRNAs. The
                            coverage of UTR-derived sRNAs should not only exceed
                            the --min_utr_TEX_coverage, --min_utr_noTEX_coverage
                            and --min_utr_fragmented_coverage, but also this
                            value. Default is 50.
      --cutoff_energy CUTOFF_ENERGY, -ce CUTOFF_ENERGY
                            If "sec_str" is included in --filter_info, please
                            assign the maximum folding energy change (normalized
                            by length of gene). Default is -0.05.
      --mountain_plot, -m   Generating mountain plot of sRNA candidate. Default is
                            False.
      --nr_format, -nf      Format nr database. Default is False.
      --srna_format, -sf    Format sRNA database. Default is False.
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
                            Similar to --decrease_intergenic_antisense. This value
                            is for UTR-derived sRNA. Default is 0.05.
      --tolerance_intergenic_antisense TOLERANCE_INTERGENIC_ANTISENSE, -ti TOLERANCE_INTERGENIC_ANTISENSE
                            The 5'ends and 3'ends of intergenic and antisense
                            sRNAs will be extended by this value (nucleotides) for
                            detecting the significant coverage decrease. (please
                            check --decrease_intergenic_antisense). For example,
                            the location of intergenic sRNA is 300-400, and
                            --tolerance_intergenic_antisense is 30. The searching
                            region is 270-430. Default is 10.
      --tolerance_utr TOLERANCE_UTR, -tu TOLERANCE_UTR
                            Similar to --tolerance_intergenic_antisense. This is
                            for UTR-derived sRNAs. Default is 10.
      --cutoff_nr_hit CUTOFF_NR_HIT, -cn CUTOFF_NR_HIT
                            The maximum hits number in nr database. Default is 0.
      --blast_e_nr BLAST_E_NR, -en BLAST_E_NR
                            The maximum e-value for searching in nr database.
                            Default is 0.0001.
      --blast_e_srna BLAST_E_SRNA, -es BLAST_E_SRNA
                            The maximum e-value for searching in sRNA database.
                            Default is 0.0001.
      --blast_score_srna BLAST_SCORE_SRNA, -bs BLAST_SCORE_SRNA
                            The minimum score for searching in sRNA database.
                            Default is None.
      --blast_score_nr BLAST_SCORE_NR, -bn BLAST_SCORE_NR
                            The minimum score for searching in nr database.
                            Default is None.
      --detect_srna_in_cds, -ds
                            Searching sRNA in CDS (e.g. the genome annotation is
                            not correct). More sRNA candidates which overlap with
                            CDS will be detected. Default is False.
      --overlap_percent_cds OVERLAP_PERCENT_CDS, -oc OVERLAP_PERCENT_CDS
                            The maximum ratio of overlapping between CDS and sRNA
                            candidates. It only works if --detect_srna_in_cds is
                            true. Default is 0.5
      --search_poly_u SEARCH_POLY_U, -sp SEARCH_POLY_U
                            The tolerance length for searching poly U tail of
                            sRNA. If this value is assigned by 0, the 3'end of
                            sRNA will not be extended by searching poly U tail.
                            Default is 15.
      --min_u_poly_u MIN_U_POLY_U, -np MIN_U_POLY_U
                            The minimum number of U that poly U tail should
                            contain. Default is 5.
      --mutation_poly_u MUTATION_POLY_U, -mp MUTATION_POLY_U
                            The minimum number of nts which are not U can be
                            tolerated. Default is 2.
      --ignore_hypothetical_protein, -ih
                            For ignoring hypothetical proteins in the genome
                            annotation file. Default is False.
      --ranking_time_promoter RANKING_TIME_PROMOTER, -rp RANKING_TIME_PROMOTER
                            If --promoter_tables is provided, the information of
                            promoter can be use for ranking sRNA candidates. The
                            ranking score is --ranking_time_promoter * average
                            coverage. For example, a sRNA candidate which is
                            associated with a promoter and its average coverage is
                            10. If --ranking_time_promoter is 2, the ranking score
                            will be 20 (2*10). For the candidate which are not
                            associated with a promoter, the
                            --ranking_time_promoter will be 1. Therefore,
                            --ranking_time_promoter can not be smaller than 1.
                            Default is 2.
      --exclude_srna_in_annotation_file, -ea
                            For excluding the sRNAs which are already annotated in
                            --annotation_files. Default is False.
    optional arguments:
      -h, --help            show this help message and exit

- **Output files**

Output files are stored in ``$ANNOgesic/output/sRNAs``. the output folders and files are following:

**sRNA_2d_$GENOME:** The secondary structures of all sRNA candidates.

**sRNA_seq_$GENOME:** The sequences of all sRNA candidates.

**blast_results_and_misc:** Stores the results of blast.

	**nr_blast_$GENOME.txt:** output of BLAST for the nr database.

	**sRNA_blast_$GENOME.txt:** output of BLAST for the sRNA database.

**figs:** Stores the figures about secondary structures of sRNAs.

	**mountain_plots:** Stores mountain plots of the sRNA candidates. Filename is as ``srna10_NC_009839.1_335339_335435_+_mountain.pdf``.
	``srna10``, ``NC_009839.1``, ``335339``, ``335435``, ``+`` are ID of sRNA gff file, genome name, starting point, end point and strand, respectively.

	**sec_plots:** Stores the secondary structure plots of sRNA candidates. 
	Filename of is as ``srna10_NC_009839.1_335339_335435_+_rss.ps``.
	``srna10``, ``NC_009839.1``, ``335339``, ``335435``, ``+`` are ID of sRNA gff file, genome name, starting point, end point and strand, respectively.

	**dot_plots:** Stores the dot plots of sRNA candidates. 
	Filename of dot plot is as ``srna10_NC_009839.1_335339_335435_+_dp.ps``.
	``srna10``, ``NC_009839.1``, ``335339``, ``335435``, ``+`` are ID of sRNA gff file, genome name, starting point, end point and strand, respectively.

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

	**Name:** sRNA Name which is shown in gff file. The sRNA name is generated from the BLAST search.

	**Start:** Starting point of this sRNA.

	**End:** End point of this sRNA.

	**Strand:** Strand of this sRNA.

	**Start_with_TSS/Cleavage_site:** This sRNA starts with which TSS or cleavage site.

	**End_with_cleavage:** If the sRNA ends with a cleavage site, the information of this cleavage site will be showed here.

	**Candidates:** Position of this sRNA.

	**Lib_type:** This sRNA is detected by TEX+/- or fragmented/conventional library.

	**Best_avg_coverage:** Based on coverage of all libraries, The best average coverage of this sRNA will be showed here.

	**Normalized_secondary_energy_change(by_length):** Secondary folding energy change (normalized by length) of this sRNA.

	**sRNA_types:** Shows the sRNA type.

	**Conflict_sORF:** If this sRNA overlaps sORF, the overlapped sORF will be showed here.

	**nr_hit_number:** The hit numbers of this sRNA in nr database.

	**sRNA_hit_number:** The hit numbers of this sRNA in sRNA database.

	**nr_hit_top3|ID|e-value|score:** The top 3 hits of this sRNA in nr database will be showed here. 
	The information includes protein name, ID, e-value, and score.

	**sRNA_hit|e-value|score:** If the homology of this sRNA can be found in sRNA database, the information will be showed here.
	The information includes sRNA name, e-value, and score.

	**Overlap_CDS_forward:** If the sRNA overlaps genomic features, the information of the overlapped features will be showed here (for forward strand).

	**Overlap_nts_forward:** If the sRNA overlaps genomic features, the length and percentage of the overlapping regions will be showed here (for forward strand).

	**Overlap_CDS_reverse:** If the sRNA overlaps genomic features, the information of the overlapped features will be showed here (for reverse strand).

	**Overlap_nts_reverse:** If the sRNA overlaps genomic features, the length and percentage of the overlapping regions will be showed here (for reverse strand).

	**End_with_terminator:** The terminator which is associated with this sRNA.

	**Associated_promoter:** The promoter which is associated with this sRNA.

	**sRNA_length:** sRNA length.

	**Avg_coverage:$LIB_NAME:** Shows the average coverage information of the libraries about this sRNA.

**gffs:** Stores gff files of the sRNA. There are also some sub-folders:

	**for_classes:** Stores the results based on different sRNA classes.

	**best_candidates:** Stores the best results of sRNAs after filtering.

	**all_candidates:** Stores all candidates without filtering.

Some useful information can be found in the tags of the attributes within the sRNA gff file.
Based on this information, we can know the details of the specific sRNA. The tags are as following:

	**Name:** The sRNA name is generated from the BLAST search.

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

**Gff files of the genome annotations containing CDSs, tRNAs, rRNAs, etc**
**Gff files of the transcripts**

**Wiggle files of TEX+/- or fragmented/conventional libraries:** Please refer to the section :ref:`The input format of libraries for running ANNOgesic`.

**fasta files of the genome sequences**

- **Optional input files**

**Gff files of the TSSs:** For checking the sORFs start from TSS or not. We strongly recommend to input this file. 

**Gff files of sRNAs:** For checking the overlap of sRNAs and sORFs.

- **Basic arguments**

::

    usage: annogesic sorf --project_path PROJECT_PATH [--utr_derived_sorf]
                          --fasta_files FASTA_FILES [FASTA_FILES ...]
                          --transcript_files TRANSCRIPT_FILES
                          [TRANSCRIPT_FILES ...] --annotation_files
                          ANNOTATION_FILES [ANNOTATION_FILES ...]
                          [--tss_files TSS_FILES [TSS_FILES ...]]
                          [--tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]]
                          [--frag_libs FRAG_LIBS [FRAG_LIBS ...]]
                          [--tex_notex {1,2}]
                          [--replicate_tex REPLICATE_TEX [REPLICATE_TEX ...]]
                          [--replicate_frag REPLICATE_FRAG [REPLICATE_FRAG ...]]
                          [--additional arguments]
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --utr_derived_sorf, -u
                            Detect UTR-derived sORF. Default is False.
      --fasta_files FASTA_FILES [FASTA_FILES ...], -f FASTA_FILES [FASTA_FILES ...]
                            Paths of the fasta files of the reference genome.
      --transcript_files TRANSCRIPT_FILES [TRANSCRIPT_FILES ...], -a TRANSCRIPT_FILES [TRANSCRIPT_FILES ...]
                            Paths of the transcript gff files.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the the genome annotation gff files
                            containing CDSs, tRNAs, rRNAs, etc. Gff files of
                            transcripts, sRNAs, and TSSs need to be separately
                            assigned to --transcript_files, --srna_files, and
                            --tss_files, respectively.
      --tss_files TSS_FILES [TSS_FILES ...], -t TSS_FILES [TSS_FILES ...]
                            Paths of TSS gff files.
      --tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...], -tl TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]
                            TEX+/- wig files. The format is:
                            wig_file_path:TEX+/-(tex or notex):condition_id(intege
                            r):replicate_id(alphabet):strand(+ or -). If multiple
                            wig files need to be assigned, please use spaces to
                            separate the wig files. For example,
                            my_lib_tex_forward.wig:tex:1:a:+
                            my_lib_tex_reverse.wig:tex:1:a:-.
      --frag_libs FRAG_LIBS [FRAG_LIBS ...], -fl FRAG_LIBS [FRAG_LIBS ...]
                            Wig files of RNA-Seq with fragmented transcripts. The
                            format is: wig_file_path:frag:condition_id(integer):re
                            plicate_id(alphabet):strand(+ or -). If multiple wig
                            files need to be assigned, please use spaces to
                            separate the wig files. For example,
                            my_lib_frag_forward.wig:frag:1:a:+
                            my_lib_frag_reverse.wig:frag:1:a:-.
      --tex_notex {1,2}, -te {1,2}
                            If the TEX+/- libraries are provided, this value is
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
                            Similar to --replicates_tex. This value is for
                            fragmented (or conventional) libraries.

- **Additional arguments**

::

    additional arguments:
      --srna_files SRNA_FILES [SRNA_FILES ...], -s SRNA_FILES [SRNA_FILES ...]
                            Paths of the sRNA gff files for comparing sORF and
                            sRNA to detect the overlapping.
      --utr_length UTR_LENGTH, -ul UTR_LENGTH
                            The utr length for comparing TSS with sORF. The
                            default number is 300.
      --min_length MIN_LENGTH, -ml MIN_LENGTH
                            The minimum nucleotide length of sORF. Default is 30.
      --max_length MAX_LENGTH, -Ml MAX_LENGTH
                            The maximum nucleotide length of sORF. Default is 300.
      --cutoff_intergenic_coverage CUTOFF_INTERGENIC_COVERAGE, -ci CUTOFF_INTERGENIC_COVERAGE
                            The minimum coverage of intergenic sORF candidates.
                            Default is 10.
      --cutoff_antisense_coverage CUTOFF_ANTISENSE_COVERAGE, -ai CUTOFF_ANTISENSE_COVERAGE
                            The minimum coverage of antisense sORF candidates.
                            Default is 10.
      --cutoff_5utr_coverage CUTOFF_5UTR_COVERAGE, -cu5 CUTOFF_5UTR_COVERAGE
                            The minimum coverage for 5'UTR derived sORF
                            candidates. This value can be assigned by percentile
                            ("p") or the amount of reads ("n"). For example, p_0.5
                            means that the coverage of sORF candidates should be
                            higher than the 50 percentile of all 5'UTR
                            transcripts. n_10 means that the coverage of sORF
                            candidates should be higher than 10 reads. Default is
                            p_0.5.
      --cutoff_3utr_coverage CUTOFF_3UTR_COVERAGE, -cu3 CUTOFF_3UTR_COVERAGE
                            Similar to --cutoff_5utr_coverage. This value is for
                            3'UTRs. Default is p_0.5.
      --cutoff_intercds_coverage CUTOFF_INTERCDS_COVERAGE, -cuf CUTOFF_INTERCDS_COVERAGE
                            Similar to --cutoff_5utr_coverage. This value is for
                            interCDS. Default is p_0.5.
      --cutoff_base_coverage CUTOFF_BASE_COVERAGE, -cub CUTOFF_BASE_COVERAGE
                            The general minimum coverage of all sORF candidates.
                            All candidates should exceed this value. Default is
                            10.
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
                            Include the sORFs which are not associated with
                            ribosome binding site to the high-confidence sORF
                            list. Default is False.
      --rbs_seq RBS_SEQ [RBS_SEQ ...], -rs RBS_SEQ [RBS_SEQ ...]
                            The sequence of ribosome binding site. If multiple
                            sequences need to be assigned, please use space to
                            split them. Default is AGGAGG.
      --tolerance_rbs TOLERANCE_RBS, -tr TOLERANCE_RBS
                            The number of nucleotides of ribosome binding sites
                            allowed to be different from AGGAGG. Default is 2.
      --tolerance_3end TOLERANCE_3END, -t3 TOLERANCE_3END
                            The number of nucleotides can be extended from the end
                            of transcript for searching stop codon. Default is 30.
      --tolerance_5end TOLERANCE_5END, -t5 TOLERANCE_5END
                            The number of nucleotides can be extended from the
                            starting point of transcript for searching start
                            codon. Default is 5.
      --print_all_combination, -pa
                            For printing all combinations of multiple start and
                            stop codons. Default is False.
      --best_no_srna, -bs   Excluding the sORFs which overlap with sRNAs to highly
                            confidence sORF list. Default is False.
      --best_no_tss, -bt    Excluding the sORFs which do not start with TSS to
                            highly confidence sORF list. Default is False.
      --ignore_hypothetical_protein IGNORE_HYPOTHETICAL_PROTEIN, -ih IGNORE_HYPOTHETICAL_PROTEIN
                            For ignoring hypothetical protein in the genome
                            annotation file. Default is False.

    optional arguments:
      -h, --help            show this help message and exit

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

	**Conflict_sRNA:** If this sORF overlaps sRNA, the overlapped sRNA will be showed here.

	**Frame_shift:** If there are sORF candidates which can be found by frame shift, 
	the number of frame shift will be showed here. ``1`` means there 
	are some candidates can be found by frame shift once. ``2`` means there are some candidates can be found by frame shift twice.

	**Lib_type:** This sORF can be detected in TEX+/- or fragmented/conventional libraries.

	**Best_avg_coverage:** Based on coverage of all libraries, The best average coverage of this sORF will be showed here.

	**Seq:** Sequence of this sORF.

	**Avg_coverage:$LIB_NAME:** Shows the average coverage information of the libraries about this sORF.

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

- **Basic arguments**

::

    usage: annogesic promoter [--program {meme,glam2,both}] --project_path
                              PROJECT_PATH --fasta_files FASTA_FILES
                              [FASTA_FILES ...] --tss_files TSS_FILES
                              [TSS_FILES ...]
                              [--motif_width MOTIF_WIDTH [MOTIF_WIDTH ...]]
                              [--num_motifs NUM_MOTIFS]
                              [--nt_before_tss NT_BEFORE_TSS]
                              [--additional arguments]
    
    basic arguments:
      --program {meme,glam2,both}, -p {meme,glam2,both}
                            Please choose the program -- meme, glam2 or both.
                            Default is both
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --fasta_files FASTA_FILES [FASTA_FILES ...], -f FASTA_FILES [FASTA_FILES ...]
                            Paths of the genome fasta files.
      --tss_files TSS_FILES [TSS_FILES ...], -t TSS_FILES [TSS_FILES ...]
                            Paths of the TSS gff files.
      --motif_width MOTIF_WIDTH [MOTIF_WIDTH ...], -w MOTIF_WIDTH [MOTIF_WIDTH ...]
                            Length of the motifs to detects. For a range insert
                            "-" between two values. Moreover, if multiple lengths
                            need to be assigned, please use spaces to separate
                            them. For an example, 50 2-10 means that the lengths
                            of motifs are 50 and within 2 to 10. The number should
                            be less or equal than --nt_before_TSS. Default is 50.
      --num_motifs NUM_MOTIFS, -n NUM_MOTIFS
                            The number of motifs. Default is 10.
      --nt_before_tss NT_BEFORE_TSS, -b NT_BEFORE_TSS
                            The number of nucleotides upstream of the TSSs for
                            promoter prediction. Default is 50.

- **Additional arguments**

::

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
                            The number of parallel jobs.
      --tss_source, -s      The TSS gff files are generated from ANNOgesic or not.
                            Default is True (from ANNOgesic)
      --tex_libs TEX_LIBS [TEX_LIBS ...], -tl TEX_LIBS [TEX_LIBS ...]
                            If --tss_source is False, please assign the name of
                            the TEX+/- library. The format is:
                            wig_file_path:TEX+/-(tex or notex):condition_id(intege
                            r):replicate_id(alphabet):strand(+ or -). If multiple
                            wig files need to be assigned, please use spaces to
                            separate the wig files. For an example,
                            $WIG_PATH_1:tex:1:a:+ $WIG_PATH_2:tex:1:a:-.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            If --tss_source is False, please assign the paths of
                            the genome annotation gff files containing CDSs,
                            tRNAs, rRNAs, etc.
      --combine_all, -c     Generate global promoter motifs across all reference
                            sequences in the TSS files. Default is False.

    optional arguments:
      -h, --help            show this help message and exit

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
	``NC_000915.1``, ``allgenome``, ``primary`` and ``45_nt`` are gff filename, genome name, TSS type and upstream nucleotides of TSS, respectively.
	If genome name is ``allgenome``, this means that the result is generated by the information of all genomes of gff files. 
	If there is only one genome in the gff file, the genome name will be assigned as ``allgenome`` as well. Several files are stored in the sub-folder:
	
		**Figures of the promoter motifs:** Contains EPS and PNG files.
	
		**Details of the promoter motifs:** Contains HTML file, XML file and TXT file. These files include the TSS information.
	
		**Promoter tables:** ``meme.csv`` or ``glam2.csv`` is the promoter table which also includes the TSS information. 
		Moreover, it can used as an input for sRNA detection (``srna``). Please check the section ``srna``.

**TSS_classes:** If the TSSs are not computed by ANNOgesic, ``TSS_classes`` will be generated for classification of TSS.
TSS gff files with TSS types will be stored here.

.. _operon:

operon (Operon detection)
-------------------------

``operon`` will search operons and sub-operons based on TSSs, transcripts, and genes. If the transcripts 
are not associated with genes, they would not be counted as operons.

- **Required files**

**Gff files of the genome annotations containing CDSs, tRNAs, rRNAs, etc**

**Gff files of the transcripts**

- **Optional input files**

**Gff files of the TSSs**: We strongly recommend to input this file for detecting sub-operon.

**Gff files of the terminators**

- **Basic arguments**

::

    usage: annogesic operon --project_path PROJECT_PATH
                            [--tss_files TSS_FILES [TSS_FILES ...]]
                            --annotation_files ANNOTATION_FILES
                            [ANNOTATION_FILES ...] --transcript_files
                            TRANSCRIPT_FILES [TRANSCRIPT_FILES ...]
                            [--terminator_files TERMINATOR_FILES [TERMINATOR_FILES ...]]
                            [--additional arguments]
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --tss_files TSS_FILES [TSS_FILES ...], -t TSS_FILES [TSS_FILES ...]
                            Paths of the TSS gff files.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the genome annotation gff files containing
                            CDSs, tRNAs, rRNAs etc. Gff files of TSSs,
                            transcripts, terminators, need to be separately
                            assigned to --tss_files, --transcript_files,
      --transcript_files TRANSCRIPT_FILES [TRANSCRIPT_FILES ...], -a TRANSCRIPT_FILES [TRANSCRIPT_FILES ...]
                            Paths of the transcript gff files.
      --terminator_files TERMINATOR_FILES [TERMINATOR_FILES ...], -e TERMINATOR_FILES [TERMINATOR_FILES ...]
                            Paths of terminator gff files.

- **Additional arguments**

::

    additional arguments:
      --tss_tolerance TSS_TOLERANCE, -tt TSS_TOLERANCE
                            The 5'ends of transcripts will be extended or withdrew
                            by this value (nucleotides) for searching the
                            associated TSSs. Default is 5.
      --terminator_tolerance TERMINATOR_TOLERANCE, -et TERMINATOR_TOLERANCE
                            The 3'ends of transcripts will be extended or withdrew
                            by this value (nucleotides) for searching the
                            associated terminators. Default is 30.
      --min_length MIN_LENGTH, -l MIN_LENGTH
                            The minimum length of operon. Default is 20.

    optional arguments:
      -h, --help            show this help message and exit

- **Output files**

Output files are stored in ``$ANNOgesic/output/operons``. The output folders are as following:

**gffs:** The gff files of operons. **associated_gene** shows the genes located in the operon.

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

**Gff files of the genome annotations containing CDSs, tRNAs, rRNAs, etc**

- **Basic Arguments**

::

    usage: annogesic circrna --project_path PROJECT_PATH
                             [--read_files READ_FILES [READ_FILES ...]]
                             [--bam_files BAM_FILES [BAM_FILES ...]] --fasta_files
                             FASTA_FILES [FASTA_FILES ...] --annotation_files
                             ANNOTATION_FILES [ANNOTATION_FILES ...]
                             [--additional arguments]
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --read_files READ_FILES [READ_FILES ...], -rp READ_FILES [READ_FILES ...]
                            Paths of read fasta or fastq files. ANNOgesic will map
                            the reads via segemehl (with -S). Required format:
                            $SET_NAME:$READ1,$READ2,... If multiple data sets need
                            to be assigned, please separated them by spaces. The
                            read files compressed by bz2 or gz files can be
                            accepted as well. For using BAM files, please check
                            --bam_files.
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
                            Paths of the genome annotation gff files containing
                            CDSs, tRNAs, rRNAs, etc.

- **Additional arguments**

::

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
                            the starting points of candidates. Default is 0.5.
      --end_ratio END_RATIO, -er END_RATIO
                            Similar to --start_ratio. This value is for the end
                            points of candidates. Default is 0.5.
      --ignore_hypothetical_protein IGNORE_HYPOTHETICAL_PROTEIN, -ih IGNORE_HYPOTHETICAL_PROTEIN
                            For ignoring hypothetical protein in the genome
                            annotation file. Default is False.

    optional arguments:
      -h, --help            show this help message and exit

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

	**Genome:** Genome name.

	**Strand:** Strand of the circular RNA.

	**Start:** Starting point of the circular RNA.

	**End:** End point of the circular RNA.

	**Annotation_overlap:** If there is a genome annotation (like a CDS) which overlap this circular RNA, the overlapped feature will be showed here.

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

**Gff files of the genome annotations containing CDSs, tRNAs, rRNAs, etc**

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
                            Paths of the genome annotation gff files containg
                            CDSs.
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
`IntaRNA <https://github.com/BackofenLab/IntaRNA/>`_.

- **Required files**

**Gff files of the genome annotations containing CDSs**

**Gff files of the sRNAs**

**Fasta files of the genomes**

- **Basic arguments**

::

    usage: annogesic srna_target --project_path PROJECT_PATH
                                 --annotation_files ANNOTATION_FILES
                                 [ANNOTATION_FILES ...] --fasta_files FASTA_FILES
                                 [FASTA_FILES ...] --srna_files SRNA_FILES
                                 [SRNA_FILES ...]
                                 [--query_srnas QUERY_SRNAS [QUERY_SRNAS ...]]
                                 --program {RNAplex,RNAup,IntaRNA}
                                 [{RNAplex,RNAup,IntaRNA} ...] [--top TOP]
                                 [--additional arguments]
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the genome annotation gff files containing
                            CDSs.
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
      --program {RNAplex,RNAup,IntaRNA} [{RNAplex,RNAup,IntaRNA} ...], -p {RNAplex,RNAup,IntaRNA} [{RNAplex,RNAup,IntaRNA} ...]
                            The program for detecting sRNA-mRNA interaction.
                            Please choose "RNAplex", "RNAup" or "IntaRNA". If
                            multiple programs need to be executed, please use
                            space to separate them.
      --top TOP, -t TOP     The ranking number of targets which will be included
                            to final output. The ranking is based on the binding
                            energy. Default is 50.

- **Additional arguments**

::
    
    additional arguments:
      --rnaplfold_path RNAPLFOLD_PATH
                            Path of RNAplfold of the Vienna package.
      --rnaplex_path RNAPLEX_PATH
                            Path of RNAplex of the Vienna package.
      --rnaup_path RNAUP_PATH
                            Path of RNAup of the Vienna package.
      --intarna_path INTARNA_PATH
                            Path of IntaRNA.
      --interaction_length INTERACTION_LENGTH, -i INTERACTION_LENGTH
                            Maximum length of an interaction. Default is 30.
      --window_size_target_rnaplex WINDOW_SIZE_TARGET_RNAPLEX, -wt WINDOW_SIZE_TARGET_RNAPLEX
                            The average of the pair probabilities over windows for
                            mRNA target. It is only applied for "RNAplex". Default
                            is 240.
      --span_target_rnaplex SPAN_TARGET_RNAPLEX, -st SPAN_TARGET_RNAPLEX
                            The maximum allowed separation of a base pair to span
                            for mRNA target. It is only applied for "RNAplex".
                            Default is 160.
      --window_size_srna_rnaplfold WINDOW_SIZE_SRNA_RNAPLFOLD, -ws WINDOW_SIZE_SRNA_RNAPLFOLD
                            Similar to --window_size_target, but for sRNA. Default
                            is 30.
      --span_srna_rnaplfold SPAN_SRNA_RNAPLFOLD, -ss SPAN_SRNA_RNAPLFOLD
                            Similar to --span_target, but for sRNA. Default is 30.
      --unstructured_region_rnaplex_target UNSTRUCTURED_REGION_RNAPLEX_TARGET, -ut UNSTRUCTURED_REGION_RNAPLEX_TARGET
                            Calculate the mean probability of the unpaired region
                            for mRNA target. It only works for "RNAplex". Default
                            is 30.
      --unstructured_region_rnaplex_srna UNSTRUCTURED_REGION_RNAPLEX_SRNA, -us UNSTRUCTURED_REGION_RNAPLEX_SRNA
                            Similar to --unstructured_region_rnaplex_target, but
                            for sRNA. Default is 30.
      --unstructured_region_rnaup UNSTRUCTURED_REGION_RNAUP, -uu UNSTRUCTURED_REGION_RNAUP
                            Compute the mean probability of unpaired region. It
                            only works for "RNAup". Default is 40.
      --energy_threshold_rnaplex ENERGY_THRESHOLD_RNAPLEX, -e ENERGY_THRESHOLD_RNAPLEX
                            The minimum energy for a duplex. It only works for
                            "RNAplex". Default is -8.
      --duplex_distance_rnaplex DUPLEX_DISTANCE_RNAPLEX, -d DUPLEX_DISTANCE_RNAPLEX
                            Distance between target 3'ends of two consecutive
                            duplexes. It works for "RNAplex". Default is 20.
      --parallels_rnaplex PARALLELS_RNAPLEX, -pp PARALLELS_RNAPLEX
                            The number of parallel jobs for running RNAplex.
                            Default is 5.
      --parallels_rnaup PARALLELS_RNAUP, -pu PARALLELS_RNAUP
                            The number of parallel jobs for running RNAup. Default
                            is 20.
      --parallels_intarna PARALLELS_INTARNA, -pi PARALLELS_INTARNA
                            The number of parallel jobs for running IntaRNA.
                            Default is 10.
      --continue_rnaup, -cr
                            For running RNAup based on the previous intermediate
                            results if the previous process stopped. Default is
                            False.
      --slide_window_size_srna_intarna SLIDE_WINDOW_SIZE_SRNA_INTARNA, -sw SLIDE_WINDOW_SIZE_SRNA_INTARNA
                            The silding window size of sRNA sequences. 0 will use
                            the full sequence to execute IntaRNA. Default is 150.
      --max_loop_length_srna_intarna MAX_LOOP_LENGTH_SRNA_INTARNA, -ls MAX_LOOP_LENGTH_SRNA_INTARNA
                            The maximal loop length of sRNA. If the value is
                            assigned by 0, --slide_window_size_srna_intarna will
                            be used for the maximal loop length of sRNA. Default
                            is 100.
      --slide_window_size_target_intarna SLIDE_WINDOW_SIZE_TARGET_INTARNA, -tw SLIDE_WINDOW_SIZE_TARGET_INTARNA
                            The silding window size of target sequences. 0 will
                            use the full sequence to execute IntaRNA. Default is
                            150.
      --max_loop_length_target_intarna MAX_LOOP_LENGTH_TARGET_INTARNA, -lt MAX_LOOP_LENGTH_TARGET_INTARNA
                            The maximal loop length of target. If the value is
                            assigned by 0, --slide_window_size_target_intarna will
                            be used for the maximal loop length of target. Default
                            is 100.
      --mode_intarna {H,E,M}, -mi {H,E,M}
                            The prediction mode of IntaRNA. 'H' is heuristic, 'M'
                            is exact with long computational time. 'E' is exact
                            with long computational time and high memory.
      --potential_target_start POTENTIAL_TARGET_START, -ps POTENTIAL_TARGET_START
                            Distance for the extraction of upstream nucleotides of
                            --target_feature. Default is 200.
      --potential_target_end POTENTIAL_TARGET_END, -pe POTENTIAL_TARGET_END
                            Distance for the extraction of downstream nucleotides
                            of --target_feature. Default is 150.
      --target_feature TARGET_FEATURE [TARGET_FEATURE ...], -tf TARGET_FEATURE [TARGET_FEATURE ...]
                            The feature name of potential targets. If multiple
                            features need to be assigned, please use spaces to
                            separate them. For example, CDS exon. Default is CDS.

    optional arguments:
      -h, --help            show this help message and exit

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

**IntaRNA_results:** Stored all results of IntaRNA. ``$GENOME_IntaRNA.txt`` is raw results of IntaRNA.
``$GENOME_RNAup_rank.csv`` is the tables with details, and the targets are
sorted by binding energy. The meaning of each column is similar to the table of RNAplex.

**merged_results:** Store the results which are merged by the results of ``RNAplex_results``, ``RNAup_results``, 
and ``IntaRNA_results``. 
``$GENOME_merge.csv`` contains all candidates of the all assigned programs. 
``$GENOME_overlap.csv`` contains the results which are top 50 (default) in the all assigned methods. 
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

**Gff files of the genome annotations containing CDSs**

- **Basic arguments**

::

    usage: annogesic ppi_network --project_path PROJECT_PATH
                                 --annotation_files ANNOTATION_FILES
                                 [ANNOTATION_FILES ...] --species_string
                                 SPECIES_STRING --query_strains QUERY_STRAINS
                                 [QUERY_STRAINS ...] [--query QUERY [QUERY ...]]
                                 [--without_strain_pubmed] [--additional arguments]
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the genome annotation gff files containing
                            genes and CDSs with proper locus_tag items in the
                            attributes.
      --species_string SPECIES_STRING, -d SPECIES_STRING
                            Path of the species table of STRING
                            (species.$VERSION.txt).
      --query_strains QUERY_STRAINS [QUERY_STRAINS ...], -s QUERY_STRAINS [QUERY_STRAINS ...]
                            The name of the input file of the query genomes.
                            Required format: $GFF_FILE:$STRAIN_IN_GFF:$STRAIN_IN_S
                            TRING:$STRAIN_FOR_PUBMED. $GFF_FILE is the name of the
                            gff file, $STRAIN_IN_GFF is the name/ID of the strain
                            in the gff file, $STRAIN_IN_STRING is the strain name
                            in species table of STRING (species.$VERSION.txt), and
                            $STRAIN_FOR_PUBMED is the strain name for searching in
                            Pubmed. If the strain is not available in STRING
                            database, it can be relaced by a related strain. For
                            example, Staphylococcus_aureus_HG003.gff:Staphylococcu
                            s_aureus_HG003:"Staphylococcus aureus NCTC
                            8325":"Staphylococcus aureus" (Staphylococcus aureus
                            NCTC 8325 is a related strain of HG003 since HG003 is
                            not available in STRING). If the assigned name is with
                            spaces, please use double quotes. For assigning
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

- **Additional arguments**

::

    additional arguments:
      --score SCORE, -ps SCORE
                            The minimum PIE score for searching literature. The
                            value is from -1 (worst) to 1 (best). Default is 0.
      --node_size NODE_SIZE, -ns NODE_SIZE
                            The size of nodes in figure, default is 4000.

    optional arguments:
      -h, --help            show this help message and exit

- **Output files**

Output files are stored in ``$ANNOgesic/output/PPI_networks``. The output folders are as following:

**best_results:** Stores the results which the scores of `PIE <http://www.ncbi.nlm.nih.gov/CBBresearch/Wilbur/IRET/PIE/>`_
for supported literature are higher than ``--score``.

**all_results:** Stores the results of all protein-protein interactions
(including the low score(`PIE <http://www.ncbi.nlm.nih.gov/CBBresearch/Wilbur/IRET/PIE/>`_) literature).

Under ``best_results`` and ``all_results``, several files and folders are generated:

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

**Gff files of the genome annotations containing CDSs**

**Fasta files of the genome sequences**

- **Optional input files**

**Gff files of the transcripts:** For detecting subcellular localization only based on expressed CDSs.

- **Basic arguments**

::

    usage: annogesic localization --project_path PROJECT_PATH
                                  --annotation_files ANNOTATION_FILES
                                  [ANNOTATION_FILES ...] --fasta_files FASTA_FILES
                                  [FASTA_FILES ...]
                                  [--transcript_files TRANSCRIPT_FILES [TRANSCRIPT_FILES ...]]
                                  --bacteria_type {positive,negative}
                                  [--additional arguments]
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of genome annotation gff files containing CDSs.
      --fasta_files FASTA_FILES [FASTA_FILES ...], -f FASTA_FILES [FASTA_FILES ...]
                            Paths of genome fasta files.
      --transcript_files TRANSCRIPT_FILES [TRANSCRIPT_FILES ...], -a TRANSCRIPT_FILES [TRANSCRIPT_FILES ...]
                            Paths of the transcript gff files for detecting the
                            subcellular localization based on expressed CDS and
                            all CDS.
      --bacteria_type {positive,negative}, -b {positive,negative}
                            The type of bacteria (Gram-positive or Gram-negative).
                            Please assign 'positive' or 'negative'.

- **Additional arguments**

::

    additional arguments:
      --psortb_path PSORTB_PATH
                            Path of Psortb.
      --difference_multi DIFFERENCE_MULTI, -d DIFFERENCE_MULTI
                            For the protein which have multiple predicted
                            locations, if the difference of psorb scores is
                            smaller than this value, all locations will be printed
                            out. Default is 0.5. The maximum value is 10.

    optional arguments:
      -h, --help            show this help message and exit

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

**Gff files of the genome annotations containing CDSs, tRNAs, rRNAs, etc**

**Gff files of the transcripts**


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

- **Optional input files**

**Gff files of the TSSs**: For checking the ribosome binding site. We strongly recommend to input this file.

- **Basic arguments**

::

    usage: annogesic riboswitch_thermometer --project_path PROJECT_PATH
                                            [--program {riboswitch,thermometer,both}]
                                            [--riboswitch_id_file RIBOSWITCH_ID_FILE]
                                            [--rna_thermometer_id_file RNA_THERMOMETER_ID_FILE]
                                            --annotation_files ANNOTATION_FILES
                                            [ANNOTATION_FILES ...]
                                            [--tss_files TSS_FILES [TSS_FILES ...]]
                                            --transcript_files TRANSCRIPT_FILES
                                            [TRANSCRIPT_FILES ...] --fasta_files
                                            FASTA_FILES [FASTA_FILES ...]
                                            --rfam_path RFAM_PATH
                                            [--additional arguments]
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --program {riboswitch,thermometer,both}, -p {riboswitch,thermometer,both}
                            Please choose the feature for the detection. The
                            options can be "riboswitch", "thermometer", "both".
                            Default is both.
      --riboswitch_id_file RIBOSWITCH_ID_FILE, -ri RIBOSWITCH_ID_FILE
                            Path of the file which contains the information of
                            riboswitches in Rfam. Required format of the file:
                            $RFAM_ID{tab}$RIBOSWITCH_NAME{tab}$DESCRIPTION. Please
                            check an example in https://github.com/Sung-Huan/ANNOg
                            esic/blob/master/database/Rfam_riboswitch_ID.csv
      --rna_thermometer_id_file RNA_THERMOMETER_ID_FILE, -ti RNA_THERMOMETER_ID_FILE
                            Same format as for -riboswitch_id_file, but for RNA
                            thermometers. Please check an example in
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

- **Additional arguments**

::

    additional arguments:
      --cmscan_path CMSCAN_PATH, -cs CMSCAN_PATH
                            Path of cmscan in Infernal package.
      --cmpress_path CMPRESS_PATH, -cp CMPRESS_PATH
                            Path of cmpress in Infernal package.
      --utr_length UTR_LENGTH, -u UTR_LENGTH
                            The UTR length. Default is 300.
      --cutoff CUTOFF, -cf CUTOFF
                            The cutoff of the infernal search. The cutoff can be
                            assigned by e value (assigned by 'e') or score
                            (assigned by 's'). For example, 'e_0.001' represents
                            using e value as a cutoff and the maximum value is
                            0.001. 's_8' represents using score as a cutoff and
                            the minimum score is 8. Default is e_0.001.
      --output_all, -o      One query sequence may fit multiple riboswitches or
                            RNA thermometers. It can print multiple riboswitches
                            or RNA thermometers. Otherwise, only the highest
                            confident one will be printed. Default is False.
      --tolerance TOLERANCE, -to TOLERANCE
                            The 5'ends and 3'ends of potential riboswitches or RNA
                            thermometers will be extended by this value
                            (nucleotides) for extracting the sequences to search
                            in Rfam. Default is 10.
      --without_rbs, -wr    Running the prediction without considering ribosome
                            binding site.Default is False.
      --rbs_seq RBS_SEQ [RBS_SEQ ...], -rs RBS_SEQ [RBS_SEQ ...]
                            The sequences of ribosome binding site. If the
                            multiple sequences needs to be assigned, please use
                            space to split them. Default is AGGAGG.
      --tolerance_rbs TOLERANCE_RBS, -tr TOLERANCE_RBS
                            The number of nucleotides of ribosome binding site
                            allow to be different with AGGAGG. Default is 2.

    optional arguments:
      -h, --help            show this help message and exit

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

	**score:** Score of searching this riboswitch/RNA_thermometer to Rfam.

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

	**Score:** Score of searching this riboswitch/RNA_thermometer to Rfam.

	**Start_align:** Position of this riboswich/RNA_thermometer can be aligned to the genome.

	**End_align:** Position this riboswich/RNA_thermometer can be aligned to the genome.

**statistics:** Stores the file which contains the riboswich/RNA_thermometer with corresponding amount.

.. _crispr:

crispr (CRISPR detection)
-------------------------
``crispr`` integrates CRISPR Recognition Tool (`CRT <http://www.room220.com/crt/>`_) which can detect the repeat 
units and spacers of CRISPR. Moreover, the false positive can be removed by comparing candidates with genome annotations.

- **Required tools**

`CRT <http://www.room220.com/crt/>`_.

- **Required files**

**Fasta files of the genome sequences**

- **Optional input files**

**Gff files of the genome annotations containing CDSs, tRNAs, rRNAs, etc:** This file can be used for removing false positive.

- **Basic arguments**

::

    usage: annogesic crispr --project_path PROJECT_PATH --fasta_files
                            FASTA_FILES [FASTA_FILES ...]
                            [--annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...]]
                            [--additional arguments]
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --fasta_files FASTA_FILES [FASTA_FILES ...], -f FASTA_FILES [FASTA_FILES ...]
                            Paths of the genome fasta files.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the genome gff files containing CDSs for
                            comparing CRISPRs and the genome annotation to remove
                            the false positives. Default is None.

- **Additional arguments**

::

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
                            To ignore hypothetical proteins. Default is inactive.

    optional arguments:
      -h, --help            show this help message and exit

- **Output files**

Output files are stored in ``$ANNOgesic/output/crisprs``. The folders which are generated by the subcommand are as following:

**CRT_results:** Stores the output of `CRT <http://www.room220.com/crt/>`_.

**gffs:** Stores CRSIPR gff files. **Please be aware that CRISPR has no strand information (shows '.' in gff files)**. Two sub-folders are under this folder:

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

**Gff files of the genome annotations containing CDSs, tRNAs, rRNAs, etc**

**Gff files of the manual-detected TSSs**

- **Basic arguments**

::

    usage: annogesic optimize_tss_ps --project_path PROJECT_PATH
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
                                     [--additional arguments]
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --program {TSS,PS}, -p {TSS,PS}
                            The feature for optimization. Please assign "TSS" or
                            "PS". Default is TSS.
      --fasta_files FASTA_FILES [FASTA_FILES ...], -f FASTA_FILES [FASTA_FILES ...]
                            Paths of the fasta file of the reference genome.
      --annotation_files ANNOTATION_FILES [ANNOTATION_FILES ...], -g ANNOTATION_FILES [ANNOTATION_FILES ...]
                            Paths of the genome annotation gff file containing
                            CDSs, tRNAs, rRNAs, etc.
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

- **Additional arguments**

::

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
                            Number of parallel threats for optimization. Default
                            is 4.
      --steps STEPS, -s STEPS
                            Number of tota runs for the optimization. Default is
                            4000 runs.

    optional arguments:
      -h, --help            show this help message and exit

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

- **Basic arguments**

::

    usage: annogesic screenshot --project_path PROJECT_PATH --fasta_file
                                FASTA_FILE --main_gff MAIN_GFF
                                [--side_gffs SIDE_GFFS [SIDE_GFFS ...]]
                                [--tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]]
                                [--frag_libs FRAG_LIBS [FRAG_LIBS ...]]
                                --output_folder OUTPUT_FOLDER
                                [--additional arguments]
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --fasta_file FASTA_FILE, -f FASTA_FILE
                            Path of the genome fasta file.
      --main_gff MAIN_GFF, -mg MAIN_GFF
                            Screenshots will be generated based on the positions
                            of genomic features in this gff file.
      --side_gffs SIDE_GFFS [SIDE_GFFS ...], -sg SIDE_GFFS [SIDE_GFFS ...]
                            The gff files of further genomic features (besides
                            --main_gff). For assigning multiple files, please use
                            spaces to separated them.
      --tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...], -tl TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]
                            TEX+/- wig files. The format is:
                            wig_file_path:TEX+/-(tex or notex):condition_id(intege
                            r):replicate_id(alphabet):strand(+ or -). If multiple
                            wig files need to be assigned, please use spaces to
                            separate the wig files.
      --frag_libs FRAG_LIBS [FRAG_LIBS ...], -fl FRAG_LIBS [FRAG_LIBS ...]
                            Wig files of RNA-Seq of fragmented transcripts. The
                            format is: wig_file_path:frag:condition_id(integer):re
                            plicate_id(alphabet):strand(+ or -). If multiple wig
                            files need to be assigned, please use spaces to
                            separate the wig files. For example,
                            my_lib_frag_forward.wig:frag:1:a:+
                            my_lib_frag_reverse.wig:frag:1:a:-.
      --output_folder OUTPUT_FOLDER, -o OUTPUT_FOLDER
                            Path of the output folder. A sub-folder "screenshots"
                            in the --output_folder will be created to store the
                            results.

- **Additional arguments**

::

    additional arguments:
      --height HEIGHT, -he HEIGHT
                            Height of the screenshot. Default is 1500.
      --present {expand,collapse,squish}, -p {expand,collapse,squish}
                            The presentation types (expand, collapse, or squish)
                            of the features in the screenshot. Default is expand.

    optional arguments:
      -h, --help            show this help message and exit

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

- **Basic arguments**

::

    usage: annogesic colorize_screenshot_tracks --project_path PROJECT_PATH
                                                --screenshot_folder
                                                SCREENSHOT_FOLDER --track_number
                                                TRACK_NUMBER
                                                [--additional arguments]
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --screenshot_folder SCREENSHOT_FOLDER, -f SCREENSHOT_FOLDER
                            The folder containing "screenshots" which are
                            generated from the subcommand "screenshot".
      --track_number TRACK_NUMBER, -t TRACK_NUMBER
                            The number of tracks.

- **Additional arguments**

::

    additional arguments:
      --imagemagick_covert_path IMAGEMAGICK_COVERT_PATH, -m IMAGEMAGICK_COVERT_PATH
                            Path of "convert" in the ImageMagick package.

    optional arguments:
      -h, --help            show this help message and exit

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
Moreover, if the genomic features from different gff files are overlapped, the module can select the data based on 
user-assigned source or keep all overlapping features to the final output.

- **Basic arguments**

::

    usage: annogesic merge_features --project_path PROJECT_PATH
                                    --output_prefix OUTPUT_PREFIX
                                    [--transcript_file TRANSCRIPT_FILE]
                                    [--other_features_files OTHER_FEATURES_FILES [OTHER_FEATURES_FILES ...]]
                                    [--source_for_overlapping SOURCE_FOR_OVERLAPPING [SOURCE_FOR_OVERLAPPING ...]]
                                    [--additional arguments]
    
    basic arguments:
      --project_path PROJECT_PATH, -pj PROJECT_PATH
                            Path of the project folder.
      --output_prefix OUTPUT_PREFIX, -op OUTPUT_PREFIX
                            The prefix name of output gff file. The file name will
                            be $OUTPUT_PREFIX_merge_features.gff.
      --transcript_file TRANSCRIPT_FILE, -a TRANSCRIPT_FILE
                            Path of transcript gff file. The parent transcripts
                            ("Parent" in gff attributes) of all features will be
                            generated.
      --other_features_files OTHER_FEATURES_FILES [OTHER_FEATURES_FILES ...], -of OTHER_FEATURES_FILES [OTHER_FEATURES_FILES ...]
                            Paths of the gff files (besides transcript gff file).
                            For assigning multiple gff files, please use spaces to
                            separate them.
      --source_for_overlapping SOURCE_FOR_OVERLAPPING [SOURCE_FOR_OVERLAPPING ...], -s SOURCE_FOR_OVERLAPPING [SOURCE_FOR_OVERLAPPING ...]
                            If the locations of features are overlapped, the
                            module will only keep the feature positions which
                            provided from --source_for_overlapping. For example,
                            if --source_for_overlapping is 'Ref' 'ANNOgesic', the
                            overlapping features from these two sources will be
                            both kept. However, if --source_for_overlapping is
                            'Ref', only the overlapping features from 'RefSeq will
                            be kept.' The value of --source_for_overlapping should
                            be the same as the second column of the gff file. if
                            all the sources of the overlapping features can not be
                            found in --source_for_overlapping, all the overlapping
                            features will be kept. Default is keeping all overlap
                            features.

- **Additional arguments**

::

    additional arguments:
      --terminator_tolerance TERMINATOR_TOLERANCE, -et TERMINATOR_TOLERANCE
                            The 3'ends of transcripts will be extended or withdrew
                            by this value (nucleotides) for searching the
                            associated terminators. Default is 30.
      --tss_tolerance TSS_TOLERANCE, -tt TSS_TOLERANCE
                            The 5'ends of transcripts will be extended or withdrew
                            by this value (nucleotides) for searching the
                            associated TSSs. Default is 5.

    optional arguments:
      -h, --help            show this help message and exit

- **Output files**

Output gff files are stored in ``$ANNOgesic/output/merge_all_features``. The tag - ``Parent`` in the attributes of 
gff file shows the parent transcript.
