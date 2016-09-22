ANNOgesic's subcommands
===============

In general, the subcommands need at least one argument - the analysis
folder. If this is not given, ANNOgesic assumes that the current
folder is the analysis folder.

The format of filename
--------------------
In order to recognize the file types and the relationship between different files,
please follow the principle of the format of filename.

Keep the filename of genome fasta file the same as annotation file. For example,
``NC_007795.fa, NC_007795.gff, NC_007795.ptt, NC_007795.rnt, NC_007705.gbk``.

All the ``.gff`` files which are input of ``ANNOgesic``, please follow the format.
``$STRAINNAME_$FEATURE.gff``. For example, ``NC_007795_TSS.gff, NC_007795_transcript.gff``.

The names of features are:

===============  ===========================
Feature name     meaning
---------------  --------------------------- 
TSS              transcription starting site
processing       processing site
transcript       transcript assembly
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

Please avoid ``|`` in the filename, strain name of Gff3 files or fasta file.

The input format of libraries for running ANNOgesic
-----------------------------

Some modules of ``ANNOgesic`` need the information of libraries.
please follow the following format:

``$LIBRARY_FILENAME:$LIBRARY_TYPE:$CONDITION:$REPLICATE:$STRAND``

``$LIBRARY_FILENAME`` means the ``.wig`` file.

``$LIBRARY_TYPE`` can be ``tex`` (TEX+) or ``notex`` (TEX-) or ``frag`` (fragmented).

``$CONDITION`` is the index of conditions. Please use 1, 2, 3, ... to represent different conditions.

``$REPLICATE`` is the index of replicates. Please use a, b, c, ... to represent different replicates.

``$STRAND`` is the strand of wiggle file. Please use + or -.

All the libraries should be separated by comma like the above.
The libraries should be listed in one line.

For example, 

for TEX +/- treated libraries:

::

  TSB_OD_0.2_TEX_reverse.wig:tex:1:a:-,
  TSB_OD_0.5_TEX_reverse.wig:tex:2:a:-,
  TSB_OD_0.2_TEX_forward.wig:tex:1:a:+,
  TSB_OD_0.5_TEX_forward.wig:tex:2:a:+,
  TSB_OD_0.2_reverse.wig:notex:1:a:-,
  TSB_OD_0.5_reverse.wig:notex:2:a:-,
  TSB_OD_0.2_forward.wig:notex:1:a:+,
  TSB_OD_0.5_forward.wig:notex:2:a:+,

And for fragmented libraries:

::

  fragmented_forward.wig:frag:1:a:+,fragmented_reverse.wig:frag:1:a:-

If only conventional RNA-seq data without fragementation or TEX treated can be provided, 
it can still be assigned to fragmented libraries and run ANNOgesic.
However, it may influnece on the results.

The format of sRNA database
-----------------------------
If searching sRNA database is assigned to sRNA detection, please follow the format. Or the output will become
nonsense. The format is 

::
  $ID|$STRAIN|$SRNANAME

for example

::
  srn_4840|S._aureus_NCTC8325|RsaOV

You can also download sRNA database `BSRD <http://www.bac-srna.org/BSRD/index.jsp>`_ from our
`Github <https://github.com/Sung-Huan/ANNOgesic/tree/master/database>`_ easily.

Definition of reference strain and target strain
------------------------------
"target strain" represents the strain which user want to annotate.
"reference strain" represents the strain which is close to the "target strain".
If user have no fasta file or genome annotation files of "target strain", 
ANNOgesic can generate them by applying the mutations between "reference strain" and "target strain".

Riboswitch and RNA thermometer dataset of Rfam
----------------------------
For riboswitch and RNA thermometer detection, it needs the information of riboswitch and RNA thermometer in Rfam. 
The input format is the following.

======== ==== ==========================
#Rfam_ID Name Description
-------- ---- --------------------------
RF00162  SAM  SAM riboswitch box leader
RF00059  TPP  TPP riboswitch THI element
======== ===  ==========================

All columns are splited by ``tab``. You can also download the data from our 
`Github <https://github.com/Sung-Huan/ANNOgesic/tree/master/database>`_.

create
-----

``create`` generates the folders for analysis. Once these folders are created, 
please put the required files into the proper folders.

BAMs: For ``.bam`` files. ``BAMs_map_reference`` 
is for the ``.bam`` files which mapped on "reference strain".
``BAMs_map_target`` is for the ``.bam`` files which mapped on "target strain".

database: For all databases

manual_TSS: If the manual detected transcription starting sites(TSSs) is provided,
it can be stored here for running ``TSS_optimization`` or merging 
the automatic predicted ones and manual detected ones. Please use gff3 format.

manual_processing_site: It is similar with ``manual_TSS``, it is for 
processing sites.

mutation_table: If the mutations table between "reference strain" and 
"target strain" is provided, please put the file here. Please refer
to the section of ``get_target_fasta`` for the format of 
mutation table.

reads: For running ``circrna`` with mapping reads by ANNOgesic,
please put the reads here. It can also deal with ``.bzip2`` and ``.gzip``.
       
reference: For annotation files and fasta files of "reference strain". 
If they can be downloaded from NCBI, the files can also be gain via running ``get_input_files``.

riboswitch_ID: For storing the file which contains all the Rfam ID of riboswitch.
For the details of format, please refer to the section of 
``Riboswitch and RNA thermometer dataset of Rfam``.

RNA_thermometer_ID: For storing the file which contains all the Rfam ID of RNA thermometer.
For the details of format, please refer to the section of
``Riboswitch and RNA thermometer dataset of Rfam``.

wigs: For wiggle files. Based on the methods of RNA-Seq, wiggle files can be stored in  
``fragment`` (fragmented libraries) or ``tex_notex`` (TEX +/- treated libraries).


- Arguments

::

  usage: ANNOgesic.py create [-h] project_path
  
  positional arguments:
    project_path  Name/path of the project.
  
  optional arguments:
    -h, --help    show this help message and exit

get_input_files
--------------

``get_input_files`` is the subcommand for downloading required files (fasta, annotation files) from NCBI. 
Therefore, user needs to assign the IP of the reference genome in NCBI. For example,
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000013425.1_ASM1342v1
Then, user can assign the file type for download.


- Pre-reqired information

``FTP source``: The IP of NCBI.

- Arguments


::

    usage: annogesic get_input_files [-h] [--FTP_path FTP_PATH] [--ref_fasta]
                                     [--ref_gff] [--ref_ptt] [--ref_rnt]
                                     [--ref_gbk] [--convert_embl] [--for_target]
                                     [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --FTP_path FTP_PATH, -F FTP_PATH
                            Path of NCBI FTP which can download the required
                            files.
      --ref_fasta, -f       Download fasta files of reference. Default is False.
      --ref_gff, -g         Download gff files of reference. Default is False.
      --ref_ptt, -p         Download ptt files of reference. Default is False.
      --ref_rnt, -r         Download rnt files of reference. Default is False.
      --ref_gbk, -k         Download genbank files of reference. Default is False.
      --convert_embl, -e    Convert gbk to embl files of reference. Default is
                            False.
      --for_target, -t      If the required files of query strain can be
                            downloaded from NCBI (you won't modify the genome),
                            The files can store in target folder in stead of
                            reference folder.

- Output files

The output files will be stored in ``$ANNOgesic_folder/input/reference`` if ``--for_target`` is False.
The output files will be stored in ``$ANNOgesic_folder/output/target`` if ``--for_target`` is True.

``fasta``: Fasta files.

``annotation``: Annotation files.

get_target_fasta
--------------

``get_target_fasta`` is the subcommand for generating fasta files of "target strain" from 
"reference strain". The format of mutation table is following:

==============  =========  ============  ========  =========  ====================  =============  ====  ============
 #reference_id  target_id  reference_nt  position  target_nt  impact_of_correction  locus_tag      gene  Description 
--------------  ---------  ------------  --------  ---------  --------------------  -------------  ----  ------------
 NC_007795.1     HG003     a             333       c                                SAOUHSC_00002  dnaA  XXXXXX      
 NC_007795.1     HG003     t             543       \-          deletion                                  YYYYYY      
 NC_007795.1     HG003     \-            600       g           insertion            SAOUHSC_00132                    
==============  =========  ============  ========  =========  ====================  =============  ====  ============

If the titles of columns is presented on the top, they need to start from ``#``. 
Each column is separated by ``tab``. If the mutation type is deletion or insertion, 
user can put ``-`` to represent them. The information of ``target_id``, ``reference_id``,
``reference_nt``, ``position``, ``target_nt`` is required. The others can be blank. 
However, please still use ``tab`` to separate all blank columns.

If no mutation information is provided, ``SNP_calling`` can be used for detecting mutations. 
(one module of ``ANNOgesic``). Please refer to the section of ``SNP_calling``.

- Pre-required files

Fasta files of reference genome.

Mutation table: it indicates the information of mutations between reference and target strain.

- Arguments

::

    usage: annogesic get_target_fasta [-h] [--ref_fasta_folder REF_FASTA_FOLDER]
                                      [--mutation_table MUTATION_TABLE]
                                      [--output_format OUTPUT_FORMAT]
                                      [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --ref_fasta_folder REF_FASTA_FOLDER, -r REF_FASTA_FOLDER
                            The path of the folder of fasta files.
      --mutation_table MUTATION_TABLE, -m MUTATION_TABLE
                            The path of mutation table.
      --output_format OUTPUT_FORMAT, -o OUTPUT_FORMAT
                            Please assign the output filename and the strain which
                            should be included in the file. For example:
                            FILE1:strain1_and_strain2,FILE2:strain3. FILE1 is a
                            output fasta file which include the information of
                            strain1 and strain2 (import multi-strains to one file
                            should be separated by "_and_".) And FILE2 is for
                            strain3. The files are splitted by comma.

- Output files

Fasta files of target genome will be stored in ``$ANNOgesic_folder/output/target/fasta``.

annotation_transfer
-----------

``annotation transfer`` is the subcommand for transfering the annotation from "reference strain" 
to target "strain". In this subcommand, `RATT <http://www.sanger.ac.uk/resources/software/pagit/>`_ 
is integrated to achieve it. The similarity of "reference strain" and "target strain" should be closed enough or 
it will influence the final results.
Be attation, before running RATT (annotation transfer), 
please run ``source $PAGIT_HOME/sourceme.pagit`` first. it will modify the path for execute RATT. 
If you use Dockerfile to execute ANNOgesic, the path modification can be skipped.

- Pre-required tools and files

`RATT <http://www.sanger.ac.uk/resources/software/pagit/>`_.

The genbank files of reference genome.

The fasta files of reference genome.

The fasta files of target genome.

- Arguments

::

    usage: annogesic annotation_transfer [-h] [--RATT_path RATT_PATH]
                                         [--compare_pair COMPARE_PAIR]
                                         [--element ELEMENT]
                                         [--transfer_type TRANSFER_TYPE]
                                         [--ref_embl REF_EMBL] [--ref_gbk REF_GBK]
                                         [--ref_fasta REF_FASTA]
                                         [--target_fasta TARGET_FASTA]
                                         [--convert_to_gff_rnt_ptt]
                                         [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --RATT_path RATT_PATH
                            Path of the start.ratt.sh file of RATT folder. Default
                            is start.ratt.sh.
      --compare_pair COMPARE_PAIR, -p COMPARE_PAIR
                            Please assign the name of strain pairs. ex.
                            NC_007795:NEW_NC_007795. The reference strain is
                            NC_007795 and the target strain is NEW_NC_007795. the
                            assigned names is the strain name in the fasta file,
                            not the filenames of fasta files. For multiple
                            strains, please use comma to separate the strains.
      --element ELEMENT, -e ELEMENT
                            It will be assigned to the prefix of all output file.
      --transfer_type TRANSFER_TYPE, -t TRANSFER_TYPE
                            The transfer type for running RATT. (For the details,
                            please refer to the manual of RATT.) Default is
                            Strain.
      --ref_embl REF_EMBL, -re REF_EMBL
                            The folder of embl files.
      --ref_gbk REF_GBK, -rg REF_GBK
                            If you have no embl file, you can assign the folder of
                            genbank files. The genbank can be ended by .gbk, .gbff
                            or .gb
      --ref_fasta REF_FASTA, -rf REF_FASTA
                            The folder of reference fasta files.
      --target_fasta TARGET_FASTA, -tf TARGET_FASTA
                            The folder of target fasta files.
      --convert_to_gff_rnt_ptt, -g
                            Convert the annotation to gff, rnt and ptt. Default is
                            False.

- Output files

All the output files from `RATT <http://www.sanger.ac.uk/resources/software/pagit/>`_
will be stored in ``$ANNOgesic_folder/output/annotation_transfer``.

All annotation files (``.gff``, ``.ptt``, ``.rnt``) will be stored in ``$ANNOgesic_folder/output/target/annotation``.

snp
-------

``snp`` can compare the alignment files and fasta files to detect the mutations by 
`Samtools <https://github.com/samtools>`_, `Bcftools <https://github.com/samtools>`_. 
There are multiple programs which can be applied (with BAQ, without BAQ and extend BAQ) and set the filters 
(QUAL, DP, DP4, etc.) to run ``snp``. Moreover, 
It can also be used for generating the fasta file of "target strain".

- Pre-required files and tools:

`Samtools <https://github.com/samtools>`_.

`Bcftools <https://github.com/samtools>`_.

BAM files for fragmented libraries or TEX +/- treated libraries.

Reference or target genome fasta files.

- Arguments

::

    usage: annogesic snp [-h] [--samtools_path SAMTOOLS_PATH]
                         [--bcftools_path BCFTOOLS_PATH] [--bam_type BAM_TYPE]
                         [--program PROGRAM] [--fasta_path FASTA_PATH]
                         [--tex_bam_path TEX_BAM_PATH]
                         [--frag_bam_path FRAG_BAM_PATH] [--quality QUALITY]
                         [--read_depth_range READ_DEPTH_RANGE] [--ploidy PLOIDY]
                         [--RG_tag] [--min_sample_number MIN_SAMPLE_NUMBER]
                         [--caller CALLER] [--DP4_cutoff DP4_CUTOFF]
                         [--indel_fraction INDEL_FRACTION]
                         [--filter_tag_info FILTER_TAG_INFO]
                         [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --samtools_path SAMTOOLS_PATH
                            If you want to assign the path of samtools, please
                            assign here.
      --bcftools_path BCFTOOLS_PATH
                            If you want to assign the path of bcftools, please
                            assign here.
      --bam_type BAM_TYPE, -t BAM_TYPE
                            Please assign the type of BAM. If the BAM files are
                            produced by mapping to reference strain, please keyin
                            "reference". It is for detecting the mutations between
                            refenece strain and target strain, and generate
                            potential target fasta file. If the BAM files are
                            produced by mapping to target strain, please keyin
                            'target'. It is for detecting the mutations of target
                            genome sequence.
      --program PROGRAM, -p PROGRAM
                            Please assign the program for detecting SNP of
                            transcript: 1: calculate with BAQ, 2: calculate
                            without BAQ, 3: calculate with extend BAQ. Multi-
                            programs can be executed at the same time. For
                            example: 1,2,3. Default is 1,2,3.
      --fasta_path FASTA_PATH, -f FASTA_PATH
                            The path of genome fasta folder.
      --tex_bam_path TEX_BAM_PATH, -tw TEX_BAM_PATH
                            The path of tex+/- wig folder.
      --frag_bam_path FRAG_BAM_PATH, -fw FRAG_BAM_PATH
                            The path of fragmented wig folder.
      --quality QUALITY, -q QUALITY
                            The minimum quality which considers a real snp.
                            Default is 40.
      --read_depth_range READ_DEPTH_RANGE, -d READ_DEPTH_RANGE
                            The range of read depth. If the read depth is higher
                            or lower than this range, it will be excluded. The
                            format is $MIN,$MAX. --read_depth_range can be
                            assigned by different types: 1. real number (r), 2.
                            the read depth based on the number of samples (n) or
                            3. times of average read depth (a). For example,
                            n_10,a_2 means the range of read depth should be
                            higher than average 10 reads of --min_sample_number
                            (if --min_sample_number is 2, DP value in output will
                            be higher than 20) and lower than 2 times of average
                            read depth.If r_10,a_2 is assigned, it means that the
                            minimum read depth becomes 10 without considering the
                            number of samples. Default is n_10,a_2.
      --ploidy PLOIDY, -pl PLOIDY
                            haploid or diploid. Default is haploid.
      --RG_tag, -R          It is opposite of --ignore-RG in samtools. (one BAM
                            file includes multi samples). Default is False.
      --min_sample_number MIN_SAMPLE_NUMBER, -ms MIN_SAMPLE_NUMBER
                            The minimum numbers of samples will effect
                            --read_depth_range, --DP4_cutoff and --indel_fraction.
      --caller CALLER, -c CALLER
                            The types of caller - consensus-caller or
                            multiallelic-caller. For details, please check
                            bcftools. "c" represents consensus-caller. "m"
                            represents multiallelic-caller. Default is m.
      --DP4_cutoff DP4_CUTOFF, -D DP4_CUTOFF
                            The cutoff of DP4. DP4 is compose of four numbers:
                            high-quality reference forward bases (number 1),
                            reference reverse bases (number 2), alternate forward
                            bases (number 3) and alternative reverse bases (number
                            4). Two cutoff values can be assigned, ex: n_10,0.8.
                            First cutoff is for (number 3 + number 4). It can be
                            assigned based on 1. real number (r), 2. the read
                            depth based on the number of samples (n) or 3. the
                            times of average read depth (a). The second cutoff is
                            for (number 3 + number 4) / (number 1 + number 2 +
                            number 3 + number 4). These two cutoff is splited by
                            comma. For example, n_10,0.8 means the sum of read
                            depth of number 3 and number 4 should be higher than
                            average 10 reads of --min_sample_number (if
                            --min_sample_number is 2, DP value in output will be
                            higher than 20). And the fraction should be higher
                            than 0.8. If r_10,0.8 is assigned, it means that the
                            sum of read depth of number 3 and number 4 become 10
                            without considering the number of samples. Default is
                            n_10,0.8.
      --indel_fraction INDEL_FRACTION, -if INDEL_FRACTION
                            The fraction of maximum read depth (IMF) and read
                            number of each sample (IDV), which supports insertion
                            of deletion. There are different types can be
                            assigned: 1. real number (r), 2. the read depth based
                            on the number of samples (n) or 3. the times of
                            average read depth (a) for IDV. For example, n_10,0.8
                            means the IDV should be higher than average 10 reads
                            of --min_sample_number (if --min_sample_number, DP
                            value in output will be higher than 20). And IMF
                            should be higher than 0.8. If r_10,0.8 is assigned, it
                            means that IDV become 10 without considering the
                            number of samples. Default is n_10,0.8.
      --filter_tag_info FILTER_TAG_INFO, -ft FILTER_TAG_INFO
                            Please assign 1. the tag, 2. bigger or samller and 3.
                            value for filters. For example, "RPB_b0.1,MQ0F_s0"
                            means RPB should bigger than 0.1 and MQ0F should
                            smaller than 0. Default is
                            RPB_b0.1,MQSB_b0.1,MQB_b0.1,BQB_b0.1.

- Output files

If ``bam_type`` is ``reference``, 
the results will be stored in ``$ANNOgesic/output/SNP_calling/compare_reference``. 
If it is ``target``, the results will be stored in ``$ANNOgesic/output/SNP_calling/validate_target``.

The raw data from `Samtools <https://github.com/samtools>`_ and `Bcftools <https://github.com/samtools>`_
will be stored in ``$ANNOgesic/output/SNP_calling/$BAM_TYPE/SNP_raw_outputs``.

The results will be stored in ``$ANNOgesic/output/SNP_calling/$BAM_TYPE/SNP_table``.

The meaning of filenames are:

``$STRAIN_$PROGRAM_best.vcf`` which is in ``$ANNOgesic/output/SNP_calling/$BAM_TYPE/SNP_table/$STRAIN``. 
It means the results after filtering by cutoff.

``$STRAIN_$PROGRAM.vcf`` which is in ``$ANNOgesic/output/SNP_calling/$BAM_TYPE/SNP_raw_output/$STRAIN``. 
It means the results match the condition of read depth and quality.

``$STRAIN_$PROGRAM_seq_reference.csv`` is the index of fasta files which generated by ``snp``.

For example,

::

  Staphylococcus_aureus_HG003     1632629 .       AaA     AA      57      .
  Staphylococcus_aureus_HG003     1632630 .       aA      a       57      .
  Staphylococcus_aureus_HG003     1499572 .       T       TT,TTTTT        43.8525 .

These mutations will cause conflict. Then, the conflict will effect the positions of other mutations.
Therefore, it will generate four different fasta files. The fisrt two lines are "position conflict", and 
the last line is "mutation conflict".
``$STRAIN_$PROGRAM_seq_reference.csv`` is the index for these four fasta files.

::

   1       1632629 1       1499572:TT      Staphylococcus_aureus_HG003
   1       1632629 2       1499572:TTTTT   Staphylococcus_aureus_HG003
   2       1632630 1       1499572:TT      Staphylococcus_aureus_HG003
   2       1632630 2       1499572:TTTTT   Staphylococcus_aureus_HG003

The first column is the index of position conflict. The second column is the position which be selected.
The third one is the index of mutations conflict. The fourth one is
the position and nucleotides of selected mutation. The last column is the name of strain.
If you refer to ``$ANNOgesic/output/SNP_calling/$BAM_TYPE/seqs``, the filename of fasta is like 
``$FILENAME_$STRIANNAME_$INDEXofPOSITIONCONNFLICT_$INDEXofMUTATIONCONFLICT.fa``. Therefore, the first line of 
``$STRAIN_$PROGRAM_seq_reference.csv`` will generate 
``Staphylococcus_aureus_HG003_Staphylococcus_aureus_HG003_1_1.fa`` 
(if the file name of genome is Staphylococcus_aureus_HG003). The second line will generate
``Staphylococcus_aureus_HG003_Staphylococcus_aureus_HG003_1_2.fa`` and so forth.

The statistics files will be stored in ``$ANNOgesic/output/SNP_calling/$BAM_TYPE/statistics``.

tsspredator(TSS and processing site prediction)
--------------

``tsspredator`` can generate the candidates of TSSs and processing sites via 
`TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_. The parameters 
of `TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_ can be assigned. For  
optimization of parameters of `TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_,
please refer to the section of ``optimize_tsspredator``.

For the information of libraries, please refer to the section 
``The format of libraries for import to ANNOgesic``.

- Pre-required tools and files

`TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_.

The libraries and wiggle files of TEX +/-. Please refer to ``The format of libraries for import to ANNOgesic``.

Fasta file of genome sequence.

Gff file of genome annotation.

If the gff file of manual detected TSSs is provided, ``tsspredator`` can merge the manual one
and predicted one.

If comparing TSSs with transcripts is needed, the gff files of transcripts were required.
For the transcripts, please refer to the section of ``transcript_assembly``.

- Arguments

::

    usage: annogesic tsspredator [-h] [--TSSpredator_path TSSPREDATOR_PATH]
                                 [--fasta_folder FASTA_FOLDER]
                                 [--annotation_folder ANNOTATION_FOLDER]
                                 [--wig_folder WIG_FOLDER] [--height HEIGHT]
                                 [--height_reduction HEIGHT_REDUCTION]
                                 [--factor FACTOR]
                                 [--factor_reduction FACTOR_REDUCTION]
                                 [--enrichment_factor ENRICHMENT_FACTOR]
                                 [--processing_factor PROCESSING_FACTOR]
                                 [--base_height BASE_HEIGHT]
                                 [--replicate_match REPLICATE_MATCH]
                                 [--utr_length UTR_LENGTH] [--lib LIB]
                                 [--output_prefix OUTPUT_PREFIX]
                                 [--merge_manual MERGE_MANUAL] [--statistics]
                                 [--validate_gene]
                                 [--compute_program COMPUTE_PROGRAM]
                                 [--compare_transcript_assembly COMPARE_TRANSCRIPT_ASSEMBLY]
                                 [--fuzzy FUZZY] [--cluster CLUSTER]
                                 [--length LENGTH] [--re_check_orphan]
                                 [--overlap_feature OVERLAP_FEATURE]
                                 [--reference_gff_folder REFERENCE_GFF_FOLDER]
                                 [--remove_low_expression REMOVE_LOW_EXPRESSION]
                                 [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --TSSpredator_path TSSPREDATOR_PATH
                            If you want to assign the path of TSSpredator, please
                            assign here. Default is /usr/local/bin/TSSpredator.jar
      --fasta_folder FASTA_FOLDER, -f FASTA_FOLDER
                            Path of the target genome fasta folder.
      --annotation_folder ANNOTATION_FOLDER, -g ANNOTATION_FOLDER
                            Path of the target genome gff folder.
      --wig_folder WIG_FOLDER, -w WIG_FOLDER
                            The folder of the TEX+/- wig folder.
      --height HEIGHT, -he HEIGHT
                            This value relates to the minimal number of read
                            starts at a certain genomic position to be considered
                            as a TSS candidate. Default is 0.3.
      --height_reduction HEIGHT_REDUCTION, -rh HEIGHT_REDUCTION
                            When comparing different strains/conditions and the
                            step height threshold is reached in at least one
                            strain/condition, the threshold is reduced for the
                            other strains/conditions by the value set here. This
                            value must be smaller than the step height threshold.
                            Default is 0.2.
      --factor FACTOR, -fa FACTOR
                            This is the minimal factor by which the TSS height has
                            to exceed the local expression background. Default is
                            2.0.
      --factor_reduction FACTOR_REDUCTION, -rf FACTOR_REDUCTION
                            When comparing different strains/conditions and the
                            step factor threshold is reached in at least one
                            strain/condition, the threshold is reduced for the
                            other strains/conditions by the value set here. This
                            value must be smaller than the step factor threshold.
                            Default is 0.5.
      --enrichment_factor ENRICHMENT_FACTOR, -ef ENRICHMENT_FACTOR
                            This is the minimal enrichment factor. Default is 2.0.
      --processing_factor PROCESSING_FACTOR, -pf PROCESSING_FACTOR
                            This is the minimal processing factor. If untreated
                            library is higher than the treated library and above
                            which the TSS candidate is considered as a processing
                            site and not annotated as detected. Default is 1.5.
      --base_height BASE_HEIGHT, -bh BASE_HEIGHT
                            This is the minimal number of reads should be mapped
                            on TSS. Default is 0.0.
      --replicate_match REPLICATE_MATCH, -rm REPLICATE_MATCH
                            The TSS candidates should be detected at least in the
                            number of the replicates. The format is
                            $NUMBERofCONDITION_$NUMBERofREPLICATED. For assigning
                            different --replicate_match to different conditions,
                            please use comma to separate them. For example,
                            1_2,2_2,3_3 means number 1 and 2 condition assign 2 to
                            --replicate_match, and number 3 condition assign 3 to
                            --replcate_match. For assigning the same
                            --replicate_match to all conditions, just assign like
                            all_1 (all condition use 1 --replicate_match). Default
                            is all_1.
      --utr_length UTR_LENGTH, -u UTR_LENGTH
                            The length of UTR. It is for Primary and Secondary
                            definition. Default is 300.
      --lib LIB, -l LIB     The libraries of TEX+/- wig files for TSSpredator. The
                            format is: wig_file_name:tex_treat_or_not(tex or notex
                            ):condition_id(integer):replicate_id(alphabet):strand(
                            + or -). For multiple wig files, please use comma to
                            separate the wig files. For example,
                            wig1:tex:1:a:+,wig2:tex:1:a:-.
      --output_prefix OUTPUT_PREFIX, -p OUTPUT_PREFIX
                            The output prefix of all conditions. For multiple
                            conditions, please use comma to separate them. For
                            example, prefix_condition1,prefix_condition2.
      --merge_manual MERGE_MANUAL, -m MERGE_MANUAL
                            If the gff file of manual checked TSS is provided, it
                            will merge manual checked ones and predicted ones.
                            please assign the path of gff file of manual checked
                            TSS.
      --statistics, -s      Doing statistics for TSS candidates. it will store in
                            statistics folder. Default is False.
      --validate_gene, -v   Using TSS candidates to validate genes in annotation
                            file. it will store in statistics folder. Default is
                            False.
      --compute_program COMPUTE_PROGRAM, -t COMPUTE_PROGRAM
                            The program for prediction (TSS or processing_site).
                            Default is TSS.
      --compare_transcript_assembly COMPARE_TRANSCRIPT_ASSEMBLY, -ta COMPARE_TRANSCRIPT_ASSEMBLY
                            For comparing with transcriptome assembly, please
                            assign the folder of gff file of transcript assembly.
                            Default is False.
      --fuzzy FUZZY, -fu FUZZY
                            The fuzzy for comparing TSS and transcript assembly.
                            Default is 5.
      --cluster CLUSTER, -c CLUSTER
                            This number is for comparing manual detected TSS and
                            prediced one. If the position between manual checked
                            one and predicted one is smaller or equal than this
                            value, It will only print one of them. Default is 2.
      --length LENGTH, -le LENGTH
                            The length of genome for comparing between predicted
                            one and manual checked one for statistics. For
                            comparing whole genome, please don't use this function
                            (Default). The default is comparing whole genome.
      --re_check_orphan, -ro
                            If the annotation file lacks information of gene or
                            locus_tag, all TSSs will be assigned to orphan TSSs.
                            The function can compare with CDS to classify the TSS.
                            Default is False.
      --overlap_feature OVERLAP_FEATURE, -of OVERLAP_FEATURE
                            If processing site and TSS are overlap, you can keep
                            "TSS" or "processing_site" or "both". Default is both.
      --reference_gff_folder REFERENCE_GFF_FOLDER, -rg REFERENCE_GFF_FOLDER
                            If --overlap_feature is "TSS" or "processing_site",
                            --reference_gff_folder need to be assigned. For
                            running TSS, please assign the folder of processing
                            site. For running processing_site, please assign the
                            folder of TSS. If --overlap_feature is "both", please
                            don't use this function (Default). Default is None
                            (for keep both).
      --remove_low_expression REMOVE_LOW_EXPRESSION, -rl REMOVE_LOW_EXPRESSION
                            If removing low expressed TSS/processing site is
                            needed, please assign the file of manual checked gff
                            file here. It will remove the low expressed ones based
                            on comparison of manual checked ones. Please Be
                            ATTENTION: this parameter may remove some True
                            positive, too. So, please make sure you want to do it.

- Output files

The output files will be stored in ``$ANNOgesic/output/TSS``.

``MasterTables``: The MasterTable from `TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_.

``statistics``: Statistics files.

The output files of processing sites are similar. Just replace ``TSS`` to ``processing_site``
like ``$ANNOgesic/output/processing_site``.

``configs``: The configuration files for running TSSpredator.

``gffs``: The gff files of TSSs.

There are some useful tags in the attributes of gff files:

``method``: The TSSs are from manual detection or `TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_.

``type``: The type of TSSs. It could be Primary, Secondary, Internal, Antisense or Orphan.

``utr_length``: The length of UTR.

``associated_gene``: Which genes are associated with this TSS.

``Parent``: Which transcript are associated with this TSS, if user has compared with transcript.

If the comparison between TSSs and genome annotation files is done, the tag - ``start_TSS`` will appear in the gff files 
of genome annotation. It represents the TSSs which associates with the CDS/tRNA/rRNA.

If the comparison between TSSs and transcripts is done, the tag - ``associated_tss`` will appear in the gff files
of transcript. It will show the associated TSSs which is in the transcript.

transcript_assembly
-------------------

``transcript_assembly`` will detect transcripts based on the coverage.

For importing the information of libraries, please refer to the section of 
``The format of libraries for import to ANNOgesic``.

- Pre-required tools and files

Wiggle files of fragmented libraries or TEX+/- treated libraries.

If user wants to compare transcripts with TSSs, it requires ``.gff`` files of TSSs.
If user wants to compare transcripts with genome anntation, it requires ``.gff`` files of genomes.

- Arguments

::

    usage: annogesic transcript_assembly [-h]
                                         [--annotation_folder ANNOTATION_FOLDER]
                                         [--length LENGTH]
                                         [--tex_wig_path TEX_WIG_PATH]
                                         [--frag_wig_path FRAG_WIG_PATH]
                                         [--height HEIGHT] [--width WIDTH]
                                         [--tolerance TOLERANCE]
                                         [--tolerance_coverage TOLERANCE_COVERAGE]
                                         [-rt REPLICATES_TEX]
                                         [--replicates_frag REPLICATES_FRAG]
                                         [--tex_notex TEX_NOTEX]
                                         [--compare_TSS COMPARE_TSS]
                                         [--compare_genome_annotation COMPARE_GENOME_ANNOTATION]
                                         [--compare_feature_genome COMPARE_FEATURE_GENOME]
                                         [--TSS_fuzzy TSS_FUZZY]
                                         [--Tex_treated_libs TEX_TREATED_LIBS]
                                         [--fragmented_libs FRAGMENTED_LIBS]
                                         [--table_best]
                                         [--terminator_folder TERMINATOR_FOLDER]
                                         [--fuzzy_term FUZZY_TERM]
                                         [--max_length_distribution MAX_LENGTH_DISTRIBUTION]
                                         [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --annotation_folder ANNOTATION_FOLDER, -g ANNOTATION_FOLDER
                            It is for comparing transcript assembly and genome
                            annotation gff file. This function can use annotation
                            gff file as reference and modify transcript assembly
                            file, If the genome annotation gff folder is provided.
      --length LENGTH, -l LENGTH
                            The minimum width of transcript. It is for comparing
                            to annotation file (--annotation_folder). If
                            --annotation_folder is assigned, it will be the final
                            output. Otherwise, --width would be minimum length for
                            the final output. The default is 20.
      --tex_wig_path TEX_WIG_PATH, -tw TEX_WIG_PATH
                            The path of TEX+/- wig folder.
      --frag_wig_path FRAG_WIG_PATH, -fw FRAG_WIG_PATH
                            The path of fragment wig folder.
      --height HEIGHT, -he HEIGHT
                            The minimum height of coverage to be a transcript. The
                            default is 10.
      --width WIDTH, -w WIDTH
                            The minimum width of transcript. It is for not
                            comparing to annotation file (--annotation_folder). It
                            will be the final output if --annotation_folder is not
                            provided. Otherwise, --length would be the minimum
                            length of transcript for the final output. The default
                            is 20.
      --tolerance TOLERANCE, -t TOLERANCE
                            It indicates the number of nucleotides which coverages
                            drop below --height can be ignore. The default is 5.
      --tolerance_coverage TOLERANCE_COVERAGE, -tc TOLERANCE_COVERAGE
                            If the coverage is lower than tolerance_coverage, even
                            the length is within --tolerance, it will still
                            terminate the current transcript. The default is 0.
      -rt REPLICATES_TEX, --replicates_tex REPLICATES_TEX
                            The transcript should be detected at least in the
                            number of the replicates. The format is
                            $NUMBERofCONDITION_$NUMBERofREPLICATED. For assigning
                            different --replicates_tex to different conditions,
                            please use comma to separate it. For example,
                            1_2,2_2,3_3 means number 1 and 2 condition assign 2 to
                            --replicate_tex, and number 3 condition assign 3 to
                            --replcate_tex. For assigning the same
                            --replicates_tex to all conditions, just assign like
                            all_1 (all condition use 1 --replicate_tex).
      --replicates_frag REPLICATES_FRAG, -rf REPLICATES_FRAG
                            The meaning and input type is the same as
                            --replicates_tex.
      --tex_notex TEX_NOTEX, -te TEX_NOTEX
                            For TEX+/- library, transcript should be detected in
                            both (TEX+ and TEX-) or can be detected in only one
                            library (TEX+ or TEX-). Please assign 1 or 2. Default
                            is 1.
      --compare_TSS COMPARE_TSS, -ct COMPARE_TSS
                            If comparing with TSS is needed, please assign TSS
                            folder.
      --compare_genome_annotation COMPARE_GENOME_ANNOTATION, -cg COMPARE_GENOME_ANNOTATION
                            If comparing with genome annotation file and searching
                            the parent transcript of gene is needed, please assign
                            annotation folder.
      --compare_feature_genome COMPARE_FEATURE_GENOME, -cf COMPARE_FEATURE_GENOME
                            If --compare_genome_annotation is provided, please
                            assign the feature which you want to compare. Default
                            is gene,CDS. For multiple features, just insert comma
                            between each feature, such as gene,CDS.
      --TSS_fuzzy TSS_FUZZY, -fu TSS_FUZZY
                            The fuzzy for comparing TSS and transcript assembly.
                            Default is 5.
      --Tex_treated_libs TEX_TREATED_LIBS, -tl TEX_TREATED_LIBS
                            Tex+/- library. The format is:
                            wig_file_name:TEX+/-(tex or notex):condition_id(intege
                            r):replicate_id(alphabet):strand(+ or -). For multiple
                            wig files, please use comma to separate the wig files.
                            For example, wig1:tex:1:a:+,wig2:tex:1:a:-.
      --fragmented_libs FRAGMENTED_LIBS, -fl FRAGMENTED_LIBS
                            Fragmented library. The format is: wig_file_name:fragm
                            ented(frag):condition_id(integer):replicate_id(alphabe
                            t):strand(+ or -). For multiple wig files, please use
                            comma to separate the wig files. For example,
                            wig1:frag:1:a:+,wig2:frag:1:a:-.
      --table_best, -tb     The output table only includes the best library.
                            Default is False.
      --terminator_folder TERMINATOR_FOLDER, -tr TERMINATOR_FOLDER
                            If comparing between transcripts and terminators is
                            needed, please assign the folder of gff files of
                            terminator here. Default is None.
      --fuzzy_term FUZZY_TERM, -fz FUZZY_TERM
                            If --terminator_folder is provided, please assign the
                            fuzzy here. Default is 30.
      --max_length_distribution MAX_LENGTH_DISTRIBUTION, -mb MAX_LENGTH_DISTRIBUTION
                            For generating the figure of distribution of
                            transcript length, please assign the maximum length
                            that you want to include. Default is 2000.

- Output files

The output files will be stored in ``$ANNOgesic/output/transcriptome_assembly``.

``gffs``: The gff files of transcript.

``tables``: The table of transcript with more details.

``statistics``: Statistics files.

There are some useful tags in gff files.

``compare_FEATURE``: The situation of overlap between transcripts and features (--compare_feature_genome)
(If --compare_genome_annotation is assigned.) 

``associated_tss``: Which TSSs are located in this transcripts. 
(If --compare_TSS is assigned.) 

``associated_$FEATURE``: It shows the feature (--compare_feature_genome) which are located in this transcripts.
(If --compare_genome_annotation is assigned.) 

``detect_lib``: The transcript is detected by tex-treated libraries or fragmented libraries.

``best_avg_coverage``: The average coverage of highest expressed library.

If --compare_genome_annotation is assigned, the associated transcript will be assigned in ``Parent``
of genome annotations.

If --compare_TSS is assigned, the tag - ``Parent`` will appear
in the gff files of TSSs. It will show which transcripts that TSSs are located.


terminator
-----------

``terminator`` will predict the rho-independent terminators. ``ANNOgesic`` combine the results of 
two methods in order to get more reliable candidates. First one is using `TranstermHP <http://transterm.cbcb.umd.edu/>`_.
The other one is detect the specific secondary structure between converging pairs  
of transcripts and CDSs. ``ANNOgesic`` can also compare with coverages in order to generate the terminators 
which has coverage significant decrease.

- Pre-required tools and files

`TranstermHP <http://transterm.cbcb.umd.edu/>`_

RNAfold of `ViennaRNA <http://www.tbi.univie.ac.at/RNA/>`_.

Gff files target genome annotation.

Fasta files of target genome sequence.

Wiggle files of TEX +/- treated libraries or fragmented libraries.

Gff files of transcript.

- Arguments

::

    usage: annogesic terminator [-h] [--TransTermHP_path TRANSTERMHP_PATH]
                                [--expterm_path EXPTERM_PATH]
                                [--RNAfold_path RNAFOLD_PATH]
                                [--fasta_folder FASTA_FOLDER]
                                [--annotation_folder ANNOTATION_FOLDER]
                                [--transcript_folder TRANSCRIPT_FOLDER]
                                [--sRNA SRNA] [--statistics]
                                [--tex_wig_folder TEX_WIG_FOLDER]
                                [--frag_wig_folder FRAG_WIG_FOLDER]
                                [--decrease DECREASE]
                                [--fuzzy_detect_coverage FUZZY_DETECT_COVERAGE]
                                [--fuzzy_within_transcript FUZZY_WITHIN_TRANSCRIPT]
                                [--fuzzy_downstream_transcript FUZZY_DOWNSTREAM_TRANSCRIPT]
                                [--fuzzy_within_gene FUZZY_WITHIN_GENE]
                                [--fuzzy_downstream_gene FUZZY_DOWNSTREAM_GENE]
                                [--highest_coverage HIGHEST_COVERAGE]
                                [-tl TEX_NOTEX_LIBS] [-fl FRAG_LIBS]
                                [-te TEX_NOTEX] [-rt REPLICATES_TEX]
                                [-rf REPLICATES_FRAG] [-tb] [-ml MIN_LOOP_LENGTH]
                                [-Ml MAX_LOOP_LENGTH] [-ms MIN_STEM_LENGTH]
                                [-Ms MAX_STEM_LENGTH] [-mr MISS_RATE]
                                [-mu MIN_U_TAIL_LENGTH] [-ru RANGE_U_TAIL] [-kp]
                                [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --TransTermHP_path TRANSTERMHP_PATH
                            Please assign the path of "transterm" in TransTermHP.
      --expterm_path EXPTERM_PATH
                            Please assign the path of expterm.dat for TransTermHP.
                            Default is /usr/local/bin/expterm.dat
      --RNAfold_path RNAFOLD_PATH
                            If you want to assign the path of "RNAfold" of Vienna
                            package, please assign here.
      --fasta_folder FASTA_FOLDER, -f FASTA_FOLDER
                            The path of genome fasta folder.
      --annotation_folder ANNOTATION_FOLDER, -g ANNOTATION_FOLDER
                            The path of genome annotation gff folder.
      --transcript_folder TRANSCRIPT_FOLDER, -a TRANSCRIPT_FOLDER
                            The path of transcript assembly gff folder.
      --sRNA SRNA, -sr SRNA
                            For including sRNA information, please assign the
                            folder of sRNA gff files.
      --statistics, -s      Doing statistics for terminator. The name of
                            statistics file is - stat_terminator_$STRAIN_NAME.csv.
                            Default is False.
      --tex_wig_folder TEX_WIG_FOLDER, -tw TEX_WIG_FOLDER
                            If TEX+/- libraries are provided, please assign TEX+/-
                            wig folder.
      --frag_wig_folder FRAG_WIG_FOLDER, -fw FRAG_WIG_FOLDER
                            If fragmented libraries are provided, please assign
                            fragmented wig folder.
      --decrease DECREASE, -d DECREASE
                            If the (lowest coverage / highest coverage) in the
                            terminator is smaller than this number, it will
                            consider this terminator have coverage dramatic
                            decreasing in it. Default is 0.5.
      --fuzzy_detect_coverage FUZZY_DETECT_COVERAGE, -fc FUZZY_DETECT_COVERAGE
                            For extending the region of detection of coverage
                            significant decreasing. Ex: the location of terminator
                            is 300-400, and --fuzzy_detect_coverage is 30. If the
                            coverage decrease is detected within 270-430, it will
                            still consider the terminator have coverage dramatic
                            decrease. Default is 30.
      --fuzzy_within_transcript FUZZY_WITHIN_TRANSCRIPT, -fut FUZZY_WITHIN_TRANSCRIPT
                            If the candidates are within transcript and the
                            distance between the end of gene/transcript and
                            terminator candidate is within this number, it will be
                            consider as terminator. Default is 30.
      --fuzzy_downstream_transcript FUZZY_DOWNSTREAM_TRANSCRIPT, -fdt FUZZY_DOWNSTREAM_TRANSCRIPT
                            The meaning is similar with --fuzzy_within_transcript.
                            It is for the candidates which are downstream of
                            transcript. Default is 30.
      --fuzzy_within_gene FUZZY_WITHIN_GENE, -fuc FUZZY_WITHIN_GENE
                            The meaning is similar with --fuzzy_within_transcript.
                            It is for gene in stead of transcript. Default is 10.
      --fuzzy_downstream_gene FUZZY_DOWNSTREAM_GENE, -fdg FUZZY_DOWNSTREAM_GENE
                            The meaning is similar with
                            --fuzzy_downstream_transcript. It is for gene in stead
                            of transcript. Default is 310.
      --highest_coverage HIGHEST_COVERAGE, -hc HIGHEST_COVERAGE
                            If the highest coverage of terminator is below to this
                            number, the terminator will be classify to non-detect
                            and not included in "best" results, but still included
                            in "all_candidates". Default is 10.
      -tl TEX_NOTEX_LIBS, --tex_notex_libs TEX_NOTEX_LIBS
                            Library name of TEX+/- library. The format is:
                            wig_file_name:TEX+/-(tex or notex):condition_id(intege
                            r):replicate_id(alphabet):strand(+ or -). For multiple
                            wig files, please use comma to separate the wig files.
                            For example, wig1:tex:1:a:+,wig2:tex:1:a:-.
      -fl FRAG_LIBS, --frag_libs FRAG_LIBS
                            Library name of fragmented library. The format is: wig
                            _file_name:fragmented(frag):condition_id(integer):repl
                            icate_id(alphabet):strand(+ or -). For multiple wig
                            files, please use comma to separate the wig files. For
                            example, wig1:frag:1:a:+,wig2:frag:1:a:-.
      -te TEX_NOTEX, --tex_notex TEX_NOTEX
                            For TEX+/- library, terminators should be detected in
                            both (TEX+ and TEX-) or can be detected in only one
                            library (TEX+ or TEX-). Please assign 1 or 2. Default
                            is 1.
      -rt REPLICATES_TEX, --replicates_tex REPLICATES_TEX
                            The Terminator candidates should be detected at least
                            in the number of the replicates. The format is
                            $NUMBERofCONDITION_$NUMBERofREPLICATED. For assigning
                            different --replicate_tex to different conditions,
                            please use comma to separate it. For example,
                            1_2,2_2,3_3 means number 1 and 2 condition assign 2 to
                            --replicate_tex, and number 3 condition assign 3 to
                            --replcate_tex. For assigning the same
                            --replicates_tex to all conditions, just assign like
                            all_1 (all condition use 1 --replicate_tex).
      -rf REPLICATES_FRAG, --replicates_frag REPLICATES_FRAG
                            The meaning and input type is the same as
                            --replicates_tex.
      -tb, --table_best     Output table only contains the library which has
                            coverage most significant decreasing. Default is
                            False.
      -ml MIN_LOOP_LENGTH, --min_loop_length MIN_LOOP_LENGTH
                            The minimum length of loop for terminator. Default is
                            3 nts.
      -Ml MAX_LOOP_LENGTH, --max_loop_length MAX_LOOP_LENGTH
                            The maximum length of loop for terminator. Default is
                            10 nts.
      -ms MIN_STEM_LENGTH, --min_stem_length MIN_STEM_LENGTH
                            The minimum length of stem for terminator. Default is
                            4 nts.
      -Ms MAX_STEM_LENGTH, --max_stem_length MAX_STEM_LENGTH
                            The maximum length of stem for terminator. Default is
                            20 nts.
      -mr MISS_RATE, --miss_rate MISS_RATE
                            The percentage of nucleotides which can be no base
                            pair in the stem. Default is 0.25.
      -mu MIN_U_TAIL_LENGTH, --min_U_tail_length MIN_U_TAIL_LENGTH
                            The minimum length of U tail for terminator. Default
                            is 3 nts.
      -ru RANGE_U_TAIL, --range_U_tail RANGE_U_TAIL
                            The range of nucleotides for detection of U tail. For
                            example, if --range_U_tail is 6 and
                            --min_U_tail_length is 3, and there are 3 U within 6
                            nts, it will be assigned to detecting U tail
                            successfully. Default is 6.
      -kp, --keep_multi_term
                            Sometimes, one gene is associated with more terminator
                            candidates. In default, it will only keep the high
                            confident one. The function can keep all terminators
                            which associated with the same gene. Default is False.

- Output files

The output files will be stored in ``$ANNOgesic/output/terminator``.

``statistics``: Statistics files.

``transtermhp``: All output of `TranstermHP <http://transterm.cbcb.umd.edu/>`_.

``gffs``: Gff files of terminator.
There are four different sub-folders to store terminators.

``all_candidates`` will store all terminators which ``ANNOgesic`` can detect.

``express`` will store the terminators which has gene expression.

``best`` will store the terminators which not only has gene expression but also
has coverage dramatic decrease.

``non_express`` will store the terminators which has no gene expression.

``tables``: The tables of terminators with more details.

The tags of gff files:

``method``: The method that this terminator be detected.

``coverage_decrease``: The coverage of the terminator has dramatic decreasing or not.

``express``: The terminator has gene expression or not.

``diff_coverage``: The highest coverage and lowest coverage of the library which expresses highest.
The numbers in parens are highest coverage and lowest coverage.

``Parent``: This tag presents the parent transcript of terminator.

utr
-----

``utr`` can compare with TSSs, CDSs/tRNAs/sRNAs, transcripts and terminators
to generate proper UTRs. 5'UTRs are based on detecting the regions between TSSs and CDSs/tRNAs/sRNAs. 
3'UTRs are based on detecting the 
regions between the end of transcripts and CDSs/tRNAs/sRNAs. If the gff files of TSSs are not computed by 
ANNOgesic, please use --TSS_source. ``utr`` would classify TSSs for the analysis.

- Pre-required files

Gff files of genome annotations, TSSs and transcripts.

If the information of terminators is needed, the gff files of terminators are required.

- Arguments

::

    usage: annogesic utr [-h] [--annotation_folder ANNOTATION_FOLDER]
                         [--TSS_folder TSS_FOLDER]
                         [--transcript_assembly_folder TRANSCRIPT_ASSEMBLY_FOLDER]
                         [--terminator_folder TERMINATOR_FOLDER] [--TSS_source]
                         [--base_5UTR BASE_5UTR] [--UTR_length UTR_LENGTH]
                         [--base_3UTR BASE_3UTR]
                         [--terminator_fuzzy TERMINATOR_FUZZY]
                         [--fuzzy_3utr FUZZY_3UTR] [--fuzzy_5utr FUZZY_5UTR]
                         [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --annotation_folder ANNOTATION_FOLDER, -g ANNOTATION_FOLDER
                            The path of genome annotation gff folder.
      --TSS_folder TSS_FOLDER, -t TSS_FOLDER
                            The path of TSS folder.
      --transcript_assembly_folder TRANSCRIPT_ASSEMBLY_FOLDER, -a TRANSCRIPT_ASSEMBLY_FOLDER
                            The path of transcriptome assembly folder.
      --terminator_folder TERMINATOR_FOLDER, -e TERMINATOR_FOLDER
                            If including the information of terminator is
                            required, please assign the path of terminator folder
                            here.
      --TSS_source, -s      If TSS file is not generated from ANNOgesic, please
                            turn it on. Default is True(ANNOgesic).
      --base_5UTR BASE_5UTR, -b5 BASE_5UTR
                            The information for detection of 5'UTR. It can be
                            "TSS" or "transcript" or "both". Default is both.
      --UTR_length UTR_LENGTH, -l UTR_LENGTH
                            The maximum length of UTR. Default is 300.
      --base_3UTR BASE_3UTR, -b3 BASE_3UTR
                            the information for detection of 3'UTR. It can be
                            "transcript" or "terminator" or "both". Default is
                            transcript.
      --terminator_fuzzy TERMINATOR_FUZZY, -f TERMINATOR_FUZZY
                            This is only for --base_3UTR which assigned by
                            "transcript" or "both", and terminator file are
                            provided. If the distance (nucleotides) between
                            terminator and the end of transcript lower than this
                            value, it will assign the terminator associated with
                            the 3'UTR. Therefore, the relationship between 3'UTRs
                            and terminator can be gain. Default is 30.
      --fuzzy_3utr FUZZY_3UTR, -f3 FUZZY_3UTR
                            If --base_3UTR includes transcript, please assign the
                            fuzzy of 3'UTR. Default is 10 nucleotides.
      --fuzzy_5utr FUZZY_5UTR, -f5 FUZZY_5UTR
                            If --base_5UTR includes transcript, please assign the
                            fuzzy of 5'UTR. Default is 5 nucleotides.

- Output files

All output of 5'UTRs will be stored in ``$ANNOgesic/output/UTR/5UTR``.

All output of 3'UTRs will be stored in ``$ANNOgesic/output/UTR/3UTR``.

``gffs``: Gff files of 5'UTR/3'UTR

The tags of gff files:

``length``: UTR length.

``associated_cds``: Which CDSs/rRNAs/tRNAs are associated with this UTR.

``associated_gene``: Which genes are associated with this UTR.

``Parent``: Which transcript is associated with this UTR.

``associated_tss``: Which TSSs are associated with this 5'UTR.

``tss_type``: What types of TSSs are associated with this 5'UTR.

``associated_term``: Which terminators are associated with this 3'UTR.

srna
-----
``srna`` can predict different types of sRNAs. For intergenic and antisense sRNA, it 
is detected via comparison of the transcripts and annotation profile. 
For UTR-derived sRNA, the detection is based on the TSSs and processing sites, 
transcript and genome annotation.

- Pre-required tools and files

Gff files of genome annotation and Transcript data.

wiggle files: The libraries and wiggle files, Please refer to the ``The format of libraries for import to ANNOgesic``.

If you want to detect the UTR-derived sRNAs, it is necessary to input
TSS information. It is for the detection of 5'UTR-derived sRNA and interCDS-derived sRNA. 
If you don't want to detect UTR-derived sRNAs,
TSS information still can be provided as a filter.

Optional input file:

processing site: It is for checking the sRNAs which end with processing sites. Moreover,
Some 3'UTR-derived and interCDS-derived sRNA candidates start
from processing sites not TSSs. If you don't want to detect UTR-derived sRNAs,
This information still can be provided to increase the accuracy, especially for some
long non-coding RNAs.

There are some filters can improve the prediction, and they need some input information.

Secondary structure: `ViennaRNA <http://www.tbi.univie.ac.at/RNA/>`_, 
`Ps2pdf14 <http://pages.cs.wisc.edu/~ghost/doc/AFPL/6.50/Ps2pdf.htm>`_ and 
Fasta files of genome sequence is required.

TSS: gff file of TSS is necessary. Based on different types of TSS, the cutoff of coverage can be 
assigned separately.

Searching sRNA database: `Blast+ <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`_ and 
sRNA database, such as `BSRD <http://www.bac-srna.org/BSRD/index.jsp>`_ are necessary.
The format of header should be ``$ID|$STRAIN|$SRNANAME``.
For example, ``>saci2813.1|Acinetobacter sp. ADP1|Aar``. 
The ID is saci403.1; the strain of this sRNA is Acinetobacter sp. ADP1 and the name of sRNA is Aar.
If the database does not follow the format, it will occur error when the user runs with ``--sRNA_blast_stat, -sb``.
Or the results will be meaningless. User can also assign ``--best_with_all_sRNAhit, -ba`` for including
all candidates which can find hits in sRNA database (even without matching other filters).

Searching nr database: `Blast+ <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`_ and
`nr database <ftp://ftp.ncbi.nih.gov/blast/db/FASTA/>`_ are required.

Terminator: gff file of terminator is needed. If sRNA must be associated with terminator, 
please add ``-bt``. Without ``-bt``, it will only provide to the relationship between terminator and sRNA 
but does not remove the the candidates which not associated with terminator.

sORF: gff file of sORF needs to be assigned. If sRNA must be not overlap with sORF,  
please add ``-bs``. Without ``-bs``, it will only provide to the relationship between sORF and sRNA
but does not remove the the candidates which are overlap with sORF.

Promoter: Table of promoters have to be input. It can compare between promoters and sRNA for ranking 
the sRNA candidate as well. Moreover, If sRNA must be associated with promoter,
please add ``-bp``. Without ``-bp``, it will only provide to the relationship between promoter and sRNA
but does not remove the the candidates which not associated with promoter.


The format should be

===========  ============  ==========  =======
strain       TSS_position  TSS_strand  Motif
-----------  ------------  ----------  -------
NC_000915.1  237118        \-          MOTIF_1
NC_000915.1  729009        \-          MOTIF_1
===========  ============  ==========  =======

First row is the header of table, the last column is the name of motif/promoter.
If subcommand ``promoter`` was used for detecting promoter, the table will be generated automatically.
Please refer to the section of ``promoter``.


- Arguments

::

    usage: annogesic srna [-h] [--Vienna_folder VIENNA_FOLDER]
                          [--Vienna_utils VIENNA_UTILS]
                          [--blast_plus_folder BLAST_PLUS_FOLDER]
                          [--ps2pdf14_path PS2PDF14_PATH] [--UTR_derived_sRNA]
                          [--import_info IMPORT_INFO]
                          [--transcript_assembly_folder TRANSCRIPT_ASSEMBLY_FOLDER]
                          [--annotation_folder ANNOTATION_FOLDER]
                          [--TSS_folder TSS_FOLDER]
                          [--processing_site_folder PROCESSING_SITE_FOLDER]
                          [--promoter_table PROMOTER_TABLE]
                          [--promoter_name PROMOTER_NAME] [--TSS_source]
                          [--TSS_intergenic_fuzzy TSS_INTERGENIC_FUZZY]
                          [--TSS_5UTR_fuzzy TSS_5UTR_FUZZY]
                          [--TSS_3UTR_fuzzy TSS_3UTR_FUZZY]
                          [--TSS_interCDS_fuzzy TSS_INTERCDS_FUZZY]
                          [--terminator_folder TERMINATOR_FOLDER]
                          [--terminator_fuzzy_in_CDS TERMINATOR_FUZZY_IN_CDS]
                          [--terminator_fuzzy_out_CDS TERMINATOR_FUZZY_OUT_CDS]
                          [--min_length MIN_LENGTH] [--max_length MAX_LENGTH]
                          [--tex_wig_folder TEX_WIG_FOLDER]
                          [--frag_wig_folder FRAG_WIG_FOLDER]
                          [--run_intergenic_TEX_coverage RUN_INTERGENIC_TEX_COVERAGE]
                          [--run_intergenic_noTEX_coverage RUN_INTERGENIC_NOTEX_COVERAGE]
                          [--run_intergenic_fragmented_coverage RUN_INTERGENIC_FRAGMENTED_COVERAGE]
                          [--run_antisense_TEX_coverage RUN_ANTISENSE_TEX_COVERAGE]
                          [--run_antisense_noTEX_coverage RUN_ANTISENSE_NOTEX_COVERAGE]
                          [--run_antisense_fragmented_coverage RUN_ANTISENSE_FRAGMENTED_COVERAGE]
                          [--intergenic_tolerance INTERGENIC_TOLERANCE]
                          [--run_utr_TEX_coverage RUN_UTR_TEX_COVERAGE]
                          [--run_utr_noTEX_coverage RUN_UTR_NOTEX_COVERAGE]
                          [--run_utr_fragmented_coverage RUN_UTR_FRAGMENTED_COVERAGE]
                          [--min_utr_coverage MIN_UTR_COVERAGE]
                          [--fasta_folder FASTA_FOLDER]
                          [--cutoff_energy CUTOFF_ENERGY] [--mountain_plot]
                          [--nr_format] [--srna_format]
                          [--sRNA_database_path SRNA_DATABASE_PATH]
                          [--nr_database_path NR_DATABASE_PATH]
                          [--tex_notex_libs TEX_NOTEX_LIBS]
                          [--frag_libs FRAG_LIBS] [--tex_notex TEX_NOTEX]
                          [-rt REPLICATES_TEX] [--replicates_frag REPLICATES_FRAG]
                          [--table_best]
                          [--decrease_intergenic DECREASE_INTERGENIC]
                          [--decrease_utr DECREASE_UTR]
                          [--fuzzy_intergenic FUZZY_INTERGENIC]
                          [--fuzzy_utr FUZZY_UTR] [--cutoff_nr_hit CUTOFF_NR_HIT]
                          [--blast_e_nr BLAST_E_NR] [--blast_e_srna BLAST_E_SRNA]
                          [--sORF SORF] [--best_with_all_sRNAhit]
                          [--best_without_sORF_candidate] [--best_with_terminator]
                          [--best_with_promoter] [--detect_sRNA_in_CDS]
                          [--overlap_percent_CDS OVERLAP_PERCENT_CDS]
                          [--ignore_hypothetical_protein]
                          [--ranking_time_promoter RANKING_TIME_PROMOTER]
                          [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --Vienna_folder VIENNA_FOLDER
                            Please assign the folder of Vienna package. It should
                            include RNAfold.
      --Vienna_utils VIENNA_UTILS
                            Please assign the folder of Utils of Vienna package.
                            It should include relplot.pl and mountain.pl.
      --blast_plus_folder BLAST_PLUS_FOLDER
                            Please assign the folder of blast+ which include
                            blastn, blastx, makeblastdb.
      --ps2pdf14_path PS2PDF14_PATH
                            Please assign the path of ps2pdf14.
      --UTR_derived_sRNA, -u
                            The function is for detecting UTR-derived sRNA.
                            Default is False.
      --import_info IMPORT_INFO, -d IMPORT_INFO
                            There are several types of information that you can
                            import to detect and filter sRNAs: tss (the sRNA
                            should start from a TSS), sec_str (free energy change
                            of secondary structure (normalized by length)),
                            blast_nr (blast to non-redundant database), blast_srna
                            (blast to sRNA), sorf (compare with sORF), term
                            (compare with terminator), promoter (compare with
                            promoter motif). ATTENTION: without filters, the
                            results may include many false positives. Please
                            assign the information that you want to import (comma
                            for separate the filters), ex: tss,sec_str,blast_nr -
                            means it used 1. TSS, 2. free energy change of
                            secondary structure and 3. blast to nr database to
                            detect sRNA. Besides these, it will also consider the
                            sequence length of sRNA. ATTENTION: if you want to
                            import sRNA database, please follow the format:
                            $ID|$STRAIN|$SRNANAME. Default is
                            tss,sec_str,blast_nr,blast_srna.
      --transcript_assembly_folder TRANSCRIPT_ASSEMBLY_FOLDER, -a TRANSCRIPT_ASSEMBLY_FOLDER
                            The path of transcriptome assembly folder.
      --annotation_folder ANNOTATION_FOLDER, -g ANNOTATION_FOLDER
                            The path of genome annotation gff folder.
      --TSS_folder TSS_FOLDER, -t TSS_FOLDER
                            If TSS information is needed, please assign the path
                            of gff folder of TSS. For detection of UTR-derived
                            sRNA, TSS information MUST be provided.
      --processing_site_folder PROCESSING_SITE_FOLDER, -p PROCESSING_SITE_FOLDER
                            If processing site information is needed, please
                            assign the path of gff folder of processing site.For
                            detection of UTR-derived sRNA, processing site
                            information can improve the results.
      --promoter_table PROMOTER_TABLE, -pt PROMOTER_TABLE
                            If promoter information is needed, please assign the
                            path of promoter table. The format of table is $STRAIN
                            $TSS_POSITION $TSS_STRAND $PROMOTER_NAME. TSS
                            information is also required.
      --promoter_name PROMOTER_NAME, -pn PROMOTER_NAME
                            If --promoter_table is provided, please assign the
                            promoter name (the last column of promoter table)
                            which you want to compare. For multiple promoters,
                            please put comma between the promoters. Default is
                            None.
      --TSS_source, -ts     If the gff file of TSS is not generated by ANNOgesic,
                            please use this function to classify TSSs and generate
                            the proper format for sRNA prediction. Default is
                            True.
      --TSS_intergenic_fuzzy TSS_INTERGENIC_FUZZY, -ft TSS_INTERGENIC_FUZZY
                            If --TSS_folder is provided, please assign the fuzzy
                            for comparing TSS and transcript. It is for intergenic
                            sRNA. Default is 3.
      --TSS_5UTR_fuzzy TSS_5UTR_FUZZY, -f5 TSS_5UTR_FUZZY
                            If --TSS_folder is provided, please assign the fuzzy
                            for comparing TSS and transcript. It is for 5'UTR of
                            UTR derived sRNA.The type of input can be percentage
                            or the real amount of reads. Ex: p_0.05 means the
                            fuzzy is 5 percent of the length of 5'UTR. n_10 means
                            the fuzzy is 10 base pair. Default is n_3.
      --TSS_3UTR_fuzzy TSS_3UTR_FUZZY, -f3 TSS_3UTR_FUZZY
                            The meaning is similar with --TSS_5UTR_fuzzy. It is
                            for 3'UTR instead of 5'UTR. Default is p_0.04.
      --TSS_interCDS_fuzzy TSS_INTERCDS_FUZZY, -fc TSS_INTERCDS_FUZZY
                            The meaning is similar with --TSS_5UTR_fuzzy. It is
                            for interCDS instead of 5'UTR. Default is p_0.04.
      --terminator_folder TERMINATOR_FOLDER, -tf TERMINATOR_FOLDER
                            If terminator information is needed, please assign the
                            path of gff folder of terminator.
      --terminator_fuzzy_in_CDS TERMINATOR_FUZZY_IN_CDS, -tfi TERMINATOR_FUZZY_IN_CDS
                            If --terminator_folder is provided, please assign the
                            fuzzy for comparing terminator and transcript. It is
                            the fuzzy for the terminator which is within CDS.
                            Default is 30.
      --terminator_fuzzy_out_CDS TERMINATOR_FUZZY_OUT_CDS, -tfo TERMINATOR_FUZZY_OUT_CDS
                            If --terminator_folder is provided, please assign the
                            fuzzy for comparing terminator and transcript. It is
                            the fuzzy for the terminator which is outside of CDS.
                            Default is 30.
      --min_length MIN_LENGTH, -lm MIN_LENGTH
                            Please assign the minimum length of sRNA. Default is
                            30.
      --max_length MAX_LENGTH, -lM MAX_LENGTH
                            Please assign the maximum length of sRNA. Default is
                            500.
      --tex_wig_folder TEX_WIG_FOLDER, -tw TEX_WIG_FOLDER
                            The path of TEX+/- wig folder.
      --frag_wig_folder FRAG_WIG_FOLDER, -fw FRAG_WIG_FOLDER
                            The path of fragment wig folder.
      --run_intergenic_TEX_coverage RUN_INTERGENIC_TEX_COVERAGE, -it RUN_INTERGENIC_TEX_COVERAGE
                            The minimum average coverage of intergenic sRNA
                            candidates for TEX+. The cutoff of coverage for sRNA
                            prediction is based on different types of TSSs. The
                            order of numbers is "Primary Secondary Internal
                            Antisense Orphan" (separated by comma). Ex: The cutoff
                            is 0,0,0,50,10, it means that antisense (cutoff
                            coverage is 50) TSS and orphan (cutoff coverage is 10)
                            TSS are used for sRNA prediction. 0 means that not use
                            it for prediction. If TSS information is not provided,
                            it will choose the lowest one as a general cutoff for
                            prediction. Ex: if the cutoff is 0,0,0,50,10 and
                            --TSS_folder is not provided, 10 will be the general
                            cutoff for prediction. Default is 0,0,0,40,20.
      --run_intergenic_noTEX_coverage RUN_INTERGENIC_NOTEX_COVERAGE, -in RUN_INTERGENIC_NOTEX_COVERAGE
                            The meaning is the same as
                            --run_intergenic_TEX_coverage. This is for TEX-
                            library. Default is 0,0,0,30,10.
      --run_intergenic_fragmented_coverage RUN_INTERGENIC_FRAGMENTED_COVERAGE, -if RUN_INTERGENIC_FRAGMENTED_COVERAGE
                            The meaning is the same as
                            --run_intergenic_TEX_coverage. This is for fragmented
                            library. Default is 400,200,0,50,20.
      --run_antisense_TEX_coverage RUN_ANTISENSE_TEX_COVERAGE, -at RUN_ANTISENSE_TEX_COVERAGE
                            The meaning is the same as
                            --run_intergenic_TEX_coverage. Just apply to
                            antisense. Default is 0,0,0,40,20.
      --run_antisense_noTEX_coverage RUN_ANTISENSE_NOTEX_COVERAGE, -an RUN_ANTISENSE_NOTEX_COVERAGE
                            The meaning is the same as
                            --run_intergenic_noTEX_coverage. Just apply to
                            antisense. Default is 0,0,0,30,10.
      --run_antisense_fragmented_coverage RUN_ANTISENSE_FRAGMENTED_COVERAGE, -af RUN_ANTISENSE_FRAGMENTED_COVERAGE
                            The meaning is the same as
                            --run_intergenic_fragmented_coverage. Just apply to
                            antisense. Default is 400,200,0,50,20.
      --intergenic_tolerance INTERGENIC_TOLERANCE, -ti INTERGENIC_TOLERANCE
                            This number indicates the tolerance of temporary drop
                            below cutoff of coverage. Default is 5.
      --run_utr_TEX_coverage RUN_UTR_TEX_COVERAGE, -ut RUN_UTR_TEX_COVERAGE
                            The minimum average coverage of UTR-derived sRNA
                            candidates for TEX+. The cutoff can be assigned by the
                            percentile or real number of coverage. The order of
                            numbers are "5'UTR, 3'UTR and interCDS" (separated by
                            comma). Ex: if the cutoff is "p_0.7,p_0.5,p_0.5", it
                            will use 70 percentile of coverage as cutoff for
                            5'UTR, median of coverage as cutoff for 3'UTR and
                            interCDS. Ex: if the cutoff is "n_30,n_10,n_20 " it
                            will use 30 as cutoff for 5'UTR and 10 as cutoff for
                            3'UTR and 20 for interCDS. Default is
                            p_0.8,p_0.6,p_0.7.
      --run_utr_noTEX_coverage RUN_UTR_NOTEX_COVERAGE, -un RUN_UTR_NOTEX_COVERAGE
                            The meaning is the same as --run_utr_TEX_coverage.
                            This is for TEX- library. Default is
                            p_0.7,p_0.5,p_0.6.
      --run_utr_fragmented_coverage RUN_UTR_FRAGMENTED_COVERAGE, -uf RUN_UTR_FRAGMENTED_COVERAGE
                            The meaning is the same as --run_utr_TEX_coverage.
                            This is for fragmented library. Default is
                            p_0.7,p_0.5,p_0.6.
      --min_utr_coverage MIN_UTR_COVERAGE, -mu MIN_UTR_COVERAGE
                            The minimum coverage of UTR-derived sRNA. The coverage
                            should not only fit the --run_utr_TEX_coverage,
                            --run_utr_noTEX_coverage and
                            --run_utr_fragmented_coverage, but also this value.
                            Defaul is 50.
      --fasta_folder FASTA_FOLDER, -f FASTA_FOLDER
                            If "sec_str" or "blast_nr" or "blast_srna" is assigned
                            to --import_info, please assign the path of genome
                            fasta folder.
      --cutoff_energy CUTOFF_ENERGY, -e CUTOFF_ENERGY
                            If secondary structure information is needed, please
                            assign the cutoff of folding energy change (normalized
                            by length of gene). Default is -0.05.
      --mountain_plot, -m   It is for generating mountain plot of sRNA candidate.
                            Default is False.
      --nr_format, -nf      If nr database is not formatted, it is for formating
                            nr database. Default is False.
      --srna_format, -sf    If sRNA databse is not formatted, it is for formating
                            sRNA database. Default is False.
      --sRNA_database_path SRNA_DATABASE_PATH, -sd SRNA_DATABASE_PATH
                            If blast results of sRNA is needed, please assign the
                            path of sRNA database.
      --nr_database_path NR_DATABASE_PATH, -nd NR_DATABASE_PATH
                            If blast results of nr is needed, please assign the
                            path of nr database.
      --tex_notex_libs TEX_NOTEX_LIBS, -tl TEX_NOTEX_LIBS
                            library name of TEX+/- libraries. The format is:
                            wig_file_name:TEX+/-(tex or notex):condition_id(intege
                            r):replicate_id(alphabet):strand(+ or -). For multiple
                            wig files, please use comma to separate the wig files.
                            For example, wig1:tex:1:a:+,wig2:tex:1:a:-.
      --frag_libs FRAG_LIBS, -fl FRAG_LIBS
                            library name of fragmented libraries. The format is: w
                            ig_file_name:fragmented(frag):condition_id(integer):re
                            plicate_id(alphabet):strand(+ or -). For multiple wig
                            files, please use comma to separate the wig files. For
                            example, wig1:frag:1:a:+,wig2:frag:1:a:-.
      --tex_notex TEX_NOTEX, -te TEX_NOTEX
                            For TEX+/- library, sRNA should be detected in both
                            (TEX+ and TEX-) or can be detected in only one library
                            (TEX+ or TEX-). Please assign 1 or 2. Default is 2.
      -rt REPLICATES_TEX, --replicates_tex REPLICATES_TEX
                            The sRNA should be detected at least in the number of
                            the replicates. The format is
                            $NUMBERofCONDITION_$NUMBERofREPLICATED. For assigning
                            different --replicates_tex to different conditions,
                            please use comma to separate it. For example,
                            1_2,2_2,3_3 means number 1 and 2 condition assign 2 to
                            --replicates_tex, and number 3 condition assign 3 to
                            --replcates_tex. For assigning the same
                            --replicates_tex to all conditions, just assign like
                            all_1 (all condition use 1 --replicates_tex).
      --replicates_frag REPLICATES_FRAG, -rf REPLICATES_FRAG
                            The meaning and input type is the same as
                            --replicates_tex.
      --table_best, -tb     The output table of sRNA candidates only prints the
                            best library. Default is False.
      --decrease_intergenic DECREASE_INTERGENIC, -di DECREASE_INTERGENIC
                            If the length of intergenic potential sRNA is longer
                            than the max_length, it will check the sRNA candidates
                            based on coverage. If (the lowest coverage / the
                            highest coverage) is smaller than this number, it will
                            consider that the spot of lowest coverage as end
                            point. If new length is suitable for a sRNA candidate,
                            this candiate is included as a sRNA. Default is 0.1.
      --decrease_utr DECREASE_UTR, -du DECREASE_UTR
                            It is similar with --decrease_intergenic. This is for
                            UTR-derived sRNAs. Default is 0.05.
      --fuzzy_intergenic FUZZY_INTERGENIC, -fi FUZZY_INTERGENIC
                            If the situation is like --decrease_intergenic
                            mentioned, This is a fuzzy value between the end of
                            sRNA. Default is 10.
      --fuzzy_utr FUZZY_UTR, -fu FUZZY_UTR
                            It is simliar with --fuzzy_intergenic. This is for
                            UTR-derived sRNAs. Default is 10.
      --cutoff_nr_hit CUTOFF_NR_HIT, -cn CUTOFF_NR_HIT
                            The cutoff of hits number in nr database. If the
                            number of nr hits more than this cutoff, it will be
                            excluded. Default is 0.
      --blast_e_nr BLAST_E_NR, -en BLAST_E_NR
                            The cutoff of blast e value for nr alignment. Default
                            is 0.0001.
      --blast_e_srna BLAST_E_SRNA, -es BLAST_E_SRNA
                            The cutoff of blast e value for sRNA alignment.
                            Default is 0.0001.
      --sORF SORF, -O SORF  If comparing sORF and sRNA is needed, please assign
                            the path of sORF gff folder.
      --best_with_all_sRNAhit, -ba
                            The sRNA candidates which have the homology in sRNA
                            database can be included in best results without
                            fitting other information (ex. TSS, blast in nr...) if
                            this parameter is True. Or it will just select the
                            best candidates based on all filter conditions.
                            Default is False.
      --best_without_sORF_candidate, -bs
                            It is for generating the best sRNA candidates without
                            including the sRNA candidates which are overlap with
                            sORFs.Default is False.
      --best_with_terminator, -bt
                            It is for generating the best sRNA candidates which
                            must be associated with terminator. It also includes
                            the sRNA which ends with processing site. Default is
                            False.
      --best_with_promoter, -bp
                            It is for generating the best sRNA candidates which is
                            must be associated with promoter.Default is False.
      --detect_sRNA_in_CDS, -ds
                            It is for searching sRNA in CDS (ex: the genome
                            annotation is not correct). It may find more sRNA
                            candidates which overlap with CDS. Default is False.
      --overlap_percent_CDS OVERLAP_PERCENT_CDS, -oc OVERLAP_PERCENT_CDS
                            If --detect_sRNA_in_CDS is True, please assign the
                            cutoff of the ratio of overlap between CDS and sRNA
                            candidates. Default is 0.5
      --ignore_hypothetical_protein, -ih
                            For ignoring hypothetical protein in genome annotation
                            file. Default is False.
      --ranking_time_promoter RANKING_TIME_PROMOTER, -rp RANKING_TIME_PROMOTER
                            If --promoter_table is provided, it will also use for
                            ranking sRNA candidates. The ranking will base on
                            --ranking_time_promoter * average coverage. For
                            example, one candidates which average coverage is 10,
                            associated with promoter and --ranking_time_promoter
                            is 2, the score for ranking will be 20 (2*10). The
                            candidates which are not associated with promoters,
                            the --ranking_time_promoter is 1. Default is 2. This
                            number can not be smaller than 1.

- Output files

All output files will be stored in ``$ANNOgesic/output/sRNA``.

``sRNA_2d_$STRAIN_NAME``: The secondary structure of all sRNA candidates.

``sRNA_seq_$STRAIN_NAME``: The sequence of all sRNA candidates.

``blast_result_and_misc``: The results of blast.

``mountain_plot``: The mountain plots of sRNA candidates.

``sec_structure``: The dot plots and secondary structure plots of sRNA candidates.

``statistics``: Statistics files. ``stat_$STRAIN_NAME_sRNA_blast.csv`` is the results of analysis of blast sRNA databases.
``stat_sRNA_class_Staphylococcus_aureus_HG003.csv`` is the results of classification of sRNA candidates.

``tables``: sRNA tables with more details. It also includes the ranking of sRNA candidates. 
``for class`` is for different classes of sRNAs.
``best`` is the best results of sRNAs after filtering. ``all_candidates`` is for all candidates without filtering.

``gffs``: Gff files of sRNAs. The meanings of ``for class``, ``best``, ``all_candidates`` are the same as ``tables``.

The useful tags of gff files:

``sRNA_type``: The sRNA is from 5'UTR, 3'UTR, interCDS, intergenic, antisense or within CDS.

``with_TSS``: Which TSSs are related to this sRNA. "NA" means the sRNA is not related to any TSSs.

``sORF``: Which sORFs are overlap with this sRNA.

``sRNA_hit``: The blast hit of sRNA database.

``nr_hit``: The blast hit of nr database.

``2d_energy``: The normalized (by the length of sRNA) free energy change of secondary structure of sRNA candidate.

``with_term``: The terminators which are associated with the sRNA candidate.

If you assigned ``--TSS_source`` for sRNA prediction, ``TSS_class`` will be generated and store the gff files of TSSs.

``promoter``: The promoters which are associated with this sRNA candidate.

``overlap_cds``: The CDSs which are overlap with sRNA.

``overlap_percent``: If there are CDSs overlap with sRNA, it will shows the percentage of overlap.

sorf
----------
``sorf`` can detect sORF based on searching ribosome binding sites, start codons and stop codons within the non-annotated transcripts.
Since non-annotated region may be sRNAs or sORFs, it also provides the function to compare sORFs and sRNAs. 
If there are some sORFs overlaped, it will merge them together. Therefore, one region may contain more than one sORF. 
BE CAREFUL, The position of start codon is assigned to the first nucleotide. The position of stop codon is the last nucleotide. 
``sorf`` provides the region which covers all possible sORFs. Thus, the region may contain different frame shift. 
Ex: (200, 202, 203) are the positions of three start codons and (241, 243) are two stop codons in 
a small transcript. Therefore, there are three possible ORFs(200-241, 203-241 and 202-243).
Please be aware this point for using the results.

- Pre-required tools and files

The gff files of genome annotation and transcripts.

The libraries and wiggle files, Please refer to the ``The format of libraries for import to ANNOgesic``.

The fasta files of genome sequence for detection of ribosome binding sites, start codons and stop codons.

Some useful information can be used to improve the prediction:

gff files of TSSs for checking the sORFs start from TSS or not. 

gff files of sRNAs for checking the overlap of sRNAs and sORFs.

- Arguments

::

    usage: annogesic sorf [-h] [--UTR_derived_sORF]
                          [--transcript_assembly_folder TRANSCRIPT_ASSEMBLY_FOLDER]
                          [--annotation_folder ANNOTATION_FOLDER]
                          [--TSS_folder TSS_FOLDER] [--utr_length UTR_LENGTH]
                          [--min_length MIN_LENGTH] [--max_length MAX_LENGTH]
                          [--tex_wig_folder TEX_WIG_FOLDER]
                          [--frag_wig_folder FRAG_WIG_FOLDER]
                          [--cutoff_intergenic_coverage CUTOFF_INTERGENIC_COVERAGE]
                          [--cutoff_antisense_coverage CUTOFF_ANTISENSE_COVERAGE]
                          [--cutoff_5utr_coverage CUTOFF_5UTR_COVERAGE]
                          [--cutoff_3utr_coverage CUTOFF_3UTR_COVERAGE]
                          [--cutoff_interCDS_coverage CUTOFF_INTERCDS_COVERAGE]
                          [--cutoff_background CUTOFF_BACKGROUND]
                          [--fasta_folder FASTA_FOLDER]
                          [--tex_notex_libs TEX_NOTEX_LIBS]
                          [--frag_libs FRAG_LIBS] [--tex_notex TEX_NOTEX]
                          [-rt REPLICATES_TEX] [--replicates_frag REPLICATES_FRAG]
                          [--table_best] [--sRNA_folder SRNA_FOLDER]
                          [--start_codon START_CODON] [--stop_codon STOP_CODON]
                          [--min_rbs_distance MIN_RBS_DISTANCE]
                          [--max_rbs_distance MAX_RBS_DISTANCE]
                          [--rbs_not_after_TSS] [--fuzzy_rbs FUZZY_RBS]
                          [--print_all_combination] [--best_no_sRNA]
                          [--best_no_TSS]
                          [--ignore_hypothetical_protein IGNORE_HYPOTHETICAL_PROTEIN]
                          [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --UTR_derived_sORF, -u
                            It is for detecting UTR-derived sORF. Default is
                            False.
      --transcript_assembly_folder TRANSCRIPT_ASSEMBLY_FOLDER, -a TRANSCRIPT_ASSEMBLY_FOLDER
                            The path of transcriptome assembly folder.
      --annotation_folder ANNOTATION_FOLDER, -g ANNOTATION_FOLDER
                            The path of genome annotation gff folder.
      --TSS_folder TSS_FOLDER, -t TSS_FOLDER
                            If TSS information is needed, please assign the path
                            of gff folder of TSS.
      --utr_length UTR_LENGTH, -ul UTR_LENGTH
                            If TSS information is needed, please assign the utr
                            length for comparing TSS and sORF. The default number
                            is 300.
      --min_length MIN_LENGTH, -lm MIN_LENGTH
                            Please assign the minimum residue length of sORF.
                            Default is 30.
      --max_length MAX_LENGTH, -lM MAX_LENGTH
                            Please assign the maximum residue length of sORF.
                            Default is 150.
      --tex_wig_folder TEX_WIG_FOLDER, -tw TEX_WIG_FOLDER
                            The path of TEX+/- wig folder.
      --frag_wig_folder FRAG_WIG_FOLDER, -fw FRAG_WIG_FOLDER
                            The path of fragment wig folder.
      --cutoff_intergenic_coverage CUTOFF_INTERGENIC_COVERAGE, -ci CUTOFF_INTERGENIC_COVERAGE
                            The cutoff of minimum coverage of intergenic sORF
                            candidates.
      --cutoff_antisense_coverage CUTOFF_ANTISENSE_COVERAGE, -ai CUTOFF_ANTISENSE_COVERAGE
                            The cutoff of minimum coverage of antisense sORF
                            candidates.
      --cutoff_5utr_coverage CUTOFF_5UTR_COVERAGE, -cu5 CUTOFF_5UTR_COVERAGE
                            The cutoff of minimum coverage of 5'UTR derived sORF
                            candidates. It can be assigned by percentage or the
                            amount of reads. p_0.05 means the coverage of sORF
                            candidates should higher than 5 percentile of all
                            5'UTR transcripts. n_10 means the coverage of sORF
                            candidates should be 10. Default is p_0.5.
      --cutoff_3utr_coverage CUTOFF_3UTR_COVERAGE, -cu3 CUTOFF_3UTR_COVERAGE
                            The meaning is the same as --cutoff_5utr_coverage.
                            This is for 3'UTR. Default is p_0.5.
      --cutoff_interCDS_coverage CUTOFF_INTERCDS_COVERAGE, -cuf CUTOFF_INTERCDS_COVERAGE
                            The meaning is the same as --cutoff_5utr_coverage.
                            This is for interCDS. Default is p_0.5.
      --cutoff_background CUTOFF_BACKGROUND, -cub CUTOFF_BACKGROUND
                            The cutoff of minimum coverage of all sORF candidates.
                            Default is 10.
      --fasta_folder FASTA_FOLDER, -f FASTA_FOLDER
                            The folder of genome fasta file.
      --tex_notex_libs TEX_NOTEX_LIBS, -tl TEX_NOTEX_LIBS
                            Library name of TEX+/- library. The format is:
                            wig_file_name:TEX+/-(tex or notex):condition_id(intege
                            r):replicate_id(alphabet):strand(+ or -). For multiple
                            wig files, please use comma to separate the wig files.
                            For example, wig1:tex:1:a:+,wig2:tex:1:a:-.
      --frag_libs FRAG_LIBS, -fl FRAG_LIBS
                            Library name of fragmented library The format is: wig_
                            file_name:fragmented(frag):condition_id(integer):repli
                            cate_id(alphabet):strand(+ or -). For multiple wig
                            files, please use comma to separate the wig files. For
                            example, wig1:frag:1:a:+,wig2:frag:1:a:-.
      --tex_notex TEX_NOTEX, -te TEX_NOTEX
                            For TEX+/- library, sORF should be detected in both
                            (TEX+ and TEX-) or can be detected in only one library
                            (TEX+ or TEX-). Please assign 1 or 2. Default is 2.
      -rt REPLICATES_TEX, --replicates_tex REPLICATES_TEX
                            The sORF should be detected at least in the number of
                            the replicates. The format is
                            $NUMBERofCONDITION_$NUMBERofREPLICATED. For assigning
                            different --replicates_tex to different conditions,
                            please use comma to separate it. For example,
                            1_2,2_2,3_3 means number 1 and 2 condition assign 2 to
                            --replicate_tex, and number 3 condition assign 3 to
                            --replcate_tex. For assigning the same
                            --replicates_tex to all conditions, just assign like
                            all_1 (all condition use 1 --replicate_tex).
      --replicates_frag REPLICATES_FRAG, -rf REPLICATES_FRAG
                            The meaning and input type is the same as
                            --replicates_tex.
      --table_best, -tb     The output table of sORF candidates only includes the
                            best library. Default is False.
      --sRNA_folder SRNA_FOLDER, -s SRNA_FOLDER
                            If comparing sORF and sRNA is needed, please assign
                            the path of sORF gff folder.
      --start_codon START_CODON, -ac START_CODON
                            The types of start coden. For assigning multiple types
                            of start codon, please use comma to separate them.
                            Default is ATG.
      --stop_codon STOP_CODON, -oc STOP_CODON
                            The types of stop coden. For assigning multiple types
                            of stop codon, please use comma to separate them.
                            Default is TTA,TAG,TGA.
      --min_rbs_distance MIN_RBS_DISTANCE, -mr MIN_RBS_DISTANCE
                            The minimum distance between the ribosome binding site
                            and start codon. Default is 3.
      --max_rbs_distance MAX_RBS_DISTANCE, -Mr MAX_RBS_DISTANCE
                            The maximum distance between the ribosome binding site
                            and start codon. Default is 15.
      --rbs_not_after_TSS, -at
                            This function can generate best results which also
                            include ribosome binding site not after TSS, Default
                            is False.
      --fuzzy_rbs FUZZY_RBS, -zr FUZZY_RBS
                            The number of nucleotides of ribosome binding site can
                            be different with AGGAGG? Default is 2.
      --print_all_combination, -pa
                            Non-annotated transcript may has many start codons and
                            stop codons. This function can print all combinations
                            of start codons and stop codons. Default is False.
      --best_no_sRNA, -bs   This function can generate best results which excluded
                            the candidate overlap with sRNA, please turn it on.
                            Default is False.
      --best_no_TSS, -bt    This function can generate best results without
                            referring to TSS. Default is False.
      --ignore_hypothetical_protein IGNORE_HYPOTHETICAL_PROTEIN, -ih IGNORE_HYPOTHETICAL_PROTEIN
                            For ignoring hypothetical protein in genome annotation
                            file. Default is False.

- Output files

All output files will be stored in ``$ANNOgesic/output/sORF``.

``statistics``: Statistics files.

``tables``: The tables of sORFs with more details. ``all_candidates`` is for all sORF candidates without filtering.
``best`` is for the best sORF candidates with filtering.

``gffs``: Gff files of sORFs. The meanings of ``all_candidates`` and ``best`` are the same as ``tables``.

The tags of gff files:

``start_TSS``: The starting TSS of this sORF.

``with_TSS``: Which TSSs are associated with this sORFs.

``sORF_type``: The type of the sORF (5'UTR, 3'UTR, interCDS, intergenic, antisense or within CDS).

``sRNA``: Which sRNAs are overlap with this sORFs.

``rbs``: The ribosome binding sites of this sORFs.

``frame_shift``: The number of frame shifts in the regions.

promoter
-----------

``promoter`` can scan the upstream of TSSs to discover the promoter motifs.
the regions of upstream TSSs for extraction can be assigned. We integrate 
`MEME <http://meme-suite.org/tools/meme>`_ to compute the promoters.
Visulization HTML files will be generated. If the gff files of TSSs is not computed by 
ANNOgesic, please use --TSS_source. ``promoter`` will classify the TSSs for computing 
promoter motifs.

- Pre-required tools and files

`MEME <http://meme-suite.org/tools/meme>`_.

`MPICH <https://http://www.mpich.org/>`_ (if parallel runs are required)

Fasta files of genome sequence.

Gff files of genome annotation.

Gff files of TSSs.

If the gff files of TSSs is not computed by ANNOgesic, the libraries and wiggle files are necessary.
Please refer to the ``The format of libraries for import to ANNOgesic`` in order to assign the correct format.

- Arguments

::

    usage: annogesic promoter [-h] [--MEME_path MEME_PATH]
                              [--fasta_folder FASTA_FOLDER]
                              [--TSS_folder TSS_FOLDER] [--num_motif NUM_MOTIF]
                              [--nt_before_TSS NT_BEFORE_TSS] [--e_value E_VALUE]
                              [--motif_width MOTIF_WIDTH] [--parallel PARALLEL]
                              [--TSS_source] [--tex_libs TEX_LIBS]
                              [--tex_wig_path TEX_WIG_PATH]
                              [--annotation_folder ANNOTATION_FOLDER]
                              [--combine_all]
                              [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --MEME_path MEME_PATH
                            path of MEME.
      --fasta_folder FASTA_FOLDER, -f FASTA_FOLDER
                            Please assign the folder of genome fasta file.
      --TSS_folder TSS_FOLDER, -t TSS_FOLDER
                            The folder of TSS gff file.
      --num_motif NUM_MOTIF, -n NUM_MOTIF
                            The number of motifs. Default is 10.
      --nt_before_TSS NT_BEFORE_TSS, -b NT_BEFORE_TSS
                            The number of upstream nucleotides of TSS for promoter
                            prediction? Default is 50.
      --e_value E_VALUE, -e E_VALUE
                            The cutoff of e value. Default is 0.05.
      --motif_width MOTIF_WIDTH, -w MOTIF_WIDTH
                            MEME will try to search motifs based on this length.
                            For detecting a range of width, please insert "-"
                            between two values. Moreover, for computing more than
                            one motif length, please use comma to separate them.
                            for example, 50,2-10. It means the width for detection
                            is 50 and within 2 to 10. The number should be less or
                            equal than --nt_before_TSS. Default is 50.
      --parallel PARALLEL, -pl PARALLEL
                            It is for the parallel running promoter. Please input
                            the number of parallel runs.
      --TSS_source, -s      It is for the TSS gff file which is generated from
                            other metho. It can classify TSS and generate the
                            proper format for promoter detection. Default is True
                            (from ANNOgesic)
      --tex_libs TEX_LIBS, -tl TEX_LIBS
                            Library name of TEX+/- library. If --TSS_source is
                            False, please must assign it.The format is:
                            wig_file_name:TEX+/-(tex or notex):condition_id(intege
                            r):replicate_id(alphabet):strand(+ or -). For multiple
                            wig files, please use comma to separate the wig files.
                            For example, wig1:tex:1:a:+,wig2:tex:1:a:-.
      --tex_wig_path TEX_WIG_PATH, -tw TEX_WIG_PATH
                            The path of TEX+/- wig folder. If --TSS_source is
                            False, please must assign it.
      --annotation_folder ANNOTATION_FOLDER, -g ANNOTATION_FOLDER
                            The path of genome annotation gff folder. If
                            --TSS_source is False, please must assign it..
      --combine_all, -c     It will combine all TSS files in "TSS_folder" to
                            generate global promoter motifs. Default is False.

- Output files

All output files will be stored in ``$ANNOgesic/output/promoter_analysis``.

``allfasta``: If ``--combine_all`` is True, it will combine all TSS files in ``--TSS_folder`` 
to generate promoter motifs. The results will be stored in this folder.

``fasta_class``: The fasta files of different types of TSSs.

Each strain will generate one folder for storing the output of promoter motifs.
The format of sub-folder is ``promoter_motifs_$FILENAME_$STRAINNAME_$TSSTYPE_$PROMOTERLEGNTH``. 
Each sub-folder stores PNG files and HTML files of promoters. It also contains a table (``meme.csv``) of 
prometer which can be used to sRNA detection (please refer to ``srna``).

If the TSSs are not computed by ANNOgesic, ``TSS_class`` will be generated. It will classify the 
TSSs and store as gff files.

operon
----------

``operon`` will group TSSs, genes/CDSs/tRNAs/rRNAs, transcripts, terminators and UTRs to operons and 
sub-operons.

- Pre-required tools or files

Gff files of TSSs, annotations, transcripts, 5'UTRs, and 3'UTRs.

``operon`` can integrate terminators as well, but it is not necessary.

- Arguments

::

    usage: annogesic operon [-h] [--TSS_folder TSS_FOLDER]
                            [--annotation_folder ANNOTATION_FOLDER]
                            [--transcript_folder TRANSCRIPT_FOLDER]
                            [--UTR5_folder UTR5_FOLDER]
                            [--UTR3_folder UTR3_FOLDER]
                            [--term_folder TERM_FOLDER] [--TSS_fuzzy TSS_FUZZY]
                            [--term_fuzzy TERM_FUZZY] [--min_length MIN_LENGTH]
                            [--statistics] [--combine_gff]
                            [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --TSS_folder TSS_FOLDER, -t TSS_FOLDER
                            The path of TSS gff folder.
      --annotation_folder ANNOTATION_FOLDER, -g ANNOTATION_FOLDER
                            The path of genome annotation gff folder.
      --transcript_folder TRANSCRIPT_FOLDER, -a TRANSCRIPT_FOLDER
                            The path of transcript gff folder.
      --UTR5_folder UTR5_FOLDER, -u5 UTR5_FOLDER
                            The path of 5'UTR gff folder.
      --UTR3_folder UTR3_FOLDER, -u3 UTR3_FOLDER
                            The path of 3'UTR gff folder.
      --term_folder TERM_FOLDER, -e TERM_FOLDER
                            If the information of terminator is needed, please
                            assign the path of terminator gff folder.
      --TSS_fuzzy TSS_FUZZY, -tf TSS_FUZZY
                            The fuzzy for comparing between TSS and transcript
                            assembly. Default is 5.
      --term_fuzzy TERM_FUZZY, -ef TERM_FUZZY
                            The fuzzy for comparing bewteen terminator and
                            transcript assembly. Default is 30.
      --min_length MIN_LENGTH, -l MIN_LENGTH
                            The minimum length of operon. Default is 20.
      --statistics, -s      Doing statistics for operon analysis. Default is
                            False. The name of statistics file is -
                            stat_operon_$STRAIN_NAME.csv.
      --combine_gff, -c     Convert all assigned features to one gff file. Default
                            is False.

- Output files

All output files will be stored in ``$ANNOgesic/output/operon``.

``gffs``: The gff files which integrate the information of TSSs, annotations, 
transcripts, 5'UTRs, and 3'UTRs and assign parent transcript to all features (presented by 
``Parent``).

``tables``: The tables of operons which store all information of operons and sub-operons.

``statistics``: Statistics files.

circrna
--------------

``circrna`` can detect the potential circular RNAs via `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_. 
Moreover, it remove the false positive by checking genome annotation files and quality of splicing site detection. 
User can assign reads for mapping and detecting circular RNAs or assign alignment files to skip mapping.
But BE CAREFUL, the alignment files must be mapped by `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_ 
with ``-S`` or ``circrna`` can't find the proper candidates.

- Pre-required tools and files

`Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_.

Fasta files of reads or alignment files (BAM or SAM file). If you want to use alignment files directly, they should be 
mapped by `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_ with ``-S``.

Fasta files and gff files of genome.

- Arguments

::

    usage: annogesic circrna [-h] [--segemehl_folder SEGEMEHL_FOLDER]
                             [--samtools_path SAMTOOLS_PATH] [--align]
                             [--tex_bam_path TEX_BAM_PATH]
                             [--fragmented_bam_path FRAGMENTED_BAM_PATH]
                             [--read_path READ_PATH] [--process PROCESS]
                             [--fasta_path FASTA_PATH]
                             [--annotation_path ANNOTATION_PATH]
                             [--support_reads SUPPORT_READS]
                             [--start_ratio START_RATIO] [--end_ratio END_RATIO]
                             [--ignore_hypothetical_protein IGNORE_HYPOTHETICAL_PROTEIN]
                             [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --segemehl_folder SEGEMEHL_FOLDER, -sg SEGEMEHL_FOLDER
                            Please assign the folder of segemehl.
      --samtools_path SAMTOOLS_PATH, -st SAMTOOLS_PATH
                            Please assign the path of samtools.
      --align, -a           Using segemehl to map reads (included splice
                            detection). If you already usd segemehl with -S to map
                            your reads, you can skip this step. Please be
                            attention, it only use default parameters of segemehl
                            to map the reads. Moreover, it will map all read files
                            in ANNOgesic/input/reads. For running some specific
                            functions of segemehl, please directly run segemehl.
                            Default is False.
      --tex_bam_path TEX_BAM_PATH, -tb TEX_BAM_PATH
                            If Bam files can be provided, Please assign the TEX+/-
                            Bam path.
      --fragmented_bam_path FRAGMENTED_BAM_PATH, -fb FRAGMENTED_BAM_PATH
                            If Bam files can be provided, Please assign the
                            fragmented Bam path.
      --read_path READ_PATH, -rp READ_PATH
                            If --align is True, please assign the path of reads.
      --process PROCESS, -p PROCESS
                            The number of parallels processes for --align. Default
                            is 10.
      --fasta_path FASTA_PATH, -f FASTA_PATH
                            The folder of genome fasta.
      --annotation_path ANNOTATION_PATH, -g ANNOTATION_PATH
                            The folder of genome annotation gff files.
      --support_reads SUPPORT_READS, -s SUPPORT_READS
                            The cutoff of supported reads. Default is 10.
      --start_ratio START_RATIO, -sr START_RATIO
                            The ratio of (read support circ / all read) at
                            starting point. The ratio of candidates should higher
                            than this cutoff. Default is 0.5.
      --end_ratio END_RATIO, -er END_RATIO
                            The ratio of (read support circ / all read) at end
                            point. The ratio of candidates should higher than this
                            cutoff. Default is 0.5.
      --ignore_hypothetical_protein IGNORE_HYPOTHETICAL_PROTEIN, -ih IGNORE_HYPOTHETICAL_PROTEIN
                            For ignoring hypothetical protein in genome annotation
                            file, please turn it on. Default is False.

- Output files

All the output files will be stored in ``$ANNOgesic/output/circRNA``.

``gffs``: Gff files of circular RNAs. ``$STRAINNAME_best.gff`` is the gff files for best result after comparing 
with genome annotation and quality of splicing. ``$STRAINNAME_all.gff`` is for all candidates without filering.

``circRNA_tables``: The tables for circular RNAs with more details.

``statistics``: Statistics files.

``segemehl_align``: If ``circrna`` starts from read mapping, the folder is for results of mapping.

``segemehl_splice``: The results of splicing detection. The information of the splicing tables, please 
refer to `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_.

go_term
----------

``go_term`` can retreive the information of Gene Ontology from Uniprot.
It also provides some analyses of Go terms.

- Pre-required tools and files

`idmapping_selected.tab from Uniprot <http://www.uniprot.org/downloads>`_.

`goslim.obo <http://geneontology.org/page/go-slim-and-subset-guide>`_.

`go.obo <http://geneontology.org/page/download-ontology>`_.

Gff files of genome annotation.

For detecting based on the expressed CDSs, the transcript gff file is required.

- Arguments

::

    usage: annogesic go_term [-h] [--annotation_path ANNOTATION_PATH]
                             [--transcript_path TRANSCRIPT_PATH]
                             [--UniProt_id UNIPROT_ID] [--go_obo GO_OBO]
                             [--goslim_obo GOSLIM_OBO]
                             [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --annotation_path ANNOTATION_PATH, -g ANNOTATION_PATH
                            The path of genome annotation gff folder.
      --transcript_path TRANSCRIPT_PATH, -a TRANSCRIPT_PATH
                            If the path of transcript gff folder is provided, it
                            can retrieve GO term based on expressed CDS and all
                            CDS.
      --UniProt_id UNIPROT_ID, -u UNIPROT_ID
                            The path of UniProt ID mapping database. Default is
                            ANNOgesic/input/database/idmapping_selected.tab.
      --go_obo GO_OBO, -go GO_OBO
                            The path of go.obo. Default is
                            ANNOgesic/input/database/go.obo.
      --goslim_obo GOSLIM_OBO, -gs GOSLIM_OBO
                            The path of goslim.obo. Default is
                            ANNOgesic/input/database/goslim_generic.obo.

- Output files

All output files will be stored in ``ANNOgesic/output/Go_term``.

``Go_term_results``: The tables of Go terms information.

``statistics``: Statistics files and figures.

srna_target
---------------

``srna_target`` will search the potential targets of sRNA. Different programs can
be assigned for running (RNAup or RNAplex or both). We recommand running with both 
programs. ``srna_target`` can compare the both results and provide the best ones.

- Pre-required tools and files

`ViennaRNA <http://www.tbi.univie.ac.at/RNA/>`_ .

Gff files of genome annotation.

Gff files of sRNAs.

Fasta files of genome.

- Arguments

::

    usage: annogesic srna_target [-h] [--Vienna_folder VIENNA_FOLDER]
                                 [--annotation_path ANNOTATION_PATH]
                                 [--fasta_path FASTA_PATH] [--sRNA_path SRNA_PATH]
                                 [--query_sRNA QUERY_SRNA] [--program PROGRAM]
                                 [--interaction_length INTERACTION_LENGTH]
                                 [--window_size_target WINDOW_SIZE_TARGET]
                                 [--span_target SPAN_TARGET]
                                 [--window_size_srna WINDOW_SIZE_SRNA]
                                 [--span_srna SPAN_SRNA]
                                 [--unstructured_region_RNAplex_target UNSTRUCTURED_REGION_RNAPLEX_TARGET]
                                 [--unstructured_region_RNAplex_srna UNSTRUCTURED_REGION_RNAPLEX_SRNA]
                                 [--unstructured_region_RNAup UNSTRUCTURED_REGION_RNAUP]
                                 [--energy_threshold ENERGY_THRESHOLD]
                                 [--duplex_distance DUPLEX_DISTANCE] [--top TOP]
                                 [--process_rnaplex PROCESS_RNAPLEX]
                                 [--process_rnaup PROCESS_RNAUP]
                                 [--continue_rnaup]
                                 [--potential_target_start POTENTIAL_TARGET_START]
                                 [--potential_target_end POTENTIAL_TARGET_END]
                                 [--target_feature TARGET_FEATURE]
                                 [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --Vienna_folder VIENNA_FOLDER
                            Please assign the folder of Vienna package. It should
                            include RNAplfold, RNAup and RNAplex.
      --annotation_path ANNOTATION_PATH, -g ANNOTATION_PATH
                            The path of genome annotation gff folder.
      --fasta_path FASTA_PATH, -f FASTA_PATH
                            The path of genome fasta folder.
      --sRNA_path SRNA_PATH, -r SRNA_PATH
                            The path of sRNA gff folder.
      --query_sRNA QUERY_SRNA, -q QUERY_SRNA
                            Please assign the query sRNA. For computing all sRNA
                            in gff file, please keyin 'all'.the input format
                            should be like, $STRAIN:$STRAND:$START:$END. For
                            assigning more than one sRNA, please use comma to
                            separate them.For example,
                            NC_007795.1:+:200:534,NC_007795.1:-:6767:6900. Default
                            is all.
      --program PROGRAM, -p PROGRAM
                            Using RNAplex, RNAup or both. Default is both.
      --interaction_length INTERACTION_LENGTH, -i INTERACTION_LENGTH
                            Maximum length of an interaction. Default is 30.
      --window_size_target WINDOW_SIZE_TARGET, -wt WINDOW_SIZE_TARGET
                            Only works when --program is RNAplex or both. Average
                            the pair probabilities over windows of given size for
                            RNAplex. Default is 240.
      --span_target SPAN_TARGET, -st SPAN_TARGET
                            Only works when --program is RNAplex or both. Set the
                            maximum allowed separation of a base pair to span for
                            RNAplex. Default is 160.
      --window_size_srna WINDOW_SIZE_SRNA, -ws WINDOW_SIZE_SRNA
                            Only works when --program is RNAplex or both. Average
                            the pair probabilities over windows of given size for
                            RNAplex. Default is 30.
      --span_srna SPAN_SRNA, -ss SPAN_SRNA
                            Only works when --program is RNAplex or both. Set the
                            maximum allowed separation of a base pair to span for
                            RNAplex. Default is 30.
      --unstructured_region_RNAplex_target UNSTRUCTURED_REGION_RNAPLEX_TARGET, -ut UNSTRUCTURED_REGION_RNAPLEX_TARGET
                            Only works when --program is RNAplex or both. Compute
                            the mean probability that regions of length 1 to a
                            given length are unpaired for RNAplex. Default is 30.
      --unstructured_region_RNAplex_srna UNSTRUCTURED_REGION_RNAPLEX_SRNA, -us UNSTRUCTURED_REGION_RNAPLEX_SRNA
                            Only works when --program is RNAplex or both. Compute
                            the mean probability that regions of length 1 to a
                            given length are unpaired for RNAplex. Default is 30.
      --unstructured_region_RNAup UNSTRUCTURED_REGION_RNAUP, -uu UNSTRUCTURED_REGION_RNAUP
                            Only works when --program is RNAup or both. Compute
                            the mean probability that regions of length 1 to a
                            given length are unpaired for RNAplex. Default is 30.
      --energy_threshold ENERGY_THRESHOLD, -e ENERGY_THRESHOLD
                            Only works when --program is RNAplex or both. Minimum
                            energy for a duplex to be returned for RNAplex.
                            Default is -8.
      --duplex_distance DUPLEX_DISTANCE, -d DUPLEX_DISTANCE
                            Only works when --program is RNAplex or both. Distance
                            between target 3' ends of two consecutive duplexes for
                            RNAplex. Default is 20.
      --top TOP, -t TOP     The output file only includes the candidates which
                            ranking number is smaller or equal than --top. Default
                            is 20.
      --process_rnaplex PROCESS_RNAPLEX, -pp PROCESS_RNAPLEX
                            The amount of parallel processes for running RNAplex
                            prediction. Default is 5.
      --process_rnaup PROCESS_RNAUP, -pu PROCESS_RNAUP
                            The amount of parallel processes for running RNAup
                            prediction. Default is 20.
      --continue_rnaup, -cr
                            The running time of RNAup is long, when numerous sRNAs
                            are assigned. It can run RNAup based on the previous
                            run if the process was crushed. Default is False.
      --potential_target_start POTENTIAL_TARGET_START, -ps POTENTIAL_TARGET_START
                            The number of upstream nucleotides of the start point
                            of --target_feature for extracting as potential
                            target. Default is 200.
      --potential_target_end POTENTIAL_TARGET_END, -pe POTENTIAL_TARGET_END
                            The number of downstream nucleotides of the start
                            point of --target_feature for extracting as potential
                            target. Default is 150.
      --target_feature TARGET_FEATURE, -tf TARGET_FEATURE
                            The features which are applied for extracting as
                            potential targets. If multi-features need to be
                            assigned, just use comma to separate them. Ex:
                            CDS,gene. Default is CDS.

- Output files

All output files will be stored in ``$ANNOgesic/output/sRNA_targets``.

``RNAplex``: All results of RNAplex. ``$STRAIN_RNAplex.txt`` is raw results of RNAplex.
It includes the information of binding situations. ``$STRAIN_RNAplex_rank.csv`` is the results 
that sorted by binding energy.

``RNAup``: All results of RNAup. ``$STRAIN_RNAup.txt`` is raw results of RNAup.
It includes the information of binding situations. ``$STRAIN_RNAup_rank.csv`` is the results
that sorted by binding energy.

``merge``: The results which merge ``RNAplex`` and ``RNAup``. ``$STRAIN_merge.csv`` merges all candidates of both programs. 
``$STRAIN_overlap.csv`` lists the results which be top 20 (default) in both methods.

``sRNA_seqs``: The fasta sequences of sRNAs.

``target_seqs``: The fasta sequences of potential targets.

ppi_network
-------------

``ppi_network`` will retrieve the data from `STRING <http://string-db.org/>`_. 
Then using `PIE <http://www.ncbi.nlm.nih.gov/CBBresearch/Wilbur/IRET/PIE/>`_ to search 
the literatures to support the protein-protein interaction networks. Therefore, 
``ppi_network`` can generate the protein-protein interaction networks with supported literatures.

- Pre-required tools and files

`species.vXXXX.txt from STRING <http://string-db.org/cgi/download.pl>`_.

Ptt files of genome annotation.

- Arguments

::

    usage: annogesic ppi_network [-h] [--gff_path GFF_PATH]
                                 [--proteinID_strains PROTEINID_STRAINS]
                                 [--without_strain_pubmed]
                                 [--species_STRING SPECIES_STRING] [--score SCORE]
                                 [--node_size NODE_SIZE] [--query QUERY]
                                 [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --gff_path GFF_PATH, -g GFF_PATH
                            The path of genome annotation gff folder. Please
                            confirm the gff files which have proper locus_tag item
                            in attributes. The locus_tag items can be real locus
                            tag like locus_tag=SAOUHSC_00003, or they can be
                            assigned by gene name, such as, locus_tag=dnaA.
      --proteinID_strains PROTEINID_STRAINS, -s PROTEINID_STRAINS
                            This is for query strain and filename. In order to
                            retrieve the data from STRING and Pubmed, The close
                            reference must be assigned. For example, the query
                            strain is Staphylococcus aureus HG003, but there is no
                            Staphylococcus aureus HG003 in STRING database.
                            However, there is Staphylococcus aureus 8325 which is
                            close to Staphylococcus aureus HG003. Therefore, it
                            can be used as reference. The format is like Staphyloc
                            occus_aureus_HG003.gff:Staphylococcus_aureus_HG003:"St
                            aphylococcus aureus 8325":"Staphylococcus aureus". or 
                            Staphylococcus_aureus_HG003.gff:Staphylococcus_aureus_
                            HG003:"93061":"Staphylococcus aureus". or Staphylococc
                            us_aureus_HG003.gff:Staphylococcus_aureus_HG003:"Staph
                            ylococcus aureus NCTC 8325":"Staphylococcus aureus". (
                            gff_filename:gff_strain_name:STRING_name:Pubmed_name).
                            First one is the gff file name. Second one is the
                            strain name of gff files. Third one is for STRING
                            database, and the fourth one is for Pubmed. Running
                            multiple strains at the same time is acceptable, just
                            put comma between these strains. Before running it,
                            please check the species table which should be located
                            in ANNOgesic/input/database .If the file was not
                            downloaded before, please download it. taxon_id,
                            STRING_name_compact or official_name_NCBI can be
                            assigned to STRING_name.BE CAREFUL, if the assigned
                            name with spaces, please put "" at two ends. The name
                            of Pubmed can be not so specific. If a specific name
                            is assigned, it may not be able to find the related
                            literatures.
      --without_strain_pubmed, -n
                            Retrieving pubmed without assigning strains, Default
                            is False.
      --species_STRING SPECIES_STRING, -d SPECIES_STRING
                            Please assign the path of species table of STRING.
      --score SCORE, -ps SCORE
                            Please assign the cutoff of text-mining score. The
                            value is from -1 to 1. Default is 0.
      --node_size NODE_SIZE, -ns NODE_SIZE
                            Please assign the size of nodes in figure, default is
                            4000.
      --query QUERY, -q QUERY
                            Please assign the query protein here. The format is
                            $STRAINOFGFF:$START_POINT:$END_POINT:$STRAND.For
                            multiple strains, please use comma to separate them.
                            For example, Staphylococcus_aureus_HG003:345:456:+,Sta
                            phylococcus_aureus_HG003:2000:3211:-. For computing
                            all protein, just type all. Default is all.

- Output files

All the output files will be stored in ``$ANNOgesic/output/PPI``.

``best_results``: The results which the score of `PIE <http://www.ncbi.nlm.nih.gov/CBBresearch/Wilbur/IRET/PIE/>`_
is high enough. 
``$STRAIN_without_strain.csv`` is the results of searching literatures without specific strain. 
``$STRAIN_with_strain.csv`` is the results of searching literatures with specific strain. 
For example, Staphylococcus_aureus_8325_without_strain.csv is search with Staphylococcus aureus; 
Staphylococcus_aureus_8325_without_strain.csv is search without Staphylococcus aureus. 
``without_strain`` stores all interaction information which search without specific strain. 
``with_strain`` stores all interaction information which search with specific strain. 

``all_results``: The results of all protein-protein interactions. 
(Even the text-mining(`PIE <http://www.ncbi.nlm.nih.gov/CBBresearch/Wilbur/IRET/PIE/>`_) score is too low)

``figures``: The thickness represents how many literatures can be found for the interactions. 
The solid line means there is some literatures which strongly support the interactions. The dash-dot line 
means the supported literatures are very weak. The dot line means there is no literatures which can support the 
interactions. The color is the best score of the literatures of the interactions.

subcellular_localization
------------------

``subcellular localization`` can predict the subcellular localization of CDSs. It also provides some 
statistics and visualization files.

- Pre-required tools and files

`Psortb <http://www.psort.org/psortb/>`_.

Gff files of genome annotation.

Fasta files of genome sequence.

For computing based on the expressed CDSs, the transcript gff file is required.

- Arguments

::

    usage: annogesic subcellular_localization [-h] [--Psortb_path PSORTB_PATH]
                                              [--gff_path GFF_PATH]
                                              [--fasta_path FASTA_PATH]
                                              [--transcript_path TRANSCRIPT_PATH]
                                              [--bacteria_type BACTERIA_TYPE]
                                              [--difference_multi DIFFERENCE_MULTI]
                                              [--merge_to_gff]
                                              [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --Psortb_path PSORTB_PATH
                            If you want to assign the path of Psortb, please
                            assign here.
      --gff_path GFF_PATH, -g GFF_PATH
                            The path of genome annotation gff folder.
      --fasta_path FASTA_PATH, -f FASTA_PATH
                            The path of genome fasta folder.
      --transcript_path TRANSCRIPT_PATH, -a TRANSCRIPT_PATH
                            If the path of transcript gff folder is provided, it
                            can get the results based on expressed CDS and all
                            CDS.
      --bacteria_type BACTERIA_TYPE, -b BACTERIA_TYPE
                            The type of bacteria (Gram-positive or Gram-negative).
                            Please assign 'positive' or 'negative'.
      --difference_multi DIFFERENCE_MULTI, -d DIFFERENCE_MULTI
                            For the protein which have multiple locations, if the
                            difference of psorb scores is smaller than
                            --difference_multi, it will be printed out as well.
                            Default is 0.5. The maximum value is 10.
      --merge_to_gff, -m    Merging the information to genome annotation gff file.
                            Default is False.

- Output files

All output files will be stored in ``$ANNOgesic/output/subcellular_localization``.

``psortb_results``: The results of Psortb.

``statistics``: Statistics files and figures.

riboswitch_thermometer
--------------

``riboswitch_thermometer`` will search ribosome binding sites in the region between 
a TSS（the starting point of transcript was assigned if no TSS was detected) 
and a downstream CDSs. Then using `Infernal <http://infernal.janelia.org/>`_ to scan 
riboswitch or RNA thermometer in `Rfam <http://rfam.xfam.org/>`_.

- Pre-required tools and files

`Infernal <http://infernal.janelia.org/>`_.

`Rfam <http://rfam.xfam.org/>`_.

Gff files of genome annotation.

Fasta files of genome sequence.

File of ``riboswitch_ID`` or ``RNA_thermometer_ID``. The file should contain Rfam ID, 
Name and Description of riboswitchs or RNA thermometer (refer to 
``Riboswitch and RNA thermometer dataset of Rfam``). You can download the file from our 
`Github <https://github.com/Sung-Huan/ANNOgesic>`_ (Rfam_riboswitch_ID.csv or Rfam_RNA_thermometer_ID.csv).

- Arguments

::

    usage: annogesic riboswitch_thermometer [-h] [--program PROGRAM]
                                            [--infernal_path INFERNAL_PATH]
                                            [--riboswitch_ID RIBOSWITCH_ID]
                                            [--RNA_thermometer_ID RNA_THERMOMETER_ID]
                                            [--gff_path GFF_PATH]
                                            [--tss_path TSS_PATH]
                                            [--UTR_length UTR_LENGTH]
                                            [--transcript_path TRANSCRIPT_PATH]
                                            [--fasta_path FASTA_PATH]
                                            [--Rfam RFAM] [--e_value E_VALUE]
                                            [--output_all] [--fuzzy FUZZY]
                                            [--fuzzy_rbs FUZZY_RBS]
                                            [--start_codon START_CODON]
                                            [--max_dist_rbs MAX_DIST_RBS]
                                            [--min_dist_rbs MIN_DIST_RBS]
                                            [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --program PROGRAM, -p PROGRAM
                            Please assign the query feature. The options are
                            "riboswitch", "thermometer", "both". Default is both.
      --infernal_path INFERNAL_PATH, -if INFERNAL_PATH
                            Please assign the folder of Infernal (where is cmscan
                            and cmsearch located).
      --riboswitch_ID RIBOSWITCH_ID, -ri RIBOSWITCH_ID
                            If --program is "riboswitch" or "both", please assigh
                            the path of the riboswitch ID of Rfam. The file should
                            include the Accession, ID and Description of
                            riboswitch in Rfam.
      --RNA_thermometer_ID RNA_THERMOMETER_ID, -ti RNA_THERMOMETER_ID
                            If --program is "thermometer" or "both", please assigh
                            the path of the RNA thermometer ID of Rfam. The file
                            should include the Accession, ID and Description of
                            RNA thermometer in Rfam.
      --gff_path GFF_PATH, -g GFF_PATH
                            The path of annotation gff folder.
      --tss_path TSS_PATH, -t TSS_PATH
                            The path of TSS gff folder.
      --UTR_length UTR_LENGTH, -u UTR_LENGTH
                            The length of UTR. Default is 300.
      --transcript_path TRANSCRIPT_PATH, -a TRANSCRIPT_PATH
                            The path of transcript gff folder.
      --fasta_path FASTA_PATH, -f FASTA_PATH
                            The path of genome fasta folder.
      --Rfam RFAM, -R RFAM  The path of Rfam CM database.
      --e_value E_VALUE, -e E_VALUE
                            The cutoff of e value. Default is 0.001.
      --output_all, -o      One query sequence may fit multiple riboswitches or
                            RNA thermometers. If --output_all is True, all results
                            all be printed out. Or it will only print the best
                            one. Default is False.
      --fuzzy FUZZY, -z FUZZY
                            The fuzzy of nucleotides for 3' or 5' end. Default is
                            10.
      --fuzzy_rbs FUZZY_RBS, -zr FUZZY_RBS
                            The number of nucleotides of ribosome binding site can
                            be different with AGGAGG. Default is 2.
      --start_codon START_CODON, -ac START_CODON
                            The types of start coden. If more than one type are
                            assigned, please use comma to separate them. Default
                            is ATG.
      --max_dist_rbs MAX_DIST_RBS, -Mr MAX_DIST_RBS
                            The maximum distance between ribosome binding site and
                            start codon. Default is 14.
      --min_dist_rbs MIN_DIST_RBS, -mr MIN_DIST_RBS
                            The minmum distance between ribosome binding site and
                            start codon. Default is 5.

- Output files

All output files of riboswitch will be stored in ``$ANNOgesic/output/riboswitch`` and 
all output files of RNA thermometer will be stored in ``$ANNOgesic/output/RNA_thermometer``.

``gffs``: The gff files of riboswitchs / RNA_thermometer.

``tables``: The tables of riboswichs / RNA_thermometer with more details.

``scan_Rfam``: The results of searching to Rfam with ``cmscan`` (`Infernal <http://infernal.janelia.org/>`_). 
User can check the results of blast.

``statistics``: The statistics files and figures.

crispr
---------------
``crispr`` integrates CRISPR Recognition Tool (`CRT <http://www.room220.com/crt/>`_) which can detect the repeat 
units and spacer of CRISPR. It can also compare to genome annotation to exclude some false positives.

- Pre-required tools and files

`CRT <http://www.room220.com/crt/>`_.

Fasta file of genome sequence.

Some useful information can be assigned to improve the prediction:

Gff files of genome annotation.

- Arguments

::

    usage: annogesic crispr [-h] [--CRT_path CRT_PATH] [--gff_path GFF_PATH]
                            [--fasta_path FASTA_PATH] [--window_size WINDOW_SIZE]
                            [--min_number_repeat MIN_NUMBER_REPEAT]
                            [--min_length_repeat MIN_LENGTH_REPEAT]
                            [--Max_length_repeat MAX_LENGTH_REPEAT]
                            [--min_length_spacer MIN_LENGTH_SPACER]
                            [--Max_length_spacer MAX_LENGTH_SPACER]
                            [--ignore_hypothetical_protein]
                            [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --CRT_path CRT_PATH   If you want to assign the path of CRT.jar (the execute
                            file), please assign here. Default is
                            /usr/local/bin/CRT.jar
      --gff_path GFF_PATH, -g GFF_PATH
                            If genome gff file is provided, it will compare with
                            CRISPR for removing the false positives. Default is
                            None.
      --fasta_path FASTA_PATH, -f FASTA_PATH
                            The path of genome fasta folder.
      --window_size WINDOW_SIZE, -w WINDOW_SIZE
                            Length of window size for searching (range: 6-9).
                            Default is 8.
      --min_number_repeat MIN_NUMBER_REPEAT, -mn MIN_NUMBER_REPEAT
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

- Output files

It will generate three folders - ``CRT_output``, ``gffs`` and ``statistics``. 

``CRT_output`` stores the output of `CRT <http://www.room220.com/crt/>`_.

``gffs`` stores the gff files of CRISPR. ``all_candidates`` is for the all CRISPRs which ANNOgesic can detect. 
``best`` is for the CRISPRs which are not overlap with genome annotations.

``statistics`` is for the statistics files.


optimize_tsspredator
---------------

``optimize_tsspredator`` can adapt the parameter set of `TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_. 
For running it, only manual checking around 200kb is required (using gff format).
We suggest the manual detected TSSs should more than 50. If there are less than 50 TSSs within 200kb, 
please continue checking until 50 TSSs detected.
Then ``optimize_tsspredator`` can scan whole genome based on the principle of manual detection to get the optimized parameters.

- Pre-required tools and files

`TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_.

The libraries and wiggle files of TEX +/-. 
Please refer to the ``The format of libraries for import to ANNOgesic``.

Fasta file of genome sequence.

Gff file of genome annotation.

Gff file of manual detection.

- Arguments

::

    usage: annogesic optimize_tsspredator [-h]
                                          [--TSSpredator_path TSSPREDATOR_PATH]
                                          [--fasta_file FASTA_FILE]
                                          [--annotation_file ANNOTATION_FILE]
                                          [--wig_folder WIG_FOLDER]
                                          [--manual MANUAL]
                                          [--strain_name STRAIN_NAME]
                                          [--max_height MAX_HEIGHT]
                                          [--max_height_reduction MAX_HEIGHT_REDUCTION]
                                          [--max_factor MAX_FACTOR]
                                          [--max_factor_reduction MAX_FACTOR_REDUCTION]
                                          [--max_base_height MAX_BASE_HEIGHT]
                                          [--max_enrichment_factor MAX_ENRICHMENT_FACTOR]
                                          [--max_processing_factor MAX_PROCESSING_FACTOR]
                                          [--utr_length UTR_LENGTH] [--lib LIB]
                                          [--output_prefix OUTPUT_PREFIX]
                                          [--cluster CLUSTER] [--length LENGTH]
                                          [--core CORE] [--program PROGRAM]
                                          [--replicate_match REPLICATE_MATCH]
                                          [--steps STEPS]
                                          [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --TSSpredator_path TSSPREDATOR_PATH
                            If you want to assign the path of TSSpredator, please
                            assign here. Default is /usr/local/bin/TSSpredator.jar
      --fasta_file FASTA_FILE, -fs FASTA_FILE
                            Path of one target genome fasta file.
      --annotation_file ANNOTATION_FILE, -g ANNOTATION_FILE
                            Path of one target genome annotation gff file.
      --wig_folder WIG_FOLDER, -w WIG_FOLDER
                            The folder of the TEX+/- wig folder.
      --manual MANUAL, -m MANUAL
                            The file of manual checked gff file.
      --strain_name STRAIN_NAME, -n STRAIN_NAME
                            The name of the strain for optimization.
      --max_height MAX_HEIGHT, -he MAX_HEIGHT
                            This value relates to the minimum number of read
                            starts at a certain genomic position to be considered
                            as a TSS candidate. During optimization --max_height
                            will be never larger than this value. Default is 2.5.
      --max_height_reduction MAX_HEIGHT_REDUCTION, -rh MAX_HEIGHT_REDUCTION
                            When comparing different strains/conditions and the
                            step height threshold is reached in at least one
                            strain/condition, the threshold is reduced for the
                            other strains/conditions by the value set here. This
                            value must be smaller than the step height threshold.
                            During optimization, --max_height_reduction will be
                            never larger than this value. Default is 2.4.
      --max_factor MAX_FACTOR, -fa MAX_FACTOR
                            This is the minimum factor by which the TSS height has
                            to exceed the local expression background. During
                            optimization, --max_factor will be never larger than
                            this value. Default is 10.
      --max_factor_reduction MAX_FACTOR_REDUCTION, -rf MAX_FACTOR_REDUCTION
                            When comparing different strains/conditions and the
                            step factor threshold is reached in at least one
                            strain/condition, the threshold is reduced for the
                            other strains/conditions by the value set here. This
                            value must be smaller than the step factor threshold.
                            During optimization, --max_factor_reduction will be
                            never larger than this value. Default is 9.9.
      --max_base_height MAX_BASE_HEIGHT, -bh MAX_BASE_HEIGHT
                            This is the minimum number of reads should be mapped
                            on TSS. During optimization, --max_base_height will be
                            never larger than this value. Default is 0.06.
      --max_enrichment_factor MAX_ENRICHMENT_FACTOR, -ef MAX_ENRICHMENT_FACTOR
                            This is the minimum enrichment factor. During
                            optimization, --max_enrichment_factor will be never
                            larger than this value. Default is 6.0.
      --max_processing_factor MAX_PROCESSING_FACTOR, -pf MAX_PROCESSING_FACTOR
                            This is the minimum processing factor. If untreated
                            library is higher than the treated library and above
                            which the TSS candidate is considered as a processing
                            site and not annotated as detected. During
                            optimization, --max_processing_factor will be never
                            larger than this value. Default is 6.0
      --utr_length UTR_LENGTH, -u UTR_LENGTH
                            The length of UTR. It is for Primary and Secondary
                            TSSs. Default is 300.
      --lib LIB, -l LIB     The libraries of TEX+/- wig files for TSSpredator. The
                            format is: wig_file_name:TEX+/-(tex or notex):conditio
                            n_id(integer):replicate_id(alphabet):strand(+ or -).
                            For multiple wig files, please use comma to separate
                            the wig files. For example,
                            wig1:tex:1:a:+,wig2:tex:1:a:-.
      --output_prefix OUTPUT_PREFIX, -p OUTPUT_PREFIX
                            The output prefix of all conditions. For multiple
                            conditions, please use comma to separate them. For
                            example, prefix_condition1,prefix_condition2.
      --cluster CLUSTER, -cu CLUSTER
                            If the position between manual one and predicted one
                            is smaller or equal than this value, it will only
                            print one of them. Default is 2.
      --length LENGTH, -le LENGTH
                            The length of genome for running optimization. If
                            running for whole genome is not the case, please
                            assign the length of partial genome. Default is
                            comparing whole genome.
      --core CORE, -c CORE  Paralle runs for optimization. Default is 4.
      --program PROGRAM, -t PROGRAM
                            The feature for optimization. Please assign "TSS" or
                            "Processint_site". Default is TSS.
      --replicate_match REPLICATE_MATCH, -rm REPLICATE_MATCH
                            The TSS candidates should be detected at lease in the
                            number of replicates. The format is
                            $NUMBERofCONDITION_$NUMBERofREPLICATED. For assigning
                            different --replicate_match to different conditions,
                            please use comma to separate it. For example,
                            1_2,2_2,3_3 means number 1 and 2 condition assign 2 to
                            --replicate_match, and number 3 condition assign 3 to
                            --replcate_match. For assigning the same
                            --replicate_match to all conditions, just assign like
                            all_1 (all condition use 1 --replicate_match). Default
                            is all_1.
      --steps STEPS, -s STEPS
                            The total runs for optimization. Default is 4000 runs.

- Output files

Based on the programs (TSS/processing site), the output files will be stored in 
``$ANNOgesic/output/TSS/optimized_TSSpredator`` or ``$ANNOgesic/output/processing_site/optimized_TSSpredator``.

``stat.csv`` stores the information of every step. The first column is the number of run.
The second column is the parameter set. ``he`` represents height; ``rh`` represents 
height reduction; ``fa`` means factor; ``rf`` means factor reduction; ``bh`` is 
base height; ``ef`` is enrichment factor; ``pf`` is processing factor. About the details 
of parameters, please refer to `TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_.
For example, ``he_2.0_rh_1.8_fa_4.4_rf_2.8_bh_0.08_ef_3.0_pf_2.6`` means height is 2.0, 
height reduction is 1.8, factor is 4.4, factor reduction is 2.8, base height is 0.08, 
enrichment factor is 3.0 and processing factor is 2.6. The third and fourth columns are 
the number of true positive. The fifth and sixth columns are true positive rate. The seventh 
and eighth columns are the number of false positive. The ninth and tenth false positive rate. 
The eleventh and twelfth columns are the number of false negative. The thirteenth and fourteenth 
columns are missing rate.

``best.csv`` stores the best parameter set. The meanings of all columns are the same as ``stat.csv``.

screenshot
-----------

``screenshot`` will generate batch files for producing screenshot of `IGV <https://www.broadinstitute.org/igv>`_.
It can reduce the time for checking the results in genome browser.
When the batch files is produced, user just needs to open `IGV <https://www.broadinstitute.org/igv>`_, then presses ``tools`` 
on the top tags and choose ``run batch script``. The program will automatically produce screenshots. 

- Pre-required tools and files

`IGV <https://www.broadinstitute.org/igv>`_.

Gff files that user want to produce screenshots. All screenshots will be produced based on the positions of ``main_gff``. 
``side_gffs`` are the gff files that user want to compare with ``main_gff``. They will also appear in the screenshots.

Fasta files of genome.

The libraries and wiggle files. Please refer to the ``The format of libraries for import to ANNOgesic``.

- Arguments

::

    usage: annogesic screenshot [-h] [--main_gff MAIN_GFF] [--side_gffs SIDE_GFFS]
                                [--fasta FASTA]
                                [--frag_wig_folder FRAG_WIG_FOLDER]
                                [--tex_wig_folder TEX_WIG_FOLDER]
                                [--height HEIGHT] [--tex_libs TEX_LIBS]
                                [--frag_libs FRAG_LIBS] [--present PRESENT]
                                [--output_folder OUTPUT_FOLDER]
                                [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --main_gff MAIN_GFF, -mg MAIN_GFF
                            Screenshot will based on the position of main_gff file
                            to generate screenshot.
      --side_gffs SIDE_GFFS, -sg SIDE_GFFS
                            If comparing between other features and --main_gff is
                            needed, please assign the other gff files here. For
                            multiple gff files, just use comma to separate them.
      --fasta FASTA, -f FASTA
                            The path of genome fasta folder.
      --frag_wig_folder FRAG_WIG_FOLDER, -fw FRAG_WIG_FOLDER
                            If the information of fragmented wig file is required,
                            please assign the folder.
      --tex_wig_folder TEX_WIG_FOLDER, -tw TEX_WIG_FOLDER
                            If the information of TEX+/- wig file is required,
                            please assign the folder.
      --height HEIGHT, -he HEIGHT
                            The height of screenshot. Default is 1500.
      --tex_libs TEX_LIBS, -tl TEX_LIBS
                            If the TEX+/- wig file is required, please also assign
                            proper format here. The format is:
                            wig_file_name:TEX+/-(tex or notex):condition_id(intege
                            r):replicate_id(alphabet):strand(+ or -). For multiple
                            wig files, please use comma to separate the wig files.
                            For example, wig1:tex:1:a:+,wig2:tex:1:a:-.
      --frag_libs FRAG_LIBS, -fl FRAG_LIBS
                            If the fragmented wig file is required, please also
                            assign proper format here. The format is: wig_file_nam
                            e:fragmented(frag):condition_id(integer):replicate_id(
                            alphabet):strand(+ or -). For multiple wig files,
                            please use comma to separate the wig files. For
                            example, wig1:frag:1:a:+,wig2:frag:1:a:-.
      --present PRESENT, -p PRESENT
                            The type of presentation in the screenshot.
                            expand/collapse/squish. Default is expand.
      --output_folder OUTPUT_FOLDER, -o OUTPUT_FOLDER
                            Please assign the output folder. It will create a sub-
                            folder "screenshots" in the output folder to store the
                            results.

- Output files

Based on the paths of ``main_gff``, ``screenshot`` will generate a folder - ``screenshots`` under the 
folder of ``main_gff``. All output files will be generated in this folder.

``forward.txt`` is the batch file of forward strand for running on IGV.

``reverse.txt`` is the batch file of reverse strand for running on IGV.

``forward`` is the folder for storing screenshots of forward strand.

``reverse`` is the folder for storing screenshots of reverse strand.

When batch files are executed on IGV, the screenshots will automatically store in ``forward`` and ``reverse``. 
The format of filenames will be ``$STRAIN:$START-$END.png``. For example, ``NC_007795:1051529-1051696.png`` 
means the strain is NC_007795. The feature's start point is 1051529 and the end point is 
1051696.

color_png
----------

``color_png`` is a following procedure of ``screenshot``. If numerous samples are included in one figure, 
it will be difficult to distinguish the tracks. ``color_png`` can color the tracks based on TEX +/- libraries 
for assisting the checking process.

- Pre-required tools and files

`ImageMagick <http://www.imagemagick.org/script/index.php>`_.

The screenshots which generated from ``screenshot``. Please make sure the folders of ``forward`` and ``reverse`` 
exist in the folder of ``screenshots``.

- Arguments

::

    usage: annogesic color_png [-h] [--screenshot_folder SCREENSHOT_FOLDER]
                               [--track_number TRACK_NUMBER]
                               [--ImageMagick_covert_path IMAGEMAGICK_COVERT_PATH]
                               [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --screenshot_folder SCREENSHOT_FOLDER, -f SCREENSHOT_FOLDER
                            The folder which stores "screenshots" (a folder
                            generated by subcommand "screenshot").
      --track_number TRACK_NUMBER, -t TRACK_NUMBER
                            The number of tracks.
      --ImageMagick_covert_path IMAGEMAGICK_COVERT_PATH, -m IMAGEMAGICK_COVERT_PATH
                            Please assign the path of "convert" in ImageMagick
                            package.

- Output files

The new screenshots will replace the previous ones automatically.

merge_features
--------------
If merging multiple features to one gff file is needed, ``merge_features`` can achieve this purpose. 
It will merge all features that user assigned and search the parent transcript to each feature.


- Arguments

::

    usage: annogesic merge_features [-h] [--transcript_path TRANSCRIPT_PATH]
                                    [--other_features_path OTHER_FEATURES_PATH]
                                    [--fuzzy_term FUZZY_TERM]
                                    [--fuzzy_TSS FUZZY_TSS]
                                    [--strain_name STRAIN_NAME]
                                    [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given, the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --transcript_path TRANSCRIPT_PATH, -a TRANSCRIPT_PATH
                            If transcript gff file is provided. The output will be
                            generated based on transcripts and assign
                            "Parent_tran" to each feature in attributes of gff
                            file. If there is no transcript, the output will just
                            simply combine all input gff files.
      --other_features_path OTHER_FEATURES_PATH, -of OTHER_FEATURES_PATH
                            Please assign the gff files of other features which
                            you want to merge.
      --fuzzy_term FUZZY_TERM, -fm FUZZY_TERM
                            For merging terminators, please assign the fuzzy
                            between transcript and terminator. ATTENTION, the
                            third column of gff file of terminator should be
                            exactly "terminator". Default is 30.
      --fuzzy_TSS FUZZY_TSS, -ft FUZZY_TSS
                            For merging TSSs, please assign the fuzzy between TSS
                            and transcript. ATTENTION, the third column of gff
                            file of terminator should be exactly "TSS". Default is
                            5.
      --strain_name STRAIN_NAME, -s STRAIN_NAME
                            Please assign the strain name of the input files. It
                            will become the prefix name of output gff file.

-Output files

The gff file will be stored in ``merge_all_features``. The tag - ``Parent`` in attributes of 
gff file shows the parent transcript.
