ANNOgesic's subcommands
===============

In general the subcommands need at least one argument - the analysis
folder. If this is not given ANNOgesic assumes that the current
folder is the analysis folder.

The format of filename
--------------------
In order to recognize the type of files and the relationship between different files,
please follow the principle of the format of filename.

Keep the filename of genome fasta file the same as annotation file. For example,
``NC_007795.fa, NC_007795.gff, NC_007795.ptt, NC_007795.rnt, NC_007705.gbk``.

All the ``.gff`` files which are import into ``ANNOgesic``, please follow the format.
``$STRAINNAME_$FEATURE.gff``. For example, ``NC_007795_TSS.gff, NC_007795_transcript.gff``.

The names of features are:

============  ===========================
Feature name  meaning
------------  --------------------------- 
TSS           transcription starting site
processing    processing site
transcript    transcript assembly
sRNA          small RNA
sORF          small open reading frame
term          terminator
5UTR          5'UTR
3UTR          3'UTR
============  ===========================

For wiggle files, please use ``forward`` and ``reverse`` to represent strands.
The format of wiggle file is ``$LIBRARYNAME_$STRAND.wig``. 
For example, ``TSB_OD_0.2_forward.wig``.

Please also notice the name of strain or genome. It should be avoid to use ``|`` in the name. Because ``|`` is one character which already used in Unix system. However, if your fasta file or gff file are download from NCBI, ANNOgesic will automatically modify the name.

The format of libraries for import to ANNOgesic
-----------------------------

Some modules of ``ANNOgesic`` need to import the information of libraries.
For import the information of libraries, please follow this format:

``$LIBRARY_FILENAME:$LIBRARY_TYPE:$CONDITION:$REPLICATE:$STRAND``

``$LIBRARY_FILENAME`` means the ``.wig`` file.

``$LIBRARY_TYPE`` can be ``tex`` (TEX+) or ``notex`` (TEX-) or ``frag`` (fragmented, not for ``tsspredator``).

``$CONDITION`` is the index of conditions. Please use 1, 2, 3, ... to represent different conditions.

``$REPLICATE`` is the index of replicates. Please use a, b, c, ... to represent different replicates.

``$STRAND`` is the strand of wiggle file. Please use + or -.

For example, 

for TEX +/- treated libraries:

::

  TSB_OD_0.2_TEX_reverse.wig:tex:1:a:- 
  TSB_OD_0.5_TEX_reverse.wig:tex:2:a:- 
  TSB_OD_0.2_TEX_forward.wig:tex:1:a:+ 
  TSB_OD_0.5_TEX_forward.wig:tex:2:a:+ 
  TSB_OD_0.2_reverse.wig:notex:1:a:- 
  TSB_OD_0.5_reverse.wig:notex:2:a:- 
  TSB_OD_0.2_forward.wig:notex:1:a:+ 
  TSB_OD_0.5_forward.wig:notex:2:a:+

The above block is one file. Each line is combine with space.

for fragmented libraries:

::

  fragmented_forward.wig:frag:1:a:+ fragmented_reverse.wig:frag:1:a:-


create
-----

``create`` generates the required folder structure for input and
output files. Once these folders are created the input files have to
be placed into the correct locations. After creating the folders,
please put the required files in the proper folders.

BAMs: For ``.bam`` files

database: For all databases. For example: BSRD, nr...

manual_TSS: If you detected transcription starting sites(TSS) manually,
you can put the results here. When you compute TSS, ANNOgesic 
will merge them  together. If you want to run ``TSS_optimization``  
for TSS. It also required.

manual_processing_site: It is similar with ``manual_TSS``, it is for 
processing sites.

mutation_table: If you detect the mutation between reference genome and 
target genome manually, please put the file here. When
you run ``get_target_fasta``. it is required. Please refer
to the section of ``get_target_fasta`` for the format of 
mutation table.

promoter_analysis: Please leave it blank. It is for store the fasta for
``promoter_analysis``.

reads: If you want to run ``circrna`` and you also want to realign read data,
please put the read data here. It can also deal with ``.bzip2`` and ``.gzip``.
       
reference: For annotation files and fasta files. If the information of 
the reference strain can be download from NCBI, you can also get
the files through running ``get_input_files``.

riboswitch_ID: For store the file which contains all the Rfam ID of riboswitch.
For the details of format, please refer to the section of 
``riboswicth``.

TSSpredator: Please leave it blank. It is for config files of ``TSSpredator``.

wigs: For wiggle files. Based on the type of RNA-Seq, you can put them to 
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

``get_input_files`` is the subcommand for download required files (fasta, annotation files)from NCBI. 
Therefore, user need to assign the IP of the reference genome. For example,
ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Staphylococcus_aureus_NCTC_8325_uid57795
Then, user can assign which kinds of files he/she wants to download.


- Pre-reqired information

``FTP source``: The IP of NCBI where store all the data of your reference strain.

- Arguments


::

    usage: annogesic get_input_files [-h] [--FTP_path FTP_PATH] [--ref_fasta]
                                     [--ref_gff] [--ref_ptt] [--ref_rnt]
                                     [--ref_gbk] [--convert_embl]
                                     [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --FTP_path FTP_PATH, -F FTP_PATH
                            Path of website where can download the required files
      --ref_fasta, -f       Download fasta files of reference. Default is False
      --ref_gff, -g         Download gff files of reference. Default is False
      --ref_ptt, -p         Download ptt files of reference. Default is False
      --ref_rnt, -r         Download rnt files of reference. Default is False
      --ref_gbk, -k         Download genbank files of reference. Default is False
      --convert_embl, -e    Convert gbk to embl files of reference. Default is
                            False

- Output files

The output files will store in ``$ANNOgesic_folder/input/reference``.

``fasta``: fasta files.

``annotation``: annotation files.

get_target_fasta
--------------

``get_target_fasta`` is the subcommand for generating target fasta file from 
reference genome files. It is based on the mutation table to modify the reference 
fasta to target fasta. Therefore, the similarity of reference genome and target genome
should be close. The example of format of mutation table is following:

===========  ============  ============  ========  =========  ====================  =============  ====  ============
 #target_id  reference_id  reference_nt  position  target_nt  impact of correction  locus tag      gene  Description 
-----------  ------------  ------------  --------  ---------  --------------------  -------------  ----  ------------
 HG003       NC_007795.1   a             333       c                                SAOUHSC_00002  dnaA  XXXXXX      
 HG003       NC_007795.1   t             543       \-          deletion                                  YYYYYY      
 HG003       NC_007795.1   \-            600       g           insertion            SAOUHSC_00132                    
===========  ============  ============  ========  =========  ====================  =============  ====  ============

If user wants to put the name of column in the top, it needs to start from ``#``. 
Each column is separated by ``tab``. If the mutation type if deletion or insertion, 
user can put - to represent them. The information of target_id, reference_id,
reference_nt, position, target_nt is required. The others can leave them blank. 
However, please still use tab to separate all blank columns.

If user has no mutation information between the reference genome and target 
genome, user can also use ``SNP_calling`` (one module of ``ANNOgesic``) to compute 
it. Please refer to the description of SNP_calling.

- Pre-required files

Fasta file of reference genome.

Mutation table which indicate the information of mutations between reference and target genome.

- Arguments

::

    usage: annogesic get_target_fasta [-h] [--ref_fasta_folder REF_FASTA_FOLDER]
                                  [--mutation_table MUTATION_TABLE]
                                  [--output_format OUTPUT_FORMAT [OUTPUT_FORMAT ...]]
                                  [project_path]

    positional arguments:
      project_path          Path of the project folder. If none is given the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --ref_fasta_folder REF_FASTA_FOLDER, -r REF_FASTA_FOLDER
                            Folder of fasta files.
      --mutation_table MUTATION_TABLE, -m MUTATION_TABLE
                            The path of mutation table.
      --output_format OUTPUT_FORMAT [OUTPUT_FORMAT ...], -o OUTPUT_FORMAT [OUTPUT_FORMAT ...]
                            Please assign the output filename and which strain
                            should be included in it. For example:
                            FILE1:strain1,strain2. FILE1 is a output fasta file
                            which include the information of strain1 and strain2.

- Output files

Fasta files of target genome will store in ``$ANNOgesic_folder/output/target/fasta``

annotation_transfer
-----------

``annotation transfer`` is the subcommand for transfer the annotation from reference genome 
to target genome. In this module, we use `PAGIT and RATT <http://www.sanger.ac.uk/resources/software/pagit/>`_ 
to achieve it. The similarity of reference genome and target genome should be closed.
Or it will influence the final results.
Please be attation, before you start to run RATT(annotation transfer), 
run ``source $PAGIT_HOME/sourceme.pagit`` first. it will modify the path for execute RATT.

- Pre-required tools and files

`PAGIT and RATT <http://www.sanger.ac.uk/resources/software/pagit/>`_

The genbank file of reference genome.

The fasta file of reference genome.

The fasta file of target genome.

- Arguments

::

    usage: annogesic annotation_transfer [-h] [--RATT_path RATT_PATH]
                                         [--compare_pair COMPARE_PAIR [COMPARE_PAIR ...]]
                                         [--element ELEMENT]
                                         [--transfer_type TRANSFER_TYPE]
                                         [--ref_gbk REF_GBK]
                                         [--ref_fasta REF_FASTA]
                                         [--target_fasta TARGET_FASTA]
                                         [--convert_to_gff_rnt_ptt]
                                         [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --RATT_path RATT_PATH
                            Path of the start.ratt.sh file of RATT folder. Default
                            is start.ratt.sh.
      --compare_pair COMPARE_PAIR [COMPARE_PAIR ...], -p COMPARE_PAIR [COMPARE_PAIR ...]
                            Please assign the name of fasta files of pairs. ex.
                            NC_007795:NEW_NC_007795. The reference fasta file is
                            NC_007795.fa and the target fasta file is
                            NEW_NC_007795.fa. ATTENTION:please make sure the ref
                            name is the same as embl file.
      --element ELEMENT, -e ELEMENT
                            It will become the prefix of all output file.
      --transfer_type TRANSFER_TYPE, -t TRANSFER_TYPE
                            The transfer type for running RATT.(details can refer
                            to the manual of RATT.) default is Strain.
      --ref_gbk REF_GBK, -re REF_GBK
                            The folder which stores every reference embl
                            folders.If you have no embl folder, you can assign the
                            folder of genbank.
      --ref_fasta REF_FASTA, -rf REF_FASTA
                            The folder of reference fasta files.
      --target_fasta TARGET_FASTA, -tf TARGET_FASTA
                            The folder which stores target fasta files.
      --convert_to_gff_rnt_ptt, -g
                            Do you want to convert to gff, rnt and ptt? Default is
                            False

- Output files

All the output files from `PAGIT and RATT <http://www.sanger.ac.uk/resources/software/pagit/>`_
will store in ``$ANNOgesic_folder/output/annotation_transfer``.

All annotation files(``.gff``, ``.ptt``, ``.rnt``) will store in ``$ANNOgesic_folder/output/target/annotation``.

snp
-------

``snp`` can detect the potential mutation positions by comparing the results of alignment and fasta file.
`Samtools <https://github.com/samtools>`_, `Bcftools <https://github.com/samtools>`_ are the main tools
for detect mutations. User can choose programs (with BAQ, without BAQ and extend BAQ) to run ``snp``.
User can also set the quality, read depth and the fraction of maximum read depth which support for indel.
User can use it for getting the target fasta file from the alignment results of reference genome.
User can also use it for checking the alignment results of target genome.

- Pre-required files and tools:

`Samtools <https://github.com/samtools>`_

`Bcftools <https://github.com/samtools>`_

BAM files for fragmented libraries or TEX +/- treated libraries.

Reference or target fasta files.

- Arguments

::

   usage: ANNOgesic.py snp [-h] [--samtools_path SAMTOOLS_PATH]
                           [--bcftools_path BCFTOOLS_PATH] [--bam_type BAM_TYPE]
                           [--program PROGRAM [PROGRAM ...]]
                           [--fasta_path FASTA_PATH]
                           [--tex_bam_path TEX_BAM_PATH]
                           [--frag_bam_path FRAG_BAM_PATH] [--quality QUALITY]
                           [--read_depth READ_DEPTH]
                           [--indel_fraction INDEL_FRACTION]
                           [project_path]
   
   positional arguments:
     project_path          Path of the project folder. If none is given the
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
                           Please assign the type of BAM. If your BAM file is
                           mapping to reference genome and you want to know the
                           difference between refenece genome and target genome,
                           plase keyin 'reference'. If your BAM file already
                           mapped to target genome and you want to check the
                           genome sequence has SNP or not, please keyin 'target'
     --program PROGRAM [PROGRAM ...], -p PROGRAM [PROGRAM ...]
                           Please assign the program for detecting SNP of
                           transcript: 1: calculate with BAQ, 2: calculate
                           without BAQ, 3: calculate with extend BAQ. You can
                           assign more than 1 program. For example: 1 2 3
     --fasta_path FASTA_PATH, -f FASTA_PATH
                           The path of fasta folder.
     --tex_bam_path TEX_BAM_PATH, -tw TEX_BAM_PATH
                           The path of tex+/- wig folder. If you want to use tex
                           treated and untreated bam files, please assign the
                           path. Or it will not combine the bam files
     --frag_bam_path FRAG_BAM_PATH, -fw FRAG_BAM_PATH
                           The path of fragmented wig folder. If you want to use
                           fragmented bam files, please assign the path. Or it
                           will not combine the bam files
     --quality QUALITY, -q QUALITY
                           The min quality which consider a real snp. Default is
                           20
     --read_depth READ_DEPTH, -d READ_DEPTH
                           The minimum read depth, below to it will be excluded.
                           default is 5 * number of BAM files,if the cutoff
                           higher than 40, it will use 40.
     --indel_fraction INDEL_FRACTION, -imf INDEL_FRACTION
                           The fraction of maximum read depth, which support
                           insertion of deletion. Default is 0.5

- Output files

The results will store according to the type of Bam files. If the Bam files are from mapping
with reference, the results will store in ``$ANNOgesic/output/SNP_calling/compare_reference``. 
If the Bam files are from mapping with target genome, the results will store in 
``$ANNOgesic/output/SNP_calling/validate_target``.

The raw data from `Samtools <https://github.com/samtools>`_ and `Bcftools <https://github.com/samtools>`_
will store in ``$ANNOgesic/output/SNP_calling/$BAM_TYPE/SNP_raw_outputs``.

The results which fit the conditions which user set will store in 
``$ANNOgesic/output/SNP_calling/$BAM_TYPE/SNP_table``.

The meaning of file names are:

``$STRAIN_$PROGRAM_depth_only.vcf`` means the results only match the condition of read depth.

``$STRAIN_$PROGRAM_depth_quality.vcf`` means the results match the condition of read depth and quality.

``$STRAIN_$PROGRAM_seq_reference.csv`` is the index of fasta files which are generated based on the results of ``snp``.

For example,

::

  Staphylococcus_aureus_HG003     1632629 .       AaA     AA      57      .
  Staphylococcus_aureus_HG003     1632630 .       aA      a       57      .
  Staphylococcus_aureus_HG003     1499572 .       T       TT,TTTTT        43.8525 .

These first two mutation will cause conflicts. Then the conflicts will effect other mutations.
Therefore, based on these information of mutations, it will generate four different fasta files.
``$STRAIN_$PROGRAM_seq_reference.csv`` is the index for these four fasta files.

::

   1       1632629 1       1499572:TT      Staphylococcus_aureus_HG003
   1       1632629 2       1499572:TTTTT   Staphylococcus_aureus_HG003
   2       1632630 1       1499572:TT      Staphylococcus_aureus_HG003
   2       1632630 2       1499572:TTTTT   Staphylococcus_aureus_HG003

The first column is the index of conflicts. The second column is the position which be selected.
The third one is the index of two potential mutations in the same position. The fourth one is
the position and nucleotides of mutation. The last column is the name of strain.
If you refer to ``$ANNOgesic/output/SNP_calling/$BAM_TYPE/seqs``, the filename of fasta is like 
``$FILENAME_$STRIANNAME_$INDEXofCONFLICT_$INDEXofTWOMUTATION.fa``. Therefore, the first line of 
``$STRAIN_$PROGRAM_seq_reference.csv`` will generate 
``Staphylococcus_aureus_HG003_Staphylococcus_aureus_HG003_1_1.fa`` 
(if the file name of genome is Staphylococcus_aureus_HG003). The second line will generate
``Staphylococcus_aureus_HG003_Staphylococcus_aureus_HG003_1_2.fa`` and so forth.

The statistics files will store in ``$ANNOgesic/output/SNP_calling/$BAM_TYPE/statistics``.

expression_analysis
--------------
``expression_analysis`` can analyze which CDSs or other features express in which libraries.
It can be used to compare between different conditions. It is also a good way to detect housekeeping gene.

- Pre-required tools and files

The gff file of annotation.

The library and wiggle file. Please refer to the ``The format of libraries for import to ANNOgesic`` in 
order to assign correct format.

- Arguments

::

    usage: ANNOgesic.py expression_analysis [-h] [-g ANNOTATION_FOLDER]
                                            [-tl TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]]
                                            [-fl FRAG_LIBS [FRAG_LIBS ...]]
                                            [-te TEX_NOTEX] [-rt REPLICATES_TEX]
                                            [-rf REPLICATES_FRAG]
                                            [--tex_wig_folder TEX_WIG_FOLDER]
                                            [--frag_wig_folder FRAG_WIG_FOLDER]
                                            [--cutoff_overlap_tex CUTOFF_OVERLAP_TEX]
                                            [--cutoff_overlap_frag CUTOFF_OVERLAP_FRAG]
                                            [--cutoff_coverage CUTOFF_COVERAGE]
                                            [--features FEATURES [FEATURES ...]]
                                            [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      -g ANNOTATION_FOLDER, --annotation_folder ANNOTATION_FOLDER
                            The folder of annotation file which you want to
                            analyze.
      -tl TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...], --tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]
                            Library name of tex and notex library. The format is:
                            wig_file_name:tex_treat_or_not(tex or notex):condition
                            _id(integer):replicate_id(alphabet):strand(+ or -).
      -fl FRAG_LIBS [FRAG_LIBS ...], --frag_libs FRAG_LIBS [FRAG_LIBS ...]
                            Library name of fragment library. The format is: wig_f
                            ile_name:fragmented(frag):condition_id(integer):replic
                            ate_id(alphabet):strand(+ or -).
      -te TEX_NOTEX, --tex_notex TEX_NOTEX
                            For tex +/- library, expressing CDS should be detected
                            by both or just one.(1 or 2)
      -rt REPLICATES_TEX, --replicates_tex REPLICATES_TEX
                            The expressing CDS of tex +/- library should be
                            detected more than this number of replicates.
      -rf REPLICATES_FRAG, --replicates_frag REPLICATES_FRAG
                            The expressing CDS of fragmented library should be
                            detected more than this number of replicates.
      --tex_wig_folder TEX_WIG_FOLDER, -tw TEX_WIG_FOLDER
                            The folder of TEX+/- wigge files.
      --frag_wig_folder FRAG_WIG_FOLDER, -fw FRAG_WIG_FOLDER
                            The folder of fragmented wigge files.
      --cutoff_overlap_tex CUTOFF_OVERLAP_TEX, -ot CUTOFF_OVERLAP_TEX
                            This value is for decision of CDS which is expressing
                            or not in TEX+/- library. If the expressing nts more
                            than this value, it will consider the CDS is
                            expressing one. You can assign by percentage or
                            nucleotide. ex: p_0.5 means the percentage of
                            expressing nts should higher 0.5. n_100 means there
                            should be 100 nts which are expressing. Default is
                            "all" which means as long as there is a nt's coverage
                            higher than cutoff_coverage, it would consider the CDS
                            which is expressing.
      --cutoff_overlap_frag CUTOFF_OVERLAP_FRAG, -of CUTOFF_OVERLAP_FRAG
                            This value is for decision of CDS which is expressing
                            or not in fragmented library. If the expressing nts
                            more than this value, it will consider the CDS is
                            expressing. You can assign by percentage or
                            nucleotide. ex: p_0.5 means the percentage of
                            expressing nts should higher 0.5. n_100 means there
                            should be 100 nts which are expressing. Default is
                            "all" which means as long as there is a nt's coverage
                            higher than cutoff_coverage, it would consider the CDS
                            which is expressing.
      --cutoff_coverage CUTOFF_COVERAGE, -c CUTOFF_COVERAGE
                            If the coverage is higher than this value, it will
                            consider the nt is expressing
      --features FEATURES [FEATURES ...], -f FEATURES [FEATURES ...]
                            The features which you want to compute, ex: CDS tRNA

- Output files

All output files will store in ``$ANNOgesic/output/target/annotation/for_libs`` of the input folder of annotation. 

``gffs``: gff files. Based on the format of libraries, the name of output annotation file would be 
``STRAIN_FEATURE_LIBTYPE.gff``. For example, ``Staphylococcus_aureus_HG003_CDS_10_texnotex.gff``.
It means the strain is Staphylococcus_aureus_HG003. The feature of this gff file is for analysis of CDS.
About the ``LIBTYPE``, ``10_texnotex`` is the number of condition of tex treated libraries.
``all_libs`` means the CDSs or other features express in all libraries. ``at_least_one_lib`` 
means the CDSs or other features express at least in one libraries. ``no_express`` means 
the CDSs or other features don't express in any libraries.

``statistics``: statistics files.

tsspredator(TSS and processing site prediction)
--------------

``tsspredator`` can generate the candidates of TSS and processing sites. The main tool is
`TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_. We can easily switch the
TEX + libraries and TEX - libraries to detect processing sites. User can assign the parameters 
for `TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_. If User want to get the 
optimized parameters of `TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_,
there is ``optimize_tsspredator`` for this purpose. Please refer to this section.

For import the information of libraries, please refer to the section 
``The format of libraries for import to ANNOgesic``.

- Pre-required tools and files

`TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_

The libraries and wiggle files of TEX +/-. Please refer to the ``The format of libraries for import to ANNOgesic``.

The fasta and anntation file of genome.

If User has gff file of TSS by manual detection, ``tsspredator`` can merge the manual one
and predicted one.

If User want to compare TSS with transcripts, it also need the gff file of transcripts.
For the transcripts, please refer to the section of ``transcript_assembly``.

- Arguments

::

   usage: ANNOgesic.py tsspredator [-h] [--TSSpredator_path TSSPREDATOR_PATH]
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
                                   [--utr_length UTR_LENGTH]
                                   [--lib LIB [LIB ...]]
                                   [--output_prefix OUTPUT_PREFIX [OUTPUT_PREFIX ...]]
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
     project_path          Path of the project folder. If none is given the
                           current directory is used.
   
   optional arguments:
     -h, --help            show this help message and exit
     --TSSpredator_path TSSPREDATOR_PATH
                           If you want to assign the path of TSSpredator, please
                           assign here.
     --fasta_folder FASTA_FOLDER, -f FASTA_FOLDER
                           Path of the target fasta folder.
     --annotation_folder ANNOTATION_FOLDER, -g ANNOTATION_FOLDER
                           Path of the target gff folder.
     --wig_folder WIG_FOLDER, -w WIG_FOLDER
                           The folder of the wig folder.
     --height HEIGHT, -he HEIGHT
                           This value relates to the minimal number of read
                           starts at a certain genomic position to be considered
                           as a TSS candidate. Default is 0.3
     --height_reduction HEIGHT_REDUCTION, -rh HEIGHT_REDUCTION
                           When comparing different strains/conditions and the
                           step height threshold is reached in at least one
                           strain/condition, the threshold is reduced for the
                           other strains/conditions by the value set here. This
                           value must be smaller than the step height threshold.
                           Default is 0.2
     --factor FACTOR, -fa FACTOR
                           This is the minimal factor by which the TSS height has
                           to exceed the local expression background. Default is
                           2.0
     --factor_reduction FACTOR_REDUCTION, -rf FACTOR_REDUCTION
                           When comparing different strains/conditions and the
                           step factor threshold is reached in at least one
                           strain/condition, the threshold is reduced for the
                           other strains/conditions by the value set here. This
                           value must be smaller than the step factor threshold.
                           Default is 0.5
     --enrichment_factor ENRICHMENT_FACTOR, -ef ENRICHMENT_FACTOR
                           This is the minimal enrichment factor. During
                           optimization will never larger than this value.
                           Default is 2.0
     --processing_factor PROCESSING_FACTOR, -pf PROCESSING_FACTOR
                           This is the minimal processing factor. If untreated
                           library is higher than the treated library and above
                           which the TSS candidate is considered as a processing
                           site and not annotated as detected. During
                           optimization will never larger than this value.
                           Default is 1.5
     --base_height BASE_HEIGHT, -bh BASE_HEIGHT
                           This is the minimal number of reads should be mapped
                           on TSS. Default is 0.0
     --replicate_match REPLICATE_MATCH, -rm REPLICATE_MATCH
                           The TSS candidates should match to how many number of
                           the replicates. Default is 1.
     --utr_length UTR_LENGTH, -u UTR_LENGTH
                           The length of UTR. It is for Primary and Secondary
                           definition. Default is 300
     --lib LIB [LIB ...], -l LIB [LIB ...]
                           The libraries of wig files for TSSpredator. The format
                           is: wig_file_name:tex_treat_or_not(tex or notex):condi
                           tion_id(integer):replicate_id(alphabet):strand(+ or
                           -).
     --output_prefix OUTPUT_PREFIX [OUTPUT_PREFIX ...], -p OUTPUT_PREFIX [OUTPUT_PREFIX ...]
                           The output prefix of all conditions.
     --merge_manual MERGE_MANUAL, -m MERGE_MANUAL
                           If you have gff file of manual checked TSS, you can
                           use this function to merge manual checked ones and
                           predicted ones.
     --statistics, -s      Doing statistics for TSS candidates. it will store in
                           statistics folder. Default is False
     --validate_gene, -v   Using TSS candidates to validate genes in annotation
                           file. it will store in statistics folder. Default is
                           False
     --compute_program COMPUTE_PROGRAM, -t COMPUTE_PROGRAM
                           Which program do you want to predict. (TSS or
                           processing_site)
     --compare_transcript_assembly COMPARE_TRANSCRIPT_ASSEMBLY, -ta COMPARE_TRANSCRIPT_ASSEMBLY
                           If you want to compare with transcriptome assembly,
                           please assign the folder of gff file of transcript
                           assembly.Default is False.
     --fuzzy FUZZY, -fu FUZZY
                           The fuzzy for comparing TSS and transcript assembly.
                           Default is 5
     --cluster CLUSTER, -c CLUSTER
                           This number is for compare manual detected TSS and
                           prediced one. If the position between manual checked
                           one and predicted one is smaller or equal than this
                           value, It will only print one of them. Default is 2
     --length LENGTH, -le LENGTH
                           The length that you want to compare with manual check
                           for statistics. If you want to compare whole genome,
                           please don't turn it on. The default is comparing
                           whole genome
     --re_check_orphan, -ro
                           If your annotation file lack information of gene or
                           locus_tag, you can turn it on. It will try to compare
                           with CDS. Default is False
     --overlap_feature OVERLAP_FEATURE, -of OVERLAP_FEATURE
                           If processing site and TSS are overlap, you can keep
                           "TSS" or "processing_site" or "both". Default is both.
     --reference_gff_folder REFERENCE_GFF_FOLDER, -rg REFERENCE_GFF_FOLDER
                           For overlap_feature, if you want to only keep "TSS" or
                           "processing_site", please assign the
                           reference_gff_folder. If you are running TSS, please
                           assign the folder of processing site. If you are
                           running processing_site, please assign the folder of
                           TSS. If you want to keep "both" at overlap position,
                           please don't turn it on.
     --remove_low_expression REMOVE_LOW_EXPRESSION, -rl REMOVE_LOW_EXPRESSION
                           If you want to remove low expressed TSS/processing
                           site, please assign the file of manual checked gff
                           file here. Please Be ATTENTION: this parameter may
                           remove some True positive, too. So, please make sure
                           you want to do it.

- Output files

The output files will be stored in ``$ANNOgesic/output/TSS``.

``MasterTables``: The output files from `TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_.

``gffs``: The gff file of TSS.
The second column of gff file shows the TSS is from manual detection or
`TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_.

There are some useful tags of in the attributes of gff file:

``type``: It represents the type of TSS. It could be Primary, Secondary, Internal, Antisense or Orphan.

``UTR_length``: It represents the length of UTR.

``associated_gene``: It shows the TSS associates which gene.

``Parent_tran``: if user has compared with transcript, it will show that the TSS is located in which transcript.

If user has compared with genome annotation file, the tag - ``start_TSS`` will appear in the gff file 
of genome annotation. It will show the TSS of CDS/rRNA/tRNA.

If user has compared with transcript, the tag - ``associated_tss`` will appear in the gff file
of transcript. It will show the associated TSS in the transcript.

``statistics``: statistics files.

The output file of processing site is similar. Just replace ``TSS`` to ``processing_site``
like ``$ANNOgesic/output/processing_site``.

transcript_assembly
-------------------

``transcript_assembly`` will base on the coverage to generate transcripts. User can 
assign the parameters of ``transcript_assembly``.

For importing the information of libraries, please refer to the section 
``The format of libraries for import to ANNOgesic``.

- Pre-required tools and files

Wiggle files of fragmented libraries or TEX +/- treated libraries. We don't recommend only
use wiggle files of TEX +/- treated libraries to generate transcripts. It will lose some reads
in 3'end.

If user wants to compare transcript with TSS, it requires ``.gff`` file of TSS.
If user wants to compare transcript with genome anntation, it requires ``.gff`` file of genome.

- Arguments

::

   usage: ANNOgesic.py transcript_assembly [-h]
                                           [--annotation_folder ANNOTATION_FOLDER]
                                           [--sort_annotation] [--length LENGTH]
                                           [--normal_wig_path NORMAL_WIG_PATH]
                                           [--frag_wig_path FRAG_WIG_PATH]
                                           [--height HEIGHT] [--width WIDTH]
                                           [--tolerance TOLERANCE]
                                           [--tolerance_coverage TOLERANCE_COVERAGE]
                                           [--replicates_tex REPLICATES_TEX]
                                           [--replicates_frag REPLICATES_FRAG]
                                           [--tex_notex TEX_NOTEX]
                                           [--compare_TSS COMPARE_TSS]
                                           [--compare_CDS COMPARE_CDS]
                                           [--TSS_fuzzy TSS_FUZZY]
                                           [--Tex_treated_libs TEX_TREATED_LIBS [TEX_TREATED_LIBS ...]]
                                           [--fragmented_libs FRAGMENTED_LIBS [FRAGMENTED_LIBS ...]]
                                           [project_path]
   
   positional arguments:
     project_path          Path of the project folder. If none is given the
                           current directory is used.
   
   optional arguments:
     -h, --help            show this help message and exit
     --annotation_folder ANNOTATION_FOLDER, -g ANNOTATION_FOLDER
                           It is for comparing transcript assembly and annotation
                           gff file. It can use annotation gff file as reference
                           and modify transcript assembly file. If you want to do
                           it, please assign the annotation gff folder.
                           Otherwise, don't turn it on.
     --sort_annotation, -s
                           The annotation gff files in annotation folder are
                           sorted or not. If they didn't be sorted, please turn
                           it on. Default is False
     --length LENGTH, -l LENGTH
                           The minimum width of region to be a transcript. It is
                           for refer to annotation file. If you want to compare
                           with annotation files, it will be the final output. If
                           you don't want to compare with annotation files,
                           --width would be the length for the final output. The
                           default is 20.
     --normal_wig_path NORMAL_WIG_PATH, -nw NORMAL_WIG_PATH
                           The path of normal wig folder.
     --frag_wig_path FRAG_WIG_PATH, -fw FRAG_WIG_PATH
                           The path of fragment wig folder.
     --height HEIGHT, -he HEIGHT
                           The minimum height of coverage to be a transcript. The
                           default is 5.
     --width WIDTH, -w WIDTH
                           The minimum width of region to be a transcript. It is
                           for without annotation to be reference. If you don't
                           want to compare with annotation files (--length), it
                           will be the final output. Otherwise, --length would be
                           the length of transcript for the final output.The
                           default is 20.
     --tolerance TOLERANCE, -t TOLERANCE
                           This number indicates how willing the algorithm is to
                           ignore a temporary drop below this number. The default
                           is 5.
     --tolerance_coverage TOLERANCE_COVERAGE, -tc TOLERANCE_COVERAGE
                           If the coverage is lower than tolerance_coverage, even
                           the range is within tolerance, it will terminate the
                           current transcript. The default is 0.
     --replicates_tex REPLICATES_TEX, -rt REPLICATES_TEX
                           The position is included in the transcript if there
                           are more than the replicate which you assign here to
                           supported it. (for tex +/- library)
     --replicates_frag REPLICATES_FRAG, -rf REPLICATES_FRAG
                           The position is included in the transcript if there
                           are more than the replicate which you assign here to
                           supported it. (for fragmented library)
     --tex_notex TEX_NOTEX, -te TEX_NOTEX
                           If you use tex +/- libraries to run transcript
                           assembly, please assign the tex +/- should both
                           consider or just one. (1 or 2). Default is 2
     --compare_TSS COMPARE_TSS, -ct COMPARE_TSS
                           If you want to compare with TSS, please assign TSS
                           folder.
     --compare_CDS COMPARE_CDS, -cg COMPARE_CDS
                           If you want to compare with annotation file, please
                           assign annotation folder.
     --TSS_fuzzy TSS_FUZZY, -fu TSS_FUZZY
                           The fuzzy for comparing TSS and transcript assembly.
                           Default is 5
     --Tex_treated_libs TEX_TREATED_LIBS [TEX_TREATED_LIBS ...], -tl TEX_TREATED_LIBS [TEX_TREATED_LIBS ...]
                           Input of tex +/- library. The format is:
                           wig_file_name:tex_treat_or_not(tex or notex):condition
                           _id(integer):replicate_id(alphabet):strand(+ or -).
     --fragmented_libs FRAGMENTED_LIBS [FRAGMENTED_LIBS ...], -fl FRAGMENTED_LIBS [FRAGMENTED_LIBS ...]
                           Input of fragmented library. The format is: wig_file_n
                           ame:fragmented(frag):condition_id(integer):replicate_i
                           d(alphabet):strand(+ or -).

- Output files

The output files will be stored in ``$ANNOgesic/output/transcriptome_assembly``.

``gffs``: The gff files of transcript.
In the second column of gff file shows the transcript is from which kinds of libraries.

There are some useful tags in gff files.

``type``: It shows the situation of overlap between transcript and CDS/tRNA/rRNA
(cover_CDS, within_CDS, not_related_CDS, left_shift_CDS or right_shift_CDS).
(If user has compared transcript with genome annotation.) 

``associated_tss``: It shows the TSS which are located in the transcript. 
(If user has compared transcript with TSS.) 

``associated_cds``: It shows the CDS/tRNA/rRNA which are located in the transcript.
(If user has compared transcript with genome annotation.) 

If user has compared transcript with genome annotation. The tag - ``Parent_tran`` will appear
in the gff file of genome annotation. It will show the CDS/tRNA/rRNA is located in which transcript.

If user has compared transcript with TSS. The tag - ``Parent_tran`` will appear
in the gff file of TSS. It will show the TSS is located in which transcript.

``statistics``: statistics files.

terminator
-----------

``terminator`` will predict the rho-independent terminator. ``ANNOgesic`` combine the results of 
two methods in order to get more reliable candidates. First one is using `TranstermHP <http://transterm.cbcb.umd.edu/>`_.
The other one is detect the specific secondary structure in intergenic region of forward and reverse strand. 
``ANNOgesic`` can also compare with coverage in order to generate the terminators which has dramatic coverage
decreasing.

- Pre-required tools and files

`TranstermHP <http://transterm.cbcb.umd.edu/>`_

RNAfold of `ViennaRNA <http://www.tbi.univie.ac.at/RNA/>`_

Annotation file and fasta file of target genome

Wiggle file of TEX +/- treated libraries or fragmented libraries. we don't 
recommand only use TEX +/- treated libraries. Because it will lose reads in 3'end.

Gff file of transcript

- Arguments

::

    usage: ANNOgesic.py terminator [-h] [--TransTermHP_path TRANSTERMHP_PATH]
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
                                   [--fuzzy_upstream_transcript FUZZY_UPSTREAM_TRANSCRIPT]
                                   [--fuzzy_downstream_transcript FUZZY_DOWNSTREAM_TRANSCRIPT]
                                   [--fuzzy_upstream_cds FUZZY_UPSTREAM_CDS]
                                   [--fuzzy_downstream_cds FUZZY_DOWNSTREAM_CDS]
                                   [--highest_coverage HIGHEST_COVERAGE]
                                   [-tl TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]]
                                   [-fl FRAG_LIBS [FRAG_LIBS ...]] [-te TEX_NOTEX]
                                   [-rt REPLICATES_TEX] [-rf REPLICATES_FRAG]
                                   [-tb]
                                   [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --TransTermHP_path TRANSTERMHP_PATH
                            Please assign the path of transterm in TransTermHP.
      --expterm_path EXPTERM_PATH
                            Please assign the path of your expterm.dat.
      --RNAfold_path RNAFOLD_PATH
                            If you want to assign the path of RNAfold of Vienna
                            package, please assign here.
      --fasta_folder FASTA_FOLDER, -f FASTA_FOLDER
                            The path of fasta folder.
      --annotation_folder ANNOTATION_FOLDER, -g ANNOTATION_FOLDER
                            The path of annotation gff folder.
      --transcript_folder TRANSCRIPT_FOLDER, -a TRANSCRIPT_FOLDER
                            The path of the folder which store gff files of
                            transcript assembly.
      --sRNA SRNA, -sr SRNA
                            If you want to include sRNA information, please assign
                            the folder of gff files of sRNA.
      --statistics, -s      Doing statistics for TransTermHP. the name of
                            statistics file is - stat_terminator_$STRAIN_NAME.csv.
      --tex_wig_folder TEX_WIG_FOLDER, -tw TEX_WIG_FOLDER
                            If you want to use tex +/- libraries, please assign
                            tex +/- wig folder.
      --frag_wig_folder FRAG_WIG_FOLDER, -fw FRAG_WIG_FOLDER
                            If you want to use fragmented libraries, please assign
                            fragmented wig folder.
      --decrease DECREASE, -d DECREASE
                            If the (lowest coverage / highest coverage) in the
                            terminator is smaller than this number, it will
                            consider this terminator have dramatic coverage
                            decrease in it.
      --fuzzy_detect_coverage FUZZY_DETECT_COVERAGE, -fc FUZZY_DETECT_COVERAGE
                            It will elongate the number of nt(you assign here)
                            from both terminal site. If it can found the coverage
                            dramatic decrease within this range, it will consider
                            the terminator have dramatic coverage decrease in it.
      --fuzzy_upstream_transcript FUZZY_UPSTREAM_TRANSCRIPT, -fut FUZZY_UPSTREAM_TRANSCRIPT
                            If the candidates are upstream of transcript and the
                            distance between the end of gene and terminator
                            candidate is within this number, it will be consider
                            as terminator.
      --fuzzy_downstream_transcript FUZZY_DOWNSTREAM_TRANSCRIPT, -fdt FUZZY_DOWNSTREAM_TRANSCRIPT
                            If the candidates are downstream of transcript and the
                            distance is within this number, it will be consider as
                            terminator.
      --fuzzy_upstream_cds FUZZY_UPSTREAM_CDS, -fuc FUZZY_UPSTREAM_CDS
                            If the candidates are upstream of CDS/tRNA/rRNA/sRNA
                            and the distance between the end of gene and
                            terminator candidate is within this number, it will be
                            consider as terminator.
      --fuzzy_downstream_cds FUZZY_DOWNSTREAM_CDS, -fdg FUZZY_DOWNSTREAM_CDS
                            If the candidates are downstream of CDS/tRNA/rRNA/sRNA
                            and the distance is within this number, it will be
                            consider as terminator.
      --highest_coverage HIGHEST_COVERAGE, -hc HIGHEST_COVERAGE
                            If the highest coverage of the region of terminator is
                            below to this number, the terminator will be classify
                            to non-detect.
      -tl TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...], --tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]
                            Library name of tex and notex library. The format is:
                            wig_file_name:tex_treat_or_not(tex or notex):condition
                            _id(integer):replicate_id(alphabet):strand(+ or -).
      -fl FRAG_LIBS [FRAG_LIBS ...], --frag_libs FRAG_LIBS [FRAG_LIBS ...]
                            Library name of fragment library. The format is: wig_f
                            ile_name:fragmented(frag):condition_id(integer):replic
                            ate_id(alphabet):strand(+ or -).
      -te TEX_NOTEX, --tex_notex TEX_NOTEX
                            For tex +/- library, terminators should be detected by
                            both or just one.(1/2)
      -rt REPLICATES_TEX, --replicates_tex REPLICATES_TEX
                            The terminator of tex +/- library should be detected
                            more than this number of replicates.
      -rf REPLICATES_FRAG, --replicates_frag REPLICATES_FRAG
                            The terminator of fragmented library should be
                            detected more than this number of replicates.
      -tb, --table_best     Output sRNA table only most decreasing track.

- Output files

The output files will stored in ``$ANNOgesic/output/terminator``.

``gffs``: gff files of terminator.
There are three different folders to store terminators.

``all_candidates`` will store all terminators which ``ANNOgesic`` can detect.

``express`` will store the terminators which has gene expression

``detect`` will store the terminators which not only has gene expression but also
the coverage has dramatic decreasing.

The tags of gff files:

``coverage_decrease``: The coverage of the terminator has dramatic decreasing or not.

``express``: The terminator has gene expression or not.

``diff_coverage``: The best library which can detect the terminator. The numbers in parens 
are highest coverage and lowest coverage.

``tables``: the table of terminator which store more details.

``statistics``: statistics files.

``transtermhp``: all output of `TranstermHP <http://transterm.cbcb.umd.edu/>`_.

utr
-----

``utr`` can compare with the information of TSS, CDS/tRNA/sRNA, transcripts and terminators
to generate proper UTRs. 5'UTR is based on detecting the intergenic region of  TSS (which 
are located in transcript) and CDS/tRNA/sRNA. 3'UTR is based on detecting the intergenic
region of the terminal of transcript and CDS/tRNA/sRNA. If the gff file of TSS is not computed by 
ANNOgesic, please use --TSS_source. ``utr`` would compute the class of TSS for analysis.

- Pre-required files

Gff files of genome annotation, TSS and transcript.

If user want to combine the information of terminator, it also need the gff file of terminator.

- Arguments

::

   usage: ANNOgesic.py utr [-h] [--annotation_folder ANNOTATION_FOLDER]
                           [--TSS_folder TSS_FOLDER]
                           [--transcript_assembly_folder TRANSCRIPT_ASSEMBLY_FOLDER]
                           [--terminator_folder TERMINATOR_FOLDER]
                           [--terminator_fuzzy TERMINATOR_FUZZY] [--TSS_source]
                           [--base_5UTR BASE_5UTR]
                           [project_path]
   
   positional arguments:
     project_path          Path of the project folder. If none is given the
                           current directory is used.
   
   optional arguments:
     -h, --help            show this help message and exit
     --annotation_folder ANNOTATION_FOLDER, -g ANNOTATION_FOLDER
                           The path of annotation gff folder.
     --TSS_folder TSS_FOLDER, -t TSS_FOLDER
                           The path of TSS folder.
     --transcript_assembly_folder TRANSCRIPT_ASSEMBLY_FOLDER, -a TRANSCRIPT_ASSEMBLY_FOLDER
                           The path of transcriptome assembly folder.
     --terminator_folder TERMINATOR_FOLDER, -e TERMINATOR_FOLDER
                           If you want to add the information of terminator, you
                           can assign the path of terminator folder here.
     --terminator_fuzzy TERMINATOR_FUZZY, -f TERMINATOR_FUZZY
                           If the distance(nt) between terminator and the end of
                           transcript assembly belows to this value, it will
                           assign the terminator associated with the 3'UTR.
     --TSS_source, -s      If you generate TSS from other method not from
                           ANNOgesic, please turn it on.
     --base_5UTR BASE_5UTR, -b BASE_5UTR
                           Which kind of information that you want to use for
                           generating 5'UTR. TSS/transcript/both. Default is
                           both.

- Output files

All output of 5'UTR will store in ``$ANNOgesic/output/UTR/5UTR``.

All output of 3'UTR will store in ``$ANNOgesic/output/UTR/3UTR``.

``gffs``: gff files of 5'UTR/3'UTR

The tags of gff file:

``length``: UTR length.

``associated_cds``: Which CDS/rRNA/tRNA are associated with this UTR.

``associated_gene``: Which genes are associated with this UTR.

``associated_tss``: Which TSS are associated with this 5'UTR.

``TSS_type``: What types of TSS are associated with this 5'UTR.

``associated_tran``: Which transcript is associated with this 3'UTR. 

``associated_term``: Which terminators are associated with this 3'UTR.

srna
-----
``srna`` can predict sRNA candidates by comparing the transcripts and annotation profile. 
The transcripts in intergenic region might be sRNA candidates. Moreover, based on 
the information of TSS and processing site, we can also predict UTR derived sRNA candidates.

- Pre-required tools and files

`ViennaRNA <http://www.tbi.univie.ac.at/RNA/>`_

`Ps2pdf14 <http://pages.cs.wisc.edu/~ghost/doc/AFPL/6.50/Ps2pdf.htm>`_

`Blast+ <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`_

`BSRD <http://www.bac-srna.org/BSRD/index.jsp>`_

`nr database <ftp://ftp.ncbi.nih.gov/blast/db/FASTA/>`_

Gff files of genome annotation and Transcript assembly data.

It also can import more information to improve the accuracy of prediction.

TSS: To filter out the sRNA candidates which have no TSS related.

Processing site: For checking the end point of sRNA candidates which are located in long transcripts.

wiggle file: To detect the intergenic transcripts which have dramatic coverage 
decreasing, it could be find the sRNA which located in long transcripts,
The libraries and wiggle files, Please refer to the ``The format of libraries for import to ANNOgesic``.

sRNA database: It could be used to search the known sRNA. The format of header should be ``$ID|$STRAIN|$sRNANAME|$OTHER_INFO|$OTHER_INFO``. 
For example, ``>saci2813.1|Acinetobacter sp. ADP1|Aar|12240|2812430|forward``.
The ID is saci403.1; the strain of this sRNA is Acinetobacter sp. ADP1 and the name of sRNA is Aar. 
``srna`` only consider the first three columns. Therefore, the format of the first three columns should be 
fit the rule. BE ATTENTION, after format the sRNA database, the information after the third columns will be removed.
If the user doesn't follow the format, it will occur some error when the user run ``--sRNA_blast_stat, -sb``. 
Of course, it also can run ``srna`` without ``--sRNA_blast_stat, -sb``.

nr database: It could be used to search the known protein in order to exclude false positive.

sORF: It could compare sORF and sRNA. User can refer these information to find the best sRNA candidates.

If user want to detect the UTR derived sRNA, it will be necessary to import 
more information.

TSS: It becomes necessary information. Because UTR derived sRNAs must to from TSS.

wiggle file: It becomes necessary information. Because the terminal of UTR derived 
sRNAs should be processing site or dramatic coverage decreasing point.
The libraries and wiggle files, Please refer to the ``The format of libraries for import to ANNOgesic``.

processing site: Some 3'UTR derived and interCDS derived sRNA candidates start 
from processing site not TSS.

- Arguments

::

    usage: ANNOgesic.py srna [-h] [--Vienna_folder VIENNA_FOLDER]
                             [--Vienna_utils VIENNA_UTILS]
                             [--blast_plus_folder BLAST_PLUS_FOLDER]
                             [--ps2pdf14_path PS2PDF14_PATH] [--UTR_derived_sRNA]
                             [--import_info IMPORT_INFO [IMPORT_INFO ...]]
                             [--transcript_assembly_folder TRANSCRIPT_ASSEMBLY_FOLDER]
                             [--annotation_folder ANNOTATION_FOLDER]
                             [--TSS_folder TSS_FOLDER]
                             [--processing_site_folder PROCESSING_SITE_FOLDER]
                             [--TSS_intergenic_fuzzy TSS_INTERGENIC_FUZZY]
                             [--TSS_5UTR_fuzzy TSS_5UTR_FUZZY]
                             [--TSS_3UTR_fuzzy TSS_3UTR_FUZZY]
                             [--TSS_interCDS_fuzzy TSS_INTERCDS_FUZZY]
                             [--min_length MIN_LENGTH] [--max_length MAX_LENGTH]
                             [--tex_wig_folder TEX_WIG_FOLDER]
                             [--frag_wig_folder FRAG_WIG_FOLDER]
                             [--cutoff_intergenic_coverage CUTOFF_INTERGENIC_COVERAGE]
                             [--cutoff_5utr_coverage CUTOFF_5UTR_COVERAGE]
                             [--cutoff_3utr_coverage CUTOFF_3UTR_COVERAGE]
                             [--cutoff_interCDS_coverage CUTOFF_INTERCDS_COVERAGE]
                             [--fasta_folder FASTA_FOLDER]
                             [--cutoff_energy CUTOFF_ENERGY] [--mountain_plot]
                             [--database_format]
                             [--sRNA_database_path SRNA_DATABASE_PATH]
                             [--nr_database_path NR_DATABASE_PATH]
                             [--sRNA_blast_stat]
                             [--tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]]
                             [--frag_libs FRAG_LIBS [FRAG_LIBS ...]]
                             [--tex_notex TEX_NOTEX]
                             [--replicates_tex REPLICATES_TEX]
                             [--replicates_frag REPLICATES_FRAG] [--table_best]
                             [--decrease_intergenic DECREASE_INTERGENIC]
                             [--decrease_utr DECREASE_UTR]
                             [--fuzzy_intergenic FUZZY_INTERGENIC]
                             [--fuzzy_utr FUZZY_UTR]
                             [--cutoff_nr_hit CUTOFF_NR_HIT]
                             [--blast_e_nr BLAST_E_NR]
                             [--blast_e_srna BLAST_E_SRNA] [--sORF SORF]
                             [--best_with_all_sRNAhit]
                             [--best_without_sORF_candidate]
                             [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given the
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
                            If you want to detect UTR derived sRNA, please turn it
                            on. Default is False.
      --import_info IMPORT_INFO [IMPORT_INFO ...], -d IMPORT_INFO [IMPORT_INFO ...]
                            There are several types of information you can import
                            to detect and filter sRNA: TSS(1), energy of secondary
                            structure(2), blast to nr(3), blast to sRNA(4),
                            sORF(5), without any information, only detect sRNA
                            (without any information) by transcriptome
                            assembly(6).Please assign the number of type you want
                            to import, i.e. 1 2 4 - means it used TSS, energy and
                            blast result to detect sRNA. Besides these
                            information, it will also consider the sequence length
                            of sRNA.
      --transcript_assembly_folder TRANSCRIPT_ASSEMBLY_FOLDER, -a TRANSCRIPT_ASSEMBLY_FOLDER
                            The path of transcriptome assembly folder.
      --annotation_folder ANNOTATION_FOLDER, -g ANNOTATION_FOLDER
                            The path of annotation gff folder.
      --TSS_folder TSS_FOLDER, -t TSS_FOLDER
                            If you want to import TSS information, please assign
                            the path of gff folder of TSS. If you want to detect
                            UTR derived sRNA, you must assign the folder of TSS.
      --processing_site_folder PROCESSING_SITE_FOLDER, -p PROCESSING_SITE_FOLDER
                            If you want to import processing site information,
                            please assign the path of gff folder of processing
                            site.If you want to detect UTR derived sRNA, you must
                            assign the folder of processing site.
      --TSS_intergenic_fuzzy TSS_INTERGENIC_FUZZY, -ft TSS_INTERGENIC_FUZZY
                            If you want to import TSS information, you need to
                            assign the fuzzy for comparing TSS and transcript
                            assembly. It is for intergenic.The default number is
                            2.
      --TSS_5UTR_fuzzy TSS_5UTR_FUZZY, -f5 TSS_5UTR_FUZZY
                            If you want to import TSS information, you need to
                            assign the fuzzy for comparing TSS and transcript
                            assembly. It is for 5'UTR of UTR derived sRNA.The
                            default number is 2.
      --TSS_3UTR_fuzzy TSS_3UTR_FUZZY, -f3 TSS_3UTR_FUZZY
                            If you want to import TSS information, you need to
                            assign the fuzzy for comparing TSS and transcript
                            assembly. It is for 3'UTR of UTR derived sRNA.The
                            default number is 10.
      --TSS_interCDS_fuzzy TSS_INTERCDS_FUZZY, -fc TSS_INTERCDS_FUZZY
                            If you want to import TSS information, you need to
                            assign the fuzzy for comparing TSS and transcript
                            assembly. It is for interCDS derived sRNA.The default
                            number is 10.
      --min_length MIN_LENGTH, -lm MIN_LENGTH
                            Please assign the minium length of sRNA. It will
                            classify sRNA candidates based on the value. Default
                            is 30.
      --max_length MAX_LENGTH, -lM MAX_LENGTH
                            Please assign the maxium length of sRNA. It will
                            classify sRNA candidates based on the value. Default
                            is 500.
      --tex_wig_folder TEX_WIG_FOLDER, -tw TEX_WIG_FOLDER
                            The path of tex+/- wig folder.
      --frag_wig_folder FRAG_WIG_FOLDER, -fw FRAG_WIG_FOLDER
                            The path of fragment wig folder.
      --cutoff_intergenic_coverage CUTOFF_INTERGENIC_COVERAGE, -ci CUTOFF_INTERGENIC_COVERAGE
                            The cutoff of minimal coverage of intergenic sRNA
                            candidates. Default is 5
      --cutoff_5utr_coverage CUTOFF_5UTR_COVERAGE, -cu5 CUTOFF_5UTR_COVERAGE
                            The cutoff of minimal coverage of 5'UTR derived sRNA
                            candidates. You can also assign median or mean.
                            Default is median
      --cutoff_3utr_coverage CUTOFF_3UTR_COVERAGE, -cu3 CUTOFF_3UTR_COVERAGE
                            The cutoff of minimal coverage of 3'UTR derived sRNA
                            candidates. You can also assign median or mean.
                            Default is median
      --cutoff_interCDS_coverage CUTOFF_INTERCDS_COVERAGE, -cuf CUTOFF_INTERCDS_COVERAGE
                            The cutoff of minimal coverage of inter CDS sRNA
                            candidates. You can also assign median or mean.
                            Default is median
      --fasta_folder FASTA_FOLDER, -f FASTA_FOLDER
                            If you want to import secondary structure information,
                            please assign the path of fasta folder.
      --cutoff_energy CUTOFF_ENERGY, -e CUTOFF_ENERGY
                            If you want to import secondary structure information,
                            please assign the cutoff of folding energy. It will
                            classify sRNA candidates based on the value. Default
                            is 0.
      --mountain_plot, -m   If you want to generate mountain plots of sRNA
                            candidates, please turn it on. Default is False.
      --database_format, -fd
                            If you already format your database, you don't need to
                            turn it on. Default is False
      --sRNA_database_path SRNA_DATABASE_PATH, -sd SRNA_DATABASE_PATH
                            If you want to import blast results of sRNA, please
                            assign the path of sRNA database.
      --nr_database_path NR_DATABASE_PATH, -nd NR_DATABASE_PATH
                            If you want to import blast results of nr, please
                            assign the path of nr database.
      --sRNA_blast_stat, -sb
                            If the sRNA database which you used are the same
                            format as our default sRNA database, you can run
                            sRNA_blast_stat for do statistics of the result of
                            sRNA blast.If your format is not the same as our
                            default database, please don't turn it on. Out default
                            format of header is ID|strain|srna_name
      --tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...], -tl TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]
                            library name of tex and notex library. The format is:
                            wig_file_name:tex_treat_or_not(tex or notex):condition
                            _id(integer):replicate_id(alphabet):strand(+ or -).
      --frag_libs FRAG_LIBS [FRAG_LIBS ...], -fl FRAG_LIBS [FRAG_LIBS ...]
                            library name of fragment library. The format is: wig_f
                            ile_name:fragmented(frag):condition_id(integer):replic
                            ate_id(alphabet):strand(+ or -).
      --tex_notex TEX_NOTEX, -te TEX_NOTEX
                            For tex +/- library, sRNA candidates should be
                            detected by both or just one.(1/2) Default is 2.
      --replicates_tex REPLICATES_TEX, -rt REPLICATES_TEX
                            The sRNA of tex +/- library should be detected more
                            than this number of replicates.
      --replicates_frag REPLICATES_FRAG, -rf REPLICATES_FRAG
                            The sRNA of fragmented library should be detected more
                            than this number of replicates.
      --table_best, -tb     The output table of sRNA candidates only print the
                            best track. Default is False
      --decrease_intergenic DECREASE_INTERGENIC, -di DECREASE_INTERGENIC
                            If the intergenic region is longer than the
                            max_length, it will based on coverage to check the
                            sRNA candidates. If the ratio of lowest coverage of
                            intergenic region and the highest coverage of
                            intergenic region is smaller than this number, it will
                            consider the the point of lowest coverage to be end of
                            sRNA. If the length of sRNA candidate is properly, it
                            also assign the transcript to be one of sRNA
                            candidates. Default is 0.5.
      --decrease_utr DECREASE_UTR, -du DECREASE_UTR
                            If the kind of utr derived is 5'UTR, you have to
                            consider the end of it's end.If the ratio of lowest
                            coverage of it and the highest coverage of it is
                            smaller than this number, it will consider the the
                            point of lowest coverage to be end of sRNA. If the
                            length of sRNA candidate is properly, it also assign
                            the transcript to be one of sRNA candidates. Default
                            is 0.5.
      --fuzzy_intergenic FUZZY_INTERGENIC, -fi FUZZY_INTERGENIC
                            If the situation is like decrease_intergenic
                            mentioned, the value would be fuzzy between the end of
                            sRNA
      --fuzzy_utr FUZZY_UTR, -fu FUZZY_UTR
                            If the situation is like decrease_utr mentioned, the
                            value would be fuzzy between the end of sRNA
      --cutoff_nr_hit CUTOFF_NR_HIT, -cn CUTOFF_NR_HIT
                            The cutoff of number of hits in nr database. If the
                            number of nr hits more than this cutoff, program will
                            exclude it during classification.
      --blast_e_nr BLAST_E_NR, -en BLAST_E_NR
                            The cutoff of blast e value for nr alignment.
      --blast_e_srna BLAST_E_SRNA, -es BLAST_E_SRNA
                            The cutoff of blast e value for sRNA alignment.
      --sORF SORF, -O SORF  If you want to compare sORF and sRNA, please assign
                            the path of sORF gff folder.
      --best_with_all_sRNAhit, -ba
                            When you want to generate the files which store the
                            best sRNA candidates, it should include all the sRNA
                            candidates which can find the homology from blast sRNA
                            database no matter other information(ex. TSS, blast in
                            nr...).Please turn it on. Or it will just select the
                            the best candidates based on all filter conditions.
                            Default is False.
      --best_without_sORF_candidate, -bs
                            If you want to generate the files which store the best
                            sRNA candidates excluded all the sRNA candidates which
                            also can be detected by sORF file.Please turn it on.
                            Or it will just select the the best candidates based
                            on all filter conditions. Default is False.

- Output files

All output files will be stored in ``$ANNOgesic/output/sRNA``

``sRNA_2d_$STRAIN_NAME``: The secondary structure of all sRNA candidates.

``sRNA_seq_$STRAIN_NAME``: The sequence of all sRNA candidates.

``blast_result_and_misc``: the results of blast.

``mountain_plot``: the mountain plot of sRNA candidates.

``sec_structure``: the dot plot and secondary structure plot of sRNA candidates.

``statistics``: statistics files. ``stat_sRNA_blast_class_$STRAIN_NAME.csv`` is the results of analysis of blast sRNA database.
``stat_sRNA_class_Staphylococcus_aureus_HG003.csv`` is the results of classification of sRNA candidates.

``tables``: sRNA tables for more details. ``for class`` is for different classes of sRNAs.
``best`` is the best results of sRNA. ``all_candidates`` is for all candidates without filtering.

``gffs``: gff files of sRNA. The meaning of ``for class``, ``best``, ``all_candidates`` is the same as ``tables``.

The tags of gff file:

Second column presents the type of sRNA - ``intergenic`` or ``UTR_derived``

``UTR_type``: The sRNA is from 5'UTR or 3'UTR or interCDS.

``with_TSS``: The sRNA is from which TSS. NA means the sRNA is not related with TSS.

``with_cleavage``: The sRNA is ended at which processing site.

``sORF``: Which sORF overlap with this sRNA.

``sRNA_hit``: The blast hit of sRNA database.

``nr_hit``: The blast hit of nr database.

``2d_energy``: The folding energy of sRNA candidate.

sorf
----------
``sorf`` can detect start codon and stop codon within the intergenic region.
User can also import some information to filter false positive. Because non-coding region 
may be sRNA or sORF, it also provide the function to compare sORF and sRNA. If start and stop 
codons are more than one in sORF region. ``sorf`` will provide the longest one and all the start 
and stop codons. User can refer to it. BE CAREFUL, The length between start codon and stop codon 
should be multiple of 3 or it will not be sORF. The position of start codon is the first nucleotide.
The position of stop codon is the last nucleotide.

- Pre-required tools and files

The gff files of CDS/tRNA/rRNA and transcripts.

The libraries and wiggle files, Please refer to the ``The format of libraries for import to ANNOgesic``.

The fasta for detect start codon and stop codon.

User can also import some useful information to improve the prediction:

gff file of TSS for checking the sORF start from TSS or not. 
gff file of sRNA for checking the conflict of sRNA and sORF.

- Arguments

::

    usage: ANNOgesic.py sorf [-h] [--UTR_derived_sORF]
                             [--transcript_assembly_folder TRANSCRIPT_ASSEMBLY_FOLDER]
                             [--annotation_folder ANNOTATION_FOLDER]
                             [--TSS_folder TSS_FOLDER] [--utr_length UTR_LENGTH]
                             [--min_length MIN_LENGTH] [--max_length MAX_LENGTH]
                             [--tex_wig_folder TEX_WIG_FOLDER]
                             [--frag_wig_folder FRAG_WIG_FOLDER]
                             [--cutoff_intergenic_coverage CUTOFF_INTERGENIC_COVERAGE]
                             [--cutoff_5utr_coverage CUTOFF_5UTR_COVERAGE]
                             [--cutoff_3utr_coverage CUTOFF_3UTR_COVERAGE]
                             [--cutoff_interCDS_coverage CUTOFF_INTERCDS_COVERAGE]
                             [--cutoff_background CUTOFF_BACKGROUND]
                             [--fasta_folder FASTA_FOLDER]
                             [--tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]]
                             [--frag_libs FRAG_LIBS [FRAG_LIBS ...]]
                             [--tex_notex TEX_NOTEX]
                             [--replicates_tex REPLICATES_TEX]
                             [--replicates_frag REPLICATES_FRAG] [--table_best]
                             [--sRNA_folder SRNA_FOLDER]
                             [--start_coden START_CODEN [START_CODEN ...]]
                             [--stop_coden STOP_CODEN [STOP_CODEN ...]]
                             [--condition_best CONDITION_BEST]
                             [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --UTR_derived_sORF, -u
                            If you want to detect UTR derived sORF, please turn it
                            on. Default is False.
      --transcript_assembly_folder TRANSCRIPT_ASSEMBLY_FOLDER, -a TRANSCRIPT_ASSEMBLY_FOLDER
                            The path of transcriptome assembly folder.
      --annotation_folder ANNOTATION_FOLDER, -g ANNOTATION_FOLDER
                            The path of annotation gff folder.
      --TSS_folder TSS_FOLDER, -t TSS_FOLDER
                            If you want to import TSS information, please assign
                            the path of gff folder of TSS.
      --utr_length UTR_LENGTH, -ul UTR_LENGTH
                            If you want to import TSS information, please assign
                            the utr length for comparing TSS and sORF. The default
                            number is 300.
      --min_length MIN_LENGTH, -lm MIN_LENGTH
                            Please assign the minium length of sORF. It will
                            classify sORF candidates based on the value. Default
                            is 30.
      --max_length MAX_LENGTH, -lM MAX_LENGTH
                            Please assign the maxium length of sORF. It will
                            classify sORF candidates based on the value. Default
                            is 500.
      --tex_wig_folder TEX_WIG_FOLDER, -tw TEX_WIG_FOLDER
                            The path of tex+/- wig folder.
      --frag_wig_folder FRAG_WIG_FOLDER, -fw FRAG_WIG_FOLDER
                            The path of fragment wig folder.
      --cutoff_intergenic_coverage CUTOFF_INTERGENIC_COVERAGE, -ci CUTOFF_INTERGENIC_COVERAGE
                            The cutoff of minimal coverage of intergenic sORF
                            candidates.
      --cutoff_5utr_coverage CUTOFF_5UTR_COVERAGE, -cu5 CUTOFF_5UTR_COVERAGE
                            The cutoff of minimal coverage of 5'UTR derived sORF
                            candidates. You also can assign median or mean.
                            Default is median.
      --cutoff_3utr_coverage CUTOFF_3UTR_COVERAGE, -cu3 CUTOFF_3UTR_COVERAGE
                            The cutoff of minimal coverage of 3'UTR derived sORF
                            candidates. You also can assign median or mean.
                            Default is median.
      --cutoff_interCDS_coverage CUTOFF_INTERCDS_COVERAGE, -cuf CUTOFF_INTERCDS_COVERAGE
                            The cutoff of minimal coverage of interCDS derived
                            sORF candidates. You also can assign median or mean.
                            Default is median.
      --cutoff_background CUTOFF_BACKGROUND, -cub CUTOFF_BACKGROUND
                            The cutoff of minimal coverage of all sORF candidates.
                            Default is 5.
      --fasta_folder FASTA_FOLDER, -f FASTA_FOLDER
                            The folder of fasta file.
      --tex_notex_libs TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...], -tl TEX_NOTEX_LIBS [TEX_NOTEX_LIBS ...]
                            Library name of tex and notex library. The format is:
                            wig_file_name:tex_treat_or_not(tex or notex):condition
                            _id(integer):replicate_id(alphabet):strand(+ or -).
      --frag_libs FRAG_LIBS [FRAG_LIBS ...], -fl FRAG_LIBS [FRAG_LIBS ...]
                            Library name of fragment library The format is: wig_fi
                            le_name:fragmented(frag):condition_id(integer):replica
                            te_id(alphabet):strand(+ or -)..
      --tex_notex TEX_NOTEX, -te TEX_NOTEX
                            For tex +/- library, sORF candidates should be
                            detected by both or just one.(1/2) Default is 2.
      --replicates_tex REPLICATES_TEX, -rt REPLICATES_TEX
                            The sORF of tex +/- library should be detected more
                            than this number of replicates.
      --replicates_frag REPLICATES_FRAG, -rf REPLICATES_FRAG
                            The sORF of fragmented library should be detected more
                            than this number of replicates.
      --table_best, -tb     The output table of sORF candidates only print the
                            best track. Default is False.
      --sRNA_folder SRNA_FOLDER, -s SRNA_FOLDER
                            If you want to compare sORF and sRNA, please assign
                            the path of sORF gff folder.
      --start_coden START_CODEN [START_CODEN ...], -ac START_CODEN [START_CODEN ...]
                            What kinds of start coden ATG/GTG/TTG you want to use.
      --stop_coden STOP_CODEN [STOP_CODEN ...], -oc STOP_CODEN [STOP_CODEN ...]
                            What kinds of stop coden TTA/TAG/TGA you want to use.
      --condition_best CONDITION_BEST, -c CONDITION_BEST
                            For generating the result of best sORF, please assign
                            which information you want to consider
                            (TSS/sRNA/both). default is TSS.

- Output files

All output files will be stored in ``$ANNOgesic/output/sORF``

``statistics``: statistics files.

``tables``: the tables of sORF for more details. ``all_candidates`` is for all sORF candidates without filtering.
``best`` is for the best sORF candidates with filtering. The table stores all the TSS, start codon and 
stop codon information.

``gffs``: gff files of sORF. The meanings of ``all_candidates`` and ``best`` are the same as ``tables``.

The tags of gff file:

``start_TSS`` is the sORF start from with TSS.

``with_TSS`` means all TSSs which are located in the region of this sORF.

``UTR_type`` is the type of UTR if the sORF is UTR derived one.

``sRNA`` means which sRNA overlaps with this sORF.

promoter
-----------

``promoter`` can scan the upstream of TSS to discover the promoter motif.
User can assign the region of upstream TSS. We integrate MEME to compute the promoter.
User can view the result very easily. If the gff file of TSS is not computed by 
ANNOgesic, please use --TSS_source. ``promoter`` will classify the TSS for computing 
promoter motif.

- Pre-required tools and files

`MEME <http://meme-suite.org/tools/meme>`_

Fasta and gff file of genome.

Gff file of TSS.

If the gff file of TSS is not computed by ANNOgesic, the libraries and wiggle files are necessary.
Please refer to the ``The format of libraries for import to ANNOgesic``.

- Arguments

::

    usage: ANNOgesic.py promoter [-h] [--MEME_path MEME_PATH]
                                 [--fasta_folder FASTA_FOLDER]
                                 [--TSS_folder TSS_FOLDER] [--num_motif NUM_MOTIF]
                                 [--motif_width MOTIF_WIDTH [MOTIF_WIDTH ...]]
                                 [--parallel PARALLEL] [--TSS_source]
                                 [--tex_libs TEX_LIBS [TEX_LIBS ...]]
                                 [--tex_wig_path TEX_WIG_PATH]
                                 [--annotation_folder ANNOTATION_FOLDER]
                                 [--combine_all]
                                 [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --MEME_path MEME_PATH
                            path of MEME
      --fasta_folder FASTA_FOLDER, -f FASTA_FOLDER
                            Please assign the folder of gemonic fasta file.
      --TSS_folder TSS_FOLDER, -t TSS_FOLDER
                            The folder of TSS gff file.
      --num_motif NUM_MOTIF, -n NUM_MOTIF
                            How many of motifs you want to produce.
      --motif_width MOTIF_WIDTH [MOTIF_WIDTH ...], -w MOTIF_WIDTH [MOTIF_WIDTH ...]
                            Motif length - it will refer the value to find the
                            motif. if you want to detect a range of width, you can
                            insert "-" between two values. for example, 50 2-10.
                            It means the range of width which you want to detect
                            is 50 and within 2 to 10.
      --parallel PARALLEL, -p PARALLEL
                            How many process you want to use to run paralle.
      --TSS_source, -s      If you generate TSS from other method, please turn it
                            on.
      --tex_libs TEX_LIBS [TEX_LIBS ...], -tl TEX_LIBS [TEX_LIBS ...]
                            Library name of tex+/- library. If your TSS is not
                            from ANNOgesic, please assign the libs of tex+/-
                            too.The format is: wig_file_name:tex_treat_or_not(tex
                            or notex):condition_id(integer):replicate_id(alphabet)
                            :strand(+ or -).
      --tex_wig_path TEX_WIG_PATH, -tw TEX_WIG_PATH
                            The path of tex+/- wig folder. If your TSS is not from
                            ANNOgesic, please assign the wig path too.
      --annotation_folder ANNOTATION_FOLDER, -g ANNOTATION_FOLDER
                            The path of annotation gff folder. If your TSS is not
                            from ANNOgesic, please assign the annotation gff path
                            too.
      --combine_all, -c     If you want to combine all TSS in TSS output folder to
                            generate a overall promoter motif, please turn it on.
                            Default is False.

- Output files

All output file will be stored in ``$ANNOgesic/output/promoter_analysis``.

``allfasta``: the promoter information of all TSS (including all strains in gff file)

Every strain will generate one folder for storing the information of promoter motif.

If the TSS is not computed by ANNOgesic, ``TSS_class`` will be generated. It will classify the 
TSS and store as gff file in it.

operon
----------

``operon`` will group TSS, gene/CDS/tRNA/rRNA, transcript, terminator and UTR to operon and 
suboperon.

- Pre-required tools or files

Gff files of TSS, annotation, transcript, 5'UTR, and 3'UTR.

If user want to import the information of terminator, ``operon`` can integrate terminator, too.

- Arguments

::

    usage: ANNOgesic.py operon [-h] [--TSS_folder TSS_FOLDER]
                               [--annotation_folder ANNOTATION_FOLDER]
                               [--transcript_folder TRANSCRIPT_FOLDER]
                               [--UTR5_folder UTR5_FOLDER]
                               [--UTR3_folder UTR3_FOLDER]
                               [--term_folder TERM_FOLDER] [--TSS_fuzzy TSS_FUZZY]
                               [--term_fuzzy TERM_FUZZY] [--min_length MIN_LENGTH]
                               [--statistics] [--combine_gff]
                               [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --TSS_folder TSS_FOLDER, -t TSS_FOLDER
                            The path of TSS gff folder.
      --annotation_folder ANNOTATION_FOLDER, -g ANNOTATION_FOLDER
                            The path of annotation gff folder.
      --transcript_folder TRANSCRIPT_FOLDER, -a TRANSCRIPT_FOLDER
                            The path of transcript assembly gff folder.
      --UTR5_folder UTR5_FOLDER, -u5 UTR5_FOLDER
                            The path of 5'UTR gff folder.
      --UTR3_folder UTR3_FOLDER, -u3 UTR3_FOLDER
                            The path of 3'UTR gff folder.
      --term_folder TERM_FOLDER, -e TERM_FOLDER
                            If you want to import the information of terminator,
                            please assign the path of terminator gff folder.
      --TSS_fuzzy TSS_FUZZY, -tf TSS_FUZZY
                            The fuzzy for comparing TSS and transcript assembly.
                            The default number is 5.
      --term_fuzzy TERM_FUZZY, -ef TERM_FUZZY
                            The fuzzy for comparing terminator and transcript
                            assembly. The default number is 30.
      --min_length MIN_LENGTH, -l MIN_LENGTH
                            The minimum length of operon. The default number is
                            20.
      --statistics, -s      Doing statistics for Operon analysis. Default is
                            False.The name of statistics file is -
                            stat_operon_$STRAIN_NAME.csv.
      --combine_gff, -c     Convert the operon and all features you assigned to
                            one gff file. Default is False.

- Output files

All output files will be stored in ``$ANNOgesic/output/operon``

``gffs``: the gff files which integrate the information of TSS, annotation, 
transcript, 5'UTR, and 3'UTR. The order of rows in gff file is according to the operon.

``tables``: the table of operon which stores all information of operon and suboperon.

``statistics``: the statistics files.

circrna
--------------

``circrna`` can detect the potential circular RNA. It uses `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_ 
to detect circular RNA. Then check annotation file and quality of splicing site detection to 
exclude false positive. User can assign reads for mapping and detect circular RNA or assign alignment files to skip mapping.
But BE CAREFUL, If user uses alignment files, they should be mapped by `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_ 
and with ``-S``. Or ``circrna`` can't find the proper candidates.

- Pre-required tools and files

`Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_

Fasta files of reads or alignment files. If you want input alignment files directly, remember they should be 
mapped by `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_ and with ``-S``.

Fasta file and gff file of genome

- Arguments

::

    usage: ANNOgesic.py circrna [-h] [--segemehl_folder SEGEMEHL_FOLDER]
                                [--samtools_path SAMTOOLS_PATH] [--align]
                                [--normal_bam_path NORMAL_BAM_PATH]
                                [--fragmented_bam_path FRAGMENTED_BAM_PATH]
                                [--process PROCESS] [--fasta_path FASTA_PATH]
                                [--annotation_path ANNOTATION_PATH]
                                [--convert_to_gff] [--support_reads SUPPORT_READS]
                                [--start_ratio START_RATIO]
                                [--end_ratio END_RATIO]
                                [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --segemehl_folder SEGEMEHL_FOLDER, -sg SEGEMEHL_FOLDER
                            Please assign the folder of segemehl.
      --samtools_path SAMTOOLS_PATH, -st SAMTOOLS_PATH
                            Please assign the path of samtools.
      --align, -a           Using segemehl to align read (included splice
                            detection). If you already usd segemehl with -S to
                            align your reads, you can skip this step, don't need
                            to turn it on. Please be attention, it only use
                            default parameters of segemehl to align your reads.
                            Moreover, it will align all read files in
                            ANNOgesic/input/reads. If you want to run some
                            specific functions of segemehl, please run it
                            seperately.
      --normal_bam_path NORMAL_BAM_PATH, -nb NORMAL_BAM_PATH
                            If you already has Bam files, Please assign the normal
                            Bam path or fragmented_Bam_path.
      --fragmented_bam_path FRAGMENTED_BAM_PATH, -fb FRAGMENTED_BAM_PATH
                            If you already has Bam files, Please assign the
                            fragmented Bam path or normal Bam path.
      --process PROCESS, -p PROCESS
                            How many parallels processes for --align.
      --fasta_path FASTA_PATH, -f FASTA_PATH
                            The folder of genome fasta.
      --annotation_path ANNOTATION_PATH, -g ANNOTATION_PATH
                            The folder of annotation gff files.
      --convert_to_gff, -cg
                            If you want to convert circRNA candidates to gff file,
                            please turnn it on.
      --support_reads SUPPORT_READS, -s SUPPORT_READS
                            If you want to convert circRNA candidates to gff file,
                            please also assign the cut off of supported reads.
                            Default is 5.
      --start_ratio START_RATIO, -sr START_RATIO
                            The ratio of (read support circ / all read) at
                            starting point. The ratio of candidates should higher
                            than this cutoff. Default is 0.25.
      --end_ratio END_RATIO, -er END_RATIO
                            The ratio of (read support circ / all read) at end
                            point. The ratio of candidates should higher than this
                            cutoff. Default is 0.25.

- Output files

All the output files will be stored in ``$ANNOgesic/output/circRNA``

``gffs``: gff files of circular RNA. ``$STRAINNAME_best.gff`` is the gff file for best result after comparing 
with annotation and quality of splicing. ``$STRAINNAME_all.gff`` is for all candidates without filering.

``circRNA_tables``: the tables for circular RNA with more details.

``statistics``: statistics files

``segemehl_align``: if ``circrna`` starts from read mapping, the folder is for results of mapping.

``segemehl_splice``: the results of splicing detection. The information of the splicing table, please 
refer to `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_.

go_term
----------

``go_term`` can compare the annotation file and the Uniprot to retreive the information of Gene Ontology.
It also provide some analysis of the classes of Go term.

- Pre-required tools and files

`idmapping_selected.tab from Uniprot <http://www.uniprot.org/downloads>`_

`goslim.obo <http://geneontology.org/page/go-slim-and-subset-guide>`_

`go.obo <http://geneontology.org/page/download-ontology>`_

Gff file of annotation of genome.

- Arguments

::

    usage: ANNOgesic.py go_term [-h] [--annotation_path ANNOTATION_PATH]
                                [--UniProt_id UNIPROT_ID] [--go_obo GO_OBO]
                                [--goslim_obo GOSLIM_OBO]
                                [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --annotation_path ANNOTATION_PATH, -g ANNOTATION_PATH
                            The path of annotation gff folder.
      --UniProt_id UNIPROT_ID, -u UNIPROT_ID
                            The path of UniProt ID mapping database. default is
                            ANNOgesic/input/database/idmapping_selected.tab.
      --go_obo GO_OBO, -go GO_OBO
                            The path of go.obo. Default is
                            ANNOgesic/input/database/go.obo.
      --goslim_obo GOSLIM_OBO, -gs GOSLIM_OBO
                            The path of goslim.obo. Default is
                            ANNOgesic/input/database/goslim.obo.

- Output files

All output files will be stored in ``ANNOgesic/output/Go_term``

``Go_term_results``: the table of Go term information.

``statistics``: statistics files and visualization files.

srna_target
---------------

``srna_target`` will search the possible target of sRNA. Users can assign which 
program they want to use (RNAup or RNAplex or both). We recommand use both of the 
programs. ``srna_target`` can compare the both results and provide the best ones.

- Pre-required tools and files

`ViennaRNA <http://www.tbi.univie.ac.at/RNA/>`_ 

Gff files of annotation of genome and sRNA

Fasta file of genome

- Arguments

::

    usage: ANNOgesic.py srna_target [-h] [--Vienna_folder VIENNA_FOLDER]
                                    [--annotation_path ANNOTATION_PATH]
                                    [--fasta_path FASTA_PATH]
                                    [--sRNA_path SRNA_PATH]
                                    [--query_sRNA QUERY_SRNA [QUERY_SRNA ...]]
                                    [--program PROGRAM]
                                    [--interaction_length INTERACTION_LENGTH]
                                    [--window_size_target WINDOW_SIZE_TARGET]
                                    [--span_target SPAN_TARGET]
                                    [--window_size_srna WINDOW_SIZE_SRNA]
                                    [--span_srna SPAN_SRNA]
                                    [--unstructured_region_RNAplex_target UNSTRUCTURED_REGION_RNAPLEX_TARGET]
                                    [--unstructured_region_RNAplex_srna UNSTRUCTURED_REGION_RNAPLEX_SRNA]
                                    [--unstructured_region_RNAup UNSTRUCTURED_REGION_RNAUP]
                                    [--energy_threshold ENERGY_THRESHOLD]
                                    [--duplex_distance DUPLEX_DISTANCE]
                                    [--top TOP]
                                    [--process_rnaplex PROCESS_RNAPLEX]
                                    [--process_rnaup PROCESS_RNAUP]
                                    [--continue_rnaup]
                                    [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --Vienna_folder VIENNA_FOLDER
                            Please assign the folder of Vienna package. It should
                            include RNAplfold, RNAup and RNAplex.
      --annotation_path ANNOTATION_PATH, -g ANNOTATION_PATH
                            The path of annotation gff folder.
      --fasta_path FASTA_PATH, -f FASTA_PATH
                            The path of genome fasta folder.
      --sRNA_path SRNA_PATH, -r SRNA_PATH
                            The path of sRNA gff folder.
      --query_sRNA QUERY_SRNA [QUERY_SRNA ...], -q QUERY_SRNA [QUERY_SRNA ...]
                            Please assign the query sRNA. If you want to compute
                            all sRNA in gff file, please keyin 'all'.the input
                            format should be like, $STRAIN:$STRAND:$START:$END.For
                            example, NC_007795:+:200:534 NC_007795:-:6767:6900
      --program PROGRAM, -p PROGRAM
                            Using RNAplex, RNAup or both. Default is both.
      --interaction_length INTERACTION_LENGTH, -i INTERACTION_LENGTH
                            Maximal length of an interaction. Default is 30.
      --window_size_target WINDOW_SIZE_TARGET, -wt WINDOW_SIZE_TARGET
                            Only work when --program is RNAplex or both. Average
                            the pair probabilities over windows of given size for
                            RNAplex. Default is 240.
      --span_target SPAN_TARGET, -st SPAN_TARGET
                            Only work when --program is RNAplex or both. Set the
                            maximum allowed separation of a base pair to span for
                            RNAplex. Default is 160.
      --window_size_srna WINDOW_SIZE_SRNA, -ws WINDOW_SIZE_SRNA
                            Only work when --program is RNAplex or both. Average
                            the pair probabilities over windows of given size for
                            RNAplex. Default is 30.
      --span_srna SPAN_SRNA, -ss SPAN_SRNA
                            Only work when --program is RNAplex or both. Set the
                            maximum allowed separation of a base pair to span for
                            RNAplex. Default is 30.
      --unstructured_region_RNAplex_target UNSTRUCTURED_REGION_RNAPLEX_TARGET, -ut UNSTRUCTURED_REGION_RNAPLEX_TARGET
                            Only work when --program is RNAplex or both. Compute
                            the mean probability that regions of length 1 to a
                            given length are unpaired for RNAplex. Default is 30.
      --unstructured_region_RNAplex_srna UNSTRUCTURED_REGION_RNAPLEX_SRNA, -us UNSTRUCTURED_REGION_RNAPLEX_SRNA
                            Only work when --program is RNAplex or both. Compute
                            the mean probability that regions of length 1 to a
                            given length are unpaired for RNAplex. Default is 30.
      --unstructured_region_RNAup UNSTRUCTURED_REGION_RNAUP, -uu UNSTRUCTURED_REGION_RNAUP
                            Only work when --program is RNAup or both. Compute the
                            mean probability that regions of length 1 to a given
                            length are unpaired for RNAplex. Default is 30.
      --energy_threshold ENERGY_THRESHOLD, -e ENERGY_THRESHOLD
                            Only work when --program is RNAplex or both. Minimal
                            energy for a duplex to be returned for RNAplex.
                            Default is -8.
      --duplex_distance DUPLEX_DISTANCE, -d DUPLEX_DISTANCE
                            Only work when --program is RNAplex or both. Distance
                            between target 3' ends of two consecutive duplexes for
                            RNAplex. Default is 20.
      --top TOP, -t TOP     The output file only include top ones(default is 20).
      --process_rnaplex PROCESS_RNAPLEX, -pp PROCESS_RNAPLEX
                            How many parallel processes for running RNAplex
                            prediction.
      --process_rnaup PROCESS_RNAUP, -pu PROCESS_RNAUP
                            How many parallel processes for running RNAup
                            prediction.
      --continue_rnaup, -cr
                            RNAup will take a long time for running if you want to
                            compute a lot of sRNA. If the process crush, you can
                            turn it on. This flag will continue running RNAup
                            based on your previous running.

- Output files

All output files will be stored in ``$ANNOgesic/output/sRNA_targets``.

``RNAplex``: the results of RNAplex. ``$STRAIN_RNAplex.txt`` is raw results from RNAplex.
It includes the information of binding situation. ``$STRAIN_RNAplex_rank.csv`` is the results 
that sort by binding energy.

``RNAup``: the results of RNAup. ``$STRAIN_RNAup.txt`` is raw results from RNAup.
It includes the information of binding situation. ``$STRAIN_RNAup_rank.csv`` is the results
that sort by binding energy.

``merge``: the results which merge ``RNAplex`` and ``RNAup``. ``$STRAIN_merge.csv`` is just 
merge the results. ``$STRAIN_overlap.csv`` is list the results which be top 20(default) in both methods.

``sRNA_seqs``: the fasta sequence of sRNA.

``target_seqs``: the fasta sequence of potential target.

ppi_network
-------------

``ppi_network`` will retrieve the data from `STRING <http://string-db.org/>`_. 
Then use `PIE <http://www.ncbi.nlm.nih.gov/CBBresearch/Wilbur/IRET/PIE/>`_ to search 
the literatures to support the Protein-protein interaction network. Therefore, 
``ppi_network`` can generate the Protein-protein interaction network with literatures.
User can refer to it and exclude false positive.

- Pre-required tools and files

`species.vXXXX.txt from STRING <http://string-db.org/newstring_cgi/show_download_page.pl?UserId=ReWbu8uLrfAN&sessionId=_FAQBbatf7RX>`_

Ptt file of genome

- Arguments

::

    usage: annogesic ppi_network [-h] [--ptt_path PTT_PATH] [--gff_path GFF_PATH]
                                 [--proteinID_strains PROTEINID_STRAINS [PROTEINID_STRAINS ...]]
                                 [--without_strain_pubmed]
                                 [--species_STRING SPECIES_STRING] [--score SCORE]
                                 [--node_size NODE_SIZE]
                                 [--query QUERY [QUERY ...]]
                                 [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --ptt_path PTT_PATH, -p PTT_PATH
                            The path of .ptt annotation folder.
      --gff_path GFF_PATH, -g GFF_PATH
                            The path of gff annotation folder.
      --proteinID_strains PROTEINID_STRAINS [PROTEINID_STRAINS ...], -s PROTEINID_STRAINS [PROTEINID_STRAINS ...]
                            This is for assigning protein Id which you want to
                            predict. In order to retrieve the data from STRING and
                            Pubmed, you also have to assign the similar reference.
                            For example, if you want to run all proteins in
                            Staphylococcus aureus HG003, you can assign the it
                            like Staphylococcus_aureus_HG003.ptt:Staphylococcus_au
                            reus_HG003:"Staphylococcus aureus
                            8325":"Staphylococcus aureus". or Staphylococcus_aureu
                            s_HG003.ptt:Staphylococcus_aureus_HG003:"93061":"Staph
                            ylococcus aureus". or Staphylococcus_aureus_HG003.ptt:
                            Staphylococcus_aureus_HG003:"Staphylococcus aureus
                            NCTC 8325":"Staphylococcus aureus".
                            (ptt_filename:header_ptt:STRING_name:Pubmed_name).
                            First one is the ptt file name. Second one is the
                            header of gff files. If first line of ptt
                            file has a comma, it will extract the string before
                            comma as the header, ex: Helicobacter pylori 26695
                            chromosome, complete genome - 1..1667867. The header
                            will be "Helicobacter pylori 26695 chromosome". If
                            there is no comma in first line, the header will be
                            the string before dash. ex: Helicobacter pylori -
                            1..1667867. The header will be "Helicobacter pylori".
                            Third one is for STRING database, and the fourth one
                            is for Pubmed. Of course, you can run the script for
                            several strains at the same time. Before running it,
                            please check the species file which located in
                            ANNOgesic/input/database .If you didn't download the
                            file, please download it. You can use taxon_id,
                            STRING_name_compact or official_name_NCBI to represent
                            STRING_name.BE CAREFUL, if the name which you assigned
                            has spaces, please put "" at two ends. For the name of
                            Pubmed, you can assign the name not so specific. If
                            you assign a specific name, it may not be able to find
                            the related literatures.
      --without_strain_pubmed, -n
                            If you want to retrieve pubmed without assign any
                            strains, please turn it on. Default is False.
      --species_STRING SPECIES_STRING, -d SPECIES_STRING
                            Please assign the path of species file of STRING.
      --score SCORE, -ps SCORE
                            Please assign the cutoff of score. The value is from
                            -1 to 1
      --node_size NODE_SIZE, -ns NODE_SIZE
                            Please size of the nodes in figure, default is 4000.
      --query QUERY [QUERY ...], -q QUERY [QUERY ...]
                            Please assign the query protein here. The format is
                            $HEADEROFPTT:$START_POINT:$END_POINT:$STRAND.For
                            example, Helicobacter pylori 26695
                            chromosome:345:456:+ Helicobacter pylori 26695
                            chromosome:2000:3211:-. If you want to compute all
                            protein, just type all.The default is all.


- Output files

All the output files will be stored in ``$ANNOgesic/output/PPI``.

``best_results``: the results which have supported literatures. ``$STRAIN_without_strain.csv`` is 
the results of searching without specific strain. ``$STRAIN_with_strain.csv`` is the results of 
searching with specific strain. For example, Staphylococcus_aureus_8325_without_strain.csv is 
search with Staphylococcus aureus; Staphylococcus_aureus_8325_without_strain.csv is search with 
Staphylococcus aureus 8325. ``without_strain`` stores all interaction information which search without 
specific strain. ``with_strain`` stores all interaction information which search with specific strain. 

``all_results``: the results of all protein-protein interaction. 
Even the text-mining(`PIE <http://www.ncbi.nlm.nih.gov/CBBresearch/Wilbur/IRET/PIE/>`_) score is too low,
They still stores in this folder.

``figures``: the thickness is represent how many literatures can be found for the interaction of two proteins. 
The solid line means there is some literatures which strongly support the interaction. The dash-dot line 
means the supported literatures are very weak. The dot line means there is no literatures support the 
interaction. The color is the best score of the literatures of the interaction. Basically, these figures are 
generated from the ``best_results``

subcellular_localization
------------------

``subcellular localization`` can predict where is the CDS located. It also provide some statistics and 
visualization files.

- Pre-required tools and files

`Psortb <http://www.psort.org/psortb/>`_

Gff file and fasta file of genome

- Arguments

::

    usage: ANNOgesic.py subcellular_localization [-h] [--Psortb_path PSORTB_PATH]
                                                 [--gff_path GFF_PATH]
                                                 [--fasta_path FASTA_PATH]
                                                 [--bacteria_type BACTERIA_TYPE]
                                                 [--difference_multi DIFFERENCE_MULTI]
                                                 [--merge_to_gff]
                                                 [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --Psortb_path PSORTB_PATH
                            If you want to assign the path of Psortb, please
                            assign here.
      --gff_path GFF_PATH, -g GFF_PATH
                            The path of annotation gff folder.
      --fasta_path FASTA_PATH, -f FASTA_PATH
                            The path of fasta folder.
      --bacteria_type BACTERIA_TYPE, -b BACTERIA_TYPE
                            Is Gram-positive or Gram-negative. Please assign
                            'positive' or 'negative'.
      --difference_multi DIFFERENCE_MULTI, -d DIFFERENCE_MULTI
                            If the protein may have multiple location, it will
                            calculte the difference of scores(psortb) between best
                            one and others. If the difference is within this
                            value, it will print it out, too. Default is 0.5. The
                            maximum value is 10.
      --merge_to_gff, -m    If you want to merge the information to annotation gff
                            file, please turn it on.

- Output files

All output files will be stored in ``$ANNOgesic/output/subcellular_localization``

``psortb_results``: the results of Psortb.

``statistics``: statistics files and visualization files

riboswitch
--------------

``riboswitch`` will search the intergenic region which has ribosome binding site. 
``riboswitch`` use `Infernal <http://infernal.janelia.org/>`_ to scan 
`Rfam <http://rfam.xfam.org/>`_ to find the homologs of known riboswitchs.

- Pre-required tools and files

`Infernal <http://infernal.janelia.org/>`_

`Rfam <http://rfam.xfam.org/>`_

Gff files and fasta files of genome

File of ``riboswitch_ID``. The file should contain Accession of Rfam, ID and Description of riboswitch.
The format is ``$ACCESSION{tab}$ID{tab}$DESCRIPTION``. You can download the file from our 
`Github <https://github.com/Sung-Huan/ANNOgesic>`_ (Rfam_riboswitch_ID.csv). Ours is based on the 
section of Rfam and literatures. You also can create your own one.

- Arguments

::

    usage: ANNOgesic.py riboswitch [-h] [--infernal_path INFERNAL_PATH]
                                   [--riboswitch_ID RIBOSWITCH_ID]
                                   [--gff_path GFF_PATH] [--fasta_path FASTA_PATH]
                                   [--Rfam RFAM] [--re_scan] [--e_value E_VALUE]
                                   [--output_all] [--fuzzy FUZZY]
                                   [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --infernal_path INFERNAL_PATH, -if INFERNAL_PATH
                            Please assign the folder of Infernal(where is cmscan
                            and cmsearch located).
      --riboswitch_ID RIBOSWITCH_ID, -i RIBOSWITCH_ID
                            The path of the riboswitch ID of Rfam.
      --gff_path GFF_PATH, -g GFF_PATH
                            The path of annotation gff folder.
      --fasta_path FASTA_PATH, -f FASTA_PATH
                            The path of fasta folder.
      --Rfam RFAM, -R RFAM  The path of Rfam CM database.
      --re_scan, -r         Based on the results of first scaning, it will modify
                            the input of sequence and re scan again. Default is
                            False.
      --e_value E_VALUE, -e E_VALUE
                            The cutoff of e value. Default si 0.001.
      --output_all, -o      One sequence may fit multiple riboswitches. If you
                            want to output all of them, please turn it on. Or it
                            will only print the best one.
      --fuzzy FUZZY, -z FUZZY
                            It will extend some nts of 3' ans 5' end. Default is
                            10.


- Output files

All output files will be stored in ``$ANNOgesic/output/riboswitch``

``gffs``: the gff file of riboswitch.

``tables``: the table of riboswich with more detail information.

``scan_Rfam``: the results of ``cmscan`` of `Infernal <http://infernal.janelia.org/>`_. 
User can find the results of blast.

``statistics``: the statistics files and visualization files.

optimize_tsspredator
---------------

``optimize_tsspredator`` can adapt the parameter set for the input genome. User need to 
manual check more than 100kb of genome sequence and the detected TSS need to be more than 50.
Then ``optimize_tsspredator`` can based on the principle of manual detection to scan whole genome 
in order to get the best results.

- Pre-required tools and files

`TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_.

The libraries and wiggle files of TEX +/-. 
Please refer to the ``The format of libraries for import to ANNOgesic``.

The fasta and gff file of genome.

The gff file of manual detection.

- Arguments

::

    usage: ANNOgesic.py optimize_tsspredator [-h]
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
                                             [--utr_length UTR_LENGTH]
                                             [--lib LIB [LIB ...]]
                                             [--output_prefix OUTPUT_PREFIX [OUTPUT_PREFIX ...]]
                                             [--cluster CLUSTER] [--length LENGTH]
                                             [--core CORE] [--program PROGRAM]
                                             [--replicate_match REPLICATE_MATCH]
                                             [--steps STEPS]
                                             [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --TSSpredator_path TSSPREDATOR_PATH
                            If you want to assign the path of TSSpredator, please
                            assign here.
      --fasta_file FASTA_FILE, -fs FASTA_FILE
                            Path of one target fasta file which you want to
                            opimize it.
      --annotation_file ANNOTATION_FILE, -g ANNOTATION_FILE
                            Path of one target gff file which you want to opimize
                            it.
      --wig_folder WIG_FOLDER, -w WIG_FOLDER
                            The folder of the wig folder.
      --manual MANUAL, -m MANUAL
                            The file of manual checked gff file.
      --strain_name STRAIN_NAME, -n STRAIN_NAME
                            The name of the strain you want to optimize.
      --max_height MAX_HEIGHT, -he MAX_HEIGHT
                            This value relates to the minimal number of read
                            starts at a certain genomic position to be considered
                            as a TSS candidate.During optimization will never
                            larger than this value.
      --max_height_reduction MAX_HEIGHT_REDUCTION, -rh MAX_HEIGHT_REDUCTION
                            When comparing different strains/conditions and the
                            step height threshold is reached in at least one
                            strain/condition, the threshold is reduced for the
                            other strains/conditions by the value set here. This
                            value must be smaller than the step height
                            threshold.During optimization will never larger than
                            this value.
      --max_factor MAX_FACTOR, -fa MAX_FACTOR
                            This is the minimal factor by which the TSS height has
                            to exceed the local expression background.During
                            optimization will never larger than this value.
      --max_factor_reduction MAX_FACTOR_REDUCTION, -rf MAX_FACTOR_REDUCTION
                            When comparing different strains/conditions and the
                            step factor threshold is reached in at least one
                            strain/condition, the threshold is reduced for the
                            other strains/conditions by the value set here. This
                            value must be smaller than the step factor
                            threshold.During optimization will never larger than
                            this value.
      --max_base_height MAX_BASE_HEIGHT, -bh MAX_BASE_HEIGHT
                            This is the minimal number of reads should be mapped
                            on TSS. During optimization will never larger than
                            this value.
      --max_enrichment_factor MAX_ENRICHMENT_FACTOR, -ef MAX_ENRICHMENT_FACTOR
                            This is the minimal enrichment factor. During
                            optimization will never larger than this value.
      --max_processing_factor MAX_PROCESSING_FACTOR, -pf MAX_PROCESSING_FACTOR
                            This is the minimal processing factor. If untreated
                            library is higher than the treated library and above
                            which the TSS candidate is considered as a processing
                            site and not annotated as detected. During
                            optimization will never larger than this value.
      --utr_length UTR_LENGTH, -u UTR_LENGTH
                            The length of UTR. It is for Primary and Secondary
                            definition.
      --lib LIB [LIB ...], -l LIB [LIB ...]
                            The libraries of wig files for TSSpredator. The format
                            is: wig_file_name:tex_treat_or_not(tex or notex):condi
                            tion_id(integer):replicate_id(alphabet):strand(+ or
                            -).
      --output_prefix OUTPUT_PREFIX [OUTPUT_PREFIX ...], -p OUTPUT_PREFIX [OUTPUT_PREFIX ...]
                            The output prefix of all conditions.
      --cluster CLUSTER, -cu CLUSTER
                            This number if for compare manual detected TSS and
                            prediced one. If the position between manual one and
                            predicted one is smaller or equal than this value, it
                            will only print one of them.
      --length LENGTH, -le LENGTH
                            The length of nts for running optimization.Default is
                            compare whole genome
      --core CORE, -c CORE  How many paralle running do you want to use. Default
                            is 4
      --program PROGRAM, -r PROGRAM
                            The type which you want to run TSSpredator (TSS or
                            Processing_site). Default is TSS
      --replicate_match REPLICATE_MATCH, -rm REPLICATE_MATCH
                            The TSS candidates should match to how many number of
                            the replicates. Default is 1
      --steps STEPS, -s STEPS
                            How many steps do you want to run. Default is 4000
                            runs.

- Output files

Based on the program (TSS/processing site), the output files will be stored in 
``$ANNOgesic/output/TSS/optimized_TSSpredator`` or ``$ANNOgesic/output/processing_site/optimized_TSSpredator``.

``stat.csv`` stores the information of every run. The first column is the number of run.
The second column is the parameter set. ``he`` represent height; ``rh`` represent 
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

``best.csv`` stores the best parameter set. The meaning of all columns are the same as ``stat.csv``.

screenshot
-----------

``screenshot`` will generate batch files based on input gff files for producing screenshot of `IGV <https://www.broadinstitute.org/igv>`_.
When the batch file produced, user just need to open `IGV <https://www.broadinstitute.org/igv>`_, then press tools 
and run batch script. The program will automatically produce screenshot. Then user can refer to them easily.

- Pre-required tools and files

`IGV <https://www.broadinstitute.org/igv>`_

Gff files that user want to produce screenshots. All screenshots will be produced based on the position of ``main_gff``. 
``side_gffs`` are the gff files that user want to compare with ``main_gff``. They will also show in the screenshots.

Fasta file of genome.

The libraries and wiggle files. Please refer to the ``The format of libraries for import to ANNOgesic``.

- Arguments

::

    usage: ANNOgesic.py screenshot [-h] [--main_gff MAIN_GFF]
                                   [--side_gffs SIDE_GFFS [SIDE_GFFS ...]]
                                   [--fasta FASTA]
                                   [--frag_wig_folder FRAG_WIG_FOLDER]
                                   [--tex_wig_folder TEX_WIG_FOLDER]
                                   [--height HEIGHT]
                                   [--tex_libs TEX_LIBS [TEX_LIBS ...]]
                                   [--frag_libs FRAG_LIBS [FRAG_LIBS ...]]
                                   [--present PRESENT]
                                   [--output_folder OUTPUT_FOLDER]
                                   [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --main_gff MAIN_GFF, -mg MAIN_GFF
                            Screenshot will based on the position of main_gff file
                            to generate screenshot.
      --side_gffs SIDE_GFFS [SIDE_GFFS ...], -sg SIDE_GFFS [SIDE_GFFS ...]
                            If you have more than one gff want to plot together,
                            please assign here.
      --fasta FASTA, -f FASTA
                            The path of genome fasta folder.
      --frag_wig_folder FRAG_WIG_FOLDER, -fw FRAG_WIG_FOLDER
                            If you want to include the information of fragmented
                            wig file, please assign the folder.
      --tex_wig_folder TEX_WIG_FOLDER, -tw TEX_WIG_FOLDER
                            If you want to include the information of tex treated
                            wig file, please assign the folder.
      --height HEIGHT, -he HEIGHT
                            You can assign the height of screenshot. Default is
                            1500.
      --tex_libs TEX_LIBS [TEX_LIBS ...], -tl TEX_LIBS [TEX_LIBS ...]
                            If you want to include the tex treated wig file,
                            please also assign proper format here. The format is:
                            wig_file_name:tex_treat_or_not(tex or notex):condition
                            _id(integer):replicate_id(alphabet):strand(+ or -).
      --frag_libs FRAG_LIBS [FRAG_LIBS ...], -fl FRAG_LIBS [FRAG_LIBS ...]
                            If you want to include the fragmented wig file, please
                            also assign proper format here. The format is: wig_fil
                            e_name:fragmented(frag):condition_id(integer):replicat
                            e_id(alphabet):strand(+ or -)..
      --present PRESENT, -p PRESENT
                            Which type you want to present in the screen shot.
                            expand/collapse/squish.
      --output_folder OUTPUT_FOLDER, -o OUTPUT_FOLDER
                            Please assign the output folder. If the folder does
                            not exist, it will generate automatically.

- Output files

Based on the path of ``main_gff``, ``screenshot`` will generate a folder - ``screenshots`` under the 
folder of ``main_gff``. All output files will be stored in this folder.

``forward.txt`` is the batch file of forward strand.

``reverse.txt`` is the batch file of reverse strand.

``forward`` is the folder for storing screenshots of forward strand.

``reverse`` is the folder for storing screenshots of reverse strand.

When user run batch files on IGV. The screenshots will automatically store in ``forward`` and ``reverse``. 
The filename will be ``$STRAIN:$START-$END.png``. For example, ``NC_007795:1051529-1051696.png`` 
means the strain of the screenshot is NC_007795. The feature's start point is 1051529 and the end point is 
1051696.

color_png
----------

``color_png`` is a following procedure of ``screenshot``. If there are many wiggle files, it will be difficult 
to distinguish the tracks. Especially, when user want to compare TEX +/- libraries to check TSS or processing sites, 
it is not convenient. ``color_png`` can color the tracks based on TEX +/- libraries. Therefore, the user can 
refer to the screenshots much easier.

- Pre-required tools and files

`ImageMagick <http://www.imagemagick.org/script/index.php>`_

The screenshots which generated from ``screenshot``. Please make sure the folder of ``forward`` and ``reverse`` 
exist in the folder of ``screenshots``.

- Arguments

::

    usage: ANNOgesic.py color_png [-h] [--screenshot_folder SCREENSHOT_FOLDER]
                                  [--track_number TRACK_NUMBER]
                                  [--ImageMagick_covert_path IMAGEMAGICK_COVERT_PATH]
                                  [project_path]
    
    positional arguments:
      project_path          Path of the project folder. If none is given the
                            current directory is used.
    
    optional arguments:
      -h, --help            show this help message and exit
      --screenshot_folder SCREENSHOT_FOLDER, -f SCREENSHOT_FOLDER
                            The folder which stores the folder of screenshots.
      --track_number TRACK_NUMBER, -t TRACK_NUMBER
                            How many number of tracks do you have.
      --ImageMagick_covert_path IMAGEMAGICK_COVERT_PATH, -m IMAGEMAGICK_COVERT_PATH
                            Please assign the path of covert in ImageMagick
                            package.

- Output files

The new screenshots will replace the previous ones automatically.
