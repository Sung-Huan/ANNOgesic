Tutorial of ANNOgesic
=====================

This tutorial guids you through a small test case to show you how to
use ANNOgesic's subcommands. It builds on a differential RNA-Seq data
set from *Campylobacter jejuni* subsp. jejuni 81116 which can be
downloaded from NCBI GEO and that was part of a work by `Dugar et
al. <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38883>`_.
As part of the tutorial several output files will be generated in
different formats including CSV files (tabular separated plain text
files), that can be opened in LibreOffice or Excel, GFF3 files, plain
text files and figures. For the viewing GFF3 files you can use genome
browsers like `IGB <http://bioviz.org/igb/index.html>`_, `IGV
<https://www.broadinstitute.org/igv/>`_ or `Jbrowse
<http://jbrowse.org/>`_.

Before we start, please check :ref:`The format of filename` and
:ref:`The input format of libraries for running ANNOgesic`. All the 
parameters and details can be found in the
section :ref:`ANNOgesic's subcommands`.
Moreover, all the requirements are listed in the section
:ref:`Required tools or databases`. 

If a subcommand requires third-party softwares (e.g. TSSpredator in
the subcommand ``tss_ps``) make sure that path of the executable file
is properly specified or is part of the environmental variable
``$PATH``. Moreover, if ``annogesic`` (executable file of ANNOgesic)
is not in your ``$PATH``, please specify its full path.

Using Singularity (optional)
----------------------------

Using `Singularity <https://singularity.lbl.gov/index.html>`_  to build an 
image of ANNOgesic can install all required tools and 
environment automatically. If you would like to use Singularity to run ANNOgesic, 
please check the following commands.


First of all, we need to build up an image of ANNOgesic via Docker image.


::

    $ singularity build \
        annogesic.img \
        docker://silasysh/annogesic:latest

Now, we can use Singularity to run ANNOgesic. You just need to add the following 
line before the command that you want to run.

::

    singularity exec -B $STORAGE_PATH annogesic.img

Please put the storage path of your home directory to ``$STORAGE_PATH``. ``df`` can be used to check the
storage system.

For example, if you want to creat the analyzing folder called ``ANNOgesic`` and your storage path is ``storage1``, you can 
use the following command to run ANNOgesic via Singularity.

::

     singularity exec -B /storage1 annogesic.img \
        annogesic create -pj ANNOgesic

Run scripts for tutorial
------------------------

We also provide a shell script which stores all commands that were used in 
this tutorial. You can run/comment the module by simply add/remove ``#`` in front of
the module. The run script (file name is run.sh) can be found `here
<https://github.com/Sung-Huan/ANNOgesic/tree/master/tutorial_data>`_. If you are 
using Singularity, the file name of run script is run_with_singularity.sh which can 
be found in the same folder.

Generating a project
--------------------

First of all, we need to create a working folder (``--project_path``) by running ``create``.

::

    $ annogesic create --project_path ANNOgesic

This will create the following folder structure:

::


   $ tree ANNOgesic
   ANNOgesic
   ├── input
   │   ├── BAMs
   │   │   ├── BAMs_map_reference_genomes
   │   │   │   ├── fragment
   │   │   │   └── tex_notex
   │   │   └── BAMs_map_related_genomes
   │   │       ├── fragment
   │   │       └── tex_notex
   │   ├── databases
   │   ├── manual_processing_sites
   │   ├── manual_TSSs
   │   ├── mutation_tables
   │   ├── reads
   │   ├── references
   │   │   ├── annotations
   │   │   └── fasta_files
   │   ├── riboswitch_ID_file
   │   ├── RNA_thermometer_ID_file
   │   └── wigs
   │       ├── fragment
   │       └── tex_notex
   ├── output
   └── used_annogesic_version.txt

   22 directories, 1 file


Retrieving the genome sequences and annotation files
----------------------------------------------------

For our test case, the genome and annotation data has to be retrieved
`NCBI <ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Campylobacter_jejuni/latest_assembly_versions/GCF_000017905.1_ASM1790v1/>`_.

ANNOgesic offers a convenient to do that and we download fasta files
(``--ref_fasta``), gff files (``--ref_gff``), gbk files
(``--ref_gbk``), ptt files (``--ref_ptt``), rnt files (``--ref_rnt``),
and convert the the gff files to embl format (``--convert_embl``).

::

    $ annogesic get_input_files \
        --ftp_path ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Campylobacter_jejuni/latest_assembly_versions/GCF_000017905.1_ASM1790v1/ \
        --ref_gff --ref_fasta --ref_gbk --ref_ptt --ref_rnt --convert_embl \
        --project_path ANNOgesic

The file will be place in the following locations:

::

    $ ls ANNOgesic/input/references/fasta_files/
    NC_009839.1.fa
    $ ls ANNOgesic/input/references/annotations/
    NC_009839.1.embl  NC_009839.1.gbk  NC_009839.1.gff

Alternatively you can manually copy the file into these subfolders.

Retrieving wiggle and read files
--------------------------------

We need to download reads in SRA format form GEO
`here <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38883>`_.

::

    $ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP013/SRP013869/SRR515254/SRR515254.sra
    $ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP013/SRP013869/SRR515255/SRR515255.sra
    $ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP013/SRP013869/SRR515256/SRR515256.sra
    $ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP013/SRP013869/SRR515257/SRR515257.sra

Then we can convert SRA files to Fasta or Fastq format for mapping by
using ``fastq-dump`` of the the `SRA toolkit
<http://www.ncbi.nlm.nih.gov/books/NBK158900/>`_.

::
  
   $ fastq-dump --fasta SRR515254.sra
   $ fastq-dump --fasta SRR515255.sra
   $ fastq-dump --fasta SRR515256.sra
   $ fastq-dump --fasta SRR515257.sra
   $ mv *.fasta ANNOgesic/input/reads
   $ rm SRR515254.sra SRR515255.sra SRR515256.sra SRR515257.sra

Then we have to download the coverage files in wiggle format.

::

   $ wget -cP ANNOgesic/input/wigs/tex_notex ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM951nnn/GSM951380/suppl/GSM951380%5FLog%5F81116%5FR1%5Fminus%5FTEX%5Fin%5FNC%5F009839%5Fminus.wig.gz
   $ wget -cP ANNOgesic/input/wigs/tex_notex ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM951nnn/GSM951380/suppl/GSM951380%5FLog%5F81116%5FR1%5Fminus%5FTEX%5Fin%5FNC%5F009839%5Fplus.wig.gz
   $ wget -cP ANNOgesic/input/wigs/tex_notex ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM951nnn/GSM951381/suppl/GSM951381%5FLog%5F81116%5FR1%5Fplus%5FTEX%5Fin%5FNC%5F009839%5Fminus.wig.gz
   $ wget -cP ANNOgesic/input/wigs/tex_notex ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM951nnn/GSM951381/suppl/GSM951381%5FLog%5F81116%5FR1%5Fplus%5FTEX%5Fin%5FNC%5F009839%5Fplus.wig.gz
   $ cd ANNOgesic/input/wigs/tex_notex
   $ gunzip GSM951380_Log_81116_R1_minus_TEX_in_NC_009839_minus.wig.gz \
            GSM951380_Log_81116_R1_minus_TEX_in_NC_009839_plus.wig.gz \
            GSM951381_Log_81116_R1_plus_TEX_in_NC_009839_minus.wig.gz \
            GSM951381_Log_81116_R1_plus_TEX_in_NC_009839_plus.wig.gz
   $ cd ../../../../

If we check the wiggle files, we will find that the fasta filename (presented by "chrom") is not the same as fasta or annotation gff file.

::

   $ head ANNOgesic/input/wigs/tex_notex/GSM951380_Log_81116_R1_minus_TEX_in_NC_009839_minus.wig 
     track type=wiggle_0 name="Log_81116_R1_minus_TEX_in_NC_009839"
     variableStep chrom=NC_009839 span=1
     7	-1.0
     8	-1.0
     9	-1.0
     10	-1.0
     11	-1.0
     12	-1.0
     13	-1.0
     14	-1.0

Our genome fasta file is NC_009839.1.fa. Thus "chrom" in wiggle file should be NC_009839.1 not NC_009839. 
We can use `replace_seq_id.py <https://github.com/Sung-Huan/ANNOgesic/tree/master/tutorial_data>`_ from our 
Git repository to replace the genome name in wiggle files properly. If the genome names in your fasta, annotation, 
wiggle files are the same, you don't need to do this step.

::

   $ wget https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/tutorial_data/replace_seq_id.py
   $ python3 replace_seq_id.py -i ANNOgesic/input/wigs/tex_notex -n NC_009839.1
   $ rm replace_seq_id.py

Since this is a tutorial, we only download one replicate to reduce the running time.

Mapping (optional)
------------------

ANNOgesic needs several input files, such as alignment files and wiggle files. In this tutorial, all 
the required files can be downloaded from public database. However, the users may have their own dataset. 
In order to generate these required files, we recommand the users to use READemption for mapping their reads. 
The following example is for  mapping the read files of NC_009839.1 that we just retrieved from SRA.

First of all, we need to create an folder for the analysis.

::

    $ reademption create READemption

Then, we need to move or copy the required input files to the corresponding folders.

::

    $ cp ANNOgesic/input/reads/* READemption/input/reads
    $ cp ANNOgesic/input/references/fasta_files/*.fa READemption/input/reference_sequences
    $ cp ANNOgesic/input/references/annotations/*.gff READemption/input/annotations

Now, the mapping can be performed..

::

    $ reademption \
        align \
        -r \
        -S \
        -a 95 \
        -l 12 \
        --poly_a_clipping \
        --progress \
        READemption

We can assign ``-p`` for running parallel, such as ``-p 20``.

When the alignment is done, the alignemnt files in BAM format can be found in 
``READemption_analysis/output/align/alignments/``.

Now, we can generate coverage files in wiggle format.

::

    $ reademption \
        coverage \
        READemption

``-p`` can be assigned for ``coverage`` as well. The raw and normalized coverage files can be 
found in ``READemption_analysis/output/coverage``.

For the details of using READemption, please check 
`READemption <http://pythonhosted.org/READemption/>`_

Improving the reference genome
------------------------------

If the fasta file of the target reference genome is not available, ``update_genome_fasta`` 
can generate it from a related genome.

Now, we assume that we need to generate fasta file of the reference genome. 
First of all, we need to find a closely related genome (fasta file and gff file can be found) of our target genome. 
Then, we need to generate a mutation table (please check the section :ref:`ANNOgesic's subcommands`) 
between these two genomes. When these files are available, 
we can run subcommand ``update_genome_fasta`` for getting fasta file of the target genome. 

A simple example of the mutation table can be found in 
`mutation table <https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/tutorial_data/mutation.csv>`_.
Each column of the table is separated by tab. 

::

     $ wget https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/tutorial_data/mutation.csv
     $ mv mutation.csv ANNOgesic/input/mutation_tables

We assume NC_009839.1.fa is a related genome of our     
target genome -- NC_test.1 and test_case2. The fasta files of the new genomes (NC_test.1 and test_case2)
will be generated in ``ANNOgesic/output/updated_references/fasta_files``.

::

     $ annogesic update_genome_fasta \
        --related_fasta_files ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        --mutation_table ANNOgesic/input/mutation_tables/mutation.csv \
        --updated_seq_name NC_test.1 \
        --project_path ANNOgesic

``--related_fasta_files`` is path of the fasta file of closely related genome. 
When the running process is done, the following information will appear.

::

    $ Updating the reference sequences
      Please use the new fasta file to remapping again.

Since the data (``ANNOgesic/output/updated_references/fasta_files``) that we generated is only a dummy data,
we can ignore the information now. Otherwise,
you have to re-map again in order to get the correct alignment and coverage files.

Now we can check the results.

::

    $ head ANNOgesic/input/references/fasta_files/NC_009839.1.fa
    >NC_009839.1
    ATGAATCCAAATCAAATACTTGAAAATTTAAAAAAAGAATTAAGTGAAAACGAATACGAAAATTATATCGCTATCTTAAA
    ATTTAACGAAAAACAAAGCAAAGCAGATTTTCTAGTCTTTAACGCTCCTAATGAGCTTTTAGCCAAATTCATACAAACAA
    AATACGGTAAAAAAATTTCACATTTTTATGAAGTACAAAGCGGAAATAAAGCGAGCGTTTTGATACAAGCACAAAGTGCT
    AAACAAAGTAGCAAAAGCACTAAAATCGATATCGCTCATATCAAGGCGCAAAGTACGATTTTAAATCCTTCTTTTACTTT
    TGAAAGCTTTGTAGTGGGGGATTCTAACAAATACGCTTATGGAGCTTGTAAAGCTATCTCACAAAAAGACAAACTGGGAA
    AACTTTATAATCCTATCTTTATCTATGGGCCTACAGGGCTTGGAAAAACGCACTTGCTTCAAGCTGTGGGAAATGCAAGT
    TTGGAAATGGGAAAAAAAGTGATTTATGCTACGAGTGAAAATTTTATCAATGATTTTACTTCAAATTTAAAAAATGGCTC
    TTTAGATAAATTTCACGAAAAATATAGAAATTGTGATGTTTTACTCATAGATGATGTGCAGTTTTTAGGAAAAACCGATA
    AAATTCAAGAAGAATTTTTCTTTATATTTAATGAAATCAAAAATAACGATGGACAAATCATCATGACTTCAGACAATCCA
    $ head ANNOgesic/output/updated_references/fasta_files/NC_test.1.fa
    >NC_test.1
    ATcAACCAAATCAAATACTTGAAAATTTAAAAAAAGAATTAAGTGAAAACGAATACGAAA
    ATTATATCGCTATCTTAAAATTTAACGAAAAACAAAGCAAAGCAGATTTTCTAGTCTTTA
    ACGCTCCTAATGAGCTTTTAGCCAAATTCATACAAACAAAATACGGTAAAAAAATTTCAC
    ATTTTTATGAAGTACAAAGCGGAAATAAAGCGAGCGTTTTGATACAAGCACAAAGTGCTA
    AACAAAGTAGCAAAAGCACTAAAATCGATATCGCTCATATCAAGGCGCAAAGTACGATTT
    TAAATCCTTCTTTTACTTTTGAAAGCTTTGTAGTGGGGGATTCTAACAAATACGCTTATG
    GAGCTTGTAAAGCTATCTCACAAAAAGACAAACTGGGAAAACTTTATAATCCTATCTTTA
    TCTATGGGCCTACAGGGCTTGGAAAAACGCACTTGCTTCAAGCTGTGGGAAATGCAAGTT
    TGGAAATGGGAAAAAAAGTGATTTATGCTACGAGTGAAAATTTTATCAATGATTTTACTT

We can see the third nucleotide of ``NC_test.1.fa`` is replaced from G to c. Moreover, The sixth nucleotide T is deleted.

If the mutation table can not be provided, we can also use subcommand ``snp`` to detect mutations and generate 
fasta files automatically. For ``snp``, we will go through it later.

Generating annotation files
---------------------------

If the genome annotations of the target reference genome in GFF format is not available, ``annotation_transfer``
can generate it from a related genome.

Like the last section, we assume NC_009839.1.fa is a related genome of our
target genome -- NC_test.1 and test_case2.

Before we running this subcommand, we have to modify environment paths of `RATT <http://ratt.sourceforge.net/>`_. 
If you run ANNOgesic in docker container, the path is already set. 
Otherwise, please check 
`RATT <http://ratt.sourceforge.net/>`_ to set your environment paths properly. Moreover, if 
the error message related to 'defined(@array)' occurs, please check :ref:`FAQ`.

After setting the environment, we can try it.

::

    anngesic annotation_transfer \
        --related_embl_files ANNOgesic/input/references/annotations/NC_009839.1.embl \
        --related_fasta_files ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        --target_fasta_files ANNOgesic/output/updated_references/fasta_files/NC_test.1.fa \
        --element chromosome \
        --transfer_type Strain \
        --compare_pair NC_009839.1:NC_test.1 \
        --project_path ANNOgesic


``--element`` is prefix name of the output embl files. 
``--transfer_type`` is a program of `RATT <http://ratt.sourceforge.net/>`_.
We use ``Strain`` because the similarity between two genomes is higher than 90% (please check 
`RATT <http://ratt.sourceforge.net/>`_). In ``--compare_pair``, the pairs of the genomes 
(NC_test.1 and test_case2) and their closely related genomes (NC_000915.1) are assigned. 
The annotation information in embl, GFF3, ptt, and rnt format will be stored in 
``ANNOgesic/output/updated_references/annotations``.

Once the transfer is done, we can see

::

    $ ls ANNOgesic/output/updated_references/annotations/
    NC_test.1.gff  NC_test.1.ptt  NC_test.1.rnt
    $ ls ANNOgesic/output/annotation_transfer/
    chromosome.NC_test.1.final.embl  log.txt  NC_test.1.gff  ratt_log.txt

In ``ANNOgesic/output/updated_references/annotations``, we can find ptt, rnt and gff files. In ``ANNOgesic/output/annotation_transfer``,
we can find the output of `RATT <http://ratt.sourceforge.net/>`_.

We already saw how to update genome fasta and annotation files. 
For the following subcommands, we will use ``ANNOgesic/input/references/annotations/NC_009839.1.gff`` 
and ``ANNOgesic/input/references/fasta_files/NC_009839.1`` as the reference genome.

TSS and processing site prediction and optimization
---------------------------------------------------

Before running following subcommands, we need to setup our libraries in a correct format.
First, we set the paths of wig files.

::

    WIG_FOLDER="ANNOgesic/input/wigs/tex_notex"

Then, we can setup our libraries -- ``$WIGFILE:$TEXorNOTEXorFRAG:CONDITION:REPLICATE:STRAND`` 
(:ref:`The input format of libraries for running ANNOgesic`).

::

    TEX_LIBS="$WIG_FOLDER/GSM951380_Log_81116_R1_minus_TEX_in_NC_009839_minus.wig:notex:1:a:- \
              $WIG_FOLDER/GSM951381_Log_81116_R1_plus_TEX_in_NC_009839_minus.wig:tex:1:a:- \
              $WIG_FOLDER/GSM951380_Log_81116_R1_minus_TEX_in_NC_009839_plus.wig:notex:1:a:+ \
              $WIG_FOLDER/GSM951381_Log_81116_R1_plus_TEX_in_NC_009839_plus.wig:tex:1:a:+"

Now, we can start to test other subcommands.
 
Before running ``tss_ps``, we can use ``optimize_tss_ps`` to optimize the parameters 
(It is an optional step, but we highly recommand you to do it). 
The optimization requires a small set of the manual-detected TSSs in GFF3 format. 
In our experience, we recommend you to detect at least 50 TSSs and check more than 200kb of the genome. 

For this test case, you can download the `manual TSS file <https://github.com/Sung-Huan/ANNOgesic/tree/master/tutorial_data>`_ 
from our git repository. 

::

    $ wget -cP ANNOgesic/input/manual_TSSs/ https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/tutorial_data/NC_009839_manual_TSS.gff

Now, we have a manual-detected TSS gff file which is stored in ``ANNOgesic/input/manual_TSSs``. 
we can try ``optimize_tss_ps`` right now (since we only check first 200kb, we set ``--genome_lengths`` 
as "NC_009839.1:200000" which means only first 200kb of NC_009839.1 is valid.).

::

    $ annogesic optimize_tss_ps \
        --fasta_files ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        --annotation_files ANNOgesic/input/references/annotations/NC_009839.1.gff \
        --tex_notex_libs $TEX_LIBS \
        --condition_names TSS --steps 25 \
        --manual_files ANNOgesic/input/manual_TSSs/NC_009839_manual_TSS.gff \
        --curated_sequence_length NC_009839.1:200000 \
        --replicate_tex all_1 \
        --project_path ANNOgesic

``optimize_tss_ps`` will compare manual-checked TSSs with predicted TSSs to search the optimized parameters. 
Results of the different parameters will be printed in the screen, and stored in ``stat_NC_009839.1.csv`` as well. 
We only set 25 runs for testing. For optimization of processing sites, we just need to change ``--program`` from TSS to PS. 
``--replicate_tex`` means the minimum replicates that a TSS can be detected. ``all_1`` means that a TSS 
should be detected in at least one replicate in all conditions. ``--replicate_tex`` can be also assigned like ``all_2`` 
(a TSS should be detected in at least two replicates in all conditions) 
or ``1_2 2_2 3_3`` (in condition 1 and 2 (based on the setting of ``--tex_notex_libs``, 
a TSS should be detected in at least two replicates, and a TSS should be predicted in three replicates in condition 3).
Once the optimization is done, you can find several files.

::

    $ ls ANNOgesic/output/TSSs/optimized_TSSpredator/
    best_NC_009839.1.csv  log.txt  results_all_steps.txt   stat_NC_009839.1.csv

``best_NC_009839.1.csv`` is for the results of the optimized parameters; ``stat_NC_009839.1.csv`` is for the results of each step.

Now, we assume the optimized parameters are following: height is 0.4, height_reduction is 0.1, factor is 1.7, factor_reduction is 0.2, 
base_height is 0.039, enrichment_factor is 1.1, processing_factor is 4.5. We can set these parameters for running  
``tss``.

::

    $ annogesic tss_ps \
        --fasta_files ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        --annotation_files ANNOgesic/input/references/annotations/NC_009839.1.gff \
        --tex_notex_libs $TEX_LIBS \
        --condition_names test \
        --height 0.4 \
        --height_reduction 0.1 \
        --factor 1.7 \
        --factor_reduction 0.2 \
        --base_height 0.039 \
        --enrichment_factor 1.1 \
        --processing_factor 4.5 \
        --validate_gene \
        --program TSS \
        --replicate_tex all_1 \
        --curated_sequence_length NC_009839.1:200000 \
        --manual_files ANNOgesic/input/manual_TSSs/NC_009839_manual_TSS.gff \
        --project_path ANNOgesic

We assigned the manual-checked TSS gff file to ``--manual_files``. Therefore,
the output gff file will contain the manual-detected TSSs and predicted TSSs.
If we didn't assign it, Only the predicted TSSs will be included in output gff file.

If you did not run ``optimize_tss_ps`` before and just want to do TSS prediction with 
default parameters, ``--height``, ``--height_reduction``, ``--factor``, ``--factor_reduction``,
``--base_height``, and ``--enrichment_factor`` do not need to be assigned.

::

    $ annogesic tss_ps \
        --fasta_files ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        --annotation_files ANNOgesic/input/references/annotations/NC_009839.1.gff \
        --tex_notex_libs $TEX_LIBS \
        --condition_names test \
        --validate_gene \
        --program TSS \
        --replicate_tex all_1 \
        --curated_sequence_length NC_009839.1:200000 \
        --manual_files ANNOgesic/input/manual_TSSs/NC_009839_manual_TSS.gff \
        --project_path ANNOgesic


The output files are gff file, MasterTable and statistic files.

::

    $ ls ANNOgesic/output/TSSs/
    configs  gffs  MasterTables  log.txt  optimized_TSSpredator  screenshots  statistics
    $ ls ANNOgesic/output/TSSs/configs/
    config_NC_009839.1.ini
    $ ls ANNOgesic/output/TSSs/gffs/
    NC_009839.1_TSS.gff
    $ ls ANNOgesic/output/TSSs/MasterTables/MasterTable_NC_009839.1/
    AlignmentStatistics.tsv  err.txt  log.txt  MasterTable.tsv  superConsensus.fa  superTSS.gff  superTSStracks.gff  test_super.fa  test_super.gff  test_TSS.gff
    $ ls ANNOgesic/output/TSSs/statistics/NC_009839.1/
    stat_compare_TSSpredator_manual_NC_009839.1.csv  stat_TSS_class_NC_009839.1.csv  TSS_class_NC_009839.1.png  TSS_venn_NC_009839.1.png
    stat_gene_vali_NC_009839.1.csv                   stat_TSS_libs_NC_009839.1.csv   TSSstatistics.tsv

If we want to predict processing sites, the procedures are the same. We just need to change the program from TSS to 
PS (``--program``) and assign the proper parameter sets. We assume the best parameter sets are following: 
height is 0.2, height_reduction is 0.1, factor is 2.0, factor_reduction is 0.5,
base_height is 0.009, enrichment_factor is 1.2, processing_factor is 1.5.

::

    $ annogesic tss_ps \
        --fasta_files ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        --annotation_files ANNOgesic/input/references/annotations/NC_009839.1.gff \
        --tex_notex_libs $TEX_LIBS \
        --condition_names test \
        --height 0.2 \
        --height_reduction 0.1 \
        --factor 2.0 \
        --factor_reduction 0.5 \ 
        --base_height 0.009 \
        --enrichment_factor 1.2 \
        --processing_factor 1.5 \ 
        --replicate_tex all_1 \
        --program PS \
        --project_path ANNOgesic

The output files are following:

::

    $ ls ANNOgesic/output/processing_sites/
    configs  gffs  MasterTables  log.txt  statistics
    $ ls ANNOgesic/output/processing_sites/configs/
    config_NC_009839.1.ini
    $ ls ANNOgesic/output/processing_sites/gffs/
    NC_009839.1_processing.gff
    $ ls ANNOgesic/output/processing_sites/MasterTables/MasterTable_NC_009839.1/
    AlignmentStatistics.tsv  err.txt  log.txt  MasterTable.tsv  superConsensus.fa  superTSS.gff  superTSStracks.gff  test_super.fa  test_super.gff  test_TSS.gff
    $ ls ANNOgesic/output/processing_sites/statistics/NC_009839.1/
    processing_class_NC_009839.1.png  processing_venn_NC_009839.1.png  stat_processing_class_NC_009839.1.csv  stat_processing_libs_NC_009839.1.csv  TSSstatistics.tsv

Since we use TSSpredator to detect processing site, the files in 
``ANNOgesic/output/processing_sites/MasterTables/MasterTable_NC_009839.1/`` are for processing site not for TSS.

Performing transcript detection
-------------------------------

Transcript detection is a basic procedure for detecting transcript boundary. 
we can use subcommand ``transcript`` to detect the transcript. Normally, we strongly 
recommend that the user should provide the libraries of RNA-Seq with transcript fragmented 
(``--frag_libs``) because dRNA-Seq focus on 5'end and usually loses some information 
of 3'end. However, we only use TEX +/- for testing since we have no fragmented libraries.

There are several options for modifying transcripts by comparing transcripts and genome annotations like CDSs (``--modify_transcript``). 
By assigning ``--modify_transcript``, transcripts can be merged or extended based on the genome annotations.
If you want to know the details, please check :ref:`transcript`. Now, we use default setting to run this module: 

::

    $ annogesic transcript \
        --annotation_files ANNOgesic/input/references/annotations/NC_009839.1.gff \
        --tex_notex_libs $TEX_LIBS \
        --replicate_tex all_1 \
        --compare_feature_genome gene CDS \
        --tss_files ANNOgesic/output/TSSs/gffs/NC_009839.1_TSS.gff \
        --project_path ANNOgesic

The output files are gff files, tables and statistic files.

::

    $ ls ANNOgesic/output/transcripts/gffs
    NC_009839.1_transcript.gff
    $ ls ANNOgesic/output/transcripts/tables
    NC_009839.1_transcript.csv
    $ ls ANNOgesic/output/transcripts/statistics
    NC_009839.1_length_all.png  NC_009839.1_length_less_2000.png  stat_compare_transcript_TSS_NC_009839.1.csv  stat_compare_transcript_genome_NC_009839.1.csv

Prediction of terminator
------------------------

We can run subcommand ``terminator`` to detect terminators. The command is like following: 

::

    $ annogesic terminator \
        --fasta_files ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        --annotation_files ANNOgesic/input/references/annotations/NC_009839.1.gff \
        --transcript_files ANNOgesic/output/transcripts/gffs/NC_009839.1_transcript.gff \
        --tex_notex_libs $TEX_LIBS \
        --replicate_tex all_1 \
        --project_path ANNOgesic

Four different kinds of gff files and tables will be generated.

::

    $ ls ANNOgesic/output/terminators/gffs/
    all_candidates  best_candidates  expressed_candidates  non_expressed_candidates
    $ ls ANNOgesic/output/terminators/tables
    all_candidates  best_candidates  expressed_candidates  non_expressed_candidates

``all_candidates`` is for all candidates; ``expressed_candidates`` is for the candidates which reveal gene expression; 
``best_candidates`` is for the candidates which reveal gene expression and their coverages show significant decreasing; 
``non_expressed_candidates`` is for the candidates which have no expression.

::

    $ ls ANNOgesic/output/terminators/gffs/best_candidates
    NC_009839.1_term.gff
    $ ls ANNOgesic/output/terminators/gffs/expressed_candidates
    NC_009839.1_term.gff
    $ ls ANNOgesic/output/terminators/gffs/all_candidates
    NC_009839.1_term.gff
    $ ls ANNOgesic/output/terminators/gffs/non_expressed_candidates
    NC_009839.1_term.gff
    $ ls ANNOgesic/output/terminators/tables/best_candidates
    NC_009839.1_term.csv
    $ ls ANNOgesic/output/terminators/tables/expressed_candidates
    NC_009839.1_term.csv
    $ ls ANNOgesic/output/terminators/tables/all_candidates
    NC_009839.1_term.csv
    $ ls ANNOgesic/output/terminators/tables/non_expressed_candidates
    NC_009839.1_term.csv

In transtermhp folder, output files of `TranstermHP <http://transterm.cbcb.umd.edu/>`_ can be found.

::

    $ ls ANNOgesic/output/terminators/transtermhp_results/NC_009839.1
    NC_009839.1_best_terminator_after_gene.bag  NC_009839.1_terminators.txt  NC_009839.1_terminators_within_robust_tail-to-tail_regions.t2t

Moreover, statistic files are stored in ``statistics``.

::

    $ ls ANNOgesic/output/terminators/statistics/
    stat_compare_terminator_transcript_NC_009839.1_all_candidates.csv   stat_compare_terminator_transcript_NC_009839.1_expressed_candidates.csv
    stat_compare_terminator_transcript_NC_009839.1_best_candidates.csv  stat_NC_009839.1.csv

Generating UTR
--------------

Now, we have the information of TSSs, transcripts and terminators. We can detect the 5'UTRs and 3'UTRs by using 
subcommand ``utr``.

::

    $ annogesic utr \
        --annotation_files ANNOgesic/input/references/annotations/NC_009839.1.gff \
        --tss_files ANNOgesic/output/TSSs/gffs/NC_009839.1_TSS.gff \
        --transcript_files ANNOgesic/output/transcripts/gffs/NC_009839.1_transcript.gff \
        --terminator_files ANNOgesic/output/terminators/gffs/best_candidates/NC_009839.1_term.gff \
        --project_path ANNOgesic

If the TSS gff file is not generated by ANNOgesic, please add ``--tss_source`` which can classify TSSs for generating UTRs.
Output gff files and statistic files will be stored in ``ANNOgesic/output/UTRs/5UTRs`` and ``ANNOgesic/output/UTRs/3UTRs``.

::

    $ ls ANNOgesic/output/UTRs/3UTRs
    gffs/       statistics/
    $ ls ANNOgesic/output/UTRs/5UTRs
    gffs/       statistics/
    $ ls ANNOgesic/output/UTRs/3UTRs/gffs
    NC_009839.1_3UTR.gff
    $ ls ANNOgesic/output/UTRs/5UTRs/gffs
    NC_009839.1_5UTR.gff
    $ ls ANNOgesic/output/UTRs/5UTRs/statistics
    NC_009839.1_all_5utr_length.png
    $ ls ANNOgesic/output/UTRs/3UTRs/statistics
    NC_009839.1_all_3utr_length.png

Now, we have all information for defining the transcript boundary.

Detecting operon and suboperon
------------------------------

We have TSSs, transcripts, terminators, CDSs, UTRs now. We can integrate all these feature to 
detect operons and suboperons by executing subcommand ``operon``.

::

    $ annogesic operon \
        --annotation_files ANNOgesic/input/references/annotations/NC_009839.1.gff \
        --tss_files ANNOgesic/output/TSSs/gffs/NC_009839.1_TSS.gff \
        --transcript_files ANNOgesic/output/transcripts/gffs/NC_009839.1_transcript.gff \
        --terminator_files ANNOgesic/output/terminators/gffs/best_candidates/NC_009839.1_term.gff \
        --project_path ANNOgesic

Three folders will be generated to store table and statistics files.

::

    $ ls ANNOgesic/output/operons/
    gffs  log.txt  statistics  tables 
    $ ls ANNOgesic/output/operons/gffs/
    NC_009839.1_operon.gff
    $ ls ANNOgesic/output/operons/tables/
    NC_009839.1_operon.csv
    $ ls ANNOgesic/output/operons/statistics/
    stat_NC_009839.1_operon.csv

Promoter motif detection
------------------------

As long as we have TSSs, we can use subcommand ``promoter`` to identify promoters. If the TSS gff files are not generated by ``ANNOgesic``,
please add ``--tss_source`` and corresponding genome annotation file containing CDSs, tRNAs, rRNAs, etc, (``--annotation_files``) in order to  
classify TSSs for detecting promoters.
Now, let's try ``promoter`` (``--program`` is assigned by "both" in default. If you want to only run 
`MEME <http://meme-suite.org/tools/meme>`_ or `GLAM2 <http://meme-suite.org/tools/glam2>`_, 
please assign "meme" or "glam2" to ``--program``), the process may take a while.

::

    $ annogesic promoter \
        --tss_files ANNOgesic/output/TSSs/gffs/NC_009839.1_TSS.gff \
        --fasta_files ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        --motif_width 45 2-10 \
        --project_path ANNOgesic

We defined the length of the motifs as ``50`` and ``2-10``. ``2-10`` means the width can be from 2 to 10.

If the software for running MEME and GLAM2 in parallels is installed, ``--parallels`` can also be assigned for
running MEME and GLAM2 in parallels.

::

    $ annogesic promoter \
        --tss_files ANNOgesic/output/TSSs/gffs/NC_009839.1_TSS.gff \
        --fasta_files ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        --motif_width 45 2-10 \
        --parallels 10 \
        --project_path ANNOgesic

Based on different types of the TSSs and the length of the motif, numerous output files will be generated.

::

    $ ls ANNOgesic/output/promoters/
    fasta_classes  NC_009839.1   log.txt
    $ ls ANNOgesic/output/promoters/fasta_classes/NC_009839.1
    NC_009839.1_allgenome_all_types.fa  NC_009839.1_allgenome_internal.fa  NC_009839.1_allgenome_primary.fa    NC_009839.1_allgenome_without_orphan.fa
    NC_009839.1_allgenome_antisense.fa  NC_009839.1_allgenome_orphan.fa    NC_009839.1_allgenome_secondary.fa
    $ ls ANNOgesic/output/promoters/NC_009839.1
    MEME GLAM2
    $ ls ANNOgesic/output/promoters/NC_009839.1/MEME
    promoter_motifs_NC_009839.1_allgenome_all_types_2-10_nt  promoter_motifs_NC_009839.1_allgenome_all_types_45_nt
    $ ls ANNOgesic/output/promoters/NC_009839.1/GLAM2
    promoter_motifs_NC_009839.1_allgenome_all_types_2-10_nt  promoter_motifs_NC_009839.1_allgenome_all_types_45_nt
    $ ls ANNOgesic/output/promoters/NC_009839.1/MEME/promoter_motifs_NC_009839.1_allgenome_all_types_45_nt/
    logo1.eps  logo1.png  logo2.eps  logo2.png  logo3.eps  logo3.png  logo_rc1.eps  logo_rc1.png  logo_rc2.eps  logo_rc2.png  logo_rc3.eps  logo_rc3.png  meme.csv  meme.html  meme.txt  meme.xml
    $ ls ANNOgesic/output/promoters/NC_009839.1/GLAM2/promoter_motifs_NC_009839.1_allgenome_all_types_45_nt/
    glam2.csv   glam2.txt   logo1.eps  logo2.png  logo4.eps  logo5.png  logo7.eps  logo8.png  logo_ssc10.eps  logo_ssc1.png  logo_ssc3.eps  logo_ssc4.png  logo_ssc6.eps  logo_ssc7.png  logo_ssc9.eps
    glam2.html  logo10.eps  logo1.png  logo3.eps  logo4.png  logo6.eps  logo7.png  logo9.eps  logo_ssc10.png  logo_ssc2.eps  logo_ssc3.png  logo_ssc5.eps  logo_ssc6.png  logo_ssc8.eps  logo_ssc9.png
    glam2.meme  logo10.png  logo2.eps  logo3.png  logo5.eps  logo6.png  logo8.eps  logo9.png  logo_ssc1.eps   logo_ssc2.png  logo_ssc4.eps  logo_ssc5.png  logo_ssc7.eps  logo_ssc8.png

Prediction of sRNA and sORF
---------------------------

Based on transcripts, genome annotations and coverage information, sRNAs can be detected. Moreover, we 
have TSSs and processing sites which can be used for detecting UTR-derived sRNAs as well. Now, we can 
get sRNAs by running subcommand ``srna``. Normally, we recommend that the user uses the libraries of RNA-Seq 
with transcript fragmented as well. Here, we only use TEX +/- for testing.

For running ``srna``, we can apply several filters to improve the detection. These filters are ``tss``, ``sec_str``,
``blast_nr``, ``blast_srna``, ``promoter``, ``term``, ``sorf``. Normally, ``tss``, ``sec_str``,
``blast_nr``, ``blast_srna`` are recommended to be used.

Please be aware, filters are strict. For examples, if your filters include ``term``, only the sRNAs which are 
associated with terminators will be included in the list of best candidates. If you want to include terminator information 
but do not use terminator as a filter, you can remove ``term`` in filters and still assign the path of terminator gff file. 
The results will include the sRNAs which are not associated with terminators, and terminator information can be shown 
and checked in the results as well. Many parameters can be used for adjustment the prediction, such as ``blast_e_srna``, 
``blast_e_nr``, ``blast_score_nr``, ``blast_score_srna``, etc. 
For details of the filters, please check the section :ref:`srna` in ANNOgesic subcommand.

Before running ``srna``, we have to download sRNA database (we can use `BSRD <http://www.bac-srna.org/BSRD/index.jsp>`_) and 
`nr database <ftp://ftp.ncbi.nih.gov/blast/db/FASTA/>`_ (if you have not downloaded before). 
We can download fasta file of `BSRD <http://www.bac-srna.org/BSRD/index.jsp>`_ from our 
`Git repository <https://github.com/Sung-Huan/ANNOgesic/tree/master/database>`_.

::

    $ wget -cP ANNOgesic/input/databases/ https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/database/sRNA_database_BSRD.fa

If you have your sRNA database in other folders, please assign your path of databases to ``--srna_database_path`` 
(please check :ref:`srna` to modify the headers of your database).
If your database was formatted before, you can remove ``--srna_format``.
In order to use the recommended filters, we have to download 
`nr database <ftp://ftp.ncbi.nih.gov/blast/db/FASTA/>`_ (takes a while). If you already downloaded it, 
you can skip this step.

::

    $ wget -cP ANNOgesic/input/databases/ ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
    $ gunzip ANNOgesic/input/databases/nr.gz
    $ mv ANNOgesic/input/databases/nr ANNOgesic/input/databases/nr.fa

If your nr database is in other folders, please assign your path to ``--nr_database_path``.
You can also remove ``--nr_format`` if your database is already formatted.
Now, we can use the recommended filters to run ``srna``, but it may takes a while.

::

    $ annogesic srna \
        --filter_info tss blast_srna sec_str blast_nr \
        --annotation_files ANNOgesic/input/references/annotations/NC_009839.1.gff \
        --tss_files ANNOgesic/output/TSSs/gffs/NC_009839.1_TSS.gff \
        --processing_site_files ANNOgesic/output/processing_sites/gffs/NC_009839.1_processing.gff \
        --transcript_files ANNOgesic/output/transcripts/gffs/NC_009839.1_transcript.gff \
        --fasta_files ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        --terminator_files ANNOgesic/output/terminators/gffs/best_candidates/NC_009839.1_term.gff \
        --promoter_tables ANNOgesic/output/promoters/NC_009839.1/MEME/promoter_motifs_NC_009839.1_allgenome_all_types_45_nt/meme.csv \
        --promoter_names MOTIF_1 \
        --mountain_plot \
        --utr_derived_srna \
        --compute_sec_structures \
        --srna_format \
        --nr_format \
        --nr_database_path ANNOgesic/input/databases/nr \
        --srna_database_path ANNOgesic/input/databases/sRNA_database_BSRD \
        --tex_notex_libs $TEX_LIBS \
        --replicate_tex all_1 \
        --project_path ANNOgesic

If you have sORF information, you can also assign the path of sORF gff file to ``--sorf_files``. 
Then, the comparison between sRNAs and sORFs can be executed.

Output files are following.

::

    $ ls ANNOgesic/output/sRNAs/
    blast_results_and_misc  figs  gffs  log.txt  sRNA_2d_NC_009839.1  sRNA_seq_NC_009839.1  statistics  tables

``blast_results_and_misc`` stores the results of blast; ``figs`` stores plots of sRNAs; 
``statistics`` stores statistic files.

``sRNA_2d_NC_009839.1`` and ``sRNA_seq_NC_009839.1`` are text files of sRNA sequences and secondary structures.

::

    $ ls ANNOgesic/output/sRNAs/blast_results_and_misc/
    nr_blast_NC_009839.1.txt  sRNA_blast_NC_009839.1.txt
    $ ls ANNOgesic/output/sRNAs/figs/
    dot_plots  mountain_plots  sec_plots
    $ ls ANNOgesic/output/sRNAs/figs/mountain_plots/NC_009839.1/
    srna0_NC_009839.1_36954_37044_-_mountain.pdf     srna25_NC_009839.1_854600_854673_-_mountain.pdf    srna40_NC_009839.1_1091155_1091251_-_mountain.pdf  srna56_NC_009839.1_1440826_1441414_+_mountain.pdf
    srna10_NC_009839.1_248098_248257_-_mountain.pdf  srna26_NC_009839.1_879881_880088_-_mountain.pdf    srna41_NC_009839.1_1097654_1097750_-_mountain.pdf  srna57_NC_009839.1_1448211_1448306_+_mountain.pdf
    ...

    $ ls ANNOgesic/output/sRNAs/figs/dot_plots/NC_009839.1/
    srna0_NC_009839.1_36954_37044_-_dp.ps     srna25_NC_009839.1_854600_854673_-_dp.ps    srna40_NC_009839.1_1091155_1091251_-_dp.ps  srna56_NC_009839.1_1440826_1441414_+_dp.ps
    srna10_NC_009839.1_248098_248257_-_dp.ps  srna26_NC_009839.1_879881_880088_-_dp.ps    srna41_NC_009839.1_1097654_1097750_-_dp.ps  srna57_NC_009839.1_1448211_1448306_+_dp.ps
    ...

    $ ls ANNOgesic/output/sRNAs/figs/sec_plots/NC_009839.1/
    rna0_NC_009839.1_36954_37044_-_rss.ps     srna25_NC_009839.1_854600_854673_-_rss.ps    srna40_NC_009839.1_1091155_1091251_-_rss.ps  srna56_NC_009839.1_1440826_1441414_+_rss.ps
    srna10_NC_009839.1_248098_248257_-_rss.ps  srna26_NC_009839.1_879881_880088_-_rss.ps    srna41_NC_009839.1_1097654_1097750_-_rss.ps  srna57_NC_009839.1_1448211_1448306_+_rss.ps
    ...

    $ ls ANNOgesic/output/sRNAs/statistics/
    stat_NC_009839.1_sRNA_blast.csv  stat_sRNA_class_NC_009839.1.csv

In ``gffs`` and ``tables``, three different folders are generated. ``all_candidates`` is for all candidates 
without filtering; ``best_candidates`` is for the candidates after filtering; 
``for_classes`` is for different sRNA types based on ``stat_sRNA_class_NC_009839.1.csv``. 

::

    $ ls ANNOgesic/output/sRNAs/gffs/
    all_candidates  best_candidates  for_classes
    $ ls ANNOgesic/output/sRNAs/tables/
    all_candidates  best_candidates  for_classes
    $ ls ANNOgesic/output/sRNAs/gffs/all_candidates/
    NC_009839.1_sRNA.gff
    $ ls ANNOgesic/output/sRNAs/tables/all_candidates/
    NC_009839.1_sRNA.csv
    $ ls ANNOgesic/output/sRNAs/gffs/best_candidates/
    NC_009839.1_sRNA.gff
    $ ls ANNOgesic/output/sRNAs/tables/best_candidates/
    NC_009839.1_sRNA.csv
    $ ls ANNOgesic/output/sRNAs/gffs/for_classes/NC_009839.1/
    class_1_all.gff                                          class_1_class_2_class_7_all.gff                  class_2_all.gff                                  class_3_all.gff
    class_1_class_2_all.gff                                  class_1_class_3_all.gff                          class_2_class_3_all.gff                          class_3_class_4_all.gff
    ...

    $ ls ANNOgesic/output/sRNAs/tables/for_classes/NC_009839.1/
    class_1_all.csv                                          class_1_class_2_class_7_all.csv                  class_2_all.csv                                  class_3_all.csv
    class_1_class_2_all.csv                                  class_1_class_3_all.csv                          class_2_class_3_all.csv                          class_3_class_4_all.csv
    ...

If the amount of sRNA candidates are too many for the user, please check the :ref:`FAQ Q9` to do the further filtering.
As we know, expressed regions without annotations may be sORFs as well. 
In order to get information of sORFs, we can use subcommand ``sorf``.

::

    $ annogesic sorf \
        --annotation_files ANNOgesic/input/references/annotations/NC_009839.1.gff \
        --tss_files ANNOgesic/output/TSSs/gffs/NC_009839.1_TSS.gff \
        --transcript_files ANNOgesic/output/transcripts/gffs/NC_009839.1_transcript.gff \
        --fasta_files ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        --srna_files ANNOgesic/output/sRNAs/gffs/best_candidates/NC_009839.1_sRNA.gff \
        --tex_notex_libs $TEX_LIBS \
        --replicate_tex all_1 -u \
        --project_path ANNOgesic

For generating best candidates, some filters can be assigned 
(ex: with ribosome binding site (Shine-Dalgarno sequence), with TSS, without overlap with sRNA, etc.).
After running ``sorf``, gff files, statistic files and tables of the sORF will be generated. ``all_candidates`` 
stores the gff files and tables without filtering; ``best_candidates`` stores the gff_files and tables with filtering.

::

    $ ls ANNOgesic/output/sORFs/gffs/all_candidates/
    NC_009839.1_sORF.gff
    $ ls ANNOgesic/output/sORFs/gffs/best_candidates/
    NC_009839.1_sORF.gff
    $ ls ANNOgesic/output/sORFs/tables/all_candidates/
    NC_009839.1_sORF.csv
    $ ls ANNOgesic/output/sORFs/tables/best_candidates/
    NC_009839.1_sORF.csv
    $ ls ANNOgesic/output/sORFs/statistics/
    stat_NC_009839.1_sORF.csv

Performing sRNA target prediction
---------------------------------

Now we have sRNA candidates. If we want to know targets of these sRNAs, we can use ``srna_target``.

::

    $annogesic srna_target \
        --annotation_files ANNOgesic/input/references/annotations/NC_009839.1.gff \
        --fasta_files ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        --srna_files ANNOgesic/output/sRNAs/gffs/best_candidates/NC_009839.1_sRNA.gff \
        --query_srnas NC_009839.1:36954:37044:- \
        --mode_intarna H \
        --program RNAup IntaRNA RNAplex \
        --project_path ANNOgesic

For testing, we only assign one sRNA to do the prediction. You can also assign several sRNAs like 
``NC_009839.1:36954:37044:- NC_009839.1:75845:75990:+``. If you want to compute all sRNAs, you 
can assign ``all`` to ``--query_srnas`` (may take several days).

Several output folders will be generated. 

::

    $ ls ANNOgesic/output/sRNA_targets/
    IntaRNA_results  merged_results  log.txt  RNAplex_results  RNAup_results  sRNA_seqs  target_seqs

``sRNA_seqs`` and ``target_seqs`` are for sequences of the sRNAs and the potential targets.

::

    $ ls ANNOgesic/output/sRNA_targets/sRNA_seqs
    NC_009839.1_sRNA.fa
    $ ls ANNOgesic/output/sRNA_targets/target_seqs
    NC_009839.1_target.fa

``IntaRNA_results``, ``RNAplex_results`` and ``RNAup_results`` are for the output of 
`IntaRNA <https://github.com/BackofenLab/IntaRNA/>`_, `RNAplex <http://www.tbi.univie.ac.at/RNA/RNAplex.1.html>`_ and
`RNAup <http://www.tbi.univie.ac.at/RNA/RNAup.1.html>`_.

::

    $ ls ANNOgesic/output/sRNA_targets/RNAplex_results/NC_009839.1/
    NC_009839.1_RNAplex_rank.csv  NC_009839.1_RNAplex.txt
    $ ls ANNOgesic/output/sRNA_targets/RNAup_results/NC_009839.1/
    NC_009839.1_RNAup.log  NC_009839.1_RNAup_rank.csv  NC_009839.1_RNAup.txt
    $ ls ANNOgesic/output/sRNA_targets/IntaRNA_results/NC_009839.1/
    NC_009839.1_IntaRNA_rank.csv  NC_009839.1_IntaRNA.txt

``merged_results`` is for the merged results of `IntaRNA <https://github.com/BackofenLab/IntaRNA/>`_, 
`RNAplex <http://www.tbi.univie.ac.at/RNA/RNAplex.1.html>`_ and
`RNAup <http://www.tbi.univie.ac.at/RNA/RNAup.1.html>`_. ``NC_009839.1_merge.csv``  contains all results of the 
assigned methods, and ``NC_009839.1_overlap.csv`` only stores candidates which are top 50 (default) in the assigned methods.

::

    $ ls ANNOgesic/output/sRNA_targets/merged_results/NC_009839.1/
    NC_009839.1_merge.csv  NC_009839.1_overlap.csv

Mapping and detecting of circular RNA
-------------------------------------

You may also be interested in circular RNAs. The subcommand ``circrna`` can help us to predict circular RNAs by  
using `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_. Since 
we didn't map reads before, we can also do mapping by running ``circrna``. If you already mapped 
the reads by `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_ with ``--splits``, you can 
add path of the bam files to ``--bam_files`` directly. However, 
if you mapped the reads by other tools or you mapped the reads by 
`Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_ without ``--splits``, Unfortunately, 
you have to re-map the reads(``--read_files``) again. You can assign the number of parallels (``--parallels``) for mapping.

Since we just want to test the subcommand. 
Thus, we can reduce the running time by selecting the subset of the reads (first 50000) for only testing.

::

     $ head -n 50000 ANNOgesic/input/reads/SRR515254.fasta > ANNOgesic/input/reads/SRR515254_50000.fasta
     $ head -n 50000 ANNOgesic/input/reads/SRR515255.fasta > ANNOgesic/input/reads/SRR515255_50000.fasta
     $ head -n 50000 ANNOgesic/input/reads/SRR515256.fasta > ANNOgesic/input/reads/SRR515256_50000.fasta
     $ head -n 50000 ANNOgesic/input/reads/SRR515257.fasta > ANNOgesic/input/reads/SRR515257_50000.fasta
     $ rm ANNOgesic/input/reads/SRR515254.fasta
     $ rm ANNOgesic/input/reads/SRR515255.fasta
     $ rm ANNOgesic/input/reads/SRR515256.fasta
     $ rm ANNOgesic/input/reads/SRR515257.fasta

Then we setup the read files.

::
    $ READ_FILES=ANNOgesic/input/reads/SRR515254_50000.fasta,\
    ANNOgesic/input/reads/SRR515255_50000.fasta,\
    ANNOgesic/input/reads/SRR515256_50000.fasta,\
    ANNOgesic/input/reads/SRR515257_50000.fasta


After that, we assign ``all_samples:$READ_FILE`` to ``--read_files``. ``all_sample`` is the set name of read files. 
The all four read files will be compute together. Now, we can try ``circrna``

::

     $ annogesic circrna \
         --fasta_files ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
         --parallels 10 \
         --annotation_files ANNOgesic/input/references/annotations/NC_009839.1.gff \
         --read_files all_samples:$READ_FILES \
         --project_path ANNOgesic

``testrealign.x`` is not available, please refer to :ref:`Required tools or databases`.

Several output folders will be generated.

::

    $ ls ANNOgesic/output/circRNAs/
    circRNA_tables  gffs  log.txt  segemehl_alignment_files  segemehl_splice_results  statistics

``segemehl_alignment_files`` and ``segemehl_splice_results`` are for output of 
`Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_. ``segemehl_alignment_files`` stores Bam files of 
the alignment and ``segemehl_splice_results`` stores results of the splice detection.

::

    $ ls ANNOgesic/output/circRNAs/segemehl_alignment_files/NC_009839.1/
    SRR515254_50000_NC_009839.1.bam  SRR515256_50000_NC_009839.1.bam
    SRR515255_50000_NC_009839.1.bam  SRR515257_50000_NC_009839.1.bam
    $ ls ANNOgesic/output/circRNAs/segemehl_splice_results/NC_009839.1/
    NC_009839.1_all_samples_splicesites.bed  NC_009839.1_all_samples_transrealigned.bed

Gff files, tables and statistic files are stored in ``gffs``, ``circRNA_tables`` and ``statistics``.

::

    $ ls ANNOgesic/output/circRNAs/gffs/NC_009839.1/
    NC_009839.1_all_samples_circRNA_all.gff  NC_009839.1_all_samples_circRNA_best.gff
    $ ls ANNOgesic/output/circRNAs/circRNA_tables/NC_009839.1/
    NC_009839.1_all_samples_circRNA_all.csv  NC_009839.1_all_samples_circRNA_best.csv
    $ ls ANNOgesic/output/circRNAs/statistics/
    stat_NC_009839.1_all_samples_circRNA.csv

``NC_009839.1_all_samples_circRNA_all.gff`` and ``NC_009839.1_all_samples_circRNA_all.csv`` store all circular RNAs without filtering. 
``NC_009839.1_all_samples_circRNA_best.gff`` and ``NC_009839.1_all_samples_circRNA_best.csv`` store
the circular RNAs after filtering. In our case, there are some circular RNAs can be detected without filtering, but no one 
can exist after filtering.

SNP calling
--------------

If we want to know SNPs or mutations based on our RNA-Seq data, we can use ``snp`` to achieve this purpose.
``snp`` is compose of two parts. One part is for obtaining the differences between our reference genome 
and the closely related genome. If we have no fasta file of our reference genome, 
this part will be very useful. We just need to map the reads on the fasta file of the closely related genome. Then 
using ``snp`` can automatically detect differences between the closely related genome and our reference genome. 
Furthermore, potential fasta files of the refernce genome can be generated automatically as well. 
The other part is for detecting SNPs or mutations of the reference genome if the fasta file of the reference genome can be provided.
In this part, you can know real mutations of our reference genonme.

Before running the subcommand, BAM files are required. Since we already generated them via 
running ``circrna``, we can just put them to the corresponding folder. Please remember that the mapping function of 
``circrna`` is only basic one.

First, we copy the bam files to ``BAMs_map_reference_genomes``.

::

    $ cp ANNOgesic/output/circRNAs/segemehl_alignment_files/NC_009839.1/SRR51525* ANNOgesic/input/BAMs/BAMs_map_reference_genomes/tex_notex

Now, we can set our bam files

::
    $ BAM_FILES=ANNOgesic/input/BAMs/BAMs_map_reference_genomes/tex_notex/SRR515254_50000_NC_009839.1.bam,\
      ANNOgesic/input/BAMs/BAMs_map_reference_genomes/tex_notex/SRR515255_50000_NC_009839.1.bam,\
      ANNOgesic/input/BAMs/BAMs_map_reference_genomes/tex_notex/SRR515256_50000_NC_009839.1.bam,\
      ANNOgesic/input/BAMs/BAMs_map_reference_genomes/tex_notex/SRR515257_50000_NC_009839.1.bam

Then we can run the subcommand with three programs -- ``extend_BAQ``, ``with_BAQ`` and ``without_BAQ``. 
``all_sample:$BAM_FILES`` for ``--bam_files`` means the set name of BAM files is "all_sample", 
and all four BAM files need to be computed together.

::

    $ annogesic snp \
        --bam_type reference_genome \
        --program with_BAQ without_BAQ extend_BAQ \
        --bam_files all_samples:$BAM_FILES \
        --fasta_files ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        --project_path ANNOgesic

Two output folders will be generated, ``compare_related_and_reference_genomes`` is for the results of comparison between closely related genome 
and reference genome, ``mutations_of_reference_genomes`` is for results of detecting mutations of the reference genome.

::

    $ ls ANNOgesic/output/SNP_calling/                                                                                                      
    compare_related_and_reference_genomes  mutations_of_reference_genomes  log.txt

Since we run ``reference_genome``,  the output folders are generated under ``mutations_of_reference_genomes``.

::

    $ ls ANNOgesic/output/SNP_calling/mutations_of_reference_genomes/
    seqs  SNP_raw_outputs  SNP_tables  statistics

The output folders are compose of three parts - ``extend_BAQ``, ``with_BAQ`` and ``without_BAQ``.

::

    $ ls ANNOgesic/output/SNP_calling/mutations_of_reference_genomes/seqs/
    extend_BAQ/  with_BAQ/    without_BAQ/

In ``seqs``, the potential sequences can be found.

::

    $ ls ANNOgesic/output/SNP_calling/mutations_of_reference_genomes/seqs/with_BAQ/NC_009839.1/
    NC_009839.1_all_samples_NC_009839.1_1_1.fa

``SNP_raw_outputs`` stores output of `Samtools and Bcftools <https://github.com/samtools>`_. 
``SNP_tables`` stores results after filtering and the indices of potential sequence 
(potential sequences are stored in ``seqs``).
``statistics`` stores the statistic files.

::

    $ ls ANNOgesic/output/SNP_calling/mutations_of_reference_genomes/SNP_raw_outputs/NC_009839.1/
    NC_009839.1_extend_BAQ_all_samples.vcf  NC_009839.1_with_BAQ_all_samples.vcf  NC_009839.1_without_BAQ_all_samples.vcf
    $ ls ANNOgesic/output/SNP_calling/mutations_of_reference_genomes/SNP_tables/NC_009839.1/
    NC_009839.1_extend_BAQ_all_samples_best.vcf           NC_009839.1_with_BAQ_all_samples_best.vcf           NC_009839.1_without_BAQ_all_samples_best.vcf
    NC_009839.1_extend_BAQ_all_samples_seq_reference.csv  NC_009839.1_with_BAQ_all_samples_seq_reference.csv  NC_009839.1_without_BAQ_all_samples_seq_reference.csv
    $ ls ANNOgesic/output/SNP_calling/mutations_of_reference_genomes/statistics/
    figs                                                  stat_NC_009839.1_with_BAQ_all_samples_SNP_best.csv     stat_NC_009839.1_without_BAQ_all_samples_SNP_raw.csv
    stat_NC_009839.1_extend_BAQ_all_samples_SNP_best.csv  stat_NC_009839.1_with_BAQ_all_samples_SNP_raw.csv
    stat_NC_009839.1_extend_BAQ_all_samples_SNP_raw.csv   stat_NC_009839.1_without_BAQ_all_samples_SNP_best.csv
    $ ls ANNOgesic/output/SNP_calling/mutations_of_reference_genomes/statistics/figs
    NC_009839.1_extend_BAQ_all_samples_NC_009839.1_SNP_QUAL_best.png  NC_009839.1_with_BAQ_all_samples_NC_009839.1_SNP_QUAL_best.png  NC_009839.1_without_BAQ_all_samples_NC_009839.1_SNP_QUAL_best.png
    NC_009839.1_extend_BAQ_all_samples_NC_009839.1_SNP_QUAL_raw.png   NC_009839.1_with_BAQ_all_samples_NC_009839.1_SNP_QUAL_raw.png   NC_009839.1_without_BAQ_all_samples_NC_009839.1_SNP_QUAL_raw.png

Mapping Gene ontology
---------------------

Gene ontology is useful for understanding functions of gene products. 
``go_term`` can search GO terms of the proteins in annotation files. Before running ``go_term``, we 
need to prepare some databases. First, please download 
`goslim.obo <http://geneontology.org/page/go-slim-and-subset-guide>`_ and 
`go.obo <http://geneontology.org/page/download-ontology>`_ and 
`idmapping_selected.tab <http://www.uniprot.org/downloads>`_.

::

    $ wget -cP ANNOgesic/input/databases http://www.geneontology.org/ontology/subsets/goslim_generic.obo
    $ wget -cP ANNOgesic/input/databases http://geneontology.org/ontology/go.obo
    $ wget -cP ANNOgesic/input/databases ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
    $ gunzip ANNOgesic/input/databases/idmapping_selected.tab.gz

Now, we have all required databases. We can also import information of the transcripts to 
generate results which only contain the expressed CDSs.

Let's try it.

::

    $ annogesic go_term \
        --annotation_files ANNOgesic/input/references/annotations/NC_009839.1.gff \
        --transcript_files ANNOgesic/output/transcripts/gffs/NC_009839.1_transcript.gff \
        --go_obo ANNOgesic/input/databases/go.obo \
        --goslim_obo ANNOgesic/input/databases/goslim_generic.obo \
        --uniprot_id ANNOgesic/input/databases/idmapping_selected.tab \
        --project_path ANNOgesic

The results of ``go_term`` are stored in ``GO_term_results``. The statistic files and 
figures are stored in ``statistics``.

::

    $ ls ANNOgesic/output/GO_terms/
    all_CDSs  expressed_CDSs  log.txt
    $ ls ANNOgesic/output/GO_terms/all_CDSs/
    GO_term_results  statistics
    $ ls ANNOgesic/output/GO_terms/all_CDSs/GO_term_results/NC_009839.1/
    all_genomes_uniprot.csv
    $ ls ANNOgesic/output/GO_terms/all_CDSs/statistics/NC_009839.1/
    figs  stat_NC_009839.1.csv
    $ ls ANNOgesic/output/GO_terms/all_CDSs/statistics/NC_009839.1/figs/
    NC_009839.1_biological_process.png  NC_009839.1_cellular_component.png  NC_009839.1_molecular_function.png  NC_009839.1_three_roots.png

Prediction of Subcellular localization
--------------------------------------

Subcellular localization is also a useful information for analysis of protein functions. For 
detecting subcellular localization, we can use the subcommand 
``localization``. We can also import 
information of the transcripts to generate results which only contain the expressed CDSs.

::

    $ annogesic localization \
        --annotation_files ANNOgesic/input/references/annotations/NC_009839.1.gff \
        --fasta_files ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        --transcript_files ANNOgesic/output/transcripts/gffs/NC_009839.1_transcript.gff \
        --bacteria_type negative \
        --project_path ANNOgesic

Two output folders will be generated. ``psortb_results`` stores the output files 
of `Psortb <http://www.psort.org/psortb/>`_. ``statistics`` stores 
statistic files and figures.

::

    $ ls ANNOgesic/output/subcellular_localization/
    all_CDSs  expressed_CDSs  log.txt
    $ ls ANNOgesic/output/subcellular_localization/all_CDSs/
    psortb_results  statistics
    $ ls ANNOgesic/output/subcellular_localization/all_CDSs/psortb_results/NC_009839.1/
    NC_009839.1_raw.txt  NC_009839.1_table.csv
    $ ls ANNOgesic/output/subcellular_localization/all_CDSs/statistics/NC_009839.1/
    NC_009839.1_NC_009839.1_sublocal.png  stat_NC_009839.1_sublocal.csv

Generating protein-protein interaction network
----------------------------------------------

``ppi_network`` can detect protein-protein interaction based on `STRING <http://string-db.org/>`_ 
(a database of protein-protein interaction) and searching the literatures by implementing 
`PIE <http://www.ncbi.nlm.nih.gov/CBBresearch/Wilbur/IRET/PIE/>`_ 
(text-mining for protein-protein interaction). Therefore, ``ppi_network`` can generate protein-protein 
interaction networks based on literatures.

Before running the subcommand, you need to download 
`species.v{$VERSIO}.txt from STRING <http://string-db.org/cgi/download.pl>`_

::

    $ wget -cP ANNOgesic/input/databases http://string-db.org/newstring_download/species.v10.5.txt

Now, we can try the subcommand.

::

    $ annogesic ppi_network \
        --query_strains NC_009839.1.gff:NC_009839.1:'Campylobacter jejuni 81176':'Campylobacter jejuni' \
        --annotation_files ANNOgesic/input/references/annotations/NC_009839.1.gff \
        --species_string ANNOgesic/input/databases/species.v10.5.txt \
        --query NC_009839.1:962231:963001:- NC_009839.1:123943:125151:+ \
        --without_strain_pubmed \
        --project_path ANNOgesic

We only detected for two proteins. If you want to detect for all proteins in gff files, 
you can easily assign ``all`` in ``--query``.

Three output folders were generated.

::

    $ ls ANNOgesic/output/PPI_networks/
    all_results/  best_results/ figures/  log.txt

``all_results`` is for all interactions without filtering. ``best_results`` is for the interactions with 
the high `PIE <http://www.ncbi.nlm.nih.gov/CBBresearch/Wilbur/IRET/PIE/>`_ score. ``figures`` is for 
figures of the protein-protein interaction networks. There are two subfolders - ``with_strain`` and ``without_strain`` in 
``figures``, ``all_results``, and ``best_results``. The two subfolders store all information of the interactions and PIE scores. 
``with_strain`` is for results with assigning specific strain name for searching literatures. 
``without_strain`` is for results without giving specific strain name for searching literatures.

::

    $ ls ANNOgesic/output/PPI_networks/all_results/PPI_NC_009839.1/
    NC_009839.1_without_strain.csv  NC_009839.1_with_strain.csv  without_strain  with_strain
    $ ls ANNOgesic/output/PPI_networks/best_results/PPI_NC_009839.1/
    NC_009839.1_without_strain.csv  NC_009839.1_with_strain.csv  without_strain  with_strain
    $ ls ANNOgesic/output/PPI_networks/figures/PPI_NC_009839.1/
    without_strain  with_strain
    $ ls ANNOgesic/output/PPI_networks/all_results/PPI_NC_009839.1/with_strain/NC_009839.1/
    atpC_atpD.csv                     Cjejjejuni_010100005380_livH.csv                     Cjejjejuni_010100005385_livF.csv  Cjejjejuni_010100005385_livM.csv  livH_livG.csv  livM_livG.csv
    Cjejjejuni_010100005380_livF.csv  Cjejjejuni_010100005380_livM.csv                     Cjejjejuni_010100005385_livG.csv  livG_livF.csv                     livH_livM.csv
    Cjejjejuni_010100005380_livG.csv  Cjejjejuni_010100005385_Cjejjejuni_010100005380.csv  Cjejjejuni_010100005385_livH.csv  livH_livF.csv                     livM_livF.csv
    $ ls ANNOgesic/output/PPI_networks/all_results/PPI_NC_009839.1/without_strain/NC_009839.1/
    atpC_atpD.csv                     Cjejjejuni_010100005380_livH.csv                     Cjejjejuni_010100005385_livF.csv  Cjejjejuni_010100005385_livM.csv  livH_livG.csv  livM_livG.csv
    Cjejjejuni_010100005380_livF.csv  Cjejjejuni_010100005380_livM.csv                     Cjejjejuni_010100005385_livG.csv  livG_livF.csv                     livH_livM.csv
    Cjejjejuni_010100005380_livG.csv  Cjejjejuni_010100005385_Cjejjejuni_010100005380.csv  Cjejjejuni_010100005385_livH.csv  livH_livF.csv                     livM_livF.csv
    $ ls ANNOgesic/output/PPI_networks/best_results/PPI_NC_009839.1/without_strain/NC_009839.1/
    Cjejjejuni_010100005380_livF.csv  Cjejjejuni_010100005380_livH.csv  Cjejjejuni_010100005385_livF.csv  Cjejjejuni_010100005385_livH.csv  livG_livF.csv  livH_livG.csv  livM_livF.csv
    Cjejjejuni_010100005380_livG.csv  Cjejjejuni_010100005380_livM.csv  Cjejjejuni_010100005385_livG.csv  Cjejjejuni_010100005385_livM.csv  livH_livF.csv  livH_livM.csv  livM_livG.csv
    $ ls ANNOgesic/output/PPI_networks/best_results/PPI_NC_009839.1/with_strain/NC_009839.1/
    Cjejjejuni_010100005385_Cjejjejuni_010100005380.csv
    $ ls ANNOgesic/output/PPI_networks/figures/PPI_NC_009839.1/with_strain/NC_009839.1/
    C8J_RS04960_livG.png
    $ ls ANNOgesic/output/PPI_networks/figures/PPI_NC_009839.1/without_strain/NC_009839.1/
    C8J_RS04960_livG.png

Generating riboswitch and RNA thermometer
-----------------------------------------

If we want to detect riboswitches and RNA thermometers, we can use subcommand ``riboswitch_thermometer``.
Before running it, we need to get the information of known riboswitches and RNA thermometers in Rfam. 
The `riboswitches and RNA thermometer files <https://github.com/Sung-Huan/ANNOgesic/tree/master/database>`_ 
can be downloaded them from our Git repository.

::

    $ wget -cP ANNOgesic/input/riboswitch_ID_file/ https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/database/Rfam_riboswitch_ID.csv
    $ wget -cP ANNOgesic/input/RNA_thermometer_ID_file/ https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/database/Rfam_RNA_thermometer_ID.csv

We also need to download `Rfam <http://rfam.xfam.org/>`_.

::

    $ wget -cP ANNOgesic/input/databases ftp://ftp.ebi.ac.uk/pub/databases/Rfam/12.0/Rfam.tar.gz
    $ cd ANNOgesic/input/databases
    $ tar -zxvf Rfam.tar.gz
    $ rm Rfam.tar.gz
    $ cd ../../../

Now we can try the subcommand.

::

    $ annogesic riboswitch_thermometer \
        --annotation_files ANNOgesic/input/references/annotations/NC_009839.1.gff \
        --fasta_files ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        --riboswitch_id_file ANNOgesic/input/riboswitch_ID_file/Rfam_riboswitch_ID.csv \
        --rna_thermometer_id_file ANNOgesic/input/RNA_thermometer_ID_file/Rfam_RNA_thermometer_ID.csv \
        --rfam_path ANNOgesic/input/databases/CMs/Rfam.cm \
        --transcript_files ANNOgesic/output/transcripts/gffs/NC_009839.1_transcript.gff \
        --tss_files ANNOgesic/output/TSSs/gffs/NC_009839.1_TSS.gff \
        --project_path ANNOgesic

Output files are following, ``gffs`` stores gff files of the riboswitchs / RNA_thermometers; 
``tables`` stores tables of the riboswitchs / RNA_thermometers; 
``scan_Rfam_results`` stores output files of scanning Rfam; ``statistics`` is for statistic files.

::

     $ ls ANNOgesic/output/riboswitches/
     gffs  log.txt  scan_Rfam_results  statistics  tables
     $ ls ANNOgesic/output/riboswitches/gffs/
     NC_009839.1_riboswitch.gff
     $ ls ANNOgesic/output/riboswitches/scan_Rfam_results/NC_009839.1/
     NC_009839.1_riboswitch_prescan.txt  NC_009839.1_riboswitch_scan.txt
     $ ls ANNOgesic/output/riboswitches/tables/
     NC_009839.1_riboswitch.csv
     $ ls ANNOgesic/output/riboswitches/statistics/
     stat_NC_009839.1_riboswitch.txt
     $ ls ANNOgesic/output/RNA_thermometers/
     gffs  log.txt  scan_Rfam_results  statistics  tables
     $ ls ANNOgesic/output/RNA_thermometers/gffs/
     NC_009839.1_RNA_thermometer.gff
     $ ls ANNOgesic/output/RNA_thermometers/scan_Rfam_results/NC_009839.1/
     NC_009839.1_RNA_thermometer_prescan.txt  NC_009839.1_RNA_thermometer_scan.txt
     $ ls ANNOgesic/output/RNA_thermometers/tables/
     NC_009839.1_RNA_thermometer.csv
     $ ls ANNOgesic/output/RNA_thermometers/statistics/
     stat_NC_009839.1_RNA_thermometer.txt

Detection of CRISPR
-------------------
CRISPR is an important features for research of immunology. ``crispr`` is a useful subcommand for CRISPR detection. 
Let's try it.

::

     $ annogesic crispr \
        --annotation_files ANNOgesic/input/references/annotations/NC_009839.1.gff \
        --fasta_files ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        --project_path ANNOgesic

Output are as following, ``CRT_results`` stores output of `CRT <http://www.room220.com/crt/>`_; 
``gffs`` stores gff files of the CRISPRs; ``statistics`` is for statistic files.

::

     $ ls ANNOgesic/output/crisprs/
     CRT_results  gffs  log.txt  statistics
     $ ls ANNOgesic/output/crisprs/CRT_results
     NC_009839.1.txt
     $ ls ANNOgesic/output/crisprs/gffs
     all_candidates  best_candidates
     $ ls ANNOgesic/output/crisprs/gffs/all_candidates
     NC_009839.1_CRISPR.gff
     $ ls ANNOgesic/output/crisprs/gffs/best_candidates
     NC_009839.1_CRISPR.gff
     $ ls ANNOgesic/output/crisprs/statistics
     NC_009839.1.csv

Merge all features to be one gff file
-------------------------------------

Now, we generated all features that ANNOgesic can provide. Sometimes, merging all features into 
one gff file is useful. ``merge_features`` is the subcommand to achieve this purpose. 
Moreover, ``merge_features`` can search parent transcript to each feature that 
we assigned.

Now let's do it. We merge all features that we have.

::

    ALL_FEATURES="ANNOgesic/output/TSSs/gffs/NC_009839.1_TSS.gff \
                  ANNOgesic/input/references/annotations/NC_009839.1.gff \
                  ANNOgesic/output/UTRs/5UTRs/gffs/NC_009839.1_5UTR.gff \
                  ANNOgesic/output/UTRs/3UTRs/gffs/NC_009839.1_3UTR.gff \
                  ANNOgesic/output/terminators/gffs/best_candidates/NC_009839.1_term.gff \
                  ANNOgesic/output/processing_sites/gffs/NC_009839.1_processing.gff \
                  ANNOgesic/output/sRNAs/gffs/best_candidates/NC_009839.1_sRNA.gff \
                  ANNOgesic/output/sORFs/gffs/best_candidates/NC_009839.1_sORF.gff \
                  ANNOgesic/output/riboswitches/gffs/NC_009839.1_riboswitch.gff \
                  ANNOgesic/output/RNA_thermometers/gffs/NC_009839.1_RNA_thermometer.gff \
                  ANNOgesic/output/crisprs/gffs/best_candidates/NC_009839.1_CRISPR.gff"

::

    $ annogesic merge_features \
       --transcript_file ANNOgesic/output/transcripts/gffs/NC_009839.1_transcript.gff \
       --other_features_files $ALL_FEATURES \
       --output_prefix NC_009839.1 \
       --source_for_overlapping RefSeq \
       --project_path ANNOgesic

In the tutorial, if duplicated features exist, only the data from RefSeq will be kept. 
Output gff file is stored in ``merge_all_features``

::

    $ ls ANNOgesic/output/merge_all_features/
    NC_009839.1_merge_features.gff  log.txt

Producing the screenshots
-------------------------

It is a good idea if we can get screenshots of our interesting features. Then we can 
check them very quickly. Therefore, ANNOgesic provides a subcommand ``screenshot`` for 
generating screenshots.

Before we running it, we have to install `IGV <https://www.broadinstitute.org/software/igv/home>`_.

For testing, we use TSSs as main feature, sRNAs and CDSs as side features.

::

    $ annogesic screenshot \
        --main_gff ANNOgesic/output/TSSs/gffs/NC_009839.1_TSS.gff \
        --side_gffs ANNOgesic/input/references/annotations/NC_009839.1.gff \
                    ANNOgesic/output/sRNAs/gffs/best_candidates/NC_009839.1_sRNA.gff \
        --fasta_file ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        --output_folder ANNOgesic/output/TSSs \
        --tex_notex_libs $TEX_LIBS \
        --project_path ANNOgesic

Two txt files and two folders will be generated.

::

    $ ls ANNOgesic/output/TSSs/screenshots/NC_009839.1/
    forward/     forward.txt  reverse/     reverse.txt

``forward.txt`` and ``reverse.txt`` are batch files for running in `IGV <https://www.broadinstitute.org/software/igv/home>`_.
``forward`` and ``reverse`` are the folders for storing screenshots.

Since there are numerous candidates, we only generate several ones for testing.

::

    head -n 30 ANNOgesic/output/TSSs/screenshots/NC_009839.1/forward.txt > ANNOgesic/output/TSSs/screenshots/NC_009839.1/forward_6_cases.txt
    head -n 30 ANNOgesic/output/TSSs/screenshots/NC_009839.1/reverse.txt > ANNOgesic/output/TSSs/screenshots/NC_009839.1/reverse_6_cases.txt


Now, please open `IGV <https://www.broadinstitute.org/software/igv/home>`_ and follow the procedures: Tools -> 
Run Batch Script -> choose ``forward_6_cases.txt``. Once it is done, please do it again for the reverse strand: Tools ->
Run Batch Script -> choose ``reverse_6_cases.txt``. If you want to generate the screenshots for all candidates, 
you can run ``forward.txt`` and ``reverse.txt``. Please be careful, if you use Docker container, the path may be not correct.

As soon as the generation of the screenshots is done, 
we can see that there are several screenshots in ``forward`` and ``reverse``.

::

    $ ls ANNOgesic/output/TSSs/screenshots/NC_009839.1/forward
    NC_009839.1:1396-1396.png  NC_009839.1:14812-14812.png  NC_009839.1:6676-6676.png  NC_009839.1:6680-6680.png  NC_009839.1:8098-8098.png  NC_009839.1:9295-9295.png
    $ ls ANNOgesic/output/TSSs/screenshots/NC_009839.1/reverse
    NC_009839.1:15670-15670.png  NC_009839.1:18053-18053.png  NC_009839.1:18360-18360.png  NC_009839.1:2199-2199.png  NC_009839.1:4463-4463.png  NC_009839.1:856-856.png

Coloring the screenshots
------------------------

If we have numerous libraries and we want to check TSSs, distinguishing the 
tracks of TEX+ and TEX- will be painful. Therefore, we provide a subcommand ``colorize_screenshot_tracks`` to colorize 
our screenshots based on the tracks.

::

    $ annogesic colorize_screenshot_tracks \
        --track_number 2 \
        --screenshot_folder ANNOgesic/output/TSSs \
        --project_path ANNOgesic

The output filenames are the same as before. However, when we open the files of figures, the tracks are colorized.

::

    $ ls ANNOgesic/output/TSSs/screenshots/NC_009839.1/forward
    NC_009839.1:1396-1396.png  NC_009839.1:14812-14812.png  NC_009839.1:6676-6676.png  NC_009839.1:6680-6680.png  NC_009839.1:8098-8098.png  NC_009839.1:9295-9295.png
    $ ls ANNOgesic/output/TSSs/screenshots/NC_009839.1/reverse
    NC_009839.1:15670-15670.png  NC_009839.1:18053-18053.png  NC_009839.1:18360-18360.png  NC_009839.1:2199-2199.png  NC_009839.1:4463-4463.png  NC_009839.1:856-856.png

Now we already finished the first wonderful trip of ANNOgesic. Hopefully, you enjoy it!
