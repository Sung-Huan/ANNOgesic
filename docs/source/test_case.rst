Tutorial of ANNOgesic
=====================

This tutorial through a small test case to show you how to use
ANNOgesic's subcommands. It build on a public differential RNA-Seq
data set from **Campylobacter jejuni** subsp. jejuni 81116 that can be
downloaded from NCBI GEO and that was part of a work by `Dugar et
al. <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38883>`_.
There will be several output files which are generated in different formats - 
The CSV (tabular separated plain text files) files (opened by LibreOffice or Excel), GFF3 files, TXT files and figures. 
For viewing GFF3 file, you can use a genome browser, for example `IGB <http://bioviz.org/igb/index.html>`_, 
`IGV <https://www.broadinstitute.org/igv/>`_ or `Jbrowse <http://jbrowse.org/>`_.

Before we start, please check :ref:`The format of filename` and 
:ref:`The input format of libraries for running ANNOgesic` in 
the section :ref:`ANNOgesic's subcommands`. All the details are also in :ref:`ANNOgesic's subcommands`. 
Moreover, all the requirements are listed in the section :ref:`Required tools or databases`.
For command lines which we will present later, please check the 
`run.sh <https://github.com/Sung-Huan/ANNOgesic/tree/master/tutorial_data>`_ in our Git repository.

If the subcommand integrates third-party software, ex: TSSpredator,
please check path of the execute file. If necessary, please assign it properly. Moreover, 
if your execute path is not in your ``$PATH``, please specify your execute path of ANNOgesic as well.

Generating a project
--------------------

First of all, we need to create a working folder by running ``create``. ``-pj`` represents project path.

::

    $ annogesic create -pj ANNOgesic

Then we will see 

::

    Created folder "ANNOgesic" and required subfolders.
    $ ls 
    ANNOgesic

Retrieving the input data
-------------------------

For our test case, the input data can be downloaded from 
`NCBI <ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Campylobacter_jejuni/latest_assembly_versions/GCF_000017905.1_ASM1790v1/>`_.
Let's setting the ``$FTP_SOURCE`` first

::

    FTP_SOURCE=ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Campylobacter_jejuni/latest_assembly_versions/GCF_000017905.1_ASM1790v1/

Then, we can download fasta files(``-f``), gff files(``-g``), gbk files(``-k``), ptt files(``-p``), 
rnt files(``-r``), and converts to embl(``-e``).

::

    $ annogesic get_input_files \
        -F $FTP_SOURCE \
        -g -f -e -k -p -r \
        -pj ANNOgesic

Then we will get following results

::

    $ ls ANNOgesic/input/references/fasta_files/
    NC_009839.1.fa
    $ ls ANNOgesic/input/references/annotations/
    NC_009839.1.embl  NC_009839.1.gbk  NC_009839.1.gff

In fact, these fasta and gff files are exactly what we want to use for the test case.
But, in order to testing ``update_genome_fasta`` and ``annotation_transfer``, we used these files to 
generate some dummy genomes..

Putting wig, and reads to proper location
---------------------------------------------------
For the test case, we can download reads from 
`here <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38883>`_.

::

    $ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP013/SRP013869/SRR515254/SRR515254.sra
    $ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP013/SRP013869/SRR515255/SRR515255.sra
    $ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP013/SRP013869/SRR515256/SRR515256.sra
    $ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP013/SRP013869/SRR515257/SRR515257.sra

Then we need to convert SRA files to Fasta or Fastq format for mapping. We can 
use `SRA toolkit <http://www.ncbi.nlm.nih.gov/books/NBK158900/>`_ for that.

::
  
   $ wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.5.2/sratoolkit.2.5.2-ubuntu64.tar.gz
   $ tar -zxvf sratoolkit.2.5.2-ubuntu64.tar.gz
   $ rm sratoolkit.2.5.2-ubuntu64.tar.gz
   $ ./sratoolkit.2.5.2-ubuntu64/bin/fastq-dump.2.5.2 --fasta SRR515254.sra
   $ ./sratoolkit.2.5.2-ubuntu64/bin/fastq-dump.2.5.2 --fasta SRR515255.sra
   $ ./sratoolkit.2.5.2-ubuntu64/bin/fastq-dump.2.5.2 --fasta SRR515256.sra
   $ ./sratoolkit.2.5.2-ubuntu64/bin/fastq-dump.2.5.2 --fasta SRR515257.sra
   $ mv *.fasta ANNOgesic/input/reads
   $ rm SRR515254.sra SRR515255.sra SRR515256.sra SRR515257.sra

Now we get the reads. Then we have to download the wiggle files.

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

Our genome fasta file is NC_009839.1.fa. Thus the fasta filename in wiggle file should be NC_009839.1 not NC_009839. 
Thus, we need to change it. We can use `replace_seq_id.py <https://github.com/Sung-Huan/ANNOgesic/tree/master/tutorial_data>`_ from our 
Git repository to replace the genome name properly. If the genome names in your fasta, annotation, 
wiggle files are the same, you don't need to do this step.

::

   $ wget https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/tutorial_data/replace_seq_id.py
   $ python3 replace_seq_id.py -i ANNOgesic/input/wigs/tex_notex -n NC_009839.1
   $ rm replace_seq_id.py

We only download one replicate to reduce the running time.

Improving the reference genome
------------------------------

Again, if the data retrieved from NCBI is exactly what you want, you can skip this step and ``annotation_transfer``. 

Although the data that we downloaded before is our real data (``ANNOgesic/input/references``),
we will generate some new dummy files via this step and ``annotation_transfer`` in order to 
show you the function of these subcommands.

Now, we assume that we need to generate fasta file of our query genome. 
First of all, we need to find a close genome (fasta file and gff file can be found) of our query genome. 
Then, we need to generate a mutation table between these two genomes. When these files are produced, 
we can run subcommand ``update_genome_fasta`` for getting fasta file of the target genome. 
For mutation table format, please check the section :ref:`ANNOgesic's subcommands`.

We use a simple example to modify our test case, please check 
`mutation table <https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/tutorial_data/mutation.csv>`_.
Every column of the table is separated by tab. The new genome will be NC_test.1 and test_case2. Therefore, two fasta files 
will be generated in ``ANNOgesic/output/updated_references/fasta_files``.

::

     $ wget -cP ANNOgesic/input/mutation_table https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/tutorial_data/mutation.csv

Now, let's try it

::

     $ annogesic update_genome_fasta \
        -c ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        -o ANNOgesic/output/updated_references/fasta_files/test_case1.fa:NC_test.1 \
           ANNOgesic/output/updated_references/fasta_files/test_case2.fa:test_case2 \
        -m ANNOgesic/input/mutation_table/mutation.csv \
        -pj ANNOgesic

``-r`` is path of the close genome fasta file. In ``-o``, we assign a pairs of output filenames and 
the genomes that we want to put into the file. In our case, "test_case1" is the first output fasta file, and "test_case2" 
is the second output fasta file. "test_case1" stores the sequence of the new genome "NC_test.1", 
and "test_case2" stores the sequence of the other new genome - "test_case2". 

When the running process is done, the following information will appear.

::

    $ Transfering to target fasta
      Please use the new fasta file to remapping again.

Since the data (``ANNOgesic/output/updated_references/fasta_files``) that we generated is not real,
we can ignore the information now. However, if the new fasta file is real query one,
you have to remap again in order to get the correct alignment and coverage files.

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
    $ head ANNOgesic/output/updated_references/fasta_files/test_case1.fa
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

We can see that the sequence of "NC_tetst.1" is stored in ``test_case1.fa``. 
Moreover, the third nucleotide is replaced from G to c. Moreover, The sixth nucleotide T is deleted.
If we check ``test_case2.fa``, the modification is also according to the mutation table and our setting.

If we have no mutation table, we can also use subcommand ``snp`` to detect mutations and generate 
fasta files automatically. For ``snp``, we will go through it later.

Generating annotation files
---------------------------

We have fasta files of our new dummy query genome now. We can use them to generate annotation files. If the annotation files 
retrieved from NCBI is exactly what you want, you can skip this step. 

Before we running this subcommand, we have to modify environment paths of `RATT <http://ratt.sourceforge.net/>`_. 
If you execute ANNOgesic by using docker container, the path is already set. 
If you setup ANNOgesic by yourself, please check 
`RATT <http://ratt.sourceforge.net/>`_ to set your environment paths properly.

After setting the environment, we can try it.

::

    anngesic annotation_transfer \
        -ce ANNOgesic/input/references/annotations/NC_009839.1.embl \
        -cf ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        -uf ANNOgesic/output/updated_references/fasta_files/test_case1.fa \
            ANNOgesic/output/updated_references/fasta_files/test_case2.fa \
        -e chromosome \
        -t Strain \
        -p NC_009839.1:NC_test.1 NC_009839.1:test_case2 \
        -g \
        -pj ANNOgesic


``-e`` is prefix name of the output embl files. ``-t`` is a program of `RATT <http://ratt.sourceforge.net/>`_.
We use ``Strain`` because the similarity is higher than 90%. For other situations, please check 
`RATT <http://ratt.sourceforge.net/>`_. In ``-p``, we assign pairs of the target genomes (NC_test.1 and test_case2) 
and their close genomes (NC_000915.1). Please be careful, the information that we assign to ``-p`` 
is genome names in gff files not fasta filenames. ``-g`` means that we want to transfer the 
output embl files to GFF3 files and store in ``ANNOgesic/output/updated_references/annotations``.

Once the transfer is done, we can see

::

    $ ls ANNOgesic/output/updated_references/annotations/
    test_case1.gff  test_case1.ptt  test_case1.rnt  test_case2.gff  test_case2.ptt  test_case2.rnt
    $ ls ANNOgesic/output/annotation_transfer/
    chromosome.NC_test.1.final.embl  chromosome.test_case2.final.embl  NC_test.1.gff  ratt_log.txt  test_case2.gff

In ``ANNOgesic/output/updated_references/annotations``, we can find ptt, rnt and gff files. In ``ANNOgesic/output/annotation_transfer``,
we can find the output of `RATT <http://ratt.sourceforge.net/>`_.

We already saw how to update genome fasta and annotation files. 
We will use ``ANNOgesic/input/references/annotations/NC_009839.1.gff`` and ``ANNOgesic/input/references/fasta_files/NC_009839.1`` 
for running the following subcommands. If the fasta files and annotation files of your genome need to be updated, 
please replace the files with the fasta and annotation files in ``ANNOgesic/output/updated_references``.

TSS and processing site prediction and optimization
---------------------------------------------------

Before running following subcommands, we need to setup our libraries as a correct format.
First, we set the path of wig file folder.

::

    WIG_FOLDER="ANNOgesic/input/wigs/tex_notex"

Then, we can setup our libraries.

::

    TEX_LIBS="$WIG_FOLDER/GSM951380_Log_81116_R1_minus_TEX_in_NC_009839_minus.wig:notex:1:a:- \
              $WIG_FOLDER/GSM951381_Log_81116_R1_plus_TEX_in_NC_009839_minus.wig:tex:1:a:- \
              $WIG_FOLDER/GSM951380_Log_81116_R1_minus_TEX_in_NC_009839_plus.wig:notex:1:a:+ \
              $WIG_FOLDER/GSM951381_Log_81116_R1_plus_TEX_in_NC_009839_plus.wig:tex:1:a:+"

Now, we can start to test other subcommands. 
Before running ``tss_ps``, if we want to use the optimized parameters, 
we need to run ``optimize_tss_ps`` first. The optimization requires a gff file of the manual-detected TSSs. 
In our experience, we recommend you to detect at least 50 TSSs and check more than 200kb of genome. 

For the test case, you can download the `manual TSS file <https://github.com/Sung-Huan/ANNOgesic/tree/master/tutorial_data>`_ 
from our git repository. 

::

    $ wget -cP ANNOgesic/input/manual_TSSs/ https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/tutorial_data/NC_009839_manual_TSS.gff

Now, we have a manual TSS gff file which is stored in ``ANNOgesic/input/manual_TSSs``. 
we can try ``optimize_tss_ps`` right now (since we only check first 200kb, we set ``-le`` as "NC_009839.1:200000" which 
means only first 200kb of NC_009839.1 is valid.).

::

    $ annogesic optimize_tss_ps \
        -f ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        -g ANNOgesic/input/references/annotations/NC_009839.1.gff \
        -tl $TEX_LIBS \
        -p TSS -s 25 \
        -m ANNOgesic/input/manual_TSSs/NC_009839_manual_TSS.gff \
        -le NC_009839.1:200000 \
        -rt all_1 \
        -pj ANNOgesic

``optimize_tss_ps`` will compare manual checked TSSs with predicted TSSs to search the best parameters. 
Results of the different parameters will be printed in the screen. We only set 25 runs for testing. 
Once the optimization is done, you can find several files.

::

    $ ls ANNOgesic/output/TSSs/optimized_TSSpredator/
    best_NC_009839.1.csv  log.txt  stat_NC_009839.1.csv

``best_NC_009839.1.csv`` is for the best parameters; ``stat_NC_009839.1.csv`` is for parameters of each step.

Now, we assume the best parameters are following: height is 0.4, height_reduction is 0.1, factor is 1.7, factor_reduction is 0.2, 
base_height is 0.039, enrichment_factor is 1.1, processing_factor is 4.5. We can set these parameters for running  
``tss``.

::

    $ annogesic tss_ps \
        -f ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        -g ANNOgesic/input/references/annotations/NC_009839.1.gff \
        -tl $TEX_LIBS \
        -p test \
        -he 0.4 \
        -rh 0.1 \
        -fa 1.7 \
        -rf 0.2 \
        -bh 0.039 \
        -ef 1.1 \
        -pf 4.5 \
        -v \
        -rt all_1 \
        -le NC_009839.1:200000 \
        -m ANNOgesic/input/manual_TSSs/NC_009839_manual_TSS.gff \
        -pj ANNOgesic

We assign the manual-checked TSS gff file to ``-m``. Therefore, the output gff file contains the manual-detected TSSs and predicted TSSs. 
If we didn't assign it, Only the predicted TSSs will be included in output gff file. 
The output files are gff file, MasterTable and statistic files.

::

    $ ls ANNOgesic/output/TSSs/
    configs  gffs  MasterTables  optimized_TSSpredator  screenshots  statistics
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
processing_site (``-t``) and assign the proper parameter sets. We assume the best parameter sets are following: 
height is 0.2, height_reduction is 0.1, factor is 2.0, factor_reduction is 0.5,
base_height is 0.009, enrichment_factor is 1.2, processing_factor is 1.5.

::

    $ annogesic tss_ps \
        -f ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        -g ANNOgesic/input/references/annotations/NC_009839.1.gff \
        -tl $TEX_LIBS \
        -p test \
        -he 0.2 \
        -rh 0.1 \
        -fa 2.0 \
        -rf 0.5 \ 
        -bh 0.009 \
        -ef 1.2 \
        -pf 1.5 \ 
        -rt all_1 \
        -t processing_site \
        -pj ANNOgesic

The output files are following:

::

    $ ls ANNOgesic/output/processing_sites/
    configs  gffs  MasterTables  statistics
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
we can use subcommand ``transcript`` to do it. Normally, we strongly 
recommend that the user should provide fragmented libraries. Because dRNA-Seq usually loses some information 
of 3'end. However, we only use TEX +/- for testing.

There are several options for modifying transcripts by comparing transcripts and genome annotations. 
If you want to know the details, please check :ref:`transcript`. Now, we use default setting run this module: 

::

    $ annogesic transcript \
        -g ANNOgesic/input/references/annotations/NC_009839.1.gff \
        -tl $TEX_LIBS \
        -rt all_1 \
        -cf gene CDS \
        -ct ANNOgesic/output/TSSs/gffs/NC_009839.1_TSS.gff \
        -pj ANNOgesic

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

We can use subcommand ``terminator`` to detect terminators. ``terminator`` integrates `RNAfold <http://www.tbi.univie.ac.at/RNA/RNAfold.1.html>`_ 
for computing secondary structure of potential terminators. Therefore, this process will take a while. The command is like following: 

::

    $ annogesic terminator \
        -f ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        -g ANNOgesic/input/references/annotations/NC_009839.1.gff \
        -a ANNOgesic/output/transcripts/gffs/NC_009839.1_transcript.gff \
        -tl $TEX_LIBS \
        -rt all_1 -tb \
        -pj ANNOgesic

Four different kinds of gff files and tables will be generated.

::

    $ ls ANNOgesic/output/terminators/gffs/
    all_candidates  best_candidates  expressed_candidates  non_expressed_candidates
    $ ls ANNOgesic/output/terminators/tables
    all_candidates  best_candidates  expressed_candidates  non_expressed_candidates

``all_candidates`` is for all candidates; ``expressed_candidates`` is for the candidates which reveal gene expression; 
``best_candidates`` is for the candidates which reveal gene expression and their coverage shows significant decreasing; 
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

Now, we have the information of TSSs, transcripts and terminators. We can detect the 5'UTRs and 3'UTRs easily by using 
subcommand ``utr``.

::

    $ annogesic utr \
        -g ANNOgesic/input/references/annotations/NC_009839.1.gff \
        -t ANNOgesic/output/TSSs/gffs/NC_009839.1_TSS.gff \
        -a ANNOgesic/output/transcripts/gffs/NC_009839.1_transcript.gff \
        -e ANNOgesic/output/terminators/gffs/best_candidates/NC_009839.1_term.gff \
        -pj ANNOgesic

If the TSS gff file is not generated by ANNOgesic, please assign ``-s``,  the TSSs can be classified for generating UTRs.
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

We already had TSSs, transcripts, terminators, CDSs, UTRs. We can integrate all these feature to 
detect operons and suboperons by executing subcommand ``operon``.

::

    $ annogesic operon \
        -g ANNOgesic/input/references/annotations/NC_009839.1.gff \
        -t ANNOgesic/output/TSSs/gffs/NC_009839.1_TSS.gff \
        -a ANNOgesic/output/transcripts/gffs/NC_009839.1_transcript.gff \
        -u5 ANNOgesic/output/UTRs/5UTRs/gffs/NC_009839.1_5UTR.gff \
        -u3 ANNOgesic/output/UTRs/3UTRs/gffs/NC_009839.1_3UTR.gff \
        -e ANNOgesic/output/terminators/gffs/best_candidates/NC_009839.1_term.gff \
        -pj ANNOgesic

Three folders will be generated to store gff files, tables and statistics files.

::

    $ ls ANNOgesic/output/operons/
    gffs  statistics  tables
    $ ls ANNOgesic/output/operons/gffs/
    NC_009839.1_operon.gff
    $ ls ANNOgesic/output/operons/tables/
    NC_009839.1_operon.csv
    $ ls ANNOgesic/output/operons/statistics/
    stat_NC_009839.1_operon.csv

Promoter motif detection
------------------------

As long as we have TSSs, we can use subcommand ``promoter`` to get promoters. The promoters can be detected 
by different types of the TSSs. Therefore, if the TSSs gff files are not generated by ``ANNOgesic``,
you need to add ``-s`` and assign corresponding genome annotation file to ``-g``.
Now, let try ``promoter`` by running MEME and GLAM2 (``-p`` is assigned by "both" in default. If you want to only run 
MEME or GLAM2, please assign "meme" or "glam2" to ``-p``), the process may take a while.

::

    $ annogesic promoter \
        -t ANNOgesic/output/TSSs/gffs/NC_009839.1_TSS.gff \
        -f ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        -w 45 2-10 \
        -pj ANNOgesic

We define the length of the motifs as ``50`` and ``2-10``. ``2-10`` means the width can be from 2 to 10.

Based on different types of the TSSs and the length of the motif, numerous output files will be generated.

::

    $ ls ANNOgesic/output/promoters/
    fasta_classes  NC_009839.1
    $ ls ANNOgesic/output/promoters/fasta_classes/NC_009839.1
    NC_009839.1_allgenome_all_types.fa  NC_009839.1_allgenome_internal.fa  NC_009839.1_allgenome_primary.fa    NC_009839.1_allgenome_without_orphan.fa
    NC_009839.1_allgenome_antisense.fa  NC_009839.1_allgenome_orphan.fa    NC_009839.1_allgenome_secondary.fa
    $ ls ANNOgesic/output/promoters/NC_009839.1
    MEME GLAM2
    $ ls ANNOgesic/output/promoters/NC_009839.1/MEME
    promoter_motifs_NC_009839.1_allgenome_all_types_2-10_nt  promoter_motifs_NC_009839.1_allgenome_internal_45_nt   promoter_motifs_NC_009839.1_allgenome_secondary_2-10_nt
    promoter_motifs_NC_009839.1_allgenome_all_types_45_nt    promoter_motifs_NC_009839.1_allgenome_orphan_2-10_nt   promoter_motifs_NC_009839.1_allgenome_secondary_45_nt
    promoter_motifs_NC_009839.1_allgenome_antisense_2-10_nt  promoter_motifs_NC_009839.1_allgenome_orphan_45_nt     promoter_motifs_NC_009839.1_allgenome_without_orphan_2-10_nt
    promoter_motifs_NC_009839.1_allgenome_antisense_45_nt    promoter_motifs_NC_009839.1_allgenome_primary_2-10_nt  promoter_motifs_NC_009839.1_allgenome_without_orphan_45_nt
    promoter_motifs_NC_009839.1_allgenome_internal_2-10_nt   promoter_motifs_NC_009839.1_allgenome_primary_45_nt
    $ ls ANNOgesic/output/promoters/NC_009839.1/GLAM2
    promoter_motifs_NC_009839.1_allgenome_all_types_2-10_nt  promoter_motifs_NC_009839.1_allgenome_internal_45_nt   promoter_motifs_NC_009839.1_allgenome_secondary_2-10_nt
    promoter_motifs_NC_009839.1_allgenome_all_types_45_nt    promoter_motifs_NC_009839.1_allgenome_orphan_2-10_nt   promoter_motifs_NC_009839.1_allgenome_secondary_45_nt
    promoter_motifs_NC_009839.1_allgenome_antisense_2-10_nt  promoter_motifs_NC_009839.1_allgenome_orphan_45_nt     promoter_motifs_NC_009839.1_allgenome_without_orphan_2-10_nt
    promoter_motifs_NC_009839.1_allgenome_antisense_45_nt    promoter_motifs_NC_009839.1_allgenome_primary_2-10_nt  promoter_motifs_NC_009839.1_allgenome_without_orphan_45_nt
    promoter_motifs_NC_009839.1_allgenome_internal_2-10_nt   promoter_motifs_NC_009839.1_allgenome_primary_45_nt
    $ ls ANNOgesic/output/promoters/NC_009839.1/MEME/promoter_motifs_NC_009839.1_allgenome_all_types_45_nt/
    logo1.eps  logo1.png  logo2.eps  logo2.png  logo3.eps  logo3.png  logo_rc1.eps  logo_rc1.png  logo_rc2.eps  logo_rc2.png  logo_rc3.eps  logo_rc3.png  meme.csv  meme.html  meme.txt  meme.xml
    $ ls ANNOgesic/output/promoters/NC_009839.1/GLAM2/promoter_motifs_NC_009839.1_allgenome_all_types_45_nt/
    glam2.csv   glam2.txt   logo1.eps  logo2.png  logo4.eps  logo5.png  logo7.eps  logo8.png  logo_ssc10.eps  logo_ssc1.png  logo_ssc3.eps  logo_ssc4.png  logo_ssc6.eps  logo_ssc7.png  logo_ssc9.eps
    glam2.html  logo10.eps  logo1.png  logo3.eps  logo4.png  logo6.eps  logo7.png  logo9.eps  logo_ssc10.png  logo_ssc2.eps  logo_ssc3.png  logo_ssc5.eps  logo_ssc6.png  logo_ssc8.eps  logo_ssc9.png
    glam2.meme  logo10.png  logo2.eps  logo3.png  logo5.eps  logo6.png  logo8.eps  logo9.png  logo_ssc1.eps   logo_ssc2.png  logo_ssc4.eps  logo_ssc5.png  logo_ssc7.eps  logo_ssc8.png

Prediction of sRNA and sORF
---------------------------

Based on transcripts, genome annotation and coverage information, sRNAs can be detected. Moreover, we 
have TSSs and processing sites which can be used for detecting UTR-derived sRNAs as well. Now, we can 
get sRNAs by running subcommand ``srna``. Normally, we recommend that the user inputs fragmented libraries as well.
Here, we only use TEX +/- for testing.

For running ``srna``, we can apply several filters to improve the detection. These filters are ``tss``, ``sec_str``,
``blast_nr``, ``blast_srna``, ``promoter``, ``term``, ``sorf``. Normally, ``tss``, ``sec_str``,
``blast_nr``, ``blast_srna`` are recommended to used.

Please be aware, filters are strict. For example, if your filters are included ``term``, only the sRNAs which are 
associated with terminators will be included in best list. If you want to include terminator information 
but not use terminator as a filter, you can remove ``term`` in filters and still assign the path of terminator gff file. 
The results will include the sRNAs which are not associated with terminators and also store terminator information.

Before running ``srna``, we have to get sRNA database (we can use `BSRD <http://www.bac-srna.org/BSRD/index.jsp>`_) and 
`nr database <ftp://ftp.ncbi.nih.gov/blast/db/FASTA/>`_ (if you have not downloaded before). 
We can download fasta file of `BSRD <http://www.bac-srna.org/BSRD/index.jsp>`_ from our 
`Git repository <https://github.com/Sung-Huan/ANNOgesic/tree/master/database>`_.

::

    $ wget -cP ANNOgesic/input/databases/ https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/database/sRNA_database_BSRD.fa



If you already had sRNA database in other folders, please assign your path of databases to ``-sd``.
If your database is formatted before, you can remove ``-sf``.
In order to use the recommended filters, we have to download 
`nr database <ftp://ftp.ncbi.nih.gov/blast/db/FASTA/>`_ (takes a while). If you already had it, 
you can skip this step.

::

    $ wget -cP ANNOgesic/input/databases/ ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
    $ gunzip ANNOgesic/input/databases/nr.gz
    $ mv ANNOgesic/input/databases/nr ANNOgesic/input/databases/nr.fa

If your nr database is in other folders, please assign your path ``-nd``.
You can also remove ``-nf`` if your database is already formatted.
Now, we can use the recommended filters to run ``srna``, but it may takes several hours.

::

    $ annogesic srna \
        -d tss blast_srna sec_str blast_nr \
        -g ANNOgesic/input/references/annotations/NC_009839.1.gff \
        -t ANNOgesic/output/TSSs/gffs/NC_009839.1_TSS.gff \
        -p ANNOgesic/output/processing_sites/gffs/NC_009839.1_processing.gff \
        -a ANNOgesic/output/transcripts/gffs/NC_009839.1_transcript.gff \
        -f ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        -tf ANNOgesic/output/terminators/gffs/best_candidates/NC_009839.1_term.gff \
        -pt ANNOgesic/output/promoters/NC_009839.1/MEME/promoter_motifs_NC_009839.1_allgenome_all_types_45_nt/meme.csv \
        -pn MOTIF_1 \
        -m \
        -u \
        -cs \
        -sf \
        -nf \
        -nd ANNOgesic/input/databases/nr \
        -sd ANNOgesic/input/databases/sRNA_database_BSRD \
        -tl $TEX_LIBS \
        -rt all_1 \
        -pj ANNOgesic

If you have sORF information, you can also assign path of the sORF gff folder to ``-O``. 
Then, comparison of sRNAs and sORFs can be done.

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
    srna0_NC_009839.1_36954_37044_-_dp.pdf     srna25_NC_009839.1_854600_854673_-_dp.pdf    srna40_NC_009839.1_1091155_1091251_-_dp.pdf  srna56_NC_009839.1_1440826_1441414_+_dp.pdf
    srna10_NC_009839.1_248098_248257_-_dp.pdf  srna26_NC_009839.1_879881_880088_-_dp.pdf    srna41_NC_009839.1_1097654_1097750_-_dp.pdf  srna57_NC_009839.1_1448211_1448306_+_dp.pdf
    ...

    $ ls ANNOgesic/output/sRNAs/figs/sec_plots/NC_009839.1/
    rna0_NC_009839.1_36954_37044_-_rss.pdf     srna25_NC_009839.1_854600_854673_-_rss.pdf    srna40_NC_009839.1_1091155_1091251_-_rss.pdf  srna56_NC_009839.1_1440826_1441414_+_rss.pdf
    srna10_NC_009839.1_248098_248257_-_rss.pdf  srna26_NC_009839.1_879881_880088_-_rss.pdf    srna41_NC_009839.1_1097654_1097750_-_rss.pdf  srna57_NC_009839.1_1448211_1448306_+_rss.pdf
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

As we know, expressed region without annotation may be sORF as well. 
In order to get information of sORFs, we can use subcommand ``sorf``.

::

    $ annogesic sorf \
        -g ANNOgesic/input/references/annotations/NC_009839.1.gff \
        -t ANNOgesic/output/TSSs/gffs/NC_009839.1_TSS.gff \
        -a ANNOgesic/output/transcripts/gffs/NC_009839.1_transcript.gff \
        -f ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        -s ANNOgesic/output/sRNAs/gffs/best_candidates/NC_009839.1_sRNA.gff \
        -tl $TEX_LIBS \
        -rt all_1 -u \
        -pj ANNOgesic

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

    $ annogesic srna_target \
        -g ANNOgesic/input/references/annotations/NC_009839.1.gff \
        -f ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        -r ANNOgesic/output/sRNAs/gffs/best_candidates/NC_009839.1_sRNA.gff \
        -q NC_009839.1:36954:37044:- \
        -p both \
        -pj ANNOgesic

For testing, we only assign one sRNA to do the prediction. You can also assign several of sRNAs like 
``NC_009839.1:36954:37044:- NC_009839.1:75845:75990:+``. If you want to compute all sRNAs, you 
can assign ``all`` to ``-q`` (may take several days).

Several output folders will be generated. 

::

    $ ls ANNOgesic/output/sRNA_targets/
    merged_results  RNAplex_results  RNAup_results  sRNA_seqs  target_seqs

``sRNA_seqs`` and ``target_seqs`` are for sequences of the sRNAs and the potential targets.

::

    $ ls ANNOgesic/output/sRNA_targets/sRNA_seqs
    NC_009839.1_sRNA.fa
    $ ls ANNOgesic/output/sRNA_targets/target_seqs
    NC_009839.1_target.fa

``RNAplex_results`` and ``RNAup_results`` are for output of `RNAplex and RNAup <http://www.tbi.univie.ac.at/RNA/>`_.

::

    $ ls ANNOgesic/output/sRNA_targets/RNAplex_results/NC_009839.1/
    NC_009839.1_RNAplex_rank.csv  NC_009839.1_RNAplex.txt
    $ ls ANNOgesic/output/sRNA_targets/RNAup_results/NC_009839.1/
    NC_009839.1_RNAup.log  NC_009839.1_RNAup_rank.csv  NC_009839.1_RNAup.txt

``merged_results`` is for the merged results of `RNAplex <http://www.tbi.univie.ac.at/RNA/RNAplex.1.html>`_ and 
`RNAup <http://www.tbi.univie.ac.at/RNA/RNAup.1.html>`_. ``NC_009839.1_merge.csv``  contains all results of the 
both methods. ``NC_009839.1_overlap.csv`` only stores candidates which are top 20 (default) in the both methods.

::

    $ ls ANNOgesic/output/sRNA_targets/merged_results/NC_009839.1/
    NC_009839.1_merge.csv  NC_009839.1_overlap.csv

Mapping and detecting of circular RNA
-------------------------------------

You may also be interested in circular RNAs. The subcommand ``circrna`` can help us to get circular RNAs by  
using `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_. Since 
we didn't map reads of the test case before, we can also do mapping by running ``circrna``. If you already mapped 
the reads by `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_ with ``-S``, then you can 
remove ``-a`` and add path of the bam files to ``-nb`` or ``-fb``. However, 
if you mapped the reads by other tools or you mapped the reads by 
`Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_ without ``-S``, Unfortunately, 
you have to re-map the reads again. You can assign parallel (``-p``) for mapping.

In normal situation, the reads should be directly given to ``circrna``. However, we just want to test the 
subcommand. Thus, we can reduce the running time by selecting the subset of reads (first 50000) for only testing.

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


After that, we assign ``all_samples:$READ_FILE`` to ``-rp``. ``all_sample`` is the name of the set of read files. 
The all four read files will be compute together. Now, we can try ``circrna``

::

     $ annogesic circrna \
         -f ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
         -p 10 \
         -g ANNOgesic/input/references/annotations/NC_009839.1.gff \
         -rp all_samples:$READ_FILES \
         -pj ANNOgesic

If you can't find testrealign.x, please refer to :ref:`Required tools or databases`.
Several output folders will be generated.

::

    $ ls ANNOgesic/output/circRNAs/
    circRNA_tables  gffs  segemehl_alignment_files  segemehl_splice_results  statistics

``segemehl_alignment_files`` and ``segemehl_splice_results`` are for results of 
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
the circular RNAs after filtering. In our case, there are some circular RNAs can be detected, but no circular RNAs 
can exist after filtering.

SNP calling
--------------

If we want to know SNPs or mutations of our RNA-seq data, we can use ``snp`` to achieve this purpose.
``snp`` is compose of two parts. One part is for obtaining the differences between our query genome 
and the close genome of our query one. If we have no fasta file of our query genome, 
this part will be very useful. We just need to map reads of our query genome on the fasta file of the closed genome. Then 
using ``snp`` can automatically detect differences between the closed genome and our query genome. 
Furthermore, potential fasta files of our query genome can be generated automatically as well. 
The other part is for detecting SNPs or mutations of our query genome if the fasta file of our query genome can be provided.
In this part, you can know real mutations of our query genonme.

Before running the subcommand, bam files are required. Since we already generated them through 
running ``circrna``, we can just put them to corresponding folder. Please remember that the mapping function of 
``circrna`` is very basic.

Now, we can try to detect mutations. Since we already got the Bam files of NC_009839.1 (our query genome) via ``circrna``, 
we can set ``-t`` as ``query_genome``. The procedures of comparing closed genome and query genome are similar, 
you just need to put Bam files, and fasta files to corresponding folders and set ``-t`` as ``closed_genome``.

First, we copy the bam files to ``BAMs_map_query_genomes``.

::

    $ cp ANNOgesic/output/circRNAs/segemehl_alignment_files/NC_009839.1/SRR51525* ANNOgesic/input/BAMs/BAMs_map_query_genomes/tex_notex

Now, we can set our bam files

::
    $ BAM_FILES=ANNOgesic/input/BAMs/BAMs_map_query_genomes/tex_notex/SRR515254_50000_NC_009839.1.bam,\
      ANNOgesic/input/BAMs/BAMs_map_query_genomes/tex_notex/SRR515255_50000_NC_009839.1.bam,\
      ANNOgesic/input/BAMs/BAMs_map_query_genomes/tex_notex/SRR515256_50000_NC_009839.1.bam,\
      ANNOgesic/input/BAMs/BAMs_map_query_genomes/tex_notex/SRR515257_50000_NC_009839.1.bam

Then we can run the subcommand with three programs -- ``extend_BAQ``, ``with_BAQ`` and ``without_BAQ``. 
``all_sample:2:$BAM_FILES`` for ``-b`` means the name of the set of bam files is all_sample, there are two 
samples in this set, and all four bam files need to be compute together.

::

    $ annogesic snp \
        -t query_genome \
        -p with_BAQ without_BAQ extend_BAQ \
        -b all_samples:2:$BAM_FILES \
        -f ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        -pj ANNOgesic

Two output folders will be generated, ``compare_closed_and_updated_references`` is for results of the comparison between closed genome 
and query genome, ``mutations_of_query_genomes`` is for results of detecting mutations of the query genome.

::

    $ ls ANNOgesic/output/SNP_calling/                                                                                                      
    compare_closed_and_updated_references  mutations_of_query_genomes

Since we run ``query_genome``,  the output folders are produced under ``mutations_of_query_genomes``.

::

    $ ls ANNOgesic/output/SNP_calling/mutations_of_query_genomes/
    seqs  SNP_raw_outputs  SNP_tables  statistics

The output folders are compose of three parts - ``extend_BAQ``, ``with_BAQ`` and ``without_BAQ``.

::

    $ ls ANNOgesic/output/SNP_calling/mutations_of_query_genomes/seqs/
    extend_BAQ/  with_BAQ/    without_BAQ/

In ``seqs``, the potential sequences can be found.

::

    $ ls ANNOgesic/output/SNP_calling/mutations_of_query_genomes/seqs/with_BAQ/NC_009839.1/
    NC_009839.1_all_samples_NC_009839.1_1_1.fa

``SNP_raw_outputs`` stores output of `Samtools and Bcftools <https://github.com/samtools>`_. 
``SNP_tables`` stores results after filtering and the indices of potential sequence 
(potential sequences are stored in ``seqs``).
``statistics`` stores the statistic files.

::

    $ ls ANNOgesic/output/SNP_calling/mutations_of_query_genomes/SNP_raw_outputs/NC_009839.1/
    NC_009839.1_extend_BAQ_all_samples.vcf  NC_009839.1_with_BAQ_all_samples.vcf  NC_009839.1_without_BAQ_all_samples.vcf
    $ ls ANNOgesic/output/SNP_calling/mutations_of_query_genomes/SNP_tables/NC_009839.1/
    NC_009839.1_extend_BAQ_all_samples_best.vcf           NC_009839.1_with_BAQ_all_samples_best.vcf           NC_009839.1_without_BAQ_all_samples_best.vcf
    NC_009839.1_extend_BAQ_all_samples_seq_reference.csv  NC_009839.1_with_BAQ_all_samples_seq_reference.csv  NC_009839.1_without_BAQ_all_samples_seq_reference.csv
    $ ls ANNOgesic/output/SNP_calling/mutations_of_query_genomes/statistics/
    figs                                                  stat_NC_009839.1_with_BAQ_all_samples_SNP_best.csv     stat_NC_009839.1_without_BAQ_all_samples_SNP_raw.csv
    stat_NC_009839.1_extend_BAQ_all_samples_SNP_best.csv  stat_NC_009839.1_with_BAQ_all_samples_SNP_raw.csv
    stat_NC_009839.1_extend_BAQ_all_samples_SNP_raw.csv   stat_NC_009839.1_without_BAQ_all_samples_SNP_best.csv
    $ ls ANNOgesic/output/SNP_calling/mutations_of_query_genomes/statistics/figs
    NC_009839.1_extend_BAQ_all_samples_NC_009839.1_SNP_QUAL_best.png  NC_009839.1_with_BAQ_all_samples_NC_009839.1_SNP_QUAL_best.png  NC_009839.1_without_BAQ_all_samples_NC_009839.1_SNP_QUAL_best.png
    NC_009839.1_extend_BAQ_all_samples_NC_009839.1_SNP_QUAL_raw.png   NC_009839.1_with_BAQ_all_samples_NC_009839.1_SNP_QUAL_raw.png   NC_009839.1_without_BAQ_all_samples_NC_009839.1_SNP_QUAL_raw.png

Mapping Gene ontology
---------------------

Gene ontology is useful for understanding function of gene products. 
Implementing ``go_term`` can map our annotations to gene ontology. Before running ``go_term``, we 
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
generate results which only included the expressed CDS.

Let's try it.

::

    $ annogesic go_term \
        -g ANNOgesic/input/references/annotations/NC_009839.1.gff \
        -a ANNOgesic/output/transcripts/gffs/NC_009839.1_transcript.gff \
        -go ANNOgesic/input/databases/go.obo \
        -gs ANNOgesic/input/databases/goslim_generic.obo \
        -u ANNOgesic/input/databases/idmapping_selected.tab \
        -pj ANNOgesic

Output of ``go_term`` are stored in ``GO_term_results``. The statistic files and 
figures are stored in ``statistics``.

::

    $ ls ANNOgesic/output/GO_terms/
    all_CDSs  expressed_CDSs
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

Subcellular localization is also a useful information for analysis of protein function. For 
detecting subcellular localization, we can use the subcommand 
``localization``. Like ``go_term``, we can also import 
information of the transcript to generate results which only included the expressed CDS.

::

    $ annogesic localization \
        -g ANNOgesic/input/references/annotations/NC_009839.1.gff \
        -f ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        -a ANNOgesic/output/transcripts/gffs/NC_009839.1_transcript.gff \
        -m -b negative \
        -pj ANNOgesic

Two output folders will be generated. ``psortb_results`` stores output 
of `Psortb <http://www.psort.org/psortb/>`_. ``statistics`` stores 
statistic files and figures.

::

    $ ls ANNOgesic/output/subcellular_localization/
    all_CDSs  expressed_CDSs
    $ ls ANNOgesic/output/subcellular_localization/all_CDSs/
    psortb_results  statistics
    $ ls ANNOgesic/output/subcellular_localization/all_CDSs/psortb_results/NC_009839.1/
    NC_009839.1_raw.txt  NC_009839.1_table.csv
    $ ls ANNOgesic/output/subcellular_localization/all_CDSs/statistics/NC_009839.1/
    NC_009839.1_NC_009839.1_sublocal.png  stat_NC_009839.1_sublocal.csv

Generating protein-protein interaction network
----------------------------------------------

``ppi_network`` can detect protein-protein interaction from `STRING <http://string-db.org/>`_ 
(database of protein-protein interaction) and searching the literatures by implementing 
`PIE <http://www.ncbi.nlm.nih.gov/CBBresearch/Wilbur/IRET/PIE/>`_ 
(text-mining for protein-protein interaction). Therefore, ``ppi_network`` can generate protein-protein 
interaction networks with supported literatures.

Before running the subcommand, you need to download 
`species.v{$VERSIO}.txt from STRING <http://string-db.org/cgi/download.pl>`_

::

    $ wget -cP ANNOgesic/input/databases http://string-db.org/newstring_download/species.v10.txt

Now, we can try the subcommand.

::

    $ annogesic ppi_network \
        -s NC_009839.1.gff:NC_009839.1:'Campylobacter jejuni 81176':'Campylobacter jejuni' \
        -g ANNOgesic/input/references/annotations/NC_009839.1.gff \
        -d ANNOgesic/input/databases/species.v10.txt \
        -q NC_009839.1:70579:71463:+ NC_009839.1:102567:103973:+ \
        -n \
        -pj ANNOgesic

We only detected for two proteins. If you want to detect for all proteins in ptt files, 
you can easily assign ``all`` in ``-q``.

Three output folders will be generated.

::

    $ ls ANNOgesic/output/PPI_networks/
    all_results/  best_results/ figures/

``all_results`` is for all interactions without filtering. ``best_results`` is for the interactions with 
the high `PIE <http://www.ncbi.nlm.nih.gov/CBBresearch/Wilbur/IRET/PIE/>`_ score. ``figures`` is for 
figures of the protein-protein interaction networks. There are two subfolders - ``with_strain`` and ``without_strain`` in 
``figures``. These two folders store all information of the interactions and literature scores. 
``with_strain`` is for results with assigning specific strain name for searching literature. 
``without_strain`` is for results without giving specific strain name for searching literature.

::

    $ ls ANNOgesic/output/PPI_networks/all_results/PPI_NC_009839.1/
    NC_009839.1_without_strain.csv  NC_009839.1_with_strain.csv  without_strain  with_strain
    $ ls ANNOgesic/output/PPI_networks/best_results/PPI_NC_009839.1/
    NC_009839.1_without_strain.csv  NC_009839.1_with_strain.csv  without_strain  with_strain
    $ ls ANNOgesic/output/PPI_networks/figures/PPI_NC_009839.1/
    without_strain  with_strain
    $ ls ANNOgesic/output/PPI_networks/all_results/PPI_NC_009839.1/with_strain/NC_009839.1/
    flgB_flgD.csv    flgE_flgD.csv  flgF_fliG.csv  flgG_fliG.csv  fliG_fliF.csv
    flgE-1_flgD.csv  flgF_flgC.csv  flgG_flgC.csv  flgI_flgH.csv  pyrB_ansA.csv
    $ ls ANNOgesic/output/PPI_networks/all_results/PPI_NC_009839.1/without_strain/NC_009839.1/
    flgB_flgD.csv    flgE_flgD.csv  flgF_fliG.csv  flgG_fliG.csv  fliG_fliF.csv
    flgE-1_flgD.csv  flgF_flgC.csv  flgG_flgC.csv  flgI_flgH.csv  pyrB_ansA.csv
    $ ls ANNOgesic/output/PPI_networks/best_results/PPI_NC_009839.1/without_strain/NC_009839.1/
    flgB_flgD.csv    flgE_flgD.csv  flgG_flgC.csv  fliG_fliF.csv
    flgE-1_flgD.csv  flgF_flgC.csv  flgI_flgH.csv
    $ ls ANNOgesic/output/PPI_networks/best_results/PPI_NC_009839.1/with_strain/NC_009839.1/
    fliG_fliF.csv
    $ ls ANNOgesic/output/PPI_networks/figures/PPI_NC_009839.1/with_strain/NC_009839.1/
    C8J_RS00250_flgD.png
    $ ls ANNOgesic/output/PPI_networks/figures/PPI_NC_009839.1/without_strain/NC_009839.1/
    C8J_RS00250_flgD.png

Generating riboswitch and RNA thermometer
-----------------------------------------

If we want to detect riboswitches and RNA thermometer, we can use subcommand ``riboswitch_thermometer``.
Before running it, we need to get information of the known riboswitches and RNA thermometer in Rfam. 
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
        -g ANNOgesic/input/references/annotations/NC_009839.1.gff \
        -f ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        -ri ANNOgesic/input/riboswitch_ID_file/Rfam_riboswitch_ID.csv \
        -ti ANNOgesic/input/RNA_thermometer_ID_file/Rfam_RNA_thermometer_ID.csv \
        -R ANNOgesic/input/databases/CMs/Rfam.cm \
        -a ANNOgesic/output/transcripts/gffs/NC_009839.1_transcript.gff \
        -t ANNOgesic/output/TSSs/gffs/NC_009839.1_TSS.gff \
        -pj ANNOgesic

Output files are following, ``gffs`` stores gff files of the riboswitchs / RNA_thermometer; 
``tables`` stores tables of the riboswitchs / RNA_thermometer; 
``scan_Rfam_results`` stores output files of scanning Rfam; ``statistics`` is for statistic files.

::

     $ ls ANNOgesic/output/riboswitches/
     gffs  scan_Rfam_results  statistics  tables
     $ ls ANNOgesic/output/riboswitches/gffs/
     NC_009839.1_riboswitch.gff
     $ ls ANNOgesic/output/riboswitches/scan_Rfam_results/NC_009839.1/
     NC_009839.1_riboswitch_prescan.txt  NC_009839.1_riboswitch_scan.txt
     $ ls ANNOgesic/output/riboswitches/tables/
     NC_009839.1_riboswitch.csv
     $ ls ANNOgesic/output/riboswitches/statistics/
     stat_NC_009839.1_riboswitch.txt
     $ ls ANNOgesic/output/RNA_thermometers/
     gffs  scan_Rfam_results  statistics  tables
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
CRISPR is an unique features for research of immunology. ``crispr`` is a useful subcommand for CRISPR detection. 
``crispr`` integrates `CRT <http://www.room220.com/crt/>`_ and compare genome 
annotation to remove false positive. Let's try it.

::

     $ annogesic crispr \
        -g ANNOgesic/input/references/annotations/NC_009839.1.gff \
        -f ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        -pj ANNOgesic

Output are as following, ``CRT_results`` stores output of `CRT <http://www.room220.com/crt/>`_; 
``gffs`` stores gff files of the CRISPRs; ``statistics`` is for statistic files.

::

     $ ls ANNOgesic/output/crisprs/
     CRT_results  gffs  statistics
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
you assigned. The relationship between all features can be revealed.

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
       -a ANNOgesic/output/transcripts/gffs/NC_009839.1_transcript.gff \
       -of $ALL_FEATURES \
       -op NC_009839.1 \
       -pj ANNOgesic

Output gff file is stored in ``merge_all_features``

::

    $ ls ANNOgesic/output/merge_all_features/
    NC_009839.1_merge_features.gff

Producing the screenshots
-------------------------

It is a good idea if we can get screenshots of our interesting features. Then we can 
check them very quickly. Therefore, ANNOgesic provides a subcommand ``screenshot`` for 
generating screenshots.

Before we running it, we have to install `IGV <https://www.broadinstitute.org/software/igv/home>`_.

For testing, we use TSSs as main feature, sRNAs and CDSs as side features.

::

    $ annogesic screenshot \
        -mg ANNOgesic/output/TSSs/gffs/NC_009839.1_TSS.gff \
        -sg ANNOgesic/input/references/annotations/NC_009839.1.gff \
            ANNOgesic/output/sRNAs/gffs/best_candidates/NC_009839.1_sRNA.gff \
        -f ANNOgesic/input/references/fasta_files/NC_009839.1.fa \
        -o ANNOgesic/output/TSSs \
        -tl $TEX_LIBS \
        -pj ANNOgesic

Two txt files and two folders will be generated.

::

    $ ls ANNOgesic/output/TSSs/screenshots/NC_009839.1/
    forward/     forward.txt  reverse/     reverse.txt

``forward.txt`` and ``reverse.txt`` are batch files for running in `IGV <https://www.broadinstitute.org/software/igv/home>`_.
``forward`` and ``reverse`` are the folders for storing screenshots.

Since there are numerous candidates, we can only generate several ones in order to reduce the running time for testing.

::

    head -n 30 ANNOgesic/output/TSSs/screenshots/NC_009839.1/forward.txt > ANNOgesic/output/TSSs/screenshots/NC_009839.1/forward_6_cases.txt
    head -n 30 ANNOgesic/output/TSSs/screenshots/NC_009839.1/reverse.txt > ANNOgesic/output/TSSs/screenshots/NC_009839.1/reverse_6_cases.txt


Now, please open `IGV <https://www.broadinstitute.org/software/igv/home>`_ and follow the procedures: Tools -> 
Run Batch Script -> choose ``forward_6_cases.txt``. Once it is done, please do it again for reverse strand: Tools ->
Run Batch Script -> choose ``reverse_6_cases.txt``. If you just want to generate the screenshots for all candidates, 
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

If we have numerous samples and we want to check TSSs, Distinguishing the 
tracks of TEX+ and TEX- will be painful. Therefore, we provide a subcommand ``color_png`` to color
our screenshots.

::

    $ annogesic color_png \
        -t 2 \
        -f ANNOgesic/output/TSSs \
        -pj ANNOgesic

We will see output filenames are the same as before. However, when we open the figures, the tracks are colored.

::

    $ ls ANNOgesic/output/TSSs/screenshots/NC_009839.1/forward
    NC_009839.1:1396-1396.png  NC_009839.1:14812-14812.png  NC_009839.1:6676-6676.png  NC_009839.1:6680-6680.png  NC_009839.1:8098-8098.png  NC_009839.1:9295-9295.png
    $ ls ANNOgesic/output/TSSs/screenshots/NC_009839.1/reverse
    NC_009839.1:15670-15670.png  NC_009839.1:18053-18053.png  NC_009839.1:18360-18360.png  NC_009839.1:2199-2199.png  NC_009839.1:4463-4463.png  NC_009839.1:856-856.png

Now we already finished the first wonderful trip of ANNOgesic. Hopefully, you enjoy it!
