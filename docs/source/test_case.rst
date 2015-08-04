Performing a test cast
===================

Here we will guide you through a small test case using ANNOgesic. 
You will see most subcommands of ANNOgesic. The test case is a publicly 
available RNA-Seq from NCBI GEO that was part of a work by
`Bischler et al. <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67564>`_.
The differential RNA-seq data is of Helicobacter pylori 26695. 
There will be several output files which are generated in different formats. 
The CSV (tabular separated plain text files) files (opened by LibreOffice or Excel). 
For GFF3 file you can use a genome browser for example `IGB <http://bioviz.org/igb/index.html>`_, 
`IGV <https://www.broadinstitute.org/igv/>`_ or `Jbrowse <http://jbrowse.org/>`_.
There are also several PNG and TXT files.

Before we start, please refer to the section of ``The format of filename`` and 
``The format of libraries for import to ANNOgesic`` in 
`subcommand <https://github.com/Sung-Huan/ANNOgesic/blob/master/docs/source/subcommands.rst>`_. 
All the details are also in `subcommand <https://github.com/Sung-Huan/ANNOgesic/blob/master/docs/source/subcommands.rst>`_.

Generating a project
--------------------

First of all, we need to create a working folder, we use ``create`` subcommand.

::

    annogesic create ANNOgesic

Then you will see 

::

    Created folder "ANNOgesic" and required subfolders.
    $ ls 
    ANNOgesic

Retrieving the input data
-------------------

For our test case, we can download from `NCBI <ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Helicobacter_pylori_26695_uid57787/>`_.
We can set the ``$FTP_SOURCE`` first

::

    FTP_SOURCE=ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Helicobacter_pylori_26695_uid57787

Then download fasta files(``-f``), gff files(``-g``), gbk files(``-k``), ptt files(``-p``), 
rnt files(``-r``), and transfer to embl(``-e``).

::

    $ annogesic get_input_files -F $FTP_SOURCE -g -f -e -k -p -r ANNOgesic

you can also write 

::

    $ annogesic get_input_files -F $FTP_SOURCE -gfekpr ANNOgesic

Then you will get the following results

::

    $ ls ANNOgesic/input/reference/fasta/
    NC_000915.fa
    $ ls ANNOgesic/input/reference/annotation/
    NC_000915.embl  NC_000915.gbk  NC_000915.gff  NC_000915.ptt  NC_000915.rnt

Putting wig, bam files and reads to proper location
------------------
For the test case, you can download files from 
`here <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67564>`_.

::

    $ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP056%2FSRP056856/SRR1951997/SRR1951997.sra
    $ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP056%2FSRP056856/SRR1951998/SRR1951998.sra

Then you need to convert the SRA files to Fasta or Fastq format for mapping. You can 
use `SRA toolkit <http://www.ncbi.nlm.nih.gov/books/NBK158900/>`_ for that.

::
  
   $ wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.5.2/sratoolkit.2.5.2-ubuntu64.tar.gz
   $ tar -zxvf sratoolkit.2.5.2-ubuntu64.tar.gz
   $ rm sratoolkit.2.5.2-ubuntu64.tar.gz
   $ ./sratoolkit.2.5.2-ubuntu64/bin/fastq-dump.2.5.2 --fasta SRR1951997.sra
   $ ./sratoolkit.2.5.2-ubuntu64/bin/fastq-dump.2.5.2 --fasta SRR1951998.sra
   $ mv *.fasta ANNOgesic/input/reads
   $ rm SRR1951997.sra
   $ rm SRR1951998.sra

Now you get the reads. You can select a part of them to reduce the running time.

Then we have to download the wiggle files.

::

    wget -cP ANNOgesic/input/wigs/tex_notex ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1649nnn/GSM1649587/suppl/GSM1649587%5FHp26695%5FML%5FB1%5FHS1%5F%2DTEX%5Fforward%2Ewig%2Egz
    wget -cP ANNOgesic/input/wigs/tex_notex ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1649nnn/GSM1649587/suppl/GSM1649587%5FHp26695%5FML%5FB1%5FHS1%5F%2DTEX%5Freverse%2Ewig%2Egz
    wget -cP ANNOgesic/input/wigs/tex_notex ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1649nnn/GSM1649588/suppl/GSM1649588%5FHp26695%5FML%5FB1%5FHS1%5F%2BTEX%5Fforward%2Ewig%2Egz
    wget -cP ANNOgesic/input/wigs/tex_notex ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1649nnn/GSM1649588/suppl/GSM1649588%5FHp26695%5FML%5FB1%5FHS1%5F%2BTEX%5Freverse%2Ewig%2Egz
    cd ANNOgesic/input/wigs/tex_notex
    gunzip GSM1649587_Hp26695_ML_B1_HS1_-TEX_forward.wig.gz \
           GSM1649587_Hp26695_ML_B1_HS1_-TEX_reverse.wig.gz \
           GSM1649588_Hp26695_ML_B1_HS1_+TEX_forward.wig.gz \
           GSM1649588_Hp26695_ML_B1_HS1_+TEX_reverse.wig.gz
    cd ../../../../

Improving the reference genome
------------------

If the input data which we retrieved is exactly what you want, you can skip this step. 
You just need to copy the data which you retrieved to ``ANNOgesic/output/target/fasta``. 

If the retrieved data is only the closed strain of your target strain, you may need to run 
subcommand ``get_target_fasta``. However, this command need the mutation table. please refer 
to the user manual. Once you have the mutation table, you can improve the fasta files.

We use a simple example to modify our test case. The mutation table is 

::

    #target_id	reference_id	reference_nt	position	target_nt	impact of correction	locus tag	gene	Description
    NC_test.1	NC_000915.1	a	3	c		SAOUHSC_00002	dnaA	XXXXXX
    NC_test.1	NC_000915.1	t	6	-	deletion			YYYYYY
    test_case2	NC_000915.1	-	6	g	insertion	SAOUHSC_00132		

You can see the new strain will be NC_test.1 and test_case2. Although there are four fasta files 
in ``ANNOgesic/input/reference/fasta``, we just modify two of them. Therefore, there will be 
only two fasta files in ``ANNOgesic/output/target/fasta``, too.

Now, let's try it

::

     $ annogesic get_target_fasta \
        -r ANNOgesic/input/reference/fasta \
        -o test_case1:NC_test.1 test_case2:test_case2 \
        -m ANNOgesic/input/mutation_table/mutation.csv \
        ANNOgesic

``-r`` is the folder of original fasta files. In ``-o`` you can assign the filename of output fasta and 
which strains you want to put in it. In our case, we call the first fasta file test_case1 and the 
second one test_case2. In test_case1 stores the fasta of NC_test.1 and test_case2 stores test_case2. 
Now we can check the retuls.

::

    $ head ANNOgesic/input/reference/fasta/NC_016810.fa
    >NC_000915.1
    TGATTAGTGATTAGTGATTAGTGATTAGTGATTAGTGATTAGTGATTAGTGATTAGTGATTAGTGATTAG
    TGATTAGTGATTAGTGATTAGTGATTAGTGATTAGTGATTAGTGATTAGTGATTAGTGATTAGTGATTAG
    TGATTAGTGATTAGTGATTAGTGATTAGTGATTAGTGATTATAGCATCATTTTTTAAATTTAGGTATAAA
    ACACCCTCAATTCAAGGGTTTTTGAGTGAGCTTTTTGCTCAAAGAATCCAAGATAGCGTTTAAAAATTTA
    GGGGTGTTAGGCTCAGCGTAGAGTTTGCCAAGCTCTATGCATTCATTGATGATGATAGGGTTTTGCGTGG
    GCGTGAAGCCAATTTCATACGCTCCTAAGCGTAAAATCGCCTTTTCCATGCTCCCTAATCGCTTGAAATC
    CCAGTCTTTTAAATGCGGCTCGATGAGGGCGTCAATTTCATTGATTTTTTCTAACACGCCATTAAAAAGG
    CTTAAAGCGAAAGCGAGTTGGTTGTTTTTAATCTTTTTTTCTTCTAACATGCTAGAAGCGATTTTTTTAA
    TTTCTTCATTACCGCTCTCAAACGCATACAACAATTCAACCACAGCCCCCCTGGCTTGAGTTCGTGTCGC
    $ head ANNOgesic/output/target/fasta/test_case1.fa
    >NC_test.1
    TGcTTGTGATTAGTGATTAGTGATTAGTGATTAGTGATTAGTGATTAGTGATTAGTGATT
    AGTGATTAGTGATTAGTGATTAGTGATTAGTGATTAGTGATTAGTGATTAGTGATTAGTG
    ATTAGTGATTAGTGATTAGTGATTAGTGATTAGTGATTAGTGATTAGTGATTAGTGATTA
    TAGCATCATTTTTTAAATTTAGGTATAAAACACCCTCAATTCAAGGGTTTTTGAGTGAGC
    TTTTTGCTCAAAGAATCCAAGATAGCGTTTAAAAATTTAGGGGTGTTAGGCTCAGCGTAG
    AGTTTGCCAAGCTCTATGCATTCATTGATGATGATAGGGTTTTGCGTGGGCGTGAAGCCA
    ATTTCATACGCTCCTAAGCGTAAAATCGCCTTTTCCATGCTCCCTAATCGCTTGAAATCC
    CAGTCTTTTAAATGCGGCTCGATGAGGGCGTCAATTTCATTGATTTTTTCTAACACGCCA
    TTAAAAAGGCTTAAAGCGAAAGCGAGTTGGTTGTTTTTAATCTTTTTTTCTTCTAACATG

In new fasta file, the third nucleotide replace from A to c. Moreover, The sixth nucleotide is deleted.
In test_case2, it also modify the fasta file based on the mutation table.

when you have the correct fasta files, you can map your reads to the correct fasta file.

If you have no mutation table, you can also use the subcommand ``snp`` to detect the mutation automatically. 
For this subcommand, we will go through it later.

Generating annotation files
-------------------

We have the fasta files now. We can use it to generate our annotation files. If the annotaion files which 
you retrieved by ``get_input_files`` is exactly what you want, you can skip this step. You just need to 
copy the annotation files to ``ANNOgesic/output/target/annotation``.

Before you running it, you have to notice path. If you are using docker container, the path is alread setup.
You don't need to worry about it. However, if you are build ANNOgesic by yourself, remember to assign the 
path by running

::

    . $PAGIT_HOME/sourceme.pagit

``$PAGIT_HOME`` is the directory of PAGIT. The better way is change the environment. Or you have to run 
this command everytime. For changing the environment, you just need to copy all the information of 
``$PAGIT_HOME/sourceme.pagit`` to ``.bashrc``.

Now, we can try it.

::

    python3 annogesic annotation_transfer \
        -re ANNOgesic/input/reference/annotation \
        -rf ANNOgesic/input/reference/fasta \
        -tf ANNOgesic/output/target/fasta \
        -e chromosome \
        -t Strain \
        -p NC_000915:NC_test.1 NC_000915:test_case2 \
        -g \
        ANNOgesic

``-e`` is the prefix of output embl files. ``-t`` is a program of `RATT <http://ratt.sourceforge.net/>`_.
We use ``Strain`` because the similarity is higher than 90%. We assign the pairs of transfer at ``-p``. 
The names are strains' name not filenames of fasta files. ``-g`` means we want to transfer the embl files 
to GFF3 files and store in ``ANNOgesic/output/target/annotation``.

When the computation is done, you can see

::

    $ ls ANNOgesic/output/target/annotation/
    test_case1.gff  test_case1.ptt  test_case1.rnt  test_case2.gff  test_case2.ptt  test_case2.rnt
    $ ls ANNOgesic/output/annotation_transfer/
    chromosome.NC_test.1.final.embl  chromosome.test_case2.final.embl  NC_test.1.gff  ratt_log.txt  test_case2.gff

In ``ANNOgesic/output/target/annotation``, you can see ptt, rnt and gff files. In ``ANNOgesic/output/annotation_transfer``,
you can see the results of `RATT <http://ratt.sourceforge.net/>`_. ``chromosome.NC_test.1.final.embl`` and 
``chromosome.test_case2.final.embl`` are generated by `RATT <http://ratt.sourceforge.net/>`_. Gff files are 
transferred from these embl files.

Gene expression analysis
-----------------

Now we already saw how to generate the fasta and annotation files. In order to running
the following subcommand, we have to replace the files of ``ANNOgesic/output/target/annotation``
and ``ANNOgesic/output/target/fasta``.

:: 

    rm -rf ANNOgesic/output/target/annotation
    cp -rf ANNOgesic/input/reference/annotation ANNOgesic/output/target
    rm -rf ANNOgesic/output/target/fasta
    cp -rf ANNOgesic/input/reference/fasta ANNOgesic/output/target

We do this because the wiggle file is generated by NC_000915. If you running by your data, you don't need to do this. 
Now, our fasta and annotation files are fit with wiggle files. We can other subcommands now.

For analyzing gene expression, we can run ``expression analysis``. Based on this subcommand, you can 
know which CDS expresses in which library or discover housekeeping gene.

Before running it, you can set the library first.

::

    tex_notex_libs="GSM1649587_Hp26695_ML_B1_HS1_-TEX_forward.wig:notex:1:a:+ \
                    GSM1649587_Hp26695_ML_B1_HS1_-TEX_reverse.wig:notex:1:a:- \
                    GSM1649588_Hp26695_ML_B1_HS1_+TEX_forward.wig:tex:1:a:+ \
                    GSM1649588_Hp26695_ML_B1_HS1_+TEX_reverse.wig:tex:1:a:-"

The command would be like

::

    python3 annogesic expression_analysis \
        -g $ANNOGESIC_FOLDER/output/target/annotation \
        -tl $tex_notex_libs \
        -tw $ANNOGESIC_FOLDER/input/wigs/tex_notex \
        -f CDS \
        -rt 1 \
        $ANNOGESIC_FOLDER

It will generate several gff files and statistics files.

::

    $ ls ANNOgesic/output/target/annotation/for_libs/gffs
    NC_000915_CDS_1_texnotex.gff  NC_000915_CDS_all_libs.gff  NC_000915_CDS_at_least_one_lib.gff  NC_000915_CDS_no_express.gff
    $ ls ANNOgesic/output/target/annotation/for_libs/statistics
    NC_000915_CDS.csv

``NC_000915_CDS_1_texnotex.gff`` stores the CDS which express in the library of 1_texnotex. In our definition 
of library, it is ``GSM1649587_Hp26695_ML_B1_HS1``. ``NC_000915_CDS_all_lib.gff`` is the CDS which express in all libraries.
``NC_000915_CDS_at_least_one_lib.gff`` stores the CDS which express at least one library. ``NC_000915_CDS_no_express.gff`` 
is the CDS which don't express in any libraries.

TSS and processing site prediction and optimization
-----------------

Before running ``tsspredator``, if you want to use the optimized parameters, 
you need to run ``optimize_TSSpredator`` first.

Then need to manual check some TSS. In our experience, we recommand you detect more than 50 TSS and longer than 100kb. 
For test case, we can just copy the default parameter as a manual detected ones. This is only for test. If you are 
running your own data and you have manual detected TSS, please use it.

::

    python3 annogesic tsspredator \
        -w ANNOgesic/input/wigs/tex_notex \
        -f ANNOgesic/output/target/fasta \
        -g ANNOgesic/output/target/annotation \
        -l $tex_notex_libs \
        -p test \
        ANNOgesic
    mv ANNOgesic/output/TSS/gffs/NC_000915_TSS.gff ANNOgesic/input/manual_TSS

Now, we have a fake manual detected TSS file. we can try optimization of TSS right now.

::

    python3 annogesic optimize_tsspredator \
        -w ANNOgesic/input/wigs/tex_notex \
        -fs ANNOgesic/output/target/fasta \
        -g ANNOgesic/output/target/annotation \
        -n NC_000915.1 \
        -l $tex_notex_libs \
        -p TSS -s 25 \
        -m ANNOgesic/input/manual_TSS/NC_000915_TSS.gff  \
        ANNOgesic
    ...
    Current Parameter:step=16_height=0.5_reduction_height=0.1_factor=8.3_reduction_factor=2.4_base_height=0.068_enrichment_factor=3.0_processing_factor=2.4
    Current:TP=665  TP_rate=0.38394919168591224     FP=23   FP_rate=1.3804403604749916e-05  FN=1067 Missing_ratio=0.6160508083140878
    Best Parameter:height=0.5       reduction_height=0.1    factor=8.3      reduction_factor=2.4    base_height=0.068       enrichment_factor=3.0   processing_factor=2.4
    Best:TP=665     TP_rate=0.38394919168591224     FP=23   FP_rate=1.3804403604749916e-05  FN=1067 Missing_ratio=0.6160508083140878
    Current Parameter:step=17_height=0.5_reduction_height=0.1_factor=9.6_reduction_factor=6.4_base_height=0.1_enrichment_factor=3.0_processing_factor=1.0
    Current:TP=539  TP_rate=0.3112009237875289      FP=13   FP_rate=7.802488993989082e-06   FN=1193 Missing_ratio=0.6887990762124712
    Best Parameter:height=0.5       reduction_height=0.1    factor=8.3      reduction_factor=2.4    base_height=0.068       enrichment_factor=3.0   processing_factor=2.4
    Best:TP=665     TP_rate=0.38394919168591224     FP=23   FP_rate=1.3804403604749916e-05  FN=1067 Missing_ratio=0.6160508083140878
    Current Parameter:step=18_height=1.4_reduction_height=0.1_factor=8.3_reduction_factor=2.8_base_height=0.1_enrichment_factor=3.0_processing_factor=1.0
    Current:TP=433  TP_rate=0.25    FP=11   FP_rate=6.602106071836916e-06   FN=1299 Missing_ratio=0.75
    Best Parameter:height=0.5       reduction_height=0.1    factor=8.3      reduction_factor=2.4    base_height=0.068       enrichment_factor=3.0   processing_factor=2.4
    Best:TP=665     TP_rate=0.38394919168591224     FP=23   FP_rate=1.3804403604749916e-05  FN=1067 Missing_ratio=0.6160508083140878
    Current Parameter:step=19_height=0.5_reduction_height=0.1_factor=8.3_reduction_factor=4.9_base_height=0.1_enrichment_factor=3.0_processing_factor=5.6
    Current:TP=582  TP_rate=0.33602771362586603     FP=16   FP_rate=9.603063377217333e-06   FN=1150 Missing_ratio=0.663972286374134
    Best Parameter:height=0.5       reduction_height=0.1    factor=8.3      reduction_factor=2.4    base_height=0.068       enrichment_factor=3.0   processing_factor=2.4
    Best:TP=665     TP_rate=0.38394919168591224     FP=23   FP_rate=1.3804403604749916e-05  FN=1067 Missing_ratio=0.6160508083140878
    ...

``optimize_TSSpredator`` will compare gff files of manual checked TSS and predicted TSS to find the best parameters. 
You can check the results and parameters of each step in screen. Because we just test it, we set the steps only 20. 
When the program finished, you can find several files.

::

    $ ls ANNOgesic/output/TSS/optimized_TSSpredator/
    best.csv  log.txt  stat.csv

``best.csv`` is for the best parameters; ``stat.csv`` is for the parameters of each step.

Let's assume the best parameters are height is 0.4, height_reduction is 0.2, factor is 2.0, factor_reduction is 0.5, 
base_height is 0.0, enrichment_factor is 2.0, processing_factor is 1.5. Now we can set the parameter set for running  
``tss``.

::

    python3 annogesic tsspredator \
        -w ANNOgesic/input/wigs/tex_notex \
        -f ANNOgesic/output/target/fasta \
        -g ANNOgesic/output/target/annotation \
        -l $tex_notex_libs \
        -p test \
        -he 0.4 \
        -rh 0.2 \
        -fa 2.0 \
        -rf 0.5 \
        -bh 0.0 \
        -ef 2.0 \
        -pf 1.5 \
        -s \
        -v \
        -m ANNOgesic/input/manual_TSS/NC_000915_TSS.gff \
        ANNOgesic

If you set the manual checked TSS ``-m``, the subcommand will merge the manual checked TSS with predicted TSS. 
If you didn't set it, the subcommand will only produce predicted TSS. You will get gff file, MasterTable and statistic file.

::

    $ ls ANNOgesic/output/TSS/gffs/
    NC_000915_TSS.gff
    $ ls ANNOgesic/output/TSS/MasterTables/MasterTable_NC_000915.1/
    AlignmentStatistics.tsv  err.txt  log.txt  MasterTable.tsv  superConsensus.fa  superTSS.gff  superTSStracks.gff  test_super.fa  test_super.gff  test_TSS.gff  TSSstatistics.tsv
    $ ls ANNOgesic/output/TSS/statistics/NC_000915/
    stat_compare_TSSpredator_manual_NC_000915.csv  stat_gene_vali_NC_000915.csv  stat_TSS_class_NC_000915.csv  stat_TSS_libs_NC_000915.csv  TSS_class_NC_000915.1.png  TSS_venn_NC_000915.1.png

If you want to predict processing site, the procedures are the same. You just need to change the program from TSS to 
processing_site (``-t``).

Performing transcript assembly
----------------

For detecting transcript boundary, transcript assembly is the basic feature. 
we can use the subcommand ``transcript_assembly`` to do it. Normally, we strongly 
recommand you to get fragmented libraries for ``transcript_assembly``. However, there is no 
fragmented libraries available in the database. Therefore, we only use TEX +/- to do it.

The command would be like the following.

::

    python3 annogesic transcript_assembly \
        -g ANNOgesic/output/target/annotation \
        -nw ANNOgesic/input/wigs/tex_notex \
        -tl $tex_notex_libs \
        -rt 1 \
        -te 2 \
        -ct ANNOgesic/output/TSS/gffs \
        -cg ANNOgesic/output/target/annotation \
        ANNOgesic

It will generate gff files. Because we also compared with TSS and annotation files, it will generate statistics files.

::

    $ ls ANNOgesic/output/transcriptome_assembly/gffs
    NC_000915_transcript.gff
    $ ls ANNOgesic/output/transcriptome_assembly/statistics
    stat_compare_Transcriptome_assembly_CDS_NC_000915.csv  stat_compare_Transcriptome_assembly_TSS_NC_000915.csv

Prediction of terminator
----------------------

For predicting terminator, we can use subcommand ``terminator``. The command is like following.

::

    python3 annogesic terminator \
        -f $ANNOGESIC_FOLDER/output/target/fasta \
        -g $ANNOGESIC_FOLDER/output/target/annotation \
        -s \
        -tw $ANNOGESIC_FOLDER/input/wigs/tex_notex \
        -a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs \
        -tl $tex_notex_libs \
        -te 2 -rt 1 -tb \
        ANNOgesic

It will generate three different kinds of gff files and tables.

::

    $ ls ANNOgesic/output/terminator/gffs/
    all_candidates  detect  express
    $ ls ANNOgesic/output/terminator/tables
    all_candidates  detect  express

``all_candidates`` is for all candidates; ``express`` is for the candidates which have expression; 
``detect`` is for the candidates which have dramatic decreasing coverage. There is a gff file or table for 
each folder.

::

    $ ls ANNOgesic/output/terminator/gffs/detect
    NC_000915_term.gff
    $ ls ANNOgesic/output/terminator/gffs/express
    NC_000915_term.gff
    $ ls ANNOgesic/output/terminator/gffs/all_candidates
    NC_000915_term.gff
    $ ls ANNOgesic/output/terminator/tables/detect
    NC_000915_term.csv
    $ ls ANNOgesic/output/terminator/tables/express
    NC_000915_term.csv
    $ ls ANNOgesic/output/terminator/tables/all_candidates
    NC_000915_term.csv

In transtermhp folder, there are the output files from `TranstermHP <http://transterm.cbcb.umd.edu/>`_.

::

    $ ls ANNOgesic/output/terminator/transtermhp/NC_000915.1
    NC_000915.1_best_terminator_after_gene.bag  NC_000915.1_terminators.txt  NC_000915.1_terminators_within_robust_tail-to-tail_regions.t2t

Moreover, the statistics files are stored in ``statistics``.

::

    $ ls ANNOgesic/output/terminator/statistics/
    stat_NC_000915.csv

Generating UTR
--------------

Now, we have the information of TSS, transcripts and terminators. We can detect the 5'UTR and 3'UTR easily by using 
the subcommand ``utr``.

::

    python3 annogesic utr \
        -g ANNOgesic/output/target/annotation \
        -t ANNOgesic/output/TSS/gffs \
        -a ANNOgesic/output/transcriptome_assembly/gffs \
        -e ANNOgesic/output/terminator/gffs/detect \
        ANNOgesic

If your TSS is not generated by ANNOgesic, please assign ``-s``, it will classify the TSS to generate UTR.
The gff files and statistics files will be stored in ``5UTR`` and ``3UTR``.

::

    $ ls ANNOgesic/output/UTR/3UTR
    gffs/       statistics/
    $ ls ANNOgesic/output/UTR/5UTR
    gffs/       statistics/
    $ ls ANNOgesic/output/UTR/3UTR/gffs
    NC_000915_3UTR.gff
    $ ls ANNOgesic/output/UTR/5UTR/gffs
    NC_000915_5UTR.gff
    $ ls ANNOgesic/output/UTR/5UTR/statistics
    NC_000915_all_5utr_length.png
    $ ls ANNOgesic/output/UTR/3UTR/statistics
    NC_000915_all_3utr_length.png

Until now, you have all information for defining the transcript boundary.

Integrating to operon and suboperon
-----------------

Now, we already had TSS, transcript, terminator, CDS, UTR. We can integrate all these information to 
detect operons and suboperons. You can use the subcommand ``operon`` to get it.

::

    python3 annogesic operon \
        -g ANNOgesic/output/target/annotation \
        -t ANNOgesic/output/TSS/gffs \
        -a ANNOgesic/output/transcriptome_assembly/gffs \
        -u5 ANNOgesic/output/UTR/5UTR/gffs \
        -u3 ANNOgesic/output/UTR/3UTR/gffs \
        -e ANNOgesic/output/terminator/gffs/detect \
        -s -c \
        ANNOgesic

``operon`` will generate three folders to store gff files, tables and statistics files.

::

    $ ls ANNOgesic/output/operons/
    gffs  statistics  tables
    $ ls ANNOgesic/output/operons/gffs/
    NC_000915_all_features.gff
    $ ls ANNOgesic/output/operons/tables/
    operon_NC_000915.csv
    $ ls ANNOgesic/output/operons/statistics/
    stat_operon_NC_000915.csv

Prediction of sRNA and sORF
-----------------

Based on comparison of trascripts and CDS information, we can detect intergenic sRNA. Moreover, we 
have the information of TSS and processing sites. We can detect UTR derived sRNA, too. You can 
get sRNA by running subcommand ``srna``. Normally, we would recommand you use fragmented libraries to 
compute ``srna``. However, we can't get it right now. We use TEX +/- for this test case.

Before ŕunning ``srna``, we have to get `BSRD <http://www.bac-srna.org/BSRD/index.jsp>`_ and 
`nr database <ftp://ftp.ncbi.nih.gov/blast/db/FASTA/>`_. Please download it and store in ``ANNOgesic/input/database``.
If you want to store to other folders, please change ``-sd`` and ``-od``. Then running the following 
command.

::

    python3 annogesic srna \
        -d 1 2 3 4 \
        -g ANNOgesic/output/target/annotation \
        -t ANNOgesic/output/TSS/gffs \
        -p ANNOgesic/output/processing_site/gffs \
        -a ANNOgesic/output/transcriptome_assembly/gffs \
        -tw ANNOgesic/input/wigs/tex_notex \
        -f ANNOgesic/output/target/fasta \
        -m -u -fd \
        -sd ANNOgesic/input/database/sRNA_database \
        -nd ANNOgesic/input/database/nr \
        -tl $tex_notex_libs \
        -te 2 \
        -rt 1 \
        -ba -sb \
        ANNOgesic

Because we have no information of sORF right now, we can't compare sRNA and sORF. If you already have 
the information of sORF, you can assign ``-d 1 2 3 4 5`` and ``-O`` for the path of sORF. It will compare 
sORF ans sRNA.

The output of ``srna`` will be the following.

::

    $ ls ANNOgesic/output/sRNA/
    blast_result_and_misc  gffs  log.txt  mountain_plot  sec_structure  sRNA_2d_NC_000915  sRNA_seq_NC_000915  statistics  tables

``blast_result_and_misc`` will store the results of blast; ``mountain_plos`` will store the mountain plots; 
``sec_structure`` will store the plots of secondary structure of sRNA; ``statistics`` will store statistics files.

``sRNA_2d_NC_000915`` and ``sRNA_seq_NC_000915`` are text file of sequence of sRNA and secondary structure of sRNA.

::

    $ ls ANNOgesic/output/sRNA/blast_result_and_misc/
    nr_blast_NC_000915.txt  sRNA_blast_NC_000915.txt
    $ ls ANNOgesic/output/sRNA/mountain_plot/NC_000915/
    srna0_NC_000915.1_16651_16765_-_mountain.pdf        srna138_NC_000915.1_1496938_1497216_-_mountain.pdf  srna25_NC_000915.1_444979_445114_+_mountain.pdf  srna63_NC_000915.1_761094_761174_+_mountain.pdf
    srna100_NC_000915.1_1164252_1164295_-_mountain.pdf  srna139_NC_000915.1_1502542_1502637_+_mountain.pdf  srna26_NC_000915.1_445012_445139_-_mountain.pdf  srna64_NC_000915.1_773130_773564_+_mountain.pdf
    ...
    ls ANNOgesic/output/sRNA/sec_structure/dot_plot/NC_000915/
    srna0_NC_000915.1_16651_16765_-_dp.pdf        srna138_NC_000915.1_1496938_1497216_-_dp.pdf  srna25_NC_000915.1_444979_445114_+_dp.pdf  srna63_NC_000915.1_761094_761174_+_dp.pdf
    srna100_NC_000915.1_1164252_1164295_-_dp.pdf  srna139_NC_000915.1_1502542_1502637_+_dp.pdf  srna26_NC_000915.1_445012_445139_-_dp.pdf  srna64_NC_000915.1_773130_773564_+_dp.pdf
    ...
    $ ls ANNOgesic/output/sRNA/sec_structure/sec_plot/NC_000915/
    srna0_NC_000915.1_16651_16765_-_rss.pdf        srna138_NC_000915.1_1496938_1497216_-_rss.pdf  srna25_NC_000915.1_444979_445114_+_rss.pdf  srna63_NC_000915.1_761094_761174_+_rss.pdf
    srna100_NC_000915.1_1164252_1164295_-_rss.pdf  srna139_NC_000915.1_1502542_1502637_+_rss.pdf  srna26_NC_000915.1_445012_445139_-_rss.pdf  srna64_NC_000915.1_773130_773564_+_rss.pdf
    ...
    $ ls ANNOgesic/output/sRNA/statistics/
    stat_sRNA_blast_class_NC_000915.csv  stat_sRNA_class_NC_000915.csv

For ``gffs`` and ``tables``, they are divided by three kinds of results. ``all_candidates`` is for all candidates 
without filtering; ``best`` is for the best candidates of sRNA; ``for_class`` is for each class of sRNA which you 
imported before. For our test case, we import TSS(class 1), folding energy(class 2), blast to nr(class 3), 
blast to sRNA database(class 4 for without homologs, class 5 for with homologs).

::

    $ ls ANNOgesic/output/sRNA/gffs/
    all_candidates  best  for_class
    $ ls ANNOgesic/output/sRNA/tables/
    all_candidates  best  for_class
    $ ls ANNOgesic/output/sRNA/gffs/all_candidates/
    NC_000915_sRNA.gff
    $ ls ANNOgesic/output/sRNA/tables/all_candidates/
    NC_000915_sRNA.csv
    $ ls ANNOgesic/output/sRNA/gffs/best/
    NC_000915_sRNA.gff
    $ ls ANNOgesic/output/sRNA/tables/best/
    NC_000915_sRNA.csv
    $ ls ANNOgesic/output/sRNA/gffs/for_class/NC_000915/
    class_1_all.gff                          class_1_class_2_class_3_class_5_all.gff  class_1_class_3_class_4_all.gff  class_2_all.gff                  class_2_class_4_all.gff  class_3_class_5_all.gff
    class_1_class_2_all.gff                  class_1_class_2_class_4_all.gff          class_1_class_3_class_5_all.gff  class_2_class_3_all.gff          class_2_class_5_all.gff  class_4_all.gff
    class_1_class_2_class_3_all.gff          class_1_class_2_class_5_all.gff          class_1_class_4_all.gff          class_2_class_3_class_4_all.gff  class_3_all.gff          class_5_all.gff
    class_1_class_2_class_3_class_4_all.gff  class_1_class_3_all.gff                  class_1_class_5_all.gff          class_2_class_3_class_5_all.gff  class_3_class_4_all.gff
    $ ls ANNOgesic/output/sRNA/tables/for_class/NC_000915/
    class_1_all.csv                          class_1_class_2_class_3_class_5_all.csv  class_1_class_3_class_4_all.csv  class_2_all.csv                  class_2_class_4_all.csv  class_3_class_5_all.csv
    class_1_class_2_all.csv                  class_1_class_2_class_4_all.csv          class_1_class_3_class_5_all.csv  class_2_class_3_all.csv          class_2_class_5_all.csv  class_4_all.csv
    class_1_class_2_class_3_all.csv          class_1_class_2_class_5_all.csv          class_1_class_4_all.csv          class_2_class_3_class_4_all.csv  class_3_all.csv          class_5_all.csv
    class_1_class_2_class_3_class_4_all.csv  class_1_class_3_all.csv                  class_1_class_5_all.csv          class_2_class_3_class_5_all.csv  class_3_class_4_all.csv    

As we know, the sORF also has no annotation information and has coverage. Therefore, the potential sRNA 
may be sORF not sRNA. Therefore, it is a good idea to have the information of sORF. You can use subcommand 
``sorf`` to get it.

::

    python3 annogesic sorf \
        -g ANNOgesic/output/target/annotation \
        -t ANNOgesic/output/TSS/gffs \
        -a ANNOgesic/output/transcriptome_assembly/gffs \
        -tw ANNOgesic/input/wigs/tex_notex \
        -f ANNOgesic/output/target/fasta \
        -s ANNOgesic/output/sRNA/gffs/best \
        -tl $tex_notex_libs \
        -te 2 -rt 1 -u \
        ANNOgesic

Then you can get the gff files, statistics files and tables. Based on the information of TSS and 
sRNA, the gff files and tables will be divided by ``all_candidates`` and ``best``.

::

    $ ls ANNOgesic/output/sORF/gffs/all_candidates/
    NC_000915_sORF.gff
    $ ls ANNOgesic/output/sORF/gffs/best/
    NC_000915_sORF.gff
    $ ls ANNOgesic/output/sORF/tables/all_candidates/
    NC_000915_sORF.csv
    $ ls ANNOgesic/output/sORF/tables/best/
    NC_000915_sORF.csv
    $ ls ANNOgesic/output/sORF/statistics/
    stat_NC_000915_sORF.csv

Performing sRNA target prediction
------------------

Now we have the sRNA candidates. If you want to know the targets of these sRNA, you can use ``srna_target`` 
to get these information.

::

    python3 annogesic srna_target \
        -g ANNOgesic/output/target/annotation \
        -f ANNOgesic/output/target/fasta \
        -r ANNOgesic/output/sRNA/gffs/best \
        -q NC_000915.1:+:13666:13701 \
        -p both \
        ANNOgesic

For testing, we just did the prediction to one sRNA. If you want to compute all sRNA, you can assign ``all`` 
to ``-q``. However, it may take several days.

``srna_target`` will generate several folders.

::

    $ ls ANNOgesic/output/sRNA_targets/
    merge  RNAplex  RNAup  sRNA_seqs  target_seqs

``sRNA_seqs`` and ``target_seqssRNA_seqs`` are for the sequences of sRNA and potential targets.

::

    $ ls ANNOgesic/output/sRNA_targets/sRNA_seqs
    NC_000915.1_sRNA.fa
    $ ls ANNOgesic/output/sRNA_targets/target_seqs
    NC_000915.1_target.fa

``RNAplex`` and ``RNAup`` are for the output of `RNAplex and RNAup <http://www.tbi.univie.ac.at/RNA/>`_.

::

    $ ls ANNOgesic/output/sRNA_targets/RNAplex/NC_000915.1/
    NC_000915.1_RNAplex_rank.csv  NC_000915.1_RNAplex.txt
    $ ls ANNOgesic/output/sRNA_targets/RNAup/NC_000915.1/
    NC_000915.1_RNAup.log  NC_000915.1_RNAup_rank.csv  NC_000915.1_RNAup.txt

``merge`` is for the merged results of `RNAplex and RNAup <http://www.tbi.univie.ac.at/RNA/>`_.

::

    $ ls ANNOgesic/output/sRNA_targets/merge/NC_000915.1/
    NC_000915.1_merge.csv  NC_000915.1_overlap.csv

Promoter motif detection
----------------

As long as you have TSS, you can use the subcommand ``promoter`` to get promoter. It will based on 
the class of TSS to generate the promoters. Therefore, if your TSS data is not computed by ``ANNOgesic``, 
you need to assign ``-s`` and ``-g`` (for annotation files). Then ``promoter`` will help you 
to classify your TSS data.

::

    python3 annogesic promoter \
        -t ANNOgesic/output/TSS/gffs \
        -f ANNOgesic/output/target/fasta \
        -w 50 2-10 \
        -p 10 \
        ANNOgesic

You can define the length of motif. In our test case, we use ``50`` and ``2-10``. ``2-10`` means the 
width is from 2 to 10.

Based on the class of TSS, it will generate different output files.

::

    $ ls ANNOgesic/output/promoter_analysis/NC_000915
    promoter_motifs_NC_000915_allstrain_all_types_2-10_nt  promoter_motifs_NC_000915_allstrain_internal_50_nt   promoter_motifs_NC_000915_allstrain_secondary_2-10_nt
    promoter_motifs_NC_000915_allstrain_all_types_50_nt    promoter_motifs_NC_000915_allstrain_orphan_2-10_nt   promoter_motifs_NC_000915_allstrain_secondary_50_nt
    promoter_motifs_NC_000915_allstrain_antisense_2-10_nt  promoter_motifs_NC_000915_allstrain_orphan_50_nt     promoter_motifs_NC_000915_allstrain_without_orphan_2-10_nt
    promoter_motifs_NC_000915_allstrain_antisense_50_nt    promoter_motifs_NC_000915_allstrain_primary_2-10_nt  promoter_motifs_NC_000915_allstrain_without_orphan_50_nt
    promoter_motifs_NC_000915_allstrain_internal_2-10_nt   promoter_motifs_NC_000915_allstrain_primary_50_nt

Mapping and detecting of circular RNA
-------------------

You may also be interested in circular RNA. The subcommand ``circrna`` can help you to get the information. 
it apply `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_ to detect circular RNA. Because 
we didn't map the reads of test case before, we can also do it by running ``circrna``. However, if 
you already mapped the reads by other tools. Or you mapped read by 
`Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_ without ``-S``. You may need re-mapping again.
If your mapping is generated by `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_ with ``-S``, 
then you can skip ``-a`` and assign the path of bam files to ``-nb`` or ``-fb``. It can reduce the 
running time.

::

     python3 annogesic circrna \
        -f ANNOgesic/output/target/fasta \
        -p 10 \
        -g ANNOgesic/output/target/annotation \
        -cg -a \
        ANNOgesic

``circrna`` will generate several folders.

::

    $ ls ANNOgesic/output/circRNA/
    circRNA_tables  gffs  segemehl_align  segemehl_splice  statistics

``segemehl_align`` and ``segemehl_splice`` are for the results of 
`Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_. ``segemehl_align`` stores the bam files of 
alignment and ``segemehl_splice`` stores the results of splice detection.

::

    $ ls ANNOgesic/output/circRNA/segemehl_align/NC_000915.1/
    SRR1951997_NC_000915.1.bam  SRR1951998_NC_000915.1.bam
    $ ls ANNOgesic/output/circRNA/segemehl_splice/NC_000915/
    splicesites_all.bed  transrealigned_all.bed    

The gff files, tables and statistics files are stored in the other folders.

::

    $ ls ANNOgesic/output/circRNA/gffs/NC_000915/
    NC_000915_circRNA_all.gff  NC_000915_circRNA_best.gff
    $ ls ANNOgesic/output/circRNA/circRNA_tables/NC_000915/
    circRNA_NC_000915.txt
    $ ls ANNOgesic/output/circRNA/statistics/
    stat_circRNA_NC_000915.csv

``NC_000915_circRNA_all.gff`` is for all circular RNA candidates and ``NC_000915_circRNA_best.gff`` is the 
circular RNA after filering by mapping ratio and comparison of CDS.

SNP calling
--------------

If you want to know the SNP or mutation of your sRNA-seq data, you can use ``snp`` to get it.
``snp`` is divided by two part. One part compares with the "reference strain" which is the
closed strain of our strain("target strain"). You can refer to the section of ``Retrieving the input data``.
Because you may not have time to check the mutation between "reference strain" and "target strain",
it is a good way to detect the mutations automatically. You just need to put your alignment files of 
mapping with "reference strain" into correct path. It will generate the potential sequence.
The other part is for detecting the mutations of the reads and "target strain". In this part, you 
can know the real mutations in "target strain". Therefore, you need to put the alignment files with 
"target strain" to the correct folder.

Before running the subcommand, we must have the bam files. Because we already generated them through 
running ``circrna``, we can just use them. However, please remember that the mapping function of 
``circrna`` is basic one. If you have other request, please do mapping by yourself.

For testing, we only discover the mutation for "target strain" because our mapping is based on the 
fasta of "target strain" - NC_000915. Therefore, we copy the bam files to ``BAMs_map_target``.

::

    cp ANNOgesic/output/circRNA/segemehl_align/NC_000915.1/SRR195199* ANNOgesic/input/BAMs/BAMs_map_target/tex_notex

Then we can run the subcommand.

::

    python3 annogesic snp \
        -t target \
        -p 1 2 3 \
        -tw ANNOgesic/input/BAMs/BAMs_map_target/tex_notex \
        -f ANNOgesic/output/target/fasta \
        ANNOgesic

If you want to compute for comparison of "reference strain" and "target strain", you just need to change 
``-t`` to reference and assign the correct bam files.

``snp`` will generate two folders, ``compare_reference`` is for comparison of "reference strain" and 
"target strain". ``validate_target`` is for real mutations of "target strain"

::

    $ ls ANNOgesic/output/SNP_calling/                                                                                                      
    compare_reference  validate_target

Becaues we run ``validate_target``, you can see there are several folders under ``validate_target``.

::

    $ ls ANNOgesic/output/SNP_calling/validate_target/
    seqs/            SNP_raw_outputs/ SNP_table/       statistics/

All folders are divided by three parts - ``extend_BAQ``, ``extend_BAQ`` and ``extend_BAQ``.

::

    $ ls ANNOgesic/output/SNP_calling/validate_target/seqs/
    extend_BAQ/  with_BAQ/    without_BAQ/

In ``seqs``, there are the potential sequences.

::

    $ ls ANNOgesic/output/SNP_calling/validate_target/seqs/with_BAQ/NC_000915.1/
    NC_000915.1_NC_000915.1_1_100.fa  NC_000915.1_NC_000915.1_1_13.fa   NC_000915.1_NC_000915.1_1_179.fa  NC_000915.1_NC_000915.1_1_217.fa  NC_000915.1_NC_000915.1_1_256.fa  NC_000915.1_NC_000915.1_1_63.fa
    NC_000915.1_NC_000915.1_1_101.fa  NC_000915.1_NC_000915.1_1_140.fa  NC_000915.1_NC_000915.1_1_17.fa   NC_000915.1_NC_000915.1_1_218.fa  NC_000915.1_NC_000915.1_1_25.fa   NC_000915.1_NC_000915.1_1_64.fa
    ....

``SNP_raw_outputs`` stores the output of `Samtools and Bcftools<https://github.com/samtools>`_. 
``SNP_table`` stores the results after filtering and the index of potential sequence.
``statistics`` stores the statistics and visualization files.

::

    $ ls ANNOgesic/output/SNP_calling/validate_target/SNP_raw_outputs/NC_000915.1/
    NC_000915.1_extend_BAQ.vcf  NC_000915.1_with_BAQ.vcf  NC_000915.1_without_BAQ.vcf
    $ ls ANNOgesic/output/SNP_calling/validate_target/SNP_table/NC_000915.1/
    NC_000915.1_extend_BAQ_depth_only.vcf     NC_000915.1_with_BAQ_depth_only.vcf     NC_000915.1_without_BAQ_depth_only.vcf
    NC_000915.1_extend_BAQ_depth_quality.vcf  NC_000915.1_with_BAQ_depth_quality.vcf  NC_000915.1_without_BAQ_depth_quality.vcf
    NC_000915.1_extend_BAQ_seq_reference.csv  NC_000915.1_with_BAQ_seq_reference.csv  NC_000915.1_without_BAQ_seq_reference.csv
    $ ls ANNOgesic/output/SNP_calling/validate_target/statistics/
    NC_000915.1_extend_BAQ_NC_000915.1_SNP_QUAL.png  NC_000915.1_without_BAQ_NC_000915.1_SNP_QUAL.png  stat_NC_000915.1_with_BAQ_SNP.csv
    NC_000915.1_with_BAQ_NC_000915.1_SNP_QUAL.png    stat_NC_000915.1_extend_BAQ_SNP.csv               stat_NC_000915.1_without_BAQ_SNP.csv

Mapping Gene ontology
------------------

Gene ontology is useful for understanding the function of gene product. ``go_term`` is the 
subcommand for mapping your annotation to gene ontology. Before running ``go_term``, we 
need to prepare some database. First, please download 
`goslim.obo <http://geneontology.org/page/go-slim-and-subset-guide>`_ and 
`go.obo <http://geneontology.org/page/download-ontology>`_ and 
`idmapping_selected.tab <http://www.uniprot.org/downloads>`_.

::

    $ wget -cP ANNOgesic/input/database http://www.geneontology.org/ontology/subsets/goslim_generic.obo
    $ wget -cP ANNOgesic/input/database http://geneontology.org/ontology/go.obo
    $ wget -cP ANNOgesic/input/database ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
    $ gunzip ANNOgesic/input/database/idmapping_selected.tab.gz

Now, we have all required database. Let's try it.

::

    python3 annogesic go_term \
        -g ANNOgesic/output/target/annotation\
        ANNOgesic

The output of ``go_term`` will store in ``Go_term_results``. The statistics files and 
visualization files will store in ``statistics``.

::

    $ ls ANNOgesic/output/Go_term/Go_term_results/NC_000915/
    all_strains_uniprot.csv
    $ ls ANNOgesic/output/Go_term/statistics/NC_000915/
    figs  stat_NC_000915.csv
    $ ls ANNOgesic/output/Go_term/statistics/NC_000915/figs/
    NC_000915.1_biological_process.png  NC_000915.1_cellular_component.png  NC_000915.1_molecular_function.png  NC_000915.1_three_roots.png

Prediction of Subcellular localization
------------------

Subcellular localization may be a useful information for analysis of protein function. For 
generating the information of subcellular localization, you can use the subcommand 
``subcellular_localization`` to get it.

::

    python3 annogesic subcellular_localization \
        -g ANNOgesic/output/target/annotation \
        -f ANNOgesic/output/target/fasta \
        -m -b positive \
        ANNOgesic

The output of ``subcellular_localization`` will generate two folder. ``psortb_results`` will 
store the output of `Psortb <http://www.psort.org/psortb/>`_. ``statistics`` will store 
the statistics files and visualization files.

::

    $ ls ANNOgesic/output/subcellular_localization/
    psortb_results  statistics
    $ ls ANNOgesic/output/subcellular_localization/psortb_results/NC_000915/
    NC_000915.1_raw.txt  NC_000915_table.csv
    $ ls ANNOgesic/output/subcellular_localization/statistics/NC_000915/
    NC_000915_NC_000915.1.png  stat_NC_000915_sublocal.csv

Generating protein-protein interaction network
-------------------

Protein-protein interaction network is an important feature for analysis of regulation. 
The subcommand ``ppi_network`` combines STRING(database of protein-protein interaction) 
and PIE(text-mining for protein-protein interaction). It can generate the protein-protein 
interaction network with supported literatures. You can refer to relevance of literatures 
and network to find your interesting candidates.

Before running the subcommand, you need to download the 
`species.vXXXX.txt from STRING <http://string-db.org/newstring_cgi/show_download_page.pl?UserId=ReWbu8uLrfAN&sessionId=_FAQBbatf7RX>`_

::

    wget -cP ANNOgesic/input/database http://string-db.org/newstring_download/species.v10.txt

Now, we can try the subcommand.

::

    python3 annogesic ppi_network \
        -s NC_000915.ptt:'Helicobacter pylori 26695 chromosome':'Helicobacter pylori 26695':'Helicobacter pylori' \
        -p $ANNOGESIC_FOLDER/output/target/annotation \
        -d $ANNOGESIC_FOLDER/input/database/species.v10.txt \
        -q 'Helicobacter pylori 26695 chromosome':217:633:- 'Helicobacter pylori 26695 chromosome':2719:3402:+ \
        -n \
        ANNOgesic

If you want to get all proteins in ptt files, you just need to assign ``all`` in ``-q``.

``ppi_network`` will generate three folders.

::

    $ ls ANNOgesic/output/PPI/
    all_results/  best_results/ figures/

``all_results`` is for all interactions without filtering. ``best_results`` is for the interactions with 
filtering of `PIE <http://www.ncbi.nlm.nih.gov/CBBresearch/Wilbur/IRET/PIE/>`_ score. ``figures`` is for 
the figures of PPI network. There are two subfolders - ``with_strain`` and ``without_strain`` in these folders. 
These two folders stores all information of interactions and literature scores. ``with_strain`` is for 
the information with giving specific strain name for searching literatures. ``without_strain`` is for the 
information without giving specific strain name for searching literatures.

::

    $ ls ANNOgesic/output/PPI/all_results/NC_000915/
    NC_000915_without_strain.csv  NC_000915_with_strain.csv     without_strain/               with_strain/
    $ ls ANNOgesic/output/PPI/best_results/NC_000915/
    NC_000915_without_strain.csv  NC_000915_with_strain.csv  without_strain  with_strain
    $ ls ANNOgesic/output/PPI/figures/NC_000915/
    without_strain  with_strain
    $ ls ANNOgesic/output/PPI/all_results/NC_000915/with_strain/Helicobacter\ pylori\ 26695\ chromosome/
    carA_C694_01345.csv  carA_guaB.csv  carB_C694_01345.csv  guaB_C694_01345.csv  pyrB_carB.csv  pyrD_C694_01345.csv  ribD_ribH.csv
    carA_carB.csv        carA_pyrB.csv  carB_guaB.csv        nusG_rpoB.csv        pyrB_guaB.csv  pyrE_pyrF.csv        rpsJ_nusB.csv
    $ ls ANNOgesic/output/PPI/all_results/NC_000915/without_strain/Helicobacter\ pylori\ 26695\ chromosome/
    carA_C694_01345.csv  carA_guaB.csv  carB_C694_01345.csv  guaB_C694_01345.csv  pyrB_carB.csv  pyrD_C694_01345.csv  ribD_ribH.csv
    carA_carB.csv        carA_pyrB.csv  carB_guaB.csv        nusG_rpoB.csv        pyrB_guaB.csv  pyrE_pyrF.csv        rpsJ_nusB.csv
    $ ls ANNOgesic/output/PPI/best_results/NC_000915/without_strain/Helicobacter\ pylori\ 26695\ chromosome/
    carA_C694_01345.csv  carA_carB.csv  carA_pyrB.csv  carB_C694_01345.csv  guaB_C694_01345.csv  pyrD_C694_01345.csv
    $ ls ANNOgesic/output/PPI/best_results/NC_000915/with_strain/Helicobacter\ pylori\ 26695\ chromosome/
      (It can't find any interactions)
    $ ls ANNOgesic/output/PPI/figures/NC_000915/with_strain/
    Helicobacter pylori 26695 chromosome
    $ ls ANNOgesic/output/PPI/figures/NC_000915/with_strain/Helicobacter\ pylori\ 26695\ chromosome/
      (It can't find any interactions. Therefore, it has no figures)
    $ ls ANNOgesic/output/PPI/figures/NC_000915/without_strain/Helicobacter\ pylori\ 26695\ chromosome/
    HP0001_nusB.png  HP0005_pyrF.png 

Generating riboswitch
-----------------

If you want to know the candidates of riboswitch, you can use the subcommand ``riboswitch``.
Before running ``riboswitch``, you need to get the information of all riboswitch in Rfam. 
You can use ours or create your own one. For testing, you can just go 
`here <https://github.com/Sung-Huan/ANNOgesic>`_ and downloan or copy ``Rfam_riboswitch_ID.csv``.
Then you can put ``Rfam_riboswitch_ID.csv`` in ``ANNOgesic/input/riboswitch_ID``

::

    $ cp /path/of/Rfam_riboswitch_ID.csv ANNOgesic/input/riboswitch_ID

You also need to download Rfam.

::

    $ wget ANNOgesic/input/database ftp://ftp.ebi.ac.uk/pub/databases/Rfam/12.0/Rfam.tar.gz
    $ cd ANNOgesic/input/database
    $ tar -zxvf Rfam.tar.gz
    $ rm Rfam.tar.gz
    $ cd ../../../

Now we can try the subcommand.

::

    python3 annogesic riboswitch \
        -g ANNOgesic/output/target/annotation \
        -f ANNOgesic/output/target/fasta \
        -r \
        -i ANNOgesic/input/riboswitch_ID/Rfam_riboswitch_ID.csv \
        -R ANNOgesic/input/database/CMs/Rfam.cm \
        ANNOgesic

The output is the following. ``gffs`` is for gff files; ``tables`` is for tables of riboswitch; 
``scan_Rfam`` is for the output of scanning Rfam; ``statistics`` is for the statistics files.

::

     $ ls ANNOgesic/output/riboswitch/
     gffs  scan_Rfam  statistics  tables
     $ ls ANNOgesic/output/riboswitch/gffs/
     NC_000915_RBS.gff
     $ ls ANNOgesic/output/riboswitch/scan_Rfam/NC_000915/
     NC_000915.1_RBS_rescan.txt  NC_000915.1_RBS.txt
     $ ls ANNOgesic/output/riboswitch/tables/
     NC_000915_RBS.csv
     $ ls ANNOgesic/output/riboswitch/statistics/
     stat_NC_000915_RBS.txt


Producing the screenshots
-----------------

It is a good idea if we can get the screenshots of our interesting features. Then we can 
check them very quickly. Therefore, ANNOgesic provide a subcommand ``screenshot`` for 
generating screenshots.

Before you running it, you have to install `IGV <https://www.broadinstitute.org/software/igv/home>`_ 
for generating screenshot.

For testing, we use TSS as main feature, sRNA and CDS information as side features.

::

    python3 annogesic screenshot \
        -mg ANNOgesic/output/TSS/gffs/NC_000915_TSS.gff \
        -sg ANNOgesic/output/target/annotation/NC_000915.gff \
            ANNOgesic/output/sRNA/gffs/best/NC_000915_sRNA.gff \
        -f ANNOgesic/output/target/fasta/NC_000915.fa \
        -o ANNOgesic/output/TSS/screenshots \
        -tl $tex_notex_libs \
        -tw ANNOgesic/input/wigs/tex_notex \
    ANNOgesic

``screenshot`` will generate two txt files and two folders.

::

    $ ls ANNOgesic/output/TSS/screenshots/NC_000915/
    forward/     forward.txt  reverse/     reverse.txt

``forward.txt`` and ``reverse.txt`` are batch files for `IGV <https://www.broadinstitute.org/software/igv/home>`_.
``forward`` and ``reverse`` are the folders for storing screenshots.

Now, please open `IGV <https://www.broadinstitute.org/software/igv/home>`_. Please follow the command: Tools -> 
Run Batch Script -> choose ``forward.txt``. When it has done, please do it again for reverse strand: Tools ->
Run Batch Script -> choose ``reverse.txt``. If you just want to test it and don't want to wait a long time for 
generating screenshots, you can delete some lines of gff files of TSS.

After that, you can see that there are several screenshots in ``forward`` and ``reverse``.

::

    $ ls ANNOgesic/output/TSS/screenshots/NC_000915/forward
    NC_000915.1:1002230-1002230.png  NC_000915.1:1245705-1245705.png  NC_000915.1:1516949-1516949.png  NC_000915.1:246369-246369.png  NC_000915.1:463673-463673.png  NC_000915.1:741002-741002.png
    NC_000915.1:1002524-1002524.png  NC_000915.1:124623-124623.png    NC_000915.1:151822-151822.png    NC_000915.1:251753-251753.png  NC_000915.1:463731-463731.png  NC_000915.1:744418-744418.png
    NC_000915.1:1002728-1002728.png  NC_000915.1:1249488-1249488.png  NC_000915.1:1520156-1520156.png  NC_000915.1:255496-255496.png  NC_000915.1:464179-464179.png  NC_000915.1:744551-744551.png
    ...
    
    $ ls ANNOgesic/output/TSS/screenshots/NC_000915/reverse
    NC_000915.1:1002215-1002215.png  NC_000915.1:1235881-1235881.png  NC_000915.1:1481470-1481470.png  NC_000915.1:179609-179609.png  NC_000915.1:467716-467716.png  NC_000915.1:767765-767765.png
    NC_000915.1:1002707-1002707.png  NC_000915.1:1238472-1238472.png  NC_000915.1:1482537-1482537.png  NC_000915.1:181416-181416.png  NC_000915.1:46780-46780.png    NC_000915.1:769891-769891.png
    NC_000915.1:100498-100498.png    NC_000915.1:1240517-1240517.png  NC_000915.1:1482926-1482926.png  NC_000915.1:181781-181781.png  NC_000915.1:468289-468289.png  NC_000915.1:770316-770316.png
    ...


Coloring the screenshots
-----------------

If your tracks are many and you want to check TSS, it will be painful for distinguish the 
tracks of TEX+ and TEX-. Therefore, we provide a subcommand ``color_png`` for coloring 
your screenshots.

::

    python3 annogesic color_png \
        -t 2 \
        -f ANNOgesic/output/TSS \
        ANNOgesic

You will see the png files are not different as before. However, when you open them, the tracks are colored.

::

    $ ls ANNOgesic/output/TSS/screenshots/NC_000915/forward
    NC_000915.1:1002230-1002230.png  NC_000915.1:1245705-1245705.png  NC_000915.1:1516949-1516949.png  NC_000915.1:246369-246369.png  NC_000915.1:463673-463673.png  NC_000915.1:741002-741002.png
    NC_000915.1:1002524-1002524.png  NC_000915.1:124623-124623.png    NC_000915.1:151822-151822.png    NC_000915.1:251753-251753.png  NC_000915.1:463731-463731.png  NC_000915.1:744418-744418.png
    NC_000915.1:1002728-1002728.png  NC_000915.1:1249488-1249488.png  NC_000915.1:1520156-1520156.png  NC_000915.1:255496-255496.png  NC_000915.1:464179-464179.png  NC_000915.1:744551-744551.png
    ...
    
    $ ls ANNOgesic/output/TSS/screenshots/NC_000915/reverse
    NC_000915.1:1002215-1002215.png  NC_000915.1:1235881-1235881.png  NC_000915.1:1481470-1481470.png  NC_000915.1:179609-179609.png  NC_000915.1:467716-467716.png  NC_000915.1:767765-767765.png
    NC_000915.1:1002707-1002707.png  NC_000915.1:1238472-1238472.png  NC_000915.1:1482537-1482537.png  NC_000915.1:181416-181416.png  NC_000915.1:46780-46780.png    NC_000915.1:769891-769891.png
    NC_000915.1:100498-100498.png    NC_000915.1:1240517-1240517.png  NC_000915.1:1482926-1482926.png  NC_000915.1:181781-181781.png  NC_000915.1:468289-468289.png  NC_000915.1:770316-770316.png
    ...


Now you already finished your first wonderful trip of annotation. Hopefully, you enjoy it!!
