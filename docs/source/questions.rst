.. _FAQ:

FAQ
===

No matter which kinds of error messages occur, please always make sure **the seq IDs or strain names of 
all Gff files, fasta files, BAM files and wiggle files must be consistent.**

- **Q1:**

I got the error about complaining no module can be found as following:

::

    Traceback (most recent call last):
      File "ANNOgesic.py", line 6, in <module>
        from annogesiclib.controller import Controller
    ImportError: No module named annogesiclib.controller

- **A1:**

Switch to Python 3.3 or higher. The newer versions can handle this.

- **Q2:**

I try to install ANNOgesic in my own directory. However, when I execute ANNOgesic, an error comes out as following:

::

    Traceback (most recent call last):
      File "bin/annogesic", line 7, in <module>
        from annogesiclib.controller import Controller
    ImportError: No module named 'annogesiclib'

- **A2:**

Please generate a soft link of ``annogesiclib`` to ``bin``.

::

    $ cd ANNOgesic/bin
    $ ln -s ../annogesiclib .

- **Q3:**

When I finished sRNA detection, there are some files in ``sRNAs/gffs/for_classes`` and ``sRNAs/tables/for_classes``. 
What do the classes mean?

- **A3:**

Please check ``stat_sRNA_class_$GENOME.csv`` in ``stat`` folder. The information of classes can be found in it. The 
classes are based on the input information.

- **Q4:**

I got a sRNA gff file which only contain ``##gff-version 3``. What does it mean?

- **A4:**

If you only got ``##gff-version 3`` in sRNA gff file, it means that ANNOgesic can not find any sRNA in this genome. 
This case can be applied to other genomic feature prediction.

- **Q5:**

How long does sRNA target prediction take?

- **A5**

It depends on the input data size. However, most of the sRNA target prediction tools need to spend around one hour for detecting 
targets of one sRNA. RNAplex is a fast one which can detect the targets of one sRNA in several minutes. If you want to reduce the 
running time, you can specify ``--program RNAplex``. However, the accuracy may be influenced.

- **Q6:**

When I run sRNA target prediction, I found many files in ``RNAplfold``. Are they necessary?

- **A6:**

These files are intermediate input for RNAplex. They will be deleted when RNAplex is done.

- **Q7:**

If I only have fragmented RNA-Seq libraries or conventional RNA-Seq libraries, what kinds of annotations can be generated? On the 
other hands, if only dRNA-Seq libraries are available, what kinds of annotations can be generated?

- **A7:**

Only fragmented RNA-Seq or conventional RNA-Seq libraries are available:

dRNA-Seq data (TEX+/-) are mainly used for TSS detection. Therefore, if the user has no dRNA-Seq data, the annotations 
which need TSS information may not be able to detect. The annotations which can be detected without dRNA-Seq libraries are: 
transcript, operon (sub-operon can not be detected), terminator, CRISPR, sRNA (the results may be influenced, 
and UTR-derived sRNA can not be detected), sRNA target, GO term, subcellular localization, circular RNA, SNP, 
protein-protein interaction, annotation transfer. For the details of sRNA detection, if TSSs can not be obtained, 
the selection of the best sRNA candidates can only considers the information like secondary structure folding energy, 
BLAST results, etc. Thus, the accuracy of sRNA detection (especially for the selection of the best sRNA candidates) 
will be influenced. For operon detection, since the TSS information is not available, the sub-operons can not be detected.

The other features (including UTR-derived sRNAs) which can not be detected precisely or at all by only fragmented or 
conventional RNA-Seq data require the information of TSS or ribosome-binding site.

Only dRNA-Seq libraries are available:

Actually, all the genomic features can be detected by only using dRNA-Seq data. However, dRNA-Seq usually does not able 
to detect transcript boundary. Thus, the genomic features which are related to transcript boundary will be influenced, such as 
transcripts, sRNAs, sORFs, operon, terminators, etc. Although the 3'end of transcripts may be not clear, the genomic feature detections 
of ANNOgesic still show high performances by comparing to the previously published data.

- **Q8:**

An error message related to 'defined(@array)' shows when I run annotation_transfer.

- **A8:**

This issue is caused by different version of Perl. RATT is written by the old version of Perl. 
Please do the following steps to fix this error.

1. Download RATT via `PAGIT <http://www.sanger.ac.uk/resources/software/pagit/>`_
2. Set environment path.
3. Run sed -i '244s/defined//' $RATT_FOLDER/main.ratt.pl && sed -i '19s/$PAGIT_HOME/\/usr/' $RATT_FOLDER/start.ratt.sh (Please specify the folder of your RATT to $RATT_FOLDER)

We also suggest the user running ANNOgesic via our `Docker image <https://hub.docker.com/r/silasysh/annogesic/>`_. This issue was fixed in our Docker image.

.. _FAQ Q9:

- **Q9:**

After running sRNA prediction of ANNOgesic, I got a huge number of sRNA candidates. The amount is too many for me, how can I remove some low-confidence ones?

- **A9:**

Although the sRNA prediction of ANNOgesic was tested by RNA-Seq data from different strains of bacteria, 
the setting may not be able to apply to some other specific RNA-Seq data. 
If the amount of candidates are too many, it may be due to the low coverage rate of RNA-Seq data.
Therefore, many short transcripts are genearted. There are several ways to further filted out the
sRNAs.

- Using promoter and terminator for ``--filter_info``: We suggest the user to put promoter information first.
If terminator information is still needed to add to ``--filter_info``, please run ``terminator`` with sRNA files
(``--srna_files``) first. Then using that terminator output file to run ``srna``.

- Using ranking number: The output table of sRNAs contains the ranking number based on coverage and promoter information.
The user can select top ranking ones to do the following analyses.

- Adjusting coverage cutoff (not recommended): There are some coverage cutoffs can be modified by the user. The
default setting is tested based on some RNA-Seq data from different bacteria. Thus, changing these parameters may
influence the accuracy of the prediction.

