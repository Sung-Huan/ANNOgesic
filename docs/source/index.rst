ANNOgesic - Transcriptome annotation pipeline
*****************************************
Table of content
================

.. toctree::
   :maxdepth: 1

   index
   prerequired
   installation
   docker
   subcommands
   test_case
   troubleshooting
   license

Introduction
=========================
ANNOgesic is a modular, command-line tool that can
integrated different types of RNA-Seq data like dRNA-Seq or RNA-Seq
generated after transcript fragmentation and generates high quality
genome annotations. It can detect gene, CDS/tRNA/rRNA, TSS and
processing sites, transcripts, terminator, Untranslated region (UTR)
as well as small RNA (sRNA), small open reading frame (sORF), circular
RNA, CRISPR related RNAs, riboswitch and RNA-thermometer. 
It can also perform RNA-RNA
and protein-protein interaction predictions. Furthermore, it groups
genes into operon and sub-operons and reveal promotor motifs. It can
also allocate GO term and subcellular localization to genes. Several
of ANNOgesic features are new implementation while others are
performed and improved by third-party tools and for some of them
adaptive parameter-optimizations were included. Additionally, numerous
visualization and statistitcs help the user quickly evaluated feature
predictions resulting from an ANNOgesic analysis. The pipeline is
modular and was heavily tested with several RNA-Seq data set from
bacterial as well as archaeal samples.

::

    usage: annogesic [-h] [--version]
                     {create,get_input_files,get_target_fasta,annotation_transfer,tss_processing,optimize_tss_processing,terminator,transcript,utr,srna,sorf,promoter,operon,circrna,go_term,srna_target,snp,ppi_network,subcellular_localization,riboswitch_thermometer,crispr,merge_features,screenshot,color_png}
                     ...
    
    positional arguments:
      {create,get_input_files,get_target_fasta,annotation_transfer,tss_processing,optimize_tss_processing,terminator,transcript,utr,srna,sorf,promoter,operon,circrna,go_term,srna_target,snp,ppi_network,subcellular_localization,riboswitch_thermometer,crispr,merge_features,screenshot,color_png}
                            commands
        create              Create a project
        get_input_files     Get required files. (i.e. annotation files, fasta
                            files)
        get_target_fasta    Get target fasta.
        annotation_transfer
                            Run RATT to transfer the annotation files from
                            reference to target.
        tss_processing      Run TSSpredator to predict TSSs or processing sites.
        optimize_tss_processing
                            Optimize TSSpredator based on (partial)manual detect
                            one.
        terminator          Detect rho-independant terminators.
        transcript          Detect transripts based on coverage file.
        utr                 Run UTR detection to detect 5'UTR and 3'UTR.
        srna                Detect intergenic, antisense and UTR-derived sRNA.
        sorf                Detect expressed sORF.
        promoter            Run MEME to dicover promoter.
        operon              Detect operon and sub-operon.
        circrna             Detect circular RNA by segemehl.
        go_term             Extract Go terms from Uniprot.
        srna_target         Detect sRNA-mRNA interaction by RNAup and RNAplex.
        snp                 Detect SNP/mutation and generate potential fasta file.
        ppi_network         Generate protein-protein interaction with literature
                            supported.
        subcellular_localization
                            Predict subcellular localization of genome CDS.
        riboswitch_thermometer
                            Predict riboswitch and RNA thermometer.
        crispr              Run CRT to predict CRISPR.
        merge_features      Merge all features to one gff file.
        screenshot          Generate screenshot for selected feature.
        color_png           Generate color screenshots of TSS or processing site.
                            It only works after running "screenshot" (after
                            running batch script).
    
    optional arguments:
      -h, --help            show this help message and exit
      --version, -v         show version

Download
========

::

    git clone git@github.com:Sung-Huan/ANNOgesic.git

or


::

    pip3 install annogesic

Source code
===========

The source code of ANNOgesic can be found at `Github <https://github.com/Sung-Huan/ANNOgesic>`_.

Cite
====

Contact
=======

For question and requests feel free to contact `Sung-Huan Yu
<sung-huan.yu@uni-wuerzburg.de>`_
