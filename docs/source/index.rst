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
integrate different types of RNA-Seq data like dRNA-Seq or RNA-Seq
generated after transcript fragmentation and generates high quality
genome annotations. It can detect genes, CDSs/tRNAs/rRNAs, transcription
starting sites (TSS) and processing sites, transcripts, terminators, 
untranslated regions (UTR) as well as small RNAs (sRNA), small open 
reading frames (sORF), circular RNAs, CRISPR related RNAs, riboswitches 
and RNA-thermometers. It can also perform RNA-RNA
and protein-protein interactions prediction. Furthermore, it groups
genes into operons and sub-operons and reveal promoter motifs. It can
also allocate GO term and subcellular localization to genes. Several
of ANNOgesic features are new implementation while others are
performed and improved by third-party tools and for some of them
adaptive parameter-optimizations were included. Additionally, numerous
visualization and statistics help the user quickly evaluated feature
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
                            Transfer the annotations from reference genome to
                            target genome.
        tss_processing      Detect TSSs or processing sites.
        optimize_tss_processing
                            Optimize TSSs or processing sites based on manual
                            detected ones.
        terminator          Detect rho-independant terminators.
        transcript          Detect transripts based on coverage file.
        utr                 Detect 5'UTRs and 3'UTRs.
        srna                Detect intergenic, antisense and UTR-derived sRNAs.
        sorf                Detect expressed sORFs.
        promoter            Discover promoter motifs.
        operon              Detect operons and sub-operons.
        circrna             Detect circular RNAs.
        go_term             Extract Go terms from Uniprot.
        srna_target         Detect sRNA-mRNA interactions.
        snp                 Detect SNP/mutation and generate potential fasta file.
        ppi_network         Detect protein-protein interactions with literature
                            supported.
        subcellular_localization
                            Predict subcellular localization of CDSs.
        riboswitch_thermometer
                            Predict riboswitches and RNA thermometers.
        crispr              Predict CRISPR related RNAs.
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

License
=======

`ISC <https://en.wikipedia.org/wiki/ISC_license>`__ (Internet Systems
Consortium license ~ simplified BSD license) - see `LICENSE <https://pythonhosted.org/ANNOgesic/license.html>`__

Cite
====

Contact
=======

For question and requests feel free to contact `Sung-Huan Yu
<sung-huan.yu@uni-wuerzburg.de>`_
