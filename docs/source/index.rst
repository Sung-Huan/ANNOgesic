ANNOgesic - THE tool for bacterial/archaeal RNA-Seq based genome annotations
********************************************************************************

Table of content
================

.. toctree::
   :maxdepth: 1

   index
   required
   installation
   docker
   subcommands
   tutorial
   questions
   license

Introduction
=========================

ANNOgesic is the swiss army knife for RNA-Seq based annotation of
bacterial/archaeal genomes.

It is a modular, command-line tool that can integrate different types
of RNA-Seq data based on dRNA-Seq (differential RNA-Seq) or RNA-Seq
protocols that inclusde transcript fragmentation to generate high
quality genome annotations. It can detect genes, CDSs/tRNAs/rRNAs,
transcription starting sites (TSS) and processing sites, transcripts,
terminators, untranslated regions (UTR) as well as small RNAs (sRNA),
small open reading frames (sORF), circular RNAs, CRISPR related RNAs,
riboswitches and RNA-thermometers. It can also perform RNA-RNA and
protein-protein interactions prediction. Furthermore, it groups genes
into operons and sub-operons and reveal promoter motifs. It can also
allocate GO term and subcellular localization to genes. Several of
ANNOgesic features are new implementations while other build on well
known third-party tools for which it offers adaptive
parameter-optimizations. Additionally, numerous visualization and
statistics help the user to quickly evaluat feature predictions
resulting from an ANNOgesic analysis. The tool was heavily tested
with several RNA-Seq data set from bacterial as well as archaeal
samples.

::

    usage: annogesic [-h] [--version]
                     {create,get_input_files,update_genome_fasta,annotation_transfer,
                      tss_ps,optimize_tss_ps,terminator,transcript,utr,srna,sorf,
                      promoter,operon,circrna,go_term,srna_target,snp,ppi_network,
                      localization,riboswitch_thermometer,crispr,merge_features,
                      screenshot,colorize_screenshot_tracks}
                     ...
    
    positional arguments:
      {create,get_input_files,update_genome_fasta,annotation_transfer,tss_ps,optimize_tss_ps,
       terminator,transcript,utr,srna,sorf,promoter,operon,circrna,go_term,srna_target,snp,
       ppi_network,localization,riboswitch_thermometer,crispr,merge_features,screenshot,
       colorize_screenshot_tracks}
                            commands
        create              Create a project
        get_input_files     Get required files. (i.e. annotation files, fasta
                            files)
        update_genome_fasta
                            Get fasta files of query genomes if the query
                            sequences do not exist.
        annotation_transfer
                            Transfer the annotations from a closely related
                            species genome to a target genome.
        tss_ps              Detect TSSs or processing sites.
        optimize_tss_ps     Optimize TSSs or processing sites based on manual
                            detected ones.
        terminator          Detect rho-independent terminators.
        transcript          Detect transcripts based on coverage file.
        utr                 Detect 5'UTRs and 3'UTRs.
        srna                Detect intergenic, antisense and UTR-derived sRNAs.
        sorf                Detect expressed sORFs.
        promoter            Discover promoter motifs.
        operon              Detect operons and sub-operons.
        circrna             Detect circular RNAs.
        go_term             Extract GO terms from Uniprot.
        srna_target         Detect sRNA-mRNA interactions.
        snp                 Detect SNP/mutation and generate fasta file if
                            mutations were found.
        ppi_network         Detect protein-protein interactions suported by
                            literature.
        localization        Predict subcellular localization of proteins.
        riboswitch_thermometer
                            Predict riboswitches and RNA thermometers.
        crispr              Predict CRISPR related RNAs.
        merge_features      Merge all features to one gff file.
        screenshot          Generate screenshots for selected features using IGV.
        colorize_screenshot_tracks
                            Add color information to screenshots (e.g. useful for
                            dRNA-Seq based TSS and PS detection. It only works
                            after running "screenshot" (after running batch
                            script).
    
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

Please check :ref:`Required tools or databases` for the required third-party tools.

Source code
===========

The source code of ANNOgesic can be found at `Github <https://github.com/Sung-Huan/ANNOgesic>`_.

License
=======

`ISC <https://en.wikipedia.org/wiki/ISC_license>`__ (Internet Systems
Consortium license ~ simplified BSD license) - see `LICENSE <http://annogesic.readthedocs.io/en/latest/license.html>`__

Cite
====

Contact
=======

For question and requests feel free to contact `Sung-Huan Yu
<sung-huan.yu@uni-wuerzburg.de>`_
