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
ANNOgesic is a bacterial transcriptome annotation pipeline based on RNA-Seq. 
ANNOgesic covers different aspects of the genome annotation. In order to get the
best results, ANNOgesic also can optimize the parameters of some tools. 
ANNOgesic can automatically generate high-quality annotation information for
query strains. Moreover, it is modular and its subcommands can be separately used.
ANNOgesic integrates six main classes of annotations. i) Reference
genome improvement: SNP/mutation calling, Sequence modification and
annotation transfer. ii) Transcript boundary: TSS, transcript,
terminator, UTR and processing site.  iii) sRNA and sORF: sRNA, sORF
and sRNA target prediction.  iv) Functional related features:
protein-protein interaction networks, Gene ontology and subcellular
localization. v) Promoter and operon: promoter motifs and operon
with sub-operon. vi) Other features: circular RNA, CRISPR and riboswitch.

::

     usage: annogesic [-h] [--version]
                      {create,get_input_files,get_target_fasta,annotation_transfer,tsspredator,optimize_tsspredator,color_png,terminator,transcript_assembly,utr,srna,sorf,promoter,operon,circrna,go_term,srna_target,snp,ppi_network,subcellular_localization,riboswitch,crispr,screenshot}
                      ...
     
     positional arguments:
       {create,get_input_files,get_target_fasta,annotation_transfer,tsspredator,optimize_tsspredator,color_png,terminator,transcript_assembly,utr,srna,sorf,promoter,operon,circrna,go_term,srna_target,snp,ppi_network,subcellular_localization,riboswitch,crispr,screenshot}
                             commands
         create              Create a project
         get_input_files     Get required files. (i.e. annotation files, fasta
                             files)
         get_target_fasta    Get target fasta.
         annotation_transfer
                             Run RATT to transfer the annotation files from
                             reference to target.
         tsspredator         Run TSSpredator to predict TSSs or processing sites.
         optimize_tsspredator
                             Optimize TSSpredator based on (partial)manual detect
                             one.
         color_png           Generating color screenshots of TSS or processing
                             site. It only works after running batch script.
         terminator          Detect Terminators.
         transcript_assembly
                             Run transcriptome assembly for detecting transcripts.
         utr                 Run UTR detection to detect 5'UTR and 3'UTR.
         srna                Run sRNA detection to detect sRNA candidates.
         sorf                Run sORF detection to detect sORF candidates which has
                             expression.
         promoter            Run MEME to dicover promoter.
         operon              Detect operon and combine features together.
         circrna             Detect circular RNA.
         go_term             Extract and find Go terms.
         srna_target         sRNA target prediction.
         snp                 Detection of SNP of transcripts.
         ppi_network         Generate protein-protein interaction with literature
                             supported.
         subcellular_localization
                             Prediction of subcellular localization of genomic CDS.
         riboswitch          Prediction of riboswitch.
         crispr              Prediction of CRISPR.
         screenshot          Generate screenshot for selected feature.
     
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
