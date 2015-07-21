ANNOgesic - Transcriptome annotation pipeline
*****************************************
Table of content
================

.. toctree::
   :maxdepth: 1

   index
   prerequired
   subcommands
   docker
   troubleshooting

Introduction
=========================
ANNOgesic is a bacterial transcriptome annotation pipeline based on RNA-Seq. 
ANNOgesic covers different aspects of the genome annotation. In order to get the
best results, ANNOgesic also can optimize the parameters of some tools. 
ANNOgesic can automatically generate high-quality annotation information for
query strains. Moreover, it is modular and its subcommands can be separately
used. Currently, the pipeline already can detect or integrate a) fasta, CDS,
tRNA, rRNA and genes of query genome, b) transcription starting sites
(TSSs), c) rho-independent terminators, d) transcript assembly,
e) untranslated region (UTRs), f) sRNA, g) promoters, h) processing sites,
i) circular RNAs, j) protein-protein interaction networks, k) potential
sRNA target, l) single-nucleotide polymorphism (SNP), m) operons,
n) GO terms, o) subcellular localization, p) riboswitch, q) potential sORF.

::

   usage: ANNOgesic.py [-h] [--version]
                    {create,get_input_files,get_target_fasta,annotation_transfer,expression_analysis,tsspredator,optimize_tsspredator,color_png,terminator,transcript_assembly,utr,srna,sorf,promoter,operon,circrna,go_term,srna_target,snp,ppi_network,subcellular_localization,riboswitch,screenshot}
                    ...

   positional arguments:
   {create,get_input_files,get_target_fasta,annotation_transfer,expression_analysis,tsspredator,optimize_tsspredator,color_png,terminator,transcript_assembly,utr,srna,sorf,promoter,operon,circrna,go_term,srna_target,snp,ppi_network,subcellular_localization,riboswitch,screenshot}
                        commands
    create              Create a project
    get_input_files     Get required files (i.e. annotation files, fasta
                        files)
    get_target_fasta    Get target fasta.
    annotation_transfer
                        Run RATT to transfer the annotation files from
                        reference to target.
    expression_analysis
                        Run gene expression analysis to compare which CDS is
                        expressing in which libraries
    tsspredator         Run TSSpredator to predict TSSs or processing sites.
    optimize_tsspredator
                        Optimize TSSpredator based on (partial)manual detect
                        one.
    color_png           Generating color screenshots of TSS or processing
                        site. It only works after running batch script.
    terminator          Detect Terminators.
    transcript_assembly
                        Run Transcript for doing transcriptome assembly.
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
                        prediction of subcellular localization of genomic CDS.
    riboswitch          prediction of riboswitch.
    screenshot          Generate screenshot for selected feature.

   optional arguments:
   -h, --help            show this help message and exit
   --version, -v         show version

Download
========

Source code
===========

The source code of ANNOgesic can be found at `<https://github.com/Sung-Huan/ANNOgesic/tree/master>`_.

Cite
====

Contact
=======

For question and requests feel free to contact `Sung-Huan Yu
<sung-huan.yu@uni-wuerzburg.de>`_
