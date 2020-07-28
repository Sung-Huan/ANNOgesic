.. image:: https://img.shields.io/pypi/v/annogesic.svg
   :target: https://pypi.python.org/pypi/ANNOgesic/
.. image:: https://img.shields.io/pypi/l/annogesic.svg
   :target: https://pypi.python.org/pypi/ANNOgesic/
.. image:: https://zenodo.org/badge/34061246.svg
   :target: https://zenodo.org/badge/latestdoi/34061246

About ANNOgesic
---------------

test

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

Documentation
-------------

Documentation can be found on
`here <http://annogesic.readthedocs.io/en/latest/index.html>`__.

Installation
------------

Pip3
^^^^

If you have all the requirements installed, installation can be done by 
using pip3.

With root permission: 

::

    $ pip3 install ANNOgesic

without root permission

::

    $ pip3 install --user ANNOgesic

If you want to know the requirement, please refer to 
`Documentation <http://annogesic.readthedocs.io/en/latest/index.html>`__.

Docker and Singularity
^^^^^^^^^^^^^^^^^^^^^^

In order to solve the issue of installing the dependencies, a Docker image of ANNOgesic is provided in 
`Docker Hub <https://hub.docker.com/r/silasysh/annogesic/>`__. 
The image can be pulled by using `Docker <https://www.docker.com/>`__ or `Singularity <https://singularity.lbl.gov/index.html>`__. 
Moreover, the users can build an Docker image from the `Dockerfile <https://github.com/Sung-Huan/ANNOgesic>`__ by themselves. 
For the details, please check the `documentation <http://annogesic.readthedocs.io/en/latest/installation.html>`__.

Github
^^^^^^

The alternative way for installing ANNOgesic is directly clone the Git repository.

::

    $ git clone https://github.com/Sung-Huan/ANNOgesic.git

or

::

    $ git clone git@github.com:Sung-Huan/ANNOgesic.git

In order to make ANNOgesic runnable, we should create a soft link of ``annogesiclib`` in ``bin``.

::

    $ cd ANNOgesic/bin
    $ ln -s ../annogesiclib .

Then, you can run ANNOgesic via specifying the installed path if all the requirements are setup properly.

Arguments
-------------

::

    usage: annogesic [-h] [--version]
                     {create,get_input_files,update_genome_fasta,annotation_transfer,
                      tss_ps,optimize_tss_ps,terminator,transcript,utr,srna,sorf,
                      promoter,operon,circrna,go_term,srna_target,snp,ppi_network,
                      localization,riboswitch_thermometer,crispr,merge_features,
                      screenshot,colorize_screenshot_tracks}
                     ...
    
    positional arguments:
      {create,get_input_files,update_genome_fasta,annotation_transfer,tss_ps,
       optimize_tss_ps,terminator,transcript,utr,srna,sorf,promoter,operon,circrna,
       go_term,srna_target,snp,ppi_network,localization,riboswitch_thermometer,
       crispr,merge_features,screenshot,colorize_screenshot_tracks}
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

Citation
--------

SH Yu, J. Vogel, KU Förstner. 2018, GigaScience,
`DOI:10.1093/gigascience/giy096 <https://academic.oup.com/gigascience/advance-article/doi/10.1093/gigascience/giy096/5087959>`_,
`PMID:30169674 <https://www.ncbi.nlm.nih.gov/pubmed/30169674>`_.


License
-------

`ISC <https://en.wikipedia.org/wiki/ISC_license>`__ (Internet Systems
Consortium license ~ simplified BSD license) - see `LICENSE <http://annogesic.readthedocs.io/en/latest/license.html>`__

Contact
-------

If you have any questions, please contact `Sung-Huan Yu <mailto:shyu@biochem.mpg.de>`_
