.. image:: https://img.shields.io/pypi/v/annogesic.svg
   :target: https://pypi.python.org/pypi/ANNOgesic/
.. image:: https://img.shields.io/pypi/l/annogesic.svg
   :target: https://pypi.python.org/pypi/ANNOgesic/
.. image:: https://zenodo.org/badge/34061246.svg
   :target: https://zenodo.org/badge/latestdoi/34061246

About ANNOgesic
---------------
ANNOgesic is a modular, command-line tool that can
integrate different types of RNA-Seq data like dRNA-Seq or RNA-Seq
generated after transcript fragmentation and generates high quality
genome annotations. It can detect genes, CDSs/tRNAs/rRNAs, 
transcription starting sites (TSS) and processing sites, transcripts, 
terminators, untranslated regions (UTR) as well as small RNAs (sRNA), 
small open reading frames (sORF), circular RNAs, CRISPR related RNAs, 
riboswitches and RNA-thermometers. It can also perform RNA-RNA
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

Documentation
-------------

Documentation can be found on
`here <http://pythonhosted.org/ANNOgesic>`__.

Installation
------------

Short version (if you have all the requirements installed):

::

    $ pip3 install ANNOgesic

If you want to know the requirement, please refer to 
`Documentation <http://pythonhosted.org/ANNOgesic/>`__.

Arguments
-------------

::

    usage: annogesic [-h] [--version]
                     {create,get_input_files,update_genome_fasta,
                      annotation_transfer,tss_ps,optimize_tss_ps,
                      terminator,transcript,utr,srna,sorf,promoter,operon,
                      circrna,go_term,srna_target,snp,ppi_network,
                      localization,riboswitch_thermometer,crispr,
                      merge_features,screenshot,color_png}
                     ...
    
    positional arguments:
      {create,get_input_files,update_genome_fasta,annotation_transfer,tss_ps,
       optimize_tss_ps,terminator,transcript,utr,srna,sorf,promoter,operon,
       circrna,go_term,srna_target,snp,ppi_network,localization,
       riboswitch_thermometer,crispr,merge_features,screenshot,color_png}
                            commands
        create              Create a project
        get_input_files     Get required files. (i.e. Annotation files, fasta
                            files)
        update_genome_fasta
                            Get fasta files of query genomes if the query
                            sequences do not exist.
        annotation_transfer
                            Transfer the annotations from closed genome to the
                            target genome.
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
        go_term             Extract Go terms from Uniprot.
        srna_target         Detect sRNA-mRNA interactions.
        snp                 Detect SNP/mutation and generate potential fasta file.
        ppi_network         Detect protein-protein interactions with literature
                            supported.
        localization        Predict subcellular localization of CDSs.
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

License
-------

`ISC <https://en.wikipedia.org/wiki/ISC_license>`__ (Internet Systems
Consortium license ~ simplified BSD license) - see `LICENSE <https://pythonhosted.org/ANNOgesic/license.html>`__

Contact
-------

If you have any questions, please contact `Sung-Huan Yu <mailto:sung-huan.yu@uni-wuerzburg.de>`_
