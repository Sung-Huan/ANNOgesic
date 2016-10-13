About ANNOgesic
---------------
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
                     {create,get_input_files,get_target_fasta,annotation_transfer,tsspredator,optimize_tsspredator,terminator,transcript_assembly,utr,srna,sorf,promoter,operon,circrna,go_term,srna_target,snp,ppi_network,subcellular_localization,riboswitch_thermometer,crispr,merge_features,screenshot,color_png}
                     ...
    
    positional arguments:
      {create,get_input_files,get_target_fasta,annotation_transfer,tsspredator,optimize_tsspredator,terminator,transcript_assembly,utr,srna,sorf,promoter,operon,circrna,go_term,srna_target,snp,ppi_network,subcellular_localization,riboswitch_thermometer,crispr,merge_features,screenshot,color_png}
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
        terminator          Detect rho-independant terminators.
        transcript_assembly
                            Run transcriptome assembly for detecting transcripts.
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

License
-------

`ICSL <https://en.wikipedia.org/wiki/ISC_license>`__ (Internet Systems
Consortium license ~ simplified BSD license) - see `LICENSE <https://pythonhosted.org/ANNOgesic/license.html>`__

Contact
-------

If you have any questions, please contact `Sung-Huan Yu <mailto:sung-huan.yu@uni-wuerzburg.de>`_
