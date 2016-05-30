Pre-required tools or databases
==============

Tools
----

`Python <https://www.python.org/>`_ : version higher or equal than 3.4

`BioPython <http://biopython.org/wiki/Main_Page>`_: version higher or equal than 1.65

`BioPerl <http://www.bioperl.org/wiki/Main_Page>`_:  version higher or equal than 1.6.1

`Blast+ <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`_ : version higher or equal than 2.2.28+ (for sRNA_detection)

`ImageMagick <http://www.imagemagick.org/script/index.php>`_ : version higher or equal than 6.9.0-0 (for generate screenshot)

`Infernal <http://infernal.janelia.org/>`_ : version higher or equal than 1.1.1 (for riboswitch detection)

`Matplotlib <http://matplotlib.org/>`_ : version higher or equal than 1.5.0 (for generate all statistics figure)

`Networkx <https://networkx.github.io/>`_ : version higher or equal than 1.10 (for generate PPI figure)

`MEME <http://meme-suite.org/tools/meme>`_ : version higher or equal than 4.10.0 (for promoter detection)

`PAGIT and RATT <http://www.sanger.ac.uk/resources/software/pagit/>`_ : version higher or equal than 1.64 (for annotation transfer)
Please be attation, before you start to run RATT(annotation transfer), run ``source $PAGIT_HOME/sourceme.pagit`` first. it will
modify the path for execute RATT.

`Psortb <http://www.psort.org/psortb/>`_ : version higher or equal than 3.0 (for subcellular localization detection)

`Samtools <https://github.com/samtools>`_ : version higher or equal than 1.0 (using htslib 1.0) (for SNP calling, CircRNA_detection)

`Bcftools <https://github.com/samtools>`_ : version higher or equal than 1.0 (using htslib 1.0) (for SNP calling)

`Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_ : version higher or equal than 0.1.9 (for CircRNA_detection)

`TranstermHP <http://transterm.cbcb.umd.edu/>`_ : version higher or equal than 2.09 (for rho-independent terminator prediction)

`TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_ : version higher or equal than 1.04 (for TSS and processing site prediction, TSS optimization)

`ViennaRNA <http://www.tbi.univie.ac.at/RNA/>`_ : version higher or equal than 2.1.7 (for rho-independent terminator prediction, sRNA detection, sRNA_target detection)

`Ps2pdf14 <http://pages.cs.wisc.edu/~ghost/doc/AFPL/6.50/Ps2pdf.htm>`_ (for sRNA detection)

`IGV <https://www.broadinstitute.org/software/igv/home>`_ (for screenshot and color_png)

Databases
---------

`BSRD <http://www.bac-srna.org/BSRD/index.jsp>`_ (for sRNA detection)

`nr database <ftp://ftp.ncbi.nih.gov/blast/db/FASTA/>`_ (for sRNA detection)

`Rfam <http://rfam.xfam.org/>`_ (for riboswitch)

`idmapping_selected.tab from Uniprot <http://www.uniprot.org/downloads>`_ (for Go term detection)

`goslim.obo <http://geneontology.org/page/go-slim-and-subset-guide>`_ (for Go term detection)

`go.obo <http://geneontology.org/page/download-ontology>`_ (for Go term detection)

`species.vXXXX.txt from STRING <http://string-db.org/cgi/download.pl>`_ (for protein-protein interaction prediction)
