.. _Required tools or databases:

Required tools or databases
===========================

Each module has different required tools or databases. **For running a module, the 
users just need to install the specific tools for it.** For examples, running 
TSS prediction only needs TSSpredator.

Please also check our `Docker image <https://hub.docker.com/r/silasysh/annogesic/>`_.
Via Docker image, all the required tools are installed.

All the used tools are compatible with the ISC license of ANNOgesic or open source.

Tools
-----

**Basic requirement:**

	`Python <https://www.python.org/>`_ : version higher or equal to 3.4.

	`BioPython <http://biopython.org/wiki/Main_Page>`_: version higher or equal to 1.65.

	`Wget <https://www.gnu.org/software/wget>`_:  version higher or equal to 1.17.1.

	`Matplotlib <http://matplotlib.org/>`_ : version higher or equal to 1.5.0.

	`NumPy <http://www.numpy.org/>`_ : version higher or equal to 1.9.2. 

**Annotation transfer:**

	`BioPerl <http://www.bioperl.org/wiki/Main_Page>`_:  version higher or equal to 1.6.1.

	`RATT <http://www.sanger.ac.uk/resources/software/pagit/>`_ : version higher or equal to 1.64. Please be attention, before you start to run RATT (annotation transfer), run ``source $PAGIT_HOME/sourceme.pagit`` first. It will modify the path for execute RATT. If you run ANNOgesic through Docker container, you can skip this step.

**SNP calling:**

	`Samtools <https://github.com/samtools>`_ : version higher or equal to 1.3.1 (using htslib 1.3.1).

	`Bcftools <https://github.com/samtools>`_ : version higher or equal to 1.3.1 (using htslib 1.3.1).

**TSS and PS prediction:**

	`TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_ : version higher or equal to 1.06.

**TSS and PS parameter optimization:**

        `TSSpredator <http://it.inf.uni-tuebingen.de/?page_id=190>`_ : version higher or equal to 1.06.

**sRNA detection:**

	`Blast+ <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`_ : version higher or equal to 2.2.28+.

	`ViennaRNA <http://www.tbi.univie.ac.at/RNA/>`_ : version higher or equal to 2.3.2. RNAfold, mountain.pl and relplot.pl are needed for sRNA prediction.

**Terminator detection:**

	`TranstermHP <http://transterm.cbcb.umd.edu/>`_ : version higher or equal to 2.09.

	`ViennaRNA <http://www.tbi.univie.ac.at/RNA/>`_ : version higher or equal to 2.3.2. RNAfold is needed for terminator prediction.

**Promoter search:**

	`MEME <http://meme-suite.org/tools/meme>`_ : version higher or equal to 4.11.1.

	`GLAM2 <http://meme-suite.org/tools/glam2>`_ : version higher or equal to 4.11.1.

	`MPICH <http://www.mpich.org/>`_ : version higher or equal to 3.2. It is for parallel version of promoter detection.

**sRNA target prediction:**

	`ViennaRNA <http://www.tbi.univie.ac.at/RNA/>`_ : version higher or equal to 2.3.2.
RNAup, RNAplex, RNAplfold are required for executing many modules of ANNOgesic.

	`IntaRNA <https://github.com/BackofenLab/IntaRNA/>`_: version higher or equal to 2.0.4.

**Circular RNA detection:**

	`Samtools <https://github.com/samtools>`_ : version higher or equal to 1.3.1 (using htslib 1.3.1).

	`Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_ : version higher or equal to 0.1.9. When you install Segemehl, please type 'make all' instead of 'make' after running configure. Otherwise, the testrealign.x won't appear.

**Riboswitch and RNA thermometer identification:**

	`Infernal <http://infernal.janelia.org/>`_ : version higher or equal to 1.1.1.

**CRISPR detection:**

	`CRT <http://www.room220.com/crt/>`_: version higher or equal to 1.2.

**Subcellular localization prediction:**

	`Psortb <http://www.psort.org/psortb/>`_ : version higher or equal to 3.0.

**Protein-protein interaction detection:**

	`Networkx <https://networkx.github.io/>`_ : version higher or equal to 1.10.

**Generating screenshots of IGV:**

	`IGV <https://www.broadinstitute.org/software/igv/home>`_

**Colorization of screenshots:**

	`ImageMagick <http://www.imagemagick.org/script/index.php>`_ : version higher or equal to 6.9.0-0.

Databases
---------

**sRNA detection:**

	`BSRD <http://www.bac-srna.org/BSRD/index.jsp>`_ (A suggestion for sRNA prediction and can be found `here <https://github.com/Sung-Huan/ANNOgesic/tree/master/database>`_.)

	`nr database <ftp://ftp.ncbi.nih.gov/blast/db/FASTA/>`_

**Riboswitch and RNA thermometer prediction:**

	`Rfam <http://rfam.xfam.org/>`_ (Only including Rfam IDs of riboswitches and RNA thermometers -- the datasets can be found `here <https://github.com/Sung-Huan/ANNOgesic/tree/master/database>`_.)

**GO term identification:**

	`idmapping_selected.tab from Uniprot <http://www.uniprot.org/downloads>`_

	`goslim.obo <http://geneontology.org/page/go-slim-and-subset-guide>`_

	`go.obo <http://geneontology.org/page/download-ontology>`_

**Protein-protein interaction detection:**

	`species.v${VERSION}.txt from STRING <http://string-db.org/cgi/download.pl>`_ (${VERSION} represents the version number.)
