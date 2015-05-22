ANNOgesic is a transcriptome annotation pipeline for RNA-seq.

1. It include many useful modules.
2. It is easy to run these modules seperately, in order to reach your specific goal.
3. It also can be used to run whole modules from RNA-seq data to any annotation information.

usage: ANNOgesic.py [-h] [--version] [--bin_path BIN_PATH]
                    {create,get_input_files,get_target_fasta,annotation_transfer,tsspredator,optimize_tsspredator,color_png,terminator,transcript_assembly,utr,srna,sorf,promoter,operon,circrna,go_term,srna_target,snp,ppi_network,subcellular_localization,riboswitch,screenshot}
                    ...

positional arguments:
  {create,get_input_files,get_target_fasta,annotation_transfer,tsspredator,optimize_tsspredator,color_png,terminator,transcript_assembly,utr,srna,sorf,promoter,operon,circrna,go_term,srna_target,snp,ppi_network,subcellular_localization,riboswitch,screenshot}
                        commands
    create              Create a project
    get_input_files     Get required files (i.e. annotation files, fasta
                        files)
    get_target_fasta    Get target fasta.
    annotation_transfer
                        Run RATT to transfer the annotation files from
                        reference to target
    tsspredator         Run TSSpredator to predict TSSs.
    optimize_tsspredator
                        Optimize TSSpredator based on (partial)manual detect
                        one.
    color_png           Generating color screenshots of TSS or processing
                        site. It only works after running batch script. If the
                        color bands are not located at the proper positions,
                        Please increase --figure_height.
    terminator          Run TransTermHP for detect Terminators.
    transcript_assembly
                        Run Transcript for doing transcriptome assembly.
    utr                 Run UTR detection to detecting 5'UTR and 3'UTR.
    srna                Run sRNA detection to detecting sRNA candidates.
    sorf                Run sORF detection to detecting sORF candidates which
                        has expression.
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
  --bin_path BIN_PATH, -B BIN_PATH
                        path of scripts
