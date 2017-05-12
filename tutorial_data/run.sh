main(){
    PATH_FILE=$(pwd)
    PYTHON_PATH=python3
    ANNOGESIC_PATH=annogesic #### specify the path of annogesic
    ANNOGESIC_FOLDER=ANNOgesic
    FTP_SOURCE="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Campylobacter_jejuni/latest_assembly_versions/GCF_000017905.1_ASM1790v1/"
    WIG_FOLDER="ANNOgesic/input/wigs/tex_notex"
    TEX_LIBS="$WIG_FOLDER/GSM951380_Log_81116_R1_minus_TEX_in_NC_009839_minus.wig:notex:1:a:- \
              $WIG_FOLDER/GSM951381_Log_81116_R1_plus_TEX_in_NC_009839_minus.wig:tex:1:a:- \
              $WIG_FOLDER/GSM951380_Log_81116_R1_minus_TEX_in_NC_009839_plus.wig:notex:1:a:+ \
	      $WIG_FOLDER/GSM951381_Log_81116_R1_plus_TEX_in_NC_009839_plus.wig:tex:1:a:+"


    set_up_analysis_folder
#    get_wig_and_read_files
#    get_input_files    
#    update_genome_fasta
#    annotation_transfer
#    Optimize_TSSpredator
#    TSS_prediction
#    processing_site_prediction
#    Transcript_detection
#    Terminator_prediction
#    utr_detection
#    operon_detection
#    promoter_detection
#    sRNA_detection
#    sORF_detection
#    sRNA_target
#    CircRNA_detection
#    SNP_calling
#    GO_term
#    Subcellular_localization
#    PPI_network
#    riboswitch_and_RNA_thermometer
#    crispr
#    merge_features
#    gen_screenshot
#    color_png
}


set_up_analysis_folder(){
    if ! [ -d $ANNOGESIC_FOLDER ]
    then
        $ANNOGESIC_PATH create \
	-pj $ANNOGESIC_FOLDER
    fi
}

get_wig_and_read_files(){
    #### This is for downloading tutorial data.
    #### If you have your own data, please put your data in corresponding folders and skip this step.

    #### Download sratoolkit
    wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.5.2/sratoolkit.2.5.2-ubuntu64.tar.gz
    tar -zxvf sratoolkit.2.5.2-ubuntu64.tar.gz
    rm sratoolkit.2.5.2-ubuntu64.tar.gz

    #### transfer SRA files to fasta files
    for SRA in SRR515254 SRR515255 SRR515256 SRR515257
    do
        wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP013/SRP013869/$SRA/$SRA.sra
        ./sratoolkit.2.5.2-ubuntu64/bin/fastq-dump.2.5.2 --fasta $SRA.sra
        rm $SRA.sra
    done
    mv *.fasta ANNOgesic/input/reads

    #### Download and unzip wiggle files
    wget -cP ANNOgesic/input/wigs/tex_notex ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM951nnn/GSM951380/suppl/GSM951380%5FLog%5F81116%5FR1%5Fminus%5FTEX%5Fin%5FNC%5F009839%5Fminus.wig.gz
    wget -cP ANNOgesic/input/wigs/tex_notex ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM951nnn/GSM951380/suppl/GSM951380%5FLog%5F81116%5FR1%5Fminus%5FTEX%5Fin%5FNC%5F009839%5Fplus.wig.gz
    wget -cP ANNOgesic/input/wigs/tex_notex ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM951nnn/GSM951381/suppl/GSM951381%5FLog%5F81116%5FR1%5Fplus%5FTEX%5Fin%5FNC%5F009839%5Fminus.wig.gz
    wget -cP ANNOgesic/input/wigs/tex_notex ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM951nnn/GSM951381/suppl/GSM951381%5FLog%5F81116%5FR1%5Fplus%5FTEX%5Fin%5FNC%5F009839%5Fplus.wig.gz
    cd ANNOgesic/input/wigs/tex_notex
    gunzip GSM951380_Log_81116_R1_minus_TEX_in_NC_009839_minus.wig.gz \
           GSM951380_Log_81116_R1_minus_TEX_in_NC_009839_plus.wig.gz \
           GSM951381_Log_81116_R1_plus_TEX_in_NC_009839_minus.wig.gz \
           GSM951381_Log_81116_R1_plus_TEX_in_NC_009839_plus.wig.gz
    cd ../../../../

    #### Chromosome names of the wiggle files are not the same as fasta files.
    #### We can use replace_seq_id.py to modify them.
    wget https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/tutorial_data/replace_seq_id.py
    python3 replace_seq_id.py -i ANNOgesic/input/wigs/tex_notex -n NC_009839.1
    rm replace_seq_id.py
}

get_input_files(){
    $ANNOGESIC_PATH \
	get_input_files \
	-F $FTP_SOURCE \
	-g \
	-f \
	-e \
	-k \
	-p \
	-r \
	-pj $ANNOGESIC_FOLDER
}

update_genome_fasta(){
    #### The mutation.csv is only for our tutorial.
    wget https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/tutorial_data/mutation.csv
    mv mutation.csv ANNOgesic/input/mutation_table 

    $ANNOGESIC_PATH \
        update_genome_fasta \
	-c $ANNOGESIC_FOLDER/input/references/fasta_files/NC_009839.1.fa \
	-m $ANNOGESIC_FOLDER/input/mutation_table/mutation.csv \
	-pj $ANNOGESIC_FOLDER
}

annotation_transfer(){
    $ANNOGESIC_PATH \
        annotation_transfer \
	-ce $ANNOGESIC_FOLDER/input/references/annotations/NC_009839.1.embl \
	-cf $ANNOGESIC_FOLDER/input/references/fasta_files/NC_009839.1.fa \
	-uf $ANNOGESIC_FOLDER/output/updated_references/fasta_files/NC_test.1.fa \
	    $ANNOGESIC_FOLDER/output/updated_references/fasta_files/test_case2.fa \
	-e chromosome \
	-t Strain \
	-p NC_009839.1:NC_test.1 NC_009839.1:test_case2 \
	-g \
	-pj $ANNOGESIC_FOLDER
}

Optimize_TSSpredator(){
#### The manual-detected TSS file is only for the tutorial
    wget -cP ANNOgesic/input/manual_TSSs/ https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/tutorial_data/NC_009839_manual_TSS.gff

    $ANNOGESIC_PATH \
        optimize_tss_ps \
        -f $ANNOGESIC_FOLDER/input/references/fasta_files/NC_009839.1.fa \
        -g $ANNOGESIC_FOLDER/input/references/annotations/NC_009839.1.gff \
        -tl $TEX_LIBS \
        -cn TSS -s 25 \
        -m $ANNOGESIC_FOLDER/input/manual_TSSs/NC_009839_manual_TSS.gff \
        -le NC_009839.1:200000 \
        -rt all_1 \
        -pj $ANNOGESIC_FOLDER
}

TSS_prediction(){
    #### The manual-detected TSS file is only for the tutorial
    wget -cP ANNOgesic/input/manual_TSS/ https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/tutorial_data/NC_009839_manual_TSS.gff

    $ANNOGESIC_PATH \
        tss_ps \
        -f $ANNOGESIC_FOLDER/input/references/fasta_files/NC_009839.1.fa \
        -g $ANNOGESIC_FOLDER/input/references/annotations/NC_009839.1.gff \
        -tl $TEX_LIBS \
        -cn test \
        -he 0.4 \
        -rh 0.1 \
        -fa 1.7 \
        -rf 0.2 \
        -bh 0.039 \
        -ef 1.1 \
        -pf 4.5 \
        -v \
        -rt all_1 \
        -le NC_009839.1:200000 \
        -m $ANNOGESIC_FOLDER/input/manual_TSSs/NC_009839_manual_TSS.gff \
        -pj $ANNOGESIC_FOLDER
}

processing_site_prediction()
{
    $ANNOGESIC_PATH \
        tss_ps \
        -f $ANNOGESIC_FOLDER/input/references/fasta_files/NC_009839.1.fa \
        -g $ANNOGESIC_FOLDER/input/references/annotations/NC_009839.1.gff \
        -tl $TEX_LIBS \
        -cn test \
        -he 0.2 \
        -rh 0.1 \
        -fa 2.0 \
        -rf 0.5 \
        -bh 0.009 \
        -ef 1.2 \
        -pf 1.5 \
        -rt all_1 \
        -p processing_site \
        -pj $ANNOGESIC_FOLDER
}

Transcript_detection(){
    $ANNOGESIC_PATH \
        transcript \
        -g $ANNOGESIC_FOLDER/input/references/annotations/NC_009839.1.gff \
        -tl $TEX_LIBS \
        -rt all_1 \
	-cf gene CDS \
        -ct $ANNOGESIC_FOLDER/output/TSSs/gffs/NC_009839.1_TSS.gff \
        -pj $ANNOGESIC_FOLDER
}

Terminator_prediction(){
    $ANNOGESIC_PATH \
        terminator \
        -f $ANNOGESIC_FOLDER/input/references/fasta_files/NC_009839.1.fa \
        -g $ANNOGESIC_FOLDER/input/references/annotations/NC_009839.1.gff \
        -a $ANNOGESIC_FOLDER/output/transcripts/gffs/NC_009839.1_transcript.gff \
        -tl $TEX_LIBS \
        -rt all_1 \
        -pj $ANNOGESIC_FOLDER
}

utr_detection(){
    $ANNOGESIC_PATH \
        utr \
        -g $ANNOGESIC_FOLDER/input/references/annotations/NC_009839.1.gff \
        -t $ANNOGESIC_FOLDER/output/TSSs/gffs/NC_009839.1_TSS.gff \
        -a $ANNOGESIC_FOLDER/output/transcripts/gffs/NC_009839.1_transcript.gff \
        -e $ANNOGESIC_FOLDER/output/terminators/gffs/best_candidates/NC_009839.1_term.gff \
        -pj $ANNOGESIC_FOLDER
}

operon_detection(){
    $ANNOGESIC_PATH \
        operon \
        -g $ANNOGESIC_FOLDER/input/references/annotations/NC_009839.1.gff \
        -t $ANNOGESIC_FOLDER/output/TSSs/gffs/NC_009839.1_TSS.gff \
        -a $ANNOGESIC_FOLDER/output/transcripts/gffs/NC_009839.1_transcript.gff \
        -u5 $ANNOGESIC_FOLDER/output/UTRs/5UTRs/gffs/NC_009839.1_5UTR.gff \
        -u3 $ANNOGESIC_FOLDER/output/UTRs/3UTRs/gffs/NC_009839.1_3UTR.gff \
        -e $ANNOGESIC_FOLDER/output/terminators/gffs/best_candidates/NC_009839.1_term.gff \
        -pj $ANNOGESIC_FOLDER
}

promoter_detection(){
    $ANNOGESIC_PATH \
        promoter \
        -t $ANNOGESIC_FOLDER/output/TSSs/gffs/NC_009839.1_TSS.gff \
        -f $ANNOGESIC_FOLDER/input/references/fasta_files/NC_009839.1.fa \
        -w 45 2-10 \
        -pj $ANNOGESIC_FOLDER
}

sRNA_detection(){
    #### If you have no sRNA database, you can use BSRD.
    wget -cP ANNOgesic/input/databases/ https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/database/sRNA_database_BSRD.fa
    #### If you have no nr databse, please download it.
    wget -cP ANNOgesic/input/databases/ ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
    gunzip ANNOgesic/input/databases/nr.gz
    mv ANNOgesic/input/databases/nr ANNOgesic/input/databases/nr.fa


    $ANNOGESIC_PATH \
        srna \
        -d tss blast_srna sec_str blast_nr \
        -g $ANNOGESIC_FOLDER/input/references/annotations/NC_009839.1.gff \
        -t $ANNOGESIC_FOLDER/output/TSSs/gffs/NC_009839.1_TSS.gff \
        -p $ANNOGESIC_FOLDER/output/processing_sites/gffs/NC_009839.1_processing.gff \
        -a $ANNOGESIC_FOLDER/output/transcripts/gffs/NC_009839.1_transcript.gff \
        -f $ANNOGESIC_FOLDER/input/references/fasta_files/NC_009839.1.fa \
        -e $ANNOGESIC_FOLDER/output/terminators/gffs/best_candidates/NC_009839.1_term.gff \
        -pt $ANNOGESIC_FOLDER/output/promoters/NC_009839.1/MEME/promoter_motifs_NC_009839.1_allgenome_all_types_45_nt/meme.csv \
        -pn MOTIF_1 \
        -m \
        -u \
	-cs \
        -sf \
        -nf \
        -nd $ANNOGESIC_FOLDER/input/databases/nr \
	-sd $ANNOGESIC_FOLDER/input/databases/sRNA_database_BSRD \
	-tl $TEX_LIBS \
        -rt all_1 \
        -pj $ANNOGESIC_FOLDER
}

sORF_detection(){
    $ANNOGESIC_PATH \
        sorf \
        -g $ANNOGESIC_FOLDER/input/references/annotations/NC_009839.1.gff \
        -t $ANNOGESIC_FOLDER/output/TSSs/gffs/NC_009839.1_TSS.gff \
        -a $ANNOGESIC_FOLDER/output/transcripts/gffs/NC_009839.1_transcript.gff \
	-f $ANNOGESIC_FOLDER/input/references/fasta_files/NC_009839.1.fa \
        -s $ANNOGESIC_FOLDER/output/sRNAs/gffs/best_candidates/NC_009839.1_sRNA.gff \
        -tl $TEX_LIBS \
        -rt all_1 -u \
        -pj $ANNOGESIC_FOLDER
}

sRNA_target(){
    $ANNOGESIC_PATH \
        srna_target \
        -g $ANNOGESIC_FOLDER/input/references/annotations/NC_009839.1.gff \
        -f $ANNOGESIC_FOLDER/input/references/fasta_files/NC_009839.1.fa \
        -s $ANNOGESIC_FOLDER/output/sRNAs/gffs/best_candidates/NC_009839.1_sRNA.gff \
        -q NC_009839.1:36989:37044:- \
        -p both \
        -pj $ANNOGESIC_FOLDER
}

CircRNA_detection(){
    #### For tutorial, we just want to test the subcommand.
    #### Thus, we can extract the subset of reads to reduce the running time.
    #### If you are running your own data, please comment the following 5 lines 
    #### and put your read files to corresponding folder.
    for SRA in SRR515254 SRR515255 SRR515256 SRR515257
    do
        head -n 50000 ANNOgesic/input/reads/$SRA.fasta > ANNOgesic/input/reads/${SRA}_50000.fasta
        rm ANNOgesic/input/reads/$SRA.fasta
    done

    READ_FILES=$ANNOGESIC_FOLDER/input/reads/SRR515254_50000.fasta,\
$ANNOGESIC_FOLDER/input/reads/SRR515255_50000.fasta,\
$ANNOGESIC_FOLDER/input/reads/SRR515256_50000.fasta,\
$ANNOGESIC_FOLDER/input/reads/SRR515257_50000.fasta

    $ANNOGESIC_PATH \
        circrna \
	-f $ANNOGESIC_FOLDER/input/references/fasta_files/NC_009839.1.fa \
        -p 10 \
        -g $ANNOGESIC_FOLDER/input/references/annotations/NC_009839.1.gff \
        -rp all_samples:$READ_FILES \
        -pj $ANNOGESIC_FOLDER
}

SNP_calling(){
    #### This is only for tutorial.
    #### Since we already got Bam via circrna, we can put the bam files to corresponding folder
    cp ANNOgesic/output/circRNAs/segemehl_alignment_files/NC_009839.1/SRR51525* ANNOgesic/input/BAMs/BAMs_map_query_genomes/tex_notex
    BAM_FILES=$ANNOGESIC_FOLDER/input/BAMs/BAMs_map_query_genomes/tex_notex/SRR515254_50000_NC_009839.1.bam,\
$ANNOGESIC_FOLDER/input/BAMs/BAMs_map_query_genomes/tex_notex/SRR515255_50000_NC_009839.1.bam,\
$ANNOGESIC_FOLDER/input/BAMs/BAMs_map_query_genomes/tex_notex/SRR515256_50000_NC_009839.1.bam,\
$ANNOGESIC_FOLDER/input/BAMs/BAMs_map_query_genomes/tex_notex/SRR515257_50000_NC_009839.1.bam 

    $ANNOGESIC_PATH \
         snp \
	-t related_genome \
	-p with_BAQ without_BAQ extend_BAQ \
	-b all_samples:2:$BAM_FILES \
	-f $ANNOGESIC_FOLDER/input/references/fasta_files/NC_009839.1.fa \
	-pj $ANNOGESIC_FOLDER
}

GO_term(){
    #### If you have no goslim_generic.obo, go.obo, and idmapping_selected.tab, please download them.
    wget -cP ANNOgesic/input/databases http://www.geneontology.org/ontology/subsets/goslim_generic.obo
    wget -cP ANNOgesic/input/databases http://geneontology.org/ontology/go.obo
    wget -cP ANNOgesic/input/databases ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
    gunzip ANNOgesic/input/databases/idmapping_selected.tab.gz

    $ANNOGESIC_PATH \
        go_term \
	-g $ANNOGESIC_FOLDER/input/references/annotations/NC_009839.1.gff \
	-a $ANNOGESIC_FOLDER/output/transcripts/gffs/NC_009839.1_transcript.gff \
        -go $ANNOGESIC_FOLDER/input/databases/go.obo \
        -gs $ANNOGESIC_FOLDER/input/databases/goslim_generic.obo \
        -u $ANNOGESIC_FOLDER/input/databases/idmapping_selected.tab \
	-pj $ANNOGESIC_FOLDER
}

Subcellular_localization(){
    $ANNOGESIC_PATH \
        localization \
        -g $ANNOGESIC_FOLDER/input/references/annotations/NC_009839.1.gff \
        -f $ANNOGESIC_FOLDER/input/references/fasta_files/NC_009839.1.fa \
        -a $ANNOGESIC_FOLDER/output/transcripts/gffs/NC_009839.1_transcript.gff \
        -m -b negative \
        -pj $ANNOGESIC_FOLDER
}

PPI_network(){
    #### If you have no species.$Version.txt of STRING database, please download it.
    wget -cP ANNOgesic/input/databases http://string-db.org/newstring_download/species.v10.txt


    $ANNOGESIC_PATH \
        ppi_network \
	-s NC_009839.1.gff:NC_009839.1:'Campylobacter jejuni 81176':'Campylobacter jejuni' \
	-g $ANNOGESIC_FOLDER/input/references/annotations/NC_009839.1.gff \
	-d $ANNOGESIC_FOLDER/input/databases/species.v10.txt \
	-q NC_009839.1:70579:71463:+ NC_009839.1:102567:103973:+ \
	-n \
	-pj $ANNOGESIC_FOLDER
}

riboswitch_and_RNA_thermometer(){
    #### You can download the riboswitch and RNA thermometer ID list of Rfam from our Git repository.
    wget -cP ANNOgesic/input/riboswitch_ID_file/ https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/database/Rfam_riboswitch_ID.csv
    wget -cP ANNOgesic/input/RNA_thermometer_ID_file/ https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/database/Rfam_RNA_thermometer_ID.csv
    #### If you have no Rfam database, please download it.
    wget -cP ANNOgesic/input/databases ftp://ftp.ebi.ac.uk/pub/databases/Rfam/12.0/Rfam.tar.gz
    cd ANNOgesic/input/databases
    tar -zxvf Rfam.tar.gz
    rm Rfam.tar.gz
    cd ../../../


    $ANNOGESIC_PATH \
        riboswitch_thermometer \
	-g $ANNOGESIC_FOLDER/input/references/annotations/NC_009839.1.gff \
	-f $ANNOGESIC_FOLDER/input/references/fasta_files/NC_009839.1.fa \
	-ri $ANNOGESIC_FOLDER/input/riboswitch_ID_file/Rfam_riboswitch_ID.csv \
	-ti $ANNOGESIC_FOLDER/input/RNA_thermometer_ID_file/Rfam_RNA_thermometer_ID.csv \
	-R $ANNOGESIC_FOLDER/input/databases/CMs/Rfam.cm \
	-a $ANNOGESIC_FOLDER/output/transcripts/gffs/NC_009839.1_transcript.gff \
	-t $ANNOGESIC_FOLDER/output/TSSs/gffs/NC_009839.1_TSS.gff \
	-pj $ANNOGESIC_FOLDER
}

crispr(){
    $ANNOGESIC_PATH \
        crispr \
        -g $ANNOGESIC_FOLDER/input/references/annotations/NC_009839.1.gff \
        -f $ANNOGESIC_FOLDER/input/references/fasta_files/NC_009839.1.fa \
        -pj $ANNOGESIC_FOLDER
}

merge_features(){
    ALL_FEATURES="$ANNOGESIC_FOLDER/output/TSSs/gffs/NC_009839.1_TSS.gff \
                  $ANNOGESIC_FOLDER/input/references/annotations/NC_009839.1.gff \
                  $ANNOGESIC_FOLDER/output/UTRs/5UTRs/gffs/NC_009839.1_5UTR.gff \
                  $ANNOGESIC_FOLDER/output/UTRs/3UTRs/gffs/NC_009839.1_3UTR.gff \
                  $ANNOGESIC_FOLDER/output/terminators/gffs/best_candidates/NC_009839.1_term.gff \
                  $ANNOGESIC_FOLDER/output/processing_sites/gffs/NC_009839.1_processing.gff \
                  $ANNOGESIC_FOLDER/output/sRNAs/gffs/best_candidates/NC_009839.1_sRNA.gff \
                  $ANNOGESIC_FOLDER/output/sORFs/gffs/best_candidates/NC_009839.1_sORF.gff \
                  $ANNOGESIC_FOLDER/output/riboswitches/gffs/NC_009839.1_riboswitch.gff \
                  $ANNOGESIC_FOLDER/output/RNA_thermometers/gffs/NC_009839.1_RNA_thermometer.gff \
                  $ANNOGESIC_FOLDER/output/crisprs/gffs/best_candidates/NC_009839.1_CRISPR.gff"

    $ANNOGESIC_PATH \
        merge_features \
        -a $ANNOGESIC_FOLDER/output/transcripts/gffs/NC_009839.1_transcript.gff \
        -of $ALL_FEATURES \
        -op NC_009839.1 \
        -pj $ANNOGESIC_FOLDER
}

gen_screenshot(){
    $ANNOGESIC_PATH \
        screenshot \
	-mg $ANNOGESIC_FOLDER/output/TSSs/gffs/NC_009839.1_TSS.gff \
	-sg $ANNOGESIC_FOLDER/input/references/annotations/NC_009839.1.gff \
	    $ANNOGESIC_FOLDER/output/sRNAs/gffs/best_candidates/NC_009839.1_sRNA.gff \
	-f $ANNOGESIC_FOLDER/input/references/fasta_files/NC_009839.1.fa \
	-o $ANNOGESIC_FOLDER/output/TSSs \
	-tl $TEX_LIBS \
	-pj $ANNOGESIC_FOLDER

    #### We only perform 6 screenshots for each strand in order to reduce the running time for tutorial.
    #### If you are running your own data, please comment the following two lines.
    head -n 30 ANNOgesic/output/TSSs/screenshots/NC_009839.1/forward.txt > ANNOgesic/output/TSSs/screenshots/NC_009839.1/forward_6_cases.txt
    head -n 30 ANNOgesic/output/TSSs/screenshots/NC_009839.1/reverse.txt > ANNOgesic/output/TSSs/screenshots/NC_009839.1/reverse_6_cases.txt

    #### Now you can use IGV to produce screenshots.
    #### Open IGV -> Tools -> Run Batch Script -> choose forward_6_cases.txt or reverse_6_cases.txt
}

color_png(){
    $ANNOGESIC_PATH \
        colorize_screenshot_tracks \
	-t 2 \
	-f $ANNOGESIC_FOLDER/output/TSSs \
	-pj $ANNOGESIC_FOLDER
}

main
