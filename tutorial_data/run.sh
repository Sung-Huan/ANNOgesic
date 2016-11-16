main(){
    PATH_FILE=$(pwd)
    PYTHON_PATH=python3
    ANNOGESIC_PATH=/home/silas/ANNOgesic/bin/annogesic
    ANNOGESIC_FOLDER=ANNOgesic
    FTP_SOURCE="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Campylobacter_jejuni/latest_assembly_versions/GCF_000017905.1_ASM1790v1/"
    WIG_FOLDER="ANNOgesic/input/wigs/tex_notex"
    TEX_LIBS="$WIG_FOLDER/GSM951380_Log_81116_R1_minus_TEX_in_NC_009839_minus.wig:notex:1:a:- \
              $WIG_FOLDER/GSM951381_Log_81116_R1_plus_TEX_in_NC_009839_minus.wig:tex:1:a:- \
              $WIG_FOLDER/GSM951380_Log_81116_R1_minus_TEX_in_NC_009839_plus.wig:notex:1:a:+ \
	      $WIG_FOLDER/GSM951381_Log_81116_R1_plus_TEX_in_NC_009839_plus.wig:tex:1:a:+"


#    set_up_analysis_folder
#    get_input_files    
#    get_target_fasta
    annotation_transfer
#    Optimize_TSSpredator
#    TSS_prediction
#    processing_site_prediction
#    Transcriptome_assembly
#    Terminator_prediction
#    utr_detection
#    operon_detection
#    promoter_detection
#    sRNA_detection
#    sORF_detection
#    sRNA_target
#    CircRNA_detection
#    SNP_calling_reference (optional)
#    SNP_calling_target
#    Go_term
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

get_target_fasta(){
    $ANNOGESIC_PATH \
        get_target_fasta \
	-r $ANNOGESIC_FOLDER/input/reference/fasta/NC_009839.1.fa \
	-o $ANNOGESIC_FOLDER/output/target/fasta/test_case1.fa:NC_test.1 \
	   $ANNOGESIC_FOLDER/output/target/fasta/test_case2.fa:test_case2 \
	-m $ANNOGESIC_FOLDER/input/mutation_table/mutation.csv \
	-pj $ANNOGESIC_FOLDER
}

annotation_transfer(){
    $ANNOGESIC_PATH \
        annotation_transfer \
	-re $ANNOGESIC_FOLDER/input/reference/annotation/NC_009839.1.embl \
	-rf $ANNOGESIC_FOLDER/input/reference/fasta/NC_009839.1.fa \
	-tf $ANNOGESIC_FOLDER/output/target/fasta/test_case1.fa \
	    $ANNOGESIC_FOLDER/output/target/fasta/test_case2.fa \
	-e chromosome \
	-t Strain \
	-p NC_009839.1:NC_test.1 NC_009839.1:test_case2 \
	-g \
	-pj $ANNOGESIC_FOLDER
}

Optimize_TSSpredator(){
    $ANNOGESIC_PATH \
        optimize_tss_processing \
        -fs $ANNOGESIC_FOLDER/output/target/fasta/NC_009839.1.fa \
        -g $ANNOGESIC_FOLDER/output/target/annotation/NC_009839.1.gff \
        -n NC_009839.1 \
        -tl $TEX_LIBS \
        -p TSS -s 25 \
        -m $ANNOGESIC_FOLDER/input/manual_TSS/NC_009839_manual_TSS.gff \
        -le 200000 \
        -pj $ANNOGESIC_FOLDER
}

TSS_prediction(){
    $ANNOGESIC_PATH \
        tss_processing \
        -f $ANNOGESIC_FOLDER/output/target/fasta/NC_009839.1.fa \
        -g $ANNOGESIC_FOLDER/output/target/annotation/NC_009839.1.gff \
        -tl $TEX_LIBS \
        -p test \
        -he 0.4 \
        -rh 0.1 \
        -fa 1.7 \
        -rf 0.2 \
        -bh 0.039 \
        -ef 1.1 \
        -pf 4.5 \
        -s \
        -v \
        -le 200000 \
        -m $ANNOGESIC_FOLDER/input/manual_TSS/NC_009839_manual_TSS.gff \
        -pj $ANNOGESIC_FOLDER
}

processing_site_prediction()
{
    $ANNOGESIC_PATH \
        tss_processing \
        -f $ANNOGESIC_FOLDER/output/target/fasta/NC_009839.1.fa \
        -g $ANNOGESIC_FOLDER/output/target/annotation/NC_009839.1.gff \
        -tl $TEX_LIBS \
        -p test \
        -he 0.2 \
        -rh 0.1 \
        -fa 2.0 \
        -rf 0.5 \
        -bh 0.009 \
        -ef 1.2 \
        -pf 1.5 \
        -s \
        -t processing_site \
        -pj $ANNOGESIC_FOLDER
}

Transcriptome_assembly(){
    $ANNOGESIC_PATH \
        transcriptome_assembly \
        -g $ANNOGESIC_FOLDER/output/target/annotation/NC_009839.1.gff \
        -tl $TEX_LIBS \
        -rt all_1 \
	-cf gene CDS \
        -ct $ANNOGESIC_FOLDER/output/TSS/gffs/NC_009839.1_TSS.gff \
        -pj $ANNOGESIC_FOLDER
}

Terminator_prediction(){
    $ANNOGESIC_PATH \
        terminator \
        -f $ANNOGESIC_FOLDER/output/target/fasta/NC_009839.1.fa \
        -g $ANNOGESIC_FOLDER/output/target/annotation/NC_009839.1.gff \
        -s \
        -a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs/NC_009839.1_transcript.gff \
        -tl $TEX_LIBS \
        -rt all_1 -tb \
        -pj $ANNOGESIC_FOLDER
}

utr_detection(){
    $ANNOGESIC_PATH \
        utr \
        -g $ANNOGESIC_FOLDER/output/target/annotation/NC_009839.1.gff \
        -t $ANNOGESIC_FOLDER/output/TSS/gffs/NC_009839.1_TSS.gff \
        -a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs/NC_009839.1_transcript.gff \
        -e $ANNOGESIC_FOLDER/output/terminator/gffs/best/NC_009839.1_term.gff \
        -pj $ANNOGESIC_FOLDER
}

operon_detection(){
    $ANNOGESIC_PATH \
        operon \
        -g $ANNOGESIC_FOLDER/output/target/annotation/NC_009839.1.gff \
        -t $ANNOGESIC_FOLDER/output/TSS/gffs/NC_009839.1_TSS.gff \
        -a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs/NC_009839.1_transcript.gff \
        -u5 $ANNOGESIC_FOLDER/output/UTR/5UTR/gffs/NC_009839.1_5UTR.gff \
        -u3 $ANNOGESIC_FOLDER/output/UTR/3UTR/gffs/NC_009839.1_3UTR.gff \
        -e $ANNOGESIC_FOLDER/output/terminator/gffs/best/NC_009839.1_term.gff \
        -s -c \
        -pj $ANNOGESIC_FOLDER
}

promoter_detection(){
    $ANNOGESIC_PATH \
        promoter \
        -t $ANNOGESIC_FOLDER/output/TSS/gffs/NC_009839.1_TSS.gff \
        -f $ANNOGESIC_FOLDER/output/target/fasta/NC_009839.1.fa \
        -w 45 2-10 \
        -pj $ANNOGESIC_FOLDER
}

sRNA_detection(){
    $ANNOGESIC_PATH \
        srna \
        -d tss blast_srna sec_str \
	--rnafold_path /home/silas/ANNOgesic/tools/ViennaRNA-2.1.7/Progs/RNAfold \
	--relplot_path /home/silas/ANNOgesic/tools/ViennaRNA-2.1.7/Utils/relplot.pl \
	--mountain_path /home/silas/ANNOgesic/tools/ViennaRNA-2.1.7/Utils/mountain.pl \
        -g $ANNOGESIC_FOLDER/output/target/annotation/NC_009839.1.gff \
        -t $ANNOGESIC_FOLDER/output/TSS/gffs/NC_009839.1_TSS.gff \
        -p $ANNOGESIC_FOLDER/output/processing_site/gffs/NC_009839.1_processing.gff \
        -a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs/NC_009839.1_transcript.gff \
        -f $ANNOGESIC_FOLDER/output/target/fasta/NC_009839.1.fa \
        -tf $ANNOGESIC_FOLDER/output/terminator/gffs/best/NC_009839.1_term.gff \
        -pt $ANNOGESIC_FOLDER/output/promoter_analysis/NC_009839.1/promoter_motifs_NC_009839.1_allstrain_all_types_45_nt/meme.csv \
        -pn MOTIF_1 \
        -m \
        -u \
	-sd $ANNOGESIC_FOLDER/input/database/sRNA_database_BSRD \
	-tl $TEX_LIBS \
        -rt all_1 \
        -pj $ANNOGESIC_FOLDER
}

#-nd $ANNOGESIC_FOLDER/input/database/nr \
#-nf \ 
#-sf \
sORF_detection(){
    $ANNOGESIC_PATH \
        sorf \
        -g $ANNOGESIC_FOLDER/output/target/annotation/NC_009839.1.gff \
        -t $ANNOGESIC_FOLDER/output/TSS/gffs/NC_009839.1_TSS.gff \
        -a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs/NC_009839.1_transcript.gff \
	-f $ANNOGESIC_FOLDER/output/target/fasta/NC_009839.1.fa \
        -s $ANNOGESIC_FOLDER/output/sRNA/gffs/best/NC_009839.1_sRNA.gff \
        -tl $TEX_LIBS \
        -rt all_1 -u \
        -pj $ANNOGESIC_FOLDER
}

sRNA_target(){
    $ANNOGESIC_PATH \
        srna_target \
        -g $ANNOGESIC_FOLDER/output/target/annotation/NC_009839.1.gff \
        -f $ANNOGESIC_FOLDER/output/target/fasta/NC_009839.1.fa \
        -r $ANNOGESIC_FOLDER/output/sRNA/gffs/best/NC_009839.1_sRNA.gff \
        -q NC_009839.1:36954:37044:- \
        -p both \
        -pj $ANNOGESIC_FOLDER
}

CircRNA_detection(){
    $ANNOGESIC_PATH \
        circrna \
	-f $ANNOGESIC_FOLDER/output/target/fasta/NC_009839.1.fa \
        -p 10 \
        -g $ANNOGESIC_FOLDER/output/target/annotation/NC_009839.1.gff \
	-a \
        -rp $ANNOGESIC_FOLDER/input/reads/SRR515254_50000.fasta \
	    $ANNOGESIC_FOLDER/input/reads/SRR515255_50000.fasta \
	    $ANNOGESIC_FOLDER/input/reads/SRR515256_50000.fasta \
	    $ANNOGESIC_FOLDER/input/reads/SRR515257_50000.fasta \
        -pj $ANNOGESIC_FOLDER
}

SNP_calling_reference(){
    $ANNOGESIC_PATH \
         snp \
	-t reference \
	-p with_BAQ without_BAQ extend_BAQ \
        -ms 1 \
	-b $ANNOGESIC_FOLDER/input/BAMs/BAMs_map_reference/tex_notex/SRR515254_50000_NC_009839.1.bam \
	   $ANNOGESIC_FOLDER/input/BAMs/BAMs_map_reference/tex_notex/SRR515255_50000_NC_009839.1.bam \
	   $ANNOGESIC_FOLDER/input/BAMs/BAMs_map_reference/tex_notex/SRR515256_50000_NC_009839.1.bam \
	   $ANNOGESIC_FOLDER/input/BAMs/BAMs_map_reference/tex_notex/SRR515257_50000_NC_009839.1.bam \
	-f $ANNOGESIC_FOLDER/output/reference/fasta/NC_009839.1.fa \
	-pj $ANNOGESIC_FOLDER
}

SNP_calling_target(){
    $ANNOGESIC_PATH \
        snp \
	-t target \
	-p with_BAQ without_BAQ extend_BAQ \
        -ms 1 \
	-b $ANNOGESIC_FOLDER/input/BAMs/BAMs_map_target/tex_notex/SRR515254_50000_NC_009839.1.bam \
	   $ANNOGESIC_FOLDER/input/BAMs/BAMs_map_target/tex_notex/SRR515255_50000_NC_009839.1.bam \
	   $ANNOGESIC_FOLDER/input/BAMs/BAMs_map_target/tex_notex/SRR515256_50000_NC_009839.1.bam \
	   $ANNOGESIC_FOLDER/input/BAMs/BAMs_map_target/tex_notex/SRR515257_50000_NC_009839.1.bam \
	-f $ANNOGESIC_FOLDER/output/target/fasta/NC_009839.1.fa \
	-pj $ANNOGESIC_FOLDER
}

Go_term(){
    $ANNOGESIC_PATH \
        go_term \
	-g $ANNOGESIC_FOLDER/output/target/annotation/NC_009839.1.gff \
	-a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs/NC_009839.1_transcript.gff \
	-pj $ANNOGESIC_FOLDER
}

Subcellular_localization(){
    $ANNOGESIC_PATH \
        subcellular_localization \
        -g $ANNOGESIC_FOLDER/output/target/annotation/NC_009839.1.gff \
        -f $ANNOGESIC_FOLDER/output/target/fasta/NC_009839.1.fa \
        -a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs/NC_009839.1_transcript.gff \
        -m -b negative \
        -pj $ANNOGESIC_FOLDER
}

PPI_network(){
    $ANNOGESIC_PATH \
        ppi_network \
	-s NC_009839.1.gff:NC_009839.1:'Campylobacter jejuni 81176':'Campylobacter jejuni' \
	-g $ANNOGESIC_FOLDER/output/target/annotation/NC_009839.1.gff \
	-d $ANNOGESIC_FOLDER/input/database/species.v10.txt \
	-q NC_009839.1:70579:71463:+ NC_009839.1:102567:103973:+ \
	-n \
	-pj $ANNOGESIC_FOLDER
}

riboswitch_and_RNA_thermometer(){
    $ANNOGESIC_PATH \
        riboswitch_thermometer \
	-g $ANNOGESIC_FOLDER/output/target/annotation/NC_009839.1.gff \
	-f $ANNOGESIC_FOLDER/output/target/fasta/NC_009839.1.fa \
	-ri $ANNOGESIC_FOLDER/input/riboswitch_ID/Rfam_riboswitch_ID.csv \
	-ti $ANNOGESIC_FOLDER/input/RNA_thermometer_ID/Rfam_RNA_thermometer_ID.csv \
	-R $ANNOGESIC_FOLDER/input/database/CMs/Rfam.cm \
	-a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs/NC_009839.1_transcript.gff \
	-t $ANNOGESIC_FOLDER/output/TSS/gffs/NC_009839.1_TSS.gff \
	-pj $ANNOGESIC_FOLDER
}

crispr(){
    $ANNOGESIC_PATH \
        crispr \
        -g $ANNOGESIC_FOLDER/output/target/annotation/NC_009839.1.gff \
        -f $ANNOGESIC_FOLDER/output/target/fasta/NC_009839.1.fa \
        -pj $ANNOGESIC_FOLDER
}

merge_features(){
    ALL_FEATURES="$ANNOGESIC_FOLDER/output/TSS/gffs/NC_009839.1_TSS.gff \
                  $ANNOGESIC_FOLDER/output/target/annotation/NC_009839.1.gff \
                  $ANNOGESIC_FOLDER/output/UTR/5UTR/gffs/NC_009839.1_5UTR.gff \
                  $ANNOGESIC_FOLDER/output/UTR/3UTR/gffs/NC_009839.1_3UTR.gff \
                  $ANNOGESIC_FOLDER/output/terminator/gffs/best/NC_009839.1_term.gff \
                  $ANNOGESIC_FOLDER/output/processing_site/gffs/NC_009839.1_processing.gff \
                  $ANNOGESIC_FOLDER/output/sRNA/gffs/best/NC_009839.1_sRNA.gff \
                  $ANNOGESIC_FOLDER/output/sORF/gffs/best/NC_009839.1_sORF.gff \
                  $ANNOGESIC_FOLDER/output/riboswitch/gffs/NC_009839.1_riboswitch.gff \
                  $ANNOGESIC_FOLDER/output/crispr/gffs/best/NC_009839.1_CRISPR.gff"

    $ANNOGESIC_PATH \
        merge_features \
        -a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs/NC_009839.1_transcript.gff \
        -of $ALL_FEATURES \
        -s NC_009839.1 \
        -pj $ANNOGESIC_FOLDER
}

gen_screenshot(){
    $ANNOGESIC_PATH \
        screenshot \
	-mg $ANNOGESIC_FOLDER/output/TSS/gffs/NC_009839.1_TSS.gff \
	-sg $ANNOGESIC_FOLDER/output/target/annotation/NC_009839.1.gff \
	    $ANNOGESIC_FOLDER/output/sRNA/gffs/best/NC_009839.1_sRNA.gff \
	-f $ANNOGESIC_FOLDER/output/target/fasta/NC_009839.1.fa \
	-o $ANNOGESIC_FOLDER/output/TSS \
	-tl $TEX_LIBS \
	-pj $ANNOGESIC_FOLDER
}

color_png(){
    $ANNOGESIC_PATH \
        color_png \
	-t 2 \
	-f $ANNOGESIC_FOLDER/output/TSS \
	-pj $ANNOGESIC_FOLDER
}

main
