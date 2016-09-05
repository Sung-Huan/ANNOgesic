main(){
    PATH_FILE=$(pwd)
    PYTHON_PATH=python3
    ANNOGESIC_PATH=annogesic
    ANNOGESIC_FOLDER=ANNOgesic
    FTP_SOURCE=ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Helicobacter_pylori_26695_uid57787/
    PAGIT_HOME=/home/silas/ANNOgesic/tools/PAGIT
    tex_notex_libs="GSM1649587_Hp26695_ML_B1_HS1_-TEX_forward.wig:notex:1:a:+,\
GSM1649587_Hp26695_ML_B1_HS1_-TEX_reverse.wig:notex:1:a:-,\
GSM1649588_Hp26695_ML_B1_HS1_+TEX_forward.wig:tex:1:a:+,\
GSM1649588_Hp26695_ML_B1_HS1_+TEX_reverse.wig:tex:1:a:-,\
GSM1649589_Hp26695_ML_B2_HS2_-TEX_forward.wig:notex:1:b:+,\
GSM1649589_Hp26695_ML_B2_HS2_-TEX_reverse.wig:notex:1:b:-,\
GSM1649590_Hp26695_ML_B2_HS2_+TEX_forward.wig:tex:1:b:+,\
GSM1649590_Hp26695_ML_B2_HS2_+TEX_reverse.wig:tex:1:b:-,\
GSM1649591_Hp26695_ML_B3_HS2_-TEX_forward.wig:notex:1:c:+,\
GSM1649591_Hp26695_ML_B3_HS2_-TEX_reverse.wig:notex:1:c:-,\
GSM1649592_Hp26695_ML_B3_HS2_+TEX_forward.wig:tex:1:c:+,\
GSM1649592_Hp26695_ML_B3_HS2_+TEX_reverse.wig:tex:1:c:-"

    frag_libs="ML-Lib_div_by_14583533.0_multi_by_14583533.0_forward.wig:frag:1:a:+,\
ML-Lib_div_by_14583533.0_multi_by_14583533.0_reverse.wig:frag:1:a:-"

#    set_up_analysis_folder
#    get_input_files    
#    get_target_fasta
#    annotation_transfer
#    expression_analysis
#    SNP_calling_reference
#    TSS_prediction
#    Transcriptome_assembly
#    Terminator_prediction
#    processing_site_prediction
#    utr_detection
#    sRNA_detection
#    sORF_detection
#    promoter_detection
#    CircRNA_detection
#    Go_term
#    sRNA_target    
#    operon_detection
#    SNP_calling_target
#    PPI_network
#    Subcellular_localization
#    riboswitch_and_RNA_thermometer
#    Optimize_TSSpredator
#    gen_screenshot
#    color_png
#    merge_features
}


create_folders(){
    for FOLDER in bin
    do
        if ! [ -d $FOLDER ]
        then
            mkdir -p $FOLDER
        fi
    done
}

set_up_analysis_folder(){
    if ! [ -d $ANNOGESIC_FOLDER ]
    then
        $ANNOGESIC_PATH create $ANNOGESIC_FOLDER
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
	$ANNOGESIC_FOLDER
}

get_target_fasta(){
    $ANNOGESIC_PATH \
        get_target_fasta \
	-r $ANNOGESIC_FOLDER/input/reference/fasta \
	-o test_case1:NC_test.1,test_case2:test_case2 \
	-m $ANNOGESIC_FOLDER/input/mutation_table/mutation.csv \
	$ANNOGESIC_FOLDER
}

annotation_transfer(){
    # instead of using "source"
    . $PAGIT_HOME/sourceme.pagit
    $ANNOGESIC_PATH \
        annotation_transfer \
	-re $ANNOGESIC_FOLDER/input/reference/annotation \
	-rf $ANNOGESIC_FOLDER/input/reference/fasta \
	-tf $ANNOGESIC_FOLDER/output/target/fasta \
	-e chromosome \
	-t Strain \
	-p NC_000915.1:NC_test.1,NC_000915.1:test_case2 \
	-g \
	$ANNOGESIC_FOLDER
}

SNP_calling_reference(){
    $ANNOGESIC_PATH \
         snp \
	-t reference \
	-p 1,2,3 \
	-tw $ANNOGESIC_FOLDER/input/BAMs/BAMs_map_reference/tex_notex \
	-f $ANNOGESIC_FOLDER/output/reference/fasta \
	$ANNOGESIC_FOLDER
}

TSS_prediction(){
    $ANNOGESIC_PATH \
        tsspredator \
	annogesic tsspredator \
	-w $ANNOGESIC_FOLDER/input/wigs/tex_notex \
	-f $ANNOGESIC_FOLDER/output/target/fasta \
	-g $ANNOGESIC_FOLDER/output/target/annotation \
	-l $tex_notex_libs \
	-p test \
	-he 0.3 \
	-rh 0.2 \
	-fa 2.0 \
	-rf 0.5 \
	-bh 0.0 \
	-ef 1.7 \
	-pf 1.5 \
	-s \
	-v \
	-le 200000 \
	-m $ANNOGESIC_FOLDER/input/manual_TSS/NC_000915_manual_TSS.gff \
	$ANNOGESIC_FOLDER
}

Transcriptome_assembly(){
    $ANNOGESIC_PATH \
        transcript_assembly \
	-g $ANNOGESIC_FOLDER/output/target/annotation \
	-tw $ANNOGESIC_FOLDER/input/wigs/tex_notex \
	-tl $tex_notex_libs \
	-rt 1 \
	-ct $ANNOGESIC_FOLDER/output/TSS/gffs \
	-cg $ANNOGESIC_FOLDER/output/target/annotation \
	$ANNOGESIC_FOLDER
}

Terminator_prediction(){
    $ANNOGESIC_PATH \
        terminator \
	-f $ANNOGESIC_FOLDER/output/target/fasta \
	-g $ANNOGESIC_FOLDER/output/target/annotation \
	-s \
	-tw $ANNOGESIC_FOLDER/input/wigs/tex_notex \
	-a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs \
	-tl $tex_notex_libs \
	-rt 1 -tb \
	$ANNOGESIC_FOLDER
}

processing_site_prediction()
{
    $ANNOGESIC_PATH \
        tsspredator \
	annogesic tsspredator \
	-w $ANNOGESIC_FOLDER/input/wigs/tex_notex \
	-f $ANNOGESIC_FOLDER/output/target/fasta \
	-g $ANNOGESIC_FOLDER/output/target/annotation \
	-l $tex_notex_libs \
	-p test \
	-he 0.3 \
	-rh 0.2 \
	-fa 2.0 \
	-rf 0.5 \
	-bh 0.0 \
	-ef 1.9 \
	-pf 5.7 \
	-s \
	-v \
	-t processing_site \
	$ANNOGESIC_FOLDER
}

utr_detection(){
    $ANNOGESIC_PATH \
        utr \
	-g $ANNOGESIC_FOLDER/output/target/annotation \
	-t $ANNOGESIC_FOLDER/output/TSS/gffs \
	-a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs \
	-e $ANNOGESIC_FOLDER/output/terminator/gffs/best \
	$ANNOGESIC_FOLDER
}

sRNA_detection(){
    $ANNOGESIC_PATH \
        srna \
	-d tss,blast_srna,promoter,term,blast_nr,sec_str \
	-g $ANNOGESIC_FOLDER/output/target/annotation \
	-t $ANNOGESIC_FOLDER/output/TSS/gffs \
	-p $ANNOGESIC_FOLDER/output/processing_site/gffs \
	-a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs \
	-tw $ANNOGESIC_FOLDER/input/wigs/tex_notex \
	-f $ANNOGESIC_FOLDER/output/target/fasta \
	-tf $ANNOGESIC_FOLDER/output/terminator/gffs/best \
	-pt $ANNOGESIC_FOLDER/output/promoter_analysis/NC_000915.1/promoter_motifs_NC_000915.1_allstrain_all_types_45_nt/meme.csv \
	-pn MOTIF_1 \
	-m \
	-u \
	-nf \
	-sf \
	-sd $ANNOGESIC_FOLDER/input/database/sRNA_database_BSRD \
	-nd $ANNOGESIC_FOLDER/input/database/nr \
	-tl $tex_notex_libs \
	-rt 1 \
	-ba \
	$ANNOGESIC_FOLDER
}

sORF_detection(){
    $ANNOGESIC_PATH \
        sorf \
	-g $ANNOGESIC_FOLDER/output/target/annotation \
	-t $ANNOGESIC_FOLDER/output/TSS/gffs \
	-a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs \
	-tw $ANNOGESIC_FOLDER/input/wigs/tex_notex \
	-f $ANNOGESIC_FOLDER/output/target/fasta \
	-s $ANNOGESIC_FOLDER/output/sRNA/gffs/best \
	-tl $tex_notex_libs \
	-rt 1 -u \
	$ANNOGESIC_FOLDER
}

promoter_detection(){
    $ANNOGESIC_PATH \
        promoter \
	-t $ANNOGESIC_FOLDER/output/TSS/gffs \
	-f $ANNOGESIC_FOLDER/output/target/fasta \
	-w 45,2-10 \
	$ANNOGESIC_FOLDER
}

CircRNA_detection(){
    $ANNOGESIC_PATH \
        circrna \
	-f $ANNOGESIC_FOLDER/output/target/fasta \
	-p 10 \
	-g $ANNOGESIC_FOLDER/output/target/annotation \
	-a \
	$ANNOGESIC_FOLDER	
}

Go_term(){
    $ANNOGESIC_PATH \
        go_term \
	-g $ANNOGESIC_FOLDER/output/target/annotation \
	-a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs \
	$ANNOGESIC_FOLDER
}

sRNA_target(){
    $ANNOGESIC_PATH \
        srna_target \
	-g $ANNOGESIC_FOLDER/output/target/annotation \
	-f $ANNOGESIC_FOLDER/output/target/fasta \
	-r $ANNOGESIC_FOLDER/output/sRNA/gffs/best \
	-q NC_000915.1:-:7249:7321 \
	-p both \
	$ANNOGESIC_FOLDER
}

operon_detection(){
    $ANNOGESIC_PATH \
        operon \
	-g $ANNOGESIC_FOLDER/output/target/annotation \
	-t $ANNOGESIC_FOLDER/output/TSS/gffs \
	-a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs \
	-u5 $ANNOGESIC_FOLDER/output/UTR/5UTR/gffs \
	-u3 $ANNOGESIC_FOLDER/output/UTR/3UTR/gffs \
	-e $ANNOGESIC_FOLDER/output/terminator/gffs/best \
	-s -c \
	$ANNOGESIC_FOLDER
}

SNP_calling_target(){
    $ANNOGESIC_PATH \
        snp \
	-t target \
	-p 1,2,3 \
	-tw $ANNOGESIC_FOLDER/input/BAMs/BAMs_map_target/tex_notex \
	-f $ANNOGESIC_FOLDER/output/target/fasta \
	$ANNOGESIC_FOLDER
}

PPI_network(){
    $ANNOGESIC_PATH \
        ppi_network \
	-s NC_000915.1.gff:NC_000915.1:'Helicobacter pylori 26695':'Helicobacter pylori' \
	-g $ANNOGESIC_FOLDER/output/target/annotation \
	-d $ANNOGESIC_FOLDER/input/database/species.v10.txt \
	-q NC_000915.1:217:633:-,NC_000915.1:2719:3402:+ \
	-n \
	$ANNOGESIC_FOLDER
}

Subcellular_localization(){
    $ANNOGESIC_PATH \
        subcellular_localization \
	-g $ANNOGESIC_FOLDER/output/target/annotation \
	-f $ANNOGESIC_FOLDER/output/target/fasta \
	-a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs \
	-m -b negative \
	$ANNOGESIC_FOLDER
}

riboswitch_and_RNA_thermometer(){
    $ANNOGESIC_PATH \
        riboswitch_thermometer \
	-g $ANNOGESIC_FOLDER/output/target/annotation \
	-f $ANNOGESIC_FOLDER/output/target/fasta \
	-ri $ANNOGESIC_FOLDER/input/riboswitch_ID/Rfam_riboswitch_ID.csv \
	-ri $ANNOGESIC_FOLDER/input/RNA_thermometer_ID/Rfam_RNA_thermometer_ID.csv \
	-R $ANNOGESIC_FOLDER/input/database/CMs/Rfam.cm \
	-a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs \
	$ANNOGESIC_FOLDER
}

Optimize_TSSpredator(){
    $ANNOGESIC_PATH \
        optimize_tsspredator \
	-w $ANNOGESIC_FOLDER/input/wigs/tex_notex \
	-fs $ANNOGESIC_FOLDER/output/target/fasta \
	-g $ANNOGESIC_FOLDER/output/target/annotation \
	-n NC_000915.1 \
	-l $tex_notex_libs \
	-p TSS -s 25 \
	-m $ANNOGESIC_FOLDER/input/manual_TSS/NC_000915_manual_TSS.gff \
	-le 200000 \
	$ANNOGESIC_FOLDER
}

gen_screenshot(){
    $ANNOGESIC_PATH \
        screenshot \
	-mg $ANNOGESIC_FOLDER/output/TSS/gffs/NC_000915.1_TSS.gff \
	-sg $ANNOGESIC_FOLDER/output/target/annotation/NC_000915.1.gff,ANNOgesic/output/sRNA/gffs/best/NC_000915.1_sRNA.gff \
	-f $ANNOGESIC_FOLDER/output/target/fasta/NC_000915.1.fa \
	-o $ANNOGESIC_FOLDER/output/TSS \
	-tl $tex_notex_libs \
	-tw $ANNOGESIC_FOLDER/input/wigs/tex_notex \
	$ANNOGESIC_FOLDER
}

color_png(){
    $ANNOGESIC_PATH \
        color_png \
	-t 2 \
	-f $ANNOGESIC_FOLDER/output/TSS \
	$ANNOGESIC_FOLDER
}

merge_features(){
    ALL_FEATURES=$ANNOGESIC_FOLDER/output/TSS/gffs/NC_000915.1_TSS.gff,\
$ANNOGESIC_FOLDER/output/target/annotation/NC_000915.1.gff,\
$ANNOGESIC_FOLDER/output/UTR/5UTR/gffs/NC_000915.1_5UTR.gff,\
$ANNOGESIC_FOLDER/output/UTR/3UTR/gffs/NC_000915.1_3UTR.gff,\
$ANNOGESIC_FOLDER/output/terminator/gffs/best/NC_000915.1_term.gff,\
$ANNOGESIC_FOLDER/output/processing_site/gffs/NC_000915.1_processing.gff,\
$ANNOGESIC_FOLDER/output/sRNA/gffs/best/NC_000915.1_sRNA.gff,\
$ANNOGESIC_FOLDER/output/sORF/gffs/best/NC_000915.1_sORF.gff,\
$ANNOGESIC_FOLDER/output/riboswitch/gffs/NC_000915.1_riboswitch.gff,\
$ANNOGESIC_FOLDER/output/crispr/gffs/best/NC_000915.1_CRISPR.gff

    $ANNOGESIC_PATH \
	merge_features \
        -a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs/NC_000915.1_transcript.gff \
        -of $ALL_FEATURES \
	-s NC_000915.1 \
	$ANNOGESIC_FOLDER
}

main
