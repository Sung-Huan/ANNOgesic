main(){
    PATH_FILE=$(pwd)
    PYTHON_PATH=python3.4
    STRAINS=Staphylococcus_aureus_HG003
    TRANSAP_PATH=/home/silas/ANNO/Transap.py
    TRANSAP_FOLDER=Transap
    FTP_SOURCE=ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Staphylococcus_aureus_NCTC_8325_uid57795
    BIN_PATH=/home/silas/ANNO/bin_bak
    ID_MAPPING_SOURCE=ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz

    #######################################################################################################################
    # if you have many GFF or FASTA files, you want to combine them to be one file.                                       #
    # you can un-mark COMBINE_GFF, COMBINE_FASTA.                                                                         #
    # Additionally, please give them a name.                                                                              #
    # The name should be absolutive path.                                                                                 #
    #######################################################################################################################
    tex_notex_libs="TSB_OD_0.2_TEX_div_by_3598556.0_multi_by_3382258.0_reverse.wig:tex:1:a:- \
                   TSB_OD_0.5_TEX_div_by_4420442.0_multi_by_3382258.0_reverse.wig:tex:2:a:- \
                   TSB_OD_1_TEX_div_by_3956047.0_multi_by_3382258.0_reverse.wig:tex:3:a:- \
                   TSB_ON_TEX_div_by_4836916.0_multi_by_3382258.0_reverse.wig:tex:4:a:- \
                   TSB_t0_TEX_div_by_4817602.0_multi_by_3382258.0_reverse.wig:tex:5:a:- \
                   TSB_t1_TEX_div_by_4957924.0_multi_by_3382258.0_reverse.wig:tex:6:a:- \
                   TSB_t2_TEX_div_by_5024086.0_multi_by_3382258.0_reverse.wig:tex:7:a:- \
                   pMEM_OD_0.2_TEX_div_by_4581966.0_multi_by_3382258.0_reverse.wig:tex:8:a:- \
                   pMEM_OD_0.5_TEX_div_by_4968210.0_multi_by_3382258.0_reverse.wig:tex:9:a:- \
                   pMEM_OD_1_TEX_div_by_4719860.0_multi_by_3382258.0_reverse.wig:tex:10:a:- \
                   pMEM_ON_TEX_div_by_5454872.0_multi_by_3382258.0_reverse.wig:tex:11:a:- \
                   pMEM_t0_TEX_div_by_4872810.0_multi_by_3382258.0_reverse.wig:tex:12:a:- \
                   pMEM_t1_TEX_div_by_4637206.0_multi_by_3382258.0_reverse.wig:tex:13:a:- \
                   pMEM_t2_TEX_div_by_4782742.0_multi_by_3382258.0_reverse.wig:tex:14:a:- \
                   TSB_OD_0.2_TEX_div_by_3598556.0_multi_by_3382258.0_forward.wig:tex:1:a:+ \
                   TSB_OD_0.5_TEX_div_by_4420442.0_multi_by_3382258.0_forward.wig:tex:2:a:+ \
                   TSB_OD_1_TEX_div_by_3956047.0_multi_by_3382258.0_forward.wig:tex:3:a:+ \
                   TSB_ON_TEX_div_by_4836916.0_multi_by_3382258.0_forward.wig:tex:4:a:+ \
                   TSB_t0_TEX_div_by_4817602.0_multi_by_3382258.0_forward.wig:tex:5:a:+ \
                   TSB_t1_TEX_div_by_4957924.0_multi_by_3382258.0_forward.wig:tex:6:a:+ \
                   TSB_t2_TEX_div_by_5024086.0_multi_by_3382258.0_forward.wig:tex:7:a:+ \
                   pMEM_OD_0.2_TEX_div_by_4581966.0_multi_by_3382258.0_forward.wig:tex:8:a:+ \
                   pMEM_OD_0.5_TEX_div_by_4968210.0_multi_by_3382258.0_forward.wig:tex:9:a:+ \
                   pMEM_OD_1_TEX_div_by_4719860.0_multi_by_3382258.0_forward.wig:tex:10:a:+ \
                   pMEM_ON_TEX_div_by_5454872.0_multi_by_3382258.0_forward.wig:tex:11:a:+ \
                   pMEM_t0_TEX_div_by_4872810.0_multi_by_3382258.0_forward.wig:tex:12:a:+ \
                   pMEM_t1_TEX_div_by_4637206.0_multi_by_3382258.0_forward.wig:tex:13:a:+ \
                   pMEM_t2_TEX_div_by_4782742.0_multi_by_3382258.0_forward.wig:tex:14:a:+ \
                   TSB_OD_0.2_div_by_5033297.0_multi_by_3382258.0_reverse.wig:notex:1:a:- \
                   TSB_OD_0.5_div_by_4665048.0_multi_by_3382258.0_reverse.wig:notex:2:a:- \
                   TSB_OD_1_div_by_4443620.0_multi_by_3382258.0_reverse.wig:notex:3:a:- \
                   TSB_ON_div_by_7674201.0_multi_by_3382258.0_reverse.wig:notex:4:a:- \
                   TSB_t0_div_by_4053535.0_multi_by_3382258.0_reverse.wig:notex:5:a:- \
                   TSB_t1_div_by_4448582.0_multi_by_3382258.0_reverse.wig:notex:6:a:- \
                   TSB_t2_div_by_6988975.0_multi_by_3382258.0_reverse.wig:notex:7:a:- \
                   pMEM_OD_0.2_div_by_4664201.0_multi_by_3382258.0_reverse.wig:notex:8:a:- \
                   pMEM_OD_0.5_div_by_5217691.0_multi_by_3382258.0_reverse.wig:notex:9:a:- \
                   pMEM_OD_1_div_by_4699153.0_multi_by_3382258.0_reverse.wig:notex:10:a:- \
                   pMEM_ON_div_by_3977549.0_multi_by_3382258.0_reverse.wig:notex:11:a:- \
                   pMEM_t0_div_by_3606941.0_multi_by_3382258.0_reverse.wig:notex:12:a:- \
                   pMEM_t1_div_by_3722754.0_multi_by_3382258.0_reverse.wig:notex:13:a:- \
                   pMEM_t2_div_by_3382258.0_multi_by_3382258.0_reverse.wig:notex:14:a:- \
                   TSB_OD_0.2_div_by_5033297.0_multi_by_3382258.0_forward.wig:notex:1:a:+ \
                   TSB_OD_0.5_div_by_4665048.0_multi_by_3382258.0_forward.wig:notex:2:a:+ \
                   TSB_OD_1_div_by_4443620.0_multi_by_3382258.0_forward.wig:notex:3:a:+ \
                   TSB_ON_div_by_7674201.0_multi_by_3382258.0_forward.wig:notex:4:a:+ \
                   TSB_t0_div_by_4053535.0_multi_by_3382258.0_forward.wig:notex:5:a:+ \
                   TSB_t1_div_by_4448582.0_multi_by_3382258.0_forward.wig:notex:6:a:+ \
                   TSB_t2_div_by_6988975.0_multi_by_3382258.0_forward.wig:notex:7:a:+ \
                   pMEM_OD_0.2_div_by_4664201.0_multi_by_3382258.0_forward.wig:notex:8:a:+ \
                   pMEM_OD_0.5_div_by_5217691.0_multi_by_3382258.0_forward.wig:notex:9:a:+ \
                   pMEM_OD_1_div_by_4699153.0_multi_by_3382258.0_forward.wig:notex:10:a:+ \
                   pMEM_ON_div_by_3977549.0_multi_by_3382258.0_forward.wig:notex:11:a:+ \
                   pMEM_t0_div_by_3606941.0_multi_by_3382258.0_forward.wig:notex:12:a:+ \
                   pMEM_t1_div_by_3722754.0_multi_by_3382258.0_forward.wig:notex:13:a:+ \
                   pMEM_t2_div_by_3382258.0_multi_by_3382258.0_forward.wig:notex:14:a:+"
    frag_libs="ID-001873-Staph_aureus_sample_mix_div_by_3162486.0_multi_by_3162486.0_forward.wig:frag:1:a:+ \
               ID-001873-Staph_aureus_sample_mix_div_by_3162486.0_multi_by_3162486.0_reverse.wig:frag:1:a:-"
    

#    COMBINE_GFF=/home/silas/Projects/2014-06-13-Sung-Huan-Yu-transcript_annotation_pipeline/Transap/output/target/annotation/all.gff
#    COMBINE_FASTA=/home/silas/Projects/2014-06-13-Sung-Huan-Yu-transcript_annotation_pipeline/Transap/output/target/fasta/all.fasta

#    set_up_analysis_folder
#    get_input_files    
#    get_target_fasta
#    annotation_transfer
    SNP_calling_reference
#    TSS_prediction
#    color_png_TSS
#    Transcriptome_assembly
#    Terminator_prediction
#    processing_site_prediction
#    color_png_processing_site
#    utr_detection
#    sRNA_detection
#    promoter_detection
#    CircRNA_detection
#    Go_term
#    sRNA_target
#    operon_detection
#    SNP_calling_target
#    PPI_network
#    Subcellular_localization
#    riboswitch
#    Optimize_TSSpredator
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
    if ! [ -d $TRANSAP_FOLDER ]
    then
        $PYTHON_PATH $TRANSAP_PATH create $TRANSAP_FOLDER
    fi
}

get_input_files(){
    $PYTHON_PATH $TRANSAP_PATH \
	get_input_files \
	-F $FTP_SOURCE \
	-g \
	-f \
	-e \
	-k \
	$TRANSAP_FOLDER
}

get_target_fasta(){
    $PYTHON_PATH $TRANSAP_PATH \
        get_target_fasta \
        -r $TRANSAP_FOLDER/input/reference/fasta \
	-o Staphylococcus_aureus_HG003:Staphylococcus_aureus_HG003 \
        -m $TRANSAP_FOLDER/input/mutation_table/Combined_table_Berscheid_and_Schuster.csv \
        $TRANSAP_FOLDER
}

annotation_transfer(){
    # instead of using "source"
#    . $PAGIT_HOME/sourceme.pagit
#    if [ -f $COMBINE_GFF ]
#    then
#        rm $COMBINE_GFF
#    fi
# --PAGIT_path /home/silas/Transap/tools/PAGIT \
    $PYTHON_PATH $TRANSAP_PATH \
        annotation_transfer \
        -tr Staphylococcus_aureus_HG003:NC_007795.1 \
        -re $TRANSAP_FOLDER/input/reference/annotation \
        -rf $TRANSAP_FOLDER/input/reference/fasta \
        -tf $TRANSAP_FOLDER/output/target/fasta \
        -e chromosome \
	-t Strain \
	-g \
	--RATT_path /home/silas/Transap_bak/tools/ratt-code \
        $TRANSAP_FOLDER
}

SNP_calling_reference(){
    $PYTHON_PATH $TRANSAP_PATH \
         snp \
        -p 1 2 3 \
        -t reference \
        -nw $TRANSAP_FOLDER/input/BAMs/BAMs_map_reference/tex_notex \
        -f $TRANSAP_FOLDER/input/reference/fasta \
        --samtools_path /home/silas/Transap/tools/samtools-bcftools-htslib-1.0_x64-linux/bin/samtools \
        --bcftools_path /home/silas/Transap/tools/samtools-bcftools-htslib-1.0_x64-linux/bin/bcftools \
        $TRANSAP_FOLDER
}

TSS_prediction(){
    $PYTHON_PATH $TRANSAP_PATH \
        tsspredator \
        --TSSpredator_path /home/silas/Transap/tools/TSSpredator_v1-04/TSSpredator.jar \
        -w $TRANSAP_FOLDER/input/wigs/tex_notex \
        -f $TRANSAP_FOLDER/output/target/fasta \
        -g $TRANSAP_FOLDER/output/target/annotation \
        -l $tex_notex_libs \
        -p TSB_OD_0.2 \
	   TSB_OD_0.5 \
	   TSB_OD_1 \
           TSB_ON \
	   TSB_t0 \
	   TSB_t1 \
	   TSB_t2 \
           pMEM_OD_0.2 \
	   pMEM_OD_0.5 \
	   pMEM_OD_1 \
	   pMEM_ON \
	   pMEM_t0 \
	   pMEM_t1 \
	   pMEM_t2 \
        -he 1.4 \
	-rh 1.3 \
	-fa 5.2 \
	-rf 1.2 \
	-bh 0.01 \
	-rm 1 \
	-m $TRANSAP_FOLDER/input/manual_TSS/Staphylococcus_aureus_HG003_manual_TSS.gff \
        -s \
        -v \
        $TRANSAP_FOLDER
}
#-ta $TRANSAP_FOLDER/output/transcriptome_assembly/gffs \

color_png_TSS(){
    $PYTHON_PATH $TRANSAP_PATH \
        color_png \
	-s TSS \
	-t 28 \
        $TRANSAP_FOLDER
}

Transcriptome_assembly(){
    $PYTHON_PATH $TRANSAP_PATH \
        Transcript_Assembly \
        -nw $TRANSAP_FOLDER/input/wigs/tex_notex \
	-fw $TRANSAP_FOLDER/input/wigs/fragment \
        -tl $tex_notex_libs \
        -fl $frag_libs \
        -r 1 \
        -te 2 \
	-ct $TRANSAP_FOLDER/output/TSS/gffs \
        -cg $TRANSAP_FOLDER/output/target/annotation \
        -g $TRANSAP_FOLDER/output/target/annotation \
        $TRANSAP_FOLDER
}
#        -fw $TRANSAP_FOLDER/input/wigs/fragment \
#        -g $TRANSAP_FOLDER/output/target/annotation \
#        -ct $TRANSAP_FOLDER/output/TSS/gffs \
#        -cg $TRANSAP_FOLDER/output/target/annotation \

Terminator_prediction(){
    $PYTHON_PATH $TRANSAP_PATH \
        Terminator \
        --TransTermHP_folder /home/silas/Transap/tools/transterm_hp_v2.09 \
        -f $TRANSAP_FOLDER/output/target/fasta \
        -g $TRANSAP_FOLDER/output/target/annotation \
        -s \
        -fw $TRANSAP_FOLDER/input/wigs/fragment \
        -tw $TRANSAP_FOLDER/input/wigs/tex_notex \
        -a $TRANSAP_FOLDER/output/transcriptome_assembly/gffs \
        -tl $tex_notex_libs \
        -fl $frag_libs \
        -te 2 \
        -r 1 \
        -tb \
	-sr srna_test \
        $TRANSAP_FOLDER
}

processing_site_prediction()
{
    $PYTHON_PATH $TRANSAP_PATH \
        tsspredator \
        --TSSpredator_path /home/silas/Transap/tools/TSSpredator_v1-04/TSSpredator.jar \
        -w $TRANSAP_FOLDER/input/wigs/tex_notex \
        -f $TRANSAP_FOLDER/output/target/fasta \
        -g $TRANSAP_FOLDER/output/target/annotation \
        -l $tex_notex_libs \
        -p TSB_OD_0.2 \
           TSB_OD_0.5 \
           TSB_OD_1 \
           TSB_ON \
           TSB_t0 \
           TSB_t1 \
           TSB_t2 \
           pMEM_OD_0.2 \
           pMEM_OD_0.5 \
           pMEM_OD_1 \
           pMEM_ON \
           pMEM_t0 \
           pMEM_t1 \
           pMEM_t2 \
        -he 0.5 \
        -rh 0.4 \
        -fa 2.1 \
        -rf 0.1 \
	-bh 0.0 \
	-rm 1 \
	-t processing_site \
	-s \
	-m $TRANSAP_FOLDER/input/manual_processing_site/Staphylococcus_aureus_HG003_manual_processing_105000.gff \
        -le 105000 \
        $TRANSAP_FOLDER
}

color_png_processing_site(){
    $PYTHON_PATH $TRANSAP_PATH \
        color_png \
        -s processing_site \
        -t 28 \
        $TRANSAP_FOLDER
}

utr_detection(){
    $PYTHON_PATH $TRANSAP_PATH \
        UTR_detection \
        -g $TRANSAP_FOLDER/output/target/annotation \
	-t $TRANSAP_FOLDER/output/TSS/gffs \
        -a $TRANSAP_FOLDER/output/transcriptome_assembly/gffs \
	-e $TRANSAP_FOLDER/output/terminator/gffs \
        $TRANSAP_FOLDER
}

sRNA_detection(){
    $PYTHON_PATH $TRANSAP_PATH \
        sRNA_detection \
        -d 1 2 3 4 \
        -g $TRANSAP_FOLDER/output/target/annotation \
        -t $TRANSAP_FOLDER/output/TSS/gffs \
        -p $TRANSAP_FOLDER/output/processing_site/gffs \
        -a $TRANSAP_FOLDER/output/transcriptome_assembly/gffs \
        -fw $TRANSAP_FOLDER/input/wigs/fragment \
        -tw $TRANSAP_FOLDER/input/wigs/tex_notex \
        -f $TRANSAP_FOLDER/output/target/fasta \
        -m \
        -u \
        -sd $TRANSAP_FOLDER/input/database/sRNA_database \
        -nd $TRANSAP_FOLDER/input/database/nr \
        -tl $tex_notex_libs \
        -fl $frag_libs \
        -te 2 \
        -r 1 \
        -ba \
        $TRANSAP_FOLDER
}

promoter_detection(){
    $PYTHON_PATH $TRANSAP_PATH \
        Promoter \
        -t $TRANSAP_FOLDER/output/TSS/gffs \
        -f $TRANSAP_FOLDER/output/target/fasta \
	-w 50 51 45 2-10 \
	-p 10 \
	-g $TRANSAP_FOLDER/output/TSS/gffs \
        $TRANSAP_FOLDER
}

CircRNA_detection(){
    $PYTHON_PATH $TRANSAP_PATH \
        CircRNA \
	-a \
        -f $TRANSAP_FOLDER/output/target/fasta \
        -p 10 \
	-g $TRANSAP_FOLDER/output/target/annotation \
	-cg \
	$TRANSAP_FOLDER	
}

Go_term(){
    $PYTHON_PATH $TRANSAP_PATH \
        Go_term \
        -g $TRANSAP_FOLDER/output/target/annotation\
        $TRANSAP_FOLDER
}

sRNA_target(){
    $PYTHON_PATH $TRANSAP_PATH \
         sRNA_target \
        -g $TRANSAP_FOLDER/output/target/annotation \
        -f $TRANSAP_FOLDER/output/target/fasta \
        -r $TRANSAP_FOLDER/output/sRNA/gffs/best \
        -p both \
        $TRANSAP_FOLDER
}

operon_detection(){
    $PYTHON_PATH $TRANSAP_PATH \
         Operon \
        -g $TRANSAP_FOLDER/output/target/annotation \
        -t $TRANSAP_FOLDER/output/TSS/gffs \
        -a $TRANSAP_FOLDER/output/transcriptome_assembly/gffs \
	-u5 $TRANSAP_FOLDER/output/utr/5UTR/gffs \
        -u3 $TRANSAP_FOLDER/output/utr/3UTR/gffs \
	-e $TRANSAP_FOLDER/output/terminator/gffs \
	-s \
	-c \
        $TRANSAP_FOLDER
}

SNP_calling_target(){
    $PYTHON_PATH $TRANSAP_PATH \
         SNP_calling \
	-t target \
	-p 1 2 3 \
	-nw $TRANSAP_FOLDER/input/BAMs/BAMs_map_target/tex_notex \
	-fw $TRANSAP_FOLDER/input/BAMs/BAMs_map_target/fragment \
	-f $TRANSAP_FOLDER/output/target/fasta \
        $TRANSAP_FOLDER
}

PPI_network(){
    $PYTHON_PATH $TRANSAP_PATH \
         PPI_network \
        -s all:Staphylococcus_aureus_HG003.ptt:Staphylococcus_aureus_HG003:'Staphylococcus aureus 8325':'Staphylococcus aureus' \
        -p $TRANSAP_FOLDER/output/target/annotation \
        -d $TRANSAP_FOLDER/input/database/species.v9.1.txt \
        -n \
        -ns 4000 \
        $TRANSAP_FOLDER
}

Subcellular_localization(){
    $PYTHON_PATH $TRANSAP_PATH \
         subcellular_localization \
        -g $TRANSAP_FOLDER/output/target/annotation \
        -f $TRANSAP_FOLDER/output/target/fasta \
        -m \
        -b positive \
        $TRANSAP_FOLDER
}

riboswitch(){
    $PYTHON_PATH $TRANSAP_PATH \
         riboswitch \
        -g $TRANSAP_FOLDER/output/target/annotation \
        -f $TRANSAP_FOLDER/output/target/fasta \
        -r \
        -i $TRANSAP_FOLDER/input/riboswitch_ID/Rfam_riboswitch_ID.csv \
        -R $TRANSAP_FOLDER/input/database/Rfam/CMs/Rfam.cm \
        $TRANSAP_FOLDER
}

Optimize_TSSpredator(){
    $PYTHON_PATH $TRANSAP_PATH \
         optimize_TSSpredator \
        --TSSpredator_path /home/silas/Transap/tools/TSSpredator_v1-04/TSSpredator.jar \
        -w $TRANSAP_FOLDER/input/wigs/tex_notex \
        -fs $TRANSAP_FOLDER/output/target/fasta/Staphylococcus_aureus_HG003.fa \
        -g $TRANSAP_FOLDER/output/target/annotation/Staphylococcus_aureus_HG003.gff \
        -n Staphylococcus_aureus_HG003 \
	-l $tex_notex_libs \
	-p TSB_OD_0.2 \
           TSB_OD_0.5 \
           TSB_OD_1 \
           TSB_ON \
           TSB_t0 \
           TSB_t1 \
           TSB_t2 \
           pMEM_OD_0.2 \
           pMEM_OD_0.5 \
           pMEM_OD_1 \
           pMEM_ON \
           pMEM_t0 \
           pMEM_t1 \
           pMEM_t2 \
        -m $TRANSAP_FOLDER/input/manual_TSS/Staphylococcus_aureus_HG003_manual_TSS.gff  \
        -c 6 \
	-rm 1 \
        $TRANSAP_FOLDER
}

main
