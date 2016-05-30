main(){
    PATH_FILE=$(pwd)
    PYTHON_PATH=python3
    STRAINS=Staphylococcus_aureus_HG003
    ANNOGESIC_PATH=/home/silas/ANNOgesic/bin/annogesic
    ANNOGESIC_FOLDER=ANNOgesic
#    FTP_SOURCE=ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Staphylococcus_aureus_NCTC_8325_uid57795
    FTP_SOURCE=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000013425.1_ASM1342v1/
    ID_MAPPING_SOURCE=ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
    PAGIT_HOME=/home/silas/ANNOgesic/tools/PAGIT
    tex_notex_libs="pMEM_OD_0.2_div_by_4664494.0_multi_by_3162518.0_forward.wig:notex:1:a:+,\
pMEM_OD_0.2_div_by_4664494.0_multi_by_3162518.0_reverse.wig:notex:1:a:-,\
pMEM_OD_0.2_TEX_div_by_4582476.0_multi_by_3162518.0_forward.wig:tex:1:a:+,\
pMEM_OD_0.2_TEX_div_by_4582476.0_multi_by_3162518.0_reverse.wig:tex:1:a:-,\
pMEM_OD_0.5_div_by_5217975.0_multi_by_3162518.0_forward.wig:notex:2:a:+,\
pMEM_OD_0.5_div_by_5217975.0_multi_by_3162518.0_reverse.wig:notex:2:a:-,\
pMEM_OD_0.5_TEX_div_by_4968963.0_multi_by_3162518.0_forward.wig:tex:2:a:+,\
pMEM_OD_0.5_TEX_div_by_4968963.0_multi_by_3162518.0_reverse.wig:tex:2:a:-,\
pMEM_OD_1_div_by_4699406.0_multi_by_3162518.0_forward.wig:notex:3:a:+,\
pMEM_OD_1_div_by_4699406.0_multi_by_3162518.0_reverse.wig:notex:3:a:-,\
pMEM_OD_1_TEX_div_by_4720486.0_multi_by_3162518.0_forward.wig:tex:3:a:+,\
pMEM_OD_1_TEX_div_by_4720486.0_multi_by_3162518.0_reverse.wig:tex:3:a:-,\
pMEM_ON_div_by_3977820.0_multi_by_3162518.0_forward.wig:notex:4:a:+,\
pMEM_ON_div_by_3977820.0_multi_by_3162518.0_reverse.wig:notex:4:a:-,\
pMEM_ON_TEX_div_by_5455288.0_multi_by_3162518.0_forward.wig:tex:4:a:+,\
pMEM_ON_TEX_div_by_5455288.0_multi_by_3162518.0_reverse.wig:tex:4:a:-,\
pMEM_t0_div_by_3607269.0_multi_by_3162518.0_forward.wig:notex:5:a:+,\
pMEM_t0_div_by_3607269.0_multi_by_3162518.0_reverse.wig:notex:5:a:-,\
pMEM_t0_TEX_div_by_4873740.0_multi_by_3162518.0_forward.wig:tex:5:a:+,\
pMEM_t0_TEX_div_by_4873740.0_multi_by_3162518.0_reverse.wig:tex:5:a:-,\
pMEM_t1_div_by_3723072.0_multi_by_3162518.0_forward.wig:notex:6:a:+,\
pMEM_t1_div_by_3723072.0_multi_by_3162518.0_reverse.wig:notex:6:a:-,\
pMEM_t1_TEX_div_by_4637645.0_multi_by_3162518.0_forward.wig:tex:6:a:+,\
pMEM_t1_TEX_div_by_4637645.0_multi_by_3162518.0_reverse.wig:tex:6:a:-,\
pMEM_t2_div_by_3382513.0_multi_by_3162518.0_forward.wig:notex:7:a:+,\
pMEM_t2_div_by_3382513.0_multi_by_3162518.0_reverse.wig:notex:7:a:-,\
pMEM_t2_TEX_div_by_4783039.0_multi_by_3162518.0_forward.wig:tex:7:a:+,\
pMEM_t2_TEX_div_by_4783039.0_multi_by_3162518.0_reverse.wig:tex:7:a:-,\
TSB_OD_0.2_div_by_5033683.0_multi_by_3162518.0_forward.wig:notex:8:a:+,\
TSB_OD_0.2_div_by_5033683.0_multi_by_3162518.0_reverse.wig:notex:8:a:-,\
TSB_OD_0.2_TEX_div_by_3598779.0_multi_by_3162518.0_forward.wig:tex:8:a:+,\
TSB_OD_0.2_TEX_div_by_3598779.0_multi_by_3162518.0_reverse.wig:tex:8:a:-,\
TSB_OD_0.5_div_by_4665441.0_multi_by_3162518.0_forward.wig:notex:9:a:+,\
TSB_OD_0.5_div_by_4665441.0_multi_by_3162518.0_reverse.wig:notex:9:a:-,\
TSB_OD_0.5_TEX_div_by_4420785.0_multi_by_3162518.0_forward.wig:tex:9:a:+,\
TSB_OD_0.5_TEX_div_by_4420785.0_multi_by_3162518.0_reverse.wig:tex:9:a:-,\
TSB_OD_1_div_by_4444037.0_multi_by_3162518.0_forward.wig:notex:10:a:+,\
TSB_OD_1_div_by_4444037.0_multi_by_3162518.0_reverse.wig:notex:10:a:-,\
TSB_OD_1_TEX_div_by_3956305.0_multi_by_3162518.0_forward.wig:tex:10:a:+,\
TSB_OD_1_TEX_div_by_3956305.0_multi_by_3162518.0_reverse.wig:tex:10:a:-,\
TSB_ON_div_by_7674925.0_multi_by_3162518.0_forward.wig:notex:11:a:+,\
TSB_ON_div_by_7674925.0_multi_by_3162518.0_reverse.wig:notex:11:a:-,\
TSB_ON_TEX_div_by_4837094.0_multi_by_3162518.0_forward.wig:tex:11:a:+,\
TSB_ON_TEX_div_by_4837094.0_multi_by_3162518.0_reverse.wig:tex:11:a:-,\
TSB_t0_div_by_4053752.0_multi_by_3162518.0_forward.wig:notex:12:a:+,\
TSB_t0_div_by_4053752.0_multi_by_3162518.0_reverse.wig:notex:12:a:-,\
TSB_t0_TEX_div_by_4817822.0_multi_by_3162518.0_forward.wig:tex:12:a:+,\
TSB_t0_TEX_div_by_4817822.0_multi_by_3162518.0_reverse.wig:tex:12:a:-,\
TSB_t1_div_by_4448852.0_multi_by_3162518.0_forward.wig:notex:13:a:+,\
TSB_t1_div_by_4448852.0_multi_by_3162518.0_reverse.wig:notex:13:a:-,\
TSB_t1_TEX_div_by_4958138.0_multi_by_3162518.0_forward.wig:tex:13:a:+,\
TSB_t1_TEX_div_by_4958138.0_multi_by_3162518.0_reverse.wig:tex:13:a:-,\
TSB_t2_div_by_6989369.0_multi_by_3162518.0_forward.wig:notex:14:a:+,\
TSB_t2_div_by_6989369.0_multi_by_3162518.0_reverse.wig:notex:14:a:-,\
TSB_t2_TEX_div_by_5024257.0_multi_by_3162518.0_forward.wig:tex:14:a:+,\
TSB_t2_TEX_div_by_5024257.0_multi_by_3162518.0_reverse.wig:tex:14:a:-"
    frag_libs="ID-001873-Staph_aureus_sample_mix_div_by_3162518.0_multi_by_3162518.0_forward.wig:frag:1:a:+,\
ID-001873-Staph_aureus_sample_mix_div_by_3162518.0_multi_by_3162518.0_reverse.wig:frag:1:a:-"

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
    sRNA_detection
#    sORF_detection
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
#    gen_screenshot
#    color_png
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
        $PYTHON_PATH $ANNOGESIC_PATH create $ANNOGESIC_FOLDER
    fi
}

get_input_files(){
    $PYTHON_PATH $ANNOGESIC_PATH \
	get_input_files \
	-F $FTP_SOURCE \
	-g \
	-f \
	-e \
	-k \
	$ANNOGESIC_FOLDER
}

get_target_fasta(){
    $PYTHON_PATH $ANNOGESIC_PATH \
        get_target_fasta \
        -r $ANNOGESIC_FOLDER/input/reference/fasta \
	-o Staphylococcus_aureus_HG003:Staphylococcus_aureus_HG003 \
        -m $ANNOGESIC_FOLDER/input/mutation_table/Combined_table_Berscheid_and_Schuster.csv \
        $ANNOGESIC_FOLDER
}

annotation_transfer(){
    # instead of using "source"
    . $PAGIT_HOME/sourceme.pagit
    $PYTHON_PATH $ANNOGESIC_PATH \
        annotation_transfer \
        -re $ANNOGESIC_FOLDER/input/reference/annotation \
        -rf $ANNOGESIC_FOLDER/input/reference/fasta \
        -tf $ANNOGESIC_FOLDER/output/target/fasta \
        -e chromosome \
	-t Strain \
	-p NC_007795.1:Staphylococcus_aureus_HG003 \
	-g \
	--RATT_path /home/silas/ANNOgesic/tools/PAGIT/RATT/start.ratt.sh \
	$ANNOGESIC_FOLDER
}

expression_analysis(){
    $PYTHON_PATH $ANNOGESIC_PATH \
         expression_analysis \
        -g $ANNOGESIC_FOLDER/output/target/annotation \
        -tl $tex_notex_libs \
        -tw $ANNOGESIC_FOLDER/input/wigs/tex_notex \
        -f CDS,tRNA,rRNA \
        -rt 1 \
	-ct high \
        $ANNOGESIC_FOLDER
}

SNP_calling_reference(){
    $PYTHON_PATH $ANNOGESIC_PATH \
         snp \
        -p 1,2,3 \
        -t reference \
        -tw $ANNOGESIC_FOLDER/input/BAMs/BAMs_map_reference/tex_notex \
        -f $ANNOGESIC_FOLDER/input/reference/fasta \
        --samtools_path /home/silas/ANNOgesic/tools/samtools-bcftools-htslib-1.0_x64-linux/bin/samtools \
        --bcftools_path /home/silas/ANNOgesic/tools/samtools-bcftools-htslib-1.0_x64-linux/bin/bcftools \
        $ANNOGESIC_FOLDER
}

TSS_prediction(){
    $PYTHON_PATH $ANNOGESIC_PATH \
        tsspredator \
        --TSSpredator_path /home/silas/ANNOgesic/tools/TSSpredator_v1-04/TSSpredator.jar \
        -w $ANNOGESIC_FOLDER/input/wigs/tex_notex \
        -f $ANNOGESIC_FOLDER/output/target/fasta \
        -g $ANNOGESIC_FOLDER/output/target/annotation \
        -l $tex_notex_libs \
        -p TSB_OD_0.2,TSB_OD_0.5,TSB_OD_1,TSB_ON,TSB_t0,TSB_t1,TSB_t2,pMEM_OD_0.2,pMEM_OD_0.5,pMEM_OD_1,pMEM_ON,pMEM_t0,pMEM_t1,pMEM_t2 \
	-he 1.4 \
        -rh 1.3 \
        -fa 5.1 \
        -rf 0.7 \
        -bh 0.055 \
        -ef 0.4 \
        -pf 1.6 \
	-m $ANNOGESIC_FOLDER/input/manual_TSS/Staphylococcus_aureus_HG003_manual_TSS.gff \
	-rm 1 \
        -s \
        -v \
        $ANNOGESIC_FOLDER
}

#-he 1.4 \
#        -rh 1.3 \
#        -fa 5.1 \
#        -rf 0.7 \
#        -bh 0.055 \
#        -ef 0.4 \
#        -pf 1.6 \
#-ta $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs \
#-m $ANNOGESIC_FOLDER/input/manual_TSS/Staphylococcus_aureus_HG003_manual_TSS.gff \
#        -v \

Transcriptome_assembly(){
    $PYTHON_PATH $ANNOGESIC_PATH \
        transcript_assembly \
        -tw $ANNOGESIC_FOLDER/input/wigs/tex_notex \
	-fw $ANNOGESIC_FOLDER/input/wigs/fragment \
        -tl $tex_notex_libs \
        -fl $frag_libs \
        -rt 1 \
	-rf 1 \
	-cf gene,CDS \
	-ct $ANNOGESIC_FOLDER/output/TSS/gffs \
        -cg $ANNOGESIC_FOLDER/output/target/annotation \
        -g $ANNOGESIC_FOLDER/output/target/annotation \
	-tr $ANNOGESIC_FOLDER/output/terminator/gffs/best \
        $ANNOGESIC_FOLDER
}

Terminator_prediction(){
    $PYTHON_PATH $ANNOGESIC_PATH \
        terminator \
        --TransTermHP_path /home/silas/ANNOgesic/tools/transterm_hp_v2.09/transterm \
        --expterm_path /home/silas/ANNOgesic/tools/transterm_hp_v2.09/expterm.dat \
        -f $ANNOGESIC_FOLDER/output/target/fasta \
        -g $ANNOGESIC_FOLDER/output/target/annotation \
        -fw $ANNOGESIC_FOLDER/input/wigs/fragment \
        -tw $ANNOGESIC_FOLDER/input/wigs/tex_notex \
        -a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs \
        -tl $tex_notex_libs \
        -fl $frag_libs \
        -rt 1 \
	-rf 1 \
        -tb \
	-s \
        $ANNOGESIC_FOLDER
}

processing_site_prediction()
{
    $PYTHON_PATH $ANNOGESIC_PATH \
        tsspredator \
        --TSSpredator_path /home/silas/ANNOgesic/tools/TSSpredator_v1-04/TSSpredator.jar \
        -w $ANNOGESIC_FOLDER/input/wigs/tex_notex \
        -f $ANNOGESIC_FOLDER/output/target/fasta \
        -g $ANNOGESIC_FOLDER/output/target/annotation \
        -l $tex_notex_libs \
        -p TSB_OD_0.2,TSB_OD_0.5,TSB_OD_1,TSB_ON,TSB_t0,TSB_t1,TSB_t2,pMEM_OD_0.2,pMEM_OD_0.5,pMEM_OD_1,pMEM_ON,pMEM_t0,pMEM_t1,pMEM_t2 \
        -he 1.0 \
        -rh 0.1 \
        -fa 2.3 \
        -rf 0.9 \
	-bh 0.0 \
	-ef 4.7 \
	-pf 0.7 \
	-rm 1 \
	-t processing_site \
	-s \
        -le 200000 \
        -m $ANNOGESIC_FOLDER/input/manual_processing_site/Staphylococcus_aureus_HG003_manual_processing_200000.gff \
        $ANNOGESIC_FOLDER
}

utr_detection(){
    $PYTHON_PATH $ANNOGESIC_PATH \
        utr \
        -g $ANNOGESIC_FOLDER/output/target/annotation \
	-t $ANNOGESIC_FOLDER/output/TSS/gffs \
        -a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs \
	-e $ANNOGESIC_FOLDER/output/terminator/gffs/best \
        $ANNOGESIC_FOLDER
}

sRNA_detection(){
    $PYTHON_PATH $ANNOGESIC_PATH \
        srna \
        -d tss,promoter,blast_nr,blast_srna,sorf,term,sec_str \
        --Vienna_folder /home/silas/ANNOgesic/tools/ViennaRNA-2.1.7/Progs \
        --Vienna_utils /home/silas/ANNOgesic/tools/ViennaRNA-2.1.7/Utils \
        --blast_plus_folder /home/silas/ANNOgesic/tools \
        -g $ANNOGESIC_FOLDER/output/target/annotation \
        -t $ANNOGESIC_FOLDER/output/TSS/gffs \
        -p $ANNOGESIC_FOLDER/output/processing_site/gffs \
        -a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs \
        -tw $ANNOGESIC_FOLDER/input/wigs/tex_notex \
        -f $ANNOGESIC_FOLDER/output/target/fasta \
        -O $ANNOGESIC_FOLDER/output/sORF/gffs/best \
        -tf $ANNOGESIC_FOLDER/output/terminator/gffs/best \
        -pt $ANNOGESIC_FOLDER/output/promoter_analysis/Staphylococcus_aureus_HG003/promoter_motifs_Staphylococcus_aureus_HG003_allstrain_all_types_45_nt/meme.csv \
        -pn MOTIF_1 \
        -m \
        -u \
        -sd $ANNOGESIC_FOLDER/input/database/sRNA_database \
        -nd $ANNOGESIC_FOLDER/input/database/nr \
        -tl $tex_notex_libs \
        -rt 1 \
        -ba \
        -fw $ANNOGESIC_FOLDER/input/wigs/fragment \
        -fl $frag_libs \
        -rf 1 \
        $ANNOGESIC_FOLDER
}

sORF_detection(){
    $PYTHON_PATH $ANNOGESIC_PATH \
        sorf \
        -g $ANNOGESIC_FOLDER/output/target/annotation \
        -t $ANNOGESIC_FOLDER/output/TSS/gffs \
        -a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs \
        -fw $ANNOGESIC_FOLDER/input/wigs/fragment \
        -tw $ANNOGESIC_FOLDER/input/wigs/tex_notex \
        -f $ANNOGESIC_FOLDER/output/target/fasta \
	-s $ANNOGESIC_FOLDER/output/sRNA/gffs/best \
        -tl $tex_notex_libs \
        -fl $frag_libs \
        -rt 1 \
	-rf 1 \
        -u \
        $ANNOGESIC_FOLDER
}
#-s $ANNOGESIC_FOLDER/output/sRNA/gffs/best \
promoter_detection(){
    $PYTHON_PATH $ANNOGESIC_PATH \
        promoter \
        -t $ANNOGESIC_FOLDER/output/TSS/gffs \
        -f $ANNOGESIC_FOLDER/output/target/fasta \
	-w 50,51,45,2-10 \
	-c \
	-tl $tex_notex_libs \
	-tw $ANNOGESIC_FOLDER/input/wigs/tex_notex \
        $ANNOGESIC_FOLDER
}

CircRNA_detection(){
    $PYTHON_PATH $ANNOGESIC_PATH \
        circrna \
        -f $ANNOGESIC_FOLDER/output/target/fasta \
        -p 10 \
	-g $ANNOGESIC_FOLDER/output/target/annotation \
	-cg \
	-tb ANNOgesic/output/circRNA/segemehl_align/Staphylococcus_aureus_HG003 \
	--samtools_path /home/silas/ANNOgesic/tools/samtools-bcftools-htslib-1.0_x64-linux/bin/samtools \
        --segemehl_folder /home/silas/ANNOgesic/tools/segemehl_0_1_9/segemehl \
	$ANNOGESIC_FOLDER	
}

Go_term(){
    $PYTHON_PATH $ANNOGESIC_PATH \
        go_term \
        -g $ANNOGESIC_FOLDER/output/target/annotation\
        -a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs \
        $ANNOGESIC_FOLDER
}

sRNA_target(){
    $PYTHON_PATH $ANNOGESIC_PATH \
         srna_target \
        --Vienna_folder /home/silas/ANNOgesic/tools/ViennaRNA-2.1.7/Progs \
        -g $ANNOGESIC_FOLDER/output/target/annotation \
        -f $ANNOGESIC_FOLDER/output/target/fasta \
        -r $ANNOGESIC_FOLDER/output/sRNA/gffs/best \
        -q all \
        -p both \
        $ANNOGESIC_FOLDER
}

operon_detection(){
    $PYTHON_PATH $ANNOGESIC_PATH \
         operon \
        -g $ANNOGESIC_FOLDER/output/target/annotation \
        -t $ANNOGESIC_FOLDER/output/TSS/gffs \
        -a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs \
	-u5 $ANNOGESIC_FOLDER/output/UTR/5UTR/gffs \
        -u3 $ANNOGESIC_FOLDER/output/UTR/3UTR/gffs \
	-e $ANNOGESIC_FOLDER/output/terminator/gffs/best \
	-s \
	-c \
        $ANNOGESIC_FOLDER
}

SNP_calling_target(){
    $PYTHON_PATH $ANNOGESIC_PATH \
         snp \
	-t target \
	-p 1,2,3 \
	-tw $ANNOGESIC_FOLDER/input/BAMs/BAMs_map_target/tex_notex \
	-f $ANNOGESIC_FOLDER/output/target/fasta \
        --samtools_path /home/silas/ANNOgesic/tools/samtools-bcftools-htslib-1.0_x64-linux/bin/samtools \
        --bcftools_path /home/silas/ANNOgesic/tools/samtools-bcftools-htslib-1.0_x64-linux/bin/bcftools \
        $ANNOGESIC_FOLDER
}

PPI_network(){
    $PYTHON_PATH $ANNOGESIC_PATH \
         ppi_network \
        -s Staphylococcus_aureus_HG003.gff:Staphylococcus_aureus_HG003:'Staphylococcus aureus 8325':'Staphylococcus aureus' \
        -g $ANNOGESIC_FOLDER/output/target/annotation \
        -d $ANNOGESIC_FOLDER/input/database/species.v9.1.txt \
        -n \
	-q all \
        -ns 4000 \
        $ANNOGESIC_FOLDER
}

Subcellular_localization(){
    $PYTHON_PATH $ANNOGESIC_PATH \
         subcellular_localization \
        -g $ANNOGESIC_FOLDER/output/target/annotation \
        -a $ANNOGESIC_FOLDER/output/transcriptome_assembly/gffs \
        -f $ANNOGESIC_FOLDER/output/target/fasta \
        -m \
        -b positive \
        --Psortb_path /home/silas/ANNOgesic/tools/psortb/bin/psort \
        $ANNOGESIC_FOLDER
}

riboswitch(){
    $PYTHON_PATH $ANNOGESIC_PATH \
         riboswitch \
        -g $ANNOGESIC_FOLDER/output/target/annotation \
        -f $ANNOGESIC_FOLDER/output/target/fasta \
        -i $ANNOGESIC_FOLDER/input/riboswitch_ID/Rfam_riboswitch_ID.csv \
        -R $ANNOGESIC_FOLDER/input/database/Rfam/CMs/Rfam.cm \
	-t $ANNOGESIC_FOLDER/output/UTR/5UTR/gffs \
	--infernal_path /home/silas/ANNOgesic/tools/infernal-1.1.1/src \
        $ANNOGESIC_FOLDER
}

Optimize_TSSpredator(){
    $PYTHON_PATH $ANNOGESIC_PATH \
         optimize_tsspredator \
        --TSSpredator_path /home/silas/ANNOgesic/tools/TSSpredator_v1-04/TSSpredator.jar \
        -w $ANNOGESIC_FOLDER/input/wigs/tex_notex \
        -fs $ANNOGESIC_FOLDER/output/target/fasta \
        -g $ANNOGESIC_FOLDER/output/target/annotation \
        -n Staphylococcus_aureus_HG003 \
	-l $tex_notex_libs \
	-p TSB_OD_0.2,TSB_OD_0.5,TSB_OD_1,TSB_ON,TSB_t0,TSB_t1,TSB_t2,pMEM_OD_0.2,pMEM_OD_0.5,pMEM_OD_1,pMEM_ON,pMEM_t0,pMEM_t1,pMEM_t2 \
        -c 5 \
	-rm 1 \
        -m $ANNOGESIC_FOLDER/input/manual_TSS/Staphylococcus_aureus_HG003_manual_TSS.gff \
        -t TSS \
        $ANNOGESIC_FOLDER
}
#        -le 200000 \
#        -m $ANNOGESIC_FOLDER/input/manual_processing_site/Staphylococcus_aureus_HG003_manual_processing_200000.gff \
#        -t Processing_site \

gen_screenshot(){
    $PYTHON_PATH $ANNOGESIC_PATH \
         screenshot \
        -mg $ANNOGESIC_FOLDER/output/sRNA/gffs/best/Staphylococcus_aureus_HG003_sRNA.gff \
        -sg $ANNOGESIC_FOLDER/output/target/annotation/Staphylococcus_aureus_HG003.gff \
            $ANNOGESIC_FOLDER/output/TSS/gffs/Staphylococcus_aureus_HG003_TSS.gff \
            $ANNOGESIC_FOLDER/output/processing_site/gffs/Staphylococcus_aureus_HG003_processing.gff \
        -f $ANNOGESIC_FOLDER/output/target/fasta/Staphylococcus_aureus_HG003.fa \
        -o $ANNOGESIC_FOLDER/output/sRNA/screenshots \
        -fl $frag_libs \
        -tl $tex_notex_libs \
        -fw $ANNOGESIC_FOLDER/input/wigs/fragment \
        -tw $ANNOGESIC_FOLDER/input/wigs/tex_notex \
    $ANNOGESIC_FOLDER
}

color_png(){
    $PYTHON_PATH $ANNOGESIC_PATH \
        color_png \
        -t 29 \
	-f $ANNOGESIC_FOLDER/output/sRNA \
	--ImageMagick_covert_path /home/silas/ANNOgesic/tools/ImageMagick-6.9.0-0/utilities/convert \
        $ANNOGESIC_FOLDER
}

main
