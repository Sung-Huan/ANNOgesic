main(){
	ANNOGESIC_FOLDER="ANNOgesic"
	LIBS="ANNOgesic/input/wigs/fragment/S1_1_reverse.wig:frag:1:a:- \
	      ANNOgesic/input/wigs/fragment/S1_1_forward.wig:frag:1:a:+ \
	      ANNOgesic/input/wigs/fragment/S1_2_reverse.wig:frag:1:b:- \
	      ANNOgesic/input/wigs/fragment/S1_2_forward.wig:frag:1:b:+ \
	      ANNOgesic/input/wigs/fragment/S2_1_reverse.wig:frag:2:a:- \
	      ANNOgesic/input/wigs/fragment/S2_1_forward.wig:frag:2:a:+ \
	      ANNOgesic/input/wigs/fragment/S2_2_reverse.wig:frag:2:b:- \
	      ANNOgesic/input/wigs/fragment/S2_2_forward.wig:frag:2:b:+ \
	      ANNOgesic/input/wigs/fragment/S3_1_reverse.wig:frag:3:a:- \
	      ANNOgesic/input/wigs/fragment/S3_1_forward.wig:frag:3:a:+ \
	      ANNOgesic/input/wigs/fragment/S3_2_reverse.wig:frag:3:b:- \
              ANNOgesic/input/wigs/fragment/S3_2_forward.wig:frag:3:b:+ \
	      ANNOgesic/input/wigs/fragment/S4_1_reverse.wig:frag:4:a:- \
	      ANNOgesic/input/wigs/fragment/S4_1_forward.wig:frag:4:a:+ \
	      ANNOgesic/input/wigs/fragment/S4_2_reverse.wig:frag:4:b:- \
              ANNOgesic/input/wigs/fragment/S4_2_forward.wig:frag:4:b:+"
	#create_folder
	#transcript
	#terminator
	#sRNA
}


create_folder(){
	annogesic create -pj $ANNOGESIC_FOLDER
}

transcript(){
	annogesic transcript \
        -g $ANNOGESIC_FOLDER/input/references/annotations/NC002929_2.gff \
        -fl $LIBS \
        -rf all_1 \
	-cf gene CDS \
        -pj $ANNOGESIC_FOLDER
}

terminator(){
	annogesic terminator \
        -f $ANNOGESIC_FOLDER/input/references/fasta_files/NC002929_2.fna \
        -g $ANNOGESIC_FOLDER/input/references/annotations/NC002929_2.gff \
        -a $ANNOGESIC_FOLDER/output/transcripts/gffs/NC002929_2_transcript.gff \
        -fl $LIBS \
        -rf all_1 \
        -pj $ANNOGESIC_FOLDER
}

sRNA(){
#    #### If you have no sRNA database, you can use BSRD.
#    wget -cP ANNOgesic/input/databases/ https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/database/sRNA_database_BSRD.fa
#    #### If you have no nr databse, please download it.
#    wget -cP ANNOgesic/input/databases/ ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
#    gunzip ANNOgesic/input/databases/nr.gz
#    mv ANNOgesic/input/databases/nr ANNOgesic/input/databases/nr.fa
	annogesic srna \
        -d blast_srna sec_str blast_nr \
        -g $ANNOGESIC_FOLDER/input/references/annotations/NC002929_2.gff \
        -a $ANNOGESIC_FOLDER/output/transcripts/gffs/NC002929_2_transcript.gff \
        -f $ANNOGESIC_FOLDER/input/references/fasta_files/NC002929_2.fna \
        -e $ANNOGESIC_FOLDER/output/terminators/gffs/all_candidates/NC002929_2_term.gff \
        -m \
	-cs \
        -sf \
        -nf \
        -nd $ANNOGESIC_FOLDER/input/databases/nr \
	-sd $ANNOGESIC_FOLDER/input/databases/sRNA_database_BSRD \
	-fl $LIBS \
        -rf all_1 \
        -pj $ANNOGESIC_FOLDER
}

main
