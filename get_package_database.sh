main(){
    ##################################################################################################
    # This is the shell script for download package from Websit.                                     #
    # Please follow the commands of every function to modify .bashrc file if it is required.         #
    # If you want to use new version of package, please change the path in .bashrc.                  #
    # If you already have these package installed, you can skip it.                                  #
    ##################################################################################################

    PATH_FILE=$(pwd)
    TOOL_PATH=$PATH_FILE/tools  ### specify the path for installing tools
    DATABASE=$PATH_FILE/databases  ### specify the path for downloading databases
    echo "Please specify the absolute path for installing tools:"
    read TOOL_PATH
    if [[ -d "$TOOL_PATH" ]]
    then
        echo "$TOOL_PATH is found."
    else
        echo "$TOOL_PATH is not found."
        exit;
    fi
    echo "Please specify the absolute path for storing databases:"
    read DATABASE
     if [[ -d "$DATABASE" ]]
    then
        echo "$DATABASE is found."
    else
        echo "$DATABASE is not found."
        exit;
    fi
    echo "Start to install/download."

    RED='\033[1;31m'
    NC='\033[0m'
    declare -A PROGs=( ["TSSpredator"]=get_TSSpredator ["CRT"]=get_CRT ["matplotlib"]=get_matplotlib \
                       ["numpy"]=get_numpy  ["networkx"]=get_networkx ["biopython"]=get_biopython \
                       ["RATT"]=get_RATT ["nr database"]=get_nr ["TransTermHP"]=get_TransTermHP \
                       ["ViennaRNA package"]=get_Vienna_package ["BLAST+"]=get_blast_plus \
                       ["MEME"]=get_MEME ["MPICH"]=get_MPICH ["IntaRNA"]=get_IntaRNA ["Segemehl"]=get_segemehl \
                       ["Samtools and Bcftools"]=get_samtools_bcftools_htslib ["ImageMagick"]=get_ImageMagick \
                       ["Infernal"]=get_infernal ["BSRD database"]=get_BSRD ["GO database"]=get_GO_data \
                       ["STRING species.txt"]=get_STRING_species ["Rfam database"]=get_Rfam )
    declare -A VERs=( ["TSSpredator"]=1.06 ["CRT"]=1.2 ["matplotlib"]=1.5.0 ["matplotlib"]=1.9.2 \
                      ["nextworkx"]=1.10 ["biopython"]=1.65 ["RATT"]=1.64 ["TransTermHP"]=2.0.9 \
                      ["ViennaRNA package"]=2.3.2 ["BLAST+"]="2.2.28+" ["MEME"]=4.11.1 ["MPICH"]=3.2 \
                      ["IntaRNA"]=2.0.4 ["Segemehl"]=0.1.9 ["Samtools and Bcftools"]=1.3.1 ["ImageMagick"]="6.9.0-0" \
                      ["Infernal"]=1.1.1 ["nr database"]="NA" ["BSRD database"]="NA" ["GO database"]="NA" \
                      ["STRING species.txt"]="NA" ["Rfam database"]="NA" )
    for PROG in "${!PROGs[@]}"
    do
        printf "${RED}$PROG${NC} : \n \tIf ${RED}$PROG${NC} already exists and work properly, the installation can be skipped.\n " 
        if [[ ${VERs[$PROG]} != "NA" ]]
        then
            printf "\tBut please make sure the version is at least ${RED}${VERs[$PROG]}${NC}\n\t"
        else
            printf "\t"
        fi
        read -p "Do you wish to download/install it (yes/no)?" yn
        case $yn in
            [Yy]* ) (${PROGs[$PROG]});;
            [Nn]* ) ;;
            * ) echo "Please answer yes or no.";;
        esac
    done
}

get_nr(){
    wget -cP $DATABASE ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
    cd $DATABASE
    gunzip nr.gz
    mv nr nr.fa
    cd $PATH_FILE
}

get_BSRD(){
    wget -cP $DATABASE https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/database/sRNA_database_BSRD.fa
}

get_GO_data(){
    wget -cP $DATABASE http://www.geneontology.org/ontology/subsets/goslim_generic.obo
    wget -cP $DATABASE http://geneontology.org/ontology/go.obo
    wget -cP $DATABASE ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
    cd $DATABASE && gunzip ANNOgesic/input/databases/idmapping_selected.tab.gz
    cd $PATH_FILE
}

get_STRING_species(){
    VERSION=v10.5
    wget -cP ANNOgesic/input/databases https://string-db.org/download/species.$STRAIN.txt
}

get_Rfam(){
    #### download the riboswitch and RNA thermometer ID list of Rfam from our Git repository.
    wget -cP $DATABASE https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/database/Rfam_riboswitch_ID.csv
    wget -cP $DATABASE https://raw.githubusercontent.com/Sung-Huan/ANNOgesic/master/database/Rfam_RNA_thermometer_ID.csv
    #### download Rfam database
    VERSION=12.0
    wget -cP $DATABASE ftp://ftp.ebi.ac.uk/pub/databases/Rfam/$VERSION/Rfam.tar.gz
    cd $DATABASE
    tar -zxvf Rfam.tar.gz
    rm Rfam.tar.gz
    cd $PATH_FILE
}

get_matplotlib(){
    pip3 install --user --upgrade matplotlib
}

get_numpy(){
    pip3 install --user --upgrade numpy
}

get_biopython(){
    pip3 install --user --upgrade biopython
}

get_networkx(){
    pip3 install --user --upgrade networkx
}

get_RATT(){
    ##################################################
    # please add                                     #
    # export RATT_HOME=$TOOL_PATH/PAGIT/RATT         #
    # and                                            #
    # export PAGIT_HOME=$TOOL_PATH/bin/PAGIT         #
    # to the .bashrc                                 #
    # Please replace $TOOL_PATH to the absolute path #
    ##################################################
    VERSION=V1
    cd $TOOL_PATH
    git clone https://github.com/bioperl/bioperl-live.git
    cd bioperl-live
    perl Build.PL && ./Build test && ./Build install
    cd $PATH_FILE
    wget ftp://ftp.sanger.ac.uk/pub/resources/software/pagit/PAGIT.V1.64bit.tgz
    cd $TOOL_PATH && tar xzf PAGIT.{$VERSION}.64bit.tgz
    rm PAGIT.{$VERSION}.64bit.tgz
    sed -i '244s/defined//' $TOOL_PATH/PAGIT/RATT/main.ratt.pl && \
    cd $PATH_FILE
    echo "Please add the following two lines to your .bashrc file."
    printf "${RED}export RATT_HOMR=$TOOL_PATH/PAGIT/RATT\n"
    printf "export PAGIT_HOME=$TOOL_PATH/bin/PAGIT\n${NC}"
}

get_TSSpredator(){
    VERSION=1.06
    wget -cP $TOOL_PATH https://lambda.informatik.uni-tuebingen.de/nexus/content/repositories/releases/org/uni-tuebingen/it/TSSpredator/$VERSION/TSSpredator-$VERSION.jar
    cd $PATH_FILE
}

get_TransTermHP(){
    VERSION=v2.09
    wget -cP $TOOL_PATH http://transterm.cbcb.umd.edu/transterm_hp_${VERSION}.zip
    cd $TOOL_PATH
    unzip transterm_hp_${VERSION}.zip
    rm transterm_hp_${VERSION}.zip
    cd $PATH_FILE
}

get_Vienna_package(){
    VERSION=2.3.2
    wget -cP $TOOL_PATH http://www.tbi.univie.ac.at/RNA/packages/source/ViennaRNA-${VERSION}.tar.gz
    cd $TOOL_PATH
    tar -zxvf ViennaRNA-${VERSION}.tar.gz
    cd ViennaRNA-${VERSION}/ && ./configure --prefix=$HOME --without-perl --without-python && make && make install && cd ..
    rm ViennaRNA-${VERSION}.tar.gz
    cd $PATH_FILE
}

get_blast_plus(){
    VERSION="2.7.1+"
    if ! [ -e bin/makeblastdb ]
    then
        wget -cP $TOOL_PATH ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-${VERSION}-x64-linux.tar.gz
	cd $TOOL_PATH
        tar xfz ncbi-blast-${VERSION}-x64-linux.tar.gz
        rm -rf ncbi-blast-${VERSION}-x64-linux.tar.gz
	cd $PATH_FILE
    fi
}

get_sRNA_database(){
    #######################################################
    # please copy the database from BSRD manually.        #
    # http://www.bac-srna.org/BSRD/downloadIndexNew.jsp   #
    # Plase store the databse in Transap/input/database   #
    #######################################################
    echo please copy the database from BSRD manually.
    echo http://www.bac-srna.org/BSRD/downloadIndexNew.jsp
    echo Plase store the databse in Transap/input/database
}

get_MPICH(){
    #### Install MPICH for parallel running of MEME
    VERSION="3.2"
    wget -cP $TOOL_PATH http://www.mpich.org/static/downloads/$VERSION/mpich-$VERSION.tar.gz
    cd $TOOL_PATH && tar -zxvf mpich-$VERSION.tar.gz
    cd mpich-$VERSION && ./configure --prefix=$HOME && make all install
    cd $PATH_FILE

}

get_MEME(){
    #### Install MEME
    VERSION="4.11.1"
    wget -cP $TOOL_PATH http://meme-suite.org/meme-software/$VERSION/meme_$VERSION.tar.gz
    cd $TOOL_PATH
    tar -zxvf meme_${VERSION}.tar.gz

    perl -MCPAN -e 'install HTML::Template; \
    install XML::Compile::SOAP11; install XML::Compile::WSDL11; \
    install XML::Compile::Transport::SOAPHTTP; install XML::RPC'

    cpan install HTML::Template && \
    cpan HTML::PullParser && \
    cpan XML::Simple && \
    cpan XML::Compile::WSDL11 && \
    cpan XML::Compile::SOAP11

    cd meme_4.11.1 && ./configure --prefix=$HOME/meme --with-url=http://meme.nbcr.net/meme --enable-build-libxml2 --enable-build-libxslt && \
    make && make test && make install
    cd $TOOL_PATH && rm *.gz
    cd $PATH_FILE

}

get_IntaRNA(){
    VERSION=2.0.4
    wget -cP $TOOL_PATH https://github.com/BackofenLab/IntaRNA/releases/download/v$VERSION/intaRNA-$VERSION.tar.gz
    cd $TOOL_PATH && tar -zxvf intaRNA-$VERSION.tar.gz && cd intaRNA-$VERSION && ./configure --prefix=$HOME && make && make install
    cd $TOOL_PATH && rm intaRNA-$VERSION.tar.gz
    cd $PATH_FILE
}

get_segemehl(){
    ############################################################################################
    # If you want to use segemehl for detecting Circular RNA,                                  #
    # you need to generate(make) "testrealign.x" from segemehl.                                #
    # The following script will install segemehl and generate testrealign at the same time.    #
    ############################################################################################
    VERSION=0_2_0
    wget -cP $TOOL_PATH http://www.bioinf.uni-leipzig.de/Software/segemehl/segemehl_${VERSION}.tar.gz
    cd $TOOL_PATH
    tar -zxvf segemehl_${VERSION}.tar.gz
    cd segemehl_${VERSION}/segemehl/
    make all
    cd $TOOL_PATH
    rm segemehl_${VERSION}.tar.gz
    cd $PATH_FILE
}

get_samtools_bcftools_htslib(){
    cd $TOOL_PATH
    git clone git://github.com/samtools/samtools.git
    git clone git://github.com/samtools/htslib.git
    git clone git://github.com/samtools/bcftools.git
    cd samtools && make && make prefix=$HOME install && cd ..
    cd bcftools && make && make prefix=$HOME install && cd ..
    cd $PATH_FILE
}

get_CRT(){
    cd $TOOL_PATH
    wget http://www.room220.com/crt/CRT1.2-CLI.jar.zip && \
    unzip CRT1.2-CLI.jar.zip
    cd $PATH_FILE
}

get_ImageMagick(){
    VERSION=7.0.7-28
    cd $TOOL_PATH
    wget https://www.imagemagick.org/download/ImageMagick.tar.gz
    tar -zxvf ImageMagick.tar.gz && cd ImageMagick-$VERSION && ./configure --prefix=$HOME && make && make install && cd ..
    rm ImageMagick.tar.gz
    cd $PATH_FILE
}

get_infernal(){
    VERSION=1.1.1
    cd $TOOL_PATH
    wget http://eddylab.org/software/infernal/infernal-$VERSION.tar.gz
    tar -zxvf infernal-$VERSION.tar.gz && cd infernal-$VERSION && ./configure --prefix=$HOME && make && make install && cd ..
    rm infernal-$VERSION.tar.gz
    cd $PATH_PATH
}

main
