FROM ubuntu
MAINTAINER Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>
#ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update --yes
RUN apt-get upgrade --yes
RUN apt-get install default-jre default-jdk python3 python3-scipy \
vim make gcc g++ gfortran libx11-dev wget zip unzip python3-biopython \
software-properties-common python3-software-properties bioperl \
ncbi-blast+ pkg-config python3-dev libfreetype6-dev libxft-dev \
libpng-dev python3-pip python-pip python3-numpy imagemagick infernal git \
openssh-client apache2 curl build-essential net-tools librpc-xml-perl \
ncbi-blast+-legacy nano libf2c2 apache2-dev libapache-singleton-perl \
libjson-rpc-perl libncurses5-dev build-essential hmmer lua5.1 blast2 \
snap cpanminus mummer exonerate mafft fasttree libsvg-perl \
libgd-svg-perl python-setuptools libc6-i386 lib32stdc++6 lib32gcc1 \
netcat genometools last-align libboost-iostreams-dev libgsl2 libgsl-dev \
libcolamd2.9.1 liblpsolve55-dev libstdc++6 aragorn tantan libstorable-perl \
libbio-perl-perl libsqlite3-dev tree --yes --fix-missing
RUN ln -fs /usr/bin/fasttree /usr/bin/FastTree
RUN apt-get update --yes && apt-get upgrade --yes 

RUN pip3 install \
matplotlib \
networkx \
ANNOgesic

RUN mkdir tools
WORKDIR tools

# vienna package
RUN wget http://www.tbi.univie.ac.at/RNA/packages/source/ViennaRNA-2.2.5.tar.gz && \
tar -zxvf ViennaRNA-2.2.5.tar.gz && cd ViennaRNA-2.2.5 && ./configure  --without-perl --without-python && make && make install && \
cp src/Utils/relplot.pl /usr/local/bin && \
cp src/Utils/mountain.pl /usr/local/bin

# TSSpredator
RUN wget https://lambda.informatik.uni-tuebingen.de/nexus/content/repositories/releases/org/uni-tuebingen/it/TSSpredator/1.06/TSSpredator-1.06.jar && \
cp TSSpredator-1.06.jar /usr/local/bin/TSSpredator.jar

# MEME
RUN wget http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz && \
tar -zxvf mpich-3.2.tar.gz && \
cd mpich-3.2 && ./configure && make all install

RUN wget http://meme-suite.org/meme-software/4.11.1/meme_4.11.1.tar.gz && \
tar -zxvf meme_4.11.1.tar.gz

RUN perl -MCPAN -e 'install HTML::Template; \
install XML::Compile::SOAP11; install XML::Compile::WSDL11; \
install XML::Compile::Transport::SOAPHTTP; install XML::RPC'

RUN cpan install HTML::Template && \
cpan HTML::PullParser && \
cpan XML::Simple && \
cpan XML::Compile::WSDL11 && \
cpan XML::Compile::SOAP11

RUN cd meme_4.11.1 && ./configure --prefix=/tools/meme \
--with-url=http://meme.nbcr.net/meme \
--enable-build-libxml2 \
--enable-build-libxslt && \
make && make test && make install && cp /tools/meme/bin/* /usr/local/bin

# segemehl
RUN wget http://www.bioinf.uni-leipzig.de/Software/segemehl/segemehl_0_2_0.tar.gz && \
tar -zxvf segemehl_0_2_0.tar.gz && cd segemehl_0_2_0/segemehl && \
make all && cp *.x /usr/local/bin

# transtermHP
RUN wget http://transterm.cbcb.umd.edu/transterm_hp_v2.09.zip && \
unzip transterm_hp_v2.09.zip && cd transterm_hp_v2.09 && \
make && cp transterm /usr/local/bin && cp expterm.dat /usr/local/bin

# CRT
RUN wget http://www.room220.com/crt/CRT1.2-CLI.jar.zip && \
unzip CRT1.2-CLI.jar.zip && cp CRT1.2-CLI.jar /usr/local/bin/CRT.jar

# Psortb
RUN DEBIAN_FRONTEND=noninteractive apt-get update && \
apt-get -y install supervisor && apt-get clean && \
rm -rf /var/lib/apt/lists/*
ENV APACHE_RUN_USER=www-data \
APACHE_RUN_GROUP=www-data \
APACHE_LOG_DIR=/var/log/apache2 \
APACHE_LOCK_DIR=/var/lock/apache2 \
APACHE_PID_FILE=/var/run/apache2.pid

WORKDIR /usr/local/src
RUN echo '/usr/local/lib64' >>/etc/ld.so.conf && \
wget http://www.psort.org/download/docker/pft2.3.4.docker64bit.tar.gz && \
tar zxvf pft2.3.4.docker64bit.tar.gz && cp pftools/pfscan /usr/local/bin/

RUN wget http://www.psort.org/download/libpsortb-1.0.tar.gz && \
tar zxvf libpsortb-1.0.tar.gz && cd libpsortb-1.0 && \
./configure && make && make install && ldconfig

RUN wget http://www.psort.org/download/bio-tools-psort-all.3.0.4.tar.gz && \
tar zxvf bio-tools-psort-all.3.0.4.tar.gz
WORKDIR /usr/local/src/bio-tools-psort-all

RUN wget http://www.psort.org/download/docker/psortb.defaults && \
perl Makefile.PL && make && make install && cp -r psort /usr/local/psortb

RUN a2enmod cgid && \
wget http://www.psort.org/download/docker/apache.conf.fragment && \
cat apache.conf.fragment >> /etc/apache2/apache2.conf

WORKDIR /usr/local/src

RUN wget http://www.psort.org/download/docker/apache-svm.tar.gz && \
tar zxvf apache-svm.tar.gz && cd apache-svm && make && \
cp svmloc.conf /etc/apache2/conf-available/

RUN wget http://www.psort.org/download/docker/startup.txt && \
mv startup.txt startup.pl && \
wget http://www.psort.org/download/docker/apache-psort.conf && \
cp apache-psort.conf /etc/apache2/conf-available/

RUN wget http://www.psort.org/download/docker/apache-psortb.tar.gz && \
tar zxvf apache-psortb.tar.gz && cd apache-psortb \
&& perl Makefile.PL && make && make install

RUN cd /etc/apache2/conf-enabled/ && \
ln -s ../conf-available/svmloc.conf && \
ln -s ../conf-available/apache-psort.conf

RUN wget http://www.psort.org/download/docker/Request.pm && \
cp Request.pm /usr/share/perl5/Apache/Singleton/Request.pm

RUN wget http://www.psort.org/download/docker/CGI-FastTemplate-1.09.tar.gz && \
tar zxvf CGI-FastTemplate-1.09.tar.gz && \
cd CGI-FastTemplate-1.09 && perl Makefile.PL && make && make install

RUN cd /var/www/html && \
wget http://www.psort.org/download/docker/psort-web.tar.gz && \
tar zxvf psort-web.tar.gz

RUN rm -r pft2.3.4.docker64bit.tar.gz \
libpsortb-1.0.tar.gz \
libpsortb-1.0 \
bio-tools-psort-all.3.0.4.tar.gz \
bio-tools-psort-all

RUN /etc/init.d/apache2 restart
EXPOSE 80
CMD ["/opt/run.sh"]

# copy psort to global execute
RUN cp /usr/local/psortb/bin/psort /usr/local/bin

WORKDIR /tools

# htslib, samtools, bcftools
RUN wget https://github.com/samtools/htslib/releases/download/1.3.1/htslib-1.3.1.tar.bz2 && \
tar -jxvf htslib-1.3.1.tar.bz2 && cd htslib-1.3.1 && make all && make install && cd ..

RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
tar -jxvf samtools-1.3.1.tar.bz2 && cd samtools-1.3.1 && make all && make install && cd ..

RUN wget https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2 && \
tar -jxvf bcftools-1.3.1.tar.bz2 && cd bcftools-1.3.1 && make all && make install && cd ..

# RATT
RUN git clone https://github.com/sanger-pathogens/rapid_annotation_transfer_tool.git && \
mv rapid_annotation_transfer_tool /opt/RATT
# patch the error of perl version and the path of mummer
RUN sed -i '244s/defined//' /opt/RATT/main.ratt.pl && \
sed -i '19s/$PAGIT_HOME/\/usr/' /opt/RATT/start.ratt.sh

ENV RATT_HOME=/opt/RATT \
PERL5LIB=/opt/RATT/:$PERL5LIB \
PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/RATT:$PATH

RUN rm meme_4.11.1.tar.gz \
segemehl_0_2_0.tar.gz \
transterm_hp_v2.09.zip \
ViennaRNA-2.2.5.tar.gz \
htslib-1.3.1.tar.bz2 \
samtools-1.3.1.tar.bz2 \
bcftools-1.3.1.tar.bz2 \
CRT1.2-CLI.jar.zip 

RUN pip3 install ANNOgesic --upgrade

WORKDIR /home
