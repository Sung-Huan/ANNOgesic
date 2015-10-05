FROM debian
MAINTAINER Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get upgrade --yes
RUN apt-get update
RUN apt-get install --yes \
default-jre \
default-jdk \
python3 \
python3-scipy \
vim \
make \
gcc \
g++ \
gfortran \
libx11-dev \
wget \
zip unzip \
python3-biopython \
software-properties-common \
python3-software-properties \
bioperl \
ncbi-blast+ \
pkg-config \
python3-dev \
libfreetype6-dev \
libpng-dev \
python3-pip \
python3-numpy \
imagemagick \
infernal \
git \
openssh-client \
apache2 \
curl \
build-essential \
net-tools \
librpc-xml-perl \
ncbi-blast+-legacy \
nano \
libf2c2 \
apache2-dev \
libapache-singleton-perl \
libjson-rpc-perl \
libncurses5-dev

RUN pip3 install \
matplotlib \
networkx \
ANNOgesic

RUN pip3 install --upgrade ANNOgesic
            
RUN mkdir tools
WORKDIR tools

# vienna package
# RUN apt-add-repository ppa:j-4/vienna-rna --yes
# RUN apt-get update --yes
# RUN apt-get install vienna-rna --yes
RUN wget http://www.tbi.univie.ac.at/RNA/packages/source/ViennaRNA-2.1.9.tar.gz && \
tar -zxvf ViennaRNA-2.1.9.tar.gz && cd ViennaRNA-2.1.9 && ./configure && make && make install && \
cp Utils/relplot.pl /usr/local/bin && \
cp Utils/mountain.pl /usr/local/bin

# TSSpredator
RUN wget http://it.informatik.uni-tuebingen.de/software/tsspredator/TSSpredator_v1-04.zip && \
unzip TSSpredator_v1-04.zip && \
cp TSSpredator_v1-04/TSSpredator.jar /usr/local/bin

# MEME
RUN wget http://ebi.edu.au/ftp/software/MEME/4.10.1/meme_4.10.1_1.tar.gz && \
tar -zxvf meme_4.10.1_1.tar.gz

RUN perl -MCPAN -e 'install HTML::Template; \
install XML::Compile::SOAP11; install XML::Compile::WSDL11; \
install XML::Compile::Transport::SOAPHTTP; install XML::RPC'

RUN cd meme_4.10.1 && ./configure --prefix=/tools/meme \
--with-url=http://meme.nbcr.net/meme \
--enable-build-libxml2 \
--enable-build-libxslt && \
make && make test && make install && cp /tools/meme/bin/meme /usr/local/bin

# RATT
RUN wget ftp://ftp.sanger.ac.uk/pub/resources/software/pagit/PAGIT.V1.64bit.tgz && \
tar xzf PAGIT.V1.64bit.tgz && ./installme.sh

ENV PAGIT_HOME /tools/PAGIT
ENV PATH /tools/PAGIT/bin/:/tools/PAGIT/bin/pileup_v0.5/:/tools/PAGIT/bin/pileup_v0.5/ssaha2:\
/tools/PAGIT/bin/pileup_v0.5/:/tools/PAGIT/IMAGE/:/tools/PAGIT/ABACAS:\
/tools/PAGIT/ICORN/:/tools/PAGIT/RATT/:$PATH

ENV PILEUP_HOME /tools/PAGIT/bin/pileup_v0.5/
ENV ICORN_HOME /tools/PAGIT/ICORN/
ENV SNPOMATIC_HOME /tools/PAGIT/bin/
ENV RATT_HOME /tools/PAGIT/RATT
ENV RATT_CONFIG $RATT_HOME/RATT.config
ENV PERL5LIB /usr/lib/perl5/:/tools/PAGIT/lib
RUN cp PAGIT/RATT/start.ratt.sh /usr/local/bin

# segemehl
RUN wget http://www.bioinf.uni-leipzig.de/Software/segemehl/segemehl_0_2_0.tar.gz && \
tar -zxvf segemehl_0_2_0.tar.gz && cd segemehl_0_2_0/segemehl && \
make all && cp *.x /usr/local/bin

# transtermHP
RUN wget http://transterm.cbcb.umd.edu/transterm_hp_v2.09.zip && \
unzip transterm_hp_v2.09.zip && cd transterm_hp_v2.09 && \
make && cp transterm /usr/local/bin && cp expterm.dat /usr/local/bin

# Psortb
RUN DEBIAN_FRONTEND=noninteractive apt-get update && \
apt-get -y install supervisor && apt-get clean && \
rm -rf /var/lib/apt/lists/*
ENV APACHE_RUN_USER www-data
ENV APACHE_RUN_GROUP www-data
ENV APACHE_LOG_DIR /var/log/apache2
ENV APACHE_LOCK_DIR /var/lock/apache2
ENV APACHE_PID_FILE /var/run/apache2.pid

WORKDIR /usr/local/src
RUN echo '/usr/local/lib64' >>/etc/ld.so.conf
RUN wget http://www.psort.org/download/docker/pft2.3.4.docker64bit.tar.gz && \
tar zxvf pft2.3.4.docker64bit.tar.gz && cp pftools/pfscan /usr/local/bin/

RUN wget http://www.psort.org/download/libpsortb-1.0.tar.gz && \
tar zxvf libpsortb-1.0.tar.gz && cd libpsortb-1.0 && \
./configure && make && make install && ldconfig

RUN wget http://www.psort.org/download/bio-tools-psort-all.3.0.4.tar.gz && \
tar zxvf bio-tools-psort-all.3.0.4.tar.gz
WORKDIR /usr/local/src/bio-tools-psort-all

RUN wget http://www.psort.org/download/docker/psortb.defaults

RUN perl Makefile.PL && make && make install && cp -r psort /usr/local/psortb

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
bio-tools-psort-all \
/tools/ViennaRNA-2.1.9

RUN /etc/init.d/apache2 restart
EXPOSE 80
CMD ["/opt/run.sh"]

# copy psort to global execute
RUN cp /usr/local/psortb/bin/psort /usr/local/bin

WORKDIR /tools

# htslib, samtools, bcftools
RUN wget https://github.com/samtools/htslib/releases/download/1.2.1/htslib-1.2.1.tar.bz2
RUN tar -jxvf htslib-1.2.1.tar.bz2 && cd htslib-1.2.1 && make all && make install && cd ..
RUN wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2
RUN tar -jxvf samtools-1.2.tar.bz2 && cd samtools-1.2 && make all && make install && cd ..
RUN wget https://github.com/samtools/bcftools/releases/download/1.2/bcftools-1.2.tar.bz2
RUN tar -jxvf bcftools-1.2.tar.bz2 && cd bcftools-1.2 && make all && make install && cd ..

# replace the old version of samtools in PAGIT
RUN cp /tools/samtools-1.2/samtools /tools/PAGIT/bin/samtools

RUN rm TSSpredator_v1-04.zip \
meme_4.10.1_1.tar.gz \
PAGIT.V1.64bit.tgz \
segemehl_0_2_0.tar.gz \
transterm_hp_v2.09.zip \
ViennaRNA-2.1.9.tar.gz \
htslib-1.2.1.tar.bz2 \
samtools-1.2.tar.bz2 \
bcftools-1.2.tar.bz2
