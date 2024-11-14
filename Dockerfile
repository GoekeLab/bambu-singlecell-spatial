FROM debian:bookworm-slim

MAINTAINER lingminhao 

# install all required CLI tools 
RUN apt-get update \ 
 && apt-get install -y wget make g++ git-all zlib1g zlib1g-dev r-base python-is-python3 python3-pip software-properties-common \ 
 ca-certificates gnupg2 libssl-dev libcurl4-gnutls-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev

RUN wget -q https://cloud.r-project.org/bin/linux/debian/key_4.1.2 -O cran-key.asc \
 && gpg --import cran-key.asc \
 && mv cran-key.asc /etc/apt/trusted.gpg.d/cran.gpg

RUN add-apt-repository -y 'deb http://cloud.r-project.org/bin/linux/debian bookworm-cran40/' \ 
 && apt update \ 
 && apt-get install -y r-base

# make a directory to store all required softwares
RUN mkdir mnt/software

# install flexiplex 
RUN cd mnt/software && wget https://github.com/DavidsonGroup/flexiplex/releases/download/v1.01/flexiplex-1.01.tar.gz \
 && tar -xvf flexiplex-1.01.tar.gz && rm flexiplex-1.01.tar.gz && cd flexiplex-1.01 && make  

# install jaffa
RUN cd mnt/software && wget https://github.com/Oshlack/JAFFA/releases/download/version-2.3/JAFFA-version-2.3.tar.gz \
 && tar -xvf JAFFA-version-2.3.tar.gz && rm JAFFA-version-2.3.tar.gz && cd JAFFA-version-2.3 $$ ./install_linux64.sh

# install minimap2
RUN cd mnt/software && git clone https://github.com/lh3/minimap2 \ 
 && cd minimap2 && make 

# install samtools 
RUN cd mnt/software && wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 \
 && tar -xvf htslib-1.18.tar.bz2 && rm htslib-1.18.tar.bz2 && cd htslib-1.18 && make && make install \ 
 && wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 \
 && tar -xvf samtools-1.18.tar.bz2 && rm samtools-1.18.tar.bz2 && cd samtools-1.18 && make && make install

# install bambu  
RUN cd mnt/software && R -e "install.packages('devtools')" 
RUN cd mnt/software && R -e "install.packages('R.utils')" 
RUN cd mnt/software && R -e "install.packages('BiocManager')" 
RUN cd mnt/software && R -e "BiocManager::install('bambu')"
RUN cd mnt/software && git clone -b messyForest --single-branch https://github.com/GoekeLab/bambu.git
RUN cd mnt/software && R -e "library('devtools'); load_all('bambu')"
RUN cd mnt/software && R -e "install.packages('Rcpp')"

RUN cd mnt/software && wget https://raw.githubusercontent.com/DavidsonGroup/flexiplex/863a5f3b182a0ee669407406b6cb39ae9f7f0c76/scripts/filter-barcodes.py \ 
 && pip3 install numpy pandas --break-system-packages

RUN cd mnt/software/bambu && chmod 777 DESCRIPTION

RUN apt-get install -y libopenblas-dev

# environment variables 
ENV PATH=$PATH:/mnt/software/flexiplex-1.01:/mnt/software/minimap2:/mnt/software/samtools-1.18:/mnt/software/htslib-1.18:/mnt/software/bambu
