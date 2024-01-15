FROM debian:bookworm-slim

MAINTAINER lingminhao 

# install all required CLI tools 
RUN apt-get update \ 
 && apt-get install -y wget make g++ git-all zlib1g zlib1g-dev r-base python-is-python3 python3-pip software-properties-common \ 
 ca-certificates gnupg2 libssl-dev libcurl4-gnutls-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \ 
 && apt-key adv --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' \ 
 && add-apt-repository -y 'deb http://cloud.r-project.org/bin/linux/debian bookworm-cran40/' \ 
 && add-apt-repository -y 'deb http://cloud.r-project.org/bin/linux/debian bookworm-cran40/' \ 
 && apt update \ 
 && apt-get install -y r-base

# make a directory to store all required softwares
RUN mkdir mnt/software

# install flexiplex 
RUN cd mnt/software && wget https://github.com/DavidsonGroup/flexiplex/releases/download/Version-0.97/flexiplex-Version-0.97.tar.gz \
 && tar -xvf flexiplex-Version-0.97.tar.gz && rm flexiplex-Version-0.97.tar.gz && cd flexiplex-Version-0.97 && make  

# install minimap2
RUN cd mnt/software && git clone https://github.com/lh3/minimap2 \ 
 && cd minimap2 && make 

# install samtools 
RUN cd mnt/software && wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 \
 && tar -xvf htslib-1.18.tar.bz2 && rm htslib-1.18.tar.bz2 && cd htslib-1.18 && make && make install \ 
 && wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 \
 && tar -xvf samtools-1.18.tar.bz2 && rm samtools-1.18.tar.bz2 && cd samtools-1.18 && make && make install

# install bambu 
RUN cd mnt/software && git clone -b messyforest --single-branch https://github.com/lingminhao/bambu.git \ 
 && R -e "install.packages('devtools')" \ 
 && R -e "install.packages('BiocManager')" && R -e "BiocManager::install('bambu')" \ 
 && R -e "install.packages('bambu', repos = NULL, type = 'source')"

RUN cd mnt/software && wget https://raw.githubusercontent.com/DavidsonGroup/flexiplex/863a5f3b182a0ee669407406b6cb39ae9f7f0c76/scripts/filter-barcodes.py \ 
 && pip3 install numpy pandas --break-system-packages

# environment variables 
ENV PATH=$PATH:/mnt/software/flexiplex-Version-0.97:/mnt/software/minimap2:/mnt/software/samtools-1.18:/mnt/software/htslib-1.18
