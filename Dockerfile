FROM debian:latest

# install R 4.4.1 and all CLI tools 
ENV R_VERSION=4.4.1 \
    DEBIAN_FRONTEND=noninteractive

RUN apt-get update \ 
 && apt-get install -y wget make g++ git-all zlib1g zlib1g-dev r-base python-is-python3 python3-pip software-properties-common \ 
 ca-certificates gnupg2 libssl-dev libcurl4-gnutls-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \ 
 && apt update \ 
 && apt-get install -y r-base

RUN apt-get update -qq && apt-get -y install --no-install-recommends \
    ca-certificates \
    build-essential \
    gfortran \
    libreadline-dev \
    xorg-dev \
    libbz2-dev \
    liblzma-dev \
    curl \
    git-all \ 
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libsodium-dev \
    wget \
    && rm -rf /var/lib/apt/lists/*

RUN wget -c https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz \
    && tar -xf R-${R_VERSION}.tar.gz \
    && cd R-${R_VERSION} \
    && ./configure \
    && make -j$(nproc) \
    && make install \
    && cd .. \
    && rm -rf R-${R_VERSION} R-${R_VERSION}.tar.gz

# make a directory to store all required softwares
RUN mkdir mnt/software
    
# install flexiplex 
RUN cd mnt/software && wget https://github.com/DavidsonGroup/flexiplex/releases/download/v1.01/flexiplex-1.01.tar.gz \
 && tar -xvf flexiplex-1.01.tar.gz && rm flexiplex-1.01.tar.gz && cd flexiplex-1.01 && make  

# install minimap2
RUN cd mnt/software && git clone https://github.com/lh3/minimap2 \ 
 && cd minimap2 && make 

# install samtools 
RUN cd mnt/software && wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 \
 && tar -xvf htslib-1.18.tar.bz2 && rm htslib-1.18.tar.bz2 && cd htslib-1.18 && make && make install \ 
 && wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 \
 && tar -xvf samtools-1.18.tar.bz2 && rm samtools-1.18.tar.bz2 && cd samtools-1.18 && make && make install

# install bambu  
# RUN R -e "install.packages('R.utils', repos = 'http://cran.us.r-project.org')" 
RUN R -e "install.packages('BiocManager', repos = 'http://cran.us.r-project.org')" 
RUN R -e "BiocManager::install('bambu')"

# install Seurat 
RUN R -e "install.packages('SeuratObject', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('Seurat', repos = 'http://cran.us.r-project.org')"

# mischellous setup
RUN R -e "install.packages('devtools', repos = 'http://cran.us.r-project.org')" 
RUN cd mnt/software && git clone -b Multiplex_Major_Patch --single-branch https://github.com/GoekeLab/bambu.git
RUN cd mnt/software && wget https://raw.githubusercontent.com/DavidsonGroup/flexiplex/main/scripts/flexiplex_filter/main.py \ 
&& pip3 install numpy==1.26.4 pandas --break-system-packages

RUN cd mnt/software/bambu && chmod 777 DESCRIPTION
RUN R -e "library(devtools); load_all('/mnt/software/bambu')"

# install jaffa 
RUN cd mnt/software && wget https://github.com/Oshlack/JAFFA/releases/download/version-2.3/JAFFA-version-2.3.tar.gz \
 && tar -xvf JAFFA-version-2.3.tar.gz && rm JAFFA-version-2.3.tar.gz && cd JAFFA-version-2.3 $$ ./install_linux64.sh

#RUN apt-get install -y libopenblas-dev

# environment variables 
ENV PATH=$PATH:/mnt/software/flexiplex-1.01:/mnt/software/minimap2:/mnt/software/samtools-1.18:/mnt/software/htslib-1.18:/mnt/software/bambu
ENV R_LIBS_USER="/usr/local/lib/R/site-library"    
    
