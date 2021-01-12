#################################################################
# Dockerfile to build bwa, MACS2, samtools, 
# picard-tools, fastQC, bedtools, cutadapt, R, multiqc
# images
# Based on Ubuntu
#  $ cd chipseq
#  $ VERSION=master
#  $ docker build -t jianhong/chipseq:$VERSION .  ## --no-cache
#  $ docker images jianhong/chipseq:$VERSION
#  $ docker push jianhong/chipseq:$VERSION
#  $ docker tag jianhong/chipseq:$VERSION jianhong/chipseq:latest
#  $ docker push jianhong/chipseq:latest
#  $ docker system prune -a
#  $ cd ~
#  $ docker pull jianhong/chipseq:latest
#  $ mkdir tmp4chipseq
#  $ docker run -it --rm -v ${PWD}/tmp4chipseq:/work jianhong/chipseq:latest
##################################################################
FROM ubuntu:20.10
LABEL authors="Jianhong Ou" \
      description="Docker image containing all software requirements for the jianhong/chipseq pipeline"

ENV DEBIAN_FRONTEND="noninteractive" TZ="America/New_York"

# Install software from source
RUN cd ~ && \
    apt-get update --fix-missing && \
    apt-get install --yes rsync wget bzip2 gcc libssl-dev libxml2-dev libncurses5-dev libbz2-dev liblzma-dev libcurl4-openssl-dev librsvg2-dev libv8-dev make cmake build-essential bedtools picard-tools python3 python3-pip pandoc fastqc multiqc bwa samtools bamtools && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
RUN ln -s python3 /usr/bin/python
RUN wget https://raw.githubusercontent.com/jianhong/chipseq/master/bin/picard -P /usr/bin/ && \
    chmod +x /usr/bin/picard

#RUN cd ~ && \
#    wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && \
#    tar -xf samtools-1.11.tar.bz2 && \
#    cd samtools-1.11 && \
#    ./configure && make && make install &&
#    cd .. && rm -rf samtool*

#RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && \
#    tar -xf bwa-0.7.17.tar.bz2 && \
#    cd bwa-0.7.17.tar.bz2 && \
#    sed -i '/const uint8_t rle_auxtab\[8\];/d' rle.h && \
#    make && cp bwa /usr/local/bin/ && \
#    cd .. && rm -rf bwa*

#RUN wget https://github.com/pezmaster31/bamtools/archive/v2.5.1.tar.gz && \
#    tar -xf v2.5.1.tar.gz && cd bamtools-2.5.1 && \
#    mkdir build && cd build && \
#    cmake -DCMAKE_INSTALL_PREFIX=/usr/local .. && \
#    make && make install && \
#    cd ../.. && rm -rf v2.5.1.tar.gz && rm -rf bamtools-2.5.1

RUN pip install pysam deeptools MACS2 cutadapt pymdown-extensions

RUN mkdir /homer && cd /homer && \
    wget http://homer.ucsd.edu/homer/configureHomer.pl && \
    perl configureHomer.pl -install
ENV PATH $PATH:/homer/bin

RUN cd ~ && wget https://raw.githubusercontent.com/gbcs-embl/Je/master/dist/je_2.0.RC.tar.gz && \
    tar -xf je_2.0.RC.tar.gz && cd je_2.0.RC && \
    sed -i "s/bin\/sh/usr\/bin\/env bash/" je && \
    cp * /usr/local/sbin/ && cd .. && rm -rf je*

RUN wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz && \
    tar -xf 0.6.6.tar.gz && cd TrimGalore-0.6.6 && \
    cp trim_galore /usr/local/sbin/ && cd .. && \
    rm 0.6.6.tar.gz && rm -rf TrimGalore-0.6.6

RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && \
    chmod +x bedGraphToBigWig && mv bedGraphToBigWig /usr/local/sbin/

## Install Bioconductor packages
RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e 'BiocManager::install(c("optparse", "rjson", "DiffBind", "ChIPpeakAnno", \
                "rtracklayer", "ggplot2", "GenomicFeatures", "DESeq2", "vsn", \
                "RColorBrewer", "pheatmap", "lattice", "BiocParallel", \
                "reshape2", "scales", "UpSetR", "caTools"))'

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron

WORKDIR /work
ENV JAVA_HOME="/usr"