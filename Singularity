From:r-base:latest
Bootstrap:docker

%labels
    MAINTAINER Jianhong Ou <jianhong.ou@duke.edu>
    DESCRIPTION Singularity image containing all requirements for the jianhong/chipseq
    VERSION latest

%environment
  DEBIAN_FRONTEND="noninteractive"
  TZ="America/New_York"
  JAVA_HOME="/usr"
  export PATH $PATH:/homer/bin:/edirect

%post
  apt-get update --fix-missing
  apt-get install --yes rsync wget bzip2 gcc libssl-dev libxml2-dev libncurses5-dev libbz2-dev liblzma-dev libcurl4-openssl-dev librsvg2-dev libv8-dev make cmake build-essential bedtools picard-tools python3 python3-pip pandoc fastqc multiqc bwa samtools bamtools subread pigz curl libxml-simple-perl
  
  ln -s python3 /usr/bin/python
  
  wget https://raw.githubusercontent.com/jianhong/chipseq/master/assets/picard -P /usr/bin/ && \
    chmod +x /usr/bin/picard
  
  pip install pysam deeptools MACS2 cutadapt pymdown-extensions trackhub
  
  mkdir /homer && cd /homer && \
    wget http://homer.ucsd.edu/homer/configureHomer.pl && \
    perl configureHomer.pl -install
  
  cd ~ && wget https://raw.githubusercontent.com/gbcs-embl/Je/master/dist/je_2.0.RC.tar.gz && \
    tar -xf je_2.0.RC.tar.gz && cd je_2.0.RC && \
    sed -i "s/bin\/sh/usr\/bin\/env bash/" je && \
    cp * /usr/local/sbin/ && cd .. && rm -rf je*
  
  wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz && \
    tar -xf 0.6.6.tar.gz && cd TrimGalore-0.6.6 && \
    cp trim_galore /usr/local/sbin/ && cd .. && \
    rm 0.6.6.tar.gz && rm -rf TrimGalore-0.6.6
  
  wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && \
    chmod +x bedGraphToBigWig && mv bedGraphToBigWig /usr/local/sbin/
  
  yes | sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)" && \
    mv $HOME/edirect /edirect && \
    wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz && \
    tar -xf sratoolkit.tar.gz && mv sratoolkit.2.10.9-ubuntu64/bin/* /usr/local/sbin/ && \
    rm -rf sratoolkit*
  
  wget https://github.com/jianhong/chipseq/archive/dev.zip && \
    unzip dev.zip && mv chipseq-dev /pipeline && rm dev.zip
  
  Rscript -e "install.packages('BiocManager')"
  
  Rscript -e 'BiocManager::install(c("optparse", "rjson", "DiffBind", "ChIPpeakAnno", \
                "rtracklayer", "ggplot2", "GenomicFeatures", "DESeq2", "vsn", \
                "RColorBrewer", "pheatmap", "lattice", "BiocParallel", \
                "reshape2", "scales", "UpSetR", "caTools", \
                "rmarkdown", "DT"))'
  
  touch .Rprofile
  touch .Renviron
  
  
