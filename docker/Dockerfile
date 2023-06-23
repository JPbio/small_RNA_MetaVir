# syntax=docker/dockerfile:1
 
# ====================================
# -- Stage: Perl --
# ====================================

# 
# TODO: 2023-06-15 - Lock versions!
# 
# Python 3.9.2
# Perl 5.36.0
# 
FROM docker.io/library/perl:5.36.0-threaded AS stage_perl

RUN cpanm --notest --quiet \
    Bio::SeqIO \
    Statistics::Basic \
    Statistics::RankCorrelation

# ====================================
# -- Stage: R --
# ====================================

# 
# TODO: 2023-06-15 - Lock versions!
# 
# R 4.0.4
# ggplot2 '3.4.2'
# reshape2 '1.4.4'
# 

FROM stage_perl AS stage_r

RUN apt-get update && \
    apt-get install -y \
        r-base && \
    Rscript -e "install.packages(c('ggplot2', 'reshape2', 'Rtsne', 'umap', 'ComplexHeatmap'), repos='https://cloud.r-project.org/')"

# ====================================
# -- Stage: Dependencies --
# ====================================

FROM stage_r AS stage_dependencies

# Add Blast
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.14.0/ncbi-blast-2.14.0+-x64-linux.tar.gz && \
        tar xvpf ncbi-blast-2.14.0+-x64-linux.tar.gz && \
        mv ncbi-blast-2.14.0+/bin/* /usr/bin/ && \
        rm -rf ncbi-blast-2.14.0+* && \
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz && \
        tar xvpf blast-2.2.26-x64-linux.tar.gz && \
        mv blast-2.2.26/bin/formatdb /usr/bin/ && \
        mv blast-2.2.26/bin/blastall /usr/bin/ && \
        rm -rf blast-2.2.26*

# 
# Add the ones available trhough apt-get
# TODO: 2023-06-23 - Lock versions
# 
RUN apt-get update && \
    apt-get install -y \
        bowtie=1.3.0+dfsg1-1 \
        velvetoptimiser=2.2.6-3 ncbi-blast+- \
        spades=3.13.1+dfsg-2+b2

# Add cap3
COPY ./cap3.linux.opteron64.tar /usr/bin/

RUN cd /usr/bin/ && \
    tar -xvf cap3.linux.opteron64.tar && \
    rm cap3.linux.opteron64.tar

ENV PATH="${PATH}:/usr/bin/CAP3"


# 
# Add FASTX-Toolkit
# TODO: 2023-06-16 - Lock fastx version
# 
ENV LD_LIBRARY_PATH="${LD_LIBRARY_PATH}/usr/local/lib"

RUN git clone https://github.com/agordon/libgtextutils.git && \
        cd libgtextutils/ && \
        ./reconf && \
        ./configure && \
        make && \
        make install && \
        cd .. && \
        rm -rf libgtextutils && \
    git clone https://github.com/agordon/fastx_toolkit.git && \
        cd fastx_toolkit && \
        wget -t 0 https://github.com/agordon/fastx_toolkit/files/1182724/fastx-toolkit-gcc7-patch.txt && \
        patch -p1 < fastx-toolkit-gcc7-patch.txt && \
        ./reconf && \
        ./configure && \
        make && \
        make install && \
        cd .. && \
        rm -rf fastx_toolkit && \
    ldconfig

# Add diamond
RUN wget http://github.com/bbuchfink/diamond/releases/download/v2.1.6/diamond-linux64.tar.gz && \
    tar xzf diamond-linux64.tar.gz && \
    mv ./diamond /usr/bin && \
    rm diamond-linux64.tar.gz

# Add fastqc
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip -d /usr/bin/fastqc && \
    rm fastqc_v0.12.1.zip

# ====================================
# -- Stage: Main --
# ====================================

FROM stage_dependencies AS stage_main

RUN mkdir /small-rna-metavir \
    && mkdir /small-rna-metavir/src \
    && mkdir /small-rna-metavir/asset \
    && chmod -R +x /small-rna-metavir/src

WORKDIR /small-rna-metavir

# ====================================
# -- Stage: Final --
# ====================================

FROM stage_main AS srna_metavir

RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# ====================================
# -- Stage: Final [dev] --
# ====================================

FROM stage_main AS small_rna_metavir_dev

RUN apt-get install -y nano
