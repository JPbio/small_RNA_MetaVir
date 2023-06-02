# syntax=docker/dockerfile:1
 
# ====================================
# -- Stage: Perl --
# ====================================

FROM docker.io/library/perl:5.36.0-threaded AS stage_perl

# RUN cpanm --notest --quiet \
#     Bio::SeqIO \
#     Statistics::Basic \
#     Statistics::RankCorrelation

RUN cpanm --notest --quiet Bio::SeqIO
RUN cpanm --notest --quiet Statistics::Basic
RUN cpanm --notest --quiet Statistics::RankCorrelation

# ====================================
# -- Stage: Python --
# ====================================

FROM stage_perl AS stage_python

# RUN apt-get update && \
#     apt-get install -y \
#         python3-pip \
#         && \
#     pip3 install \
#         numpy pandas biopython

# ====================================
# -- Stage: R --
# ====================================

FROM stage_python AS stage_r

RUN apt-get update && \
    apt-get install -y \
        r-base && \
    Rscript -e "install.packages(c('ggplot2', 'reshape2'), repos='https://cloud.r-project.org/')"

# ====================================
# -- Stage: Dependencies --
# ====================================

FROM stage_r AS stage_dependencies

# Add Blast
# RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz && \
#         tar xvpf ncbi-blast-2.13.0+-x64-linux.tar.gz && \
#         mv ncbi-blast-2.13.0+/bin/* /usr/bin/ && \
#         rm -rf ncbi-blast-2.13.0+* && \
#     wget https://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz
#         tar xvpf blast-2.2.26-x64-linux.tar.gz && \
#         mv blast-2.2.26/bin/formatdb /usr/bin/ && \
#         mv blast-2.2.26/bin/blastall /usr/bin/ && \
#         rm -rf blast-2.2.26*

# 
# TODO: 2023-06-01 - Use blast 2014
# 

RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz && \
    tar xvpf ncbi-blast-2.13.0+-x64-linux.tar.gz && \
    mv ncbi-blast-2.13.0+/bin/* /usr/bin/ && \
    rm -rf ncbi-blast-2.13.0+*

RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz && \
        tar xvpf blast-2.2.26-x64-linux.tar.gz && \
        mv blast-2.2.26/bin/formatdb /usr/bin/ && \
        mv blast-2.2.26/bin/blastall /usr/bin/ && \
        rm -rf blast-2.2.26*

# Add the ones available trhough apt-get
# RUN apt-get update && \
#     apt-get install -y \
#         bowtie=1.3.0+dfsg1-1 \
#         velvetoptimiser=2.2.6-3 ncbi-blast+- \
#         spades=3.13.1+dfsg-2+b2
RUN apt-get update
RUN apt-get install -y bowtie=1.3.0+dfsg1-1
RUN apt-get install -y velvetoptimiser=2.2.6-3 ncbi-blast+-
RUN apt-get install -y spades=3.13.1+dfsg-2+b2

# Add cap3
COPY ./cap3.linux.opteron64.tar /usr/bin/

RUN cd /usr/bin/ && \
    tar -xvf cap3.linux.opteron64.tar && \
    rm cap3.linux.opteron64.tar

ENV PATH="${PATH}:/usr/bin/CAP3"

# Add FASTX-Toolkit
ENV LD_LIBRARY_PATH="${LD_LIBRARY_PATH}/usr/local/lib"

# RUN git clone https://github.com/agordon/libgtextutils.git && \
#         cd libgtextutils/ && \
#         ./reconf && \
#         ./configure && \
#         make && \
#         make install && \
#         cd .. && \
#         rm -rf libgtextutils && \
#     git clone https://github.com/agordon/fastx_toolkit.git && \
#         cd fastx_toolkit && \
#         wget -t 0 https://github.com/agordon/fastx_toolkit/files/1182724/fastx-toolkit-gcc7-patch.txt && \
#         patch -p1 < fastx-toolkit-gcc7-patch.txt && \
#         ./reconf && \
#         ./configure && \
#         make && \
#         make install && \
#         cd .. && \
#         rm -rf fastx_toolkit && \
#     ldconfig

RUN mkdir /temp && \
    cd /temp && \
        git clone https://github.com/agordon/libgtextutils.git && \
            cd /temp/libgtextutils/ && \
            ./reconf && \
            ./configure && \
            make && \
            make install && \
    cd /temp && \
        git clone https://github.com/agordon/fastx_toolkit.git && \
            cd /temp/fastx_toolkit && \
            wget -t 0 https://github.com/agordon/fastx_toolkit/files/1182724/fastx-toolkit-gcc7-patch.txt && \
            patch -p1 < fastx-toolkit-gcc7-patch.txt && \
            ./reconf && \
            ./configure && \
            make && \
            make install && \
    cd / && \
    rm -rf /temp && \
    ldconfig

# Add diamond
RUN wget http://github.com/bbuchfink/diamond/releases/download/v2.1.6/diamond-linux64.tar.gz && \
    tar xzf diamond-linux64.tar.gz && \
    mv ./diamond /usr/bin && \
    rm diamond-linux64.tar.gz

# ====================================
# -- Stage: Main --
# ====================================

FROM stage_dependencies AS stage_main

RUN mkdir /srna_metavir \
    && mkdir /srna_metavir/src \
    && mkdir /srna_metavir/asset \
    && mkdir /srna_metavir/runs \
    && chmod -R +x /srna_metavir/src

WORKDIR /srna_metavir

# ====================================
# -- Stage: Final --
# ====================================

FROM stage_main AS srna_metavir

RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# ====================================
# -- Stage: Final [dev] --
# ====================================

FROM stage_main AS srna_metavir_dev

RUN apt-get install -y nano
