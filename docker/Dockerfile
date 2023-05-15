# syntax=docker/dockerfile:1
 
# ====================================
# -- Stage: Perl --
# ====================================

FROM docker.io/library/perl:5.36.0-threaded AS stage_perl

RUN cpanm --notest --quiet Bio::SeqIO

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
    r-base

RUN Rscript -e "install.packages(c('ggplot2', 'reshape2'), repos='https://cloud.r-project.org/')"

# ====================================
# -- Stage: Dependencies --
# ====================================

FROM stage_r AS stage_dependencies

RUN apt-get update && \
    apt-get install -y \
    bowtie \
    velvetoptimiser

# ====================================
# -- Stage: Main --
# ====================================

FROM stage_dependencies AS stage_main

RUN mkdir /srna_metavir \
    && mkdir /srna_metavir/src \
    && mkdir /srna_metavir/asset

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
