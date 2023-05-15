# syntax=docker/dockerfile:1
# FROM debian:stable-slim as srna_metavir
# FROM stable-20230502-slim as srna_metavir

# RUN apt-get update && \
#     apt-get install -y \
#         build-essential \
#         curl
 
# ====================================
# -- Stage: Perl --
# ====================================

FROM docker.io/library/perl:5.36.0-threaded AS stage_perl

# # Download and extract the thread-enabled Perl source code
# RUN curl -L https://www.cpan.org/src/perl-5.32.1.tar.gz -o perl-5.32.1.tar.gz && \
#     tar xzf perl-5.32.1.tar.gz && \
#     cd perl-5.32.1 && \
#     ./Configure -des -Dusethreads && \
#     make && \
#     make install && \
#     cd .. && \
#     rm -rf perl-5.32.1.tar.gz perl-5.32.1

# # Install CPANM
# RUN curl -L https://cpanmin.us -o /usr/local/bin/cpanm && \
#     chmod +x /usr/local/bin/cpanm

# RUN apt-get install -y bioperl
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

# =============================================
# DESC 01

# In this Dockerfile, each software package is installed in a separate RUN instruction to keep the Docker image layers as small as possible.

# The Perl installation includes some required packages for database access. The Python installation includes numpy, pandas, and biopython. The R installation includes ggplot2, reshape2, dplyr, and tidyr.

# Bowtie and VelvetOptimiser are installed in their own RUN instructions. Finally, the default command is set to run Bash.

# Note that you can further optimize the image size by combining multiple RUN instructions into a single one and removing any unnecessary files and directories.

# =============================================
# DESC 02

# In this updated Dockerfile, we first install build-essential to be able to compile Perl from source. Then, we download and extract the thread-enabled Perl source code. We configure Perl with the -Dusethreads flag to enable threading support and compile it. Finally, we install Perl using make install.

# The rest of the Dockerfile remains the same as in the previous example. Note that you can adjust the version of Perl as needed by changing the URL in the curl command and the directory name in the cd command.