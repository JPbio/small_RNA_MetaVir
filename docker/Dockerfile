# syntax=docker/dockerfile:1

# ====================================
# -- Stage: Perl --
# ====================================

FROM docker.io/library/perl:5.36.0-threaded AS stage_perl
# docker pull perl:5.36.0-threaded-bullseye
RUN cpanm --notest --quiet Bio::SeqIO

# ====================================
# -- Stage: Python --
# ====================================

FROM docker.io/library/python:3.10-slim-bullseye AS stage_python
# COPY requirements.txt .
# RUN pip install --no-cache-dir -r requirements.txt

# ====================================
# -- Stage: R --
# ====================================

FROM docker.io/library/r-base:4.2.2 as stage_r
RUN R -e "install.packages('remotes')" \
    && installGithub.r tidyverse/ggplot2@v3.4.1 \
    && installGithub.r hadley/reshape@v1.4.1

# ====================================
# -- Stage: Dependencies --
# ====================================

# FROM debian:bullseye-slim
# FROM docker.io/library/debian:11-slim AS stage_dependencies
# RUN apt-get update && \
#     apt-get install --no-install-recommends --yes \
#         bowtie=1.3.0+dfsg1-1

# FROM docker.io/library/ubuntu:20.04 AS stage_dependencies
# RUN apt-get update && \
#     apt-get install --no-install-recommends -y \
#         bowtie=1.2.3+dfsg-4build1

FROM docker.io/library/ubuntu:22.04 AS stage_dependencies
RUN apt-get update && \
    apt-get install --no-install-recommends -y \
        # wget \
        bowtie=1.3.1-1

# ====================================
# -- Stage: Final --
# ====================================

# FROM docker.io/library/debian:11-slim AS srna_metavir
# FROM docker.io/library/ubuntu:20.04 AS srna_metavir
FROM docker.io/library/ubuntu:22.04 AS srna_metavir

# Include perl
COPY ./remove_perl.sh /
RUN chmod +x /remove_perl.sh && \
    /remove_perl.sh && \
    rm /remove_perl.sh

COPY --from=stage_perl /etc/perl /etc/perl
COPY --from=stage_perl /usr/bin/perl* /usr/bin/perl/
# COPY --from=stage_perl /usr/lib/perl* /usr/lib/perl/
COPY --from=stage_perl /usr/lib/x86_64-linux-gnu/libperl* /usr/lib/x86_64-linux-gnu/libperl/
COPY --from=stage_perl /usr/lib/x86_64-linux-gnu/perl-base /usr/lib/x86_64-linux-gnu/perl-base
COPY --from=stage_perl /usr/lib/x86_64-linux-gnu/perl* /usr/lib/x86_64-linux-gnu/perl/
# COPY --from=stage_perl /usr/lib/x86_64-linux-gnu/perl/5.34 /usr/lib/x86_64-linux-gnu/perl/5.34
# COPY --from=stage_perl /usr/lib/x86_64-linux-gnu/perl5/5.34 /usr/lib/x86_64-linux-gnu/perl5/5.34
COPY --from=stage_perl /usr/local/bin/perl /usr/local/bin/perl
COPY --from=stage_perl /usr/local/lib/perl* /usr/local/lib/perl/
COPY --from=stage_perl /usr/local/lib/perl5 /usr/local/lib/perl5
# COPY --from=stage_perl /usr/local/lib/site_perl /usr/local/lib/site_perl
# COPY --from=stage_perl /usr/local/lib/x86_64-linux-gnu/perl/5.34.0 /usr/local/lib/x86_64-linux-gnu/perl/5.34.0
# COPY --from=stage_perl /usr/local/share/perl* /usr/local/share/perl/
# COPY --from=stage_perl /usr/local/share/perl/5.34.0 /usr/local/share/perl/5.34.0
COPY --from=stage_perl /usr/share/perl* /usr/share/perl/

# COPY --from=stage_perl /usr/bin/perl*
# COPY --from=stage_perl /usr/lib/perl*
# COPY --from=stage_perl /usr/lib/x86_64-linux-gnu/libperl*
# COPY --from=stage_perl /usr/lib/x86_64-linux-gnu/perl-base
# COPY --from=stage_perl /usr/lib/x86_64-linux-gnu/perl*
# COPY --from=stage_perl /usr/lib/x86_64-linux-gnu/perl/5.32
# COPY --from=stage_perl /usr/lib/x86_64-linux-gnu/perl/5.34
# COPY --from=stage_perl /usr/lib/x86_64-linux-gnu/perl5/5.32
# COPY --from=stage_perl /usr/lib/x86_64-linux-gnu/perl5/5.34
# COPY --from=stage_perl /usr/local/lib/perl*
# COPY --from=stage_perl /usr/local/lib/site_perl
# COPY --from=stage_perl /usr/local/lib/x86_64-linux-gnu/perl/5.32.1
# COPY --from=stage_perl /usr/local/lib/x86_64-linux-gnu/perl/5.34.0
# COPY --from=stage_perl /usr/local/share/perl*
# COPY --from=stage_perl /usr/local/share/perl/5.32.1
# COPY --from=stage_perl /usr/local/share/perl/5.34.0
# COPY --from=stage_perl /usr/share/perl*



# Include python
COPY --from=stage_python /usr/local/bin/python3 /usr/local/bin/python3
COPY --from=stage_python /usr/local/lib /usr/local/lib

# Include R
COPY --from=stage_r /etc/R /etc/R
COPY --from=stage_r /usr/bin/R /usr/bin/R
COPY --from=stage_r /usr/lib/R /usr/lib/R
COPY --from=stage_r /usr/lib/x86_64-linux-gnu/* /usr/lib/x86_64-linux-gnu/
COPY --from=stage_r /usr/local/lib/R /usr/local/lib/R

# Include bowtie
COPY --from=stage_dependencies /usr/bin/bowtie* /usr/bin/

# Handle source files
RUN mkdir /srna_metavir \
    && mkdir /srna_metavir/src \
    && mkdir /srna_metavir/asset
    # REVIEW: 2023-03-09 - Mind this permission problem
    # && chmod +x -R src/utils

WORKDIR /srna_metavir





RUN apt-get update && \
    apt-get install --no-install-recommends -y \
        wget
# RUN wget -O velvetoptimiser_2.2.6.tar.gz https://github.com/tseemann/VelvetOptimiser/releases/tag/2.2.6 && \
RUN wget -O velvetoptimiser_2.2.6.tar.gz https://codeload.github.com/tseemann/VelvetOptimiser/tar.gz/refs/tags/2.2.6 --no-check-certificate && \
    mkdir /src && \
    tar zxvf velvetoptimiser_2.2.6.tar.gz && \
    mv VelvetOptimiser-2.2.6/ /src/velvetOptimiser-2.2.6/
    # && \
    # echo "export PATH=$PATH:/src/velvetOptimiser-2.2.6/"
ENV PATH="$PATH:/src/velvetOptimiser-2.2.6/"


# ====================================
# -- Stage: Final [dev] --
# ====================================

FROM srna_metavir AS srna_metavir_dev

# RUN apt-get update && \
#     apt-get install --yes \
#         nano

RUN apt-get update -y \
    && apt-get install -y \
        nano
