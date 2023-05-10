# syntax=docker/dockerfile:1

# ====================================
# -- Stage: Perl --
# ====================================

FROM docker.io/library/perl:5.36.0 AS stage_perl
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
        bowtie=1.3.1-1

# ====================================
# -- Stage: Final --
# ====================================

# FROM docker.io/library/debian:11-slim AS srna_metavir
# FROM docker.io/library/ubuntu:20.04 AS srna_metavir
FROM docker.io/library/ubuntu:22.04 AS srna_metavir

# 
# REVIEW: 2023-03-09 - Can we include less stuff?
# REVIEW: 2023-03-10 - /usr/local/bin/perl seems to be the right one
# REVIEW: 2023-03-13 - /usr/lib/x86_64-linux-gnu seems to be critical
# 

# Include perl
COPY --from=stage_perl /etc/perl /etc/perl
COPY --from=stage_perl /usr/bin/perl /usr/bin/perl
COPY --from=stage_perl /usr/bin/perl5.32-x86_64-linux-gnu /usr/bin/perl5.32-x86_64-linux-gnu
COPY --from=stage_perl /usr/lib/x86_64-linux-gnu/perl /usr/lib/x86_64-linux-gnu/perl
COPY --from=stage_perl /usr/local/bin/perl /usr/local/bin/perl
COPY --from=stage_perl /usr/share/perl /usr/share/perl

COPY --from=stage_perl /usr/local/lib/perl5 /usr/local/lib/perl5
# /usr/share/man/man1/perl.1.gz

# Include R
COPY --from=stage_r /etc/R /etc/R
COPY --from=stage_r /usr/bin/R /usr/bin/R
COPY --from=stage_r /usr/lib/R /usr/lib/R
COPY --from=stage_r /usr/lib/x86_64-linux-gnu/* /usr/lib/x86_64-linux-gnu/
# COPY --from=stage_r /usr/local/lib/R /usr/local/lib/R
# COPY --from=stage_r /usr/share/R /usr/share/R

# COPY --from=stage_r /usr/lib/R/library /usr/lib/R/library
# COPY --from=stage_r /usr/lib/R/site-library /usr/lib/R/site-library
# REVIEW: 2023-03-13 - Can we get rid of these following two?
# COPY --from=stage_r /usr/local/lib/R/site-library /usr/local/lib/R/site-library
# COPY --from=stage_r /usr/lib/x86_64-linux-gnu /usr/lib/x86_64-linux-gnu

# COPY --from=stage_r /usr/bin/R /usr/bin/R
# COPY --from=stage_r /etc/R /etc/R
# COPY --from=stage_r /usr/share/R /usr/share/R

# COPY --from=stage_r /usr/lib/R /usr/lib/R
# COPY --from=stage_r /usr/lib/libR.so /usr/lib/libR.so
# COPY --from=stage_r /usr/local/lib/R /usr/local/lib/R
# COPY --from=stage_r /usr/lib64/ /usr/lib64/
# # COPY --from=stage_r /usr/lib/x86_64-linux-gnu /usr/lib/x86_64-linux-gnu
# COPY --from=stage_r /usr/lib/x86_64-linux-gnu/libblas.so.3 /usr/lib/x86_64-linux-gnu/libblas.so.3
# COPY --from=stage_r /usr/lib/x86_64-linux-gnu/libreadline.so.8 /usr/lib/x86_64-linux-gnu/libreadline.so.8
# COPY --from=stage_r /usr/lib/x86_64-linux-gnu/libicuuc.so /usr/lib/x86_64-linux-gnu/libicuuc.so
# COPY --from=stage_r /usr/lib/x86_64-linux-gnu/libtirpc.so /usr/lib/x86_64-linux-gnu/libtirpc.so

# # COPY --from=stage_r /usr/lib64/R/bin/exec/R /usr/lib64/R/bin/exec/R
# # /usr/lib/R
# # /usr/lib/R/bin/exec/R
# # /usr/lib/R/etc/ldpaths
# # /usr/lib64/R/bin/exec/R
# # /usr/share/R

# Include bowtie
COPY --from=stage_dependencies /usr/bin/bowtie* /usr/bin/

# Handle source files
RUN mkdir /srna_metavir \
    && mkdir /srna_metavir/src \
    && mkdir /srna_metavir/asset
    # REVIEW: 2023-03-09 - Mind this permission problem
    # && chmod +x -R src/utils

WORKDIR /srna_metavir

# ====================================
# -- Stage: Final [dev] --
# ====================================

FROM srna_metavir AS srna_metavir_dev
RUN apt-get update -y && apt-get install -y nano
