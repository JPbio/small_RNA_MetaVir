# syntax=docker/dockerfile:1

# ====================================
# -- Stage: Perl --
# ====================================

FROM perl:5.36.0 AS stage_perl
RUN cpan Bio::SeqIO


# ====================================
# -- Stage: R --
# ====================================

FROM r-base:4.2.2 as stage_r
RUN R -e "install.packages('remotes')" \
    && installGithub.r tidyverse/ggplot2@v3.4.1


# ====================================
# -- Stage: Dependencies --
# ====================================

FROM ubuntu:22.04 AS stage_dependencies
RUN apt-get update -y && apt-get install -y \
    # TODO: 2023-03-09 - Lock bowtie version
    bowtie


# ====================================
# -- Stage: Final --
# ====================================

# 
# REVIEW: 2023-03-09 - Find a smaller working base image
# 

FROM ubuntu:22.04 AS srna_metavir

# 
# REVIEW: 2023-03-09 - Can we include less stuff
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

# Include bowtie
# COPY --from=stage_dependencies /usr/bin/bowtie* /usr/bin/

# Habdle source files
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

# RUN apk add --no-cache bash vim
RUN apt-get update -y && apt-get install -y nano
