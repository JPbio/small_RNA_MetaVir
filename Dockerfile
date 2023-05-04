# syntax=docker/dockerfile:1

# ====================================
# -- Stage: Perl --
# ====================================

FROM docker.io/library/perl:5.36.0 AS stage_perl
RUN cpan Bio::SeqIO


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

FROM docker.io/library/ubuntu:22.04 AS stage_dependencies
RUN apt-get update -y && apt-get install -y \
    python3=3.10.6-1~22.04 \
    bowtie=1.3.1-1


# ====================================
# -- Stage: Final --
# ====================================

# 
# REVIEW: 2023-03-09 - Find a smaller working base image
# 

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
COPY --from=stage_r /usr/local/lib/R /usr/local/lib/R
COPY --from=stage_r /usr/share/R /usr/share/R

COPY --from=stage_r /usr/lib/R/library /usr/lib/R/library
COPY --from=stage_r /usr/lib/R/site-library /usr/lib/R/site-library
# REVIEW: 2023-03-13 - Can we get rid of these following two?
COPY --from=stage_r /usr/local/lib/R/site-library /usr/local/lib/R/site-library
COPY --from=stage_r /usr/lib/x86_64-linux-gnu/* /usr/lib/x86_64-linux-gnu/

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

# RUN apk add --no-cache bash vim
RUN apt-get update -y && apt-get install -y nano
