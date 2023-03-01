FROM perl:5.36.0

RUN mkdir /srna_metavir
WORKDIR /srna_metavir

RUN mkdir src
RUN mkdir asset

RUN apt-get update -y && apt-get install -y \
    bowtie \
    # vim \