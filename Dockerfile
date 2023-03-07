FROM perl:5.36.0

RUN mkdir /srna_metavir
WORKDIR /srna_metavir

RUN mkdir asset
RUN mkdir src
# RUN chmod +x -R src/utils # TODO: Mind this problem

RUN apt-get update -y \
    && apt-get install -y \
    bowtie
    # vim \

RUN cpan Bio::SeqIO