use strict;
use warnings;

use Getopt::Long;

my $usage = "

$0 - contig1 < contig 1, the biggest contig > -contig2 < contig 2, smaller contig > -output < output name >
    $0 - h

    -
    contig1 < contig1 >: The biggest contig to merge -
    contig2 < contig2 >: The smaller contig to merge -
    output < contigs merged >: Ouput file with merged contigs -
    h: Help message

";

$ | = 1;
# forces immediate prints into files rather than the buffer.

my $c1;
my $c2;
my $out;
my $help;

GetOptions("contig1=s" => \$c1,
    "contig2=s" => \$c2,
    "output=s" => \$out,
    "h=s" => \$help,
);

if (not(defined($c1))) {
    die "\nGive a valid format! \n", $usage;
}

if (not(defined($c2))) {
    die "\nGive a valid format! \n", $usage;
}

if (not(defined($out))) {
    die "\nGive a valid format! \n", $usage;
}

my $path_utils = "/srna_metavir/src/utils";
my $path_inner = "$path_utils/merge-contigs/inner";

#formatando contig 1
`formatdb -i $c2 -p F`;

#rodando blast
`blastall -p blastn -i $c1 -d $c2 -o step1.blastn -e 1e-5 -a 10`;

#Filtrando e gerando relatorio(80 % do tamanho com 80 % de identidade)
`$path_utils/filterblast.pl -b step1.blastn -pid 60 -plen 60 -evalue 1e-5 --best > step2.blastn.report`;
`$path_inner/selectBigContigFilterBlast.pl -i step2.blastn.report  -f1 $c1 -f2 $c2 > $out`;

#`rm -rf step*`;
print "\nDONE!\n";

#`rm -rf step*`;