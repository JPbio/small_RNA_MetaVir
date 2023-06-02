use strict;
use warnings;

use Getopt::Long;

#######################################################################
### PARSE INPUTS ------------------------------------------------------
#######################################################################

my $usage = "

$0 - contig1 < contig 1, the biggest contig > -contig2 < contig 2, smaller contig > -output < output name >
    $0 - h

    -
    contig1 < contig1 >: The biggest contig to merge -
    contig2 < contig2 >: The smaller contig to merge -
    output < contigs merged >: Ouput file with merged contigs -
    h: Help message

";

$| = 1;
# forces immediate prints into files rather than the buffer.

my $c1;
my $c2;
my $out;
my $help;

GetOptions(
    "contig1=s" => \$c1,
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

#######################################################################
### Main... -----------------------------------------------------------
#######################################################################

# Set file paths
my $path_utils = "/srna_metavir/src/utils";
my $path_inner = "$path_utils/merge-contigs/inner";

my $path_filter_blast = "$path_utils/filterblast.pl";
my $path_select_big_contig = "$path_inner/selectBigContigFilterBlast.pl";

#formatando contig 1
`formatdb -i $c2 -p F`;

#rodando blast
`blastall -p blastn -i $c1 -d $c2 -o step1.blastn -e 1e-5 -a 10`;

# Filtering & generaring report (80% of the size with 80% of identity)
`perl $path_filter_blast -b step1.blastn -pid 60 -plen 60 -evalue 1e-5 --best > step2.blastn.report`;
`perl $path_select_big_contig -i step2.blastn.report  -f1 $c1 -f2 $c2 > $out`;

#`rm -rf step*`;
print "\nDONE!\n";

#`rm -rf step*`;