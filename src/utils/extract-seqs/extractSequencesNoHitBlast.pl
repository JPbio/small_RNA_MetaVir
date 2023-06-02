use strict;
use warnings;

use Getopt::Long;

my $usage = "

$0 -seq <original sequences> -out <out file> -blast <blast report after filterblast.pl>
$0 -h

-seq <reads .csfasta>         : Original sequences query on blast
-blast <reads .qual>          : report generated after filterblast.pl
-out <outputput name> 	      :  output file
-h                            : Help message

";

$| = 1;  #    forces immediate prints into files rather than the buffer.

my $blast;
my $seq;
my $out;
my $help;

GetOptions(
    "seq=s" => \$seq,
    "blast=s" => \$blast,
    "out=s" => \$out,
    "h=s" => \$help,
);

if (not(defined($blast))) {
    die "\nGive a valid format! \n", $usage;
}

if (not(defined($seq))) {
    die "\nGive a valid format! \n", $usage;
}

#######################################################################
### Main... -----------------------------------------------------------
#######################################################################

my $path_utils = "/srna_metavir/src/utils";
my $path_extract_seqs_fasta = "$path_utils/extract-seqs/inner/extract_seqs_from_fasta_file.pl";

#getIdsWithHit
`grep -oP 'Query=\\S+' $blast  | cut -d '=' -f 2 | sort -u > idsWithHit`;

#pegando todos os ids dos contigs1
`grep -oP '^>\\S+' $seq | cut -d '>' -f 2 | sort -u > idsOriginalSequences`;

#diff
`diff -u idsOriginalSequences idsWithHit > difHitAndSeq.diff`;

#Fazendo diff para ver quais contigs não estão nos contigs2
`grep -oP '^\\-\\S+' difHitAndSeq.diff  | cut -d '-' -f 2 > idsSeqNoHit`;

#extraindo contigs
`perl $path_extract_seqs_fasta idsSeqNoHit $seq $out`;

#` rm -rf idsWithHit idsOriginalSequences difHitAndSeq.diff idsSeqNoHit`;