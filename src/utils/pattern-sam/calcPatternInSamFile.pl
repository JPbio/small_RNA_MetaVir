use strict;
use warnings;

use Bio::SeqIO;  
use Bio::Seq::Quality;
use Getopt::Long;
use Statistics::Basic qw(:all);
use Statistics::RankCorrelation;

#
# TODO: 2023-05-15 - Use a decent logger
# 
use constant PATH_LOG_MAIN => "srna_metavir.main.log";

#######################################################################
### PARSE INPUTS ------------------------------------------------------
#######################################################################

my $usage = "

$0  -s <sam file> -o <output> 
$0 -h


-s  <sam file>			: Sam file from mapping on putative sequences
-o 	<output file>        : output file name
-h                      : Help message

";

$| = 1;     # forces immediate prints into files rather than the buffer.


use Bio::SeqIO;
use Bio::Seq::Quality;
use Getopt::Long;
use Statistics::Basic qw(:all);
use Statistics::RankCorrelation;

my $search;
my $output;
my $help;

GetOptions(
    "s=s"     => \$search,
    "o=s"     => \$output,
    "h!"      => \$help
);

if ($help) {
    die $usage;
}

if (not(defined($search))) {
    die "\nGive a search file name", $usage;
}

#######################################################################
### Configure logging -------------------------------------------------
#######################################################################

# 
# TODO: 2023-02-27 - Find a better way to do this...
# TODO: 2023-02-27 - Restablish the custom log file(s) option
# 

# open(metrics, ">$step8/full_metrics.txt");
# open(interest, ">$step8/metrics_of_interest.txt");
# open(LOG, ">$log");

# 
# NOTE: From here on all printed stuff will be sent to the log file
# 

# open filehandle log.txt
my $LOG_FH;
open($LOG_FH, ">>", PATH_LOG_MAIN) or die "Couldn't open: $!"; # $! is a special variable holding the error
# select $LOG_FH;

#######################################################################
### Main... -----------------------------------------------------------
#######################################################################

my $path_utils = "/srna_metavir/src/utils";
my $path_inner = "$path_utils/pattern-sam/inner";

my $path_plot_pattern = "$path_inner/plotPattern_publication.R";
my $path_plot_pattern_both_strands = "$path_inner/plotPattern_publication_both_strands.R";

my %tamanhos_p;
my %tamanhos_n;

# my $cp = 0;
# my $cn = 0;
my %ch;
my %ch2;

print $LOG_FH "\nAnalyzing $search...\n";

open(IN, "<$search");
while (<IN>) {

    if (/^\w+/) {
        
        my @campos = split(/\t/, $_);
        my $size   = length($campos[9]);
        my $c      = $campos[2];

        if (not exists $tamanhos_p{$size}) {
            $tamanhos_p{$size} = 0;
        }
        if (not exists $ch{$c}{"pos"}) {
            $ch{$c}{"pos"} = 0;
        }

        if (not exists $ch{$c}{"neg"}) {
            $ch{$c}{"neg"} = 0;
        }

        if (not exists $tamanhos_n{$size}) {
            $tamanhos_n{$size} = 0;
        }

        if (not exists $ch2{$c}{"pos"}{$size}) {
            $ch2{$c}{"pos"}{$size} = 0;
        }

        if (not exists $ch2{$c}{"neg"}{$size}) {
            $ch2{$c}{"neg"}{$size} = 0;
        }

        if ($campos[1] eq "0") {
            # $cp++;
            $ch{$c}{"pos"} = $ch{$c}{"pos"} + 1;
            $tamanhos_p{$size} = $tamanhos_p{$size} + 1;
            $ch2{$c}{"pos"}{"$size"} = $ch2{$c}{"pos"}{"$size"} + 1;

        } elsif ($campos[1] eq "16") {
            # $cn++;
            $ch{$c}{"neg"} = $ch{$c}{"neg"} + 1;
            $tamanhos_n{$size} = $tamanhos_n{$size} + 1;
            $ch2{$c}{"neg"}{$size} = $ch2{$c}{"neg"}{$size} + 1;
            #print $LOG_FH "NEGATIVE in [$c] fita [pos]  size [$size] amount: ".$ch{$c}{"neg"}{$size} ."\n";
        }
    }
}
# $t = $cp + $cn;

print $LOG_FH "Calculating Z-score...\n";

open(ALL, ">$output");

print ALL "Chromosome#";

for (my $i = 0; $i <= 41; $i++) {
    print ALL "$i#";
}

print ALL "\n";

foreach my $chr (keys %ch2) {    #for each sequence
    
    my @tot;

    # Initializing array
    for (my $i = 0; $i <= 41; $i++) {
        $tot[$i] = 0;
    }

    #  open(O,">score.$chr.tab");
    #  print O "$chr\n";
    for (my $i = 0; $i <= 41; $i++) {
        
        my $p = 0;
        my $n = 0;
        
        if (exists $ch2{$chr}{"pos"}{$i}) {
            $p = $ch2{$chr}{"pos"}{$i};
        }
        if (exists $ch2{$chr}{"neg"}{$i}) {
            $n = $ch2{$chr}{"neg"}{$i};
        }

        $tot[$i - 15] = $p;
        $tot[$i + 6]  = $n;
        # print  "$i\t$p\t$n\t$t\n";
    }

    ##z = (x - mean)/sd

    my $stddev_t = stddev(@tot);
    my $mean_t   = mean(@tot);

    if ($stddev_t > 0 and $mean_t) {

        my %zscore;
        my @zvt;
        
        for (my $i = 0; $i <= 41; $i++) {
            $zscore{$i}{"both"} = ($tot[$i] - $mean_t) / $stddev_t;
            push(@zvt, ($tot[$i] - $mean_t) / $stddev_t);
        }
        
        print ALL "$chr#";
        
        for (my $i = 0; $i <= 41; $i++) {
            # print O $zscore{$i}{"both"}."\n";
            print ALL $zscore{$i}{"both"} . "#";
        }
        
        print ALL "\n";
    }
    
    # close(O);
}
close(ALL);

`R --no-save $output $output < $path_plot_pattern`;
`R --no-save $output $output < $path_plot_pattern_both_strands`;
