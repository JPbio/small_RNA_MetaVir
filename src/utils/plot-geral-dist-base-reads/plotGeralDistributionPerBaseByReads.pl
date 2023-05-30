use strict;
use warnings;

use Bio::SeqIO;
use Bio::Seq::Quality;
use Getopt::Long;

# 
# TODO: 2023-03-07 - Use a decent logger
# 
use constant PATH_LOG_MAIN => "srna_metavir.main.log";

#######################################################################
### PARSE INPUTS ------------------------------------------------------
#######################################################################

my $usage = "

$0 -sam <sam file> -s <start> -e <end> -p prefix [-norm <size library>] [--plot] [-m min]
$0 -h

-i <input file>         : Input file with columns delimited by a specific delimiter
-h                      : Help message

";

$| = 1;     # forces immediate prints into files rather than the buffer.


# 
# TODO: 2023-03-07 - Find a better way to do this
# 
my $path_inner = "/srna_metavir/src/utils/plot-geral-dist-base-reads/inner";

my $sam;
my $delimiter;
my $index;
my $prefix;
my $profile;
my $fasta;
my $keep;
my $pattern;
my $start;
my $end;
my $plot;
my $antisense;
my $revcomp;

# 
# REVIEW: 2023-03-13 - Should these variables be 'my'
# 

our $norm;
our $minimo;

my $base;
my $help;

GetOptions(
    "sam=s"      => \$sam,
    "pattern!"   => \$pattern,
    "profile!"   => \$profile,
    "p=s"        => \$prefix,
    "s=i"        => \$start,
    "e=i"        => \$end,
    "norm=i"     => \$norm,
    "m=i"        => \$minimo,
    "h!"         => \$help,
    "keep!"      => \$keep,
    "plot!"      => \$plot,
    "antisense!" => \$antisense,
    "revcomp!"   => \$revcomp,
    "base!"      => \$base
);

if ($help) {
    die $usage;
}

if (not(defined($sam))) {
    die "\nGive a sam input file name", $usage;
}
if (not(defined($prefix))) {
    die "\nGive an prefix file name", $usage;
}

if (not defined($start) or $start < 0) {
    warn("\n\n\t\t[Bad value for start, start setted to 15]\n");
    $start = 15;
}

if (not defined($end) or ($end < 0) or ($end < $start)) {
    warn("\n\t\t[Bad value for end, end setted to 30]\n");
    $end = 30;
}
if (not defined($minimo) or ($minimo < 0)) {
    warn("\n\t\t[Bad value for min, min setted to 30]\n");
    $minimo = 30;
} else {
    warn("Minimo setted to $minimo\n");
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
select $LOG_FH;

#######################################################################
### Main... -----------------------------------------------------------
#######################################################################

open(IN, "<$sam");

# 
# REVIEW: 2023-03-10 - Should all these variables be 'my'?
# 

our %tamanhos;
our %tamanhos_p;
our %tamanhos_n;
our $cp = 0;
our $cn = 0;
our %ch;
our %ch2;
our %base;
our %base_total;
our %base_per_size;
my $min = 1000;
my $max = 0;

our %base_si;
our %base_pi;
# our %base;
our %total_reads;
our %total_si;
our %total_pi;

sub norm {
    my $t = shift;
    my $a = ($t / $norm) * 1000000;
    return $a;

}

my $spos = "0";
my $sneg = "16";

if (defined($antisense)) {
    $spos = "16";
    $sneg = "0";

}
for (my $i = 15 ; $i <= 35 ; $i++) {
    $base_per_size{$i}{"+"}{"A"} = 0;
    $base_per_size{$i}{"+"}{"C"} = 0;
    $base_per_size{$i}{"+"}{"G"} = 0;
    $base_per_size{$i}{"+"}{"T"} = 0;
    $base_per_size{$i}{"-"}{"A"} = 0;
    $base_per_size{$i}{"-"}{"C"} = 0;
    $base_per_size{$i}{"-"}{"G"} = 0;
    $base_per_size{$i}{"-"}{"T"} = 0;

}

warn("Loading SAM file...\n");

my %read;
my $nmappings = 0;
my $nreads = 0;

while (<IN>) {
    if (/^\w+/) {
        
		$nmappings += 1;
        my @campos = split(/\t/, $_);
        my $size   = length($campos[9]);
        my $c      = $campos[2];
        
		if (not exists $read{ $campos[0] }) {
			
			$read{ $campos[0] } = 1;
            $nreads += 1;
            
			if (not exists $tamanhos_p{$size}) {
                $tamanhos_p{$size}   = 0;
                $base_total{$c}{"A"} = 0;
                $base_total{$c}{"C"} = 0;
                $base_total{$c}{"G"} = 0;
                $base_total{$c}{"T"} = 0;

                $base_per_size{$size}{"+"}{"A"} = 0;
                $base_per_size{$size}{"+"}{"C"} = 0;
                $base_per_size{$size}{"+"}{"G"} = 0;
                $base_per_size{$size}{"+"}{"T"} = 0;
                $base_per_size{$size}{"-"}{"A"} = 0;
                $base_per_size{$size}{"-"}{"C"} = 0;
                $base_per_size{$size}{"-"}{"G"} = 0;
                $base_per_size{$size}{"-"}{"T"} = 0;

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
                $ch2{$c}{"pos"}{$size}       = 0;
                $base{$c}{"pos"}{$size}{"A"} = 0;
                $base{$c}{"pos"}{$size}{"C"} = 0;
                $base{$c}{"pos"}{$size}{"G"} = 0;
                $base{$c}{"pos"}{$size}{"T"} = 0;
            }

            if (not exists $ch2{$c}{"neg"}{$size}) {
                $ch2{$c}{"neg"}{$size}       = 0;
                $base{$c}{"neg"}{$size}{"A"} = 0;
                $base{$c}{"neg"}{$size}{"C"} = 0;
                $base{$c}{"neg"}{$size}{"G"} = 0;
                $base{$c}{"neg"}{$size}{"T"} = 0;
            }

            if ($size < $min) {
                $min = $size;
            }
            if ($size > $max) {
                $max = $size;
            }

            my $first = uc(substr($campos[9], 0, 1));

            if ($campos[1] eq 16) {
                my $seq = uc($campos[9]);
                $seq =~ tr /ATCG/TAGC/;
                my $revseq = reverse $seq;
                $first = uc(substr($revseq, 0, 1));
            }

            ###count miRNA preference

            if (not exists($base_total{$first})) {
                $base_total{$first} = 0;
            }
            $base_total{$first} = $base_total{$first} + 1;

            if (not exists($base{$c}{$first})) {
                $base{$c}{$first} = 0;
            }
            $base{$c}{$first} = $base{$c}{$first} + 1;

            if ($size >= 20 && $size <= 23) {
                if (not exists($base_si{$c}{$first})) {
                    $base_si{$c}{$first} = 0;
                }
                $base_si{$c}{$first} = $base_si{$c}{$first} + 1;

                if (not exists $total_si{$c}) { $total_si{$c} = 0; }
                $total_si{$c} = $total_si{$c} + 1;
            }
            if ($size >= 24 && $size <= 30) {
                if (not exists($base_pi{$c}{$first})) {
                    $base_pi{$c}{$first} = 0;
                }
                $base_pi{$c}{$first} = $base_pi{$c}{$first} + 1;

                if (not exists $total_pi{$c}) { $total_pi{$c} = 0; }
                $total_pi{$c} = $total_pi{$c} + 1;
            }

            if ($campos[1] eq $spos) {
                $cp++;
                $ch{$c}{"pos"}           = $ch{$c}{"pos"} + 1;
                $tamanhos_p{$size}       = $tamanhos_p{$size} + 1;
                $ch2{$c}{"pos"}{"$size"} = $ch2{$c}{"pos"}{"$size"} + 1;

                $base{$c}{"pos"}{$size}{$first} = $base{$c}{"pos"}{$size}{$first} + 1;

				#print "POSITIVE in [$c] fita [pos]  size [$size] amount: ".$ch{$c}{"pos"}{$size} ."\n";
                $base_per_size{$size}{"+"}{$first} = $base_per_size{$size}{"+"}{$first} + 1;
            } elsif ($campos[1] eq $sneg) {
                $cn++;
                $ch{$c}{"neg"}         = $ch{$c}{"neg"} + 1;
                $tamanhos_n{$size}     = $tamanhos_n{$size} + 1;
                $ch2{$c}{"neg"}{$size} = $ch2{$c}{"neg"}{$size} + 1;

                $base{$c}{"neg"}{$size}{$first} = $base{$c}{"neg"}{$size}{$first} + 1;

				#print "NEGATIVE in [$c] fita [pos]  size [$size] amount: ".$ch{$c}{"neg"}{$size} ."\n";
                $base_per_size{$size}{"-"}{$first} = $base_per_size{$size}{"-"}{$first} + 1;
            }
            if (not exists $total_reads{$c}) {
                $total_reads{$c} = 0;
            }
            $total_reads{$c} = $total_reads{$c} + 1;
        }
    }

}

my $t = $cp + $cn;

print "Tn: $cn \t Tp: $cp \t Tt: $t \n";
print "si: %total_si \npi: %total_pi \ntotal: %total_reads\n";

#
#Load fasta file
#

my %contig;

#
#writing file with total mappings
#

if (defined($norm)) {
    open(O7, ">$prefix.norm.profile1");
    open(O8, ">$prefix.norm.profile2");

} else {
    open(O7, ">$prefix.profile1");
    open(O8, ">$prefix.profile2");
}

print O8 "$prefix";

for (my $i = $start ; $i <= $end ; $i++) {
    my $p = 0;
    my $n = 0;
    if (exists($tamanhos_p{$i})) {
        $p = $tamanhos_p{$i};
    }
    if (exists($tamanhos_n{$i})) {
        $n = $tamanhos_n{$i};
    }
    $t = $p + $n;
    if (defined($norm)) {
        my $a = ($t / $norm) * 1000000;
        my $pnorm = ($p / $norm) * 1000000;
        my $nnorm = ($p / $norm) * 1000000;

        $p = $pnorm;
        $n = $nnorm;

        $t = $a;
    }
    print O7 "$i\t$p\t$n\t$t\n";
    print O8 "#$t";

}
print O8 "\n";
close(O7);
close(O8);

open(O5, ">$prefix.profile_per_size");

for (my $i = $start ; $i <= $end ; $i++) {
    my $p = 0;
    my $n = 0;
    if (exists($tamanhos_p{$i})) {
        $p = $tamanhos_p{$i};
    }
    if ($i == $start) {
        print O5 "$prefix#$p";
    } else {
        print O5 "#$p";
    }

}
for (my $i = $start ; $i <= $end ; $i++) {
    my $n = 0;
    if (exists($tamanhos_n{$i})) {
        $n = $tamanhos_n{$i};
    }
    print O5 "#$n";

}
print O5 "\n";
close(O5);
###################

#
#Print total of reads per chromosome
#

##############

###miRNA preference print
open(O3, ">$prefix.miRNA_total_base_distribution");
print O3 "\tA\tC\tG\tT\n";
foreach my $c (keys %ch2) {
    if (defined($base_total{$c})) {
        my $ta = $base_total{$c}{"A"};
        my $tc = $base_total{$c}{"C"};
        my $tg = $base_total{$c}{"G"};
        my $tt = $base_total{$c}{"T"};
        print O3 $c . "\t$ta\t$tc\t$tg\t$tt\n";
     }
}
close(O3);

############## GERAL preference print
open(O4, ">$prefix.Geral_base_distribution");
print O4 "size\tA\tC\tG\tT\n";
for (my $i = $start; $i <= $end ; $i++) {
    my $ta = $base_per_size{$i}{"+"}{"A"};
    my $tc = $base_per_size{$i}{"+"}{"C"};
    my $tg = $base_per_size{$i}{"+"}{"G"};
    my $tt = $base_per_size{$i}{"+"}{"T"};

    my $nta = $base_per_size{$i}{"-"}{"A"} * -1;
    my $ntc = $base_per_size{$i}{"-"}{"C"} * -1;
    my $ntg = $base_per_size{$i}{"-"}{"G"} * -1;
    my $ntt = $base_per_size{$i}{"-"}{"T"} * -1;

    print O4 $i . "\t$ta\t$tc\t$tg\t$tt\n";
    print O4 $i . "\t$nta\t$ntc\t$ntg\t$ntt\n";
}
close(O4);

############## GERAL preference print
if (defined($norm)) {
    open(O4, ">$prefix.Geral_base_distribution_norm");
    print O4 "size\tA\tC\tG\tT\n";
    for (my $i = $start ; $i <= $end ; $i++) {
        my $ta = norm($base_per_size{$i}{"+"}{"A"});
        my $tc = norm($base_per_size{$i}{"+"}{"C"});
        my $tg = norm($base_per_size{$i}{"+"}{"G"});
        my $tt = norm($base_per_size{$i}{"+"}{"T"});

        my $nta = norm($base_per_size{$i}{"-"}{"A"} * -1);
        my $ntc = norm($base_per_size{$i}{"-"}{"C"} * -1);
        my $ntg = norm($base_per_size{$i}{"-"}{"G"} * -1);
        my $ntt = norm($base_per_size{$i}{"-"}{"T"} * -1);

        print O4 $i . "\t$ta\t$tc\t$tg\t$tt\n";
        print O4 $i . "\t$nta\t$ntc\t$ntg\t$ntt\n";
    }
    close(O4);
}

if (defined($plot)) {
    if (defined($norm)) {
        `R --no-save $prefix.Geral_base_distribution_norm $prefix.NORMALIZED_Geral_base_distribution_by_reads RPM < $path_inner/calcDistributionPerBase_publication.R`;
        `R --no-save $prefix.Geral_base_distribution_norm $prefix.NORMALIZED_Geral_base_distribution_by_reads_ZOOM_piRNAS RPM < $path_inner/calcDistributionPerBase_publication_zoom_piRNAs.R`;
        `R --no-save $prefix.Geral_base_distribution_norm $prefix.NORMALIZED_Geral_base_distribution_sumstrands_by_reads RPM < $path_inner/calcDistributionPerBaseSumStrand_publication.R`;

    }

    `R --no-save $prefix.Geral_base_distribution $prefix.Geral_base_distribution_by_reads < $path_inner/calcDistributionPerBase_publication.R`;
    `R --no-save $prefix.Geral_base_distribution $prefix.Geral_base_distribution_by_reads_ZOOM_piRNAs < $path_inner/calcDistributionPerBase_publication_zoom_piRNAs.R`;
    `R --no-save $prefix.Geral_base_distribution $prefix.Geral_base_distribution_sumstrands_by_reads < $path_inner/calcDistributionPerBaseSumStrand_publication.R`;

}

if (defined($base)) {
    print "\n\nBase preference \n\n";

    foreach my $c (keys %ch2) {
        if ($total_reads{$c} > $minimo) {
            open(B, ">base_preference.$c.tab");

            my $sta = 0;
            my $stc = 0;
            my $stg = 0;
            my $stt = 0;

            if ($total_si{$c} != 0) {
                $sta = $base_si{$c}{"A"} * 100 / $total_si{$c};
                $stc = $base_si{$c}{"C"} * 100 / $total_si{$c};
                $stg = $base_si{$c}{"G"} * 100 / $total_si{$c};
                $stt = $base_si{$c}{"T"} * 100 / $total_si{$c};
            }

            my $pta = 0;
            my $ptc = 0;
            my $ptg = 0;
            my $ptt = 0;

            if ($total_pi{$c} != 0) {
                $pta = $base_pi{$c}{"A"} * 100 / $total_pi{$c};
                $ptc = $base_pi{$c}{"C"} * 100 / $total_pi{$c};
                $ptg = $base_pi{$c}{"G"} * 100 / $total_pi{$c};
                $ptt = $base_pi{$c}{"T"} * 100 / $total_pi{$c};
            }
			
            my $ta = 0;
            my $tc = 0;
            my $tg = 0;
            my $tt = 0;

            if ($total_reads{$c} != 0) {
                $ta = $base{$c}{"A"} * 100 / $total_reads{$c};
                $tc = $base{$c}{"C"} * 100 / $total_reads{$c};
                $tg = $base{$c}{"G"} * 100 / $total_reads{$c};
                $tt = $base{$c}{"T"} * 100 / $total_reads{$c};
            }

            print B "$sta\n$stc\n$stg\n$stt\n";
            print B "$pta\n$ptc\n$ptg\n$ptt\n";
            print B "$ta\n$tc\n$tg\n$tt";
            close(B);

        } else {
            print "\t\tDiscarting segment $c, ".$total_reads{$c}." <  $minimo \n";
        }

    }

}
print "Total number of mappings:\t$nmappings\nTotal number of reads:\t$nreads\n";