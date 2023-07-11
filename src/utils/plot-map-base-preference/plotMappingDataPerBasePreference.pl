use strict;
use warnings;

use Bio::SeqIO;  
use Bio::Seq::Quality;
use Getopt::Long;

#
# TODO: 2023-05-15 - Use a decent logger
# 
use constant PATH_LOG_MAIN => "small_rna_metavir.main.log";

#######################################################################
### PARSE INPUTS ------------------------------------------------------
#######################################################################

my $usage = "

$0 -sam <sam file> -s <start> -e <end> -fa <fasta> -pace <num> -p prefix --profile --keep --pattern
$0 -h

-i <input file>         : Input file with columns delimited by a specific delimiter
-h                      : Help message

";

$| = 1;     # forces immediate prints into files rather than the buffer.


my $sam;
my $profile;
my $pattern;
my $prefix;
my $start;
my $pace;
my $norm;
my $end;
my $fasta;

# 
# TODO: 2023-05-30 - Make it a 'my' variable
# 

our $minimo;
my $keep;
my $antisense;
my $help;

GetOptions(
    "sam=s" => \$sam,
    "pattern!" => \$pattern,
    "profile!" => \$profile,
    "p=s" => \$prefix,
    "s=i" => \$start,
    "pace=i" => \$pace,
    "norm=i" => \$norm,
    "e=i" => \$end,
    "fa=s" => \$fasta,
    "m=i" => \$minimo,
    "keep!" => \$keep,
    "antisense!" => \$antisense,
    "h!" => \$help,
);

if ($help) {
    die $usage;
}

if (not(defined($sam))) {
    die "\nGive a sam input file name", $usage;
}
if (not(defined($fasta))) {
    die "\nGive a fasta input file name", $usage;
}
if (not(defined($prefix))) {
    die "\nGive an prefix file name", $usage;
}

if (not(defined($pace))) {
    $pace = 1;
}

if (not defined($start) or $start < 0) {
    warn("\n\n\t\t[Bad value for start (".($start // "undef")."), start setted to 15]\n");
    $start = 15;
}

if (not defined($end) or ($end < 0) or ($end < $start)) {
    warn("\n\t\t[Bad value for end (".($end // "undef")."), end setted to 30]\n");
    $end = 30;
}

if (not defined($minimo) or ($minimo < 0)) {
    warn("\n\t\t[Bad value for min (".($minimo // "undef")."), min setted to 30]\n");
    $minimo = 30;
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

*STDERR = $LOG_FH;

#######################################################################
### Main... -----------------------------------------------------------
#######################################################################

my $path_utils = "/small-rna-metavir/src/utils";
my $path_inner = "$path_utils/plot-map-base-preference/inner";

my $path_filter_sam_bowtie = "$path_inner/filterSamBowtieBySize.pl";
my $path_split_chromosome = "$path_inner/splitSamByChromosome.pl";
my $path_calc_density_base = "$path_inner/calcDensityPerBase.pl";
my $path_plot_dist = "$path_inner/plotDistribution_ggplot2.R";
my $path_plot_density = "$path_inner/plotDensity_ggplot2.R";
my $path_calc_dist_base = "$path_inner/calcDistributionPerBase.R";
my $path_calc_dist_base_publication = "$path_inner/calcDistributionPerBase_publication2.R";

open(IN, "<$sam");

# 
# REVIEW: 2023-05-30 - Should all these 'our' variables be 'my'?
# 

our %tamanhos;
our %tamanhos_p;
our %tamanhos_n;
# our $cp = 0;
# our $cn = 0;
our %ch;
our %ch2;
our %base;
our %base_total;
our %base_per_size;

my $min = 1000;
my $max = 0;

for (my $i = 15; $i <= 35; $i++) {
    $base_per_size{$i}{"+"}{"A"} = 0;
    $base_per_size{$i}{"+"}{"C"} = 0;
    $base_per_size{$i}{"+"}{"G"} = 0;
    $base_per_size{$i}{"+"}{"T"} = 0;
    $base_per_size{$i}{"-"}{"A"} = 0;
    $base_per_size{$i}{"-"}{"C"} = 0;
    $base_per_size{$i}{"-"}{"G"} = 0;
    $base_per_size{$i}{"-"}{"T"} = 0;
}

my $spos = "0";
my $sneg = "16";

if (defined($antisense)) {
    $spos = "16";
    $sneg = "0";
}

print $LOG_FH "[plot mapping data per base preference] Loading SAM file...\n";

while (<IN>) {

    if (/^\w+/) {
        
        my @campos = split(/\t/, $_);
        my $size = length($campos[9]);
        my $c = $campos[2];

        if (not exists $tamanhos_p{$size}) {
            
            $tamanhos_p{$size} = 0;
            
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
            $ch2{$c}{"pos"}{$size} = 0;
            $base{$c}{"pos"}{$size}{"A"} = 0;
            $base{$c}{"pos"}{$size}{"C"} = 0;
            $base{$c}{"pos"}{$size}{"G"} = 0;
            $base{$c}{"pos"}{$size}{"T"} = 0;
        }
        
        if (not exists $ch2{$c}{"neg"}{$size}) {
            $ch2{$c}{"neg"}{$size} = 0;
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

        # Count miRNA preference
        if ($size >= 20 && $size <= 24) {
            $base_total{$c}{$first} = ($base_total{$c}{$first} // 0) + 1;
        }
        
        if ($campos[1] eq $spos) {
            # $cp++;
                
            $ch{$c}{"pos"} = ($ch{$c}{"pos"} // 0) + 1;
            $tamanhos_p{$size} = ($tamanhos_p{$size} // 0) + 1;
            $ch2{$c}{"pos"}{"$size"} = ($ch2{$c}{"pos"}{"$size"} // 0) + 1;
            $base{$c}{"pos"}{$size}{$first} = ($base{$c}{"pos"}{$size}{$first} // 0) + 1;
            
            # print $LOG_FH "POSITIVE in [$c] fita [pos]  size [$size] amount: ".$ch{$c}{"pos"}{$size} ."\n";
            $base_per_size{$size}{"+"}{$first} = $base_per_size{$size}{"+"}{$first} + 1;

        } elsif ($campos[1] eq $sneg) {
            # $cn++;
        
            $ch{$c}{"neg"} = ($ch{$c}{"neg"} // 0) + 1;
            $tamanhos_n{$size} = ($tamanhos_n{$size} // 0) + 1;
            $ch2{$c}{"neg"}{$size} = ($ch2{$c}{"neg"}{$size} // 0) + 1;
            $base{$c}{"neg"}{$size}{$first} = ($base{$c}{"neg"}{$size}{$first} // 0) + 1;
            
            #print $LOG_FH "NEGATIVE in [$c] fita [pos]  size [$size] amount: ".$ch{$c}{"neg"}{$size} ."\n";
            $base_per_size{$size}{"-"}{$first} = ($base_per_size{$size}{"-"}{$first} // 0) + 1;
        }
    }
}
# $t = $cp + $cn;

# Load fasta file
our %contig;

print $LOG_FH "[plot mapping data per base preference] Loading FASTA file...\n";
my $in = Bio::SeqIO->new(-format => 'fasta', -file => $fasta);

while (my $data = $in->next_seq) {
    my $size = $data->length;
	$contig{$data->primary_id} = $size;
}

# Writing file with total mappings
if (not defined($profile)) {
    
    print $LOG_FH "\n\n[Total read distribution]\n";
    
    for (my $i = $start; $i <= $end; $i++) {
        
        my $p = 0;
        my $n = 0; 
        
        if (exists ($tamanhos_p{$i})) {
            $p = $tamanhos_p{$i} ;
        }
        if (exists ($tamanhos_n{$i})) {
            $n = $tamanhos_n{$i} ;
        }

        my $t = $p + $n;
        print $LOG_FH "$i\t$p\t$n\t$t\n";
    }
}

###################

# Print total of reads per chromosome
if (defined($profile)) {

    print $LOG_FH "\n\nReads per chromosoes\n";

    open(O, ">$prefix.readsPerChromossome");

    foreach my $c (keys %ch) {
        
        my $pa = $ch{$c}{"pos"};
        my $na = $ch{$c}{"neg"};
        my $ta = $pa + $na;

        if ($pa <= 0) {
            $pa = 0;
        }
        if ($na <= 0) {
            $na = 0;
        }

        print O $c."\t". $ta . "\t" . $na .  "\t". $ta ."\n";
    }

    close(O);

    # `R --no-save --args $prefix.$c.distribution < file.R`;

    ##############
    
    print $LOG_FH "\n\nReads per chromosoes\n";

    foreach my $c (keys %ch2) {
        
        open(O, ">$prefix.$c.distribution");
        open(O2, ">$prefix.$c.base_distribution");

        print O2 "\tA\tC\tG\tT\n";

        # print $LOG_FH "\n\nChromossome [$c]\n";
        
        for (my $s = $start; $s <= $end; $s++) {
            
            my $pa = $ch2{$c}{"pos"}{$s} // 0;
            my $na = $ch2{$c}{"neg"}{$s} // 0;
            
            my $bpa = $base{$c}{"pos"}{$s}{"A"} // 0;
            my $bpc = $base{$c}{"pos"}{$s}{"C"} // 0;
            my $bpg = $base{$c}{"pos"}{$s}{"G"} // 0;
            my $bpt = $base{$c}{"pos"}{$s}{"T"} // 0;
            
            my $bna = ($base{$c}{"neg"}{$s}{"A"} // 0) * -1;
            my $bnc = ($base{$c}{"neg"}{$s}{"C"} // 0) * -1;
            my $bng = ($base{$c}{"neg"}{$s}{"G"} // 0) * -1;
            my $bnt = ($base{$c}{"neg"}{$s}{"T"} // 0) * -1;

            if ($pa <= 0) {
                $pa = 0;
            }
            if ($na <= 0) {
                $na = 0;
            }

            my $ta = $pa + $na;
            print O $s."\t$pa\t$na\t$ta\n";
            print O2 $s."\t$bpa\t$bpc\t$bpg\t$bpt\n";
            print O2 $s."\t$bna\t$bnc\t$bng\t$bnt\n";

        }
        close(O);
        # drawing distribution plot
        #`R --no-save --args $prefix.$c.distribution < file.R`;
    }

    # Creating sam files
    my $split = "perl $path_split_chromosome -in $sam -prefix $prefix";
    `$split`;

    # CalcDensityPerBase
    foreach my $c (keys %contig) {

        # my $refsize = $contig{$c};
        my $cp = "$prefix.$c";
        my $samFile = quotemeta("$cp.sam");
        
        # my $countingCmd = "grep -v \"@\" ". quotemeta($samFile) ." | wc -l  | cut -f 1 -d ' '";
        my $countingCmd = "grep -v \"@\" $samFile | wc -l  | cut -f 1 -d ' '";
        my $counting = `$countingCmd`;
        $counting =~ s/\n//g;

        if ($counting >= $minimo) {
            
            print $LOG_FH "####### $cp was USED, $counting > $minimo \n";

            my $f0 = "perl $path_filter_sam_bowtie -i $samFile -si 20 -se 23 -o $cp.20to23.sam";
            my $f1 = "perl $path_filter_sam_bowtie -i $samFile -si 21 -se 21 -o $cp.21.sam";
            my $f2 = "perl $path_filter_sam_bowtie -i $samFile -si 24 -se 30 -o $cp.24to29.sam";
            my $f3 = "perl $path_filter_sam_bowtie -i $samFile -si 22 -se 22 -o $cp.22to22.sam";
            my $f4 = "perl $path_filter_sam_bowtie -i $samFile -si 15 -se 19 -o $cp.15to19.sam";

            `$f0`;
            `$f1`;
            `$f2`;
            `$f3`;
            `$f4`;

            # print $LOG_FH "$f1\n$f2\n";
            my $comm0 = "perl $path_calc_density_base -sam $cp.20to23.sam -pace $pace -ref ".$contig{$c}."  -out  $cp.20to23.density";
            my $comm = "perl $path_calc_density_base -sam $samFile -pace $pace -ref ".$contig{$c}."  -out  $cp.density";
            my $comm2 = "perl $path_calc_density_base -sam $cp.21.sam -pace $pace -ref ".$contig{$c}."  -out $cp.21.density";
            my $comm3 = "perl $path_calc_density_base -sam $cp.24to29.sam -pace $pace -ref ".$contig{$c}."  -out $cp.24to29.density";
            my $comm4 = "perl $path_calc_density_base -sam $cp.22to22.sam -pace $pace -ref ".$contig{$c}."  -out $cp.22to22.density";
            my $comm5 = "perl $path_calc_density_base -sam $cp.15to19.sam -pace $pace -ref ".$contig{$c}."  -out $cp.15to19.density";

            # print $LOG_FH "$comm \n $comm2 \n $comm3 \t $comm0 \n";
            `$comm0`;
            `$comm`;
            `$comm2`;
            `$comm3`;
            `$comm4`;
            `$comm5`;

            # print $LOG_FH "Plotting graphs [$c]...\n";
            # `/opt/R/4.1.2/bin/R --no-save $prefix.$c.distribution $cp.distribution < $path_plot_dist 2>/dev/null `;
            # `/opt/R/4.1.2/bin/R --no-save $cp.density $cp.all.density < $path_plot_density 2>/dev/null`;
            # `/opt/R/4.1.2/bin/R --no-save $cp.21.density $cp.21.density < $path_plot_density 2>/dev/null`;       
            # `/opt/R/4.1.2/bin/R --no-save $cp.24to29.density $cp.24to29.density < $path_plot_density 2>/dev/null`;

            `R --no-save $prefix.$c.distribution $cp.distribution < $path_plot_dist  2>/dev/null`;
            `R --no-save $cp.density $cp.all.density < $path_plot_density 2>/dev/null`;
            `R --no-save $cp.21.density $cp.21.density < $path_plot_density 2>/dev/null`;
            `R --no-save $cp.24to29.density $cp.24to29.density < $path_plot_density 2>/dev/null`;
            `R --no-save $cp.20to23.density $cp.20to23.density < $path_plot_density 2>/dev/null`;
            `R --no-save $cp.22to22.density $cp.22to22.density < $path_plot_density 2>/dev/null`;
            `R --no-save $cp.15to19.density $cp.15to19.density < $path_plot_density 2>/dev/null`;
            `R --no-save $prefix.$c.base_distribution $prefix.$c.base_distribution < $path_calc_dist_base 2>/dev/null`;
            `R --no-save $prefix.$c.base_distribution $prefix.$c.publication < $path_calc_dist_base_publication 2>/dev/null`;

        } else {
            print $LOG_FH "--> $cp was [DISCARDED], $counting < $minimo <--\n";
        }
        
        # `rm -rf  $cp.21.sam $cp.24to29.sam $cp.sam`;
        if (not defined($keep)) {
            `rm -rf $cp.20to23.density $cp.15to19.density $prefix.$c.distribution $cp.density $cp.21.density $cp.24to29.density $cp.22to22.density`;
        }

        # `rm -rf $cp.density $cp.21.density $cp.24to29.density $cp.21.sam $cp.24to29.sam $cp.sam`;

        # $exec = "R --no-save --args $prefix.$c.distribution $cp.density $cp.21.density $cp.24to29.density $refsize $c < /Users/eric/Desktop/projetos/scripts_servers/R/plotMappingDensity.R ";
        # `$exec`;

        # print $LOG_FH "$comm \n";
    }
    
    # calcDensityPerBase.pl -sam -pace -ref ";
}

# miRNA preference print
open(O3, ">$prefix.miRNA_total_base_distribution");

print O3 "\tA\tC\tG\tT\n";

# ^[[CUse of uninitialized value $tc in concatenation (.) or string at /small-rna-metavir/src/utils/plot-map-base-preference/plotMappingDataPerBasePreference.pl line 481, <GEN0> line 25.


foreach my $c (keys %ch2) {
    my $ta = $base_total{$c}{"A"} // 0;
	my $tc = $base_total{$c}{"C"} // 0;
	my $tg = $base_total{$c}{"G"} // 0;
	my $tt = $base_total{$c}{"T"} // 0;
    print O3 $c."\t$ta\t$tc\t$tg\t$tt\n";
}
close(O3);

# GERAL preference print
open(O4, ">$prefix.Geral_base_distribution");

print O4 "size#A#C#G#T\n";

for (my $i = 15; $i <= 35; $i++) {

    my $ta = $base_per_size{$i}{"+"}{"A"};
	my $tc = $base_per_size{$i}{"+"}{"C"};
	my $tg = $base_per_size{$i}{"+"}{"G"};
	my $tt = $base_per_size{$i}{"+"}{"T"};
	
	my $nta = $base_per_size{$i}{"-"}{"A"} * -1;
	my $ntc = $base_per_size{$i}{"-"}{"C"} * -1;
	my $ntg = $base_per_size{$i}{"-"}{"G"} * -1;
	my $ntt = $base_per_size{$i}{"-"}{"T"} * -1;
	
    print O4 $i."#$ta#$tc#$tg#$tt\n";
    print O4 $i."#$nta#$ntc#$ntg#$ntt\n";
}
close(O4);

print $LOG_FH "Plotting graphs...\n";

#`montage -mode concatenate -tile 3x $prefix*density.png $prefix.merge.pdf`;
#`convert $prefix*density.png +append $prefix.merge.2.pdf`;
#`convert $prefix*density.png -append $prefix.merge.3.pdf`;
       
#`rm -rf *.png`;
#`rm -rf $prefix*.distribution  $prefix*.sam $prefix*.density`;
# $exec="R --no-save --args $prefix.readsPerChromossome < /Users/eric/Desktop/projetos/scripts_servers/R/plotMappingDistribution.R ";
#`$exec`;
