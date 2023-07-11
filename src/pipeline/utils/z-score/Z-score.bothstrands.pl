use strict;
use warnings;

use FindBin;
use Bio::SeqIO;  
use Bio::Seq::Quality;
use Getopt::Long;

use Statistics::Basic qw(mean stddev);
use Statistics::RankCorrelation;

#
# TODO: 2023-03-07 - Use a decent logger
# 
use constant PATH_LOG_MAIN => "small_rna_metavir.main.log";

#######################################################################
### PARSE INPUTS ------------------------------------------------------
#######################################################################

my $usage = "

$0  -sam <sam file> -p <prefix name>
$0 -h


-sam <sam file>			: Sam file from mapping on putative sequences
-fa 	<input file>        : Input fasta file of putative sequences
-h                      : Help message

";

$| = 1; # Forces immediate prints into files rather than the buffer.


my $ref;
my $search;
my $delimiter;
my $index;
my $profile;
my $fasta;
my $siRNA;
my $piRNA;
my $prefix;
my $help;

our %zscore;
our $z_stddev_p;
our $z_stddev_n;
our $z_stddev_t;
our $z_sumt=0;
our @zvp;
our @zvn;
our @zvt;

GetOptions (
	"sam=s" => \$search,
	"p=s" => \$prefix,
	"h!" => \$help
);

if ($help) {
	die $usage;
}

if (not(defined($search))) {
	die "\nGive a search file name", $usage;
}

if (not(defined($prefix))) {
	die "\nGive a prefix file name", $usage;
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
*STDERR = $LOG_FH;
# select $LOG_FH;

#######################################################################
### Main... -----------------------------------------------------------
#######################################################################

# 
# TODO: 2023-06-08 - This 2nd parameter doesn't seem to be used
# 
analyzeProfile($search, 0);

`touch $prefix.zscore.tab`;
`for f in .SCORE.*tab; do cat $prefix.zscore.tab | paste -d '#' - \$f > .temp; cp .temp $prefix.zscore.tab; rm -rf \$f; done; rm .temp`;

sub analyzeProfile {
	
	my $file = shift;
	my $compare = shift;
	
	my %scores;
	my %tamanhos;
	my %tamanhos_p;
	my %tamanhos_n;
	my $cp = 0;
	my $cn = 0;
	my %ch;
	my %ch2;

	print $LOG_FH "\nAnalyzing $file...\n";
	
	open(IN, "<$file");
	while (<IN>) {
		if (/^\w\S+\t\d+.+/) {
			
			my @campos = split(/\t/, $_);
			my $size = length($campos[9]);
			my $min;
			my $max;
			
			if ($size gt 14) {
				my $c = $campos[2];
				
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
				
				if (not defined($min) or $size < $min) {
					$min = $size;
				}
				if (not defined($max) or $size > $max) {
					$max = $size;
				}
				
				if ($campos[1] eq "0") {   
					$cp++;
					$ch{$c}{"pos"} = $ch{$c}{"pos"} + 1;
					$tamanhos_p{$size} = $tamanhos_p{$size} + 1;
					$ch2{$c}{"pos"}{"$size"} = $ch2{$c}{"pos"}{"$size"} + 1;

				} elsif ($campos[1] eq "16") {
					$cn++;
					$ch{$c}{"neg"} = $ch{$c}{"neg"} + 1; 
					$tamanhos_n{$size} = $tamanhos_n{$size} + 1;
					$ch2{$c}{"neg"}{$size} = $ch2{$c}{"neg"}{$size} + 1;
				}
			}
		}
	}
	
	my $t = $cp + $cn;

	my @pos;
	my @neg;
	my @tot;

	print $LOG_FH "\nCalculating Z-score...\n";

	foreach my $chr (keys %ch2) { # for each sequence
		
		my @neg;
		my @pos;
		my @tot;
		my @vn;
		my @vp;
		my @vt;
		my @zscore;
		
		my $sumt = 0;
		my $tot;
			
		# initializing array
		for (my $i = 0; $i <= 41; $i++) {
			$tot[$i] = 0;
		}
			
		# 
		# TODO: 2023-05-26 - Check what to do with these 'parallel' logging files
		# 
		open(O, ">.SCORE.$chr.tab");
		print O "$chr\n";
		
		for (my $i = 0; $i <= 41; $i++) {
			
			my $p = 0;
			my $n = 0;
			
			if (exists $ch2{$chr}{"pos"}{$i}) {
				$p = $ch2{$chr}{"pos"}{$i};
			}
			if (exists $ch2{$chr}{"neg"}{$i}) {
				$n = $ch2{$chr}{"neg"}{$i};
			}
			
			my $t = $p + $n;
			$tot[$i - 15] = $p;
			$tot[$i + 6] = $n;
		}

		# z = (x - mean)/sd
		my $stddev_t = stddev(@tot);
		my $mean_t = mean(@tot);

		my @zp;
		my @zn;
		my @zt;
		for (my $i = 0; $i <= 41; $i++) {
			if ($stddev_t != 0 && $mean_t != 0) {
				$zscore{$i}{"both"} = ($tot[$i] - $mean_t) / $stddev_t;
				push(@zvt, ($tot[$i] - $mean_t) / $stddev_t);
			}
		}
		
		if ($stddev_t != 0 && $mean_t != 0) {
			for (my $i = 0; $i <= 41; $i++) {
				print O $zscore{$i}{"both"} . "\n";
			}
		}

		undef(@zscore);
		undef @pos;
		undef @neg;
		undef @tot;
		undef @vp;
		undef @vn;
		undef @vt;
		close(O);
	}
}
