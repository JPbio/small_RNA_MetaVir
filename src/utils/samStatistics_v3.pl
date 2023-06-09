# use strict;
use warnings;

use Bio::SeqIO;
use Bio::Seq::Quality;
use Getopt::Long;

#
# TODO: 2023-05-15 - Use a decent logger
# 
use constant PATH_LOG_MAIN => "srna_metavir.main.log";

#######################################################################
### PARSE INPUTS ------------------------------------------------------
#######################################################################

my $usage = "
$0 -sam <sam file> -fa <fasta file> -p <prefix name> [--profile] [--counts]
$0 -h

-i <input file>	: Input file with columns delimited by a specific delimiter
-h											: Help message

#
# calc stattistics by reads
#
";

$| = 1;		# forces immediate prints into files rather than the buffer.

my $sam;
my $delimiter;
my $index;
my $profile;
my $fasta;
my $counts;
my $prefix;
my $help;

GetOptions(
	"sam=s" => \$sam,
	"fa=s" => \$fasta,
	"p=s" => \$prefix,
	"profile!" => \$profile,
	"counts!" => \$counts,
	"h!" => \$help
);

if ($help) {
	die $usage;
}

if (not(defined($sam))) {
	die "\nGive an input file name", $usage;
}
if (not(defined($prefix))) {
	die "\nGive a prefix name", $usage;
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

my %tamanhos;
my %tamanhos_p;
my %tamanhos_n;
my $cp = 0;
my $cn = 0;
my %ch;
my %ch2;

my $min = 1000;
my $max = 0;

my %contig;

if (defined($fasta)) {

	print "[sam statistics] Loading FASTA file...\n";
	
    my $in = Bio::SeqIO->new(-format => 'fasta', -file => $fasta);
	while (my $data = $in->next_seq) {
		my $size = $data->length;
		$contig{$data->primary_id} = $data->seq;
		$contig{$data->primary_id} = $data->seq;
		$contig{$data->primary_id}{"desc"} = $data->desc;
		#print " ---> ".$contig{$data->primary_id}{"desc"} ."\n";
	}
}

my %count_piRNAs;

my %reads;
while (<IN>) {
    if (/^\w+/) {
        
		my @campos = split(/\t/,$_);
        my $size = length($campos[9]);
        my $c = $campos[2];
        
		if ($size >= 15 and $size<= 35) {
            if (not exists($reads{$campos[2]}{$campos[0]})) {
                
				$reads{$campos[2]}{$campos[0]} = 1;
                
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
                
				if ($size < $min) {
					$min = $size;
				}
                
				if ($size > $max) {
					$max = $size;
				}
                
				if ($campos[1] eq "0") {   
                    $cp++;
                    $ch{$c}{"pos"} = $ch{$c}{"pos"} + 1;
                    $tamanhos_p{$size} = $tamanhos_p{$size} + 1;
                    $ch2{$c}{"pos"}{"$size"} = $ch2{$c}{"pos"}{"$size"} +1;
                    #print "POSITIVE in [$c] fita [pos]  size [$size] amount: ".$ch{$c}{"pos"}{$size} ."\n";

                } elsif ($campos[1] eq "16") {
                    $cn++;
                    $ch{$c}{"neg"} = $ch{$c}{"neg"} + 1; 
                    $tamanhos_n{$size} = $tamanhos_n{$size} + 1;
                    $ch2{$c}{"neg"}{$size} = $ch2{$c}{"neg"}{$size} + 1;
                    #print "NEGATIVE in [$c] fita [pos]  size [$size] amount: ".$ch{$c}{"neg"}{$size} ."\n";
                }
            }
        }
    }
}
my $t = $cp + $cn;

# print "Total in both strand: $t\nTotal in positve strand: $cp\nTotal in negative strand: $cn\n";
# print "\nSize\tSENSE\tANTISENSE\tTOTAL\n";
# foreach (sort { $a <=> $b } keys(%tamanhos_p)) {
# 	$s = $tamanhos_p{$_}  + $tamanhos_n{$_};
# 	if ($s ne "0") {
# 		print "$_\t".$tamanhos_p{$_}."\t".$tamanhos_n{$_}."\t".$s. "\n";
# 	}
# }

my %signal;
foreach my $c (keys %ch2) {
    
	$signal{$c}{"si"} = "no";
    $signal{$c}{"pi"} = "no";
    $signal{$c}{"primary_piRNA"} = "no";

    my $scorePI = 0;
    
    my $sumT1P = 0;
    my $sumT1N = 0;
    my $sumT2P = 0;
    my $sumT2N = 0;
    
    my $allT1 = 0;
    my $allT2 = 0;
    my $allT3 = 0;
    my $allT4 = 0;
    my $allT4_pos = 0;
    my $allT4_neg = 0;

    for (my $s = $min; $s <= $max; $s++) {
        
        my $pa = $ch2{$c}{"pos"}{$s} // 0;
        my $na = $ch2{$c}{"neg"}{$s} // 0;

        if ($pa <= 0) {
			$pa = 0;
		}
        
		if ($na <= 0) {
			$na = 0;
		}

        my $ta = $pa + $na;

        if ($s >= 15 and $s <= 19) {
            $sumT1P += $pa;
            # $sumT1N += $pn;
            $allT1 = $sumT1P + $sumT1N;

        } elsif ($s >= 20 and $s <= 22) {
            
			$allT4 = $allT4 + $pa + $na;
            $allT4_pos = $allT4_pos + $pa;
            $allT4_neg = $allT4_neg + $na;

            my $r;
            if ($s == 22) {
                if (($allT4_pos > 50 or $allT4_neg > 50) and ($allT4 > $allT1)) {
                    
					if ($na == 0) {
						$na = 1;
					}
                    
					if ($pa == 0) {
						$pa = 1;
					}
                    
					$r = $pa / $na;
                    $r = $na / $pa if $na > $pa;

                    if ($r > 0.3 and $r < 10 and $r != 1) {
                        $signal{$c}{"si"} = "yes";
                    }
                }
            }

            if ($pa > 0 or $na > 0) {
                if (($pa > 50 or $na > 50) and ($ta > $allT1)) {
                    if ($pa > 0 and $na > 0) {
                        $r = $pa / $na;
                        $r = $na / $pa if $na > $pa;
                    } else {
                        $r = 0;
                    }

                    if ($r > 0.3 and $r < 10) {
                        $signal{$c}{"si"} = "yes";
                    }
                }

                $allT4 = $allT4 + $pa + $na;
            }

        } elsif ($s == 22 or $s == 23) {
            $sumT2P += $pa;
            $sumT2N += $na;
            $allT2 = $sumT2P + $sumT2N;

        } elsif ($s >= 24 and $s <= 29) {
            if ($pa > 0 or $na > 0) {

                my $r2;

                if ($pa > 0 and $na > 0) {
                    $r2 = $pa / $na;
                    $r2 = $na / $pa if $na > $pa;

                    if ($r2 > 1000) {
                        $signal{$c}{"primary_piRNA"} = "yes";
                    }

                } else {
                    $r2 = 0;
                    $signal{$c}{"primary_piRNA"} = "yes";
                }

                if (($pa > 50 and $na > 50)) {
                    $scorePI += 1;
                }

                $allT3 += $ta;
            }

        } elsif ($s == 30) {
            if (($allT3 > $allT2 + $allT1) and ($scorePI >= 4)) {
                $signal{$c}{"pi"} = "yes";
            }
        }

        #print $s."\t$pa\t$na\t$ta\n";
    }

    if (($allT4 > $allT3) && ($allT4 > $allT1)) {
        $signal{$c}{"final"} = "siRNA";
    }

    if (($allT3 > $allT3) && ($allT4 > $allT1)) {
        if (exists($signal{$c}{"primary_piRNA"})) {
            $signal{$c}{"final"} = "primary_piRNA";

        } else {
            $signal{$c}{"final"} = "piRNA";
        }
    }
    
	if (($allT1 > $allT3) && ($allT1 > $allT4)) {
        $signal{$c}{"final"} = "degradation";
    }
}

#print " \n###############SIGNAL of siRNA or piRNA #############\n\n";

open(SIG, ">$prefix.signature");
print SIG "Contig\tsiRNA\tpiRNA\tFINAL\n";

foreach my $c (keys %signal) {

    my $si = $signal{$c}{"si"} // "";
    my $pi = $signal{$c}{"pi"} // "";
    my $final = $signal{$c}{"final"} // "";
    my $primaryDna = $signal{$c}{"primary_piRNA"} // "";

    print SIG "$c\t" . $si . "\t" . $pi . "\t" . $primaryDna . "\t" . $final . " \n";
}

close(SIG);

if (defined($profile)) {

	open(S, ">$prefix.profile");
    open(SI, ">$prefix.withSiRNA.fasta");
    open(PI, ">$prefix.withPiRNA.fasta");
    open(primaryPI, ">$prefix.with_primaryPiRNA.fasta");
    open(SIPI, ">$prefix.withSiRNA_and_PiRNA.fasta");

    foreach my $c (keys %signal) {
        print S " $c\t" . $signal{$c}{"si"} . " \t" . $signal{$c}{"pi"} . " \n";
        
		if ($signal{$c}{"si"} eq "yes" and $signal{$c}{"pi"} eq "no") {
            print SI ">$c " . $contig{$c}{"desc"} . "\n" . $contig{$c} . "\n";
        }
        
		if ($signal{$c}{"pi"} eq "yes" and $signal{$c}{"si"} eq "no" and $signal{$c}{"primary_piRNA"} eq "no") {
            print PI ">$c " . $contig{$c}{"desc"} . "\n" . $contig{$c} . "\n";
        }
        
		if ($signal{$c}{"primary_piRNA"} eq "yes" and $signal{$c}{"si"} eq "no" and $signal{$c}{"pi"} eq "no") {
            print primaryPI ">$c " . $contig{$c}{"desc"} . "\n" . $contig{$c} . "\n";
        }
        if ($signal{$c}{"pi"} eq "yes" and $signal{$c}{"si"} eq "yes") {
            print SIPI ">$c " . $contig{$c}{"desc"} . "\n" . $contig{$c} . "\n";
        }
    }

    close(S);
    close(SI);
    close(PI);
    close(primaryPI);
    close(SIPI);
}

# print "\n\nReads per chromosoes\n";

open(CHR, ">$prefix.readsPerChr");

if ($min > 5) {
    $min = 5;
}
if ($max < 30) {
    $max = 30;
}

foreach my $c (keys %ch2) {

    for (my $s = $min; $s <= $max; $s++) {
        
		my $pa = $ch2{$c}{"pos"}{$s} // 0;
        my $na = $ch2{$c}{"neg"}{$s} // 0;

        if ($pa <= 0) {
            $pa = 0;
        }
        
		if ($na <= 0) {
            $na = 0;
        }

        my $ta = $pa + $na;
        print CHR $c . "\t" . $s . "\t$pa\t$na\t$ta\n";
    }
}

close(CHR);
