#use strict;

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
### Parse inputs ------------------------------------------------------
#######################################################################

my $usage = "

$0  -sam <sam file> -p <prefix name>
$0 -h


-sam <sam file>			: Sam file from mapping on putative sequences
-fa 	<input file>        : Input fasta file of putative sequences
-h                      : Help message

";

$| = 1;     # forces immediate prints into files rather than the buffer.


my $ref;
my $search;
my $delimiter;
my $index;
my $profile;
my $fasta;
my $siRNA;
my $piRNA;
my $prefix;

our %zscore;
our $z_stddev_p ;
our $z_stddev_n ;
our $z_stddev_t ;
our $z_sumt=0;
our @zvp, @zvn,@zvt;
our $fasta;


GetOptions ("sam=s" => \$search,"fa=s" => \$fasta,"p=s" => \$prefix,"h!" => \$help);

if ($help) {
        die $usage;
}

if (not(defined($search))) {
        die "\nGive an search file name",$usage;
}

if (not(defined($prefix))) {
        die "\nGive an prefix file name",$usage;
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
### Main --------------------------------------------------------------
#######################################################################

#Loading fasta files
if(defined($fasta)){
	our %contig;

	warn("Loading FASTA file...\n");
	my $in = Bio::SeqIO->new(-format    => 'fasta',
							   -file      => $fasta);
	while (my $data = $in->next_seq) {

		$contig{$data->primary_id}=$data->seq;
		$contig{$data->primary_id}{"desc"}=$data->desc;
		$contig{$data->primary_id}{"length"}=$data->length;


	}


}


analyzeProfile($search,0);

`touch $prefix.zscore.tab`;
`for f in .SCORE.*tab; do cat $prefix.zscore.tab | paste -d '#' - \$f >.temp; cp .temp $prefix.zscore.tab; rm -rf \$f; done; rm .temp`;
 

#Calculating Z-score both strands
sub analyzeProfile{
	$file =shift;
	$compare=shift;
	my %scores;
	my %tamanhos;
	my %tamanhos_p;
	my %tamanhos_n;
	my $cp=0;
	my $cn=0;
	my %ch;
	my %ch2;
	my %countings;

	print $LOG_FH "Analyzing $file...\n";
	open(IN,"<$file");
	while (<IN>){
#		print "entrou.... ($_)\n";
		#if(/^\w+\t\d+\t\w+\t\d+\t\S+\t\S+\t\S+\t\S+\t\S+.+/){
		if(/^\w\S+\t\d+.+/){
#			print ".... ($_)\n";
			@campos = split(/\t/,$_);
			$size = length($campos[9]);
			
			if($size gt 14){
			
				$c = $campos[2];	
	
				  if (not exists $countings{$c}{"18to35"} ){
							$countings{$c}{'21'} = 0;
							$countings{$c}{"24to30"} = 0;
							$countings{$c}{"20to22"} = 0;
							$countings{$c}{"25to29"} = 0;
							$countings{$c}{"20to23"} = 0;
							$countings{$c}{"18to35"} = 0;
							$countings{$c}{"15to18"} = 0;
							
				  }
				 if (not exists $tamanhos_p{$size} ){
							$tamanhos_p{$size} = 0;
					}
				   if (not exists $ch{$c}{"pos"} ){
							$ch{$c}{"pos"} = 0;
					}
			
					if (not exists $ch{$c}{"neg"} ){
							$ch{$c}{"neg"} = 0;
					}
	
				 if (not exists $tamanhos_n{$size} ){
							$tamanhos_n{$size} = 0;
					}
		
				if(not exists $ch2{$c}{"pos"}{$size}){
					$ch2{$c}{"pos"}{$size}=0;
				}
		
				if(not exists $ch2{$c}{"neg"}{$size}){
					$ch2{$c}{"neg"}{$size}=0;
				}
		
				if($size < $min) { $min=$size;}
				if($size > $max) { $max=$size;}
		
				#computing number of reads (21) (24-30) (18-35)
				#print "analyzing ($c) read ($size)  - 21: ".$countings{$c}{'21'}." \n";
				if (exists $countings{$c} ){
						if($size == 21 ){
							$countings{$c}{'21'}++;
				#			print "entrei 21 ($c - $size)".$countings{$c}{'21'}."\n";
						}
						if($size > 23 and $size < 31 ){
							$countings{$c}{"24to30"} = $countings{$c}{"24to30"} + 1;
				#			print "entrei 24-30 ($c - $size)".$countings{$c}{"24to30"}."\n";
						}
						if($size > 24 and $size < 30 ){
							$countings{$c}{"25to29"} = $countings{$c}{"25to29"} + 1;
				#			print "entrei 24-30 ($c - $size)".$countings{$c}{"24to30"}."\n";
						}
						if($size > 19 and $size < 24 ){
							$countings{$c}{"20to23"} = $countings{$c}{"20to23"} + 1;
				#			print "entrei 20-23 ($c - $size)".$countings{$c}{"20to23"}."\n";
						}
						if($size > 19 and $size < 23 ){
							$countings{$c}{"20to22"} = $countings{$c}{"20to22"} + 1;
				#			print "entrei 20-23 ($c - $size)".$countings{$c}{"20to23"}."\n";
						}
						if($size > 17 and $size < 36 ){
							$countings{$c}{"18to35"} = $countings{$c}{"18to35"} + 1;
				#			print "entrei 18-35 ($c - $size)".$countings{$c}{"18to35"}."\n";
						}
						if($size > 14 and $size < 19 ){
							$countings{$c}{"15to18"} = $countings{$c}{"15to18"} + 1;
				#			print "entrei 18-35 ($c - $size)".$countings{$c}{"18to35"}."\n";
						}
							
				}
				
		
				if ($campos[1] eq "0"){   
					$cp++;
					$ch{$c}{"pos"} = $ch{$c}{"pos"} + 1;
					$tamanhos_p{$size}= $tamanhos_p{$size} + 1;
					$ch2{$c}{"pos"}{"$size"}= $ch2{$c}{"pos"}{"$size"} +1;
	
				}elsif($campos[1] eq "16"){
					$cn++;
					$ch{$c}{"neg"} = $ch{$c}{"neg"} + 1; 
					$tamanhos_n{$size}= $tamanhos_n{$size} + 1;
					$ch2{$c}{"neg"}{$size}= $ch2{$c}{"neg"}{$size} + 1;
					#print "NEGATIVE in [$c] fita [pos]  size [$size] amount: ".$ch{$c}{"neg"}{$size} ."\n";
			
				}
			}
		}

	}
	
	$t = $cp+$cn;
	
	my @pos;
	my @neg;
	my @tot;

	
	
	print $LOG_FH "Calculating Z-score...\n";


    foreach $chr (keys %ch2){   #for each sequence
            my @pos;
            my @neg;
            my $tot;
            my $sumt=0;
            
            #initializing array
            for (my $i=0;$i<=41;$i++){
	        	$tot[$i]=0;	
        	}
            
            open(O,">.SCORE.$chr.tab");
            print O "$chr\n";
            for (my $i=0;$i<=41;$i++){
                    $p=0;
                    $n=0;
                    if(exists $ch2{$chr}{"pos"}{$i}){ $p = $ch2{$chr}{"pos"}{$i}; }
                    if(exists $ch2{$chr}{"neg"}{$i}){ $n = $ch2{$chr}{"neg"}{$i}; }
                   
    				$t= $p+$n;
                 	$tot[$i-15]=$p;
                    $tot[$i+6]=$n;
                    #print  "$i\t$p\t$n\t$t\n";
            }

        ##z = (x - mean)/sd

        my $stddev_t   = stddev(@tot);
        my $mean_t   = mean(@tot);

        my @zp;
        my @zn;
        my @zt;
        for (my $i=0;$i<=41;$i++){
        	if($stddev_t != 0 && $mean_t  !=0){

            	$zscore{$i}{"both"}     = ($tot[$i] - $mean_t) / $stddev_t;
	            push(@zvt,($tot[$i] - $mean_t) / $stddev_t);
            }

        }
		if($stddev_t != 0 && $mean_t  !=0){
      	  for (my $i=0;$i<=41;$i++){
      	          print O $zscore{$i}{"both"}."\n";
       	 }
       	}



		###computing densities
    	#print O "(".$countings{$chr}{"21"}.")\n";
		#print O "(".$countings{$chr}{"20to23"}.")\n";
		#print O "(".$countings{$chr}{"24to30"}.")\n";
		#print O "(".$countings{$chr}{"18to35"}.")\n";
		#print O "(".$contig{$chr}{"length"}.")\n";
    	
    #	if ($countings{$chr}{"21"} gt 0     ) { $countings{$chr}{"21"} =     ($countings{$chr}{"21"}     / $contig{$chr}{"length"} ) ;} else { $countings{$chr}{"21"}=0; }
	#	if ($countings{$chr}{"24to30"} gt 0 ) { $countings{$chr}{"24to30"} = ($countings{$chr}{"24to30"} / $contig{$chr}{"length"} ) ;} else { $countings{$chr}{"24to30"}=0; }
	#	if($countings{$chr}{"20to23"} gt 0  ) { $countings{$chr}{"20to23"} = ($countings{$chr}{"20to23"} / $contig{$chr}{"length"} ) ;} else { $countings{$chr}{"20to23"}=0; }
#		if($countings{$chr}{"18to35"} gt 0  ) { $countings{$chr}{"18to35"} = ($countings{$chr}{"18to35"} / $contig{$chr}{"length"} ) ;} else { $countings{$chr}{"18to35"}=0; }
	
	###calculating ratios
	
		if($countings{$chr}{"20to22"} eq 0 or $countings{$chr}{"25to29"} eq 0){
			$ratio_si_pi=0;
		}else{
	    	$ratio_si_pi = $countings{$chr}{"20to22"} / $countings{$chr}{"25to29"};
		}
		
		$tot_20_22_pos = $ch2{$chr}{"pos"}{"20"} + $ch2{$chr}{"pos"}{"21"} +$ch2{$chr}{"pos"}{"22"};
		$tot_20_22_neg = $ch2{$chr}{"neg"}{"20"} + $ch2{$chr}{"neg"}{"21"} +$ch2{$chr}{"neg"}{"22"};
		
		if($tot_20_22_pos eq 0 or $tot_20_22_neg eq 0){
			$ratio_20_22=0;
		}else{
			$ratio_20_22=$tot_20_22_pos  / $tot_20_22_neg;
			
			if($tot_20_22_pos <$tot_20_22_neg){
				$ratio_20_22= $tot_20_22_neg  / $tot_20_22_pos;
			}
			
		}
		 
		
		
		print $LOG_FH "printing counts ($chr)  \n";
#		print O $countings{$chr}{"21"}."\n";
		print O $countings{$chr}{"15to18"}."\n";
		print O $countings{$chr}{"20to22"}."\n";
#		print O $countings{$chr}{"20to23"}."\n";
		print O $countings{$chr}{"25to29"}."\n";
#		print O $countings{$chr}{"24to30"}."\n";
		print O $ratio_si_pi."\n";
		print O $ratio_20_22."\n";
		print O $countings{$chr}{"18to35"}."\n";
		print O $contig{$chr}{"length"}."\n";

		

		

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
