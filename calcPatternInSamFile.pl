#!/usr/bin/perl



use Bio::SeqIO;  
use Bio::Seq::Quality;
use Getopt::Long;
use Statistics::Basic qw(:all);
use Statistics::RankCorrelation;


my $usage = "

$0  -s <sam file> -o <output> 
$0 -h


-s  <sam file>			: Sam file from mapping on putative sequences
-o 	<output file>        : output file name
-h                      : Help message

";

$| = 1;     # forces immediate prints into files rather than the buffer.


my $ref;
my $search;
my $delimiter;
my $index;
my $profile;
my $output;
my $siRNA;
my $piRNA;
my $all;

our %zscore;
our $z_stddev_p ;
our $z_stddev_n ;
our $z_stddev_t ;
our $z_sumt=0;
our @zvp, @zvn,@zvt;

GetOptions ("s=s" => \$search,"o=s" => \$output,"siRNA!" => \$siRNA,"piRNA!" => \$piRNA,"all!" => \$all,
	"h!" => \$help);

if ($help) {
        die $usage;
}

if (not(defined($search))) {
        die "\nGive an search file name",$usage;
}


analyzeProfile($search,0);

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

	print "Analyzing $file...\n";
	open(IN,"<$file");
	while (<IN>){
	
		if(/^\w+/){
			@campos = split(/\t/,$_);
			$size = length($campos[9]);
			$c = $campos[2];	

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
	$t = $cp+$cn;
	
	my @pos;
	my @neg;
	my @tot;

	
	
	print "Calculating Z-score...\n";


    open(ALL,">$output");
    print ALL "Chromosome#";
    for (my $i=0;$i<=41;$i++){
        print ALL "$i#";
    }
    print ALL "\n";
    foreach $chr (keys %ch2){   #for each sequence
            my @pos;
            my @neg;
            my $tot;
            my $sumt=0;
            
            #initializing array
            for (my $i=0;$i<=41;$i++){
	        	$tot[$i]=0;	
        	}
            
          #  open(O,">score.$chr.tab");
          #  print O "$chr\n";
            for (my $i=0;$i<=41;$i++){
                    $p=0;
                    $n=0;
                    if(exists $ch2{$chr}{"pos"}{$i}){ $p = $ch2{$chr}{"pos"}{$i}; }
                    if(exists $ch2{$chr}{"neg"}{$i}){ $n = $ch2{$chr}{"neg"}{$i}; }
                   
    
                 $tot[$i-15]=$p;
                    $tot[$i+6]=$n;
           #         print  "$i\t$p\t$n\t$t\n";
            }

        ##z = (x - mean)/sd

        my $stddev_t   = stddev(@tot);
        my $mean_t   = mean(@tot);

        my @zp;
        my @zn;
        my @zt;
        
        if($stddev_t > 0 and $mean_t){
            for (my $i=0;$i<=41;$i++){
                $zscore{$i}{"both"}     = ($tot[$i] - $mean_t) / $stddev_t;
                push(@zvt,($tot[$i] - $mean_t) / $stddev_t);
            }
            print ALL "$chr#";
            for (my $i=0;$i<=41;$i++){
            #        print O $zscore{$i}{"both"}."\n";
                    print ALL $zscore{$i}{"both"}."#";
            }
            print ALL "\n";
        }
    undef(@zscore);
    undef @pos;
    undef @neg;
    undef @tot;
    undef @vp;
    undef @vn;
    undef @vt;
   # close(O);
    }
    close(ALL);
}



`R --no-save $output $output < /home/ubuntu/bin/plotPattern_publication.R`;
`R --no-save $output $output < /home/ubuntu/bin/plotPattern_publication_both_strands.R`;








