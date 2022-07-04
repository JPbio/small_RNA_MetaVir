#!/usr/bin/perl



use Bio::SeqIO;  
use Bio::Seq::Quality;
use Getopt::Long;

my $usage = "

$0 -sam <sam file> -fa <fasta file> [--profile] [--counts]
$0 -h

-i <input file>         : Input file with columns delimited by a specific delimiter
-h                      : Help message

#
# calc stattistics by reads
#
";

$| = 1;     # forces immediate prints into files rather than the buffer.


my $sam;
my $delimiter;
my $index;
my $profile;
my $fasta;
my $counts;

GetOptions ("sam=s" => \$sam,"fa=s" => \$fasta,"profile!" => \$profile,"counts!" => \$counts,
	"h!" => \$help);

if ($help) {
        die $usage;
}

if (not(defined($sam))) {
        die "\nGive an input file name",$usage;
}


open(IN,"<$sam");

our %tamanhos;
our %tamanhos_p;
our %tamanhos_n;
our $cp=0;
our $cn=0;
our %ch;
our %ch2;

my $min=1000;
my $max=0;

if(defined($fasta)){
	our %contig;

	warn("Loading FASTA file...\n");
	my $in = Bio::SeqIO->new(-format    => 'fasta',
							   -file      => $fasta);
	while (my $data = $in->next_seq) {
		$size = $data->length;
		$contig{$data->primary_id}=$data->seq;
	}



}

my %count_piRNAs;

my %reads;
while (<IN>){
	
  if(/^\w+/){
	@campos = split(/\t/,$_);
	$size = length($campos[9]);
	$c = $campos[2];	
	if(not exists($reads{$campos[2]}{$campos[0]})){
	
	$reads{$campos[2]}{$campos[0]}=1;
	 
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
		#print "POSITIVE in [$c] fita [pos]  size [$size] amount: ".$ch{$c}{"pos"}{$size} ."\n";
		
	}elsif($campos[1] eq "16"){
	
		$cn++;
		$ch{$c}{"neg"} = $ch{$c}{"neg"} + 1; 
		$tamanhos_n{$size}= $tamanhos_n{$size} + 1;
		$ch2{$c}{"neg"}{$size}= $ch2{$c}{"neg"}{$size} + 1;
		#print "NEGATIVE in [$c] fita [pos]  size [$size] amount: ".$ch{$c}{"neg"}{$size} ."\n";
		
	}

  	
  }else{
  	#read exists
  }
  
  }
}
$t = $cp+$cn;


print "Total in both strand: $t\nTotal in positve strand: $cp\nTotal in negative strand: $cn\n";
print "\nSize\tSENSE\tANTISENSE\tTOTAL\n";
foreach (sort { $a <=> $b } keys(%tamanhos_p) )
{
	$s = $tamanhos_p{$_}  + $tamanhos_n{$_};
        if ($s ne "0"){
                print "$_\t".$tamanhos_p{$_}."\t".$tamanhos_n{$_}."\t".$s. "\n";
        }

}

foreach $c (keys %ch){
        $pa = $ch{$c}{"pos"};
        $na=  $ch{$c}{"neg"};
        $ta=$pa+$na;
        
        if($pa <=0){ $pa=0;} 
        if($na <=0){ $na=0;} 
            
        print O $c."\t". $ta . "\t" . $na .  "\t". $ta ."\n";
    }
if (defined($profile)){
	open(S,">$sam.total.distribution");
	print S "\nSize\tSENSE\tANTISENSE\tTOTAL\n";
	for (my $i=15;$i<30;$i++)
	{
		
		$p=0;
		$n=0;
		if(exists $tamanhos_p{$i}){ $p = $tamanhos_p{$i}; }
		if(exists $tamanhos_n{$i}){ $n = $tamanhos_n{$i}; }
		
		$t=$p+$n;
		
        print S "$i\t$p\t$n\t$t\n";

	}

	close(S);
}



print "\n\nReads per chromosoes\n";

if(defined($counts)){
    open(COUNT,">$sam.counts");
}
foreach $c (keys %ch){
	$s = $ch{$c}{"pos"} + $ch{$c}{"neg"} ;
	print $c."\t". $ch{$c}{"pos"} . "\t" . $ch{$c}{"neg"} .  "\t". $s ."\n";
	
	if(defined($counts)){
	    print COUNT $c."\t". $ch{$c}{"pos"} . "\t" . $ch{$c}{"neg"} .  "\t". $s ."\n";
	}
}

if(defined($counts)){
    close(COUNT);
}


##
# Analyzing profile of small RNAs
##

my %signal;

foreach $c (keys %ch2){
     $signal{$c}{"si" } = "no" ;
     $signal{$c}{"pi" } = "no" ;
     
	$sumT1P=0;
	$sumT1N=0;
	$allT1=0;
	$sumT2P=0;
	$sumT2N=0;
	$allT2=0;
	$scorePI=0;
	for (my $s=$min;$s<=$max;$s++){
        $pa = $ch2{$c}{"pos"}{$s};
        $na=  $ch2{$c}{"neg"}{$s};
        
        if($pa <=0){ $pa=0;} 
        if($na <=0){ $na=0;} 
       
        $ta=$pa+$na;
        
        if($s > 15 and $s<21){
        	$sumT1P+=$pa;
        	$sumT1N+=$pn;
        	$allT1= $sumT1P + $sumT1N;
        }elsif($s==21){
        	
        	if($pa > 0 and $na>0){
        		if(($pa > 50 or $na > 50) and ($ta > $allT1) ){
        			$r=$pa / $na;
        			$r=$na/$pa if $na>$pa;
        			
        			if($r > 0.3 and $r < 10){
        				$signal{$c}{"si" } = "yes" ;
        			} 
        		}
        	}
        }elsif($s==22 or $s==23){
        	$sumT2P+=$pa;
        	$sumT2N+=$na;
        	$allT2= $sumT2P +$sumT2N;
        
        }elsif($s>24 and $s<=29){
        	
        	if($pa > 0 and $na>0){
        			$r2=$pa / $na;
        			$r2=$na/$pa if $na>$pa;
				if(($pa > 50 and $ta > 50) ){
					$scorePI +=1;
				}
		
			    $allT3+=$pa+$ta;  
		    }  
        }elsif($s==30){
        	
        	if( ($allT3 > $allT2 +$allT1 ) and ($scorePI == 5) ) {
        		 $signal{$c}{"pi" } ="yes" ;
        	}
        
        }
        
        
        #print $s."\t$pa\t$na\t$ta\n";
    }
}

print " \n###############SIGNAL of siRNA or piRNA #############\n\n";
print "Contig\tsiRNA\tpiRNA\n" ;
foreach $c (keys %signal){

	print "$c\t".$signal{$c}{"si" } ." \t" . $signal{$c}{"pi" }." \n" ;
}

if (defined($profile)){
	open(S,">$sam.profile");
	open(SI,">$sam.withSiRNA.fasta");
	open(PI,">$sam.withPiRNA.fasta");
	open(SIPI,">$sam.withSiRNA_and_PiRNA.fasta");
	
	foreach $c (keys %signal){
		print S " $c\t".$signal{$c}{"si" } ." \t" . $signal{$c}{"pi" }." \n" ;
		if($signal{$c}{"si"} eq "yes" and $signal{$c}{"pi"} eq "no" ){
			print SI ">$c\n".$contig{$c}."\n";
		}
		if($signal{$c}{"pi"} eq "yes" and $signal{$c}{"si"} eq "no"){
			print PI ">$c\n".$contig{$c}."\n";
		}
		if($signal{$c}{"pi"} eq "yes" and $signal{$c}{"si"} eq "yes"){
			print SIPI ">$c\n".$contig{$c}."\n";
		}
	}
	
	close(S);
	close(SI);
	close(PI);
	close(SIPI);
}

print "\n\nReads per chromosoes\n";

if($min > 5) { $min=5;} 
if($max < 30) { $max=30;} 


foreach $c (keys %ch2){
	print "\n\nChromossome [$c]\n";
	for (my $s=$min;$s<=$max;$s++){
        $pa = $ch2{$c}{"pos"}{$s};
        $na=  $ch2{$c}{"neg"}{$s};
        
        if($pa <=0){ $pa=0;} 
        if($na <=0){ $na=0;} 
       
        $ta=$pa+$na;
        print $s."\t$pa\t$na\t$ta\n";
    }
}


