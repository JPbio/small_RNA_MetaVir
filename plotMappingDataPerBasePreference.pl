#!/usr/bin/perl


use Bio::SeqIO;  
use Bio::Seq::Quality;
use Getopt::Long;

my $usage = "

$0 -sam <sam file> -s <start> -e <end> -fa <fasta> -pace <num> -p prefix --profile --keep --pattern
$0 -h

-i <input file>         : Input file with columns delimited by a specific delimiter
-h                      : Help message

";

$| = 1;     # forces immediate prints into files rather than the buffer.

our $PATH="/home/ubuntu/bin/";
my $sam;
my $delimiter;
my $index;
my $prefix;
my $profile;
my $fasta;
my $keep;
my $pattern;
our $minimo;
my $end;
my $antisense;
my $pace;


GetOptions ("sam=s" => \$sam,"pattern!" => \$pattern,"profile!" => \$profile,"p=s" => \$prefix,"s=i" => \$start,"pace=i" => \$pace,"norm=i" => \$norm,"e=i" => \$end,"fa=s" => \$fasta,"m=i" => \$minimo,
	"h!" => \$help,"keep!" => \$keep,"antisense!" => \$antisense);

if ($help) {
        die $usage;
}

if (not(defined($sam))) {
        die "\nGive a sam input file name",$usage;
}
if (not(defined($fasta))) {
        die "\nGive a fasta input file name",$usage;
}
if (not(defined($prefix))) {
        die "\nGive an prefix file name",$usage;
}

if (not(defined($pace))) {
       $pace=1;
}

if(not defined($start) or $start < 0){
    warn("\n\n\t\t[Bad value for start, start setted to 15]\n");
    $start=15;
}

if(not defined($end) or ($end < 0) or ($end < $start)){
   warn("\n\t\t[Bad value for end, end setted to 30]\n");
   $end=30;
}

if(not defined($minimo) or ($minimo < 0)){
   warn("\n\t\t[Bad value for min, min setted to 30]\n");
   $minimo=30;
}
open(IN,"<$sam");

our %tamanhos;
our %tamanhos_p;
our %tamanhos_n;
our $cp=0;
our $cn=0;
our %ch;
our %ch2;
our %base;
our %base_total;
our %base_per_size;
my $min=1000;
my $max=0;


for(my $i=15;$i<=35;$i++){   
            $base_per_size{$i}{"+"}{"A"}=0;
            $base_per_size{$i}{"+"}{"C"}=0;
            $base_per_size{$i}{"+"}{"G"}=0;
            $base_per_size{$i}{"+"}{"T"}=0;
            $base_per_size{$i}{"-"}{"A"}=0;
            $base_per_size{$i}{"-"}{"C"}=0;
            $base_per_size{$i}{"-"}{"G"}=0;
            $base_per_size{$i}{"-"}{"T"}=0;

}


my $spos ="0";
my $sneg="16";

if(defined($antisense)){
	 $spos ="16";
	 $sneg="0";
	

}


warn("Loading SAM file...\n");
while (<IN>){
	
  if(/^\w+/){
	@campos = split(/\t/,$_);
	$size = length($campos[9]);
	$c = $campos[2];	

	 if (not exists $tamanhos_p{$size} ){
            $tamanhos_p{$size} = 0;
            $base_total{$c}{"A"}=0;
	        $base_total{$c}{"C"}=0;
	        $base_total{$c}{"G"}=0;
	        $base_total{$c}{"T"}=0;
	        
	        $base_per_size{$size}{"+"}{"A"}=0;
            $base_per_size{$size}{"+"}{"C"}=0;
            $base_per_size{$size}{"+"}{"G"}=0;
            $base_per_size{$size}{"+"}{"T"}=0;
            $base_per_size{$size}{"-"}{"A"}=0;
            $base_per_size{$size}{"-"}{"C"}=0;
            $base_per_size{$size}{"-"}{"G"}=0;
            $base_per_size{$size}{"-"}{"T"}=0;
	        
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
	    $base{$c}{"pos"}{$size}{"A"}=0;
	    $base{$c}{"pos"}{$size}{"C"}=0;
	    $base{$c}{"pos"}{$size}{"G"}=0;
	    $base{$c}{"pos"}{$size}{"T"}=0;
	}
	
	if(not exists $ch2{$c}{"neg"}{$size}){
	    $ch2{$c}{"neg"}{$size}=0;
	    $base{$c}{"neg"}{$size}{"A"}=0;
	    $base{$c}{"neg"}{$size}{"C"}=0;
	    $base{$c}{"neg"}{$size}{"G"}=0;
	    $base{$c}{"neg"}{$size}{"T"}=0;
	}
	
	if($size < $min) { $min=$size;}
	if($size > $max) { $max=$size;}
	
	$first = uc(substr($campos[9],0,1));
	if($campos[1] eq 16){
     		$seq = uc($campos[9]);
    		$seq =~ tr /ATCG/TAGC/;
    		
    		$revseq=reverse $seq;
    		
    	    $first = uc(substr($revseq,0,1)); 
    }

 
    ###count miRNA preference
    if($size>=20 && $size <=24){
        $base_total{$c}{$first}= $base_total{$c}{$first} +1;
    }
    
	if ($campos[1]  eq $spos){   
	
        	
		$cp++;
		$ch{$c}{"pos"} = $ch{$c}{"pos"} + 1;
		$tamanhos_p{$size}= $tamanhos_p{$size} + 1;
		
		$ch2{$c}{"pos"}{"$size"}= $ch2{$c}{"pos"}{"$size"} +1;
		
		
		$base{$c}{"pos"}{$size}{$first}= $base{$c}{"pos"}{$size}{$first}+1;
		#print "POSITIVE in [$c] fita [pos]  size [$size] amount: ".$ch{$c}{"pos"}{$size} ."\n";
		$base_per_size{$size}{"+"}{$first}=		$base_per_size{$size}{"+"}{$first} +1;
	}elsif($campos[1]  eq $sneg){
	
		$cn++;
		$ch{$c}{"neg"} = $ch{$c}{"neg"} + 1; 
		$tamanhos_n{$size}= $tamanhos_n{$size} + 1;
		$ch2{$c}{"neg"}{$size}= $ch2{$c}{"neg"}{$size} + 1;
        
        $base{$c}{"neg"}{$size}{$first}= $base{$c}{"neg"}{$size}{$first}+1;
		#print "NEGATIVE in [$c] fita [pos]  size [$size] amount: ".$ch{$c}{"neg"}{$size} ."\n";
		$base_per_size{$size}{"-"}{$first}=	$base_per_size{$size}{"-"}{$first} +1;
	}

  }
}
$t = $cp+$cn;

#
#Load fasta file
#

our %contig;

warn("Loading FASTA file...\n");
my $in = Bio::SeqIO->new(-format    => 'fasta',
                           -file      => $fasta);
while (my $data = $in->next_seq) {
    $size = $data->length;
	$contig{$data->primary_id}=$size;
}


#
#writing file with total mappings
#

if(not defined($profile)){
    print "\n\n[Total read distribution]\n";
    for(my $i=$start; $i<=$end; $i++){
        $p=0; $n=0; 
        if(exists ($tamanhos_p{$i}) ){
            $p = $tamanhos_p{$i} ;
        }
        if(exists ($tamanhos_n{$i}) ){
            $n = $tamanhos_n{$i} ;
        }
        $t= $p+$n;
        print "$i\t$p\t$n\t$t\n";

    }
}
###################

#
#Print total of reads per chromosome
#


if(defined($profile)){

    print "\n\nReads per chromosoes\n";
    
    open(O,">$prefix.readsPerChromossome");
    
    foreach $c (keys %ch){
        $pa = $ch{$c}{"pos"};
        $na=  $ch{$c}{"neg"};
        $ta=$pa+$na;
        
        if($pa <=0){ $pa=0;} 
        if($na <=0){ $na=0;} 
            
        print O $c."\t". $ta . "\t" . $na .  "\t". $ta ."\n";
    }
    
    close(O);
    
    #`R --no-save --args $prefix.$c.distribution < file.R`;
    
	
##############

    print "\n\nReads per chromosoes\n";

    foreach $c (keys %ch2){
        open(O,">$prefix.$c.distribution");
        open(O2,">$prefix.$c.base_distribution");
    
        
        print O2 "\tA\tC\tG\tT\n"; 
    
       #print "\n\nChromossome [$c]\n";
        for (my $s=$start;$s<=$end;$s++){
            $pa = $ch2{$c}{"pos"}{$s};
            $na=  $ch2{$c}{"neg"}{$s};
            
            $bpa = $base{$c}{"pos"}{$s}{"A"};
            $bpc = $base{$c}{"pos"}{$s}{"C"};
            $bpg = $base{$c}{"pos"}{$s}{"G"};
            $bpt = $base{$c}{"pos"}{$s}{"T"};
            
            $bna = $base{$c}{"neg"}{$s}{"A"} *-1;
            $bnc = $base{$c}{"neg"}{$s}{"C"}*-1;
            $bng = $base{$c}{"neg"}{$s}{"G"}*-1;
            $bnt = $base{$c}{"neg"}{$s}{"T"}*-1;
            
        
            if($pa <=0){ $pa=0;} 
            if($na <=0){ $na=0;} 
       
            $ta=$pa+$na;
            print O $s."\t$pa\t$na\t$ta\n";
            print O2 $s."\t$bpa\t$bpc\t$bpg\t$bpt\n";
            print O2 $s."\t$bna\t$bnc\t$bng\t$bnt\n";

        }
        close(O);
        #drawing distribution plot
        #`R --no-save --args $prefix.$c.distribution < file.R`;
    }
    
    ###creating sam files
    $split="/home/ubuntu/bin/splitSamByChromosome.pl -in $sam -prefix $prefix";
    `$split`;
   
    #calcDensityPerBase
    foreach $c (keys %contig){
        $refsize= $contig{$c};
        $cp="$prefix.$c";
        
        $counting=`grep -v "@" $cp.sam | wc -l  | cut -f 1 -d ' '`;
    	$counting=~ s/\n//g;
    	
        if($counting >= $minimo) {
        		print "####### $cp was USED, $counting > $minimo \n";
        
				$f0="/home/ubuntu/bin/filterSamBowtieBySize.pl -i $cp.sam -si 20 -se 23 -o $cp.20to23.sam";
				$f1="/home/ubuntu/bin/filterSamBowtieBySize.pl -i $cp.sam -si 21 -se 21 -o $cp.21.sam";
				$f2="/home/ubuntu/bin/filterSamBowtieBySize.pl -i $cp.sam -si 24 -se 30 -o $cp.24to29.sam";
				$f3="/home/ubuntu/bin/filterSamBowtieBySize.pl -i $cp.sam -si 22 -se 22 -o $cp.22to22.sam";
				$f4="/home/ubuntu/bin/filterSamBowtieBySize.pl -i $cp.sam -si 15 -se 19 -o $cp.15to19.sam";
	  
	  
				`$f0`;
				`$f1`;
				`$f2`;
				`$f3`;
				`$f4`;
	   
				#print "$f1\n$f2\n";
 
				$comm0="/home/ubuntu/bin/calcDensityPerBase.pl -sam $cp.20to23.sam -pace $pace -ref ".$contig{$c}."  -out  $cp.20to23.density";
				$comm="/home/ubuntu/bin/calcDensityPerBase.pl -sam $cp.sam -pace $pace -ref ".$contig{$c}."  -out  $cp.density";
				$comm2="/home/ubuntu/bin/calcDensityPerBase.pl -sam $cp.21.sam -pace $pace -ref ".$contig{$c}."  -out $cp.21.density";
				$comm3="/home/ubuntu/bin/calcDensityPerBase.pl -sam $cp.24to29.sam -pace $pace -ref ".$contig{$c}."  -out $cp.24to29.density";
				$comm4="/home/ubuntu/bin/calcDensityPerBase.pl -sam $cp.22to22.sam -pace $pace -ref ".$contig{$c}."  -out $cp.22to22.density";
				$comm5="/home/ubuntu/bin/calcDensityPerBase.pl -sam $cp.15to19.sam -pace $pace -ref ".$contig{$c}."  -out $cp.15to19.density";
		
				# print "$comm \n $comm2 \n $comm3 \t $comm0 \n";
				`$comm0`;
				`$comm`;
				`$comm2`;
				`$comm3`;
				`$comm4`;
				`$comm5`;
		
			#	print "Plotting graphs [$c]...\n";
			   # `/opt/R/4.1.2/bin/R --no-save $prefix.$c.distribution $cp.distribution < $PATH/plorDistribution_ggplot2.R 2>/dev/null `;
			#`/opt/R/4.1.2/bin/R --no-save $cp.density $cp.all.density < $PATH/plotDensity_ggplot2.R 2>/dev/null`;
				#`/opt/R/4.1.2/bin/R --no-save $cp.21.density $cp.21.density < $PATH/plotDensity_ggplot2.R 2>/dev/null`;       
				#`/opt/R/4.1.2/bin/R --no-save $cp.24to29.density $cp.24to29.density < $PATH/plotDensity_ggplot2.R 2>/dev/null`;


				`R --no-save $prefix.$c.distribution $cp.distribution < $PATH/plorDistribution_ggplot2.R  2>/dev/null`;
				`R --no-save $cp.density $cp.all.density < $PATH/plotDensity_ggplot2.R 2>/dev/null`;
				`R --no-save $cp.21.density $cp.21.density < $PATH/plotDensity_ggplot2.R 2>/dev/null`;       
				`R --no-save $cp.24to29.density $cp.24to29.density < $PATH/plotDensity_ggplot2.R 2>/dev/null`;
				`R --no-save $cp.20to23.density $cp.20to23.density < $PATH/plotDensity_ggplot2.R 2>/dev/null`;
				`R --no-save $cp.22to22.density $cp.22to22.density < $PATH/plotDensity_ggplot2.R 2>/dev/null`;
				`R --no-save $cp.15to19.density $cp.15to19.density < $PATH/plotDensity_ggplot2.R 2>/dev/null`;
				`R --no-save $prefix.$c.base_distribution $prefix.$c.base_distribution < $PATH/calcDistributionPerBase.R 2>/dev/null`;
				`R --no-save $prefix.$c.base_distribution $prefix.$c.publication < $PATH/calcDistributionPerBase_publication2.R 2>/dev/null`;
				
				
        }else{
        	print "--> $cp was [DISCARDED], $counting < $minimo <--\n";
        
        }
         #`rm -rf  $cp.21.sam $cp.24to29.sam $cp.sam`;
        if(not defined($keep)){
        	`rm -rf $cp.20to23.density $cp.15to19.density $prefix.$c.distribution $cp.density $cp.21.density $cp.24to29.density $cp.22to22.density`;
        }
       
        
        #`rm -rf $cp.density $cp.21.density $cp.24to29.density $cp.21.sam $cp.24to29.sam $cp.sam`;
        
      #  $exec="R --no-save --args $prefix.$c.distribution $cp.density $cp.21.density $cp.24to29.density $refsize $c < /Users/eric/Desktop/projetos/scripts_servers/R/plotMappingDensity.R ";
      #  `$exec`;
        

        #print "$comm \n";
        
    }
    #calcDensityPerBase.pl -sam -pace -ref ";

}

###miRNA preference print
open(O3,">$prefix.miRNA_total_base_distribution");  
print O3 "\tA\tC\tG\tT\n"; 
foreach $c (keys %ch2){    

    $ta=$base_total{$c}{"A"};
	$tc=$base_total{$c}{"C"};
	$tg=$base_total{$c}{"G"};
	$tt=$base_total{$c}{"T"};
    print O3 $c."\t$ta\t$tc\t$tg\t$tt\n";
}
close(O3);	


############## GERAL preference print
open(O4,">$prefix.Geral_base_distribution");  
print O4 "size#A#C#G#T\n"; 
for(my $i=15;$i<=35;$i++){   

    $ta=$base_per_size{$i}{"+"}{"A"};
	$tc=$base_per_size{$i}{"+"}{"C"};
	$tg=$base_per_size{$i}{"+"}{"G"};
	$tt=$base_per_size{$i}{"+"}{"T"};
	
	$nta=$base_per_size{$i}{"-"}{"A"} * -1;
	$ntc=$base_per_size{$i}{"-"}{"C"} * -1;
	$ntg=$base_per_size{$i}{"-"}{"G"} * -1;
	$ntt=$base_per_size{$i}{"-"}{"T"} * -1;
	
    print O4 $i."#$ta#$tc#$tg#$tt\n";
    print O4 $i."#$nta#$ntc#$ntg#$ntt\n";
}
close(O4);
print "Plotting graphs...\n";
#`montage -mode concatenate -tile 3x $prefix*density.png $prefix.merge.pdf`;
#`convert $prefix*density.png +append $prefix.merge.2.pdf`;
#`convert $prefix*density.png -append $prefix.merge.3.pdf`;
 
       
#`rm -rf *.png`;
#`rm -rf $prefix*.distribution  $prefix*.sam $prefix*.density`;
 
 #$exec="R --no-save --args $prefix.readsPerChromossome < /Users/eric/Desktop/projetos/scripts_servers/R/plotMappingDistribution.R ";
#`$exec`;


















