#!/usr/bin/perl
########################################################################
#########################################################################
### Copyright 2013 Universidade Federal de Minas Gerais
### Author:Francisco Lobo, Eric Aguiar
### E-mail: ericgdp@gmail.com
### This script is a free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 of the License.
###
### template.pl is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY;
### without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
### See the GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License along with template.pl (file: COPYING).
### If not, see <http://www.gnu.org/licenses/>.
###
#########################################################################
#########################################################################

use Getopt::Long;

my $usage = "

$0 -sam <sam input file> -pace <pace> -ref <genome size> -out <output name>
$0 -h

-sam  <input file>         	: Sam Input file
-pace <pace>         		: pace (window to analyze coverage)
-ref  <genome size>		: Reference genome size
-out <output name>      : Output name
-h                      	: Help message

";

$| = 1;     # forces immediate prints into files rather than the buffer.

my $infile;
my $pace;
my $length;
my %hash;
my $out;

GetOptions ("sam=s" => \$infile,"out=s" => \$out,
		"pace=i" => \$pace,
		"ref=i" => \$length,
                        "h!" => \$help);

if ($help) {
        die $usage;
}

if (not(defined($infile))) {
        die "\nGive an input file name",$usage;
}

if (not(defined($out))) {
        die "\nGive an output file name",$usage;
}
if (not(defined($length)) or $length <0 ) {
        die "\nGive a valid reference genome length",$usage;
}

if (not(defined($pace))  or ($pace <0 or $pace > $length)) {
        die "\nGive a valid pace",$usage;
}



open (IN, "$infile");

my @positions_foward;

my @positions_reverse;

$positions_foward[0] = "NULL";
$positions_reverse[0] = "NULL";

#give 0 density to all bases
for (my $i = 1; $i <= $length; $i++) {
  $positions_foward[$i] = 0;
  $positions_reverse[$i] = 0;
}

my $total = 0;

my $start;
my $end;
my $reference;


while (my $line = <IN>) {
 if($line =~/^\w.+/){ 
  chomp $line;
  my @aux = split (/\t/, $line);
  my $id = $aux[0];
  my $location = $aux[1];
  $a = $aux[5];	
  $a =~ /(\d+)M.*/;
  $seq_size = length($aux[9]);
  $reference = $aux[2];
  
  my $strand = $aux[1]; 
  my $strand_bn = $aux[1];  
  
  if($strand != 0 && $strand !=16){
  
  	if($strand & hex(0x10)){
  		$strand ="16";
  	}else{
  		$strand = "0";
  	}
  }
  
  if($strand eq "+" or $strand eq "0"){
       $start = $aux[3];
       $end = $start + $seq_size - 1;
	#print "$strand_bn----$id--$location--$seq_size--$strand-------------------------positivo\n";
  }else{
	
		$start = $aux[3] - $seq_size +1;;
        $end = $aux[3];
        #print "$strand_bn---$id--$location--$seq_size--$strand----------negativo\n";
  }
 # print "$strand_bn---$id--$location--$seq_size--$strand------st: $start --se: $end --\n";
  my $center = ($end - (($end - $start)/2));
  
my $round = int($center);
  if ($round > $length) {
    die "die  ($center) ($start) ($end) [$round > $length] \n";
  }
  my $sequence = $aux[9];
  my $length = $aux[4];
  if ($strand eq "+" or $strand eq "0" ) {
  	
    for (my $e = $start; $e <= $end; $e++) {
      $positions_foward[$e]++;
    }
  }
  elsif ($strand eq "-" or $strand eq "16") {
	
    for (my $e = $start; $e <= $end; $e++) {
    	#print "($e)position $start until $end \n"; 
	  $positions_reverse[$e]++;
    }
  }
 }
}

close IN;

open(O,">$out");

print O ("#REFERENCE\tSTART\tEND\tPOSITIVE\tNEGATIVE\tTOTAL\n");
my $sum_aux_p=0;
my $sum_aux_n=0;
my $sum_aux_t=0;

for (my $e = 1; $e <= $#positions_foward; $e = $e + $pace) {
  my $start = $e;
  my $end = 0;
  my $count = 0;
  my $total = 0;
  my $totalNeg = 0;
  my $sum=0;
  until ($count == $pace) {
    $total = $total + $positions_foward[($e + $count)];
    $totalNeg = $totalNeg + $positions_reverse[($e + $count)];
    $sum = $sum + $total + $totalNeg;
    $count++;
  }
  $end = ($e + $count - 1);
  #print "foward----- \n";
  print O ("$reference\t$start\t$end\t$total\t$totalNeg\t$sum\n");
  if($total gt 0){
    $sum_aux_p++;
  }
  if($totalNeg gt 0){
    $sum_aux_n++;
  }
  if($sum gt 0){
    $sum_aux_t++;
  }
}
close(O);
$p = ($sum_aux_p / $#positions_foward) *100;
$n = ($sum_aux_n / $#positions_foward) *100;
$t = ($sum_aux_t / $#positions_foward) *100;
print "$p#$n#$t\n";

#for (my $e = 1; $e <= $#positions_reverse; $e = $e + $pace) {
#  my $start = $e;
#  my $end = 0;
#  my $count = 0;
#  my $total = 0;
#  until ($count == $pace) {
#    $total = $total + $positions_reverse[($e + $count)];
#    $count++;
#  }
#  $end = ($e + $count - 1);
#  #print "reverse-> \n";
#  print ("NC_003038\tnt\tbase_coverage\t$start\t$end\t$total\t.\t.\tName=-\n");
#}


