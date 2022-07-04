#!/usr/bin/perl
######################################################################
#######################################################################
# Copyright 2010 Fundação Oswaldo Cruz
# Author: Eric Aguiar
# template.pl is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 of the License.
#
# template.pl is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with template.pl (file: COPYING).
# If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################
#######################################################################



use Getopt::Long;
use Bio::SeqIO;


my $usage = "

$0 -i <filterblastfile> -f1 <fasta 1> -f2 <fasta 2> 
$0 -h

-i <input file>         : filterblast input file
-f1 <input file>        : fasta1 input file
-f2 <input file>        : fasta2 input file
-h                      : Help message

";

$| = 1;     # forces immediate prints into files rather than the buffer.

my $inputFile;
my $f1;
my $f2;

GetOptions ("i=s" => \$inputFile, "f1=s" => \$f1, "f2=s" => \$f2,"h!" => \$help);

if ($help) {
        die $usage;
}

if (not(defined($inputFile))) {
        die "\nGive an input file name",$usage;
}

if (not(defined($f1))) {
        die "\nGive an input file name",$usage;
}
if (not(defined($f2))) {
        die "\nGive an input file name",$usage;
}

open(F,"<$inputFile");
while(<F>){

	/Query=(.+)\tHit=(.+)\tSta.+/;
	$q = $1;
	$h=$2;
	$q =~ /.+length_(\d+).+/;
	$ql = $1;
	$h =~ /.+length_(\d+).+/;
	$hl=$1;
	if($ql>$hl){
		$remove{$h} =1;
	}else{
		 $remove{$q} =1;
	}
#	print "$q ;$ql  - $h;$hl \n";

}
close(F);

`cat $f1 $f2> aux_merge_contig`;
open(A,"<aux_merge_contig");

while(<A>){
	chomp;
	if(/^>(\S+)/){
		$print=0;
		if(not exists $remove{$1}){
			print "\n>$1\n";
			$print=1;
		}
	}else{
		#print "$print - $_";
		if($print eq 1){
			print 	$_;
		}
	}


}
close(A);
`rm -rf aux_merge_contig`;
