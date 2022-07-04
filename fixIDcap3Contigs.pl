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

$0 -i <input file> -p <prefix> -size <size>
$0 -h

-i <input file>		: Fasta Input file
-p <prefix>		: Prefix name
-s <size>		: Size to select contigs in other file
-h			: Help message

";

$| = 1;     # forces immediate prints into files rather than the buffer.

my $inputFile;
my $prefix;
my $size;

GetOptions ("i=s" => \$inputFile,"p=s" => \$prefix,"s=s" => \$size,
			"h!" => \$help);

if ($help) {
	die $usage;
}

if (not(defined($inputFile))) {
	die "\nGive an input file name",$usage;
}

if (not(defined($prefix))) {
        die "\nGive an prefix file name",$usage;
}

if (not(defined($size))) {
        die "\nGive a valid size",$usage;
}

#print scalar keys %hash;

open(C,">$prefix.contigs");
open(C2,">$prefix.gt".$size.".fasta");

my $count=0;
#print "\nLoading Contigs... \n";
my $seqio_object = Bio::SeqIO->new(-file => $inputFile);
while($seq_object = $seqio_object->next_seq){
	$id=$seq_object->id."_$count" ;
  	$s=length($seq_object->seq);
	$seq=$seq_object->seq;

		
  	if (length($seq_object->seq)> $size){
   		 print C2 ">$id\t$s\n$seq\n";
	}else{
		print C ">$id\t$s\n$seq\n";
	}	
	$count++;
}

close(C);
close(C2);
