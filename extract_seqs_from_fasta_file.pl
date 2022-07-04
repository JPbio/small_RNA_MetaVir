#!/usr/bin/perl


if (@ARGV < 1) {
        print "usage: 

extract_seqs_from_fasta.pl <input list of IDs file> <input fasta file> <output fasta file>

<input list of IDs file>  :  list of IDs, one per line
<input fasta file>        :  fasta file name
<output fasta file>       :  fasta file name\n\n";
}


$listfile = shift or die "Must provide a list of list of IDs file name\n";
$infasta = shift or die "Must provide an input fasta file name\n";
$outfasta = shift or die "Must provide an output fasta file name \n";

$| = 1;                         # forces immediate prints into files rather than the buffer.

open (LIST, "<$listfile") or die "Could not open $listfile";
open (INFASTA, "<$infasta") or die "Could not open $infasta";
open (OUTFASTA, ">$outfasta") or die "Could not open $outfasta";

while (<LIST>){
  chomp;
  if (/(\S+)/){
	   print 	"ID:".$1."\n";
	  push @IDs,$1;
  } 
} 

while (<INFASTA>){
  if (/^\>(\S+).*$/){
 	print " (fasta) $1 \n"; 
   $key = "$1";
    $sequences{$key} .= $_;
  } else {
    $sequences{$key} .= $_;
  }
}


foreach $key (keys %sequences) {
  foreach $id (@IDs) {
    if($key eq $id) {
      print OUTFASTA "$sequences{$key}\n";
    }
  }
}
