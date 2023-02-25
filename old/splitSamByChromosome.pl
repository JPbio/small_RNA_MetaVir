#!/usr/bin/perl

######################################################################
########################################################################
## Copyright 2012 Universidade Federal de Minas Gerais
## RNA interference Laboratory
## Author: Eric Aguiar
## it's a free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 of the License.
##
## template.pl is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY;
## without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
## See the GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License along with this script (file: COPYING).
## If not, see <http://www.gnu.org/licenses/>.
##
########################################################################
########################################################################


use Getopt::Long;

my $usage = "

$0 -in <input file> -prefix <output prefix name> [-c chromosome name>]
$0 -h

-in <input file>           : Sam input file
-prefix <prefix output>    : Output prefix name
-c <chromosome name>    : get only specific chromosome
-h                         : Help message

";

$| = 1;     # forces immediate prints into files rather than the buffer.

my $in;
my $out;
my $chr;
my %hash;
GetOptions ("in=s"  => \$in,
	    "prefix=s" => \$out,
	    "c=s" => \$chr,
            "h!"    => \$help);

if ($help) {
        die $usage;
}

if (not(defined($in)) or not -s $in) {
        die "\nGive an input file valid",$usage;
}

if (not(defined($out))) {
        die "\nGive an output name valid",$usage;
}

open(IN,"<$in") or die ("\nImpossible to open the input file ($in) :( \n");
print "Reading sam file...\n";
$c1=0;
$c2=0;
my $header = "";
my %chrom;

while(<IN>){

#     print $_ ."\n";	
     
    if(/^[^@](\S+)\t(\d+)/){ 
	@c = split(/\t/,$_);
	if (not exists $chrom{$c[2]}{"lines"} ) {
		$chrom{$c[2]}{"lines"} ="";
	}
	$chrom{$c[2]}{"lines"} .= $_; 	
	

    }else{
    		$header .=$_;
	}
}

print "\nwriting files...\n";
foreach $c (keys %chrom){
	if(defined($chr) and $c eq $chr ){
		open(O,">$out.$c.sam");
                print O $header;
                print O $chrom{$c}{"lines"};
                close(O);
	}else{
		open(O,">$out.$c.sam");
		print O $header;
		print O $chrom{$c}{"lines"};
		close(O);
	}
}



