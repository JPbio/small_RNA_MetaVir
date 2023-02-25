#!/usr/bin/perl -w


#######################################################################
########################################################################
## Copyright 2010 Fundação Oswaldo Cruz
## Author: Adhemar Zerlotini Neto
## filterblast.pl is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 of the License.
##
## filterblast.pl is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY;
## without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
## See the GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License along with filterblast.pl (file: COPYING).
## If not, see <http://www.gnu.org/licenses/>.
##
########################################################################
########################################################################


use strict;
use Getopt::Long;
use Bio::SearchIO;


my $usage = "

$0 -b <blast output file> -pid <% identical> -plen <% length> -evalue <evalue> --best --desc --help
$0 -h

-b <blast output file>           : Output file resulting from blast
-pid <% identical>               : % identical
-plen <% length>                 : % query length aligned
-evalue <E-value>                : E-value
-best                            : Shows only the best hit
-desc                            : Shows only the subject description
-h                               : Help message

";


my $blastfile;
my $help = 0;
my $blastline;
my $pid = 0;
my $plen = 0;
my $evalue = 100;
my $best = 0;
my $desc = 0;
my $hit;

GetOptions ("b=s" => \$blastfile,
            "pid=s" => \$pid,
            "plen=s" => \$plen,
            "evalue=s" => \$evalue,
            "best!" => \$best,
            "desc!" => \$desc,
            "h!" => \$help);

if ($help) {
        die $usage;
}

if (not(defined($blastfile))) {
        die "\nGive a blast file name",$usage;
}


my $blast = new Bio::SearchIO(-format => "blast", -file   => "$blastfile");

while( my $result = $blast->next_result ) {

  while(my $hit = $result->next_hit) {
      while(my $hsp = $hit->next_hsp) {

         if( ($hsp->percent_identity >= $pid) && ($hsp->length('query')/$result->query_length*100 >= $plen) && ($hsp->evalue < $evalue) ) {

           print "Query=",   $result->query_name ,
                 "\tHit=",        $hit->name ; 
           print "\tDesc=",        $hit->description  if($desc);
           print "\tStart=",        $hsp->start('hit') ,
                 "\tEnd=",        $hsp->end('hit') ,
                 "\tE-value=",     $hsp->evalue ,
                 "\tQuery Length=",     sprintf("%.2f",$hsp->length('query') / $result->query_length * 100), '%' ,
                "\tIdentical=", sprintf("%.2f",$hsp->percent_identity) , "%\n";


         }
      }
      last if($best);
  }

}



