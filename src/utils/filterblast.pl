
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
use warnings;

use Getopt::Long;
use Bio::SearchIO;

# 
# TODO: 2023-05-15 - Use a decent logger
# 
use constant PATH_LOG_MAIN => "srna_metavir.main.log";

#######################################################################
### PARSE INPUTS ------------------------------------------------------
#######################################################################

my $usage = "

$0 - b < blast output file > -pid < % identical > -plen < % length > -evalue < evalue > --best--desc--help
$0 - h

    -
    b < blast output file >: Output file resulting from blast -
    pid < % identical >: % identical -
    plen < % length >: % query length aligned -
    evalue < E - value >: E - value -
    best: Shows only the best hit -
    desc: Shows only the subject description -
    h: Help message

";

my $blastfile;
my $help = 0;
my $pid = 0;
my $plen = 0;
my $evalue = 100;
my $best = 0;
my $desc = 0;

GetOptions("b=s" => \$blastfile,
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
    die "\nGive a blast file name", $usage;
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
select $LOG_FH;

#######################################################################
### Main... -----------------------------------------------------------
#######################################################################

my $blast = new Bio::SearchIO(-format => "blast", -file => "$blastfile");

my $hit;
my $result;
my $hsp;

while ($result = $blast -> next_result) {
    
    while ($hit = $result -> next_hit) {
        
        while ($hsp = $hit -> next_hsp) {

            if (($hsp -> percent_identity >= $pid)
                && ($hsp -> length('query') / $result -> query_length * 100 >= $plen)
                && ($hsp -> evalue < $evalue)
            ) {

                print "Query=", $result -> query_name,
                    "\tHit=", $hit -> name;
                
                print "\tDesc=", $hit -> description if ($desc);
                
                print "\tStart=", $hsp -> start('hit'),
                    "\tEnd=", $hsp -> end('hit'),
                    "\tE-value=", $hsp -> evalue,
                    "\tQuery Length=", sprintf("%.2f", $hsp -> length('query') / $result -> query_length * 100), '%',
                    "\tIdentical=", sprintf("%.2f", $hsp -> percent_identity), "%\n";
            }
        }
        last if ($best);
    }
}