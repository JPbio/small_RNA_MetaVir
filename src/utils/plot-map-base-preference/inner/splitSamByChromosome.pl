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

use strict;

use warnings;
use Getopt::Long;

#
# TODO: 2023-05-15 - Use a decent logger
# 
use constant PATH_LOG_MAIN => "srna_metavir.main.log";

#######################################################################
### PARSE INPUTS ------------------------------------------------------
#######################################################################

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
my $help;

GetOptions(
    "in=s"     => \$in,
    "prefix=s" => \$out,
    "c=s"      => \$chr,
    "h!"       => \$help
);

if ($help) {
    die $usage;
}

if (not(defined($in)) or not -s $in) {
    die "\nGive a valid input file", $usage;
}

if (not(defined($out))) {
    die "\nGive a valid output name", $usage;
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
*STDERR = $LOG_FH;
select $LOG_FH;

#######################################################################
### Main... -----------------------------------------------------------
#######################################################################

open(IN, "<$in") or die ("\nImpossible to open the input file ($in) :( \n");

print $LOG_FH "Reading sam file...\n";

my $c1 = 0;
my $c2 = 0;
my $header = "";
my %chrom;

while (<IN>) {
    if (/^[^@](\S+)\t(\d+)/) {
        
		my @c = split(/\t/, $_);
        
		if (not exists $chrom{$c[2]}{"lines"}) {
            $chrom{$c[2]}{"lines"} = "";
        }
        
		$chrom{$c[2]}{"lines"} .= $_;

    } else {
        $header .= $_;
    }
}

print $LOG_FH "\nWriting files...\n";

foreach my $c (keys %chrom) {
    my $path_out_sam = "$out.$c.sam";
    open(O, ">$path_out_sam");
    print O $header;
    print O $chrom{$c}{"lines"};
    close(O);
}
