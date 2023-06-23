use Bio::SeqIO;
use Bio::Seq::Quality;
use Getopt::Long;

#
# TODO: 2023-05-15 - Use a decent logger
# 
use constant PATH_LOG_MAIN => "small_rna_metavir.main.log";

#######################################################################
### Parse inputs ------------------------------------------------------
#######################################################################

my $usage = "

$0 -sam <sam file> -fa <fasta file> -p <prefix name> -r <reads file> [--fastq]
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
my $prefix;
my $reads;
my $norm;
my $fastq;

GetOptions ("sam=s" => \$sam, "fa=s" => \$fasta, "p=s" => \$prefix, "r=s" => \$reads, "fastq!" => \$fastq, "profile!" => \$profile, "h!" => \$help);

if ($help) {
    die $usage;
}

if (not(defined($sam))) {
    die "\nGive an input file name", $usage;
}
if (not(defined($prefix))) {
    die "\nGive an prefix file name", $usage;
}
if (defined($profile) and not(defined($fasta))) {
    die "\nGive a fasta file ", $usage;
}

if (not(defined($reads))) {
    die "\nGive a reads input file ", $usage;
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
# select $LOG_FH;

#######################################################################
### Main --------------------------------------------------------------
#######################################################################

my $path_utils = "/small-rna-metavir/src/utils";
my $path_inner = "$path_utils/sam-sam-stranded-counts/inner";

my $path_sam_stats = "$path_inner/samStatistics_v2.pl";
my $path_define_strands = "$path_inner/defineStrandBasedOnCounts.pl";

if (defined($fastq)) {
    $type = "-q";
} else {
    $type = "-f";
}

`perl $path_sam_stats -sam $sam --counts > $prefix.stats`;

`perl $path_define_strands -c $sam.counts -f $fasta -p $prefix.stranded.fasta`;

`bowtie-build $prefix.stranded.fasta $prefix.stranded.fasta`;

`bowtie $type -S -p 40 -k 10 -v 1 $prefix.stranded.fasta $reads | awk -F '\t' '{if( \$2 != 4) print \$0}' > $prefix.v1.k10.sam`;

`rm -rf $prefix.stats $sam.counts`;