use Bio::SeqIO;  
use Bio::Seq::Quality;
use Getopt::Long;

#
# TODO: 2023-05-15 - Use a decent logger
# 
use constant PATH_LOG_MAIN => "srna_metavir.main.log";

#######################################################################
### Parse inputs ------------------------------------------------------
#######################################################################

my $usage = "

$0 -c <counts file> -f <fasta> -p <prefix> 
$0 -h

-i <input file>         : Input file with columns delimited by a specific delimiter
-h                      : Help message

";

$| = 1;     # forces immediate prints into files rather than the buffer.

my $fasta;
my $counts ;
my $prefix;

GetOptions ("c=s" => \$counts,"f=s" => \$fasta,"p=s" => \$prefix, "h!" => \$help);

if ($help) {
        die $usage;
}

if (not(defined($counts))) {
        die "\nGive an counts file name",$usage;
}

if (not(defined($fasta))) {
        die "\nGive an input fasta file name",$usage;
}
if (not(defined($prefix))) {
        die "\nGive an input prefix name",$usage;
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

my @corr;
my %l;

####Loading fasta file
our %contig;
our %sense;
warn("Loading FASTA file...\n");
my $in = Bio::SeqIO->new(-format    => 'fasta',
                           -file      => $fasta);
while (my $data = $in->next_seq) {
#	print "--> ". $data->primary_id ."\n";
    $size = $data->length;
    $contig{$data->primary_id}{"sense"}=$data->seq;
    $contig{$data->primary_id}{"antisense"}=$data->revcom()->seq;
    
  # print $data->seq."\n";
  # print $data->revcom()->seq ."\n";
}



open(IN,"< $counts");

while (<IN>){
	
  if(/^\w+/){
  
	@campos = split(/\t/,$_);
	$c = $campos[0];	
    $p = $campos[1];
    $n = $campos[2];	
	if($p >$n ){
        $sense{$c} = 1;
	}else{
	        $sense{$c} = 0;
	}
  }
}


close(IN);

open(O2,">$prefix");


foreach $cand (keys %contig){
    print O2 ">$cand\n";
    if($sense{$cand} == 1){
        print O2 $contig{$cand}{"sense"}."\n";
    }else{
        print O2 $contig{$cand}{"antisense"}."\n";
    }
     
}
