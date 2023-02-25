use strict;
use warnings;
# use YAML qw(LoadFile);
# use Data::Dumper qw(Dumper);
use Getopt::Long;

$| = 1;     # forces immediate prints into files rather than the buffer.

#######################################################################
### PARSE INPUTS ------------------------------------------------------

my $usage = "
------------------------------------------------------
USAGE:

$0 -fasta <reads in letter space without quality> -fastq <reads in fastq format with quality>  -qual <quality file> -fastqgz <reads in fastq format compressed gzip> [-nohostfilter]  -prefix <prefix>   -hostgenome <reference Host genome for reads filtering> -process <number of process to use> -size <expected size>  -log <log file name>
$0 -h

------------------------------------------------------
Args Description:

-fasta <reads are in fasta>    : Input reads are in fasta format without quality
-fastq <reads are in fasta>    : Input reads are in fastq format
-qual <reads .qual>            : Input quality file
-fastqgz <reads are in fastq.gz>	: Input reads are in fastq compressed as .gz format

-hostgenome <reference .fasta> : host reference genome
-size <expected size of genome>: The expected size of the genome to detect
-nohostfilter                      : Run directly the assembly without mapping to the host genome previously 

-si <size reads>	       		: Init range to use
-se <size reads>	       	    : End range to use
-hash <hash length>	       	    : hash to be used in velvet assembly

-process <number of process>   : number of process to use

-prefix <Output prefix>         : Output prefix folder name
-log <log file name>	       : log file name to storage the execution lines
-clean					: Clean large intermediate files

-h                             : Help message
";

# Required args
my $qual; # Quality file name
my $hostgenome;
my $process; # TODO: 2023-02-25 - Does 'proccess' stand for 'processors'?
my $size;
my $prefix; # TODO: 2023-02-25 - Should we call it simply 'output dir' (or something like it)?
my $fasta;
my $se;
my $si;

# 
# TODO: 2023-02-25 - Check optional args
# 

# Optional args
my $fastqgz;
my $log;
my $fastq;
my $hash;
my $clean;
my $nohostfilter;
my $large_index;
my $deg; # TODO: 2023-02-25 - Can we please call it 'degradation'?
my $help;

GetOptions("qual=s" => \$qual,
    "hostgenome=s" => \$hostgenome,
    "fasta=s" => \$fasta,
    "fastqgz=s" => \$fastqgz,
    "fastq=s" => \$fastq,
    "prefix=s" => \$prefix,
    "size=s" => \$size,
    "hash=s" => \$hash,
    "log=s" => \$log,
    "si=s" => \$si,
    "se=s" => \$se,
    "process=s" => \$process,
    "clean!" => \$clean,
    "nohostfilter!" => \$nohostfilter, #opstions must be lowcase and without "_"
    "degradation!" => \$deg,
    "largeindex!" => \$large_index,
    "h!" => \$help);

if ($help) {
    die $usage;
}

if (not defined($nohostfilter) and not defined($fasta) and not defined($fastqgz)) {
    if (not(defined $fastq)) {
        if ((not(defined($fastqgz))) and(not defined($fastq)) and(not defined($fasta))) {
            die "\nGive a valid input file!\n ", $usage;
        }

        if (not(defined($qual)) and(not defined($fastq)) and(not defined($fasta))) {
            die "\nGive a valid input quality file name! \n", $usage;
        }
    }
}

if (defined($large_index)) {
    $large_index = " --large-index ";
} else {
    $large_index = " ";
}

if (not(defined($hostgenome))) {
    die "\nGive a valid input reference file! \n", $usage;
}

if (not(defined($process))) {
    die "\nGive a number of process to use! \n", $usage;
}

if (not(defined($si)) or not defined($se)) {
    die "\nGive a reads range size! e.g.: -si 21 -se 23 \n", $usage;
}

if (not(defined($size))) {
    die "\nGive the genome's expected size! \n", $usage;
}

if (not(defined($prefix))) {
    die "\nGive an output folder prefix name! \n", $usage;
}

if (not(defined($hash))) {
    $hash = 15;
}


#######################################################################

print "\n\n\n#########################################################\n";
print "# Reads: $fasta\n";
print "# Reference: $hostgenome\n";
print "# si: $si\tse:$se\thash:$hash\n";
print "# prefix: $prefix\n";

if (defined($deg)) {
    print "# Searching for degradation: TRUE\n";
} else {
    print "# Searching for degradation: FALSE\n";
}



#######################################################################
print "\n\n-- THE END --\n\n"

    
# my $filename = shift or die "Usage: $0 YAML-FILE\n";
# my $data = LoadFile($filename);
# print Dumper $data;