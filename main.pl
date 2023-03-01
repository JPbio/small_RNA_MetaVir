use strict;
use warnings;
# use autodie;

use Time::HiRes;
# use Time::HiRes qw(time usleep nanosleep);
use POSIX qw(strftime);
use Getopt::Long;
# use YAML qw(LoadFile);


use constant PATH_LOG_MAIN => "srna_metavir.main.log";

# $| = 1;     # forces immediate prints into files rather than the buffer.
my $timeStart = $^T;
my $timeStartStr = strftime("%Y-%m-%d %H:%M:%S", localtime($timeStart));

#######################################################################
### PARSE INPUTS ------------------------------------------------------
#######################################################################

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
# my $log; # TODO: 2023-02-27 - Restablish the custom log file(s) option
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
    # "log=s" => \$log,
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
### Print details -----------------------------------------------------
#######################################################################

my $runDetails = "

-------------------------------------------
> Details:

Start Time: $timeStartStr

> Reads: '$fasta';
> Reference: '$hostgenome';
> si: '$si';
> se: '$se';
> si: hash:'$hash';
> prefix: '$prefix';
";

if (defined($deg)) {
    $runDetails .= "> Searching for degradation: TRUE\n";
} else {
    $runDetails .= "> Searching for degradation: FALSE\n";
}

$runDetails .= "\n-------------------------------------------\n";

print $runDetails . "\n";

#######################################################################
### Create step folders -----------------------------------------------
#######################################################################

my $step0		="$prefix/$prefix"."_00_saet";
my $step1		="$prefix/$prefix"."_01_trimming";
my $step2		="$prefix/$prefix"."_02_filter_size_gaps_convertion";
my $step3		="$prefix/$prefix"."_03_mapping_vector";
my $step4		="$prefix/$prefix"."_04_getUnmapped";
my $step5_fix	="$prefix/$prefix"."_05_assembleUnmapped_fix";
my $step5_opt	="$prefix/$prefix"."_05_assembleUnmapped_opt";
my $step5_opt_fix="$prefix/$prefix"."_05_assembleUnmapped_opt_fix";
my $step5_opt_20to23="$prefix/$prefix"."_05_assembleUnmapped_opt_20to23";
my $step5_opt_24to30="$prefix/$prefix"."_05_assembleUnmapped_opt_24to30";
my $step5_contigs= "$prefix/$prefix"."_05_assembleUnmapped_final";
my $step5_cap3  = "$prefix/$prefix"."_05_cap3";
my $step6		="$prefix/$prefix"."_06_blast";
my $step7		="$prefix/$prefix"."_07_reportBlast";
my $step8		="$prefix/$prefix"."_08_completeReport";
my $step9		="$prefix/$prefix"."_09_contigs_no_hit";
my $step10		="$prefix/$prefix"."_10_pattern";
my $step_virus	="$prefix/$prefix"."_virus";


print "Creating folders...\n\n";
if (not -e $prefix) {
    `mkdir $prefix`;
}

if (not -e $step0) {
    `mkdir $step0`;
}
if (not -e $step1) {
    `mkdir $step1`;
}
if (not -e $step2) {
    `mkdir $step2`;
}
if (not -e $step3) {
    `mkdir $step3`;
}
if (not -e $step4) {
    `mkdir $step4`;
}
if (not -e $step5_fix) {
    `mkdir $step5_fix`
}
if (not -e $step5_opt) {
    `mkdir $step5_opt`;
}
if (not -e $step5_opt_fix) {
    `mkdir $step5_opt_fix`;
}
if (not -e $step5_contigs) {
    `mkdir $step5_contigs`;
}
if (not -e $step5_cap3) {
    `mkdir $step5_cap3`;
}
if (not -e $step5_opt_20to23) {
    `mkdir $step5_opt_20to23`;
}
if (not -e $step5_opt_24to30) {
    `mkdir $step5_opt_24to30`;
}
if (not -e $step6) {
    `mkdir $step6`;
}
if (not -e $step7) {
    `mkdir $step7`;
}
if (not -e $step8) {
    `mkdir $step8`;
}
if (not -e $step9) {
    `mkdir $step9`;
}
if (not -e $step10) {
    `mkdir $step10`;
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

print "\n\n";
print "#######################################################################\n";
print ">> New execution";
print $runDetails;


# # 
# # TODO: 2023-02-25 - Should we stop using this 'binary' thing?
# # 
# # our $binary = "/home/bioinfo/eric_bin";


# #######################################################################
# ### Handle FASTQ sequences --------------------------------------------
# #######################################################################

# # 
# # TODO: 2023-02-27 - Handle FastQ sequences
# # 

# #######################################################################
# ### Handle FASTA sequences --------------------------------------------
# #######################################################################

# my $nReadsUnmapHostBac;
# my $mappedbac;
# my $nReadsUnmapHost;

# if (defined($fasta)) {
#     print "#Loading FASTA file ... \n";

#     if (not defined($nohostfilter)) {

#         # 
#         # REVIEW: 2023-02-27 - Find a better way to manage file paths
#         # 

#         print "[COPING $fasta TO $step2/trimmed_filtered_gt15.fasta]\n";
#         `cp $fasta $step2/trimmed_filtered_gt15.fasta `;
#         print "\nSTEP3\n\t cp $fasta $step2/trimmed_filtered_gt15.fasta \n";

#     } else {

#         # Mapping Host - unfiltered reads against bacters reference
#         print "[MAPPING HOST-UNFILTERED READS AGAINST BACTERIAL GENOMES]... \n";
        
#         my $exec5_1 = "bowtie -f -S -v 1 --un $step4/unmappedVectorBacters.fasta -k 1 -p $process --large-index /media/data/reference/bacterial_genomes/all_bacters.fasta $fasta > /dev/null 2>> $step4/reads_mapped_to_bacteria.log ";
        
#         print LOG "\nSTEP5_1\n\t $exec5_1\n";
#         `$exec5_1`;

#         $nReadsUnmapHostBac = `grep -c '>' $step4/unmappedVectorBacters.fasta`;
#         chomp($nReadsUnmapHostBac);
#         $mappedbac = $nReadsUnmapHost - $nReadsUnmapHostBac;

#         print metrics "#reads mapped bacter\t".$mappedbac."\n";
        
#         # Mapped reads Bacterial genomes
#         print metrics "#preprocessed reads\t".$nReadsUnmapHostBac."\n";
        
#         # pre - processed reads
#         print "\n  PRE-PROCESSING FINISHED \n";
#     }
# }


# #######################################################################

# 
# TODO: 2023-03-01 - Encapsulate these date / time formatting stuff
# 

# Calculate elapsed time
my $timeElapsed = Time::HiRes::tv_interval([$timeStart]);
my $time_elapsed_ms = 1000 * ($timeElapsed - int($timeElapsed));
my $time_elapsed_str = strftime("%H:%M:%S", localtime($timeElapsed)) . sprintf ":%03d", ($time_elapsed_ms);

print "\n-- THE END --\n";
print "Time elapsed: $time_elapsed_str", "\n";

close(LOG_FH);


    
# my $filename = shift or die "Usage: $0 YAML-FILE\n";
# my $data = LoadFile($filename);
# print Dumper $data;
