use strict;
use warnings;
# use autodie;

use Time::HiRes;
# use Time::HiRes qw(time usleep nanosleep);
use POSIX qw(strftime);
use Getopt::Long;
# use YAML qw(LoadFile);


#
# TODO: 2023-03-07 - Use a decent logger
# 
use constant PATH_LOG_MAIN => "srna_metavir.main.log";
use constant PATH_REF_BACTERIA => "/srna_metavir/asset/bacterial_genomes/all_bacters.fasta";

# 
# TODO: 2023-03-01 - Standardize naming as snake_case
# 

# $| = 1;     # forces immediate prints into files rather than the buffer.
my $time_start = $^T;
my $time_start_str = strftime("%Y-%m-%d %H:%M:%S", localtime($time_start));

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

# 
# TODO: 2023-02-25 - Should we call it simply 'output dir' (or something like it)?
# TODO: 2023-03-06 - Reenable custom naming for this folder
# 
# my $prefix = strftime("exec_%Y%m%d_%H%M%S", localtime($time_start));
my $prefix = "exec_20230307_120529";

my $se;
my $si;

# Conditionally required args

#
# REVIEW: 2023-03-01 - Shall we standardize 'path' variables (prefix 'path'?)
#

my $fasta;
my $fastq;

# 
# TODO: 2023-02-25 - Check optional args
# 

# Optional args
my $fastqgz;
# my $log; # TODO: 2023-02-27 - Restablish the custom log file(s) option
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
    # "prefix=s" => \$prefix, # TODO: 2023-03-06 - Reenable custom naming for this folder
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

Start Time: $time_start_str

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

#######################################################################
### Create step folders -----------------------------------------------
#######################################################################

# 
# REVIEW: 2023-02-25 - Should we stop using this 'binary' thing?
# 
my $path_utils = "/srna_metavir/src/utils";
# our $binary = "/home/bioinfo/eric_bin";

# 
# REVIEW: 2023-03-01 - Shall we standardize these directories as all other 'path' variables too? ('path' prefix?)
# 

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


# #######################################################################
# ### Handle FASTQ sequences --------------------------------------------
# #######################################################################

# 
# TODO: 2023-02-27 - Handle FastQ sequences
# 

# #######################################################################
# ### Handle FASTA sequences --------------------------------------------
# #######################################################################

my $path_bacterial_genome_all = "/srna_metavir/asset/bacterial_genomes/all_bacters.fasta";

# 
# REVIEW: 2023-03-01 - Shall we think of better names?
# 

my $mappedbac;
my $nReadsUnmapHostBac;
# my $nReadsUnmapHost;

if (defined($fasta)) {
    print "#Loading FASTA file ... \n";

    if (not defined($nohostfilter)) {

        # 
        # REVIEW: 2023-02-27 - Find a better way to manage file paths
        #

        my $path_trimmed_filtered_gt15 = "$step2/trimmed_filtered_gt15.fasta";
        my $cmd = "cp $fasta $path_trimmed_filtered_gt15";

        print "[COPING $fasta TO $path_trimmed_filtered_gt15]\n";
        `$cmd`;
        
        print "\nSTEP3\n\t $cmd\n";

    } else {

        # 
        # REVIEW: 2023-03-01 - Check this step numbering
        # TODO: 2023-03-01 - Test it
        # 

        # Mapping Host - unfiltered reads against bacters reference
        print "[MAPPING HOST-UNFILTERED READS AGAINST BACTERIAL GENOMES]... \n";
        

        my $path_unmapped_vector_bacters = "$step4/unmappedVectorBacters.fasta";

        my $exec5_1 = "bowtie -f -S -v 1 --un $path_unmapped_vector_bacters -k 1 -p $process --large-index $path_bacterial_genome_all $fasta > /dev/null 2>> $step4/reads_mapped_to_bacteria.log ";
        
        print "\nSTEP5_1\n\t $exec5_1\n";
        `$exec5_1`;

        $nReadsUnmapHostBac = `grep -c '>' $path_unmapped_vector_bacters`;
        chomp($nReadsUnmapHostBac);

        # Mapped reads Bacterial genomes
        # print metrics "#preprocessed reads\t".$nReadsUnmapHostBac."\n"; # TODO: 2023-03-01 - Check this 'metrics' log
        print "#preprocessed reads\t $nReadsUnmapHostBac \n"; # REVIEW: 2023-03-01 - Does this message really make sense?
        
        # $mappedbac = $nReadsUnmapHost - $nReadsUnmapHostBac; # REVIEW: 2023-03-01 - It doesn't look like it works
        # print metrics "#reads mapped bacter\t".$mappedbac."\n"; # TODO: 2023-03-01 - Check this 'metrics' log
        # print "#reads mapped bacter\t".$mappedbac."\n";
        
        print "\n  PRE-PROCESSING FINISHED \n"; # REVIEW: 2023-03-01 - Does this message really make sense?
    }
}

my $path_plot_dist_per_base_by_reads = "$path_utils/plot-geral-dist-base-reads/plotGeralDistributionPerBaseByReads.pl";

if (not defined($nohostfilter)) {
    if (defined($fastqgz)) {

        # 
        # TODO: 2023-03-06 - Test it
        # 

        # # # dealing with FASTQ TRIMMED files compacted with GZIP

        print "\n\nRunning step 2 [ converting fastq to fasta - fastq_to_fasta ]\n";

        #converting fastq to fasta
        my $exec_fq_2 = "gunzip -dc $fastqgz | fastq_to_fasta -Q 33 -o $step2/trimmed_filtered_gt15.fasta";
        print "\nSTEP2\n\t $exec_fq_2\n";
        `$exec_fq_2`;
    }

    print "[MAPPING SEQUENCE AGAINST VECTOR]\n";
    my $exec3 = "bowtie $large_index -f -S -k 1 -p $process -v 1 --un $step4/unmappedVectorReads.fasta $hostgenome $step2/trimmed_filtered_gt15.fasta | awk -F'\\t' '{if( \$2 != 4) print \$0}' > $step3/mapped_host.v1.sam  2>mapping_host.stats  ";
    print "\nSTEP3\n\t $exec3\n";
    `$exec3`;

    # #count total reads
    my $nReads = `grep -c '>' $step2/trimmed_filtered_gt15.fasta`;
    chomp($nReads);

    # print metrics "#total reads\t".$nReads. # REVIEW: 2023-03-06 - What should we do with these 'special' logs?
    print "#total reads\t".$nReads.
    "\n";
    # print interest "#total reads\t".$nReads. # REVIEW: 2023-03-06 - What should we do with these 'special' logs?
    print "#total reads\t".$nReads.
    "\n";


    print "[PLOT GERAL DISTRIBUTION READS MAPPED HOST ]\n";
    my $exec3_1 = "$path_plot_dist_per_base_by_reads -sam $step3/mapped_host.v1.sam  -s $si -e $se -p $step3/$prefix.mapped_host.v1.$si.$se -norm $nReads --plot";
    print "\nSTEP3\n\t $exec3_1\n";
    `$exec3_1`;

    print "[PLOT GERAL DISTRIBUTION READS MAPPED HOST - 15-35nt ]\n";
    my $exec3_11 = "$path_plot_dist_per_base_by_reads -sam $step3/mapped_host.v1.sam  -s 15 -e 35 -p $step3/$prefix.mapped_host.v1 -norm $nReads --plot";
    print "\nSTEP3\n\t $exec3_11\n";
    `$exec3_11`;

    my $nReadsUnmapHost = `grep -c '>' $step4/unmappedVectorReads.fasta`;
    chomp($nReadsUnmapHost);
    my $mapped = $nReads - $nReadsUnmapHost;

    # print metrics "#reads mapped host\t".$mapped. # REVIEW: 2023-03-06 - What should we do with these 'special' logs?
    print "#reads mapped host\t".$mapped.
    "\n";
    # mapped reads HOST
    # print metrics "#reads unmapped host\t".$nReadsUnmapHost. # REVIEW: 2023-03-06 - What should we do with these 'special' logs?
    print "#reads unmapped host\t".$nReadsUnmapHost.
    "\n";
    # unmapped reads HOST

    # print interest "#reads mapped host\t".$mapped. # REVIEW: 2023-03-06 - What should we do with these 'special' logs?
    print "#reads mapped host\t".$mapped.
    "\n";
    # mapped reads HOST
    # print interest "#reads unmapped host\t".$nReadsUnmapHost. # REVIEW: 2023-03-06 - What should we do with these 'special' logs?
    print "#reads unmapped host\t".$nReadsUnmapHost.
    "\n";
    # unmapped reads HOST

    `rm -rf $step3/mapped_host.v1.sam`;
    # deleting sam file mapped reads on host genome

    #Mapping Host - filtered reads against bacters reference
    print "[MAPPING HOST-FILTERED READS AGAINST BACTERIAL GENOMES]... \n";
    my $exec5_1 = "bowtie -f -S -v 1 --un $step4/unmappedVectorBacters.fasta -k 1 -p $process --large-index $path_bacterial_genome_all $step4/unmappedVectorReads.fasta > /dev/null 2>>$prefix.warn ";
    print "\nSTEP5_1\n\t $exec5_1\n";
    `$exec5_1`;

    $nReadsUnmapHostBac = `grep -c '>' $step4/unmappedVectorBacters.fasta`;
    chomp($nReadsUnmapHostBac);
    $mappedbac = $nReadsUnmapHost - $nReadsUnmapHostBac;

    # print metrics "#reads mapped bacter\t".$mappedbac. # REVIEW: 2023-03-06 - What should we do with these 'special' logs?
    print "#reads mapped bacter\t".$mappedbac.
    "\n";
    
    # mapped reads Bacterial genomes
    # print metrics "#preprocessed reads\t".$nReadsUnmapHostBac. # REVIEW: 2023-03-06 - What should we do with these 'special' logs?
    print "#preprocessed reads\t".$nReadsUnmapHostBac.
    "\n";
    # pre - processed reads

    print "\n  PRE-PROCESSING FINISHED \n";
}

# #######################################################################

# 
# TODO: 2023-03-01 - Encapsulate these date / time formatting stuff
# 

# Calculate elapsed time
my $time_elapsed = Time::HiRes::tv_interval([$time_start]);
my $time_elapsed_ms = 1000 * ($time_elapsed - int($time_elapsed));
my $time_elapsed_str = strftime("%H:%M:%S", localtime($time_elapsed)) . sprintf ":%03d", ($time_elapsed_ms);

my $msg_finish = "
-- THE END --
Time elapsed: $time_elapsed_str
";

print $msg_finish;
close($LOG_FH);

select STDOUT;
print "$msg_finish \n";



    
# my $filename = shift or die "Usage: $0 YAML-FILE\n";
# my $data = LoadFile($filename);
# print Dumper $data;
