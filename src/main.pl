use strict;
use warnings;

use Time::HiRes;
use POSIX qw(strftime);
use Getopt::Long;

#
# TODO: 2023-03-07 - Use a decent logger
# 
# use constant PATH_LOG_MAIN => "srna_metavir.main.log";
# use constant EXEC_ROOT_DIR => "/srna_metavir/runs";
use constant EXEC_ROOT_DIR => "./runs";
use constant REF_BACTERIA_GENOMES => "/srna_metavir/asset/refs/bacterial_genomes/all_bacters.fasta";
use constant REF_BLAST_DB_NT => "/srna_metavir/asset/blastdb/nt";
use constant REF_DIAMOND_NR => "/srna_metavir/asset/diamond/nr.dmnd";

# $| = 1;     # forces immediate prints into files rather than the buffer.

#######################################################################
### TIME HANDLERS -----------------------------------------------------
#######################################################################

sub getTimeDiff {
    my $t0 = $_[0] or die "Must provide at least one date for calculation!";
    my $t1 = $_[1];
    my $time_diff = Time::HiRes::tv_interval([$t0], [$t1]);
    return $time_diff;
}

sub getTimeStr {
    my $time = $_[0] // Time::HiRes::gettimeofday();
    my $time_ms = 1000 * ($time - int($time));
    my $time_str = strftime("%H:%M:%S", localtime($time)) . sprintf ":%03d", ($time_ms);
    return $time_str;
}

sub getStepTimebBeginMsg {

    my $title = $_[0] or die "Must provide step title!";
    my $time = getTimeStr(Time::HiRes::gettimeofday());
    
    return "
------------------------------------------------------
>> Begin of step '$title'

At $time...
";
}

sub getStepTimeEndMsg {

    my $title = $_[0] or die "Must provide step title!";
    my $t0 = $_[1] or die "Must provide time 00!";
    my $t1 = $_[2] or die "Must provide time 01!";

    my $t0_str = getTimeStr($t0);
    my $t1_str = getTimeStr($t1);
    my $time_diff_str = getTimeStr(getTimeDiff($t0, $t1));

    return "

>> End of step '$title'...

From '$t0_str' to '$t1_str'
Time elapsed: $time_diff_str

------------------------------------------------------
";
}

my $time_start = $^T;
my $current_time = $time_start;
my $last_time = $time_start;

my $time_diff = 0;
my $time_msg = "";

my $step_name = "";

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

my $exec_id = "exec_test_06";
my $prefix = EXEC_ROOT_DIR . "/$exec_id";

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
my $temp_prefix;

GetOptions("qual=s" => \$qual,
    "hostgenome=s" => \$hostgenome,
    "fasta=s" => \$fasta,
    "fastqgz=s" => \$fastqgz,
    "fastq=s" => \$fastq,
    # "prefix=s" => \$prefix, # TODO: 2023-03-06 - Reenable custom naming for this folder
    "tempprefix=s" => \$temp_prefix, # TODO: 2023-03-06 - Reenable custom naming for this folder
    "size=s" => \$size,
    "hash=s" => \$hash,
    # "log=s" => \$log,
    "si=s" => \$si,
    "se=s" => \$se,
    "process=s" => \$process,
    "clean!" => \$clean,
    "nohostfilter!" => \$nohostfilter, # Options must be lowcase and without "_"
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
### Set running details -----------------------------------------------
#######################################################################

my $time_start_str = getTimeStr($time_start);
my $runningDetails = "

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
    $runningDetails .= "> Searching for degradation: TRUE\n";
} else {
    $runningDetails .= "> Searching for degradation: FALSE\n";
}

$runningDetails .= "\n-------------------------------------------\n";

print STDOUT $runningDetails . "\n";

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
# open($LOG_FH, ">>", PATH_LOG_MAIN) or die "Couldn't open: $!"; # $! is a special variable holding the error
open($LOG_FH, ">>", "$prefix/$exec_id.log") or die "Couldn't open: $!"; # $! is a special variable holding the error
select $LOG_FH;

print "\n\n";
print "#######################################################################\n";
print ">> New execution";
print $runningDetails;

#######################################################################
### Define paths ------------------------------------------------------
#######################################################################

# 
# REVIEW: 2023-02-25 - Should we stop using this 'binary' thing?
# 
# our $binary = "/home/bioinfo/eric_bin";

# 
# REVIEW: 2023-02-27 - Find a better way to manage file paths
# REVIEW: 2023-03-01 - Shall we standardize these directories as all other 'path' variables too? ('path' prefix?)
# 

# Step folders
my $step0		="$prefix/00_saet";
my $step1		="$prefix/01_trimming";
my $step2		="$prefix/02_filter_size_gaps_convertion";
my $step3		="$prefix/03_mapping_vector";
my $step4		="$prefix/04_getUnmapped";
my $step5_fix	="$prefix/05_assembleUnmapped_fix";
my $step5_opt	="$prefix/05_assembleUnmapped_opt";
my $step5_opt_fix="$prefix/05_assembleUnmapped_opt_fix";
my $step5_opt_20to23="$prefix/05_assembleUnmapped_opt_20to23";
my $step5_opt_24to30="$prefix/05_assembleUnmapped_opt_24to30";
my $step5_contigs= "$prefix/05_assembleUnmapped_final";
my $step5_cap3  = "$prefix/05_cap3";
my $step6		="$prefix/06_blast";
my $step7		="$prefix/07_reportBlast";
my $step8		="$prefix/08_completeReport";
my $step9		="$prefix/09_contigs_no_hit";
my $step10		="$prefix/10_pattern";
my $step_virus	="$prefix/virus";

# Utils scripts
my $path_utils = "/srna_metavir/src/utils";

my $path_warns = "$prefix/$exec_id.warn";

my $path_filter_diamond = "$path_utils/filter_diamond.sh";
my $path_filter_fasta_by_size = "$path_utils/filter_fasta_by_size.py";
my $path_plot_dist_per_base_by_reads = "$path_utils/plot-geral-dist-base-reads/plotGeralDistributionPerBaseByReads.pl";
my $path_merge_contigs = "$path_utils/merge-contigs/mergeContigsNew.pl";
my $path_fix_cap3_contigs = "$path_utils/fixIdCap3Contigs.pl";
my $path_filter_blast = "$path_utils/filterblast.pl";
my $path_analyse_contigs_blastn = "$path_utils/analyse-contigs-blastn/analyzeContigsFilterBlastn.pl";
my $path_extract_seqs_no_hit_blast = "$path_utils/extract-seqs/extractSequencesNoHitBlast.pl";

my $path_sam_stats = "$path_utils/samStatistics_v3.pl";
my $path_z_score = "$path_utils/Z-score.bothstrands.pl";
my $path_heatmap_corr = "$path_utils/heatmap_correlation_VISA.R";
my $path_plot_map_data_base_preference = "$path_utils/plot-map-base-preference/plotMappingDataPerBasePreference.pl";
my $path_calc_pattern_sam = "$path_utils/pattern-sam/calcPatternInSamFile.pl";

my $path_trimmed_filtered_gt15 = "$step2/trimmed_filtered_gt15.fasta";

my $path_unmapped_vector_bacters = "$step4/unmappedVectorBacters.fasta";
my $path_unmapped_trimmed_filtered = "$step4/unmapped_trimmed_filtered.fasta";
my $path_unmapped_trimmed_filtered_20_23 = "$step4/unmapped_trimmed_filtered.20-23.fasta";
my $path_unmapped_trimmed_filtered_24_30 = "$step4/unmapped_trimmed_filtered.24-30.fasta";

#######################################################################
### Create step folders -----------------------------------------------
#######################################################################

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

    if (not -e "$step10/plots") {
        `mkdir $step10/plots`;
    }
}

#######################################################################
### Handle FASTQ sequences --------------------------------------------
#######################################################################

$step_name = "Handle FASTQ sequences";

$time_msg = getStepTimebBeginMsg($step_name);
print STDOUT $time_msg;
print $time_msg;

# -----------------------------------------------------------------------

my $path_00_trim_quality = "$step0/trimming.quality.fastq";
my $path_02_trim_filtered_gt15 = "$step2/trimmed_filtered_gt15.fasta";
my $path_02_trim_quality_gt15 = "$step2/trimmed.quality.gt15.fastq";

# 
# TODO: 2023-02-27 - Handle FastQ sequences
# 

# -----------------------------------------------------------------------

$current_time = Time::HiRes::gettimeofday();
$time_msg = getStepTimeEndMsg($step_name, $last_time, $current_time);
$last_time = $current_time;
print STDOUT $time_msg;

#######################################################################
### Handle FASTA sequences --------------------------------------------
#######################################################################

$step_name = "Handle FASTA sequences";

$time_msg = getStepTimebBeginMsg($step_name);
print STDOUT $time_msg;
print $time_msg;

# -----------------------------------------------------------------------

# 
# REVIEW: 2023-03-01 - Shall we think of better names?
# 

my $mappedbac;
my $nReadsUnmapHostBac;
# my $nReadsUnmapHost;

if (defined($fasta)) {
    print "#Loading FASTA file ... \n";

    if (not defined($nohostfilter)) {

        my $cmd = "cp $fasta $path_trimmed_filtered_gt15";

        print "[COPING $fasta TO $path_trimmed_filtered_gt15]\n";
        `$cmd`;
        
        print "\n[STEP 03]\n\t $cmd\n";

    } else {

        # 
        # REVIEW: 2023-03-01 - Check this step numbering
        # TODO: 2023-03-01 - Test it
        # 

        # Mapping Host - unfiltered reads against bacters reference
        print "[MAPPING HOST-UNFILTERED READS AGAINST BACTERIAL GENOMES]... \n";
        
        my $exec5_1 = "bowtie -f -S -v 1 --un $path_unmapped_vector_bacters -k 1 -p $process --large-index ".REF_BACTERIA_GENOMES." $fasta > /dev/null 2>> $step4/reads_mapped_to_bacteria.log";
        
        print "\n[STEP 05.1]\n\t $exec5_1\n";
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

if (not defined($nohostfilter)) {
    if (defined($fastqgz)) {

        # 
        # TODO: 2023-03-06 - Test it
        # 

        # # # dealing with FASTQ TRIMMED files compacted with GZIP

        print "\n\nRunning step 2 [ converting fastq to fasta - fastq_to_fasta ]\n";

        #converting fastq to fasta
        my $exec_fq_2 = "gunzip -dc $fastqgz | fastq_to_fasta -Q 33 -o $step2/trimmed_filtered_gt15.fasta";
        
        print "\n[STEP 02]\n\t $exec_fq_2\n";
        `$exec_fq_2`;
    }

    print "[MAPPING SEQUENCE AGAINST VECTOR]\n";
    my $exec3 = "bowtie $large_index -f -S -k 1 -p $process -v 1 --un $step4/unmappedVectorReads.fasta $hostgenome $step2/trimmed_filtered_gt15.fasta | awk -F'\\t' '{if( \$2 != 4) print \$0}' > $step3/mapped_host.v1.sam  2>mapping_host.stats  ";
    
    print "\n[STEP 03]\n\t $exec3\n";
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
    my $exec3_1 = "perl $path_plot_dist_per_base_by_reads -sam $step3/mapped_host.v1.sam  -s $si -e $se -p $step3/mapped_host.v1.$si.$se -norm $nReads --plot";
    
    print "\n[STEP 03]\n\t $exec3_1\n";
    `$exec3_1`;

    print "[PLOT GERAL DISTRIBUTION READS MAPPED HOST - 15-35nt ]\n";
    my $exec3_11 = "perl $path_plot_dist_per_base_by_reads -sam $step3/mapped_host.v1.sam  -s 15 -e 35 -p $step3/mapped_host.v1 -norm $nReads --plot";
    
    print "\n[STEP 03]\n\t $exec3_11\n";
    `$exec3_11`;

    my $nReadsUnmapHost = `grep -c '>' $step4/unmappedVectorReads.fasta`;
    chomp($nReadsUnmapHost);
    my $mapped = $nReads - $nReadsUnmapHost;

    # print metrics "#reads mapped host\t".$mapped. # REVIEW: 2023-03-06 - What should we do with these 'special' logs?
    print "# Reads mapped host\t".$mapped.
    "\n";
    # mapped reads HOST
    # print metrics "#reads unmapped host\t".$nReadsUnmapHost. # REVIEW: 2023-03-06 - What should we do with these 'special' logs?
    print "# Reads unmapped host\t".$nReadsUnmapHost.
    "\n";
    # unmapped reads HOST

    # print interest "#reads mapped host\t".$mapped. # REVIEW: 2023-03-06 - What should we do with these 'special' logs?
    print "# Reads mapped host\t".$mapped.
    "\n";
    # mapped reads HOST
    # print interest "#reads unmapped host\t".$nReadsUnmapHost. # REVIEW: 2023-03-06 - What should we do with these 'special' logs?
    print "# Reads unmapped host\t".$nReadsUnmapHost.
    "\n";
    # unmapped reads HOST

    # Deleting sam file mapped reads on host genome
    `rm -rf $step3/mapped_host.v1.sam`;

    # Mapping Host - filtered reads against bacters reference
    print "[MAPPING HOST-FILTERED READS AGAINST BACTERIAL GENOMES]... \n";
    my $exec5_1 = "bowtie -f -S -v 1 --un $path_unmapped_vector_bacters -k 1 -p $process --large-index ".REF_BACTERIA_GENOMES." $step4/unmappedVectorReads.fasta > /dev/null 2>>$path_warns ";
    
    print "\n[STEP 05.1]\n\t $exec5_1\n";
    # `$exec5_1`;

    $nReadsUnmapHostBac = `grep -c '>' $path_unmapped_vector_bacters`;
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

# -----------------------------------------------------------------------

$current_time = Time::HiRes::gettimeofday();
$time_msg = getStepTimeEndMsg($step_name, $last_time, $current_time);
$last_time = $current_time;
print STDOUT $time_msg;

#######################################################################
### Select filtered sequences by size ---------------------------------
#######################################################################

$step_name = "Select filtered sequences by size";

$time_msg = getStepTimebBeginMsg($step_name);
print STDOUT $time_msg;
print $time_msg;

# -----------------------------------------------------------------------

# By default we've being using 18 ~ 30 nt as arguments

# 
# TODO: 2023-05-11 - Use 15 ~ 35nt as default
# 

print "[FILTER UNMAPPED SEQUENCES BY SIZE (variable size $si to $se)]\n";

my $exec5 = "python3 $path_filter_fasta_by_size $path_unmapped_vector_bacters $si $se $path_unmapped_trimmed_filtered -t F ";

print "\n[STEP 05]\n\t $exec5\n";
`$exec5`;

print "[FILTER UNMAPPED SEQUENCES BY SIZE (20-23NT)]\n";
my $exec5_1 = "python3 $path_filter_fasta_by_size $path_unmapped_vector_bacters 20 23 $path_unmapped_trimmed_filtered_20_23";
print "\n[STEP 05.1]\n\t $exec5_1\n";
`$exec5_1`;

print "[FILTER UNMAPPED SEQUENCES BY SIZE (24-30NT)]\n";
my $exec5_2 = "python3 $path_filter_fasta_by_size $path_unmapped_vector_bacters 24 30 $path_unmapped_trimmed_filtered_24_30";
print "\n[STEP 05.2]\n\t $exec5_2\n";
`$exec5_2`;

$current_time = Time::HiRes::gettimeofday();
$time_msg = getStepTimeEndMsg($step_name, $last_time, $current_time);
$last_time = $current_time;
print STDOUT $time_msg;

#######################################################################
### Run Velvet optmiser (automatically defined hash) ------------------
#######################################################################

$step_name = "Run Velvet optmiser (automatically defined hash)";

$time_msg = getStepTimebBeginMsg($step_name);
print STDOUT $time_msg;
print $time_msg;

# -----------------------------------------------------------------------

print "\n#[RUNNING VELVET OPTIMIZER]\n";
print "\t#Running step 6 [ Assemble unmapped 21 nt - velvetOptimser.pl ]\n";

my $exec6_1 = "velvetoptimiser --d $step5_opt/run1 --t $process --s 13 --e 19 --f '-short -fasta $path_unmapped_trimmed_filtered' --a $process 2>>$path_warns";

# -----------------------------------------------------------------------

$current_time = Time::HiRes::gettimeofday();
$time_msg = getStepTimeEndMsg($step_name, $last_time, $current_time);
$last_time = $current_time;
print STDOUT $time_msg;
print "\n[STEP 06.1]\n\t $exec6_1\n";
`$exec6_1`;

print "\t#Running step 6_4 [ SPADES ] \n";
my $exec6_4 = "spades -s $path_unmapped_trimmed_filtered --careful --only-assembler -t $process -k 13,15,17,19 -o $step5_opt/run2";

print "\n[STEP 06.4]\n\t $exec6_4\n";
`$exec6_4`;

print "\t#Running step 6_5 [ merge assemblies - mergeContigs.pl ] \n";
my $exec6_5 = "perl $path_merge_contigs -contig1 $step5_opt/run1/contigs.fa -contig2 $step5_opt/run2/scaffolds.fasta -output $step5_opt/contigs.final.fasta ";

print "\n[STEP 06.5]\n\t $exec6_5\n";
`$exec6_5`;

# -----------------------------------------------------------------------

$current_time = Time::HiRes::gettimeofday();
$time_msg = getStepTimeEndMsg($step_name, $last_time, $current_time);
$last_time = $current_time;

print STDOUT $time_msg;
print $time_msg;

#######################################################################
### Running velvet (fixed hash) ---------------------------------------
#######################################################################

$step_name = "Running velvet (fixed hash)";

$time_msg = getStepTimebBeginMsg($step_name);
print STDOUT $time_msg;
print $time_msg;

# -----------------------------------------------------------------------

print "\n[RUNNING DEFAULT VELVET]\n";
print "\t#Running step 6 [ Assemble unmapped 21 nt - velvet hash $hash ]\n";

my $exec6 = "velveth $step5_fix/run1 $hash -fasta -short $path_unmapped_trimmed_filtered 2>>$path_warns";

print "\n[STEP 06]\n\t $exec6\n";
`$exec6`;

`velvetg $step5_fix/run1 -exp_cov auto -cov_cutoff auto 2>$path_warns`;

`mkdir $step5_fix/run2`;

print "\t#Running SPADES fixed hash  [ SPADES ] \n";
my $exec6_2 = "spades -s $path_unmapped_trimmed_filtered --careful --only-assembler -t $process  -k $hash -o $step5_fix/run2 ";

print "\n[STEP 06.2]\n\t $exec6_2\n";
`$exec6_2`;

print "\t#Running step 6_5 [ merge assemblies - mergeContigs.pl ] \n";
$exec6_5 = "perl $path_merge_contigs -contig1 $step5_fix/run1/contigs.fa -contig2 $step5_fix/run2/scaffolds.fasta -output $step5_fix/contigs.final.fasta ";

print "\n[STEP 06.5]\n\t $exec6_5\n";
`$exec6_5`;

# -----------------------------------------------------------------------

$current_time = Time::HiRes::gettimeofday();
$time_msg = getStepTimeEndMsg($step_name, $last_time, $current_time);
$last_time = $current_time;

print STDOUT $time_msg;
print $time_msg;

#######################################################################
### Running velvet optmiser (FIXED hash) ------------------------------
#######################################################################

$step_name = "Running velvet optmiser (FIXED hash)";

$time_msg = getStepTimebBeginMsg($step_name);
print STDOUT $time_msg;
print $time_msg;

# -----------------------------------------------------------------------

print "\n[VELVET OPTIMISER HASH ONLY 15]\n";
print "\t#Running step 6_6 [ Assemble unmapped 21 nt - velvetOptimser.pl ]\n";
my $exec6_6 = "velvetoptimiser --d $step5_opt_fix/run1 --t $process --s $hash --e $hash --f '-short -fasta $path_unmapped_trimmed_filtered' --a 2>>$path_warns";

print "\n[STEP 06.1]\n\t $exec6_6\n";
`$exec6_6`;

`mkdir $step5_opt_fix/run2`;
print "\t#Running SPADES fixed hash  [ SPADES ] \n";
my $exec6_6_2 = "spades -s $path_unmapped_trimmed_filtered --careful --only-assembler -t $process  -k 15 -o $step5_opt_fix/run2 ";

print "\n[STEP 06.2]\n\t $exec6_2\n";
`$exec6_6_2`;

print "\t#Running step 6_9 [ merge assemblies - mergeContigs.pl ] \n";
my $exec6_9 = "perl $path_merge_contigs -contig1 $step5_opt_fix/run1/contigs.fa -contig2 $step5_opt_fix/run2/scaffolds.fasta -output $step5_opt_fix/contigs.final.fasta ";

print "\n[STEP 06.5]\n\t $exec6_5\n";
`$exec6_9`;

# -----------------------------------------------------------------------

$current_time = Time::HiRes::gettimeofday();
$time_msg = getStepTimeEndMsg($step_name, $last_time, $current_time);
$last_time = $current_time;

print STDOUT $time_msg;
print $time_msg;

#######################################################################
### Running velvet optmiser (FIXED hash) 20-23 ------------------------
#######################################################################

$step_name = "Running velvet optmiser (FIXED hash) 20-23";

$time_msg = getStepTimebBeginMsg($step_name);
print STDOUT $time_msg;
print $time_msg;

# -----------------------------------------------------------------------

print "\n[VELVET OPTIMISER HASH ONLY 15 - 20-23nt]\n";
print "\t#Running step 6_10 [ Assemble unmapped 20-23nt nt - velvetOptimser.pl ]\n";
my $exec6_10 = "velvetoptimiser --d $step5_opt_20to23/run1 --t $process --s $hash --e $hash --f '-short -fasta $path_unmapped_trimmed_filtered_20_23' --a 2>>$path_warns";

print "\n[STEP 06.1]\n\t $exec6_10\n";
`$exec6_10`;

`mkdir $step5_opt_20to23/run2`;
print "\t#Running SPADES fixed hash  [ SPADES ] \n";
my $exec6_10_2 = "spades -s $path_unmapped_trimmed_filtered_20_23  --careful --only-assembler -t $process  -k $hash -o $step5_opt_20to23/run2 ";

print "\n[STEP 06.10.2]\n\t $exec6_2\n";
`$exec6_10_2`;

print "\t#Running step 6_13 [ merge assemblies - mergeContigs.pl ] \n";
my $exec6_13 = "perl $path_merge_contigs -contig1 $step5_opt_20to23/run1/contigs.fa -contig2 $step5_opt_20to23/run2/scaffolds.fasta -output $step5_opt_20to23/contigs.final.fasta ";

print "\n[STEP 06.13]\n\t $exec6_13\n";
`$exec6_13`;

# -----------------------------------------------------------------------

$current_time = Time::HiRes::gettimeofday();
$time_msg = getStepTimeEndMsg($step_name, $last_time, $current_time);
$last_time = $current_time;

print STDOUT $time_msg;
print $time_msg;

#######################################################################
### Running velvet optmiser (hash 17) 24-30 ---------------------------
#######################################################################

$step_name = "Running velvet optmiser (hash 17) 24-30";

$time_msg = getStepTimebBeginMsg($step_name);
print STDOUT $time_msg;
print $time_msg;

# -----------------------------------------------------------------------

#
# TODO: 2023-05-31 - Test this condition
# 

if (defined($deg)) {
    
    print "\n[VELVET OPTIMISER - 24-30nt]\n";
    print "\t#Running step 6_10 [ Assemble unmapped 24-30nt nt - velvetOptimser.pl ]\n";
    
    my $exec62_10 = "velvetoptimiser --d $step5_opt_24to30/run1 --t $process --s 15 --e 17 --f '-short -fasta $path_unmapped_trimmed_filtered_24_30' --a 2>>$path_warns";
    
    print "\nSTEP62_10\n\t $exec62_10\n";
    `$exec62_10`;

    `mkdir $step5_opt_24to30/run2 `;

    print "\t#Running SPADES fixed hash  [ SPADES ] \n";
    my $exec6_10_3 = "spades -s $path_unmapped_trimmed_filtered_24_30 --careful --only-assembler -t $process  -k 15,17 -o $step5_opt_24to30/run2 ";
    
    print "\n[STEP 06.10_3\n\t $exec6_2\n";
    `$exec6_10_3`;

    print "\t#Running step 6_13 [ merge assemblies - mergeContigs.pl ] \n";
    my $exec62_13 = "perl $path_merge_contigs -contig1 $step5_opt_24to30/run1/contigs.fa -contig2 $step5_opt_24to30/run2/scaffolds.fasta -output $step5_opt_24to30/contigs.final.fasta ";
    print "\nSTEP62_13\n\t $exec62_13\n";
    `$exec62_13`;
}

# -----------------------------------------------------------------------

$current_time = Time::HiRes::gettimeofday();
$time_msg = getStepTimeEndMsg($step_name, $last_time, $current_time);
$last_time = $current_time;

print STDOUT $time_msg;
print $time_msg;

#######################################################################
### Merging assemblies ------------------------------------------------
#######################################################################

$step_name = "Merging assemblies";

$time_msg = getStepTimebBeginMsg($step_name);
print STDOUT $time_msg;
print $time_msg;

# -----------------------------------------------------------------------

my $path_05_cap3_gt200 = "$step5_cap3/contigs_merged.final.gt200.fasta";

if (not defined($deg)) {
    `touch $step5_opt_24to30/contigs.final.fasta`;
}

print "\n[MERGING CONTIGS AND RUNNING CAP3]\n";
print "\t#Running STEP [CAT] Concatenating contigs...\n";
my $exec_cat = "cat $step5_opt_20to23/contigs.final.fasta $step5_opt_fix/contigs.final.fasta  $step5_fix/contigs.final.fasta $step5_opt/contigs.final.fasta $step5_opt_24to30/contigs.final.fasta > $step5_cap3/all_contigs.fasta";

print "\n[STEP 06.CAT] \n\t $exec_cat\n";
`$exec_cat`;

print "\t#Running step CAP3 [ assembly contigs from different assemblies ] \n";
my $exec_cap3 = "cap3 $step5_cap3/all_contigs.fasta  > $step5_cap3/log.cap3";
print "\nSTEP CAP3\n\t $exec_cap3\n";
`$exec_cap3`;

print "\t#Running [ Concatenning contigs and singlets from CAP3]\n";
my $exec_cap3_1 = "cat $step5_cap3/all_contigs.fasta.cap.contigs $step5_cap3/all_contigs.fasta.cap.singlets > $step5_cap3/contigs_merged.final.fasta";

print "\n[STEP 06.CAT]\n\t $exec_cap3_1\n";
`$exec_cap3_1`;

print "\t#Running step 6_7 [ selecting contigs larger than n50 - calcN50.pl ] \n";
my $step6_7 = "perl $path_fix_cap3_contigs -i $step5_cap3/contigs_merged.final.fasta -s 50  -p $step5_cap3/contigs_merged.final";

print "\n[STEP 06.7]\n\t $step6_7\n";
`$step6_7`;

my $countAssembledContigs = `grep -c '>' $step5_cap3/contigs_merged.final.gt50.fasta`;
chomp($countAssembledContigs);
# print metrics "#total assembled contigs\t" . $countAssembledContigs. "\n"; # Assembled Contigs
print "# Total assembled contigs\t" . $countAssembledContigs. "\n"; # Assembled Contigs
# print interest "#total assembled contigs\t" . $countAssembledContigs . "\n"; # Assembled Contigs
print "# Total assembled contigs\t" . $countAssembledContigs . "\n"; # Assembled Contigs

print "[FILTER CONTIGS gt 200 nt]\n";
my $exec_FC2 = "python3 $path_filter_fasta_by_size $step5_cap3/contigs_merged.final.gt50.fasta 200 1000000 $path_05_cap3_gt200";

print "\n[STEP 05.2]\n\t $exec_FC2\n";
`$exec_FC2`;

$countAssembledContigs = `grep -c '>' $path_05_cap3_gt200`;
chomp($countAssembledContigs);
# print metrics "#contigs gt200\t".$countAssembledContigs. "\n"; # Assembled Contigs
print "# Contigs gt200\t".$countAssembledContigs. "\n"; # Assembled Contigs
# print interest "#contigs gt200\t".$countAssembledContigs . "\n"; # Assembled Contigs
print "# Contigs gt200\t".$countAssembledContigs . "\n"; # Assembled Contigs

# -----------------------------------------------------------------------

$current_time = Time::HiRes::gettimeofday();
$time_msg = getStepTimeEndMsg($step_name, $last_time, $current_time);
$last_time = $current_time;

print STDOUT $time_msg;
print $time_msg;

#######################################################################
### Blastn ------------------------------------------------------------
#######################################################################

$step_name = "Blastn";

$time_msg = getStepTimebBeginMsg($step_name);
print STDOUT $time_msg;
print $time_msg;

# -----------------------------------------------------------------------

my $path_06_blastn_no_hit_1e5 = "$step6/seqNoHit.blastN.1e5.fasta";
my $path_06_blastn_merge_gt200_1e5 = "$step6/contigs_merged.final.gt200.1e5.blastn";

my $path_07_bN_analyse_blastn_prefix = "$step7/contigs.bN.blastn.analyze";
my $path_07_bN_analyze_prefix = "$step7/contigs.bN.analyze";
my $path_07_bN_analize = "$path_07_bN_analyze_prefix..contigs.fasta"; # "$step7/contigs.bN.analyze..contigs.fasta"
my $path_07_bN_analize_viral = "$path_07_bN_analyze_prefix.virus.contigs.fasta"; # "$step7/contigs.bN.blastn.analyze.virus.contigs.fasta"
my $path_07_bN_analize_non_viral = "$path_07_bN_analyze_prefix.nonviral.contigs.fasta"; # "$step7/contigs.bN.analyze.nonviral.contigs.fasta";

my $path_07_blastn_analyse = "$step7/contigs.blastN.analyze";
my $path_07_blastn_analize_virus = "$step7/contigs.blastN.virus.analyze";
my $path_07_blastn_1e5_report = "$step7/contigs.blastn.1e5.report";
my $path_07_blastn_virus = "$step7/contigs.virus.blastN.formatted.fasta";
my $path_07_blastn_virus_header = "$step7/contigs.virus.blastN.formatted.header.fasta";

my $path_07_dmnd_log = "$step7/diamond.log";
my $path_07_dmnd_hits = "$step7/diamond_blastx_Hits.fasta";
my $path_07_dmnd_out = "$step7/diamond_blastx.out";
my $path_07_dmnd_no_hits = "$step7/diamond_blastx_NoHits.fasta";
my $path_07_dmnd_no_hits_linear = "$step7/diamond_blastx_NoHits_linear.fasta";
my $path_07_dmnd_no_hits_header = "$step7/diamond_blastx_NoHits_linear.header.fasta";
my $path_07_dmnd_viral_header = "$step7/diamond_blastx_Viral.header.fasta";

my $path_07_aux_fasta = "$step7/aux.fasta";
my $path_07_aux_non_viral = "$step7/aux_nonviral";

print "\n[BlastN contigs gt 200]\n";
print "\t#Running step 10_111 [ Blast against NT - blast+ blastn ]\n";
my $exec10_111 = "blastn -query $path_05_cap3_gt200 -db ".REF_BLAST_DB_NT."/nt -num_descriptions 5 -num_alignments 5 -evalue 1e-5 -out $path_06_blastn_merge_gt200_1e5  -num_threads $process";

print "\nSTEP10_111\n\t $exec10_111\n";
`$exec10_111`;

print "\t#Running filterblast...\n";
my $exec10_112 = "perl $path_filter_blast -b $path_06_blastn_merge_gt200_1e5 -evalue 1e-5 --best --desc > $path_07_blastn_1e5_report";

print "\nSTEP10_112\n\t $exec10_112\n";
`$exec10_112`;

print "\n[Extracting contigs all Hits blastn 1e-5]\n";
`perl $path_analyse_contigs_blastn -i $path_07_blastn_1e5_report  -f $path_05_cap3_gt200 -q "" -p  $path_07_bN_analyze_prefix --fasta > $path_07_blastn_analyse`;

print "\n[Extracting contigs viral blastn 1e-5]\n";
`perl $path_analyse_contigs_blastn -i $path_07_blastn_1e5_report  -f $path_05_cap3_gt200 -q "virus" -p  $path_07_bN_analyse_blastn_prefix --fasta > $path_07_blastn_analize_virus`;

print "\n[Extracting contigs nonviral blastn 1e-5]\n";
`grep -v -i "virus" $path_07_bN_analize | grep '>' | cut -f1 -d " " >  $path_07_aux_non_viral`;

`fasta_formatter -i $path_07_bN_analize > $path_07_aux_fasta`;
`while read p; do grep -A1 \${p} $path_07_aux_fasta >> $path_07_bN_analize_non_viral;done < $path_07_aux_non_viral`;
`rm -f $path_07_aux_fasta`;
`rm -f $path_07_aux_non_viral`;

my $hitsBlastn = `grep -c '>' $path_07_bN_analize`;
chomp($hitsBlastn);

# 
# TODO: 2023-05-23 - Check what to do with these 'parallel' logging files
# 

# Assembled Contigs
# print metrics "#contigs hit blastN\t".$hitsBlastn."\n";
# print interest "#contigs hit blastN\t".$hitsBlastn."\n";
print "# Contigs hit blastN\t".$hitsBlastn."\n";

# Assembled Contigs
my $hitsVirusBlastn = `grep -c '>' $path_07_bN_analize_viral`;
chomp($hitsVirusBlastn);

# print metrics "#contigs hit VIRUS blastN\t".$hitsVirusBlastn."\n";
# print interest "#contigs hit VIRUS blastN\t".$hitsVirusBlastn."\n";
print "# Contigs hit VIRUS blastN\t".$hitsVirusBlastn."\n";

# Assembled Contigs
print "\n[Extracting contigs no hit blastn 1e-5]\n";
my $exec10_113 = "perl $path_extract_seqs_no_hit_blast -seq $path_05_cap3_gt200 -blast $path_07_blastn_1e5_report -out $path_06_blastn_no_hit_1e5";

print "\nSTEP10_113\n\t $exec10_113\n";
`$exec10_113`;

my $seqsNoHitBlastn = `grep -c '>' $path_06_blastn_no_hit_1e5`;
chomp($seqsNoHitBlastn);

# 
# TODO: 2023-05-23 - Check what to do with these 'parallel' logging files
# 

# print metrics "#contigs not hit blastN\t".$seqsNoHitBlastn."\n";
print "#contigs not hit blastN\t".$seqsNoHitBlastn."\n";
# Assembled Contigs
`cat $path_07_bN_analize_viral | perl -pi -e 's/>(\\S+) (\\S+) (\\S+) (\\S+).+/>blastN_\$1_\$2_\$3_\$4/g' > $path_07_blastn_virus`;

# -----------------------------------------------------------------------

$current_time = Time::HiRes::gettimeofday();
$time_msg = getStepTimeEndMsg($step_name, $last_time, $current_time);
$last_time = $current_time;

print STDOUT $time_msg;
print $time_msg;

#######################################################################
### DIAMOND (Blastx) --------------------------------------------------
#######################################################################

$step_name = "DIAMOND (Blastx)";

$time_msg = getStepTimebBeginMsg($step_name);
print STDOUT $time_msg;
print $time_msg;

# -----------------------------------------------------------------------

print "\n[Diamond (BlastX) contigs gt 200]\n";
print "\t# Running step 9 [ Diamond-Blast against NR ]\n";

my $exec9 = "diamond blastx -q $path_06_blastn_no_hit_1e5 -d ".REF_DIAMOND_NR." -k 5 -p $process -e 0.001 -f 0 -c 1 -b 20 --very-sensitive -o $path_07_dmnd_out --un $path_07_dmnd_no_hits --unfmt fasta --al $path_07_dmnd_hits --alfmt fasta 2>  $path_07_dmnd_log";

print "\nSTEP9\n\t $exec9\n";
`$exec9`;

print "\t# Filtering Diamond results... ]\n";   
my $exec9_11 = "$path_filter_diamond -f $path_07_dmnd_hits -o $path_07_dmnd_out -d $step7";

print "\nSTEP9_11\n\t $exec9_11\n";
`$exec9_11`;

# Replace '\t' characters with actual tabs
`sed -i 's/\\\\t/\\t/g' $path_07_dmnd_no_hits`;

# Make a 'linear' version of 'no hits' fasta (join all pieces for each sequence)
`fasta_formatter -i $path_07_dmnd_no_hits -o $path_07_dmnd_no_hits_linear`;
# `rm $path_07_dmnd_no_hits`;

# Add prefix to all contig names
`for file in $step7/*.fasta; do
    if [ -f "\$file" ]; then
        filename=\$(basename "\$file")
        filename_without_ext="\${filename%.*}"
        sed 's/^>/>${temp_prefix}_/' $step7/\$filename > $step7/\${filename_without_ext}.header.fasta
    fi
done
`;

# -----------------------------------------------------------------------

$current_time = Time::HiRes::gettimeofday();
$time_msg = getStepTimeEndMsg($step_name, $last_time, $current_time);
$last_time = $current_time;

print STDOUT $time_msg;
print $time_msg;

#######################################################################
### Merge sequences viral hits (blastn & diamond) & nohits ------------
#######################################################################

$step_name = "Merge sequences viral hits (blastn & diamond) & nohits";

$time_msg = getStepTimebBeginMsg($step_name);
print STDOUT $time_msg;
print $time_msg;

# -----------------------------------------------------------------------

`cat $path_07_blastn_virus_header $path_07_dmnd_viral_header $path_07_dmnd_no_hits_header > $step9/seq_ViralHits_and_NoHits.fasta`;

`bowtie-build $step9/seq_ViralHits_and_NoHits.fasta $step9/seq_ViralHits_and_NoHits.fasta`;

print "\t#Mapping reads against viral hits(blastn and diamond) and nohits\n";
my $exec10_13 = "bowtie -f -S -p $process -v 1 $step9/seq_ViralHits_and_NoHits.fasta $step4/unmappedVectorBacters.fasta 2>>$path_warns | awk -F'\t' '{if( \$2 != 4) print \$0}' > $step9/reads_vs_contigsHitNoHit.v1.sam ";

print"\nSTEP10_13\n\t $exec10_13\n";
`$exec10_13`;

`perl $path_sam_stats -sam $step9/reads_vs_contigsHitNoHit.v1.sam -fa $step9/seq_ViralHits_and_NoHits.fasta -p $step9/seq_ViralHits_and_NoHits --profile --counts`;

`perl $path_z_score -sam $step9/reads_vs_contigsHitNoHit.v1.sam -p $step9/reads_vs_contigsHitNoHit`;

`R --no-save $step9/reads_vs_contigsHitNoHit.zscore.tab $step9/reads_vs_contigsHitNoHit.zscore  < $path_heatmap_corr 2>/dev/null`;

# Stats
$hitsBlastn = `grep '>' $step9/seq_ViralHits_and_NoHits.withSiRNA.fasta | grep blastN | wc -l `;
chomp($hitsBlastn);

# 
# TODO: 2023-05-25 - Check what to do with all these 'parallel' logging files
# 

# print metrics "#contigs hit VIRUS blastN with siRNA\t".$hitsBlastn."\n";
# print interest "#contigs hit VIRUS blastN with siRNA\t".$hitsBlastn."\n";
print "# Contigs hit VIRUS blastN with siRNA\t".$hitsBlastn."\n";
$hitsBlastn = `grep '>' $step9/seq_ViralHits_and_NoHits.withSiRNA_and_PiRNA.fasta | grep blastN | wc -l `;
chomp($hitsBlastn);

# print metrics "#contigs hit VIRUS blastN with siRNA and piRNA\t".$hitsBlastn."\n";
# print interest "#contigs hit VIRUS blastN with siRNA and piRNA\t".$hitsBlastn."\n";
print "# Contigs hit VIRUS blastN with siRNA and piRNA\t".$hitsBlastn."\n";
$hitsBlastn = `grep '>' $step9/seq_ViralHits_and_NoHits.withPiRNA.fasta | grep blastN | wc -l `;
chomp($hitsBlastn);

# print metrics "#contigs hit VIRUS blastN with piRNA\t".$hitsBlastn."\n";
# print interest "#contigs hit VIRUS blastN with piRNA\t".$hitsBlastn. "\n";
print "# Contigs hit VIRUS blastN with piRNA\t".$hitsBlastn."\n";
$hitsVirusBlastn = `grep  '>' $step9/seq_ViralHits_and_NoHits.withSiRNA.fasta | grep blastX | wc -l `;
chomp($hitsVirusBlastn);

# print metrics "#contigs hit VIRUS BlastX with siRNA \t".$hitsVirusBlastn."\n";
# print interest "#contigs hit VIRUS BlastX with siRNA \t".$hitsVirusBlastn."\n";
print "# Contigs hit VIRUS BlastX with siRNA \t".$hitsVirusBlastn."\n";
$hitsVirusBlastn = `grep  '>' $step9/seq_ViralHits_and_NoHits.withSiRNA_and_PiRNA.fasta| grep blastX | wc -l `;
chomp($hitsVirusBlastn);

# print metrics "#contigs hit VIRUS BlastX with siRNA and piRNA \t".$hitsVirusBlastn."\n";
# print interest "#contigs hit VIRUS BlastX with siRNA and piRNA \t".$hitsVirusBlastn."\n";
print "# Contigs hit VIRUS BlastX with siRNA and piRNA \t".$hitsVirusBlastn."\n";
$hitsVirusBlastn = `grep  '>' $step9/seq_ViralHits_and_NoHits.withPiRNA.fasta | grep blastX | wc -l `;
chomp($hitsVirusBlastn);

# print metrics "#contigs hit VIRUS BlastX with piRNA \t".$hitsVirusBlastn."\n";
# print interest "#contigs hit VIRUS BlastX with piRNA \t".$hitsVirusBlastn."\n";
print "# Contigs hit VIRUS BlastX with piRNA \t".$hitsVirusBlastn."\n";
$hitsVirusBlastn = `grep  '>' $step9/seq_ViralHits_and_NoHits.withSiRNA.fasta | grep -v blastX | grep -v blastN | wc -l `;
chomp($hitsVirusBlastn);

# print metrics "#contigs NOHIT blast with siRNA \t".$hitsVirusBlastn."\n";
# print interest "#contigs NOHIT blast with siRNA \t".$hitsVirusBlastn."\n";
print "# Contigs NOHIT blast with siRNA \t".$hitsVirusBlastn."\n";
$hitsVirusBlastn = `grep  '>' $step9/seq_ViralHits_and_NoHits.withSiRNA_and_PiRNA.fasta | grep -v blastN | -v grep blastX | wc -l `;
chomp($hitsVirusBlastn);

# print metrics "#contigs NOHIT blast with siRNA and piRNA \t".$hitsVirusBlastn."\n";
# print interest "#contigs NOHIT blast with siRNA and piRNA \t".$hitsVirusBlastn."\n";
print "# Contigs NOHIT blast with siRNA and piRNA \t".$hitsVirusBlastn."\n";
$hitsVirusBlastn = `grep  '>' $step9/seq_ViralHits_and_NoHits.withPiRNA.fasta | grep -v blastN | -v grep blastX | wc -l `;
chomp($hitsVirusBlastn);

# print metrics "#contigs NOHIT blast with piRNA \t".$hitsVirusBlastn."\n";
# print interest "#contigs NOHIT blast with piRNA \t".$hitsVirusBlastn."\n";
print "# Contigs NOHIT blast with piRNA \t".$hitsVirusBlastn."\n";

# close(metrics);

# -----------------------------------------------------------------------

$current_time = Time::HiRes::gettimeofday();
$time_msg = getStepTimeEndMsg($step_name, $last_time, $current_time);
$last_time = $current_time;

print STDOUT $time_msg;
print $time_msg;

#######################################################################
### Pattern based analysis --------------------------------------------
#######################################################################

$step_name = "Pattern based analysis";

$time_msg = getStepTimebBeginMsg($step_name);
print STDOUT $time_msg;
print $time_msg;

# -----------------------------------------------------------------------

`mkdir $step_virus`;

# Merge 02 files with contigs hit virus
my $path_contig_viral_hits_all = "$step10/all_contigs_hit_virus.fasta";

`cat $path_07_blastn_virus_header $path_07_dmnd_viral_header > $path_contig_viral_hits_all`;

my $count_hit_lines = `cat $path_contig_viral_hits_all | wc -l`;
chomp($count_hit_lines);

# if ($count_hit_lines > 0) {
#     `bowtie-build $path_contig_viral_hits_all $path_contig_viral_hits_all`;

# } else {
#     print "\nNo contig viral hits found...\n";
# }

`bowtie-build $path_contig_viral_hits_all $path_contig_viral_hits_all`;

`bowtie -f -S -k 1 -p $process -v 1  $path_contig_viral_hits_all $step2/trimmed_filtered_gt15.fasta | awk -F '\t' '{if( \$2 != 4) print \$0}'  > $step10/reads.VS.contigs_hit_virus_blast.sam`;

`perl $path_sam_stats -sam $step10/reads.VS.contigs_hit_virus_blast.sam -fa $path_contig_viral_hits_all -p $step10/reads.VS.contigs_hit_virus_blast --profile`;

# Creating reference table with identified contigs hit virus by sequence similarity
`grep ">" $path_contig_viral_hits_all | cut -f 2 -d '>'  > $step10/all_contigs_hit_virus_sequence_similarity.tab`;

print "\n\n Calculating pattern viral contigs and candidates - HEATMAP \n\n";

# Bowtie

# 
# TODO: 2023-05-30 - Check the deprecation warning for passing index via positional args
# 

`bowtie -f -S -k 1 -p $process -v 1  $step9/seq_ViralHits_and_NoHits.fasta $step2/trimmed_filtered_gt15.fasta | awk -F '\t' '{if( \$2 != 4) print \$0}'  > $step10/reads.VS.contigs_virus_and_nohit.sam`;

`perl $path_z_score -sam $step10/reads.VS.contigs_virus_and_nohit.sam -p $step10/reads.VS.contigs_virus_and_nohit`;

`R --no-save $step10/reads.VS.contigs_virus_and_nohit.zscore.tab $step10/plots/reads.VS.contigs_virus_and_nohit.all < $path_heatmap_corr 2>/dev/null`;

# Identifying contigs with siRNA and piRNA signature
`perl $path_sam_stats -sam $step10/reads.VS.contigs_virus_and_nohit.sam -fa $step9/seq_ViralHits_and_NoHits.fasta -p $step10/reads.VS.contigs_virus_and_nohit --profile`;

# Formatting candidate contigs for bowtie
`bowtie-build $step10/reads.VS.contigs_virus_and_nohit.withSiRNA_and_PiRNA.fasta $step10/reads.VS.contigs_virus_and_nohit.withSiRNA_and_PiRNA.fasta > /dev/null `;

`bowtie-build $step10/reads.VS.contigs_virus_and_nohit.withSiRNA.fasta $step10/reads.VS.contigs_virus_and_nohit.withSiRNA.fasta > /dev/null `;

# Bowtie

# 
# TODO: 2023-05-30 - Check the deprecation warning for passing index via positional args
# 

`bowtie -f -S -k 1 -p $process -v 1 $step10/reads.VS.contigs_virus_and_nohit.withSiRNA.fasta $step2/trimmed_filtered_gt15.fasta | awk -F '\t' '{if( \$2 != 4) print \$0}'  > $step10/reads.VS.contigs_virus_and_nohit.siRNAs.sam`;

`bowtie -f -S -k 1 -p $process -v 1  $step10/reads.VS.contigs_virus_and_nohit.withSiRNA_and_PiRNA.fasta $step2/trimmed_filtered_gt15.fasta | awk -F '\t' '{if( \$2 != 4) print \$0}'  > $step10/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs.sam`;

print "\n\n Creating plots \n\n";

# Plotting distribution and density plots
`perl $path_plot_map_data_base_preference -sam $step10/reads.VS.contigs_virus_and_nohit.siRNAs.sam -s $si -e $se -fa $step10/reads.VS.contigs_virus_and_nohit.withSiRNA.fasta -pace 1 -p $step10/plots/reads.VS.contigs_virus_and_nohit.withSiRNA --profile --pattern`;

`perl $path_plot_map_data_base_preference -sam $step10/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs.sam -s $si -e $se -fa $step10/reads.VS.contigs_virus_and_nohit.withSiRNA_and_PiRNA.fasta -pace 1 -p $step10/plots/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs --profile --pattern`;

# Calculating mean pattern
`perl $path_calc_pattern_sam -s $step10/reads.VS.contigs_virus_and_nohit.siRNAs.sam -o $step10/plots/reads.VS.contigs_virus_and_nohit.siRNAs.pattern`;

`perl $path_calc_pattern_sam -s $step10/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs.sam -o $step10/plots/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs.pattern`;

# Plot geral distribution of reads in contigs with siRNA and siRNA+piRNAs
`perl $path_plot_dist_per_base_by_reads -sam $step10/reads.VS.contigs_virus_and_nohit.siRNAs.sam -s $si -e $se -p $step10/plots/reads.VS.contigs_virus_and_nohit.siRNAs.geral_distribution --plot`;

`perl $path_plot_dist_per_base_by_reads -sam $step10/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs.sam -s $si -e $se -p $step10/plots/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs.geral_distribution --plot`;

# Generating clusterized heatmaps
`perl $path_z_score -sam $step10/reads.VS.contigs_virus_and_nohit.siRNAs.sam -p $step10/reads.VS.contigs_virus_and_nohit.siRNAs`;

`perl $path_z_score -sam $step10/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs.sam -p $step10/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs`;

`R --no-save $step10/reads.VS.contigs_virus_and_nohit.siRNAs.zscore.tab $step10/plots/reads.VS.contigs_virus_and_nohit.siRNAs < $path_heatmap_corr 2>/dev/null`;

`R --no-save $step10/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs.zscore.tab $step10/plots/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs < $path_heatmap_corr 2>/dev/null`;

# -----------------------------------------------------------------------

$current_time = Time::HiRes::gettimeofday();
$time_msg = getStepTimeEndMsg($step_name, $last_time, $current_time);
$last_time = $current_time;

print STDOUT $time_msg;
print $time_msg;

#######################################################################

$current_time = Time::HiRes::gettimeofday();
$time_diff = getTimeDiff($time_start, $current_time);
my $time_elapsed_str = getTimeStr($time_diff);

my $msg_finish = "

-- THE END --
Time elapsed: $time_elapsed_str

";

print $msg_finish;
close($LOG_FH);

print STDOUT "$msg_finish \n";
