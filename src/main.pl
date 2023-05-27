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
my $prefix = "exec_20230511_120529";

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

# Utils scripts
my $path_utils = "/srna_metavir/src/utils";
my $path_assets = "/srna_metavir/asset";

my $path_bacterial_genome_all = "$path_assets/refs/bacterial_genomes/all_bacters.fasta";
my $path_blast_db_nt = "$path_assets/blastdb/nt";
my $path_diamond_nr = "$path_assets/diamond/nr.dmnd";

my $path_filter_diamond = "$path_utils/filter_diamond.sh";
my $path_filter_fasta_by_size = "$path_utils/filter_fasta_by_size.py";
my $path_plot_dist_per_base_by_reads = "$path_utils/plot-geral-dist-base-reads/plotGeralDistributionPerBaseByReads.pl";
my $path_merge_contigs = "$path_utils/merge-contigs/mergeContigsNew.pl";
my $path_fix_cap3_contigs = "$path_utils/fixIdCap3Contigs.pl";
my $path_filter_blast = "$path_utils/filterblast.pl";
my $path_analyse_contigs_blastn = "$path_utils/analyse-contigs-blastn/analyzeContigsFilterBlastn.pl";
my $path_extract_seqs_no_hit_blast = "$path_utils/extract-seqs/extractSequencesNoHitBlast.pl";

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
}

#######################################################################
### Handle FASTQ sequences --------------------------------------------
#######################################################################

# 
# TODO: 2023-02-27 - Handle FastQ sequences
# 

#######################################################################
### Handle FASTA sequences --------------------------------------------
#######################################################################

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
        
        print "\nSTEP3\n\t $cmd\n";

    } else {

        # 
        # REVIEW: 2023-03-01 - Check this step numbering
        # TODO: 2023-03-01 - Test it
        # 

        # Mapping Host - unfiltered reads against bacters reference
        print "[MAPPING HOST-UNFILTERED READS AGAINST BACTERIAL GENOMES]... \n";
        

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
    my $exec5_1 = "bowtie -f -S -v 1 --un $path_unmapped_vector_bacters -k 1 -p $process --large-index $path_bacterial_genome_all $step4/unmappedVectorReads.fasta > /dev/null 2>>$prefix.warn ";
    print "\nSTEP5_1\n\t $exec5_1\n";
    `$exec5_1`;

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

#######################################################################
### Select filtered sequences by size ---------------------------------
#######################################################################

# By default we've being using 18 ~ 30 nt as arguments

# 
# TODO: 2023-05-11 - Use 15 ~ 35nt as default
# 
print "[FILTER UNMAPPED SEQUENCES BY SIZE (variable size $si to $se)]\n";


my $exec5 = "python3 $path_filter_fasta_by_size $path_unmapped_vector_bacters $si $se $path_unmapped_trimmed_filtered -t F ";

print "\nSTEP5\n\t $exec5\n";
`$exec5`;

print "[FILTER UNMAPPED SEQUENCES BY SIZE (20-23NT)]\n";
my $exec5_1 = "python3 $path_filter_fasta_by_size $path_unmapped_vector_bacters 20 23 $path_unmapped_trimmed_filtered_20_23";
print "\nSTEP5_1\n\t $exec5_1\n";
`$exec5_1`;

print "[FILTER UNMAPPED SEQUENCES BY SIZE (24-30NT)]\n";
my $exec5_2 = "python3 $path_filter_fasta_by_size $path_unmapped_vector_bacters 24 30 $path_unmapped_trimmed_filtered_24_30";
print "\nSTEP5_2\n\t $exec5_2\n";
`$exec5_2`;

#######################################################################
### Run Velvet optmiser (automatically defined hash) ------------------
#######################################################################

print "\n#[RUNNING VELVET OPTIMIZER]\n";
print "\t#Running step 6 [ Assemble unmapped 21 nt - velvetOptimser.pl ]\n";
my $exec6_1 = "velvetoptimiser --d $step5_opt/run1 --t $process --s 13 --e 19 --f '-short -fasta $path_unmapped_trimmed_filtered' --a $process 2>>$prefix.warn";
print "\nSTEP6_1\n\t $exec6_1\n";
`$exec6_1`;

print "\t#Running step 6_4 [ SPADES ] \n";
my $exec6_4 = "spades -s $path_unmapped_trimmed_filtered --careful --only-assembler -t $process -k 13,15,17,19 -o $step5_opt/run2";
print "\nSTEP6_4\n\t $exec6_4\n";
`$exec6_4`;

print "\t#Running step 6_5 [ merge assemblies - mergeContigs.pl ] \n";
my $exec6_5 = "$path_merge_contigs -contig1 $step5_opt/run1/contigs.fa -contig2 $step5_opt/run2/scaffolds.fasta -output $step5_opt/contigs.final.fasta ";
print "\nSTEP6_5\n\t $exec6_5\n";
`$exec6_5`;

#######################################################################
### Running velvet (fixed hash) ---------------------------------------
#######################################################################

print "\n[RUNNING DEFAULT VELVET]\n";
print "\t#Running step 6 [ Assemble unmapped 21 nt - velvet hash $hash ]\n";
my $exec6 = "velveth $step5_fix/run1 $hash -fasta -short $path_unmapped_trimmed_filtered 2>>$prefix.warn";
print "\nSTEP6\n\t $exec6\n";
`$exec6`;

`velvetg $step5_fix/run1 -exp_cov auto -cov_cutoff auto 2>$prefix.warn`;

`mkdir $step5_fix/run2`;

print "\t#Running SPADES fixed hash  [ SPADES ] \n";
my $exec6_2 = "spades -s $path_unmapped_trimmed_filtered --careful --only-assembler -t $process  -k $hash -o $step5_fix/run2 ";
print "\nSTEP6_2\n\t $exec6_2\n";
`$exec6_2`;

print "\t#Running step 6_5 [ merge assemblies - mergeContigs.pl ] \n";
$exec6_5 = "$path_merge_contigs -contig1 $step5_fix/run1/contigs.fa -contig2 $step5_fix/run2/scaffolds.fasta -output $step5_fix/contigs.final.fasta ";
print "\nSTEP6_5\n\t $exec6_5\n";
`$exec6_5`;

# #######################################################################
# ### Running velvet optmiser (FIXED hash) ------------------------------
# #######################################################################

print "\n[VELVET OPTIMISER HASH ONLY 15]\n";
print "\t#Running step 6_6 [ Assemble unmapped 21 nt - velvetOptimser.pl ]\n";
my $exec6_6 = "velvetoptimiser --d $step5_opt_fix/run1 --t $process --s $hash --e $hash --f '-short -fasta $path_unmapped_trimmed_filtered' --a 2>>$prefix.warn";
print "\nSTEP6_1\n\t $exec6_6\n";
`$exec6_6`;

`mkdir $step5_opt_fix/run2`;
print "\t#Running SPADES fixed hash  [ SPADES ] \n";
my $exec6_6_2 = "spades -s $path_unmapped_trimmed_filtered --careful --only-assembler -t $process  -k 15 -o $step5_opt_fix/run2 ";
print "\nSTEP6_2\n\t $exec6_2\n";
`$exec6_6_2`;

print "\t#Running step 6_9 [ merge assemblies - mergeContigs.pl ] \n";
my $exec6_9 = "$path_merge_contigs -contig1 $step5_opt_fix/run1/contigs.fa -contig2 $step5_opt_fix/run2/scaffolds.fasta -output $step5_opt_fix/contigs.final.fasta ";
print "\nSTEP6_5\n\t $exec6_5\n";
`$exec6_9`;

#######################################################################
### Running velvet optmiser (FIXED hash) 20-23 ------------------------
#######################################################################

print "\n[VELVET OPTIMISER HASH ONLY 15 - 20-23nt]\n";
print "\t#Running step 6_10 [ Assemble unmapped 20-23nt nt - velvetOptimser.pl ]\n";
my $exec6_10 = "velvetoptimiser --d $step5_opt_20to23/run1 --t $process --s $hash --e $hash --f '-short -fasta $path_unmapped_trimmed_filtered_20_23' --a 2>>$prefix.warn";
print "\nSTEP6_1\n\t $exec6_10\n";
`$exec6_10`;

`mkdir $step5_opt_20to23/run2`;
print "\t#Running SPADES fixed hash  [ SPADES ] \n";
my $exec6_10_2 = "spades -s $path_unmapped_trimmed_filtered_20_23  --careful --only-assembler -t $process  -k $hash -o $step5_opt_20to23/run2 ";
print "\nSTEP6_10_2\n\t $exec6_2\n";
`$exec6_10_2`;

print "\t#Running step 6_13 [ merge assemblies - mergeContigs.pl ] \n";
my $exec6_13 = "$path_merge_contigs -contig1 $step5_opt_20to23/run1/contigs.fa -contig2 $step5_opt_20to23/run2/scaffolds.fasta -output $step5_opt_20to23/contigs.final.fasta ";
print "\nSTEP6_13\n\t $exec6_13\n";
`$exec6_13`;

   
#######################################################################
### Running velvet optmiser (hash 17) 24-30 ---------------------------
#######################################################################

if (defined($deg)) {
    
    print "\n[VELVET OPTIMISER - 24-30nt]\n";
    print "\t#Running step 6_10 [ Assemble unmapped 24-30nt nt - velvetOptimser.pl ]\n";
    
    my $exec62_10 = "velvetoptimiser --d $step5_opt_24to30/run1 --t $process --s 15 --e 17 --f '-short -fasta $path_unmapped_trimmed_filtered_24_30' --a 2>>$prefix.warn";
    print "\nSTEP62_10\n\t $exec62_10\n";
    `$exec62_10`;

    `mkdir $step5_opt_24to30/run2 `;

    print "\t#Running SPADES fixed hash  [ SPADES ] \n";
    my $exec6_10_3 = "spades -s $path_unmapped_trimmed_filtered_24_30 --careful --only-assembler -t $process  -k 15,17 -o $step5_opt_24to30/run2 ";
    print "\nSTEP6_10_3\n\t $exec6_2\n";
    `$exec6_10_3`;

    print "\t#Running step 6_13 [ merge assemblies - mergeContigs.pl ] \n";
    my $exec62_13 = "$path_merge_contigs -contig1 $step5_opt_24to30/run1/contigs.fa -contig2 $step5_opt_24to30/run2/scaffolds.fasta -output $step5_opt_24to30/contigs.final.fasta ";
    print "\nSTEP62_13\n\t $exec62_13\n";
    `$exec62_13`;
}

#######################################################################
### Merging assemblies ------------------------------------------------
#######################################################################

if (not defined($deg)) {
    `touch $step5_opt_24to30/contigs.final.fasta`;
}

print "\n[MERGING CONTIGS AND RUNNING CAP3]\n";
print "\t#Running STEP [CAT] Concatenating contigs...\n";
my $exec_cat = "cat $step5_opt_20to23/contigs.final.fasta $step5_opt_fix/contigs.final.fasta  $step5_fix/contigs.final.fasta $step5_opt/contigs.final.fasta $step5_opt_24to30/contigs.final.fasta > $step5_cap3/all_contigs.fasta";
print "\nSTEP6_CAT\n\t $exec_cat\n";
`$exec_cat`;

print "\t#Running step CAP3 [ assembly contigs from different assemblies ] \n";
my $exec_cap3 = "cap3 $step5_cap3/all_contigs.fasta  > $step5_cap3/log.cap3";
print "\nSTEP CAP3\n\t $exec_cap3\n";
`$exec_cap3`;

print "\t#Running [ Concatenning contigs and singlets from CAP3]\n";
my $exec_cap3_1 = "cat $step5_cap3/all_contigs.fasta.cap.contigs $step5_cap3/all_contigs.fasta.cap.singlets > $step5_cap3/contigs_merged.final.fasta";
print "\nSTEP6_CAT\n\t $exec_cap3_1\n";
`$exec_cap3_1`;

print "\t#Running step 6_7 [ selecting contigs larger than n50 - calcN50.pl ] \n";
my $step6_7 = "$path_fix_cap3_contigs -i $step5_cap3/contigs_merged.final.fasta -s 50  -p $step5_cap3/contigs_merged.final";
print "\nSTEP6_7\n\t $step6_7\n";
`$step6_7`;

my $countAssembledContigs = `grep -c '>' $step5_cap3/contigs_merged.final.gt50.fasta`;
chomp($countAssembledContigs);
# print metrics "#total assembled contigs\t" . $countAssembledContigs. "\n"; # Assembled Contigs
print "#total assembled contigs\t" . $countAssembledContigs. "\n"; # Assembled Contigs
# print interest "#total assembled contigs\t" . $countAssembledContigs . "\n"; # Assembled Contigs
print "#total assembled contigs\t" . $countAssembledContigs . "\n"; # Assembled Contigs

print "[FILTER CONTIGS gt 200 nt]\n";
my $exec_FC2 = "python3 $path_filter_fasta_by_size $step5_cap3/contigs_merged.final.gt50.fasta 200 1000000 $step5_cap3/contigs_merged.final.gt200.fasta";
print "\nSTEP5_2\n\t $exec_FC2\n";
`$exec_FC2`;

$countAssembledContigs = `grep -c '>' $step5_cap3/contigs_merged.final.gt200.fasta`;
chomp($countAssembledContigs);
# print metrics "#contigs gt200\t".$countAssembledContigs. "\n"; # Assembled Contigs
print "#contigs gt200\t".$countAssembledContigs. "\n"; # Assembled Contigs
# print interest "#contigs gt200\t".$countAssembledContigs . "\n"; # Assembled Contigs
print "#contigs gt200\t".$countAssembledContigs . "\n"; # Assembled Contigs

#######################################################################
### Blastn ------------------------------------------------------------
#######################################################################

print "\n[BlastN contigs gt 200]\n";
print "\t#Running step 10_111 [ Blast against NT - blast+ blastn ]\n";
my $exec10_111 = "blastn -query $step5_cap3/contigs_merged.final.gt200.fasta -db $path_blast_db_nt/nt -num_descriptions 5 -num_alignments 5 -evalue 1e-5 -out $step6/contigs_merged.final.gt200.1e5.blastn  -num_threads $process";

print "\nSTEP10_111\n\t $exec10_111\n";
`$exec10_111`;

print "\t#Running filterblast...\n";
my $exec10_112 = "$path_filter_blast -b $step6/contigs_merged.final.gt200.1e5.blastn -evalue 1e-5 --best --desc > $step7/contigs.blastn.1e5.report";

print "\nSTEP10_112\n\t $exec10_112\n";
`$exec10_112`;

print "\n[Extracting contigs all Hits blastn 1e-5]\n";
`$path_analyse_contigs_blastn -i $step7/contigs.blastn.1e5.report  -f $step5_cap3/contigs_merged.final.gt200.fasta -q "" -p  $step7/contigs.bN.analyze --fasta > $step7/contigs.blastN.analyze`;

print "\n[Extracting contigs viral blastn 1e-5]\n";
`$path_analyse_contigs_blastn -i $step7/contigs.blastn.1e5.report  -f $step5_cap3/contigs_merged.final.gt200.fasta -q "virus" -p  $step7/contigs.bN.blastn.analyze --fasta > $step7/contigs.blastN.virus.analyze`;

print "\n[Extracting contigs nonviral blastn 1e-5]\n";
`grep -v -i "virus" $step7/contigs.bN.analyze..contigs.fasta | grep '>' | cut -f1 -d " " >  $step7/aux_nonviral`;

`fasta_formatter -i $step7/contigs.bN.analyze..contigs.fasta > $step7/aux.fasta`;
`while read p; do grep -A1 \${p} $step7/aux.fasta >> $step7/contigs.bN.analyze.nonviral.contigs.fasta;done < $step7/aux_nonviral`;
`rm -f $step7/aux.fasta`;
`rm -f $step7/aux_nonviral`;

my $hitsBlastn = `grep -c '>' $step7/contigs.bN.analyze..contigs.fasta`;
chomp($hitsBlastn);

# 
# TODO: 2023-05-23 - Check what to do with these 'parallel' logging files
# 
# print metrics "#contigs hit blastN\t".$hitsBlastn."\n";
print "#contigs hit blastN\t".$hitsBlastn."\n";
# Assembled Contigs
# print interest "#contigs hit blastN\t".$hitsBlastn."\n";
print "#contigs hit blastN\t".$hitsBlastn."\n";
# Assembled Contigs

my $hitsVirusBlastn = `grep -c '>' $step7/contigs.bN.blastn.analyze.virus.contigs.fasta`;
chomp($hitsVirusBlastn);
# print metrics "#contigs hit VIRUS blastN\t".$hitsVirusBlastn."\n";
print "#contigs hit VIRUS blastN\t".$hitsVirusBlastn."\n";
# Assembled Contigs
# print interest "#contigs hit VIRUS blastN\t".$hitsVirusBlastn."\n";
print "#contigs hit VIRUS blastN\t".$hitsVirusBlastn."\n";
# Assembled Contigs

print "\n[Extracting contigs no hit blastn 1e-5]\n";
my $exec10_113 = "$path_extract_seqs_no_hit_blast -seq $step5_cap3/contigs_merged.final.gt200.fasta -blast $step7/contigs.blastn.1e5.report -out $step6/seqNoHit.blastN.1e5.fasta";
print "\nSTEP10_113\n\t $exec10_113\n";
`$exec10_113`;

my $seqsNoHitBlastn = `grep -c '>' $step6/seqNoHit.blastN.1e5.fasta`;
chomp($seqsNoHitBlastn);

# 
# TODO: 2023-05-23 - Check what to do with these 'parallel' logging files
# 

# print metrics "#contigs not hit blastN\t".$seqsNoHitBlastn."\n";
print "#contigs not hit blastN\t".$seqsNoHitBlastn."\n";
# Assembled Contigs

`cat $step7/contigs.bN.blastn.analyze.virus.contigs.fasta | perl -pi -e 's/>(\\S+) (\\S+) (\\S+) (\\S+).+/>blastN_\$1_\$2_\$3_\$4/g' > $step7/contigs.virus.blastN.formatted.fasta`;

#######################################################################
### DIAMOND (Blastx) --------------------------------------------------
#######################################################################

print "\n[Diamond (BlastX) contigs gt 200]\n";
print "\t#Running step 9 [ Diamond-Blast against NR ]\n";

my $exec9 = "diamond blastx -q $step6/seqNoHit.blastN.1e5.fasta -d $path_diamond_nr -k 5 -p $process -e 0.001 -f 0 -c 1 -b 20 --very-sensitive -o $step7/diamond_blastx.out --un $step7/diamond_blastx_NoHits.fasta --unfmt fasta --al $step7/diamond_blastx_Hits.fasta --alfmt fasta 2>  $step7/diamond.log ";

print "\nSTEP9\n\t $exec9\n";
`$exec9`;

print "\t#Filtering Diamond results... ]\n";   
my $exec9_11 = "$path_filter_diamond -f $step7/diamond_blastx_Hits.fasta -o $step7/diamond_blastx.out -d $step7";

print "\nSTEP9_11\n\t $exec9_11\n";
`$exec9_11`;

#######################################################################

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
