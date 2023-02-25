#!/usr/bin/perl

#### ericgdp & jpbio #### 
###############################################################################################################
#v2.3
#Modificacoes jp:
#1-blastn+ novo implementado
#2-Apenas reads >= 200 sao considerados para alinhamentos e analises de smal RNAs
#3-Diamond implementado
#4-ranges de montagem sao 20-23; 24-30 e recomendo o usuario passar 15-35 na linha de comando...
#5-spades esta atualizado e roda em modo normal (sem --meta) com a opcao --careful adicionada 
#6- nao tem mais opcao de colorspace
#7- adaptacao do -nohostfilter, os reads sempre serao filtrados contra bacterias
#################################################################################################################

use Getopt::Long;

my $usage = "

$0 -fasta <reads in letter space without quality> -fastq <reads in fastq format with quality>  -qual <quality file> -fastqgz <reads in fastq format compressed gzip> [-nohostfilter]  -prefix <prefix>   -hostgenome <reference Host genome for reads filtering> -process <number of process to use> -size <expected size>  -log <log file name>
$0 -h

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

$| = 1;     # forces immediate prints into files rather than the buffer.

my $fastqgz;
my $qual;
my $hostgenome;
my $process;
my $size;
my $prefix;
my $log;
my $fasta;
my $fastq;
my $hash;
my $se;
my $si;
my $clean;
my $nohostfilter;
my $large_index;

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
    print "\t $exec1\n";
    die $usage;
}
if (not defined($nohostfilter) and not defined($fasta) and not defined($fastqgz)) {
    print "aqui0";
    if (not(defined $fastq)) {
        if ((not(defined($fastqgz))) and (not defined($fastq)) and (not defined($fasta))) {
            die "\nGive an input file valid!\n ", $usage;
        }

        if (not(defined($qual)) and (not defined($fastq)) and (not defined($fasta))) {
            die "\nGive an input quality file name valid! \n", $usage;
        }
    }
}

if (defined($large_index)) {
    $large_index = " --large-index ";
} else {
    $large_index = " ";
}

if (not(defined($hostgenome))) {
    die "\nGive an input reference file valid! \n", $usage;
}

if (not(defined($process))) {
    die "\nGive an number of process to use! \n", $usage;
}

if (not(defined($si)) or not defined($se)) {
    die "\nGive a range size of reads to use! ex: -si 21 -se 23 \n", $usage;
}

if (not(defined($size))) {
    die "\nGive the expected size of the genome! \n", $usage;
}

if (not(defined($prefix))) {
    die "\nGivea prefix name! \n", $usage;
}

if (not(defined($hash))) {
    $hash = 15;
}


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

print"#########################################################\n\n\n";

my $step0		="$prefix/$prefix"."_0_saet";
my $step1		="$prefix/$prefix"."_1_trimming";
my $step2		="$prefix/$prefix"."_2_filter_size_gaps_convertion";
my $step3		="$prefix/$prefix"."_3_mapping_vector";
my $step4		="$prefix/$prefix"."_4_getUnmapped";
my $step5_fix	="$prefix/$prefix"."_5_assembleUnmapped_fix";
my $step5_opt	="$prefix/$prefix"."_5_assembleUnmapped_opt";
my $step5_opt_fix="$prefix/$prefix"."_5_assembleUnmapped_opt_fix";
my $step5_opt_20to23="$prefix/$prefix"."_5_assembleUnmapped_opt_20to23";
my $step5_opt_24to30="$prefix/$prefix"."_5_assembleUnmapped_opt_24to30";
my $step5_contigs= "$prefix/$prefix"."_5_assembleUnmapped_final";
my $step5_cap3  = "$prefix/$prefix"."_5_cap3";
my $step6		="$prefix/$prefix"."_6_blast";
my $step7		="$prefix/$prefix"."_7_reportBlast";
my $step8		="$prefix/$prefix"."_8_completeReport";
my $step9		="$prefix/$prefix"."_9_contigs_no_hit";
my $step10		="$prefix/$prefix"."_10_pattern";
my $step_virus	="$prefix/$prefix"."_virus";


##########################
   #Creating step folders	     
##########################

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
if (not -e $step5) {
    `mkdir $step5`;
}
if (not -e $step5) {
    `mkdir $step5_fix`;
    `mkdir $step5_opt`;
    `mkdir $step5_opt_fix`;
    `mkdir $step5_contigs`;
    `mkdir $step5_cap3`;
    `mkdir $step5_opt_20to23`;
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

open(metrics, ">$step8/full_metrics.txt");
open(interest, ">$step8/metrics_of_interest.txt");
open(LOG, ">$log");

our $binary = "/home/bioinfo/eric_bin";

###########################
   #handle FASTQ sequences       
###########################

if ( defined($fastq)) {

	print "Processing fastq sequences...\n";
	title("Total Reads on fastq file");
	txt(`countFastq $fastq`);
	print "\n\nRunning step 0 [ quality filter - fastq_quality_filter ]\n";
	my $exec_fq_0="fastq_quality_filter -Q 33 -q 15 -p 60 -i $fastq -o $step0/trimming.quality.fastq ";
	print LOG "\nSTEP0\n\t $exec_fq_0\n";
	`$exec_fq_0`;


	#trimming adapters
	print "\n\nRunning step 1 [ trimming adapter - fastq_clipper ]\n";
	my $exec_fq_1="fastx_clipper -Q 33 -l 17 -c -i $step0/trimming.quality.fastq -a $adapter -o $step2/trimmed.quality.gt15.fastq";
	print LOG "\nSTEP1\n\t $exec_fq_1\n";
	`$exec_fq_1`;

	print "\n\nRunning step 2 [ converting fastq to fasta - fastq_to_fasta ]\n";
	#converting fastq to fasta
	my $exec_fq_2="fastq_to_fasta -Q 33 -i $step2/trimmed.quality.gt15.fastq -o $step2/trimmed_filtered_gt15.fasta";
	print LOG "\nSTEP2\n\t $exec_fq_2\n";
	`$exec_fq_2`;
	title("Number of trimmed reads\n");
	title(`grep -c ">" $step2/trimmed_filtered_gt15.fasta`);
	
	
}

###########################
   #handle FASTA sequences        
###########################

elsif (defined($fasta)) {
     print "#Loading FASTA file ... \n";
     if (not defined ($nohostfilter)) {
         print "[COPING $fasta TO $step2/trimmed_filtered_gt15.fasta]\n";
        `cp $fasta $step2/trimmed_filtered_gt15.fasta `;
          print LOG "\nSTEP3\n\t cp $fasta $step2/trimmed_filtered_gt15.fasta \n";
     } else {         
            #Mapping Host-unfiltered reads against bacters reference    
            print "[MAPPING HOST-UNFILTERED READS AGAINST BACTERIAL GENOMES]... \n";
            my $exec5_1="bowtie -f -S -v 1 --un $step4/unmappedVectorBacters.fasta -k 1 -p $process --large-index /media/data/reference/bacterial_genomes/all_bacters.fasta $fasta > /dev/null 2>> $step4/reads_mapped_to_bacteria.log ";
            print LOG "\nSTEP5_1\n\t $exec5_1\n";
           `$exec5_1`;
           
            $nReadsUnmapHostBac = `grep -c '>' $step4/unmappedVectorBacters.fasta`;
            chomp($nReadsUnmapHostBac);
            $mappedbac = $nReadsUnmapHost - $nReadsUnmapHostBac;
			
			print metrics "#reads mapped bacter\t".$mappedbac  ."\n"; # mapped reads Bacterial genomes
			print metrics "#preprocessed reads\t".$nReadsUnmapHostBac  ."\n"; # pre-processed reads
			
			print "\n  PRE-PROCESSING FINISHED \n";
    }
    
}
    if (not defined($nohostfilter)) {
		if (defined($fastqgz)) { ### dealing with FASTQ TRIMMED files compacted with GZIP
			           
			print "\n\nRunning step 2 [ converting fastq to fasta - fastq_to_fasta ]\n";   
	        
			#converting fastq to fasta
			my $exec_fq_2="gunzip -dc $fastqgz | fastq_to_fasta -Q 33 -o $step2/trimmed_filtered_gt15.fasta";
			print LOG "\nSTEP2\n\t $exec_fq_2\n";
			`$exec_fq_2`;
         } 

	
            print "[MAPPING SEQUENCE AGAINST VECTOR]\n";
            my $exec3="bowtie $large_index -f -S -k 1 -p $process -v 1 --un $step4/unmappedVectorReads.fasta  $hostgenome $step2/trimmed_filtered_gt15.fasta | awk -F'\\t' '{if( \$2 != 4) print \$0}' > $step3/mapped_host.v1.sam  2>mapping_host.stats  ";
            print  "\nSTEP3\n\t $exec3\n";
            `$exec3  `;
                      
           ##count total reads
            $nReads = `grep -c '>' $step2/trimmed_filtered_gt15.fasta`;
            chomp($nReads);
            print metrics "#total reads\t".$nReads ."\n";
            print interest "#total reads\t".$nReads ."\n";

           
            print "[PLOT GERAL DISTRIBUTION READS MAPPED HOST ]\n";
            my $exec3_1="$binary/plotGeralDistributionPerBaseByReads.pl -sam $step3/mapped_host.v1.sam  -s $si -e $se -p $step3/$prefix.mapped_host.v1.$si.$se -norm $nReads --plot";
			print  "\nSTEP3\n\t $exec3_1\n";
            print LOG "\nSTEP3\n\t $exec3_1\n";
            `$exec3_1`; 

			    print "[PLOT GERAL DISTRIBUTION READS MAPPED HOST - 15-35nt ]\n";
            my $exec3_11="$binary/plotGeralDistributionPerBaseByReads.pl -sam $step3/mapped_host.v1.sam  -s 15 -e 35 -p $step3/$prefix.mapped_host.v1 -norm $nReads --plot";
			print  "\nSTEP3\n\t $exec3_11\n";
            print LOG "\nSTEP3\n\t $exec3_11\n";
            `$exec3_11`; 


         
			 $nReadsUnmapHost = `grep -c '>' $step4/unmappedVectorReads.fasta`;
            chomp($nReadsUnmapHost);
            $mapped = $nReads - $nReadsUnmapHost;
            
            print metrics "#reads mapped host\t".$mapped  ."\n"; # mapped reads HOST
			print metrics "#reads unmapped host\t".$nReadsUnmapHost  ."\n"; # unmapped reads HOST
			
			print interest "#reads mapped host\t".$mapped  ."\n"; # mapped reads HOST
			print interest "#reads unmapped host\t".$nReadsUnmapHost  ."\n"; # unmapped reads HOST
			
			`rm -rf $step3/mapped_host.v1.sam`; # deleting sam file mapped reads on host genome

          
            #Mapping Host-filtered reads against bacters reference
            print "[MAPPING HOST-FILTERED READS AGAINST BACTERIAL GENOMES]... \n";
            my $exec5_1="bowtie -f -S -v 1 --un $step4/unmappedVectorBacters.fasta -k 1 -p $process --large-index /media/data/reference/bacterial_genomes/all_bacters.fasta $step4/unmappedVectorReads.fasta > /dev/null 2>>$prefix.warn ";
            print LOG "\nSTEP5_1\n\t $exec5_1\n";
           `$exec5_1`;
           
            $nReadsUnmapHostBac = `grep -c '>' $step4/unmappedVectorBacters.fasta`;
            chomp($nReadsUnmapHostBac);
            $mappedbac = $nReadsUnmapHost - $nReadsUnmapHostBac;
			
			print metrics "#reads mapped bacter\t".$mappedbac  ."\n"; # mapped reads Bacterial genomes
			print metrics "#preprocessed reads\t".$nReadsUnmapHostBac  ."\n"; # pre-processed reads
			


			print "\n  PRE-PROCESSING FINISHED \n";
    }
	
#########################################
   #selecting filtered sequences by size    
#########################################
        #como default temos passado 18 - 30 nt na linha de comando na hora de executar (mas eu vou passar 15-35)
        print "[FILTER UNMAPPED SEQUENCES BY SIZE (variable size $si to $se)]\n";
        my $exec5= "python /home/bioinfo/eric_bin/filter_fasta_by_size.py $step4/unmappedVectorBacters.fasta $si $se $step4/unmapped_trimmed_filtered.fasta -t F ";
        print LOG "\nSTEP5\n\t $exec5\n";
        `$exec5`;
        
        print "[FILTER UNMAPPED SEQUENCES BY SIZE (20-23NT)]\n";
        my $exec5_1= "python /home/bioinfo/eric_bin/filter_fasta_by_size.py $step4/unmappedVectorBacters.fasta 20 23 $step4/unmapped_trimmed_filtered.20-23.fasta";
        print LOG "\nSTEP5_1\n\t $exec5_1\n";
        `$exec5_1`;
               
        
        print "[FILTER UNMAPPED SEQUENCES BY SIZE (24-30NT)]\n";
        my $exec5_2= "python /home/bioinfo/eric_bin/filter_fasta_by_size.py $step4/unmappedVectorBacters.fasta 24 30 $step4/unmapped_trimmed_filtered.24-30.fasta";
        print LOG "\nSTEP5_2\n\t $exec5_2\n";
        `$exec5_2`;
        
           
	
#########################################################
   #Running velvet optmiser (automatically defined hash)
#########################################################
	    print "\n#[RUNNING VELVET OPTIMIZER]\n";
        print "\t#Running step 6 [ Assemble unmapped 21 nt - velvetOptimser.pl ]\n";
        my $exec6_1="/home/bioinfo/source_programs_eric/VelvetOptimiser-2.2.5/VelvetOptimiser.pl --d $step5_opt/run1 --t $process --s 13 --e 19 --f '-short -fasta $step4/unmapped_trimmed_filtered.fasta' --a $process 2>>$prefix.warn";
        print LOG "\nSTEP6_1\n\t $exec6_1\n";
        `$exec6_1`;

      
        print "\t#Running step 6_4 [ SPADES ] \n";
		my $exec6_4="/media/group1_data/joaopaulo/programs_jp/SPAdes-3.15.4-Linux/bin/./spades.py -s $step4/unmapped_trimmed_filtered.fasta --careful --only-assembler -t $process  -k 13,15,17,19 -o $step5_opt/run2";
        print LOG "\nSTEP6_4\n\t $exec6_4\n";
        `$exec6_4`;

        print "\t#Running step 6_5 [ merge assemblies - mergeContigs.pl ] \n";
        my $exec6_5="$binary/mergeContigsNew.pl -contig1 $step5_opt/run1/contigs.fa -contig2 $step5_opt/run2/scaffolds.fasta -output $step5_opt/contigs.final.fasta ";
        print LOG "\nSTEP6_5\n\t $exec6_5\n";
        `$exec6_5`;
        

################################
   #Running velvet (fixed hash)                  
################################
        print "\n[RUNNING DEFAULT VELVET]\n";
        print "\t#Running step 6 [ Assemble unmapped 21 nt - velvet hash $hash ]\n";
        my $exec6 = "velveth $step5_fix/run1 $hash -fasta -short $step4/unmapped_trimmed_filtered.fasta 2>>$prefix.warn";
        print LOG "\nSTEP6\n\t $exec6\n";
        `$exec6`;

        `velvetg $step5_fix/run1 -exp_cov auto -cov_cutoff auto 2>$prefix.warn`;
		
		`mkdir $step5_fix/run2`;
		
        print "\t#Running SPADES fixed hash  [ SPADES ] \n";
		my $exec6_2="/media/group1_data/joaopaulo/programs_jp/SPAdes-3.15.4-Linux/bin/./spades.py -s $step4/unmapped_trimmed_filtered.fasta --careful --only-assembler -t $process  -k $hash -o $step5_fix/run2 ";
        print LOG "\nSTEP6_2\n\t $exec6_2\n";
        `$exec6_2`;

        print "\t#Running step 6_5 [ merge assemblies - mergeContigs.pl ] \n";
        my $exec6_5="$binary/mergeContigsNew.pl -contig1 $step5_fix/run1/contigs.fa -contig2 $step5_fix/run2/scaffolds.fasta -output $step5_fix/contigs.final.fasta ";
        print LOG "\nSTEP6_5\n\t $exec6_5\n";
        `$exec6_5`;
        
		
#########################################
   #Running velvet optmiser (FIXED hash)         
#########################################
       
        print "\n[VELVET OPTIMISER HASH ONLY 15]\n";
        print "\t#Running step 6_6 [ Assemble unmapped 21 nt - velvetOptimser.pl ]\n";
        my $exec6_6="/home/bioinfo/source_programs_eric/VelvetOptimiser-2.2.5/VelvetOptimiser.pl --d $step5_opt_fix/run1 --t $process --s $hash --e $hash --f '-short -fasta $step4/unmapped_trimmed_filtered.fasta' --a 2>>$prefix.warn";
        print LOG "\nSTEP6_1\n\t $exec6_6\n";
        `$exec6_6`;
        
        
		`mkdir $step5_opt_fix/run2`;
		print "\t#Running SPADES fixed hash  [ SPADES ] \n";
		my $exec6_6_2="/media/group1_data/joaopaulo/programs_jp/SPAdes-3.15.4-Linux/bin/./spades.py -s $step4/unmapped_trimmed_filtered.fasta --careful --only-assembler -t $process  -k 15 -o $step5_opt_fix/run2 ";
        print LOG "\nSTEP6_2\n\t $exec6_2\n";
        `$exec6_6_2`;

        print "\t#Running step 6_9 [ merge assemblies - mergeContigs.pl ] \n";
        my $exec6_9="$binary/mergeContigsNew.pl -contig1 $step5_opt_fix/run1/contigs.fa -contig2 $step5_opt_fix/run2/scaffolds.fasta -output $step5_opt_fix/contigs.final.fasta ";
        print LOG "\nSTEP6_5\n\t $exec6_5\n";
        `$exec6_9`;

###############################################
   #Running velvet optmiser (FIXED hash) 20-23  
###############################################
       
        print "\n[VELVET OPTIMISER HASH ONLY 15 - 20-23nt]\n";
        print "\t#Running step 6_10 [ Assemble unmapped 20-23nt nt - velvetOptimser.pl ]\n";
        my $exec6_10="/home/bioinfo/source_programs_eric/VelvetOptimiser-2.2.5/VelvetOptimiser.pl --d $step5_opt_20to23/run1 --t $process --s $hash --e $hash --f '-short -fasta $step4/unmapped_trimmed_filtered.20-23.fasta' --a 2>>$prefix.warn";
        print LOG "\nSTEP6_1\n\t $exec6_10\n";
        `$exec6_10`;

        `mkdir $step5_opt_20to23/run2`;
		print "\t#Running SPADES fixed hash  [ SPADES ] \n";
		my $exec6_10_2="/media/group1_data/joaopaulo/programs_jp/SPAdes-3.15.4-Linux/bin/./spades.py -s $step4/unmapped_trimmed_filtered.20-23.fasta  --careful --only-assembler -t $process  -k $hash -o $step5_opt_20to23/run2 ";
        print LOG "\nSTEP6_10_2\n\t $exec6_2\n";
        `$exec6_10_2`;

        print "\t#Running step 6_13 [ merge assemblies - mergeContigs.pl ] \n";
        my $exec6_13="$binary/mergeContigsNew.pl -contig1 $step5_opt_20to23/run1/contigs.fa -contig2 $step5_opt_20to23/run2/scaffolds.fasta -output $step5_opt_20to23/contigs.final.fasta ";
        print LOG "\nSTEP6_13\n\t $exec6_13\n";
        `$exec6_13`;

    if (defined($deg)) {
   
############################################
   #Running velvet optmiser (hash 17) 24-30  
############################################
	    print "\n[VELVET OPTIMISER - 24-30nt]\n";
        print "\t#Running step 6_10 [ Assemble unmapped 24-30nt nt - velvetOptimser.pl ]\n";
        my $exec62_10="/home/bioinfo/source_programs_eric/VelvetOptimiser-2.2.5/VelvetOptimiser.pl --d $step5_opt_24to30/run1 --t $process --s 15 --e 17 --f '-short -fasta $step4/unmapped_trimmed_filtered.24-30.fasta' --a 2>>$prefix.warn";
        print LOG "\nSTEP62_10\n\t $exec62_10\n";
        `$exec62_10`;
        
		`mkdir $step5_opt_24to30/run2 `;
		
       print "\t#Running SPADES fixed hash  [ SPADES ] \n";
		my $exec6_10_3="/media/group1_data/joaopaulo/programs_jp/SPAdes-3.15.4-Linux/bin/./spades.py -s $step4/unmapped_trimmed_filtered.24-30.fasta --careful --only-assembler -t $process  -k 15,17 -o $step5_opt_24to30/run2 ";
        print LOG "\nSTEP6_10_3\n\t $exec6_2\n";
        `$exec6_10_3`;


        print "\t#Running step 6_13 [ merge assemblies - mergeContigs.pl ] \n";
        my $exec62_13="$binary/mergeContigsNew.pl -contig1 $step5_opt_24to30/run1/contigs.fa -contig2 $step5_opt_24to30/run2/scaffolds.fasta -output $step5_opt_24to30/contigs.final.fasta ";
        print LOG "\nSTEP62_13\n\t $exec62_13\n";
        `$exec62_13`;
    }
	
	
	
#######################
   #Merging assemblies           
#######################
	    if (not defined($deg)) {
	        `touch  $step5_opt_24to30/contigs.final.fasta`;
	       }
	    
        print "\n[MERGING CONTIGS AND RUNNING CAP3]\n";
        print "\t#Running STEP [CAT] Concatenating contigs...\n";
        my $exec_cat="cat $step5_opt_20to23/contigs.final.fasta $step5_opt_fix/contigs.final.fasta  $step5_fix/contigs.final.fasta $step5_opt/contigs.final.fasta $step5_opt_24to30/contigs.final.fasta > $step5_cap3/all_contigs.fasta";
        print LOG "\nSTEP6_CAT\n\t $exec6_cat\n";
        `$exec_cat`;
 
        print "\t#Running step CAP3 [ assembly contigs from different assemblies ] \n";
        my $exec_cap3="$binary/cap3 $step5_cap3/all_contigs.fasta  > $step5_cap3/log.cap3";
        print LOG "\nSTEP CAP3\n\t $exec_cap3\n";
        `$exec_cap3`;
    
        print "\t#Running [ Concatenning contigs and singlets from CAP3]\n";
        my $exec_cap3_1="cat $step5_cap3/all_contigs.fasta.cap.contigs $step5_cap3/all_contigs.fasta.cap.singlets > $step5_cap3/contigs_merged.final.fasta";
        print LOG "\nSTEP6_CAT\n\t $exec_cap3_1\n";
        `$exec_cap3_1`;
	
  	    print "\t#Running step 6_7 [ selecting contigs larger than n50 - calcN50.pl ] \n";
	    my $step6_7="$binary/fixIDcap3Contigs.pl -i $step5_cap3/contigs_merged.final.fasta -s 50  -p $step5_cap3/contigs_merged.final";
        print LOG "\nSTEP6_7\n\t $step6_7\n";
        `$step6_7`;
      
      	$countAssembledContigs = `grep -c '>' $step5_cap3/contigs_merged.final.gt50.fasta`;
        chomp($countAssembledContigs);
        print metrics "#total assembled contigs\t".$countAssembledContigs  ."\n"; # Assembled Contigs
        print interest "#total assembled contigs\t".$countAssembledContigs  ."\n"; # Assembled Contigs
            
      
        print "[FILTER CONTIGS gt 200 nt]\n";
        my $exec_FC2= "python /home/bioinfo/eric_bin/filter_fasta_by_size.py $step5_cap3/contigs_merged.final.gt50.fasta 200 1000000 $step5_cap3/contigs_merged.final.gt200.fasta";
        print LOG "\nSTEP5_2\n\t $exec_FC2\n";
        `$exec_FC2`;
      
		$countAssembledContigs = `grep -c '>' $step5_cap3/contigs_merged.final.gt200.fasta`;
        chomp($countAssembledContigs);
        print metrics "#contigs gt200\t".$countAssembledContigs  ."\n"; # Assembled Contigs
        print interest "#contigs gt200\t".$countAssembledContigs  ."\n"; # Assembled Contigs

###########
   #Blastn
###########
	  print "\n[BlastN contigs gt 200]\n";
        print "\t#Running step 10_111 [ Blast against NT - blast+ blastn ]\n";
        my $exec10_111 = "/media/group3_data/ncbi-blast-2.13.0+/bin/blastn -query $step5_cap3/contigs_merged.final.gt200.fasta -db /media/group3_data/blastdb_2022/nt/nt -num_descriptions 5 -num_alignments 5 -evalue 1e-5 -out $step6/contigs_merged.final.gt200.1e5.blastn  -num_threads $process";
        print LOG "\nSTEP10_111\n\t $exec10_111\n";
        `$exec10_111`;
		
		print "\t#Running filterblast...\n";
        my $exec10_112 = "$binary/filterblast.pl -b $step6/contigs_merged.final.gt200.1e5.blastn -evalue 1e-5 --best --desc > $step7/contigs.blastn.1e5.report";
        print LOG "\nSTEP10_112\n\t $exec10_112\n";
        `$exec10_112`;
        print "\n[Extracting contigs all Hits blastn 1e-5]\n";   
        `$binary/analyzeContigsFilterBlastn.pl -i $step7/contigs.blastn.1e5.report  -f $step5_cap3/contigs_merged.final.gt200.fasta -q "" -p  $step7/contigs.bN.analyze --fasta > $step7/contigs.blastN.analyze`;
        print "\n[Extracting contigs viral blastn 1e-5]\n";
        `$binary/analyzeContigsFilterBlastn.pl -i $step7/contigs.blastn.1e5.report  -f $step5_cap3/contigs_merged.final.gt200.fasta -q "virus" -p  $step7/contigs.bN.blastn.analyze --fasta > $step7/contigs.blastN.virus.analyze`;
       
        print "\n[Extracting contigs nonviral blastn 1e-5]\n";
        `grep -v -i "virus" $step7/contigs.bN.analyze..contigs.fasta | grep '>' | cut -f1 -d " " >  $step7/aux_nonviral`;
        `fasta_formatter -i $step7/contigs.bN.analyze..contigs.fasta > $step7/aux.fasta`;
        `while read p; do grep -A1 \${p} $step7/aux.fasta >> $step7/contigs.bN.analyze.nonviral.contigs.fasta;done < $step7/aux_nonviral`;
        `rm -f $step7/aux.fasta`;
        `rm -f $step7/aux_nonviral`;        
        
        $hitsBlastn = `grep -c '>' $step7/contigs.bN.analyze..contigs.fasta`;
        chomp($hitsBlastn);
        print metrics "#contigs hit blastN\t".$hitsBlastn  ."\n"; # Assembled Contigs
        print interest "#contigs hit blastN\t".$hitsBlastn  ."\n"; # Assembled Contigs        
        
        $hitsVirusBlastn = `grep -c '>' $step7/contigs.bN.blastn.analyze.virus.contigs.fasta`;
        chomp($hitsVirusBlastn);
        print metrics "#contigs hit VIRUS blastN\t".$hitsVirusBlastn  ."\n"; # Assembled Contigs
        print interest "#contigs hit VIRUS blastN\t".$hitsVirusBlastn  ."\n"; # Assembled Contigs        
        
		print "\n[Extracting contigs no hit blastn 1e-5]\n";
        my $exec10_113="$binary/extractSequencesNoHitBlast.pl -seq $step5_cap3/contigs_merged.final.gt200.fasta -blast $step7/contigs.blastn.1e5.report -out $step6/seqNoHit.blastN.1e5.fasta";
        print LOG "\nSTEP10_113\n\t $exec10_113\n";
        `$exec10_113`;
     
        $seqsNoHitBlastn = `grep -c '>' $step6/seqNoHit.blastN.1e5.fasta`;
        chomp($seqsNoHitBlastn);
        print metrics "#contigs not hit blastN\t".$seqsNoHitBlastn  ."\n"; # Assembled Contigs
     
     	`cat $step7/contigs.bN.blastn.analyze.virus.contigs.fasta | perl -pi -e 's/>(\\S+) (\\S+) (\\S+) (\\S+).+/>blastN_\$1_\$2_\$3_\$4/g' > $step7/contigs.virus.blastN.formatted.fasta`;
     
       
#####################
   #DIAMOND (Blastx)
#####################
        
        print "\n[Diamond (BlastX) contigs gt 200]\n";
        print "\t#Running step 9 [ Diamond-Blast against NR ]\n";
        my $exec9 = "/media/group3_data/diamond blastx -q $step6/seqNoHit.blastN.1e5.fasta -d /media/group3_data/nr.dmnd -k 5 -p $process -e 0.001 -f 0 -c 1 -b 20 --very-sensitive -o $step7/diamond_blastx.out --un $step7/diamond_blastx_NoHits.fasta --unfmt fasta --al $step7/diamond_blastx_Hits.fasta --alfmt fasta 2>  $step7/diamond.log ";
        print LOG "\nSTEP9\n\t $exec9\n";
        `$exec9`;
        
        print "\t#Filtering Diamond results... ]\n";   
        my $exec9_11 = "$binary/./filter_diamond.sh -f $step7/diamond_blastx_Hits.fasta -o $step7/diamond_blastx.out -d $step7";
        print LOG "\nSTEP9_11\n\t $exec9_11\n";
        `$exec9_11`;

       
##############################################################         
   #Merge sequences viral hits(blastn and diamond) and nohits   
##############################################################         
         `cat $step7/contigs.virus.blastN.formatted.fasta $step7/diamond_blastx_Viral.fasta $step7/diamond_blastx_NoHits.fasta > $step9/seq_ViralHits_and_NoHits.fasta`;
         `bowtie-build $step9/seq_ViralHits_and_NoHits.fasta $step9/seq_ViralHits_and_NoHits.fasta`;
          
          print "\t#Mapping reads against viral hits(blastn and diamond) and nohits\n";
          my $exec10_13="bowtie -f -S -p $process -v 1 $step9/seq_ViralHits_and_NoHits.fasta $step4/unmappedVectorBacters.fasta 2>>$prefix.warn | awk -F'\t' '{if( \$2 != 4) print \$0}' > $step9/reads_vs_contigsHitNoHit.v1.sam ";
          print LOG "\nSTEP10_13\n\t $exec10_13\n";
          `$exec10_13`;
          
          `$binary/samStatistics_v3.pl -sam $step9/reads_vs_contigsHitNoHit.v1.sam -fa $step9/seq_ViralHits_and_NoHits.fasta -p $step9/seq_ViralHits_and_NoHits --profile --counts`;

          `$binary/Z-score.bothstrands.pl -sam $step9/reads_vs_contigsHitNoHit.v1.sam -p $step9/reads_vs_contigsHitNoHit `;
          
          `R --no-save $step9/reads_vs_contigsHitNoHit.zscore.tab $step9/reads_vs_contigsHitNoHit.zscore  < /home/bioinfo/eric_bin/R/heatmap_correlation_VISA.R 2>/dev/null`;
          
        ####stats
        $hitsBlastn = `grep '>' $step9/seq_ViralHits_and_NoHits.withSiRNA.fasta | grep blastN | wc -l `;
        chomp($hitsBlastn);
        print metrics "#contigs hit VIRUS blastN with siRNA\t".$hitsBlastn  ."\n"; 
        print interest "#contigs hit VIRUS blastN with siRNA\t".$hitsBlastn  ."\n";
        
        $hitsBlastn = `grep '>' $step9/seq_ViralHits_and_NoHits.withSiRNA_and_PiRNA.fasta | grep blastN | wc -l `;
        chomp($hitsBlastn);
        print metrics "#contigs hit VIRUS blastN with siRNA and piRNA\t".$hitsBlastn  ."\n"; 
        print interest "#contigs hit VIRUS blastN with siRNA and piRNA\t".$hitsBlastn  ."\n"; 
        
        $hitsBlastn = `grep '>' $step9/seq_ViralHits_and_NoHits.withPiRNA.fasta | grep blastN | wc -l `;
        chomp($hitsBlastn);
        print metrics "#contigs hit VIRUS blastN with piRNA\t".$hitsBlastn  ."\n";
        print interest "#contigs hit VIRUS blastN with piRNA\t".$hitsBlastn  ."\n";         
        
        $hitsVirusBlastn = `grep  '>' $step9/seq_ViralHits_and_NoHits.withSiRNA.fasta | grep blastX | wc -l `;
        chomp($hitsVirusBlastn);
        print metrics "#contigs hit VIRUS BlastX with siRNA \t".$hitsVirusBlastn  ."\n"; 
        print interest "#contigs hit VIRUS BlastX with siRNA \t".$hitsVirusBlastn  ."\n";         
        
        $hitsVirusBlastn = `grep  '>' $step9/seq_ViralHits_and_NoHits.withSiRNA_and_PiRNA.fasta| grep blastX | wc -l `;
        chomp($hitsVirusBlastn);
        print metrics "#contigs hit VIRUS BlastX with siRNA and piRNA \t".$hitsVirusBlastn  ."\n";
        print interest "#contigs hit VIRUS BlastX with siRNA and piRNA \t".$hitsVirusBlastn  ."\n";         
        
        $hitsVirusBlastn = `grep  '>' $step9/seq_ViralHits_and_NoHits.withPiRNA.fasta | grep blastX | wc -l `;
        chomp($hitsVirusBlastn);
        print metrics "#contigs hit VIRUS BlastX with piRNA \t".$hitsVirusBlastn  ."\n";
        print interest "#contigs hit VIRUS BlastX with piRNA \t".$hitsVirusBlastn  ."\n";         
        
        $hitsVirusBlastn = `grep  '>' $step9/seq_ViralHits_and_NoHits.withSiRNA.fasta | grep -v blastX | grep -v blastN | wc -l `;
        chomp($hitsVirusBlastn);
        print metrics "#contigs NOHIT blast with siRNA \t".$hitsVirusBlastn  ."\n";
        print interest "#contigs NOHIT blast with siRNA \t".$hitsVirusBlastn  ."\n";         
        
        $hitsVirusBlastn = `grep  '>' $step9/seq_ViralHits_and_NoHits.withSiRNA_and_PiRNA.fasta | grep -v blastN | -v grep blastX | wc -l `;
        chomp($hitsVirusBlastn);
        print metrics "#contigs NOHIT blast with siRNA and piRNA \t".$hitsVirusBlastn  ."\n"; 
        print interest "#contigs NOHIT blast with siRNA and piRNA \t".$hitsVirusBlastn  ."\n";        
        
        $hitsVirusBlastn = `grep  '>' $step9/seq_ViralHits_and_NoHits.withPiRNA.fasta | grep -v blastN | -v grep blastX | wc -l `;
        chomp($hitsVirusBlastn);
        print metrics "#contigs NOHIT blast with piRNA \t".$hitsVirusBlastn  ."\n";
        print interest "#contigs NOHIT blast with piRNA \t".$hitsVirusBlastn  ."\n";        
        
     close(metrics);
	
###########################
   #Pattern based analysis				
###########################

`mkdir $step_virus`;

####merge two files with contigs hit virus
`cat $step7/contigs.virus.blastN.formatted.fasta $step7/diamond_blastx_Viral.fasta > $step10/all_contigs_hit_virus.fasta`;
`bowtie-build $step10/all_contigs_hit_virus.fasta $step10/all_contigs_hit_virus.fasta`;
`bowtie -f -S -k 1 -p $process -v 1  $step10/all_contigs_hit_virus.fasta $step2/trimmed_filtered_gt15.fasta | awk -F '\t' '{if( \$2 != 4) print \$0}'  > $step10/reads.VS.contigs_hit_virus_blast.sam`;
`$binary/samStatistics_v3.pl -sam $step10/reads.VS.contigs_hit_virus_blast.sam -fa $step10/all_contigs_hit_virus.fasta -p $step10/reads.VS.contigs_hit_virus_blast --profile`;

#creating reference table with identified contigs hit virus by sequence similarity
`grep ">" $step10/all_contigs_hit_virus.fasta | cut -f 2 -d '>'  > $step10/all_contigs_hit_virus_sequence_similarity.tab`;

print "\n\n Calculating pattern viral contigs and candidates - HEATMAP \n\n";

###bowtie 
`bowtie -f -S -k 1 -p $process -v 1  $step9/seq_ViralHits_and_NoHits.fasta $step2/trimmed_filtered_gt15.fasta | awk -F '\t' '{if( \$2 != 4) print \$0}'  > $step10/reads.VS.contigs_virus_and_nohit.sam`;
`$binary/Z-score.bothstrands.pl -sam $step10/reads.VS.contigs_virus_and_nohit.sam -p $step10/reads.VS.contigs_virus_and_nohit`;
`R --no-save $step10/reads.VS.contigs_virus_and_nohit.zscore.tab $step10/plots/reads.VS.contigs_virus_and_nohit.all < /home/bioinfo/eric_bin/R/heatmap_correlation_VISA.R 2>/dev/null`;


### identifying contigs with siRNA and piRNA aignature
`$binary/samStatistics_v3.pl -sam $step10/reads.VS.contigs_virus_and_nohit.sam -fa $step9/seq_ViralHits_and_NoHits.fasta -p $step10/reads.VS.contigs_virus_and_nohit --profile`;


#formatting candidate contigs for bowtie
`bowtie-build $step10/reads.VS.contigs_virus_and_nohit.withSiRNA_and_PiRNA.fasta $step10/reads.VS.contigs_virus_and_nohit.withSiRNA_and_PiRNA.fasta > /dev/null `;

`bowtie-build $step10/reads.VS.contigs_virus_and_nohit.withSiRNA.fasta $step10/reads.VS.contigs_virus_and_nohit.withSiRNA.fasta > /dev/null `;


###bowtie 
`bowtie -f -S -k 1 -p $process -v 1  $step10/reads.VS.contigs_virus_and_nohit.withSiRNA.fasta $step2/trimmed_filtered_gt15.fasta | awk -F '\t' '{if( \$2 != 4) print \$0}'  > $step10/reads.VS.contigs_virus_and_nohit.siRNAs.sam`;

`bowtie -f -S -k 1 -p $process -v 1  $step10/reads.VS.contigs_virus_and_nohit.withSiRNA_and_PiRNA.fasta $step2/trimmed_filtered_gt15.fasta | awk -F '\t' '{if( \$2 != 4) print \$0}'  > $step10/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs.sam`;


print "\n\n Creating plots \n\n";
###creating folder to store plots
`mkdir $step10/plots`;

#plotting distribution and density plots
`$binary/plotMappingDataPerBasePreference.pl -sam $step10/reads.VS.contigs_virus_and_nohit.siRNAs.sam -s $si -e $se -fa $step10/reads.VS.contigs_virus_and_nohit.withSiRNA.fasta -pace 1 -p $step10/plots/reads.VS.contigs_virus_and_nohit.withSiRNA --profile --pattern`;

`$binary/plotMappingDataPerBasePreference.pl -sam $step10/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs.sam -s $si -e $se -fa $step10/reads.VS.contigs_virus_and_nohit.withSiRNA_and_PiRNA.fasta -pace 1 -p $step10/plots/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs --profile --pattern`;


##calculating mean pattern
`$binary/calcPatternInSamFile.pl -s $step10/reads.VS.contigs_virus_and_nohit.siRNAs.sam -o $step10/plots/reads.VS.contigs_virus_and_nohit.siRNAs.pattern`;

`$binary/calcPatternInSamFile.pl -s $step10/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs.sam -o $step10/plots/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs.pattern`;


#Plot geral distribution of reads in contigs with siRNA and siRNA+piRNAs
`$binary/plotGeralDistributionPerBaseByReads.pl -sam $step10/reads.VS.contigs_virus_and_nohit.siRNAs.sam -s $si -e $se -p $step10/plots/reads.VS.contigs_virus_and_nohit.siRNAs.geral_distribution --plot`;

`$binary/plotGeralDistributionPerBaseByReads.pl -sam $step10/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs.sam -s $si -e $se -p $step10/plots/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs.geral_distribution --plot`;


#Generating clusterized heatmaps 

`$binary/Z-score.bothstrands.pl -sam $step10/reads.VS.contigs_virus_and_nohit.siRNAs.sam -p $step10/reads.VS.contigs_virus_and_nohit.siRNAs`;

`$binary/Z-score.bothstrands.pl -sam $step10/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs.sam -p $step10/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs`;

`R --no-save $step10/reads.VS.contigs_virus_and_nohit.siRNAs.zscore.tab $step10/plots/reads.VS.contigs_virus_and_nohit.siRNAs < /home/bioinfo/eric_bin/R/heatmap_correlation_VISA.R 2>/dev/null`;

`R --no-save $step10/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs.zscore.tab $step10/plots/reads.VS.contigs_virus_and_nohit.siRNAs_and_piRNAs < /home/bioinfo/eric_bin/R/heatmap_correlation_VISA.R 2>/dev/null`;

##############################################################################################################

if (defined($clean)) {
	print "\n\n[ CLEAN LARGE INTERMEDIATE FILES ] \n\n";
	
    `$binary/samStatistics.pl -sam $step4/trimmed_filtered_gt15_X_vector.sam > $step4/trimmed_filtered_gt15_X_vector.sam.stats `;
	`$binary/samStatistics.pl -sam $step4/readsXbacters.sam  > $step4/readsXbacters.sam.stats`;

	`rm -rf $step3/*.sam  */run* $step3_1/*.sam $step4/* $prefix*/*assemble*/run* $step0/* $step4/unmappedVectorReads.fasta $step2/*.fasta `;

}


close(LOG);
close(REPORT);
close(metrics);
close(interest);

print "\n\n[[ Small RNA metagenomics pipeline is finished! ]] \n\n";
print "\n\n[[ Wake up...The Matrix has you...Follow the white rabbit... ]] \n\n";
