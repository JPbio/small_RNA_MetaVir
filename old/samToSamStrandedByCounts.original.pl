use Bio::SeqIO;  
use Bio::Seq::Quality;
use Getopt::Long;

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

GetOptions ("sam=s" => \$sam,"fa=s" => \$fasta,"p=s" => \$prefix,"r=s" => \$reads,"fastq!" => \$fastq, "profile!" => \$profile,"h!" => \$help);

if ($help) {
        die $usage;
}

if (not(defined($sam))) {
        die "\nGive an input file name",$usage;
}
if (not(defined($prefix))) {
        die "\nGive an prefix file name",$usage;
}
if ( defined($profile) and not(defined($fasta))) {
        die "\nGive an fasta file ",$usage;
}

if ( not(defined($reads))) {
        die "\nGive an reads input file ",$usage;
}

if(defined($fastq)){
	$type="-q";
}else{
		$type="-f";
}

`samStatistics_v2.pl -sam $sam --counts > $prefix.stats`;

`defineStrandBasedOnCounts.pl -c $sam.counts -f $fasta -p $prefix.stranded.fasta`;

`bowtie-build $prefix.stranded.fasta $prefix.stranded.fasta`;

`bowtie $type -S  -p 40 -k 10  -v 1  $prefix.stranded.fasta $reads   | awk -F '\t' '{if( \$2 != 4) print \$0}' > $prefix.v1.k10.sam`;

`rm -rf  $prefix.stats $sam.counts`;
