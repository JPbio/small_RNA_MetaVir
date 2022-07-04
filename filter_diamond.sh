#!/bin/bash

#   -f: fasta with diamond Hits 
#   -t: diamond output   
#   -d: dir output; 
while getopts f:o:d: option                 
do
case "${option}"
in
f) fasta=${OPTARG};;
o) output=${OPTARG};;
d) dir=${OPTARG};;
esac
done   


cp ${fasta} temp_diamond_blastx_Hits.fasta
grep '>' temp_diamond_blastx_Hits.fasta  > temp_headers

while read p
do 
 h_query=`echo $p| cut -f1 -d " "| cut -f2 -d ">"`
 h_subject=`grep $h_query -A5 $output | grep ">"|tr ">" "|"| tr "/" "-"`
 sed -i "s/$h_query/$h_query $h_subject/g" temp_diamond_blastx_Hits.fasta
done < temp_headers

fasta_formatter -i temp_diamond_blastx_Hits.fasta  -o ${dir}/diamond_blastx_AllHits.fasta

grep -A1 -i "virus" ${dir}/diamond_blastx_AllHits.fasta | grep -v ^-- > ${dir}/diamond_blastx_Viral.fasta

grep '>' ${dir}/diamond_blastx_AllHits.fasta | grep -v -i "virus" > temp_nonviral

rm -f ${dir}/diamond_blastx_NonViral.fasta
while read p
do 	
a=`echo $p |cut -f1 -d " "`
grep -A1 $a ${dir}/diamond_blastx_AllHits.fasta
done < temp_nonviral >> ${dir}/diamond_blastx_NonViral.fasta

rm -f temp_diamond_blastx_Hits.fasta temp_headers temp_nonviral

echo "Thats all Folks!"
