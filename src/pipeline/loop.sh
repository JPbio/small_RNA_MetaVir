# ls /small-rna-metavir/asset/libs/SRA_smallRNA/batch-01/*.fasta

skip_libs=("")
dir_root_libs=""
dir_root_runs=""

for file in `ls ${dir_root_libs}/*.fasta`
# for file in `ls ${dir_root_libs}/`
do

lib_id=`echo $file | cut -f 7 -d "/" | cut -f 1 -d "_"`
# echo "lib_id [cmd]: '${file} | cut -f 7 -d \"/\" | cut -f 1 -d \"_\"'"
# echo "lib_id: '${lib_id}'"

# Check if the current item should be skipped
skip=false

for skip_lib in "${skip_libs[@]}"; do
    if [[ "$lib_id" == "$skip_lib" ]]; then
        skip=true
        break
    fi
done

if [[ "$skip" == "true" ]]; then
    continue
fi

# Create results directory
exec_id=val-${lib_id}
dir_run=${dir_root_runs}/${exec_id}/

if [ ! -d "$dir_run" ]; then
    mkdir $dir_run
fi

# Run pipeline for the current lib
# echo "perl main.pl -fasta ${file} -hostgenome /small-rna-metavir/asset/refs/ref_hosts/refs_Aae_Aalb_Agam_culex.fasta -process 40 -si 15 -se 35 -size 20000 -exec-id val-${lib_id} -prefix ${lib_id} 2> ./runs/val-${lib_id}/val-${lib_id}.err.log"
# echo "perl main.pl -fasta ${file} -hostgenome /small-rna-metavir/asset/refs/ref_hosts/refs_Aae_Aalb_Agam_culex.fasta -process 40 -si 15 -se 35 -size 20000 -exec-id ${exec_id} -prefix ${lib_id} 2> ${dir_run}/${exec_id}.err.log"

perl main.pl -fasta ${file} -hostgenome /small-rna-metavir/asset/refs/ref_hosts/refs_Aae_Aalb_Agam_culex.fasta -process 40 -si 15 -se 35 -size 20000 -exec-id ${exec_id} -prefix ${lib_id} 2> ${dir_run}/${exec_id}.err.log
done
