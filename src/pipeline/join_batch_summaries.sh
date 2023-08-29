dir_root_libs="/small-rna-metavir/src/runs"
dir_summaries="/small-rna-metavir/src/run-summaries"

min_n=1
max_n=6

# Loop through batches
for ((i=min_n; i<=max_n; i++))
do
    # Check if batch exists
    batch_name="batch-0$i"
    batch_dir="${dir_root_libs}/$batch_name"

    if [ ! -d "$batch_dir" ]; then
        continue
    fi
    # echo $batch_dir

    # Loop through batch executions
    for exec_id in `ls ${batch_dir}`
    do

        # Chck if it's a valid execution path
        path_exec="${batch_dir}/$exec_id"
        
        if [ ! -d "$path_exec" ]; then
            continue
        fi

        # Create sumary dir
        path_summary_dir="$dir_summaries/$batch_name-$exec_id"
        # echo "path_summary_dir: '$path_summary_dir'"

        if [ ! -d "$path_summary_dir" ]; then
            mkdir $path_summary_dir
        fi

        # Copy summary files
        path_exec_csv="$path_exec/13_virus_eve_classif/$exec_id-viral-eve.csv"
        path_exec_summary_lines="$path_exec/$exec_id.summary.lines"
        path_exec_summary="$path_exec/$exec_id.summary.yml"

        cp $path_exec_csv $path_summary_dir/
        cp $path_exec_summary $path_summary_dir/
        cp $path_exec_summary_lines $path_summary_dir/
        # echo "cp $path_exec_summary $path_summary_dir/"
        # echo "cp $path_exec_summary_lines $path_summary_dir/"
        
    done
done