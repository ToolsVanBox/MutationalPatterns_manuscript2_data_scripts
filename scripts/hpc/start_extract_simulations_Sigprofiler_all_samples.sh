cd /hpc/pmc_vanboxtel/projects/Freek_sigprofiler/

for DIR in experiments/*; do

    context_type="SBS96"
    if [[ ${DIR} =~ "indel" ]]; then
        context_type="ID83"
    fi

    output_dir="/hpc/pmc_vanboxtel/projects/Freek_sigprofiler/extract_simulations/"${DIR}/"sigprofiler_output/"
    if [ ! -d ${output_dir} ]; then
        mkdir ${output_dir}
    fi


    for mut_mat in ${DIR}/simulated_data/*.txt; do
        name=$(basename $mut_mat)
        name=${name#"mut_mat_"}
        name=${name%".txt"}
        full_name_mut_mat=$PWD/$mut_mat

        for penalty in 0.01; do
            output_dir_final=${output_dir}"/sigprofiler_output_"${name}"_nnls_remove_"${penalty}
            time=$(date +'%H%M%S')
            output_file=${output_dir_final}/${context_type}/Suggested_Solution/COSMIC_${context_type}_Decomposed_Solution/Activities/COSMIC_${context_type}_Activities_refit.txt
            if [[ ! -f "${output_file}" ]]; then
                sbatch -e "${output_dir_final}/run_Sigprofiler_error.txt" -o "${output_dir_final}/run_Sigprofiler_output.txt" extract_simulations/run_extract_simulation_Sigprofiler.sh ${full_name_mut_mat} ${output_dir_final} ${context_type} ${penalty}
                echo "Started run for experiment: ${name} in dir: ${DIR} with context type: ${context_type} and nnls_remove_penalty: ${penalty}"
                sleep 1
            fi
        done
    done
done

###ADD for loop for stringency!!!!!
