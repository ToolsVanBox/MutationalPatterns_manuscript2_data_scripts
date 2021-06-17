#! /bin/bash

#SBATCH --job-name=sigprofiler
#SBATCH -t 96:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task 4
#SBATCH --mail-user=f.m.manders@prinsesmaximacentrum.nl
#SBATCH --mail-type=ALL
#SBATCH -e /hpc/pmc_vanboxtel/projects/Freek_sigprofiler/extract_realdata/extract_realdata_sigprofiler_error.txt
#SBATCH	-o /hpc/pmc_vanboxtel/projects/Freek_sigprofiler/extract_realdata/extract_realdata_sigprofiler_output.txt

# Activate venv
source /hpc/pmc_vanboxtel/projects/Freek_sigprofiler/sigprofiler/bin/activate

input=/hpc/pmc_vanboxtel/projects/Freek_sigprofiler/extract_realdata/AHH1_blokzijl2016_mut_mat.txt
output=/hpc/pmc_vanboxtel/projects/Freek_sigprofiler/extract_realdata/

# Run sigprofiler
python3 /hpc/pmc_vanboxtel/projects/Freek_sigprofiler/SigProfiler.py --input ${input} --output_dir ${output} --context_type SBS96 --penalty 0.001 --min_nr_signatures 3 --max_nr_signatures 3
