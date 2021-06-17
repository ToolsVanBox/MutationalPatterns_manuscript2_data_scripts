#! /bin/bash

#SBATCH --job-name=sigprofiler_AHH1
#SBATCH -t 96:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task 4
#SBATCH --mail-user=f.m.manders@prinsesmaximacentrum.nl
#SBATCH --mail-type=ALL
#SBATCH -e /hpc/pmc_vanboxtel/projects/Freek_sigprofiler/AHH1/AHH1_sigprofiler_error.txt
#SBATCH -o /hpc/pmc_vanboxtel/projects/Freek_sigprofiler/AHH1/AHH1_sigprofiler_output.txt

# Activate venv
source /hpc/pmc_vanboxtel/projects/Freek_sigprofiler/sigprofiler2/bin/activate

export PYTHONDONTWRITEBYTECODE=1

# Run sigprofiler
input=/hpc/pmc_vanboxtel/projects/Freek_sigprofiler/AHH1/mut_mat_sigprofiler.txt
output=/hpc/pmc_vanboxtel/projects/Freek_sigprofiler/AHH1/

# Run sigprofiler
python3 /hpc/pmc_vanboxtel/projects/Freek_sigprofiler/SigProfiler.py --input ${input} --output_dir ${output} --context_type SBS96 --penalty 0.05
