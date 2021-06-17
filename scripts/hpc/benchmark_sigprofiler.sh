#! /bin/bash

#SBATCH --job-name=sigprofiler
#SBATCH -t 48:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task 4
#SBATCH --mail-user=f.m.manders@prinsesmaximacentrum.nl
#SBATCH --mail-type=ALL
#SBATCH -e /hpc/pmc_vanboxtel/projects/Freek_sigprofiler/bench/benchmark_sigprofiler_error.txt
#SBATCH	-o /hpc/pmc_vanboxtel/projects/Freek_sigprofiler/bench/benchmark_sigprofiler_output.txt

# Activate venv
source /hpc/pmc_vanboxtel/projects/Freek_sigprofiler/sigprofiler/bin/activate

input=/hpc/pmc_vanboxtel/projects/Freek_sigprofiler/bench/bench_sample.txt
output=/hpc/pmc_vanboxtel/projects/Freek_sigprofiler/bench/

# Run sigprofiler
python3 /hpc/pmc_vanboxtel/projects/Freek_sigprofiler/SigProfiler.py --input ${input} --output_dir ${output} --context_type SBS96 --penalty 0.001
