#! /bin/bash

#SBATCH --job-name=sigprofiler
#SBATCH -t 144:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task 4
#SBATCH --mail-user=f.m.manders@prinsesmaximacentrum.nl
#SBATCH --mail-type=ALL

# Activate venv
source /hpc/pmc_vanboxtel/projects/Freek_sigprofiler/sigprofiler2/bin/activate

export PYTHONDONTWRITEBYTECODE=1

# Run sigprofiler
python3 /hpc/pmc_vanboxtel/projects/Freek_sigprofiler/SigProfiler.py --input $1 --output_dir $2 --context_type $3 --penalty $4
