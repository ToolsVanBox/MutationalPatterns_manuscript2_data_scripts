#! /bin/bash

#SBATCH --job-name=sigprofiler_extraction
#SBATCH -t 168:00:00
#SBATCH --mem=25G
#SBATCH --cpus-per-task 4
#SBATCH --mail-user=f.m.manders@prinsesmaximacentrum.nl
#SBATCH --mail-type=ALL

# Activate venv
source /hpc/pmc_vanboxtel/projects/Freek_sigprofiler/sigprofiler2/bin/activate

export PYTHONDONTWRITEBYTECODE=1

# Run sigprofiler
python3 /hpc/pmc_vanboxtel/projects/Freek_sigprofiler/SigProfiler.py --input $1 --output_dir $2 --context_type $3 --penalty $4  --min_nr_signatures 30 --max_nr_signatures 30
