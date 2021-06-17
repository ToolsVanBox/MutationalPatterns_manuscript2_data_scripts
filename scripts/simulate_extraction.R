library(devtools)
library(tidyverse)
library(SomaticSignatures)
library(signeR)
setwd("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/Freek_MutationalPatterns/MutationalPatterns/")
load_all()
source("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/MutationalPatterns_manuscript2_data_scripts/scripts/simulation_extraction_functions.R")
source("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/MutationalPatterns_manuscript2_data_scripts/scripts/simulation_extraction_plot_functions.R")


# Set up main output directory
set.seed(42)
out_dir = "~/surfdrive/Shared/Projects/Freek/mutpatterns/simulations_extraction"
if (!dir.exists(out_dir)){
    dir.create(out_dir)
}
setwd(out_dir)

# Get context names
feature_names <- rownames(readRDS(system.file("states/mut_mat_data.rds",
                                              package = "MutationalPatterns"
)))

# Get signatures
signatures = get_known_signatures(source = "COSMIC_v3.1")
signatures = signatures[,c(1:30)]


## Run sigprofiler on the hpc.
## Copy data from hpc to local.
######rsync -am --include='*De-Novo_Signatures.txt*' --include='*/' --exclude='*' ~/hpc/pmc_vanboxtel/projects/Freek_sigprofiler/extract_simulations/experiments/* ~/surfdrive/Shared/Projects/Freek/mutpatterns/simulations_extraction


# Mutations per sample
do_simulation_extraction_mutpatterns(signatures, "muts_per_sample")
do_simulation_extraction_somaticsignatures(signatures, "muts_per_sample")
do_simulation_extraction_signer(signatures, "muts_per_sample", feature_names)
do_simulation_extraction_sigprofiler(signatures, "muts_per_sample")
plot_simulation_extraction("muts_per_sample")

# Signature with less contribution than others
do_simulation_extraction_mutpatterns(signatures, "minor_sig")
do_simulation_extraction_somaticsignatures(signatures, "minor_sig")
do_simulation_extraction_signer(signatures, "minor_sig", feature_names)
do_simulation_extraction_sigprofiler(signatures, "minor_sig")
plot_simulation_extraction("minor_sig")

# Number of signatures
do_simulation_extraction_mutpatterns(signatures, "nrsigs")
do_simulation_extraction_somaticsignatures(signatures, "nrsigs")
do_simulation_extraction_signer(signatures, "nrsigs", feature_names)
do_simulation_extraction_sigprofiler(signatures, "nrsigs")
plot_simulation_extraction("nrsigs")


