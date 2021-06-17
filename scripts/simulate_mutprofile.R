library(devtools)
library(tidyverse)
library(deconstructSigs)
library(decompTumor2Sig)
setwd("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/Freek_MutationalPatterns/MutationalPatterns/")
load_all()
source("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/MutationalPatterns_manuscript2_data_scripts/scripts/simulation_data_functions.R")
source("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/MutationalPatterns_manuscript2_data_scripts/scripts/simulation_plot_functions.R")
source("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/MutationalPatterns_manuscript2_data_scripts/scripts/simulation_costs_functions.R")
source("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/MutationalPatterns_manuscript2_data_scripts/scripts/simulation_deconstructsigs_functions.R")
source("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/MutationalPatterns_manuscript2_data_scripts/scripts/simulation_mutpatterns_functions.R")
source("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/MutationalPatterns_manuscript2_data_scripts/scripts/simulation_sigprofiler_functions.R")
source("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/MutationalPatterns_manuscript2_data_scripts/scripts/simulation_decomptumor2sig_functions.R")

# Set up main output directory
set.seed(42)
out_dir = "~/surfdrive/Shared/Projects/Freek/mutpatterns/simulations"
if (!dir.exists(out_dir)){
    dir.create(out_dir)
}
setwd(out_dir)

# Get context names
feature_names <- rownames(readRDS(system.file("states/mut_mat_data.rds",
                                              package = "MutationalPatterns"
)))


#Copies all the simulated data to the hpc, where sigprofiler could be used.
######rsync -am --include='mut_mat*.txt' --include='*/' --exclude='*' ~/surfdrive/Shared/Projects/Freek/mutpatterns/simulations/* ~/hpc/pmc_vanboxtel/projects/Freek_sigprofiler/experiments/
#On hpc: bash start_Sigprofiler_all_samples.sh 
######rsync -am --include='*sigprofiler_output*' --include='*/' --exclude='*' ~/hpc/pmc_vanboxtel/projects/Freek_sigprofiler/experiments/* ~/surfdrive/Shared/Projects/Freek/mutpatterns/simulations/
######rsync -am --include='*Activities_refit.txt*' --include='*/' --exclude='*' ~/hpc/pmc_vanboxtel/projects/Freek_sigprofiler/experiments/* ~/surfdrive/Shared/Projects/Freek/mutpatterns/simulations/
#### How long did hpc run?: sacct --format=Jobname,Elapsed --starttime 2021-04-01


# Get signatures
signatures = get_known_signatures(source = "COSMIC_v3.1")
signatures = signatures[,c(1:30)]

# Set max_delta range and number of muts per signature
max_deltas = 1/2^c(seq(3,20))
nr_muts_l = list("50" = rep(50, 4),
                 "100" = rep(100, 4),
                 "500" = rep(500, 4),
                 "1000" = rep(1000, 4))
create_write_simulated_data(nr_muts_l, signatures, nr_samples = 300, "muts_per_sample", feature_names)
do_simulation_analysis_mutpatterns(signatures, "muts_per_sample", max_deltas)
analyze_cost_sigprofiler("muts_per_sample")
analyze_cost_deconstructsigs(signatures, "muts_per_sample", feature_names)
analyze_cost_decomptumor2sigs(signatures, "muts_per_sample", feature_names)
plot_simulation("muts_per_sample")

# Different amounts of minor signature
# Set max_delta range and number of muts per signature
max_deltas = 1/2^c(seq(3,20))
nr_muts_l = list("1000" = c(1000, 1000, 1000, 1000),
                 "minor_500" = c(1000, 1000, 1000, 500),
                 "minor_100" = c(1000, 1000, 1000, 100),
                 "minor_50" = c(1000, 1000, 1000, 50))
create_write_simulated_data(nr_muts_l, signatures, nr_samples = 300, "minor_sig", feature_names)
do_simulation_analysis_mutpatterns(signatures, "minor_sig", max_deltas)
analyze_cost_sigprofiler("minor_sig")
analyze_cost_deconstructsigs(signatures, "minor_sig", feature_names)
analyze_cost_decomptumor2sigs(signatures, "minor_sig", feature_names)
plot_simulation("minor_sig")


#Different numbers of signatures
# Set max_delta range and number of muts per signature
max_deltas = 1/2^c(seq(3,20))
nr_muts_l = list("2" = rep(500, 2),
                 "4" = rep(500, 4),
                 "6" = rep(500, 6),
                 "8" = rep(500, 8))
create_write_simulated_data(nr_muts_l, signatures, nr_samples = 300, "nrsigs", feature_names)
do_simulation_analysis_mutpatterns(signatures, "nrsigs", max_deltas)
analyze_cost_sigprofiler("nrsigs")
analyze_cost_deconstructsigs(signatures, "nrsigs", feature_names)
analyze_cost_decomptumor2sigs(signatures, "nrsigs", feature_names)
plot_simulation("nrsigs")

# Get indel signatures
signatures = get_known_signatures("indel")
feature_names = read_tsv("~/Downloads/COSMIC_v3.2_ID_GRCh37.txt")$Type
# feature_names <- rownames(readRDS(system.file("states/blood_indel_counts.rds",
#                                               package = "MutationalPatterns"
# )))


# Set max_delta range and number of muts per signature
max_deltas = 1/2^c(seq(3,20))
nr_muts_l = list("50" = rep(50, 4),
                 "100" = rep(100, 4),
                 "500" = rep(500, 4),
                 "1000" = rep(1000, 4))
create_write_simulated_data(nr_muts_l, signatures, nr_samples = 300, "indel_muts_per_sample", feature_names)
do_simulation_analysis_mutpatterns(signatures, "indel_muts_per_sample", max_deltas)
analyze_cost_sigprofiler("indel_muts_per_sample")
plot_simulation("indel_muts_per_sample")

# Different amounts of minor signature
# Set max_delta range and number of muts per signature
max_deltas = 1/2^c(seq(3,20))
nr_muts_l = list("1000" = c(1000, 1000, 1000, 1000),
                 "minor_500" = c(1000, 1000, 1000, 500),
                 "minor_100" = c(1000, 1000, 1000, 100),
                 "minor_50" = c(1000, 1000, 1000, 50))
create_write_simulated_data(nr_muts_l, signatures, nr_samples = 300, "indel_minor_sig", feature_names)
do_simulation_analysis_mutpatterns(signatures, "indel_minor_sig", max_deltas)
analyze_cost_sigprofiler("indel_minor_sig")
plot_simulation("indel_minor_sig")

#Different numbers of signatures
# Set max_delta range and number of muts per signature
max_deltas = 1/2^c(seq(3,20))
nr_muts_l = list("2" = rep(500, 2),
                 "4" = rep(500, 4),
                 "6" = rep(500, 6),
                 "8" = rep(500, 8))
create_write_simulated_data(nr_muts_l, signatures, nr_samples = 300, "indel_nrsigs", feature_names)
do_simulation_analysis_mutpatterns(signatures, "indel_nrsigs", max_deltas)
analyze_cost_sigprofiler("indel_nrsigs")
plot_simulation("indel_nrsigs")


# variables: delta, nr_muts_per_sig, nr_signatures total, ratio_sigs.
