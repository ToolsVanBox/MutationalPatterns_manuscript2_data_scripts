library(devtools)
library(tidyverse)
library(deconstructSigs)
library(decompTumor2Sig)
library(microbenchmark)
source("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/MutationalPatterns_manuscript2_data_scripts/scripts/simulation_deconstructsigs_functions.R")
source("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/MutationalPatterns_manuscript2_data_scripts/scripts/simulation_decomptumor2sig_functions.R")
setwd("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/Freek_MutationalPatterns/MutationalPatterns/")
load_all()



summarize_benchmark = function(benchmark){
    tb = benchmark %>% 
        as.data.frame() %>% 
        dplyr::mutate(expr = as.character(expr),
                      time = time / 1e6) %>% 
        dplyr::group_by(expr) %>% 
        dplyr::summarise(mean = mean(time), 
                         sdev = sd(time),
                         n = dplyr::n()) %>% 
        dplyr::mutate(error = qt(0.975, df = n-1)*sdev/sqrt(n),
                      lower = mean - error,
                      upper = mean + error)
    return(tb)
}




setwd("~/surfdrive/Shared/Projects/Freek/mutpatterns/")

# Get mut_mat
res_l = readRDS("simulations/muts_per_sample/simulated_data/res_l.rds")
mut_mat = res_l[[2]]$mut_mat[,1:10]

# Get signatures
signatures = get_known_signatures(source = "COSMIC_v3.1")
signatures = signatures[,c(1:30)]

# Get context names
feature_names <- rownames(readRDS(system.file("states/mut_mat_data.rds",
                                              package = "MutationalPatterns"
)))






# Run deconstructsigs
decon_input_l = get_correct_format_deconstructsigs(mut_mat, signatures, feature_names)
run_deconstructsigs = function(){
    contri_rel = purrr::map(seq_len(nrow(decon_input_l$input)), ~whichSignatures(tumor.ref = decon_input_l$input, 
                                                                        signatures.ref = decon_input_l$signatures, 
                                                                        sample.id = .x,
                                                                        signature.cutoff = 0.06,
                                                                        contexts.needed = TRUE))
}

# Run mutationalPatterns
run_mutationalpatterns = function(){
    fit_res = fit_to_signatures_strict(mut_mat, signatures)
}


# Run decomptumor2sig
decomp_input_l = get_correct_format_decomptumor2sig(mut_mat, signatures, feature_names)
run_decomptumor2sig = function(){
    contri_rel = decomposeTumorGenomes(decomp_input_l$decomp_input, 
                                       decomp_input_l$signatures_l, 
                                       minExplainedVariance = 0.95,
                                       greedySearch = TRUE)
}

bench = microbenchmark("MutationalPatterns" = {run_mutationalpatterns()},
               "deconstructSigs" = {run_deconstructsigs()},
               "decompTumor2Sig" = {run_decomptumor2sig()}, 
               times = 10)
bench_tb = summarize_benchmark(bench)

sigprofiler_tb = tibble("expr" = "SigProfiler", "mean" = 107640, "sdev" = NA, "n" = 1, "error" = NA, "lower" = NA, "upper" = NA)
bench_tb = rbind(bench_tb, sigprofiler_tb)
write_tsv(bench_tb, "benchmarks/refitting.txt")

refitting_bench_fig = ggplot(bench_tb, aes(x = expr, y = mean, fill = expr)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymax = upper, ymin = lower), width = 0.25) +
    scale_y_log10() +
    labs(y = "Run time (milli seconds)", x = "") +
    theme_classic() +
    theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90))
ggsave("benchmarks/refitting_fig.pdf", refitting_bench_fig)

#### How long did hpc run?: sacct --format=Jobname,Elapsed --j 3802017
