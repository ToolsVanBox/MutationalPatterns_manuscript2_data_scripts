library(devtools)
library(tidyverse)
library(SomaticSignatures)
library(signeR)
library(microbenchmark)
library(scales)
setwd("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/Freek_MutationalPatterns/MutationalPatterns/")
load_all()



summarize_benchmark = function(benchmark){
    tb = benchmark %>% 
        as.data.frame() %>% 
        dplyr::mutate(expr = as.character(expr),
                      time = time / 1e9) %>% 
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

# Get context names
feature_names <- rownames(readRDS(system.file("states/mut_mat_data.rds",
                                              package = "MutationalPatterns"
)))

# Get mut_mat
res_l = readRDS("simulations/muts_per_sample/simulated_data/res_l.rds")
mut_mat = res_l[[2]]$mut_mat
colnames(mut_mat) = as.character(seq_len(ncol(mut_mat)))
rownames(mut_mat) = feature_names


# Create functions to run with benchmark
run_extraction_mutationalpatterns_regular = function(){
    nmf_res = extract_signatures(mut_mat, rank = 30, nmf_type = "regular")
}

run_extraction_mutationalpatterns_bayes = function(){
    nmf_res = extract_signatures(mut_mat, rank = 30, nmf_type = "variational_bayes", fudge = 0.005)
}

run_extraction_somaticsignatures_nmf = function(){
    # Remove 0 row to prevent error.
    nmf_res = identifySignatures(mut_mat[!rowSums(mut_mat) == 0,], 30, nmfDecomposition, nrun = 200)
}

run_extraction_somaticsignatures_pca = function(){
    nmf_res = identifySignatures(mut_mat, 30, pcaDecomposition, nrun = 200)
}

run_extraction_signer = function(){
    nmf_res = signeR(mut_mat, nsig = 30, samples = "columns")
}

#Run benchmark
bench = microbenchmark("MutationalPatterns_regular" = {run_extraction_mutationalpatterns_regular()},
                       "MutationalPatters_bayes" = {run_extraction_mutationalpatterns_bayes()},
                       "SomaticSignatures_nmf" = {run_extraction_somaticsignatures_nmf()},
                       "SomaticSignatures_pca" = {run_extraction_somaticsignatures_pca()},
                       "signeR" = {run_extraction_signer()}, 
                       times = 1)
bench_tb = summarize_benchmark(bench)

# Add sigprofiler
sigprofiler_tb = tibble("expr" = "SigProfiler", "mean" = 172815, "sdev" = NA, "n" = 1, "error" = NA, "lower" = NA, "upper" = NA)
bench_tb = rbind(bench_tb, sigprofiler_tb)

# Save results
write_tsv(bench_tb, "benchmarks/extraction.txt")

#### How long did hpc run?: sacct --format=Jobname,Elapsed --j 4122732
### 2-00:00:15: is 172815 seconds
# Plot the benchmark
extraction_bench_fig = ggplot(bench_tb, aes(x = expr, y = mean, fill = expr)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymax = upper, ymin = lower), width = 0.25) +
    scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000), label = format_format(scientific = FALSE)) +
    labs(y = "Run time (seconds)", x = "") +
    theme_classic() +
    theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90),
          panel.grid.major.y = element_line(colour = "black"))
ggsave("benchmarks/extraction_fig.pdf", extraction_bench_fig)
