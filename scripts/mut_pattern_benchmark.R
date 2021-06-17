library(microbenchmark)
library(tidyverse)
out_dir = "~/surfdrive/Shared/Projects/Freek/mutpatterns/benchmarks_newvsold/"
if (!dir.exists(out_dir)){
    dir.create(out_dir)
}
setwd(out_dir)


ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
 
#Load old version of mutational patterns
library(MutationalPatterns)

read_vcfs_as_granges_old = read_vcfs_as_granges
mut_matrix_old = mut_matrix
mut_matrix_stranded_old = mut_matrix_stranded

#Load new version of mutational patterns
library(devtools)
load_all("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/Freek_MutationalPatterns/MutationalPatterns/")

#Load test data
vcf_files <- list.files(system.file("extdata", package = "MutationalPatterns"),
                        pattern = "sample.vcf", full.names = TRUE
)

sample_names <- c(
    "colon1", "colon2", "colon3",
    "intestine1", "intestine2", "intestine3",
    "liver1", "liver2", "liver3"
)

grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)

#Perform benchmarks
bench_reading = microbenchmark("Old" = {read_vcfs_as_granges_old(vcf_files, sample_names, ref_genome)},
                                   "New" = {read_vcfs_as_granges(vcf_files, sample_names, ref_genome, type = "all")})
bench_mut_mat = microbenchmark("Old" = {mut_matrix_old(grl, ref_genome)},
                                   "New" = {mut_matrix(grl, ref_genome)})
bench_mut_mat_strand = microbenchmark("Old" = {mut_matrix_stranded_old(grl, ref_genome, genes_hg19)},
                                          "New" = {mut_matrix_stranded(grl, ref_genome, genes_hg19)})

#Write results
write_tsv(bench_reading, "benchmark_reading.txt")
write_tsv(bench_mut_mat, "benchmark_mut_mat.txt")
write_tsv(bench_mut_mat_strand, "benchmark_mut_mat_strand.txt")


mut_mat_means = bench_mut_mat %>% 
    dplyr::group_by(expr) %>% 
    dplyr::summarise(mean = mean(time))
mut_mat_means$mean[2] / mut_mat_means$mean[1]

mut_mat_means_strand = bench_mut_mat_strand %>% 
    dplyr::group_by(expr) %>% 
    dplyr::summarise(mean = mean(time))
mut_mat_means_strand$mean[2] / mut_mat_means_strand$mean[1]


#Plot
reading_fig = autoplot(bench_reading)
ggsave("benchmark_reading.pdf", reading_fig)
mut_mat_fig = autoplot(bench_mut_mat)
ggsave("benchmark_mut_mat.pdf", mut_mat_fig)
mut_mat_strand_fig = autoplot(bench_mut_mat_strand)
ggsave("benchmark_mut_mat_strand.pdf", mut_mat_strand_fig)



#Mut matrix scaling
split_in_samples = function(gr, nr_muts, nr_samples = 100){
    
    #Sample nr_muts mutations
    gr = gr[sample.int(length(gr), nr_muts, replace = TRUE)]
    
    #Randomly reorder
    order = sample(nr_muts, nr_muts)
    gr = gr[order]
    
    #Split in nr_samples samples
    sample_i = rep(seq_len(nr_samples), each = nr_muts/nr_samples)
    grl = split(gr, sample_i)
    return(grl)
}

summarize_benchmark = function(benchmark){
    tb = benchmark %>% 
        as.data.frame() %>% 
        dplyr::mutate(expr = as.numeric(as.character(expr)),
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

get_linear_scale = function(scaling_tb){
    min_mean = scaling_tb %>% 
        dplyr::filter(expr == min(expr)) %>% 
        dplyr::pull(mean)
    
    min_expr = min(scaling_tb$expr)
    max_expr =  max(scaling_tb$expr)
    ratio = max_expr / min_expr
    max_mean = min_mean * ratio
    tb = tibble::tibble("expr" = c(min_expr, max_expr), "mean" = c(min_mean, max_mean))
    return(tb)
}

#Determine file names
hmf_overview_all = suppressWarnings(read_tsv("~/surfdrive/Shared/Projects/Freek/mutpatterns/hmf/hmf_metadata.tsv"))
hmf_overview = dplyr::filter(hmf_overview_all, primaryTumorLocation == "Skin" & cancerSubtype == "Melanoma")
hmf_overview = hmf_overview[1:200,]
hmf_dir = "~/hpc/pmc_vanboxtel/projects/Axel_GenoEcoli/HMF_analysis/somatics"
vcf_fnames = file.path(hmf_dir, hmf_overview$setName, paste0(hmf_overview$sampleId, ".purple.somatic.vcf"))

#Set sample names
sample_names = hmf_overview$sampleId

#Read data
grl = read_vcfs_as_granges(vcf_fnames, sample_names, ref_genome, type = "all")
grl = purrr::map(as.list(grl), function(gr){gr = gr[gr$FILTER == "PASS"]}) %>% 
    purrr::map(.remove_multi_alts_variants) %>% 
    GRangesList()
gr = unlist(grl)

#SNVs (Maybe don't do the split in samples. It's probably only usefull when comparing to other packages.)
snv_gr = gr[.find_snv(gr) & width(gr$REF) == 1]


snv_grl_100 = split_in_samples(snv_gr, 1e2)
snv_grl_1000 = split_in_samples(snv_gr, 1e3)
snv_grl_10000 = split_in_samples(snv_gr, 1e4)
snv_grl_100000 = split_in_samples(snv_gr, 1e5)
snv_grl_1000000 = split_in_samples(snv_gr, 1e6)
snv_grl_10000000 = split_in_samples(snv_gr, 1e7)


snv_scaling_bench = microbenchmark(
               "100" = {mut_matrix(snv_grl_100, ref_genome)},
               "1000" = {mut_matrix(snv_grl_1000, ref_genome)},
               "10000" = {mut_matrix(snv_grl_10000, ref_genome)},
               "100000" = {mut_matrix(snv_grl_100000, ref_genome)},
               "1000000" = {mut_matrix(snv_grl_1000000, ref_genome)},
               "10000000" = {mut_matrix(snv_grl_10000000, ref_genome)},
               times = 10)
write_tsv(snv_scaling_bench, "snv_scaling_bench.txt")
snv_scaling_bench %>% 
    dplyr::group_by(expr) %>% 
    dplyr::summarise(mean = mean(time))

snv_scaling_tb = summarize_benchmark(snv_scaling_bench)
linear_scale = get_linear_scale(snv_scaling_tb)

snv_scaling_fig = ggplot(snv_scaling_tb, aes(y = mean, x = expr)) +
    geom_line() +
    geom_errorbar(width = 0.2, aes(ymin = lower, ymax = upper)) +
    geom_line(data = linear_scale, linetype = "dashed", colour = "red") +
    scale_y_log10() +
    scale_x_log10() +
    labs(y = "seconds", x = "Nr SNVs (100 samples)") +
    theme_bw() +
    coord_cartesian(xlim = c(100, 10000000), ylim = c(0.1, 3000))
ggsave("snv_scaling.pdf", snv_scaling_fig)


# Look at scaling x mutations between a variable number of samples
snv_grl_1s = split_in_samples(snv_gr, 1e5, 1)
snv_grl_10s = split_in_samples(snv_gr, 1e5, 10)
snv_grl_100s = split_in_samples(snv_gr, 1e5, 100)
snv_grl_1000s = split_in_samples(snv_gr, 1e5, 1000)


snv_sample_scaling_bench = microbenchmark(
    "1" = {mut_matrix(snv_grl_1s, ref_genome)},
    "10" = {mut_matrix(snv_grl_10s, ref_genome)},
    "100" = {mut_matrix(snv_grl_100s, ref_genome)},
    "1000" = {mut_matrix(snv_grl_1000s, ref_genome)},
    times = 10)
write_tsv(snv_sample_scaling_bench, "snv_sample_scaling_bench.txt")

snv_sample_scaling_tb = summarize_benchmark(snv_sample_scaling_bench)
linear_scale = get_linear_scale(snv_sample_scaling_tb)

snv_sample_scaling_fig = ggplot(snv_sample_scaling_tb, aes(y = mean, x = expr)) +
    geom_line() +
    geom_errorbar(width = 0.2, aes(ymin = lower, ymax = upper)) +
    geom_line(data = linear_scale, linetype = "dashed", colour = "red") +
    scale_y_log10() +
    scale_x_log10() +
    labs(y = "seconds", x = "Nr samples (1e5 mutations)") +
    theme_bw() +
    coord_cartesian(xlim = c(1, 1000), ylim = c(1, 30))
ggsave("snv_sample_scaling.pdf", snv_sample_scaling_fig)






#Indels
indel_gr = gr[!.find_snv(gr)]
indel_grl_100 = split_in_samples(indel_gr, 1e2)
indel_grl_1000 = split_in_samples(indel_gr, 1e3)
indel_grl_10000 = split_in_samples(indel_gr, 1e4)
indel_grl_100000 = split_in_samples(indel_gr, 1e5)


indel_scaling_bench = microbenchmark(
    "100" = {count_indel_contexts(get_indel_context(indel_grl_100, ref_genome))},
    "1000" = {count_indel_contexts(get_indel_context(indel_grl_1000, ref_genome))},
    "10000" = {count_indel_contexts(get_indel_context(indel_grl_10000, ref_genome))},
    "100000" = {count_indel_contexts(get_indel_context(indel_grl_100000, ref_genome))},
    times = 10)
write_tsv(indel_scaling_bench, "indel_scaling_bench.txt")

indel_scaling_tb = summarize_benchmark(indel_scaling_bench)
linear_scale = get_linear_scale(indel_scaling_tb)


indel_scaling_fig = ggplot(indel_scaling_tb, aes(y = mean, x = expr)) +
    geom_line() +
    geom_errorbar(width = 0.2, aes(ymin = lower, ymax = upper)) +
    geom_line(data = linear_scale, linetype = "dashed", colour = "red") +
    scale_y_log10() +
    scale_x_log10() +
    labs(y = "seconds", x = "Nr indels (100 samples)") +
    theme_bw() +
    coord_cartesian(xlim = c(100, 100000), ylim = c(1, 30))
ggsave("indel_scaling.pdf", indel_scaling_fig)

#DBS
dbs_gr = gr[.find_snv(gr) & width(gr$REF) == 2]
dbs_grl_100 = split_in_samples(dbs_gr, 1e2)
dbs_grl_1000 = split_in_samples(dbs_gr, 1e3)
dbs_grl_10000 = split_in_samples(dbs_gr, 1e4)
dbs_grl_100000 = split_in_samples(dbs_gr, 1e5)



dbs_scaling_bench = microbenchmark(
    "100" = {count_dbs_contexts(get_dbs_context(dbs_grl_100))},
    "1000" = {count_dbs_contexts(get_dbs_context(dbs_grl_1000))},
    "10000" = {count_dbs_contexts(get_dbs_context(dbs_grl_10000))},
    "100000" = {count_dbs_contexts(get_dbs_context(dbs_grl_100000))},
    times = 10)
write_tsv(dbs_scaling_bench, "dbs_scaling_bench.txt")

dbs_scaling_tb = summarize_benchmark(dbs_scaling_bench)
linear_scale = get_linear_scale(dbs_scaling_tb)

dbs_scaling_fig = ggplot(dbs_scaling_tb, aes(y = mean, x = expr)) +
    geom_line() +
    geom_errorbar(width = 0.2, aes(ymin = lower, ymax = upper)) +
    geom_line(data = linear_scale, linetype = "dashed", colour = "red") +
    scale_y_log10() +
    scale_x_log10() +
    labs(y = "seconds", x = "Nr DBSs (100 samples)") +
    theme_bw() +
    coord_cartesian(xlim = c(100, 100000), ylim = c(1, 30))
ggsave("dbs_scaling.pdf", dbs_scaling_fig)
