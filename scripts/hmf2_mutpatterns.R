library(devtools)
setwd("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Freek/MutationalPatterns/")
load_all()


ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(biomaRt)
library(readr)
library(tidyverse)
library("ccfindR")


get_regions_grl = function(){
    region <- getBM(
        attributes = c(
            "chromosome_name", "chromosome_start",
            "chromosome_end", "feature_type_name"
        ),
        mart = regulatory
    )
    region_l = split(region, region$feature_type_name)
    grl = purrr::map(region_l, function(x){GRanges(
        x$chromosome_name, IRanges(
            x$chromosome_start,
            x$chromosome_end))
        }) %>% 
        purrr::map(GenomicRanges::reduce) %>%
        GRangesList()
    seqlevels(grl, pruning.mode = "tidy") <- c(1:22, "X", "Y")
    grl <- sort(grl)
    seqlevelsStyle(grl) = "UCSC"
    return(grl)
}

# get_region_gr = function(region_name){
#     region <- getBM(
#         attributes = c(
#             "chromosome_name", "chromosome_start",
#             "chromosome_end", "feature_type_name"
#         ),
#         filters = "regulatory_feature_type_name",
#         values = region_name,
#         mart = regulatory
#     )
#     region_gr <- reduce(GRanges(
#         region$chromosome_name,
#         IRanges(
#             region$chromosome_start,
#             region$chromosome_end
#         )
#     ))
#     seqlevels(region_gr) <- c(1:22, "X", "Y")
#     region_gr <- sort(region_gr)
#     seqlevelsStyle(region_gr) = "UCSC"
#     region_grl = GRangesList(region_gr)
#     names(region_grl) = region_name
#     return(region_grl)
# }

create_matrices_sample = function(vcf_fname, sample_name, ref_genome, regions_grl){
    
    rds_fname = paste0("mutation_matrices/", "mut_matrices_", sample_name, ".rds")
    
    if (file.exists(rds_fname)){
        return(0)
    }
    
    grl = read_vcfs_as_granges(vcf_fname, sample_name, ref_genome, type = "all")
    grl = purrr::map(as.list(grl), function(gr){gr = gr[gr$FILTER == "PASS"]}) %>% 
        purrr::map(.remove_multi_alts_variants) %>% 
        GRangesList()
    
    snv_grl = purrr::map(as.list(grl), function(gr){gr = gr[.find_substitution(gr) & width(gr$REF) == 1]}) %>% 
        GRangesList()
    
    # Create regular and extended mutation matrices
    mut_mat = mut_matrix(snv_grl, ref_genome, extension = 1)
    mut_mat_ext_context = mut_matrix(snv_grl, ref_genome, extension = 2)
    mut_mat_ext_context2 = mut_matrix(snv_grl, ref_genome, extension = 3)
    
    # # Create promoter mutation matrix
    # grl_region = split_muts_region(snv_grl, promoter_grl, include_other = TRUE)
    # mut_mat_region_promoter <- mut_matrix(grl_region, ref_genome)
    # mut_mat_long_promoter <- lengthen_mut_matrix(mut_mat_region_promoter)
    # 
    # # Create ctcf mutation matrix
    # grl_region = split_muts_region(snv_grl, ctcf_grl, include_other = TRUE)
    # mut_mat_region_ctcf <- mut_matrix(grl_region, ref_genome)
    # mut_mat_long_ctcf <- lengthen_mut_matrix(mut_mat_region_ctcf)
    
    grl_incl_region = split_muts_region(snv_grl, regions_grl, include_other = TRUE)
    mut_mat_region <- mut_matrix(grl_incl_region, ref_genome)
    mut_mat_region_extended = mut_matrix(grl_incl_region, ref_genome, extension = 2)
    
    # Create indel mutation matrix
    indel_grl = get_mut_type(grl, "indel")
    indel_grl = get_indel_context(indel_grl, ref_genome)
    indel_counts = count_indel_contexts(indel_grl)
    
    res = list("mut_mat" = mut_mat, 
               "mut_mat_ext_context" = mut_mat_ext_context, 
               "mut_mat_ext_context2" = mut_mat_ext_context2, 
               # "mut_mat_region_promoter" = mut_mat_region_promoter,
               # "mut_mat_long_promoter" = mut_mat_long_promoter,
               # "mut_mat_region_ctcf" = mut_mat_region_ctcf,
               # "mut_mat_long_ctcf" = mut_mat_long_ctcf,
               "mut_mat_region" = mut_mat_region,
               "mut_mat_region_extended" = mut_mat_region_extended,
               "indel_counts" = indel_counts)
    saveRDS(res, rds_fname)
    message("Created mutation matrices for one sample.")
    invisible(0)
}





out_dir = "~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Freek/mutpatterns/hmf2/"
if (!dir.exists(out_dir)){
    dir.create(out_dir)
}
setwd(out_dir)

# Determine regulatory regions
regulatory <- useEnsembl(biomart="regulation", 
                         dataset="hsapiens_regulatory_feature",
                         GRCh = 37)



# promoter_grl = get_region_gr("Promoter")
# ctcf_grl = get_region_gr("CTCF Binding Site")
region_grl = get_regions_grl()


#Determine file names
hmf_overview_all = suppressWarnings(read_tsv("hmf_metadata.tsv"))
hmf_overview = dplyr::filter(hmf_overview_all, primaryTumorLocation == "Skin" & cancerSubtype == "Melanoma")
#hmf_overview = hmf_overview[1:5,]
hmf_dir = "~/hpc/pmc_vanboxtel/_old_structure/projects/Axel_GenoEcoli/HMF_analysis/somatics"
vcf_fnames = file.path(hmf_dir, hmf_overview$setName, paste0(hmf_overview$sampleId, ".purple.somatic.vcf"))

#Set sample names
sample_names = hmf_overview$sampleId

#purrr::map2(vcf_fnames, sample_names, create_matrices_sample, ref_genome, region_grl)


matrices_l_l = purrr::map(sample_names, ~readRDS(paste0("mutation_matrices/", "mut_matrices_", .x, ".rds")))
mut_mat = purrr::map(matrices_l_l, 1) %>% 
    do.call(cbind, .)
mut_mat_extended = purrr::map(matrices_l_l, 2) %>% 
    do.call(cbind, .)
mut_mat_extended2 = purrr::map(matrices_l_l, 3) %>% 
    do.call(cbind, .)
# mut_mat_region_promoter = purrr::map(matrices_l_l, 4) %>% 
#     do.call(cbind, .)
# mut_mat_long_promoter = purrr::map(matrices_l_l, 5) %>% 
#     do.call(cbind, .)
# mut_mat_region_ctcf = purrr::map(matrices_l_l, 6) %>% 
#     do.call(cbind, .)
# mut_mat_long_ctcf = purrr::map(matrices_l_l, 7) %>% 
#     do.call(cbind, .)
mut_mat_region = purrr::map(matrices_l_l, 4) %>%
        do.call(cbind, .)
mut_mat_region_extended = purrr::map(matrices_l_l, 5) %>% 
        do.call(cbind, .)
indel_mut_mat = purrr::map(matrices_l_l, 6) %>%
    do.call(cbind, .)

# Look at extended contexts
ext_heat_fig = plot_profile_heatmap(mut_mat_extended, max = 0.1)
ggsave("extended_profile_heat.pdf", ext_heat_fig)
ext_heat_fig2 = plot_profile_heatmap(mut_mat_extended2, max = 0.025)
ggsave("extended_profile_heat2.pdf", ext_heat_fig2)

# Pool the regions across samples and compare them
mut_mat_region_pool = pool_mut_mat(mut_mat_region, grouping = str_remove(colnames(mut_mat_region), ".*\\."))
region_profile_fig = plot_96_profile(mut_mat_region_pool, ymax = 0.3, condensed = TRUE)
ggsave("region_profile.pdf", region_profile_fig)

# Test differences between regions
chisq.test(mut_mat_region_pool, simulate.p.value = TRUE)
cos_sim_m = cos_sim_matrix(mut_mat_region_pool, mut_mat_region_pool)
cos_heat_fig = plot_cosine_heatmap(cos_sim_m, cluster_rows = TRUE, cluster_cols = TRUE)
ggsave("cos_region_heat.pdf", cos_heat_fig)

# Pool the extended regions across samples and compare them
mut_mat_region_pool = pool_mut_mat(mut_mat_region_extended, grouping = str_remove(colnames(mut_mat_region_extended), ".*\\."))
ext_region_heat_fig = plot_profile_heatmap(mut_mat_region_pool[,c(4,5)], max = 0.1, by = colnames(mut_mat_region_pool[,c(4,5)]))
ggsave("region_extended_profile_heat.pdf")

# Create signatures with regional features.
mut_mat_long = lengthen_mut_matrix(mut_mat_region)
sc <- scNMFSet(count = mut_mat_long)
set.seed(4)
estimate_bayes <- vb_factorize(sc, ranks = 9:14, nrun = 20, 
                               progress.bar = FALSE, verbose = 0)
plot(estimate_bayes)

rank = 13
nmf_res = extract_signatures(mut_mat_long, rank = rank, nmf_type = "variational_bayes", nrun = 50)
nmf_res$signatures = prop.table(nmf_res$signatures, 2)
saveRDS(nmf_res, "regional_nmf_res.rds")

ctcf_sig = nmf_res$signatures[1:96,]
enhancer_sig = nmf_res$signatures[97:192,]
open_chrom_sig = nmf_res$signatures[193:288,]
promoter_sig = nmf_res$signatures[289:384,]
promoter_flanking_sig = nmf_res$signatures[385:480,]
tf_binding_sig = nmf_res$signatures[481:576,]
other_sig = nmf_res$signatures[577:672,]

plot_96_profile(ctcf_sig, ymax = 0.5, condensed = TRUE)
plot_96_profile(enhancer_sig, ymax = 0.5, condensed = TRUE)
plot_96_profile(promoter_sig, ymax = 0.5, condensed = TRUE)
plot_96_profile(promoter_flanking_sig, ymax = 0.5, condensed = TRUE)
plot_96_profile(other_sig, ymax = 0.5, condensed = TRUE)


sig_l = list(ctcf_sig, enhancer_sig, open_chrom_sig, promoter_sig, promoter_flanking_sig, tf_binding_sig, other_sig)
names(sig_l) = str_remove(colnames(mut_mat_region), ".*\\.")[1:7]
fig_l_l = plot_sims_region_sigs(sig_l)

fig_l_l[2]
fig_l_l[12]
ggsave("region_signature_cosheat.pdf", fig_l_l[[12]][[1]])
ggsave("region_signature_profile.pdf", fig_l_l[[12]][[2]])

region_weights = purrr::map(fig_l_l, 3) %>% 
    bind_rows()
write_tsv(region_weights, "region_signatures_weights.txt")
## THERE IS A FACTOR 5 DIFFERENCE between the sigs in their promoter weight.
max(region_weights$Promoter) / min(region_weights$Promoter)


region_weights_long = region_weights %>% 
    tibble::rownames_to_column("signature") %>% 
    tidyr::pivot_longer(-signature, names_to = "Region", values_to = "Contribution")
ggplot(region_weights_long, aes(x = Region, y = Contribution)) +
    geom_point()


plot_sims_region_sigs = function(sig_l){
    rank = ncol(sig_l[[1]])
    fig_l_l = purrr::map(seq_len(rank), plot_sims_region_singlesig, sig_l)
    return(fig_l_l)
}

plot_sims_region_singlesig = function(i, sig_l){
    m = purrr::map(sig_l, function(x){return(x[,i, drop = FALSE])}) %>% 
        do.call(cbind, .)
    colnames(m) =  names(sig_l)
    rownames(m) = rownames(mut_mat)
    cos_sim_m = cos_sim_matrix(m, m)
    heat_fig = plot_cosine_heatmap(cos_sim_m, cluster_rows = FALSE, plot_values = TRUE)
    prof_fig = plot_96_profile(m)
    region_weights = colSums(m)
    return(list(heat_fig, prof_fig, region_weights))
}
