library(devtools)
load_all("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/Freek_MutationalPatterns/MutationalPatterns/")
library(tidyverse)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(biomaRt)
library("ccfindR")

out_dir = "~/surfdrive/Shared/Projects/Freek/mutpatterns/AHH1/"
if (!dir.exists(out_dir)){
    dir.create(out_dir)
}
setwd(out_dir)

#Determine file names
root = "~/surfdrive/Shared/Projects/Markus/AHH1/hg19/"
fnames = c("HPRT1/filter_HPRT1_C_ABSENT_SC_CLONAL_SNV_AHH1HPRT1CG6SCB8.vcf",
  "MSH23/filter_MSH23_C_ABSENT_SC_CLONAL_SNV_AHH1MSH23H5SCE3.vcf",
  "UNG4/filter_UNG4_C_ABSENT_SC_CLONAL_SNV_AHH1UNG4F2SCG9.vcf",
  "WT/filter_WT_C_ABSENT_SC_CLONAL_SNV_AHH1WTG2SCG6.vcf",
  "XPC3/filter_XPC3_C_ABSENT_SC_CLONAL_SNV_AHH1XPC3G2SC1G4.vcf")
vcf_files = paste0(root, fnames)

fnames_indel = c("HPRT1/filter_HPRT1_C_ABSENT_SC_CLONAL_INDEL_AHH1HPRT1CG6SCB8.vcf",
                 "MSH23/filter_MSH23_C_ABSENT_SC_CLONAL_INDEL_AHH1MSH23H5SCE3.vcf",
                 "UNG4/filter_UNG4_C_ABSENT_SC_CLONAL_INDEL_AHH1UNG4F2SCG9.vcf",
                 "WT/filter_WT_C_ABSENT_SC_CLONAL_INDEL_AHH1WTG2SCG6.vcf",
                 "XPC3/filter_XPC3_C_ABSENT_SC_CLONAL_INDEL_AHH1XPC3G2SC1G4.vcf")

vcf_files_indel = paste0(root, fnames_indel)

#Set sample names
sample_names = c("HPRT", "MSH", "UNG", "WT", "XPC")

#Read data
grl = read_vcfs_as_granges(vcf_files, sample_names, ref_genome, type = "all")
snv_grl = get_mut_type(grl, "snv")
dbs_grl = get_mut_type(grl, "dbs")
mbs_grl = get_mut_type(grl, "mbs")

indel_grl = read_vcfs_as_granges(vcf_files_indel, sample_names, ref_genome, type = "indel")

#Create spectra
type_occurrences = mut_type_occurrences(snv_grl, ref_genome)
spec_all_fig = plot_spectrum(type_occurrences, indv_points = TRUE)
ggsave("spec_all.pdf", spec_all_fig)
spectra_fig = plot_spectrum(type_occurrences, by = sample_names, error_bars = "none")
ggsave("spectra.pdf", spectra_fig)

#Create profile
mut_mat = mut_matrix(snv_grl, ref_genome)
write.table(mut_mat, "mut_mat.txt", sep = "\t", 
            quote = FALSE, row.names = TRUE)
profile_fig = plot_96_profile(mut_mat, ymax = 0.1, condensed = TRUE)
ggsave("profile.pdf", profile_fig)

#Create profile comparison of old vs new
mut_mat_dirty = mut_matrix(grl, ref_genome)
profile_comp_fig = plot_compare_profiles(mut_mat_dirty[,4],
                                         mut_mat[,4], 
                                         profile_names = c("With DBS/MBS", "Only SNV"), 
                                         condensed = TRUE)
ggsave("profile_comp.pdf", profile_comp_fig)

#Look at extended context
mut_mat_ext_context = mut_matrix(snv_grl, ref_genome, extension = 2)
profile_ext_fig = plot_profile_heatmap(mut_mat_ext_context, 
                                       by = sample_names,
                                       condensed = TRUE)
ggsave("profile_ext.pdf", profile_ext_fig)

river_fig <- plot_river(mut_mat_ext_context, condensed = TRUE)
ggsave("river.pdf", river_fig)

#Look at indel profile
indel_grl = get_indel_context(indel_grl, ref_genome)
indel_counts = count_indel_contexts(indel_grl)
indel_counts = indel_counts[,c("WT", "MSH")]
write.table(indel_counts, "indel_counts.txt", sep = "\t", 
            quote = FALSE, row.names = TRUE)
indel_context_fig = plot_indel_contexts(indel_counts, condensed = TRUE, same_y = TRUE)
ggsave("indel_profile.pdf", indel_context_fig)
indel_main_context_fig = plot_main_indel_contexts(indel_counts, same_y = TRUE)
ggsave("indel_main_profile.pdf", indel_main_context_fig)

#Look at DBSs
dbs_grl = get_dbs_context(dbs_grl)
dbs_counts = count_dbs_contexts(dbs_grl)
dbs_counts = dbs_counts[,c("WT", "XPC")]
dbs_context_fig = plot_dbs_contexts(dbs_counts, same_y = TRUE, condensed = TRUE)
ggsave("dbs_profile.pdf", dbs_context_fig)
dbs_main_context_fig = plot_main_dbs_contexts(dbs_counts, same_y = TRUE)
ggsave("dbs_main_context.pdf", dbs_main_context_fig)


#Look at mbs
#mbs_counts = count_mbs_contexts(mbs_grl)
#plot_mbs_contexts(mbs_counts)


#Look at snv signatures
signatures = get_known_signatures(source = "COSMIC_v3.1")
fit_res = fit_to_signatures(mut_mat, signatures)
sig_contri_fig = plot_contribution(fit_res$contribution, coord_flip = FALSE, mode = "absolute")
ggsave("sig_contri.pdf", sig_contri_fig)
reconstruction_fig = plot_original_vs_reconstructed(mut_mat, fit_res$reconstructed)
ggsave("reconstruction.pdf", reconstruction_fig)

strict_refit = fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.015)
fit_res_strict = strict_refit$fit_res
sig_contri_strict_fig = plot_contribution(fit_res_strict$contribution, coord_flip = FALSE, mode = "absolute")
ggsave("sig_contri_strict.pdf", sig_contri_strict_fig)

decay_fig = strict_refit$sim_decay_fig[[2]]
ggsave("MSH2_decay.pdf", decay_fig)

reconstruction_strict_fig = plot_original_vs_reconstructed(mut_mat, fit_res_strict$reconstructed)
ggsave("reconstruction_strict.pdf", reconstruction_strict_fig)


contri_boots <- fit_to_signatures_bootstrapped(mut_mat,
                                               signatures,
                                               n_boots = 100,
                                               method = "strict",
                                               max_delta = 0.015
)
sig_contri_boots_fig = plot_bootstrapped_contribution(contri_boots)
ggsave("sig_contri_boots.pdf", sig_contri_boots_fig)

sig_contri_boots_fig2 = plot_bootstrapped_contribution(contri_boots, plot_type = "dotplot")
ggsave("sig_contri_boots2.pdf", sig_contri_boots_fig2)

fig_list = plot_correlation_bootstrap(contri_boots)
ggsave("UNG_bootstrap_correlation.pdf", fig_list[[3]])

#Look at indel signatures
indel_signatures = get_known_signatures("indel", source = "COSMIC_v3.1")
fit_res_indel = fit_to_signatures(indel_counts, indel_signatures)
indel_sig_contri_fig = plot_contribution(fit_res_indel$contribution, coord_flip = FALSE, mode = "absolute")
ggsave("indel_sig_contri.pdf", indel_sig_contri_fig)

indel_reconstruction_fig = plot_original_vs_reconstructed(indel_counts, fit_res_indel$reconstructed)
ggsave("indel_reconstruction.pdf", indel_reconstruction_fig)

#fit_res_indel$contribution = fit_res_indel$contribution[rowSums(fit_res_indel$contribution) >= 10,]
#plot_contribution(fit_res_indel$contribution, coord_flip = FALSE, mode = "absolute")

strict_refit_indel = fit_to_signatures_strict(indel_counts, indel_signatures, max_delta = 0.0001)
fit_res_indel_strict = strict_refit_indel$fit_res
indel_sig_contri_strict_fig = plot_contribution(fit_res_indel_strict$contribution, coord_flip = FALSE, mode = "absolute")
ggsave("indel_sig_contri_strict.pdf", indel_sig_contri_strict_fig)

indel_reconstruction_strict_fig = plot_original_vs_reconstructed(indel_counts, fit_res_indel_strict$reconstructed)
ggsave("indel_reconstruction_strict.pdf", indel_reconstruction_strict_fig)


contri_boots_indel <- fit_to_signatures_bootstrapped(indel_counts,
                                               indel_signatures,
                                               n_boots = 100,
                                               method = "strict",
                                               max_delta = 0.0001
)
indel_sig_contri_boots_fig = plot_bootstrapped_contribution(contri_boots_indel, mode = "relative")
ggsave("indel_sig_contri_boots.pdf", indel_sig_contri_boots_fig)

indel_sig_contri_boots_fig2 = plot_bootstrapped_contribution(contri_boots_indel, mode = "absolute",
                                                            plot_type = "dotplot")
ggsave("indel_sig_contri_boots2.pdf", indel_sig_contri_boots_fig2)

fig_list = plot_correlation_bootstrap(contri_boots_indel)
ggsave("indel_XPC_bootstrap_correlation.pdf", fig_list[[2]])


#Look at dbs signatures
dbs_signatures = get_known_signatures("dbs", source = "COSMIC_v3.1")
fit_res_dbs = fit_to_signatures(dbs_counts, dbs_signatures)
dbs_sig_contri_fig = plot_contribution(fit_res_dbs$contribution, coord_flip = FALSE, mode = "absolute")
ggsave("dbs_sig_contri.pdf", dbs_sig_contri_fig)

dbs_reconstruction_fig = plot_original_vs_reconstructed(dbs_counts, fit_res_dbs$reconstructed)
ggsave("dbs_reconstruction.pdf", dbs_reconstruction_fig)

#fit_res_dbs$contribution = fit_res_dbs$contribution[rowSums(fit_res_dbs$contribution) >= 2,]
#plot_contribution(fit_res_dbs$contribution, coord_flip = FALSE, mode = "absolute")

strict_refit_dbs = fit_to_signatures_strict(dbs_counts, dbs_signatures, max_delta = 0.02)
fit_res_dbs_strict = strict_refit_dbs$fit_res
dbs_sig_contri_strict_fig = plot_contribution(fit_res_dbs_strict$contribution, coord_flip = FALSE, mode = "absolute")
ggsave("dbs_sig_contri_strict.pdf", dbs_sig_contri_strict_fig)

dbs_reconstruction_strict_fig = plot_original_vs_reconstructed(dbs_counts, fit_res_dbs_strict$reconstructed)
ggsave("dbs_reconstruction_strict.pdf", dbs_reconstruction_srict_fig)

contri_boots_dbs <- fit_to_signatures_bootstrapped(dbs_counts,
                                                     dbs_signatures,
                                                     n_boots = 100,
                                                     method = "strict",
                                                     max_delta = 0.02
)
dbs_sig_contri_boots_fig = plot_bootstrapped_contribution(contri_boots_dbs, mode = "relative")
ggsave("dbs_sig_contri_boots.pdf", dbs_sig_contri_boots_fig)

dbs_sig_contri_boots_fig2 = plot_bootstrapped_contribution(contri_boots_dbs, mode = "absolute",
                                                          plot_type = "dotplot")
ggsave("dbs_sig_contri_boots2.pdf", dbs_sig_contri_boots_fig2)

# Region specific analyses

do_chisquare = function(type_occurrences){
  m = type_occurrences %>% 
    dplyr::select(-Sample, -Region, -`C>T`) %>% 
    as.matrix()
  
  res = chisq.test(m, simulate.p.value = TRUE)
  res_tb = broom::tidy(res)
  return(res_tb)
}

create_region_profiles = function(grl, region_grl, name, include_other = TRUE){
  
  # Set out dir
  out_dir = file.path("region_analyses", name)
  if (!dir.exists(out_dir)){
    dir.create(out_dir, recursive = TRUE)
  }

  old_dir = setwd(out_dir)
  on.exit(setwd(old_dir), add = TRUE)
  
  # Split mutations
  grl_region = split_muts_region(grl, region_grl, include_other = include_other)
  
  # Create spectra
  type_occurrences_region <- mut_type_occurrences(grl_region, ref_genome)
  type_occrurences_l = type_occurrences_region %>% 
    tibble::rownames_to_column("Sample") %>% 
    write_tsv(paste0(name, "_type_occurrences.txt")) %>% 
    tidyr::separate("Sample", c("Sample", "Region"), sep = "\\.") %>% 
    split(., .$Sample)
  res_tb = purrr::map(type_occrurences_l, do_chisquare) %>% 
    dplyr::bind_rows() %>% 
    dplyr::mutate(fdr = p.adjust(p.value, method = "fdr"),
                  method = stringr::str_remove_all(method, "\\t|\\n")) %>% 
    write_tsv(paste0(name, "_chisq_type_occur.txt"))
  
  
  region_spectra_fig = plot_spectrum_region(type_occurrences_region, 
                                            by = sample_names,
                                            error_bars = "none")
  ggsave(paste0(name, "_region_spectra.pdf"), region_spectra_fig)
  region_spectra2_fig = plot_spectrum_region(type_occurrences_region, 
                                             by = sample_names, 
                                             mode = "relative_sample",
                                             error_bars = "none")
  ggsave(paste0(name, "_region_spectra2.pdf"), region_spectra2_fig)
  
  
  #Create profile
  mut_mat_region <- mut_matrix(grl_region, ref_genome)
  mut_mat_long <- lengthen_mut_matrix(mut_mat_region)
  mut_mat_long %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("Sample") %>% 
    write_tsv(paste0(name, "_mut_mat.txt"))
  
  region_profile_fig = plot_profile_region(mut_mat_long, ymax = 0.1)
  ggsave(paste0(name, "_region_profile.pdf"), region_profile_fig)
  
  invisible(0)
}


subsample_grl = function(grl){
  nr_muts = min(elementNROWS(grl))
  subsampled_grl = purrr::map(as.list(grl), function(gr) gr[sample.int(length(gr), nr_muts)]) %>% 
    GRangesList()
  return(subsampled_grl)
}


regulatory <- useEnsembl(biomart="regulation", 
                         dataset="hsapiens_regulatory_feature",
                         GRCh = 37)
promoter <- getBM(
    attributes = c(
        "chromosome_name", "chromosome_start",
        "chromosome_end", "feature_type_name"
    ),
    filters = "regulatory_feature_type_name",
    values = "Promoter",
    mart = regulatory
)
promoter_g <- reduce(GRanges(
    promoter$chromosome_name,
    IRanges(
        promoter$chromosome_start,
        promoter$chromosome_end
    )
))
seqlevels(promoter_g) <- c(1:22, "X", "Y")
promoter_g <- sort(promoter_g)
seqlevelsStyle(promoter_g) = "UCSC"
promoter_g = GRangesList("promoter" = promoter_g)
create_region_profiles(snv_grl, promoter_g, "promoter")

exons_gr = exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
seqlevels(exons_gr, pruning.mode = "coarse") <- paste0("chr" ,c(1:22, "X", "Y"))
exons_gr = GRangesList("exons" = exons_gr)
create_region_profiles(snv_grl, exons_gr, "exons")

genes_gr = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
seqlevels(genes_gr, pruning.mode = "coarse") <- paste0("chr" ,c(1:22, "X", "Y"))
genes_gr = GRangesList("genes" = genes_gr)
create_region_profiles(snv_grl, genes_gr, "genes")

repli_fname = "~/surfdrive/Shared/Boxtel_General/Data/Replication_timing/ENCODE/blood/B_lymphocytes.bed"
repli_bed = read.table(repli_fname, header = F, sep = "\t", stringsAsFactors = F)
colnames(repli_bed) = c("chrom", "start", "end", "repli_time")
repli_gr = makeGRangesFromDataFrame(repli_bed, keep.extra.columns = T, starts.in.df.are.0based = T)
repli_gr$repli_time_cat = dplyr::case_when(
  repli_gr$repli_time >= 60 ~ "Early",
  repli_gr$repli_time <= 33 ~ "Late",
  TRUE ~ "Intermediate"
)
repli_grl = split(repli_gr, repli_gr$repli_time_cat)
create_region_profiles(snv_grl, repli_grl, "replication_time", include_other = FALSE)


# Subsampled analyses.
set.seed(42)
snv_grl_subsampled = subsample_grl(snv_grl)
create_region_profiles(snv_grl_subsampled, promoter_g, "promoter_subsampled")
create_region_profiles(snv_grl_subsampled, exons_gr, "exons_subsampled")
create_region_profiles(snv_grl_subsampled, genes_gr, "genes_subsampled")
create_region_profiles(snv_grl_subsampled, repli_grl, "replication_time_subsampled", include_other = FALSE)

# Do lesion segregation
lesion_tb = calculate_lesion_segregation(snv_grl, sample_names)
write_tsv(lesion_tb, "lesion_segregation.txt")
lesion_fig = plot_lesion_segregation(snv_grl[[2]], per_chrom = TRUE)
ggsave("lesion_segregation.pdf", lesion_fig[[1]])

lesion_fig2 = plot_lesion_segregation(snv_grl)
ggsave("lesion_segregation2.pdf", lesion_fig2)

# Do bayes NMF estimate
sc <- scNMFSet(count = mut_mat)
set.seed(4)
estimate_bayes <- vb_factorize(sc, ranks = 1:3, nrun = 1, 
                               progress.bar = FALSE, verbose = 0)
pdf("nmf_bayes_estimate.pdf")
plot(estimate_bayes)
dev.off()
