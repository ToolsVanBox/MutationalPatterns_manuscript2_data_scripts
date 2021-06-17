library(devtools)
library(tidyverse)
library(deconstructSigs)
library(decompTumor2Sig)
source("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/MutationalPatterns_manuscript2_data_scripts/scripts/simulation_deconstructsigs_functions.R")
source("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/MutationalPatterns_manuscript2_data_scripts/scripts/simulation_decomptumor2sig_functions.R")
setwd("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/Freek_MutationalPatterns/MutationalPatterns/")
load_all()


out_dir = "~/surfdrive/Shared/Projects/Freek/mutpatterns/refit_comparison/"
if (!dir.exists(out_dir)){
    dir.create(out_dir)
}
setwd(out_dir)


mut_mat = as.matrix(read.table("../AHH1/mut_mat.txt", sep = "\t", header = TRUE))

# Get signatures
signatures = get_known_signatures(source = "COSMIC_v3.1")

# Get context names
feature_names <- rownames(readRDS(system.file("states/mut_mat_data.rds",
                                              package = "MutationalPatterns"
)))

# Run deconstructsigs (this is the best cutoff when using 100 muts per sig)
decon_input_l = get_correct_format_deconstructsigs(mut_mat, signatures, feature_names)
contri_rel = purrr::map(seq_len(nrow(decon_input_l$input)), ~whichSignatures(tumor.ref = decon_input_l$input, 
                                                                             signatures.ref = decon_input_l$signatures, 
                                                                             sample.id = .x,
                                                                             signature.cutoff = 0.1,
                                                                             contexts.needed = TRUE)) %>% 
    purrr::map("weights") %>% 
    bind_rows()
contri = t(contri_rel) %*% diag(rowSums(decon_input_l$input))

colnames(contri) = colnames(mut_mat)
write.table(contri, "deconstructsigs_contri.txt", quote = F, sep = "\t")
decon_contri_fig = plot_contribution(contri, mode = "absolute")
ggsave("deconstructsigs_contri.pdf", decon_contri_fig)


# Run decomptumor2sigs (this is the best cutoff when using 100 muts per sig)
signatures_renorm = prop.table(signatures, 2) # The package won't work without this renomalization.
decomp_input_l = get_correct_format_decomptumor2sig(mut_mat, signatures_renorm, feature_names)
contri_rel = decomposeTumorGenomes(decomp_input_l$decomp_input, 
                                   decomp_input_l$signatures_l, 
                                   minExplainedVariance = 0.85,
                                   greedySearch = TRUE)

fail_f = elementNROWS(contri_rel) == 0
if (sum(fail_f)){
    contri_res_failed = decomposeTumorGenomes(decomp_input_l$decomp_input[fail_f], 
                                              decomp_input_l$signatures_l)
    contri_rel = do.call(cbind, c(contri_rel, contri_res_failed))
} else{
    contri_rel = do.call(cbind, contri_rel)
    
}

# Format contribution to correct output
contri_rel = contri_rel[,names(decomp_input_l$decomp_input)]
contri_rel[is.na(contri_rel)] = 0
contri = contri_rel %*% diag(colSums(mut_mat))

colnames(contri) = colnames(mut_mat)
write.table(contri, "decomptumor2sig_contri.txt", quote = F, sep = "\t")
decomp_contri_fig = plot_contribution(contri, mode = "absolute")
ggsave("decomptumor2sig_contri.pdf", decomp_contri_fig)


mut_mat = as.data.frame(mut_mat)
mut_mat$MutationType = feature_names
mut_mat = mut_mat[,c(ncol(mut_mat), seq(1, ncol(mut_mat)-1))]

# Write mut_mat and contribution
write.table(mut_mat,
            "mut_mat_sigprofiler.txt",
            sep = "\t",
            quote = F,
            row.names = F)

## Run sigprofiler on hpc.
contri = read_tsv("~/surfdrive/Shared/Projects/Freek/mutpatterns/refit_comparison/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities_refit.txt")
                         
# Transpose sigprofiler contribution and prepare for adding absent signatures
contri = contri %>%
    #dplyr::select(str_which(colnames(.), "SBS96[A-Z]", negate = TRUE)) %>%
    column_to_rownames("Samples") %>%
    t()




write.table(contri, "SigProfiler_contri.txt", quote = F, sep = "\t")
sigprofiler_contri_fig = plot_contribution(contri, mode = "absolute")
ggsave("SigProfiler_contri.pdf", sigprofiler_contri_fig)



