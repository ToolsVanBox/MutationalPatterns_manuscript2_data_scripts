library(devtools)
library(tidyverse)
library(SomaticSignatures)
library(signeR)
setwd("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/Freek_MutationalPatterns/MutationalPatterns/")
load_all()


calc_avg_cossim_pair = function(pair){
    sim_tb = cos_sim_matrix(pair[[1]], pair[[2]])
    rowmax = sort(apply(sim_tb, 1, max))
    colmax = sort(apply(sim_tb, 2, max))
    if (!identical(rowmax, colmax)){
        stop("The row and column max similarities should be the same.")
    }
    avg_sim = mean(rowmax)
    avg_sim_tb = tibble("method_1" = names(pair), "method_2" = rev(names(pair)), "sim" = rep(avg_sim, 2))
    return(avg_sim_tb)
}

rename_nmf_signatures <- function(nmf_signatures, 
                                  signatures, 
                                  cutoff = 0.85, 
                                  base_name = "SBS", 
                                  suffix = "-like") {
    colnames(nmf_signatures) <- seq_len(ncol(nmf_signatures))
    sim_matrix <- cos_sim_matrix(signatures, nmf_signatures)
    
    # The number of signatures with no similarity based match.
    j <- 0
    
    # Iterate over signatures in nmf_res
    for (i in seq_len(ncol(sim_matrix))) {
        sig_sim <- sim_matrix[, i]
        
        # Check if there is a highly similar signature
        cossim <- max(sig_sim)
        if (cossim > cutoff) {
            
            # Determine most similar signature
            row <- which.max(sig_sim)
            sig <- names(sig_sim)[row]
            
            # Add suffix
            sig <- paste0(sig, suffix)
            
            # Change name of nmf_res signature
            colnames(nmf_signatures)[i] <- sig
        } else {
            # If there is no similar signature, use letters for the signature name.
            j <- j + 1
            colnames(nmf_signatures)[i] <- paste0(base_name, LETTERS[j])
        }
    }
    
    # Check if there are any signatures that have been given the same name.
    dupli_sigs <- colnames(nmf_signatures) %>%
        duplicated() %>%
        any()
    if (dupli_sigs) {
        stop("You have multiple NMF signatures that are linked to the same existing signature.\n
                Please use a lower rank in the NMF or increase the cutoff at which a NMF and \n
                existing signature are considered identical.", call. = FALSE)
    }
    
    # Return the nmf_res with the updated signature names.
    return(nmf_signatures)
}



plot_save_nmf_res = function(nmf_signatures, signatures, name){
    nmf_signatures = rename_nmf_signatures(nmf_signatures, signatures)
    p96_fig = plot_96_profile(nmf_signatures, condensed = TRUE)
    ggsave(paste0(name, "_96_profile.pdf"))
}

get_max_sims = function(nmf_signatures, signatures, name){
    sim_m = cos_sim_matrix(nmf_signatures, signatures)
    max_sims = apply(sim_m, 1, max)
    sim_tb = tibble("sim" = max_sims, "name" = name, "sig" = names(max_sims))
    return(sim_tb)
}

out_dir = "~/surfdrive/Shared/Projects/Freek/mutpatterns/extraction_comparison/"
if (!dir.exists(out_dir)){
    dir.create(out_dir)
}
setwd(out_dir)

# Get signatures
signatures = get_known_signatures(source = "COSMIC_v3.1")

# Get mutation matrixes
mut_mat_AHH1 = as.matrix(read.table("../AHH1/mut_mat.txt", sep = "\t", header = TRUE))
mut_mat_axel = as.matrix(read.table("~/surfdrive/Shared/Boxtel_General/Data/MutationMatrix/hg19/Axel_HealthyBonemarrow.txt", sep = "\t", header = TRUE))
mut_mat_blokzijl = as.matrix(read.table("~/surfdrive/Shared/Boxtel_General/Data/MutationMatrix/hg19/Blokzijl_2016.txt", sep = "\t", header = TRUE))
mut_mat_fetus = as.matrix(read.table("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/Mutation_accumulation_T21/data/Mutation_matrixes/all_muts_per_cell_divstri_mut_mat.txt", sep = "\t", header = TRUE))
mut_mat = cbind(mut_mat_AHH1, mut_mat_blokzijl)

# Test number of signatures to extract.
library("ccfindR")
sc <- scNMFSet(count = mut_mat)
set.seed(4)
estimate_bayes <- vb_factorize(sc, ranks = 1:10, nrun = 2, 
                               progress.bar = FALSE, verbose = 0)
plot(estimate_bayes)


#Mutationalpatterns bayes
nmf_res = extract_signatures(mut_mat, rank = 3, nmf_type = "variational_bayes")
plot_save_nmf_res(nmf_res$signatures, signatures, "MutationalPatterns_bayes")
signatures_mutpatterns_bayes = rename_nmf_signatures(nmf_res$signatures, signatures, cutoff = 0.75)
mutationalpatterns_bayes_sim_tb = get_max_sims(signatures_mutpatterns_bayes, signatures, "MutationalPatterns_bayes")

#MutationalPatterns nmf
nmf_res = extract_signatures(mut_mat, rank = 3, nmf_type = "regular")
plot_save_nmf_res(nmf_res$signatures, signatures, "MutationalPatterns_nmf")
signatures_mutpatterns_nmf = rename_nmf_signatures(nmf_res$signatures, signatures, cutoff = 0.75)
mutationalpatterns_nmf_sim_tb = get_max_sims(signatures_mutpatterns_nmf, signatures, "MutationalPatterns_nmf")

# Somatic signatures nmf
nmf_res = identifySignatures(mut_mat, 3, nmfDecomposition, nrun = 200)
plot_save_nmf_res(signatures(nmf_res), signatures, "SomaticSignatures_nmf")
signatures_somatic_signatures_nfm = rename_nmf_signatures(signatures(nmf_res), signatures, cutoff = 0.75)
somaticsignatures_nmf_sim_tb = get_max_sims(signatures_somatic_signatures_nfm, signatures, "SomaticSignatures_nmf")

# Somatic signatures pca
nmf_res = identifySignatures(mut_mat, 3, pcaDecomposition, nrun = 200)
plot_save_nmf_res(signatures(nmf_res), signatures, "SomaticSignatures_pca")
signatures_somatic_signatures_pca = rename_nmf_signatures(signatures(nmf_res), signatures, cutoff = 0.75)
somaticsignatures_pca_sim_tb = get_max_sims(signatures_somatic_signatures_pca, signatures, "SomaticSignatures_pca")

# signeR
nmf_res = signeR(mut_mat, nsig = 3, samples = "columns")
nmf_signatures = nmf_res$Phat
rownames(nmf_signatures) = rownames(mut_mat)
plot_save_nmf_res(nmf_signatures, signatures, "signeR")
signatures_signer = rename_nmf_signatures(nmf_signatures, signatures, cutoff = 0.75)
signer_sim_tb = get_max_sims(signatures_signer, signatures, "signeR")

# Sigprofiler
# Write mut mat for sigprofiler
# mut_mat_sigprofiler = as.data.frame(mut_mat)
# mut_mat_sigprofiler$MutationType = feature_names
# mut_mat_sigprofiler = mut_mat_sigprofiler[,c(ncol(mut_mat_sigprofiler), seq(1, ncol(mut_mat_sigprofiler)-1))]
# write.table(mut_mat_sigprofiler, "AHH1_blokzijl2016_mut_mat.txt", quote = F, sep = "\t")

# Run on hpc

# Get results
sigprofiler_res = read_tsv("sigprofiler_output/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt") %>% 
    column_to_rownames("MutationsType") %>% 
    as.matrix()
plot_save_nmf_res(sigprofiler_res, signatures, "SigProfiler")
signatures_sigprofiler = rename_nmf_signatures(sigprofiler_res, signatures, cutoff = 0.75)
sigprofiler_sim_tb = get_max_sims(signatures_sigprofiler, signatures, "SigProfiler")


# Determine how similar the extracted signatures are to the canonical signatures.
sim_tb = rbind(mutationalpatterns_bayes_sim_tb, 
               mutationalpatterns_nmf_sim_tb, 
               somaticsignatures_nmf_sim_tb, 
               signer_sim_tb, 
               sigprofiler_sim_tb) %>% 
    dplyr::mutate(sig = factor(sig, levels = c("SBS1-like", "SBS5-like", "SBS18-like")))
sim_to_canon_fig = ggplot(sim_tb, aes(x = sig, y = sim, fill = name)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    theme(text = element_text(size = 20))
ggsave("sim_to_canon.pdf", sim_to_canon_fig)


# Determine the average signature similarity between the different methods.
methods_l = list("MutationalPatterns_bayes" = signatures_mutpatterns_bayes,
     "MutationalPatterns_nmf" = signatures_mutpatterns_nmf,
     "somaticSignatures_nmf" = signatures_somatic_signatures_nfm,
     "signeR" = signatures_signer,
     "SigProfiler" = signatures_sigprofiler)

methods_pairs = combn(methods_l, 2, simplify = FALSE)
avg_sim_tb = purrr::map(methods_pairs, calc_avg_cossim_pair) %>% 
    bind_rows()
diag_sim_tb = tibble("method_1" = names(methods_l), "method_2" = names(methods_l), "sim" = 1)
avg_sim_tb = rbind(avg_sim_tb, diag_sim_tb) %>% 
    dplyr::mutate(method_1 = factor(method_1, levels = names(methods_l)),
                  method_2 = factor(method_2, levels = rev(names(methods_l))))

sim_heat_fig = ggplot(avg_sim_tb, aes(x = method_1, y = method_2, fill = sim)) +
    geom_raster() +
    geom_text(label = round(avg_sim_tb$sim, 4)) +
    scale_fill_distiller(palette = "RdYlBu", limits = c(0,1)) +
    labs(fill = "Mean cosine\nsimilarity", x = "", y = "") +
    theme_classic() +
    theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90))
ggsave("heatmap_comparison.pdf", sim_heat_fig)



