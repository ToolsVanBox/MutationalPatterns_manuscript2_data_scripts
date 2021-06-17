library(ggbeeswarm)
plot_simulation_extraction = function(out_dir){
    
    # Got to output/input directory
    old_dir = setwd(out_dir)
    on.exit(setwd(old_dir), add = TRUE)
    
    # Read cost files
    cost_mutpatterns_tb = read_tsv("costs/cost_mutpatterns.txt")
    cost_somaticsignatures = read_tsv("costs/cost_somaticsignatures.txt")
    cost_signer = read_tsv("costs/cost_signer.txt")
    cost_sigprofiler = read_tsv("costs/cost_sigprofiler.txt")
    cost_tb = rbind(cost_mutpatterns_tb, 
                    cost_signer,
                    cost_somaticsignatures, 
                    cost_sigprofiler)
    
    # Correct format
    cost_tb = cost_tb %>% 
        dplyr::mutate(matrix_name = factor(matrix_name, levels = unique(matrix_name)),
                      signature = factor(signature, levels = unique(signature)),
                      method = factor(method, levels = unique(method)))
    
    #Create figures
    cossim_extract_fig = ggplot(cost_tb, aes(x = matrix_name, y = cosine, colour = method)) +
        geom_boxplot(outlier.shape = NA) +
        geom_quasirandom(dodge.width = 0.8, size = 1) +
        theme_classic() +
        theme(text = element_text(size = 20))
    ggsave("cossim_extract.pdf", cossim_extract_fig, width = 10)
    
    cost_tb_bayes = dplyr::filter(cost_tb, method == "MutationalPatterns_variational_bayes")
    per_sig_fig = ggplot(cost_tb_bayes, aes(x = signature, y = cosine, shape = method)) +
        facet_grid(matrix_name ~ .) +
        geom_bar(stat = "identity") +
        theme_classic() +
        theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90))
    ggsave("cossim_extract_per_sig_variational_bayes.pdf", per_sig_fig)
    
    invisible(0)
}
