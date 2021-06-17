# Calculate the costs of sigprofiler
analyze_cost_sigprofiler = function(out_dir){
    
    # Got to output directory
    old_dir = setwd(out_dir)
    on.exit(setwd(old_dir), add = TRUE)
    
    # Read simulated data
    res_l = readRDS("simulated_data/res_l.rds")
    
    cost_tb = purrr::imap(res_l, single_analyze_cost_sigprofiler) %>% 
        bind_rows(.id = "matrix_name")
    
    # Create output directory
    if (!dir.exists("costs")){
        dir.create("costs")
    }
    write_tsv(cost_tb, "costs/cost_sigprofiler.txt")
    invisible(0)
}

# Calculate the costs of a single matrix with sigprofiler
single_analyze_cost_sigprofiler = function(res, name, out_dir){
    cutoffs = c("0.1", "0.05", "0.01", "0.005", "0.001", "0.0005")
    cost_tb = purrr::map(cutoffs, single_strictness_analyze_cost_sigprofiler, res, name) %>% 
        bind_rows()
    return(cost_tb)
}

# Calculate the costs of a single matrix and a single strictness with sigprofiler.
single_strictness_analyze_cost_sigprofiler = function(cutoff, res, name){
    
    # Read SigProfiler results
    fname = paste0("sigprofiler_output/sigprofiler_output_",
                   name,
                   "_nnls_remove_",
                   cutoff,
                   "/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities_refit.txt")
    if (!file.exists(fname)){
        return(NULL)
    }
    contri = read_tsv(fname)
    
    
    # Transpose sigprofiler contribution and prepare for adding absent signatures
    contri = contri %>%
        #dplyr::select(str_which(colnames(.), "SBS96[A-Z]", negate = TRUE)) %>%
        column_to_rownames("Samples") %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("rowname")
    
    # prepare for adding absent signatures. (Created by NMF from sigprofiler.)
    contri_res = res$contribution %>%
        as.data.frame() %>%
        rownames_to_column("rowname")
    
    # Add absent signatures for both the real contribution and the estimated one.
    sig_ref <- tibble::tibble("rowname" = unique(c(contri_res$rowname, contri$rowname)))
    contri_all_sigs = dplyr::left_join(sig_ref, contri, by ="rowname") %>%
        column_to_rownames("rowname") %>%
        as.matrix()
    contri_all_sigs[is.na(contri_all_sigs)] <- 0
    
    contri_res_all_sigs = dplyr::left_join(sig_ref, contri_res, by ="rowname") %>%
        column_to_rownames("rowname") %>%
        as.matrix()
    contri_res_all_sigs[is.na(contri_res_all_sigs)] <- 0
    
    
    
    # Calculate cost
    cost_tb = cost_contribution(contri_res_all_sigs, contri_all_sigs)
    #cost_tb$max_delta = max_delta
    cost_tb$cutoff = cutoff
    cost_tb$method = "SigProfiler"
    
    message("Finished with one Sigprofiler run.")
    return(cost_tb)
    
}
