

analyze_cost_deconstructsigs = function(signatures, out_dir, feature_names){
    # Got to output directory
    old_dir = setwd(out_dir)
    on.exit(setwd(old_dir), add = TRUE)
    
    # Read simulated data
    res_l = readRDS("simulated_data/res_l.rds")
    
    cost_tb = purrr::map(res_l, single_analyze_cost_deconstructsigs, signatures, feature_names) %>% 
        bind_rows(.id = "matrix_name")
    
    # Create output directory
    if (!dir.exists("costs")){
        dir.create("costs")
    }
    write_tsv(cost_tb, "costs/cost_deconstructsigs.txt")
    invisible(0)
    
}

single_analyze_cost_deconstructsigs = function(res, signatures, feature_names){
    cutoffs = seq(0.00, 0.18, by = 0.02)
    cost_tb = purrr::map(cutoffs, single_strictness_analyze_cost_deconstructsigs, res, signatures, feature_names) %>% 
        bind_rows()
    message("Finished with one matrix.")
    return(cost_tb)
}

single_strictness_analyze_cost_deconstructsigs = function(cutoff, res, signatures, feature_names){
    
    decon_input_l = get_correct_format_deconstructsigs(res$mut_mat, signatures, feature_names)
    
    
    # Run whichSignatures for each sample.
    contri_rel = purrr::map(seq_len(nrow(decon_input_l$input)), ~whichSignatures(tumor.ref = decon_input_l$input, 
                                                                        signatures.ref = decon_input_l$signatures, 
                                                                        sample.id = .x,
                                                                        signature.cutoff = cutoff,
                                                                        contexts.needed = TRUE)) %>% 
        purrr::map("weights") %>% 
        bind_rows()
    contri = t(contri_rel) %*% diag(rowSums(decon_input_l$input))

    # Calculate cost
    cost_tb = cost_contribution(res$contribution, contri)
    cost_tb$cutoff = cutoff
    cost_tb$method = "deconstructSigs"
    message("Finished with one simulation.")
    return(cost_tb)
}


get_correct_format_deconstructsigs = function(mut_mat, signatures, feature_names){
    # Get correct format for deconstructsigs
    deconstruct_input = mut_mat %>% 
        t() %>% 
        as.data.frame()
    colnames(deconstruct_input) = feature_names
    rownames(deconstruct_input) = seq_len(nrow(deconstruct_input))
    
    # Get signatures in correct format
    deconstruct_signatures = signatures %>% 
        t() %>%
        as.data.frame()
    colnames(deconstruct_signatures) = feature_names
    
    return(list("input" = deconstruct_input, "signatures" = deconstruct_signatures))
}
