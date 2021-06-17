
analyze_cost_decomptumor2sigs = function(signatures, out_dir, feature_names){
    # Got to output directory
    old_dir = setwd(out_dir)
    on.exit(setwd(old_dir), add = TRUE)
    
    # Read simulated data
    res_l = readRDS("simulated_data/res_l.rds")
    
    cost_tb = purrr::map(res_l, single_analyze_cost_decomptumor2sigs, signatures, feature_names) %>% 
        bind_rows(.id = "matrix_name")
    
    # Create output directory
    if (!dir.exists("costs")){
        dir.create("costs")
    }
    write_tsv(cost_tb, "costs/cost_decomptumor2sigs.txt")
    invisible(0)
    
}

single_analyze_cost_decomptumor2sigs = function(res, signatures, feature_names){
    cutoffs = seq(0.65, 1, by = 0.025)
    cost_tb = purrr::map(cutoffs, single_strictness_analyze_cost_decomptumor2sigs, res, signatures, feature_names) %>% 
        bind_rows()
    message("Finished with one matrix.")
    return(cost_tb)
}

single_strictness_analyze_cost_decomptumor2sigs = function(cutoff, res, signatures, feature_names){
    
    
    decomp_input_l = get_correct_format_decomptumor2sig(res$mut_mat, signatures, feature_names)
    
    # Run decomposeTumorGenomes for each sample.
    contri_rel = decomposeTumorGenomes(decomp_input_l$decomp_input, 
                                      decomp_input_l$signatures_l, 
                                      minExplainedVariance = cutoff,
                                      greedySearch = TRUE)
    
    # For samples with too low variance explained even when using all signatures, 
    # run again with all signatures.
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
    contri = contri_rel %*% diag(colSums(res$mut_mat))

    # Calculate cost
    cost_tb = cost_contribution(res$contribution, contri)
    cost_tb$cutoff = cutoff
    cost_tb$method = "decomptumor2sigs"
    message("Finished with one simulation.")
    return(cost_tb)
}

# Function to get the mut matrix and signatures into the correct format.
get_correct_format_decomptumor2sig = function(mut_mat, signatures, feature_names){
    # Get correct format for decomptumor2sig
    decomp_input = mut_mat %>% 
        prop.table(2) %>% 
        as.data.frame() %>% 
        as.list()
    decomp_input = purrr::map(decomp_input, function(.x){
        names(.x) = feature_names
        return(.x)
    })
    names(decomp_input) = seq_len(ncol(mut_mat))
    
    # Get signatures in correct format
    signatures_l = signatures %>% 
        as.data.frame() %>% 
        as.list()
    signatures_l = purrr::map(signatures_l, function(.x){
        names(.x) = feature_names
        return(.x)
    })
    
    return(list("decomp_input" = decomp_input, "signatures_l" = signatures_l))
}
