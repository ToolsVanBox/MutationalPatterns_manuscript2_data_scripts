do_simulation_extraction_signer = function(signatures, out_dir, feature_names){
    # Got to output directory
    if (!dir.exists(out_dir)){
        dir.create(out_dir)
    }
    old_dir = setwd(out_dir)
    on.exit(setwd(old_dir), add = TRUE)
    
    # Read simulated data
    res_l = readRDS(paste0("../../simulations/", out_dir, "/simulated_data/res_l.rds"))
    
    sim_tb = purrr::map(res_l, single_extract_cost_signer, signatures, feature_names) %>% 
        bind_rows(.id = "matrix_name")
    
    # sim_pca_tb = purrr::imap(res_l, single_extract_cost_somaticsignatures, signatures, "pca") %>% 
    #     bind_rows(.id = "matrix_name")
    # 
    # sim_tb = rbind(sim_nmf_tb, sim_pca_tb)
    
    # Create output directory
    if (!dir.exists("costs")){
        dir.create("costs")
    }
    write_tsv(sim_tb, "costs/cost_signer.txt")
    invisible(0)
}


single_extract_cost_signer = function(res, signatures, feature_names){
    
    # Extract signatures
    mut_mat = res$mut_mat
    colnames(mut_mat) = as.character(seq_len(ncol(mut_mat)))
    rownames(mut_mat) = feature_names
    nmf_res = signeR(mut_mat, nsig = ncol(signatures), samples = "columns")
    
    # Check how well each signature has been extracted.
    # For each input signature the best cosine similarity with an extracted signature is calculated.
    cos_sim_m = cos_sim_matrix(nmf_res$Phat, signatures)
    best_sims = apply(cos_sim_m, 2, max)
    sim_tb = enframe(best_sims, name = "signature", value = "cosine") %>% 
        dplyr::mutate(method = "signeR")
    
    message("Finished with one signeR run.")
    return(sim_tb)
}


do_simulation_extraction_somaticsignatures = function(signatures, out_dir){
    # Got to output directory
    if (!dir.exists(out_dir)){
        dir.create(out_dir)
    }
    old_dir = setwd(out_dir)
    on.exit(setwd(old_dir), add = TRUE)
    
    # Read simulated data
    res_l = readRDS(paste0("../../simulations/", out_dir, "/simulated_data/res_l.rds"))
    
    sim_nmf_tb = purrr::map(res_l, single_extract_cost_somaticsignatures, signatures, "nmf") %>% 
        bind_rows(.id = "matrix_name")
    
    sim_pca_tb = purrr::map(res_l, single_extract_cost_somaticsignatures, signatures, "pca") %>% 
        bind_rows(.id = "matrix_name")
    
    sim_tb = rbind(sim_nmf_tb, sim_pca_tb)
    
    # Create output directory
    if (!dir.exists("costs")){
        dir.create("costs")
    }
    write_tsv(sim_tb, "costs/cost_somaticsignatures.txt")
    invisible(0)
}

single_extract_cost_somaticsignatures = function(res, signatures, extract_type){
    
    # Extract signatures
    if (extract_type == "nmf"){
        nmf_res = identifySignatures(res$mut_mat, ncol(signatures), nmfDecomposition, nrun = 200)
    } else if (extract_type == "pca"){
        nmf_res = identifySignatures(res$mut_mat, ncol(signatures), pcaDecomposition)
    } else{
        stop("extract_type should be either nmf or pca.")
    }
    
    # Check how well each signature has been extracted.
    # For each input signature the best cosine similarity with an extracted signature is calculated.
    cos_sim_m = cos_sim_matrix(signatures(nmf_res), signatures)
    best_sims = apply(cos_sim_m, 2, max)
    sim_tb = enframe(best_sims, name = "signature", value = "cosine") %>% 
        dplyr::mutate(method = paste0("SomaticSignatures_", extract_type))
    
    message("Finished with one SomaticSignatures run.")
    return(sim_tb)
}


do_simulation_extraction_mutpatterns = function(signatures, out_dir){
    
    # Got to output directory
    if (!dir.exists(out_dir)){
        dir.create(out_dir)
    }
    old_dir = setwd(out_dir)
    on.exit(setwd(old_dir), add = TRUE)
    
    # Read simulated data
    res_l = readRDS(paste0("../../simulations/", out_dir, "/simulated_data/res_l.rds"))
    
    sim_nmf_tb = purrr::map(res_l, single_extract_cost_mutationalpatterns, signatures, "regular", fudge = NULL) %>% 
        bind_rows(.id = "matrix_name")
    
    sim_bayes_tb = purrr::map(res_l, single_extract_cost_mutationalpatterns, signatures, "variational_bayes", fudge = 0.005) %>% 
        bind_rows(.id = "matrix_name")
    
    sim_tb = rbind(sim_nmf_tb, sim_bayes_tb)
    
    # Create output directory
    if (!dir.exists("costs")){
        dir.create("costs")
    }
    write_tsv(sim_tb, "costs/cost_mutpatterns.txt")
    invisible(0)
}

single_extract_cost_mutationalpatterns = function(res, signatures, nmf_type, fudge){
    
    # Extract signatures
    nmf_res = extract_signatures(res$mut_mat, rank = ncol(signatures), nmf_type = nmf_type, fudge = fudge)
    
    # Check how well each signature has been extracted.
    # For each input signature the best cosine similarity with an extracted signature is calculated.
    cos_sim_m = cos_sim_matrix(nmf_res$signatures, signatures)
    best_sims = apply(cos_sim_m, 2, max)
    sim_tb = enframe(best_sims, name = "signature", value = "cosine") %>% 
        dplyr::mutate(method = paste0("MutationalPatterns_", nmf_type))
    
    message("Finished with one MutationalPatterns run.")
    return(sim_tb)
}



do_simulation_extraction_sigprofiler = function(signatures, out_dir){
    # Got to output directory
    if (!dir.exists(out_dir)){
        dir.create(out_dir)
    }
    old_dir = setwd(out_dir)
    on.exit(setwd(old_dir), add = TRUE)
    
    dir_l = list.files("sigprofiler_output")
    matrix_names = str_replace(dir_l, "sigprofiler_output_(.*)_nnls_remove_.*", "\\1")
    names(dir_l) = matrix_names
    
    sim_tb = purrr::map(dir_l, single_extract_cost_sigprofiler, signatures) %>% 
        bind_rows(.id = "matrix_name")
    
    # Create output directory
    if (!dir.exists("costs")){
        dir.create("costs")
    }
    write_tsv(sim_tb, "costs/cost_sigprofiler.txt")
    invisible(0)
}

single_extract_cost_sigprofiler = function(dir, signatures){
    fname = paste0("sigprofiler_output/", dir, "/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt")
    if (!file.exists(fname)){
        fname = paste0("sigprofiler_output/", dir, "/ID83/Suggested_Solution/ID83_De-Novo_Solution/Signatures/ID83_De-Novo_Signatures.txt")
    }
    
    nmf_sigs = read_tsv(fname) %>% 
        dplyr::select(-1) %>% 
        as.matrix()
    cos_sim_m = cos_sim_matrix(nmf_sigs, signatures)
    best_sims = apply(cos_sim_m, 2, max)
    sim_tb = enframe(best_sims, name = "signature", value = "cosine") %>% 
        dplyr::mutate(method = "SigProfiler")
    return(sim_tb)
}
