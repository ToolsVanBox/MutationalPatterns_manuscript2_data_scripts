create_write_simulated_data = function(nr_muts_l, signatures, nr_samples = 300, out_dir, feature_names){
    
    #Go to out dir
    if (!dir.exists(out_dir)){
        dir.create(out_dir)
    }
    old_dir = setwd(out_dir)
    on.exit(setwd(old_dir), add = TRUE)
    
    # Create lists of simulated matrixes
    res_l = purrr::map(nr_muts_l, create_simulated_profiles, nr_samples, signatures)
    write_simulated_data(res_l, feature_names)
    saveRDS(res_l, "simulated_data/res_l.rds")
    message("Created simulated matrixes")
    invisible(0)
}

write_simulated_data = function(res_l, feature_names){
    
    # Change working directory
    if (!dir.exists("simulated_data")){
        dir.create("simulated_data")
    }
    old_dir = setwd("simulated_data")
    on.exit(setwd(old_dir), add = F)
    
    # Write the simulated mut_mat and contribution files
    purrr::iwalk(res_l, write_single_simulated_data, feature_names)
    invisible(0)
}

write_single_simulated_data = function(res, name, feature_names){
    
    # Add column names
    mut_mat = as.data.frame(res$mut_mat)
    colnames(mut_mat) = as.character(seq_len(ncol(mut_mat)))
    
    # Add mutationType column
    mut_mat$MutationType = feature_names
    mut_mat = mut_mat[,c(ncol(mut_mat), seq(1, ncol(mut_mat)-1))]
    
    # Write mut_mat and contribution
    write.table(mut_mat,
                paste0("mut_mat_", name, ".txt"),
                sep = "\t",
                quote = F,
                row.names = F)
    write.table(res$contribution, paste0("contribution_", name, ".txt"),
                sep = "\t",
                quote = F,
                row.names = T)
    invisible(0)
}


# Function to create simulated mutational profiles.
create_simulated_profiles = function(nr_muts, nr_samples, signatures){
    
    #Determine number signatures to simulate
    nr_sigs = length(nr_muts)
    
    #Create profiles and real contris
    mut_mat_contri_l = purrr::map(seq_len(nr_samples), create_simulated_profile, nr_sigs, signatures, nr_muts)
    
    #Extract mutation matrix
    mut_mat = purrr::map(mut_mat_contri_l, "mut_mat") %>%
        do.call(cbind, .)
    
    #Create contribution matrix
    contribution = purrr::map(mut_mat_contri_l, "real_contri") %>%
        purrr::reduce(dplyr::full_join, by = "sigs")
    
    sig_ref = tibble::tibble("sigs" = colnames(signatures))
    contribution <- dplyr::left_join(sig_ref, contribution, by ="sigs") %>%
        as.data.frame()
    
    # Turn contribution into matrix and remove NAs
    rownames(contribution) <- contribution$sigs
    contribution <- contribution %>%
        dplyr::select(-sigs) %>%
        as.matrix()
    contribution[is.na(contribution)] <- 0
    
    return(list("mut_mat" = mut_mat, "contribution" = contribution))
}


create_simulated_profile = function(i, nr_sigs, signatures, nr_muts){
    
    sigs_i = sample.int(ncol(signatures), nr_sigs)
    sigs = signatures[, sigs_i, drop = FALSE]
    
    #Combine sigs based on weights from nr_muts
    weighted_sig = sigs %*% nr_muts
    weighted_sig = prop.table(weighted_sig)
    
    #Sample muts
    total_muts = sum(nr_muts)
    muts = sample(nrow(signatures), total_muts, prob = weighted_sig, replace = TRUE)
    
    #Create count table
    categories = tibble::tibble("feature" = seq_len(nrow(signatures)))
    counts = muts %>%
        table %>%
        tibble::enframe(name = "feature", value = "nr_muts") %>%
        dplyr::mutate(feature = as.integer(feature),
                      nr_muts = as.vector(nr_muts))
    counts = dplyr::left_join(categories, counts, by = "feature")
    counts$nr_muts[is.na(counts$nr_muts)] = 0
    
    #Turn into matrix
    mut_mat = counts %>%
        dplyr::select(nr_muts) %>%
        as.matrix()
    colnames(mut_mat) = NULL
    
    real_contri = tibble::tibble("sigs" = colnames(sigs), "nr_muts" = nr_muts)
    colnames(real_contri)[2] = as.character(i)
    
    return(list("mut_mat" = mut_mat, "real_contri" = real_contri))
}

