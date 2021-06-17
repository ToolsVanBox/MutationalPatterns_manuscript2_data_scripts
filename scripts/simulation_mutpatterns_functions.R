do_simulation_analysis_mutpatterns = function(signatures, out_dir, max_deltas){
    
    # Got to output directory
    old_dir = setwd(out_dir)
    on.exit(setwd(old_dir), add = TRUE)
    
    # Read simulated data
    res_l = readRDS("simulated_data/res_l.rds")
    
    #Calculate costs
    cost_tb_strict = simulated_costs(max_deltas, res_l, signatures, "strict")
    cost_tb_regular = simulated_costs(0, res_l, signatures, "regular")
    cost_tb_regular10 = simulated_costs(0, res_l, signatures, "regular_10+")
    cost_tb = rbind(cost_tb_strict, cost_tb_regular, cost_tb_regular10)
    message("Calculated costs")
    
    # Create output directory
    if (!dir.exists("costs")){
        dir.create("costs")
    }
    write_tsv(cost_tb, "costs/cost_mutpatterns.txt")
    invisible(0)
}

#Returns the cost for multiple max_deltas and multiple simulated matrixes
simulated_costs = function(max_deltas, res_l, signatures, method){
    cost_tb = purrr::map(res_l, deltas_simulation_cost, max_deltas, signatures, method) %>%
        dplyr::bind_rows(.id = "matrix_name") %>%
        dplyr::mutate(matrix_name = factor(matrix_name, levels = unique(matrix_name)))
    return(cost_tb)
}

#Returns the cost for multiple max_deltas and a single simulated maxtrix
deltas_simulation_cost = function(res, max_deltas, signatures, method){
    cost_tb = purrr::map(max_deltas, single_simulation_cost, res, signatures, method) %>%
        dplyr::bind_rows()
    message("Finished with one matrix")
    return(cost_tb)
}

#Returns the cost for a single max_delta and mut_mat
single_simulation_cost = function(max_delta, res, signatures, method = c("strict", "regular", "regular_10+")){
    
    method = match.arg(method)
    
    #Fit to signatures
    if (method == "strict"){
        fit_res_strict_l = fit_to_signatures_strict(res$mut_mat, signatures, max_delta = max_delta)
        fit_res = fit_res_strict_l$fit_res
    } else if (method == "regular"){
        fit_res = fit_to_signatures(res$mut_mat, signatures)
    } else if (method == "regular_10+"){
        fit_res = fit_to_signatures(res$mut_mat, signatures)
        fit_res$contribution[fit_res$contribution < 10] = 0
    }
    
    
    #Calculate cost
    cost_tb = cost_contribution(res$contribution, fit_res$contribution)
    cost_tb$cutoff = max_delta
    cost_tb$method = method
    
    message("Finished with one refit")
    return(cost_tb)
}
