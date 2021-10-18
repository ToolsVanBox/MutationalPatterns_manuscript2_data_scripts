plot_simulation = function(out_dir){
    
    # Got to output/input directory
    old_dir = setwd(out_dir)
    on.exit(setwd(old_dir), add = TRUE)
    
    # Read cost files
    cost_tb = read_tsv("costs/cost_mutpatterns.txt")
    
    # Correct format
    cost_tb = cost_tb %>% 
        dplyr::mutate(matrix_name = factor(matrix_name, levels = unique(matrix_name)),
                      method = factor(method, levels = unique(method)))
    
    #Create figures
    fraction_correct_fig = plot_fraction_correct(cost_tb)
    fraction_mean_correct_fig = plot_mean_fraction_correct(cost_tb)
    prec_recal_fig = plot_precision_recall(cost_tb)
    delta_effect_fig = plot_delta_effect(cost_tb)
    f1_fig = plot_f1(cost_tb)
    message("Created figures")
    
    # Write out output
    readr::write_tsv(cost_tb, "cost_all.txt")
    ggsave("fraction_correct_fig.pdf", fraction_correct_fig)
    ggsave("fraction_mean_correct_fig.pdf", fraction_mean_correct_fig)
    ggsave("prec_recal_fig.pdf", prec_recal_fig)
    ggsave("delta_effect_fig.pdf", delta_effect_fig)
    ggsave("f1_fig.pdf", f1_fig)
    message("Wrote out output")
    
    invisible(0)
}


#Function to plot the fraction of correctly attributed mutations
plot_fraction_correct = function(cost_tb){
    cost_tb = cost_tb %>%
        dplyr::group_by(method, matrix_name) %>%
        dplyr::summarise(frac = max(frac))
    fig = ggplot(cost_tb, aes(x = matrix_name, y = frac, fill = method)) +
        geom_bar(position = "dodge", stat = "identity") +
        labs(y = "fraction correct") +
        theme_classic() +
        theme(text = element_text(size = 20))
    return(fig)
}

#Function to calculate the area under the curve.
calculate_auc = function(cost_tb){
    cost_tb = dplyr::select(cost_tb, precision, recall)
    max_x = max(cost_tb$precision)
    curve = approxfun(cost_tb, rule = 2, ties = max)
    auc = integrate(curve, 0, max_x)$value
    return(auc)
}

#Function to create a precision recall figure
plot_precision_recall = function(cost_tb){
    
    cost_tb_strict = dplyr::filter(cost_tb, method == "strict")
    cost_tb_regular = dplyr::filter(cost_tb, method == "regular")
    cost_tb_regular10 = dplyr::filter(cost_tb, method == "regular_10+")

    #Modify x values a little, to prevent duplicates when plotting for strict
    nr_per_group = cost_tb_strict$cutoff %>%
        unique() %>%
        length()
    small_numbers = seq(0, 1e-10, length.out = nr_per_group)
    cost_tb_strict_mod = cost_tb_strict %>%
        dplyr::group_by(matrix_name) %>%
        dplyr::mutate(precision = precision - small_numbers) %>%
        dplyr::ungroup()
    
    #Bind data back together
    cost_tb = rbind(cost_tb_strict_mod, cost_tb_regular, cost_tb_regular10)
    
    #determine aucs
    label_tb = create_label_tb(cost_tb_strict, "strict", 0.1)

    # Set size scale
    size_scale = scale_size_manual(breaks = c("strict", "regular", "regular_10+"), 
                                   values = c(1.5, 3, 3, 1.5), 
                                   guide = FALSE)
    
    #Create plot
    fig = ggplot(cost_tb, aes(x = precision, y = recall, colour = matrix_name, shape = method)) +
        geom_point(aes(size = method)) +
        geom_path(aes(linetype = method)) +
        geom_text(data = label_tb, label = label_tb$aucs_label, size = 5, show.legend = FALSE) +
        size_scale +
        coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
        theme_bw() +
        theme(text = element_text(size = 20))
    
    return(fig)
}

# Helper function for plot_precision_recall
create_label_tb = function(cost_tb, name, x_pos){
    aucs = cost_tb %>%
        split(., .$matrix_name) %>%
        purrr::map_dbl(calculate_auc)
    
    aucs_label = paste0("AUC: ", round(aucs, 3))
    label_tb = tibble::tibble("aucs_label" = aucs_label,
                              "matrix_name" = names(aucs),
                              "recall" = c(0.5, 0.4, 0.3, 0.2),
                              "precision" = rep(x_pos, 4),
                              "method" = name)
    return(label_tb)
}

#Function to plot the effect of the delta.
plot_delta_effect = function(cost_tb){
    
    cost_tb_strict = cost_tb %>% 
        dplyr::filter(method == "strict")
        
    
    opt_max_delta = cost_tb_strict %>%
        dplyr::group_by(matrix_name) %>%
        dplyr::filter(F1 == max(F1)) %>%
        dplyr::summarise(cutoff = paste0(round(cutoff, 5), collapse = ";")) %>%
        dplyr::mutate(cutoff_label = paste0("opt max_delta: ", cutoff),
                      cutoff = min(cost_tb_strict$cutoff)*4,
                      F1 = c(1, 0.95, 0.9, 0.85))
    
    fig = ggplot(cost_tb_strict, aes(x = log2(cutoff), y = F1, colour = matrix_name)) +
        geom_point() +
        geom_line() +
        geom_text(data = opt_max_delta,
                  label = opt_max_delta$cutoff_label,
                  show.legend = FALSE,
                  size = 5) +
        theme_bw() +
        theme(text = element_text(size = 20))
    
    return(fig)
}

#Function to plot the fraction of correctly attributed mutations
plot_fraction_correct = function(cost_tb){
    cost_tb = cost_tb %>%
        dplyr::group_by(method, matrix_name) %>%
        dplyr::summarise(frac = max(frac))
    fig = ggplot(cost_tb, aes(x = matrix_name, y = frac, fill = method)) +
        geom_bar(position = "dodge", stat = "identity") +
        labs(y = "fraction correct") +
        theme_classic() +
        theme(text = element_text(size = 20))
    return(fig)
}

#Function to plot the fraction of correctly attributed mutations
plot_mean_fraction_correct = function(cost_tb){
    cost_tb = cost_tb %>%
        dplyr::group_by(method, matrix_name) %>%
        dplyr::filter(frac_mean_sample == max(frac_mean_sample)) %>%
        dplyr::ungroup()
    fig = ggplot(cost_tb, aes(x = matrix_name,
                              y = frac_mean_sample,
                              fill = method)) +
        geom_bar(position = "dodge", stat = "identity") +
        geom_errorbar(aes(ymin = frac_mean_sample - frac_error_sample,
                          ymax = frac_mean_sample + frac_error_sample),
                      width = 0.3, position = position_dodge(width = 0.9)) +
        labs(y = "Mean fraction correct") +
        theme_classic() +
        theme(text = element_text(size = 20))
    return(fig)
}
plot_f1 = function(cost_tb){
    cost_tb = cost_tb %>%
        dplyr::group_by(method, matrix_name) %>%
        dplyr::summarise(F1 = max(F1))
    fig = ggplot(cost_tb, aes(x = matrix_name, y = F1, fill = method)) +
        geom_bar(position = "dodge", stat = "identity") +
        labs(y = "F1") +
        theme_classic() +
        theme(text = element_text(size = 20))
    return(fig)
}
