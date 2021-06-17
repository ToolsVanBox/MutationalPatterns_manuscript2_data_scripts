cost_contribution = function(real_contri, esti_contri){
    
    diff = abs(real_contri - esti_contri)
    
    # Fraction right
    frac = 1 - (sum(diff) / (sum(real_contri) + sum(esti_contri)))
    
    # diff = esti_contri - real_contri
    # diff[diff < 0] = 0
    # frac = 1 - (sum(diff) / sum(esti_contri))
    
    # Per sample
    frac_sample = 1 - (colSums(diff) / (colSums(real_contri) + colSums(esti_contri)))
    avg_frac_sample = mean(frac_sample)
    n = length(frac_sample)
    sd = sd(frac_sample)
    error_frac_sample = qt(0.975, df = n-1) * sd / sqrt(n)
    
    
    # Nr signatures misidentified
    id_real = real_contri != 0
    id_esti = esti_contri != 0
    
    TP = sum(id_esti & id_real)
    FP = sum(id_esti & !id_real)
    FN = sum(!id_esti & id_real)
    TN = sum(!id_esti & !id_real)
    
    precision = TP / (TP + FP)
    recall = TP / (TP + FN)
    specificity = TN / (TN + FP)
    accuracy = (TP + TN) / (TP + FP + FN + TN)
    F1 = 2 * precision * recall / (precision + recall)
    
    
    return(tibble::tibble("frac" = frac,
                          "frac_mean_sample" = avg_frac_sample,
                          "frac_error_sample" = error_frac_sample,
                          "precision" = precision,
                          "recall" = recall,
                          "specificity" = specificity,
                          "accuracy" = accuracy,
                          "F1" = F1))
}
