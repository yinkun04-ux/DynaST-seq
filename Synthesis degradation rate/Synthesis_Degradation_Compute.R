library(dplyr)     
library(tidyverse) 
library(data.table)

compute_degradation_rate <- function(half_life_matrix) {
    gamma <- log(2) / half_life_matrix
    return(gamma)
}

# Step 1: Calculate RNA degradation rate (γ), unit in 1/h
gamma_matrix <- compute_degradation_rate(half_life_impute_result)

# Step 2: Calculate RNA synthesis rate (α), unit in cp10k molecules/h
# Using the formula: α = (n * γ) / (1 - exp(-γ * t))

calculate_synthesis_rate <- function(new_rna_mean, gamma_matrix, labeling_time) {
    synthesis_rate <- (new_rna_mean * gamma_matrix) / (1 - exp(-gamma_matrix * labeling_time))
    return(synthesis_rate)
}

new_rna_mean_result_sub <- new_rna_mean_result[rownames(gamma_matrix),]

labeling_time <- 90  # Time in minutes
synthesis_rate_matrix <- calculate_synthesis_rate(new_rna_mean_result_sub, gamma_matrix, labeling_time)