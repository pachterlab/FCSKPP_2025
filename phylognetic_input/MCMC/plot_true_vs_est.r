#!/usr/bin/env Rscript

# Fixed Parameter Recovery Analysis for OU-Gaussian Mixture Model
# Matches the direct parameterization from the fitting script

# Load required libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(gridExtra)
    library(scales)
    library(rlang)
    library(tidyr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_dir <- if(length(args) > 0) args[1] else "simulations"
output_prefix <- if(length(args) > 1) args[2] else "parameter_recovery"
use_existing <- if(length(args) > 2) {
    arg3 <- toupper(trimws(args[3]))
    arg3 %in% c("TRUE", "T", "1", "YES")
} else TRUE

cat("Analyzing parameter recovery in directory:", input_dir, "\n")
cat("Use existing parameter data:", use_existing, "\n")

# Check if we should use existing parameter data
param_file <- file.path(input_dir, paste0(output_prefix, "_data.csv"))

if (use_existing && file.exists(param_file)) {
    cat("Loading existing parameter data from:", param_file, "\n")
    param_data <- read.csv(param_file, stringsAsFactors = FALSE)
    cat("Loaded parameters from", nrow(param_data), "files\n")
} else {
    # Extract parameters from fit files (existing code)
    cat("Extracting parameters from fit files...\n")
    
    # Find all fit files
    fit_files <- list.files(input_dir, pattern = "bayesian_fit_ou_gaussian_.*\\.rds$", full.names = TRUE)
    cat("Found", length(fit_files), "fit files\n")

    if (length(fit_files) == 0) {
        stop("No fit files found. Make sure the fitting script has been run.")
    }

    safe1 <- function(x) {
        if (length(x) == 1 && !is.na(x)) x else NA_real_
    }

    # Extract parameters and uncertainties from all files
    extract_params <- function(fit_file) {
        tryCatch({
            result <- readRDS(fit_file)
            posterior <- result$posterior
            true_params <- result$true_parameters
            
            cat("Processing:", basename(fit_file), "\n")
            cat("True parameter names:", names(true_params), "\n")
            cat("Posterior parameter names:", names(posterior)[1:min(10, length(names(posterior)))], "...\n")
            
            # Extract posterior summaries
            extract_summary <- function(param_name) {
                if (param_name %in% names(posterior)) {
                    values <- posterior[[param_name]]
                    # Remove infinite and extremely large values
                    values <- values[is.finite(values) & abs(values) < 1e6]
                    if (length(values) == 0) {
                        return(data.frame(mean = NA, median = NA, sd = NA, q025 = NA, q975 = NA))
                    }
                    data.frame(
                        mean = mean(values, na.rm = TRUE),
                        median = median(values, na.rm = TRUE),
                        sd = sd(values, na.rm = TRUE),
                        q025 = quantile(values, 0.025, na.rm = TRUE),
                        q975 = quantile(values, 0.975, na.rm = TRUE)
                    )
                } else {
                    data.frame(mean = NA, median = NA, sd = NA, q025 = NA, q975 = NA)
                }
            }
            
            # FIXED: Use the correct parameter names from the direct parameterization
            alpha <- extract_summary("alpha_out")           # Direct alpha parameter
            sigma_x <- extract_summary("sigma_x_out")       # Direct sigma_x parameter  
            sigma_noise <- extract_summary("sigma_noise_out") # Direct sigma_noise parameter

            ## Change this - currently post-hoc adding the correction.
            # # Correct approach: add theta_prior_var to each posterior sample, then summarize
            # if ("sigma_noise_out" %in% names(posterior) && !is.null(result$stan_data$theta_prior_sd)) {
            #     corrected_samples <- sqrt(posterior[["sigma_noise_out"]]^2 + result$stan_data$theta_prior_sd^2)
            #     sigma_noise <- data.frame(
            #         mean = mean(corrected_samples, na.rm = TRUE),
            #         median = median(corrected_samples, na.rm = TRUE), 
            #         sd = sd(corrected_samples, na.rm = TRUE),
            #         q025 = quantile(corrected_samples, 0.025, na.rm = TRUE),
            #         q975 = quantile(corrected_samples, 0.975, na.rm = TRUE)
            #     )
            # } else {
            #     sigma_noise <- data.frame(mean = NA, median = NA, sd = NA, q025 = NA, q975 = NA)
            # }

            p_noise <- extract_summary("p_noise_out")       # Direct p_noise parameter
            
            # Optional derived quantities
            phi <- extract_summary("phi")                   # Derived phylogenetic signal
            sigma_sq <- extract_summary("sigma_sq")         # Derived sigma_sq = sigma_x^2
            stationary_var <- extract_summary("stationary_variance") # Derived V = sigma_x^2/(2*alpha)
            
            # Compile results matching true parameter structure
            data.frame(
                file = basename(fit_file),

                # True values using actual parameter names from simulation
                true_H = safe1(true_params$H),                    # This is alpha in the model
                true_Sigma_x = safe1(true_params$Sigma_x),        # Diffusion parameter  
                true_Sigma_noise = safe1(true_params$Sigma_noise), # Noise standard deviation
                true_P_noise = safe1(true_params$P_noise),        # Noise proportion
                true_Theta = safe1(true_params$Theta),            # Optimum value
                true_X0 = safe1(true_params$X0),                  # Root state
                
                # Convert true parameters to derived quantities for comparison
                true_phi = safe1({
                    if (!is.null(true_params$H) && !is.null(result$stan_data$tree_height)) {
                        alpha_val <- true_params$H
                        tree_height <- result$stan_data$tree_height
                        x <- 2 * alpha_val * tree_height
                        if (x < 1e-10) {
                            alpha_val * tree_height  # Linear approximation for small x
                        } else {
                            1 - (1 - exp(-x)) / x
                        }
                    } else NA_real_
                }),
                true_V = safe1({
                    if (!is.null(true_params$Sigma_x) && !is.null(true_params$H)) {
                        true_params$Sigma_x^2 / (2 * true_params$H)
                    } else NA_real_
                }),

                # FIXED: Use direct parameter estimates (no more inf/-inf issues)
                est_alpha = safe1(alpha$mean),
                est_alpha_sd = safe1(alpha$sd),
                est_alpha_lower = safe1(alpha$q025),
                est_alpha_upper = safe1(alpha$q975),

                est_sigma_x = safe1(sigma_x$mean),
                est_sigma_x_sd = safe1(sigma_x$sd),
                est_sigma_x_lower = safe1(sigma_x$q025),
                est_sigma_x_upper = safe1(sigma_x$q975),

                est_sigma_noise = safe1(sigma_noise$mean),
                est_sigma_noise_sd = safe1(sigma_noise$sd),
                est_sigma_noise_lower = safe1(sigma_noise$q025),
                est_sigma_noise_upper = safe1(sigma_noise$q975),

                est_p_noise = safe1(p_noise$mean),
                est_p_noise_sd = safe1(p_noise$sd),
                est_p_noise_lower = safe1(p_noise$q025),
                est_p_noise_upper = safe1(p_noise$q975),
                
                # Optional derived quantities
                est_phi = safe1(phi$mean),
                est_phi_sd = safe1(phi$sd),
                est_phi_lower = safe1(phi$q025),
                est_phi_upper = safe1(phi$q975),
                
                est_V = safe1(stationary_var$mean),
                est_V_sd = safe1(stationary_var$sd),
                est_V_lower = safe1(stationary_var$q025),
                est_V_upper = safe1(stationary_var$q975),

                # Convergence diagnostics
                max_rhat = safe1(result$diagnostics$max_rhat),
                min_ess = safe1(result$diagnostics$min_ess),
                divergences = safe1(result$diagnostics$divergences)
            )
        }, error = function(e) {
            cat("Error processing", basename(fit_file), ":", e$message, "\n")
            return(NULL)
        })
    }

    # Process all files
    param_list <- lapply(fit_files, extract_params)
    param_data <- do.call(rbind, param_list[!sapply(param_list, is.null)])

    if (nrow(param_data) == 0) {
        stop("No valid parameter data extracted from fit files.")
    }

    cat("Successfully extracted parameters from", nrow(param_data), "files\n")

    # Remove any remaining infinite values
    param_data <- param_data %>%
        mutate_if(is.numeric, function(x) ifelse(is.infinite(x), NA, x))

    # Save parameter data
    write.csv(param_data, param_file, row.names = FALSE)
    cat("Parameter data saved to:", param_file, "\n")
}

# Create recovery plots
cat("Creating parameter recovery plots...\n")

# Print diagnostics
cat("Checking parameter ranges...\n")
cat("Alpha - True range:", paste(round(range(param_data$true_H, na.rm = TRUE), 3), collapse = " to "), 
    ", Estimated range:", paste(round(range(param_data$est_alpha, na.rm = TRUE), 3), collapse = " to "), "\n")
cat("Sigma_x - True range:", paste(round(range(param_data$true_Sigma_x, na.rm = TRUE), 3), collapse = " to "), 
    ", Estimated range:", paste(round(range(param_data$est_sigma_x, na.rm = TRUE), 3), collapse = " to "), "\n")
cat("Sigma_noise - True range:", paste(round(range(param_data$true_Sigma_noise, na.rm = TRUE), 3), collapse = " to "), 
    ", Estimated range:", paste(round(range(param_data$est_sigma_noise, na.rm = TRUE), 3), collapse = " to "), "\n")
cat("P_noise - True range:", paste(round(range(param_data$true_P_noise, na.rm = TRUE), 3), collapse = " to "), 
    ", Estimated range:", paste(round(range(param_data$est_p_noise, na.rm = TRUE), 3), collapse = " to "), "\n")

# Custom plotting function with log-log scale support
plot_recovery <- function(data, true_col, est_col, lower_col, upper_col, 
                         title, xlab, ylab, log_scale = TRUE) {
    
    # Filter out NA values, infinite values, and non-positive values for log scale
    plot_data <- data[!is.na(data[[true_col]]) & !is.na(data[[est_col]]) & 
                     is.finite(data[[true_col]]) & is.finite(data[[est_col]]), ]
    
    if (log_scale) {
        # For log scale, also remove non-positive values
        plot_data <- plot_data[plot_data[[true_col]] > 0 & plot_data[[est_col]] > 0 & 
                              plot_data[[lower_col]] > 0 & plot_data[[upper_col]] > 0, ]
    }
    
    if (nrow(plot_data) == 0) {
        return(ggplot() + 
               labs(title = paste(title, "(No Valid Data)")) + 
               theme_minimal())
    }
    
    # Remove extreme outliers (more conservative approach)
    if (nrow(plot_data) > 5) {
        # Use quantile-based filtering
        if (log_scale) {
            # For log scale, work in log space for outlier detection
            log_true <- log10(plot_data[[true_col]])
            log_est <- log10(plot_data[[est_col]])
            
            true_q01 <- quantile(log_true, 0.01, na.rm = TRUE)
            true_q99 <- quantile(log_true, 0.99, na.rm = TRUE)
            est_q01 <- quantile(log_est, 0.01, na.rm = TRUE)
            est_q99 <- quantile(log_est, 0.99, na.rm = TRUE)
            
            n_before <- nrow(plot_data)
            plot_data <- plot_data[
                log_true >= true_q01 & log_true <= true_q99 &
                log_est >= est_q01 & log_est <= est_q99, ]
        } else {
            # Linear space outlier detection
            true_q01 <- quantile(plot_data[[true_col]], 0.01, na.rm = TRUE)
            true_q99 <- quantile(plot_data[[true_col]], 0.99, na.rm = TRUE)
            est_q01 <- quantile(plot_data[[est_col]], 0.01, na.rm = TRUE)
            est_q99 <- quantile(plot_data[[est_col]], 0.99, na.rm = TRUE)
            
            n_before <- nrow(plot_data)
            plot_data <- plot_data[
                plot_data[[true_col]] >= true_q01 & plot_data[[true_col]] <= true_q99 &
                plot_data[[est_col]] >= est_q01 & plot_data[[est_col]] <= est_q99, ]
        }
        
        n_removed <- n_before - nrow(plot_data)
        if (n_removed > 0) {
            title <- paste0(title, " (", n_removed, " extreme outliers removed)")
        }
    }
    
    if (nrow(plot_data) == 0) {
        return(ggplot() + 
               labs(title = paste(title, "(All data were outliers)")) + 
               theme_minimal())
    }
    
    p <- ggplot(plot_data, aes(x = .data[[true_col]], y = .data[[est_col]])) +
        geom_errorbar(aes(ymin = .data[[lower_col]], ymax = .data[[upper_col]]), 
                     alpha = 0.3, width = 0, color = "gray50") +
        geom_point(aes(color = max_rhat > 1.1), alpha = 0.7, size = 2) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
        scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                          labels = c("TRUE" = "Poor convergence", "FALSE" = "Good convergence"),
                          name = "Rhat > 1.1") +
        labs(title = title, x = xlab, y = ylab) +
        theme_minimal() +
        theme(legend.position = "bottom",
              plot.title = element_text(size = 11, hjust = 0.5))
    
    if (log_scale) {
        p <- p + 
            scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
            scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
            annotation_logticks()
    }
    
    return(p)
}

# Main parameter recovery plots (4 key parameters) - ALL with log-log scale
p1 <- plot_recovery(param_data, "true_H", "est_alpha", 
                   "est_alpha_lower", "est_alpha_upper",
                   "Alpha Parameter Recovery (Log-Log)", 
                   "True H (α)", "Estimated α", 
                   log_scale = TRUE)

p2 <- plot_recovery(param_data, "true_Sigma_x", "est_sigma_x", 
                   "est_sigma_x_lower", "est_sigma_x_upper",
                   "Sigma_x Recovery (Log-Log)", 
                   "True Sigma_x", "Estimated Sigma_x", 
                   log_scale = TRUE)

p3 <- plot_recovery(param_data, "true_Sigma_noise", "est_sigma_noise", 
                   "est_sigma_noise_lower", "est_sigma_noise_upper",
                   "Sigma_noise Recovery (Log-Log)", 
                   "True Sigma_noise", "Estimated Sigma_noise", 
                   log_scale = TRUE)

p4 <- plot_recovery(param_data, "true_P_noise", "est_p_noise", 
                   "est_p_noise_lower", "est_p_noise_upper",
                   "Noise Proportion Recovery (Linear)", 
                   "True P_noise", "Estimated P_noise", 
                   log_scale = FALSE)  # P_noise is a proportion, so linear scale is more appropriate

# Combine main plots
main_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

# Save main recovery plot
main_plot_file <- file.path(input_dir, paste0(output_prefix, "_main_loglog.jpeg"))
ggsave(main_plot_file, main_plot, width = 14, height = 12, dpi = 2000)
cat("Log-log scale recovery plot saved to:", main_plot_file, "\n")

# Calculate and print recovery statistics
param_data <- param_data %>%
    mutate(
        alpha_rel_error = abs(est_alpha - true_H) / abs(true_H),
        sigma_x_rel_error = abs(est_sigma_x - true_Sigma_x) / abs(true_Sigma_x),
        sigma_noise_rel_error = abs(est_sigma_noise - true_Sigma_noise) / abs(true_Sigma_noise),
        p_noise_rel_error = abs(est_p_noise - true_P_noise) / abs(true_P_noise),
        good_convergence = max_rhat <= 1.1
    )

# Print summary statistics
cat("\n=== PARAMETER RECOVERY SUMMARY (LOG-LOG SCALE) ===\n")
cat("Files processed:", nrow(param_data), "\n")
cat("Good convergence (Rhat <= 1.1):", sum(param_data$good_convergence, na.rm = TRUE), "/", nrow(param_data), "\n")

cat("\nMean relative errors:\n")
cat("  Alpha (H):", sprintf("%.1f%%", mean(param_data$alpha_rel_error, na.rm = TRUE) * 100), "\n")
cat("  Sigma_x:", sprintf("%.1f%%", mean(param_data$sigma_x_rel_error, na.rm = TRUE) * 100), "\n")
cat("  Sigma_noise:", sprintf("%.1f%%", mean(param_data$sigma_noise_rel_error, na.rm = TRUE) * 100), "\n")
cat("  P_noise:", sprintf("%.1f%%", mean(param_data$p_noise_rel_error, na.rm = TRUE) * 100), "\n")

cat("\nCorrelations between true and estimated:\n")
cat("  Alpha (H):", sprintf("%.3f", cor(param_data$true_H, param_data$est_alpha, use = "complete.obs")), "\n")
cat("  Sigma_x:", sprintf("%.3f", cor(param_data$true_Sigma_x, param_data$est_sigma_x, use = "complete.obs")), "\n")
cat("  Sigma_noise:", sprintf("%.3f", cor(param_data$true_Sigma_noise, param_data$est_sigma_noise, use = "complete.obs")), "\n")
cat("  P_noise:", sprintf("%.3f", cor(param_data$true_P_noise, param_data$est_p_noise, use = "complete.obs")), "\n")

cat("\nAnalysis complete with log-log scale plots!\n")
cat("Usage: Rscript script.R [input_dir] [output_prefix] [use_existing_data]\n")
cat("  use_existing_data: TRUE to use existing CSV, FALSE to re-extract (default: FALSE)\n")
