#!/usr/bin/env Rscript

# Compare MCMC results across 3 fit files
# Creates box plots and overlapping posterior distributions

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
    library(gridExtra)
    library(dplyr)
    library(tidyr)
    library(RColorBrewer)
})

# Parse command line arguments
option_list <- list(
    make_option(c("-f", "--fit_files"), 
                help="Comma-separated list of 3 RDS fit files"),
    make_option(c("-o", "--output_prefix"), default="fit_comparison",
                help="Prefix for output files [default: %default]"),
    make_option(c("--plot_width"), type="numeric", default=12,
                help="Plot width in inches [default: %default]"),
    make_option(c("--plot_height"), type="numeric", default=8,
                help="Plot height in inches [default: %default]"),
    make_option(c("--dpi"), type="integer", default=1000,
                help="Plot resolution [default: %default]")
)
args <- parse_args(OptionParser(option_list=option_list))

if (is.null(args$fit_files)) {
    stop("Must specify --fit_files with 3 comma-separated file paths")
}

fit_files <- trimws(strsplit(args$fit_files, ",")[[1]])
if (length(fit_files) != 3) {
    stop("Must provide exactly 3 fit files")
}

# Custom theme
theme_comparison <- function() {
    theme_minimal() +
    theme(
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 13, face = "bold"),
        panel.grid.minor = element_blank()
    )
}

# Load all fit files
cat("Loading fit files...\n")
results <- list()
for (i in 1:3) {
    results[[i]] <- readRDS(fit_files[i])
    cat("Loaded:", fit_files[i], "\n")
}

# Extract key parameters from all fits
key_params <- c("alpha_out", "sigma_x_out", "sigma_noise_out", "p_noise_out", 
                "theta_pop_mean_out", "theta_pop_sd_out")

# Create combined data for box plots
box_data <- data.frame()
density_data <- data.frame()

for (i in 1:3) {
    posterior <- results[[i]]$posterior
    fit_name <- paste0("Fit_", i)
    
    for (param in key_params) {
        if (param %in% names(posterior)) {
            values <- posterior[[param]]
            
            # For box plots
            param_clean <- gsub("_out$", "", param)
            box_df <- data.frame(
                parameter = param_clean,
                value = values,
                fit = fit_name
            )
            box_data <- rbind(box_data, box_df)
            
            # For density plots
            density_df <- data.frame(
                parameter = param_clean,
                value = values,
                fit = fit_name
            )
            density_data <- rbind(density_data, density_df)
        }
    }
}

# Create box plots
cat("Creating box plots...\n")
p_box <- ggplot(box_data, aes(x = fit, y = value, fill = fit)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
    facet_wrap(~parameter, scales = "free_y", ncol = 2) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    labs(
        title = "Parameter Estimates Comparison (Box Plots)",
        x = "Fit",
        y = "Parameter Value",
        fill = "Fit"
    ) +
    theme_comparison()

# Create overlapping posterior density plots
cat("Creating posterior density plots...\n")
p_density <- ggplot(density_data, aes(x = value, fill = fit)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~parameter, scales = "free", ncol = 2) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    labs(
        title = "Overlapping Posterior Distributions",
        x = "Parameter Value",
        y = "Density",
        fill = "Fit"
    ) +
    theme_comparison()

# Create summary statistics table
cat("Creating summary statistics...\n")
summary_stats <- box_data %>%
    group_by(parameter, fit) %>%
    summarise(
        mean = mean(value),
        median = median(value),
        sd = sd(value),
        q025 = quantile(value, 0.025),
        q975 = quantile(value, 0.975),
        .groups = 'drop'
    )

# Save plots
box_file <- paste0(args$output_prefix, "_boxplots.png")
ggsave(box_file, p_box, 
       width = args$plot_width, height = args$plot_height, 
       dpi = args$dpi)
cat("Saved box plots:", box_file, "\n")

density_file <- paste0(args$output_prefix, "_posteriors.png")
ggsave(density_file, p_density, 
       width = args$plot_width, height = args$plot_height, 
       dpi = args$dpi)
cat("Saved posterior plots:", density_file, "\n")

# Create combined plot
combined_plot <- grid.arrange(p_box, p_density, ncol = 1, nrow = 2)
combined_file <- paste0(args$output_prefix, "_combined.png")
ggsave(combined_file, combined_plot, 
       width = args$plot_width, height = args$plot_height * 2, 
       dpi = args$dpi)
cat("Saved combined plot:", combined_file, "\n")

# Save summary table
summary_file <- paste0(args$output_prefix, "_summary.csv")
write.csv(summary_stats, summary_file, row.names = FALSE)
cat("Saved summary statistics:", summary_file, "\n")

# Print summary to console
cat("\n=== PARAMETER COMPARISON SUMMARY ===\n")
summary_display <- summary_stats
numeric_cols <- sapply(summary_display, is.numeric)
summary_display[numeric_cols] <- lapply(summary_display[numeric_cols], function(x) round(x, 4))
print(summary_display)

# Calculate pairwise differences
cat("\n=== PAIRWISE DIFFERENCES (Fit_1 vs Fit_2, Fit_1 vs Fit_3, Fit_2 vs Fit_3) ===\n")
wide_summary <- summary_stats %>%
    select(parameter, fit, mean) %>%
    pivot_wider(names_from = fit, values_from = mean)

if (all(c("Fit_1", "Fit_2", "Fit_3") %in% colnames(wide_summary))) {
    differences <- wide_summary %>%
        mutate(
            Fit1_vs_Fit2 = Fit_1 - Fit_2,
            Fit1_vs_Fit3 = Fit_1 - Fit_3,
            Fit2_vs_Fit3 = Fit_2 - Fit_3
        ) %>%
        select(parameter, Fit1_vs_Fit2, Fit1_vs_Fit3, Fit2_vs_Fit3)
    
    differences[2:4] <- lapply(differences[2:4], function(x) round(x, 4))
    print(differences)
}

cat("\nComparison complete!\n")
cat("Files created:\n")
cat("- ", box_file, "\n")
cat("- ", density_file, "\n")
cat("- ", combined_file, "\n")
cat("- ", summary_file, "\n")