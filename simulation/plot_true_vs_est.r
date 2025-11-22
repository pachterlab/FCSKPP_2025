#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i", "--input_dir"), default="sim3")
)

args <- parse_args(OptionParser(option_list=option_list))

## Choose input pattern based on fit file format.
fit_files <- list.files(args$input_dir, pattern="^ml_fit_.*\\.rds$", full.names=TRUE)
# fit_files <- list.files(args$input_dir, pattern="^fit_.*\\.rds$", full.names=TRUE)

cat("Found", length(fit_files), "fit files in", args$input_dir, "\n")

records <- list()

## Adjust mapping between fitted and simulated parameters based on models.
# --- Mapping between fit and sim parameter names ---
mapping <- list(
  "alpha"          = "alpha",
  "beta"           = "beta",
  "sigma_x1"       = "sigma_x1",
  "sigma_x2"       = "sigma_x2",
  "sigma_noise1"   = "sigma_noise1",
  "sigma_noise2"   = "sigma_noise2",
  "p_noise"        = "p_noise",
  "theta_b_pop_mean"  = "theta_b_pop_mean",
  "theta_b_pop_sd"    = "theta_b_pop_sd",
  "theta_gamma_pop_mean" = "theta_gamma_pop_mean",
  "theta_gamma_pop_sd"   = "theta_gamma_pop_sd",
  "theta_mu_pop_mean"  = "theta_mu_pop_mean",
  "theta_mu_pop_sd"    = "theta_mu_pop_sd",
  "q"                 = "q"
  # "theta_gamma_pop_mean" = "theta_gamma_pop_mean",
  # "theta_gamma_pop_sd" = "theta_gamma_pop_sd"
  # "theta_gamma_pop_mean" = "theta_b_pop_mean",
  # "theta_gamma_pop_sd" = "theta_b_pop_sd"
  # "theta_b_pop_mean" = "theta_mu_pop_mean",
  # "theta_b_pop_sd" = "theta_mu_pop_sd"
  # "theta_b_pop_mean" = "theta_gamma_pop_mean",
  # "theta_b_pop_sd" = "theta_gamma_pop_sd"
)
cat("Defined mapping\n")

# --- For each fit file ---
for (ff in fit_files) {
  fit <- readRDS(ff)
  est <- fit$estimates
  base <- sub("^ml_fit_", "", tools::file_path_sans_ext(basename(ff)))
  # base <- sub("^fit_model2", "", tools::file_path_sans_ext(basename(ff)))
  sim_file <- file.path(args$input_dir, paste0(base, ".rds"))
  sim_data <- readRDS(sim_file)
  true <- sim_data$parameters

  # Apply mapping to align names
  est_names <- names(est)
  true_names <- sapply(est_names, function(n) {
    if (!is.null(mapping[[n]])) mapping[[n]] else n
  }, USE.NAMES = TRUE)

  rec <- data.frame(
    file = basename(ff),
    param = est_names,
    true = unlist(true[true_names]),
    est = unlist(est)
  )
  records[[length(records) + 1]] <- rec
}

cat("Compiled records from fit files\n")

df <- bind_rows(records)

df <- df %>%
  mutate(param_label = recode(param,
    "alpha"                 = "alpha[b]",
    "beta"                  = "alpha[gamma]",
    "q"                     = "q",
    "p_noise"               = "p[wn]",
    "sigma_noise1"          = "sigma[wn]^b",
    "sigma_noise2"          = "sigma[wn]^gamma",
    "sigma_x1"              = "sigma[b]",
    "sigma_x2"              = "sigma[gamma]",
    "theta_b_pop_mean"      = "bar(theta)[b]^P",
    "theta_b_pop_sd"        = "tau[b]^P",
    "theta_gamma_pop_mean"  = "bar(theta)[gamma]^P",
    "theta_gamma_pop_sd"    = "tau[gamma]^P",
    # "theta_mu_pop_mean"     = "bar(theta)[mu]^P",
    # "theta_mu_pop_sd"       = "tau[mu]^P"
    "theta_mu_pop_mean"     = "bar(theta)[b]^P",
    "theta_mu_pop_sd"       = "tau[b]^P"
  )) %>%
  mutate(param_label = factor(param_label, levels = unique(param_label)))

# --- Plot ---
output <- file.path(args$input_dir, paste0("fit_vs_true", args$input_dir, ".pdf"))
pdf(output, width=10, height=8)

p <- ggplot(df, aes(x=true, y=est)) +
  geom_point(size=1, alpha=0.7) +
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed") +
  facet_wrap(~param_label, scales="free", labeller = label_parsed) +
  labs(x="True", y="Estimated") +
  theme_bw(base_size=14)

print(p)
dev.off()

cat("Saved plot to", output, "\n")
