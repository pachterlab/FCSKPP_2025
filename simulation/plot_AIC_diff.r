#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
})

# ---- Command line arguments ----
option_list <- list(
  make_option(c("-a", "--dir1"), default="sim1",
              help="First directory containing fit .rds files"),
  make_option(c("-b", "--dir2"), default="sim1_fit2",
              help="Second directory containing fit .rds files"),
  make_option(c( "--m1"), default='model 2', # remember you might need to switch these!
              help="Model 1 name for labeling"),
  make_option(c( "--m2"), default='model 1',
              help="Model 2 name for labeling")
)
args <- parse_args(OptionParser(option_list=option_list))

# If not names are given for the models, use directory names
if (is.null(args$m1)) {
  m1 <- args$dir1
} else {
  m1 <- args$m1
}
if (is.null(args$m2)) {
  m2 <- args$dir2
} else {
  m2 <- args$m2
}

## Adjust these based on the fit-file names in the two directories being compared.
pattern1 <- "^ml_fit_.*\\.rds$"
# pattern1 <- "^fit_model2"
pattern2 <- "^ml_fit_.*\\.rds$"
# pattern2 <- "^fit_model2"

fit_files1 <- list.files(args$dir1, pattern=pattern1, full.names=TRUE)
fit_files2 <- list.files(args$dir2, pattern=pattern2, full.names=TRUE)

base_name_fit1 <- function(path) sub(pattern1, "", tools::file_path_sans_ext(basename(path)))
base_name_fit2 <- function(path) sub(pattern2, "", tools::file_path_sans_ext(basename(path)))

bases1 <- sapply(fit_files1, base_name_fit1)
bases2 <- sapply(fit_files2, base_name_fit2)
common <- intersect(bases1, bases2)

records <- lapply(common, function(b) {
  f1 <- fit_files1[bases1 == b]
  f2 <- fit_files2[bases2 == b]
  fit1 <- readRDS(f1)
  fit2 <- readRDS(f2)
  data.frame(
    file = b,
    AIC_dir1 = fit1$optimization$aic,
    AIC_dir2 = fit2$optimization$aic,
    delta_AIC = fit2$optimization$aic - fit1$optimization$aic
  )
})
df <- bind_rows(records)

cat("Computed AIC differences for", nrow(df), "fits\n")

# ---- Histogram of ΔAIC ----
df$delta_AIC_capped <- pmin(df$delta_AIC, 500)

breaks <- c(seq(min(df$delta_AIC_capped), 500, by=20))
labels <- ifelse(breaks == 500, ">500", as.character(breaks))

p1 <- ggplot(df, aes(x=delta_AIC_capped)) +
  geom_histogram(
    binwidth=20, boundary=0,
    color="black", fill="skyblue", alpha=0.7
  ) +
  geom_vline(xintercept=0, color="red", linetype="dashed") +
  scale_x_continuous(
    breaks = c(seq(0, 500, by=100)),
    labels = c("0", "100", "200", "300", "400", ">500")
  ) +
  labs( # format the labels with expressions, with the given model names
    x = bquote(Delta*AIC == AIC(.(m2)) - AIC(.(m1))),
    y = "Count", # also format title with model names
    title = bquote(Delta*AIC~" for simulation under "*.(m1)*", with three-parameter H")
  ) +
  theme_bw(base_size=14)

# ---- % Correct vs cutoff ---- 
cutoffs <- seq(0, 100, by=5)
accuracy <- sapply(cutoffs, function(cut) {
  num_greater <- sum(df$delta_AIC > cut)
  num_less <- sum(df$delta_AIC < -cut)
  (num_greater / (num_greater + num_less)) * 100
})
acc_df <- data.frame(cutoff=cutoffs, correct=accuracy)

p2 <- ggplot(acc_df, aes(x=cutoff, y=correct)) +
  geom_point(size=1.5, color="darkblue") +
  labs(
    x = expression(AIC~"cutoff"),
    y = "% correctly identified",
    title = bquote("Correct identification vs AIC cutoff: Simulated "*.(m1))
  ) +
  theme_bw(base_size=14)

# ---- Save plots ----
output1 <- paste0("AIC_difference_", args$dir1, "_vs_", args$dir2, ".pdf")
output2 <- paste0("AIC_accuracy_", args$dir1, "_vs_", args$dir2, ".pdf")

pdf(output1, width=7, height=5); print(p1); dev.off()
pdf(output2, width=7, height=5); print(p2); dev.off()

cat("Saved:\n", output1, "\n", output2, "\n")
