#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
    library(optparse)
    library(rstan)
    library(ape)
    library(parallel)
    library(digest)
})

# Set up optimized BLAS/LAPACK for matrix operations
Sys.setenv("OMP_NUM_THREADS" = "1")
Sys.setenv("MKL_NUM_THREADS" = "1")
Sys.setenv("OPENBLAS_NUM_THREADS" = "1")

# RStan configuration for stability
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
rstan_options(threads_per_chain = 1)

# Set up logging
LOG_FILE <- paste0("bayesian_ou_gaussian_mixture_poptheta_realdata_", Sys.Date(), ".log")
log_msg <- function(msg) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    entry <- paste(timestamp, "-", msg, "\n")
    cat(entry)
    cat(entry, file = LOG_FILE, append = TRUE)
}

log_msg("=== REAL DATA POPULATION THETA SCRIPT STARTED ===")

# Parse arguments
option_list <- list(
    make_option(c("-t", "--trait_file"), default="trait_data.csv", 
                help="CSV file with trait data (rows=species, cols=traits)"),
    make_option(c("-e", "--error_file"), default="measurement_errors.csv",
                help="CSV file with measurement errors (same structure as trait file)"),
    make_option(c("-n", "--newick_file"), default="phylogeny.tre",
                help="Newick format phylogenetic tree file"),
    make_option(c("-o", "--output_prefix"), default="realdata_fit",
                help="Prefix for output files"),
    make_option(c("--chains"), type="integer", default=4),
    make_option(c("--iter"), type="integer", default=2000),
    make_option(c("--warmup"), type="integer", default=1000),
    make_option(c("--cores"), type="integer", default=parallel::detectCores()),
    make_option(c("--recompile_model"), type="logical", default=FALSE),
    make_option(c("--thin"), type="integer", default=1),
    make_option(c("--species_col"), default="species", 
                help="Column name for species identifiers"),
    make_option(c("--skip_cols"), default="", 
                help="Comma-separated list of columns to skip (besides species column)")
)
args <- parse_args(OptionParser(option_list=option_list))

log_msg(paste("=== INPUT FILES ==="))
log_msg(paste("Trait data file:", args$trait_file))
log_msg(paste("Error data file:", args$error_file))
log_msg(paste("Phylogeny file:", args$newick_file))
log_msg(paste("Output prefix:", args$output_prefix))
log_msg(paste("Species column:", args$species_col))

# Parse skip columns
skip_cols <- character(0)
if (args$skip_cols != "") {
    skip_cols <- trimws(strsplit(args$skip_cols, ",")[[1]])
    log_msg(paste("Skipping columns:", paste(skip_cols, collapse=", ")))
}

# Stan model code (same as before)
stan_model_code <- "
functions {
    // Efficient OU covariance matrix computation
    matrix build_ou_cov_matrix(matrix vcv_bm, real alpha, real sigma_x, int N) {
        matrix[N, N] vcv_ou;
        real sigma_sq = sigma_x^2;
        
        for (i in 1:N) {
            vcv_ou[i, i] = sigma_sq * (1 - exp(-2 * alpha * vcv_bm[i, i])) / (2 * alpha);
            for (j in 1:(i-1)) {
                real t_mrca = vcv_bm[i, j];
                real t_i = vcv_bm[i, i];
                real t_j = vcv_bm[j, j];
                real cov_val = sigma_sq * 
                              (exp(-alpha * (t_i + t_j - 2 * t_mrca)) - 
                               exp(-alpha * (t_i + t_j))) / (2 * alpha);
                vcv_ou[i, j] = cov_val;
                vcv_ou[j, i] = cov_val;
            }
        }
        return vcv_ou;
    }
    
    // Marginal likelihood with theta integrated out for OU process
    real ou_marginal_likelihood_lpdf(vector y, vector se, matrix vcv_bm, real alpha, 
                                     real sigma_x, real theta_pop_mean, real theta_pop_var,
                                     int N) {
        matrix[N, N] vcv_ou = build_ou_cov_matrix(vcv_bm, alpha, sigma_x, N);
        matrix[N, N] Sigma_obs;
        vector[N] mu_theta;
        
        // Build observed covariance (OU + measurement error)
        Sigma_obs = vcv_ou;
        for (i in 1:N) {
            Sigma_obs[i, i] += se[i]^2;
        }
        
        // Mean vector coefficient (how theta affects each tip)
        for (i in 1:N) {
            mu_theta[i] = 1 - exp(-alpha * vcv_bm[i, i]);
        }
        
        // Marginal likelihood (integrating over theta ~ N(theta_pop_mean, theta_pop_var))
        matrix[N, N] Sigma_marginal = Sigma_obs + theta_pop_var * (mu_theta * mu_theta');
        vector[N] mu_marginal = mu_theta * theta_pop_mean;
        
        return multi_normal_lpdf(y | mu_marginal, Sigma_marginal);
    }
    
    // Log mixture likelihood for a single trait
    real mixture_likelihood_lpdf(vector y, vector se, matrix vcv_bm, 
                                real alpha, real sigma_x, real sigma_noise,
                                real p_noise, real theta_pop_mean, real theta_pop_var,
                                int N) {
        
        // Log-likelihood under OU process
        real ll_ou = ou_marginal_likelihood_lpdf(y | se, vcv_bm, alpha, sigma_x, 
                                                theta_pop_mean, theta_pop_var, N);
        
        // Log-likelihood under Gaussian noise (integrate over theta)
        vector[N] mu_noise = rep_vector(theta_pop_mean, N);
        vector[N] sigma_noise_total;
        for (i in 1:N) {
            sigma_noise_total[i] = sqrt(sigma_noise^2 + theta_pop_var + se[i]^2);
        }
        real ll_noise = normal_lpdf(y | mu_noise, sigma_noise_total);
        
        // Log mixture
        return log_mix(p_noise, ll_noise, ll_ou);
    }
}

data {
    int N; int K;
    vector[N] y[K];
    matrix[N, N] vcv_bm;
    vector[N] se[K];
    real tree_height;
    
    // Initial estimates for theta population parameters
    real theta_init_mean;
    real theta_init_sd;
}

transformed data {
    real y_grand_mean;
    real y_grand_var;
    
    {
        vector[N*K] all_y;
        int idx = 1;
        for (k in 1:K) {
            for (i in 1:N) {
                all_y[idx] = y[k][i];
                idx += 1;
            }
        }
        y_grand_mean = mean(all_y);
        y_grand_var = variance(all_y);
    }
}

parameters {
    real<lower=0.01> alpha;                // Prevent alpha from going to 0
    real<lower=0> sigma_x;                 
    real<lower=0> sigma_noise;             
    real<lower=0, upper=1> p_noise;
    
    // Population theta parameters
    real theta_pop_mean;                   // Population mean of theta
    real<lower=0> theta_pop_sd;            // Population SD of theta
}

transformed parameters {
    real theta_pop_var = theta_pop_sd^2;
}

model {
    // Priors on main parameters
    alpha ~ lognormal(log(2.0/tree_height), 0.5);
    sigma_x ~ lognormal(log(sqrt(y_grand_var)), 0.5);
    sigma_noise ~ lognormal(log(sqrt(y_grand_var/2)), 0.5);
    p_noise ~ beta(2, 2);
    
    // Priors on population theta parameters
    theta_pop_mean ~ normal(theta_init_mean, theta_init_sd * 2);
    theta_pop_sd ~ exponential(2.0);  // Regularizing prior
    
    // Mixture likelihood using population theta parameters
    for (k in 1:K) {
        y[k] ~ mixture_likelihood(se[k], vcv_bm, 
                                 alpha, sigma_x, sigma_noise,
                                 p_noise, theta_pop_mean, theta_pop_var, N);
    }
}

generated quantities {
    // Primary outputs
    real alpha_out = alpha;
    real sigma_x_out = sigma_x;
    real sigma_noise_out = sigma_noise;
    real p_noise_out = p_noise;
    real theta_pop_mean_out = theta_pop_mean;
    real theta_pop_sd_out = theta_pop_sd;
    
    // Derived quantities
    real sigma_sq = sigma_x^2;
    real phylogenetic_half_life = log(2) / alpha;
    real stationary_variance = sigma_sq / (2 * alpha);
    real phi = 1 - (1 - exp(-2 * alpha * tree_height)) / (2 * alpha * tree_height);
    real p_ou_out = 1 - p_noise;
    
    // Posterior classification: probability each trait comes from noise process
    vector[K] trait_noise_prob;
    vector[K] trait_ou_prob;
    
    // Posterior theta estimates for each process and trait
    vector[K] theta_ou_posterior_mean;
    vector[K] theta_ou_posterior_sd;
    vector[K] theta_noise_posterior_mean;
    vector[K] theta_noise_posterior_sd;
    
    // Best-guess theta (weighted by process probabilities)
    vector[K] theta_mixture_mean;
    vector[K] theta_mixture_sd;
    
    {
        for (k in 1:K) {
            // Compute log-likelihoods for classification
            real ll_ou = ou_marginal_likelihood_lpdf(y[k] | se[k], vcv_bm, alpha, sigma_x,
                                                    theta_pop_mean, theta_pop_var, N);
            
            // Gaussian noise marginal likelihood
            vector[N] mu_noise = rep_vector(theta_pop_mean, N);
            vector[N] sigma_noise_total;
            for (i in 1:N) {
                sigma_noise_total[i] = sqrt(sigma_noise^2 + theta_pop_var + se[k][i]^2);
            }
            real ll_noise = normal_lpdf(y[k] | mu_noise, sigma_noise_total);
            
            // Posterior probabilities (Bayes rule)
            real log_p_ou_post = log(1 - p_noise) + ll_ou;
            real log_p_noise_post = log(p_noise) + ll_noise;
            real log_normalizer = log_sum_exp([log_p_ou_post, log_p_noise_post]);
            
            trait_ou_prob[k] = exp(log_p_ou_post - log_normalizer);
            trait_noise_prob[k] = exp(log_p_noise_post - log_normalizer);
            
            // OU process theta posterior
            {
                matrix[N, N] vcv_ou = build_ou_cov_matrix(vcv_bm, alpha, sigma_x, N);
                matrix[N, N] Sigma_obs = vcv_ou;
                for (i in 1:N) {
                    Sigma_obs[i, i] += se[k][i]^2;
                }
                
                vector[N] mu_theta_coeff;
                for (i in 1:N) {
                    mu_theta_coeff[i] = 1 - exp(-alpha * vcv_bm[i, i]);
                }
                
                matrix[N, N] Sigma_inv = inverse(Sigma_obs);
                real posterior_precision = 1.0/theta_pop_var + mu_theta_coeff' * Sigma_inv * mu_theta_coeff;
                real posterior_var = 1.0 / posterior_precision;
                real posterior_mean = posterior_var * (theta_pop_mean/theta_pop_var + 
                                                      mu_theta_coeff' * Sigma_inv * y[k]);
                
                theta_ou_posterior_mean[k] = posterior_mean;
                theta_ou_posterior_sd[k] = sqrt(posterior_var);
            }
            
            // Gaussian noise theta posterior
            {
                real sum_precision = 0;
                real sum_weighted_y = 0;
                
                for (i in 1:N) {
                    real obs_precision = 1.0 / (sigma_noise^2 + se[k][i]^2);
                    sum_precision += obs_precision;
                    sum_weighted_y += obs_precision * y[k][i];
                }
                
                real posterior_precision = 1.0/theta_pop_var + sum_precision;
                real posterior_var = 1.0 / posterior_precision;
                real posterior_mean = posterior_var * (theta_pop_mean/theta_pop_var + sum_weighted_y);
                
                theta_noise_posterior_mean[k] = posterior_mean;
                theta_noise_posterior_sd[k] = sqrt(posterior_var);
            }
            
            // Mixture estimates (probability-weighted)
            theta_mixture_mean[k] = trait_ou_prob[k] * theta_ou_posterior_mean[k] + 
                                   trait_noise_prob[k] * theta_noise_posterior_mean[k];
            theta_mixture_sd[k] = sqrt(trait_ou_prob[k] * (theta_ou_posterior_sd[k]^2 + theta_ou_posterior_mean[k]^2) +
                                      trait_noise_prob[k] * (theta_noise_posterior_sd[k]^2 + theta_noise_posterior_mean[k]^2) -
                                      theta_mixture_mean[k]^2);
        }
    }
    
    // Summary statistics
    real mean_noise_prob = mean(trait_noise_prob);
    int n_likely_noise = 0;
    real max_noise_prob = max(trait_noise_prob);
    
    // Count traits with noise probability > 0.5
    for (k in 1:K) {
        if (trait_noise_prob[k] > 0.5) {
            n_likely_noise += 1;
        }
    }
}"

# Phylogenetic mean function with error handling
compute_phylo_mean <- function(trait_values, tree) {
    tryCatch({
        vcov_matrix <- ape::vcv(tree)
        # Check for singularity
        det_vcov <- det(vcov_matrix)
        if (abs(det_vcov) < 1e-10) {
            warning("VCV matrix near singular, using arithmetic mean")
            return(mean(trait_values, na.rm = TRUE))
        }
        vcov_inv <- solve(vcov_matrix)
        ones <- rep(1, length(trait_values))
        as.numeric((t(ones) %*% vcov_inv %*% trait_values) / (t(ones) %*% vcov_inv %*% ones))
    }, error = function(e) {
        warning("Phylogenetic mean calculation failed, using arithmetic mean: ", e$message)
        return(mean(trait_values, na.rm = TRUE))
    })
}

# Load and process data files
log_msg("=== LOADING DATA FILES ===")

# Load phylogenetic tree
log_msg(paste("Loading phylogeny from:", args$newick_file))
if (!file.exists(args$newick_file)) {
    stop("Phylogeny file not found: ", args$newick_file)
}
tree <- read.tree(args$newick_file)
log_msg(paste("Tree loaded with", Ntip(tree), "tips"))

# Normalize tree to height 1 for numerical stability
original_height <- max(node.depth.edgelength(tree))
log_msg(paste("Original tree height:", round(original_height, 6)))

tree$edge.length <- tree$edge.length / original_height
normalized_height <- max(node.depth.edgelength(tree))
log_msg(paste("Normalized tree height:", round(normalized_height, 6)))

if (abs(normalized_height - 1.0) > 1e-10) {
    warning("Tree normalization may have failed - height is ", normalized_height)
}

# Function to load CSV with special header format
load_trait_csv <- function(file_path, file_type = "trait") {
    log_msg(paste("Loading", file_type, "data from:", file_path))
    
    if (!file.exists(file_path)) {
        stop(file_type, " data file not found: ", file_path)
    }
    
    # Read first few lines to check format
    first_line <- readLines(file_path, n = 1)
    log_msg(paste("First line of", file_type, "file:", first_line))
    
    # Split first line to analyze structure
    first_line_clean <- trimws(first_line)
    if (grepl(",", first_line_clean)) {
        # CSV format - split by comma
        first_parts <- strsplit(first_line_clean, ",")[[1]]
    } else {
        # Space/tab delimited - split by whitespace
        first_parts <- strsplit(first_line_clean, "\\s+")[[1]]
    }
    
    # Remove # symbols for analysis
    first_parts_clean <- gsub("^#", "", first_parts)
    
    # Check if first line contains only numbers (dimension info)
    is_dimension_line <- length(first_parts_clean) == 2 && 
                        all(grepl("^\\d+$", trimws(first_parts_clean)))
    
    if (is_dimension_line) {
        # Parse dimensions from first line
        n_species <- as.numeric(trimws(first_parts_clean[1]))
        n_traits <- as.numeric(trimws(first_parts_clean[2]))
        
        log_msg(paste("Detected dimension header:", n_species, "species,", n_traits, "traits"))
        
        # Read the actual data starting from line 2
        if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
            data_matrix <- read.csv(file_path, skip = 1, header = FALSE, stringsAsFactors = FALSE)
        } else {
            data_matrix <- read.table(file_path, skip = 1, header = FALSE, stringsAsFactors = FALSE)
        }
        
        # Validate dimensions
        if (nrow(data_matrix) != n_species) {
            warning("Number of data rows (", nrow(data_matrix), ") doesn't match header (", n_species, ")")
        }
        if (ncol(data_matrix) != (n_traits + 1)) {  # +1 for species column
            warning("Number of data columns (", ncol(data_matrix), ") doesn't match header (", n_traits + 1, ")")
        }
        
        # Create column names
        col_names <- c("species", paste0("trait", 1:n_traits))
        colnames(data_matrix) <- col_names[1:ncol(data_matrix)]
        
        log_msg(paste("Generated column names:", paste(colnames(data_matrix), collapse = ", ")))
        
    } else {
        # Standard CSV with header row
        log_msg("Standard CSV format detected with header row")
        if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
            data_matrix <- read.csv(file_path, stringsAsFactors = FALSE)
        } else {
            data_matrix <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
        }
    }
    
    return(data_matrix)
}

# Load trait data
trait_data <- load_trait_csv(args$trait_file, "trait")

# Load measurement error data  
error_data <- load_trait_csv(args$error_file, "error")

log_msg(paste("Trait data dimensions:", nrow(trait_data), "x", ncol(trait_data)))
log_msg(paste("Error data dimensions:", nrow(error_data), "x", ncol(error_data)))

# Validate data structure
if (!args$species_col %in% colnames(trait_data)) {
    # If species_col not found, assume first column is species
    if (args$species_col == "species" && ncol(trait_data) > 1) {
        log_msg("Species column 'species' not found, assuming first column contains species names")
        colnames(trait_data)[1] <- "species"
        args$species_col <- "species"
    } else {
        stop("Species column '", args$species_col, "' not found in trait data")
    }
}

if (!args$species_col %in% colnames(error_data)) {
    # If species_col not found, assume first column is species  
    if (args$species_col == "species" && ncol(error_data) > 1) {
        log_msg("Species column 'species' not found, assuming first column contains species names")
        colnames(error_data)[1] <- "species"
    } else {
        stop("Species column '", args$species_col, "' not found in error data")
    }
}

# Extract species names
trait_species <- trait_data[[args$species_col]]
error_species <- error_data[[args$species_col]]

# Check species consistency
if (!identical(sort(trait_species), sort(error_species))) {
    stop("Species names do not match between trait and error data files")
}

# Match species with tree
tree_species <- tree$tip.label
missing_in_tree <- setdiff(trait_species, tree_species)
missing_in_data <- setdiff(tree_species, trait_species)

if (length(missing_in_tree) > 0) {
    log_msg(paste("Warning: Species in data but not in tree:", paste(missing_in_tree, collapse=", ")))
}
if (length(missing_in_data) > 0) {
    log_msg(paste("Warning: Species in tree but not in data:", paste(missing_in_data, collapse=", ")))
}

# Keep only species present in both tree and data
common_species <- intersect(tree_species, trait_species)
log_msg(paste("Using", length(common_species), "species present in both tree and data"))

# Subset tree and data
tree_subset <- keep.tip(tree, common_species)
trait_subset <- trait_data[trait_data[[args$species_col]] %in% common_species, ]
error_subset <- error_data[error_data[[args$species_col]] %in% common_species, ]

# Reorder data to match tree tip order
trait_subset <- trait_subset[match(tree_subset$tip.label, trait_subset[[args$species_col]]), ]
error_subset <- error_subset[match(tree_subset$tip.label, error_subset[[args$species_col]]), ]

# Identify trait columns (exclude species column and skip columns)
all_cols <- colnames(trait_data)
trait_cols <- setdiff(all_cols, c(args$species_col, skip_cols))
log_msg(paste("Found", length(trait_cols), "trait columns:", paste(trait_cols, collapse=", ")))

# Extract trait and error matrices
trait_matrix <- as.matrix(trait_subset[, trait_cols, drop = FALSE])
error_matrix <- as.matrix(error_subset[, trait_cols, drop = FALSE])

# Check for missing values
if (any(is.na(trait_matrix))) {
    stop("Missing values found in trait data - please handle before analysis")
}
if (any(is.na(error_matrix))) {
    stop("Missing values found in error data - please handle before analysis")
}

# Check for non-positive errors
if (any(error_matrix <= 0)) {
    stop("Non-positive measurement errors found - all errors must be > 0")
}

log_msg(paste("Final data dimensions:", nrow(trait_matrix), "species x", ncol(trait_matrix), "traits"))

# Prepare data for Stan
log_msg("=== PREPARING DATA FOR STAN ===")

n_tips <- nrow(trait_matrix)
n_traits <- ncol(trait_matrix)

# Convert to lists for Stan
y_list <- lapply(1:n_traits, function(k) trait_matrix[, k])
se_list <- lapply(1:n_traits, function(k) error_matrix[, k])

# Build phylogenetic covariance matrix
vcv_bm <- vcv(tree_subset)
tree_height <- max(node.depth.edgelength(tree_subset))

# Compute initial theta estimates using phylogenetic means
log_msg("Computing initial theta estimates using phylogenetic means...")

trait_phylo_means <- sapply(1:n_traits, function(k) {
    trait_values <- trait_matrix[, k]
    compute_phylo_mean(trait_values, tree_subset)
})

theta_init_mean <- mean(trait_phylo_means, na.rm = TRUE)
theta_init_sd <- max(sd(trait_phylo_means, na.rm = TRUE), 0.1)

log_msg(paste("Initial theta estimates: mean =", round(theta_init_mean, 3), ", sd =", round(theta_init_sd, 3)))
log_msg(paste("Tree height:", round(tree_height, 3)))

# Create Stan data list
stan_data <- list(
    N = n_tips, 
    K = n_traits,
    y = y_list, 
    vcv_bm = vcv_bm, 
    se = se_list,
    tree_height = tree_height,
    theta_init_mean = theta_init_mean,
    theta_init_sd = theta_init_sd
)

# Compile model
log_msg("=== COMPILING STAN MODEL ===")

compile_or_load_model <- function(recompile = FALSE) {
    model_code <- stan_model_code
    model_hash <- digest::digest(model_code, algo = "md5")
    compiled_model_file <- paste0("compiled_model_ou_gaussian_poptheta_", substr(model_hash, 1, 8), ".rds")
    
    log_msg(paste("Model cache file:", compiled_model_file))
    
    if (!recompile && file.exists(compiled_model_file)) {
        log_msg("Found cached model file, attempting to load...")
        tryCatch({
            model <- readRDS(compiled_model_file)
            if (!inherits(model, "stanmodel")) {
                stop("Cached object is not a stanmodel")
            }
            log_msg("SUCCESS: Cached model loaded - skipping compilation")
            return(model)
        }, error = function(e) {
            log_msg(paste("FAILED to load cached model:", e$message))
            log_msg("Will recompile...")
            if (file.exists(compiled_model_file)) {
                file.remove(compiled_model_file)
            }
        })
    }
    
    # Compile new model
    log_msg("Compiling Stan model (60-90 seconds)...")
    start_time <- Sys.time()
    
    tryCatch({
        model <- stan_model(model_code = model_code, verbose = FALSE, auto_write = FALSE)
        compile_time <- difftime(Sys.time(), start_time, units="secs")
        log_msg(paste("Model compiled successfully in", round(compile_time, 2), "seconds"))
        
        tryCatch({
            saveRDS(model, compiled_model_file)
            log_msg(paste("Model cached to:", compiled_model_file))
        }, error = function(e) {
            log_msg(paste("Warning: Failed to cache model:", e$message))
        })
        
        return(model)
        
    }, error = function(e) {
        log_msg(paste("COMPILATION FAILED:", e$message))
        stop("Stan model compilation failed: ", e$message)
    })
}

model <- compile_or_load_model(args$recompile_model)

# Fit model
log_msg("=== FITTING MODEL ===")
log_msg(paste("Starting MCMC with", args$chains, "chains,", args$iter, "iterations"))

start_time <- Sys.time()

# Simple initialization based on data scale
all_trait_values <- as.vector(trait_matrix)
y_var <- var(all_trait_values)
y_sd <- sqrt(y_var)

control_settings <- list(
    adapt_delta = 0.85,
    max_treedepth = 10,
    stepsize = 0.9,
    metric = "diag_e"
)

fit <- sampling(
    model, data = stan_data,
    chains = args$chains, 
    iter = args$iter, 
    warmup = args$warmup,
    thin = args$thin,
    cores = args$cores,
    verbose = TRUE,
    refresh = max(25, args$iter %/% 16),
    control = control_settings,
    algorithm = "NUTS",
    open_progress = FALSE,
    init = function() {
        list(
            alpha = 2.0 / tree_height,
            sigma_x = y_sd * 0.8,
            sigma_noise = y_sd * 0.3,
            p_noise = 0.3,
            theta_pop_mean = theta_init_mean,
            theta_pop_sd = theta_init_sd
        )
    }
)

mcmc_time <- difftime(Sys.time(), start_time, units="secs")
log_msg(paste("MCMC completed in", round(mcmc_time, 2), "seconds"))

# Extract results and check convergence
log_msg("=== EXTRACTING RESULTS ===")

# Check if sampling was successful
if (inherits(fit, "stanfit")) {
    # Try to extract posterior first to check if samples exist
    tryCatch({
        posterior <- rstan::extract(fit)
        
        # Check if we have valid results
        if (is.null(posterior) || length(posterior) == 0) {
            stop("Failed to extract posterior samples")
        }
        
        # Check if we have the expected parameters
        expected_params <- c("alpha_out", "sigma_x_out", "p_noise_out", "theta_pop_mean_out", "theta_pop_sd_out")
        missing_params <- setdiff(expected_params, names(posterior))
        
        if (length(missing_params) > 0) {
            warning("Missing expected parameters: ", paste(missing_params, collapse=", "))
        }
        
        log_msg("Successfully extracted posterior samples")
        log_msg(paste("Posterior contains", length(posterior), "parameter groups"))
        log_msg(paste("Sample dimensions for alpha_out:", paste(dim(posterior$alpha_out), collapse=" x ")))
        
        # Extract summary statistics
        summary_stats <- summary(fit)$summary
        log_msg("Successfully extracted summary statistics")
        
    }, error = function(e) {
        log_msg(paste("ERROR extracting results:", e$message))
        stop("Failed to extract results from Stan fit: ", e$message)
    })
    
} else {
    stop("Stan fitting returned invalid object")
}

# Convergence diagnostics with error handling
tryCatch({
    max_rhat <- max(summary_stats[, "Rhat"], na.rm = TRUE)
    min_ess <- min(summary_stats[, "n_eff"], na.rm = TRUE)
    main_params <- c("alpha", "sigma_x", "sigma_noise", "p_noise", "theta_pop_mean", "theta_pop_sd")
    
    # Check if main parameters exist in summary
    available_params <- intersect(main_params, rownames(summary_stats))
    if (length(available_params) == 0) {
        warning("No main parameters found in summary statistics")
        main_max_rhat <- NA
    } else {
        main_params_rhat <- summary_stats[available_params, "Rhat"]
        main_max_rhat <- max(main_params_rhat, na.rm = TRUE)
    }
    
}, error = function(e) {
    log_msg(paste("ERROR computing convergence diagnostics:", e$message))
    max_rhat <- NA
    min_ess <- NA
    main_max_rhat <- NA
})

# Divergences with error handling
tryCatch({
    divergences <- get_num_divergent(fit)
    n_divergences <- sum(divergences)
}, error = function(e) {
    log_msg(paste("ERROR getting divergences:", e$message))
    divergences <- rep(0, args$chains)
    n_divergences <- 0
})

log_msg(paste("Convergence diagnostics:"))
log_msg(paste("  Max Rhat =", ifelse(is.na(max_rhat), "NA", round(max_rhat, 3))))
log_msg(paste("  Min ESS =", ifelse(is.na(min_ess), "NA", round(min_ess, 0))))
log_msg(paste("  Main params Max Rhat =", ifelse(is.na(main_max_rhat), "NA", round(main_max_rhat, 3))))
log_msg(paste("  Divergences =", n_divergences))

# Parameter estimates with error handling
tryCatch({
    p_noise_mean <- mean(posterior$p_noise_out)
    n_likely_noise_mean <- mean(posterior$n_likely_noise)
    alpha_mean <- mean(posterior$alpha_out)
    sigma_x_mean <- mean(posterior$sigma_x_out)
    theta_pop_mean_est <- mean(posterior$theta_pop_mean_out)
    theta_pop_sd_est <- mean(posterior$theta_pop_sd_out)
    
    log_msg(paste("Parameter estimates:"))
    log_msg(paste("  alpha =", round(alpha_mean, 3)))
    log_msg(paste("  sigma_x =", round(sigma_x_mean, 3)))
    log_msg(paste("  p_noise =", round(p_noise_mean, 3)))
    log_msg(paste("  theta_pop_mean =", round(theta_pop_mean_est, 3)))
    log_msg(paste("  theta_pop_sd =", round(theta_pop_sd_est, 3)))
    log_msg(paste("  Expected # noise traits =", round(n_likely_noise_mean, 1)))
    
}, error = function(e) {
    log_msg(paste("ERROR computing parameter estimates:", e$message))
    # Set default values
    p_noise_mean <- NA
    n_likely_noise_mean <- NA
    alpha_mean <- NA
    sigma_x_mean <- NA
    theta_pop_mean_est <- NA
    theta_pop_sd_est <- NA
})

# Save results
log_msg("=== SAVING RESULTS ===")

output_file <- paste0(args$output_prefix, "_bayesian_fit_ou_gaussian_poptheta.rds")

result <- list(
    posterior = posterior,
    summary = summary_stats,
    stan_data = stan_data,
    trait_data = list(
        trait_matrix = trait_matrix,
        error_matrix = error_matrix,
        species_names = tree_subset$tip.label,
        trait_names = trait_cols
    ),
    tree = tree_subset,
    theta_init_info = list(
        mean = theta_init_mean,
        sd = theta_init_sd,
        trait_phylo_means = trait_phylo_means
    ),
            timing = list(mcmc = as.numeric(mcmc_time)),
    diagnostics = list(
        max_rhat = max_rhat, 
        min_ess = min_ess, 
        divergences = sum(divergences),
        main_max_rhat = main_max_rhat
    ),
    mixture_summary = list(
        p_noise_mean = p_noise_mean,
        n_likely_noise_mean = n_likely_noise_mean,
        alpha_mean = alpha_mean,
        sigma_x_mean = sigma_x_mean,
        theta_pop_mean = theta_pop_mean_est,
        theta_pop_sd = theta_pop_sd_est
    ),
    settings = list(
        mixture_model = "OU-Gaussian-PopTheta",
        n_parameters = 6,
        thinning = args$thin,
        parameterization = "population_theta",
        input_files = list(
            trait_file = args$trait_file,
            error_file = args$error_file,
            newick_file = args$newick_file
        ),
        tree_scaling = list(
            original_height = original_height,
            normalized_height = normalized_height,
            scaling_factor = original_height
        )
    ),
    fit = fit
)

saveRDS(result, output_file)
log_msg(paste("Results saved to:", output_file))

# Create summary CSV - only if we have valid results
if (!is.na(alpha_mean)) {
    summary_csv_file <- paste0(args$output_prefix, "_summary.csv")
    
    summary_data <- data.frame(
        dataset = args$output_prefix,
        n_traits = n_traits,
        n_species = n_tips,
        mcmc_time_seconds = as.numeric(mcmc_time),
        original_tree_height = original_height,
        # Initial theta estimates
        theta_init_mean = theta_init_mean,
        theta_init_sd = theta_init_sd,
        # Estimated population theta parameters
        theta_pop_mean = theta_pop_mean_est,
        theta_pop_mean_sd = ifelse(is.null(posterior$theta_pop_mean_out), NA, sd(posterior$theta_pop_mean_out)),
        theta_pop_sd = theta_pop_sd_est,
        theta_pop_sd_sd = ifelse(is.null(posterior$theta_pop_sd_out), NA, sd(posterior$theta_pop_sd_out)),
        # Main parameters
        alpha_mean = alpha_mean,
        alpha_sd = ifelse(is.null(posterior$alpha_out), NA, sd(posterior$alpha_out)),
        sigma_x_mean = sigma_x_mean,
        sigma_x_sd = ifelse(is.null(posterior$sigma_x_out), NA, sd(posterior$sigma_x_out)),
        sigma_noise_mean = ifelse(is.null(posterior$sigma_noise_out), NA, mean(posterior$sigma_noise_out)),
        sigma_noise_sd = ifelse(is.null(posterior$sigma_noise_out), NA, sd(posterior$sigma_noise_out)),
        p_noise_mean = p_noise_mean,
        p_noise_sd = ifelse(is.null(posterior$p_noise_out), NA, sd(posterior$p_noise_out)),
        # Derived quantities
        sigma_sq_mean = ifelse(is.null(posterior$sigma_sq), NA, mean(posterior$sigma_sq)),
        half_life_mean = ifelse(is.null(posterior$phylogenetic_half_life), NA, mean(posterior$phylogenetic_half_life)),
        stationary_var_mean = ifelse(is.null(posterior$stationary_variance), NA, mean(posterior$stationary_variance)),
        phi_derived_mean = ifelse(is.null(posterior$phi), NA, mean(posterior$phi)),
        n_likely_noise_mean = n_likely_noise_mean,
        max_noise_prob_mean = ifelse(is.null(posterior$max_noise_prob), NA, mean(posterior$max_noise_prob)),
        # Convergence diagnostics
        max_rhat = ifelse(is.na(max_rhat), -999, max_rhat),
        min_ess = ifelse(is.na(min_ess), -999, min_ess),
        main_max_rhat = ifelse(is.na(main_max_rhat), -999, main_max_rhat),
        divergences = n_divergences,
        mixture_model = "OU-Gaussian-PopTheta",
        parameterization = "population_theta"
    )
    
    write.csv(summary_data, summary_csv_file, row.names = FALSE)
    log_msg(paste("Summary saved to:", summary_csv_file))
} else {
    log_msg("WARNING: Skipping summary CSV creation due to failed sampling")
}

# Create trait-level results CSV - only if we have valid results
if (!is.na(alpha_mean) && !is.null(posterior$trait_noise_prob)) {
    trait_results_file <- paste0(args$output_prefix, "_trait_results.csv")
    
    trait_results <- data.frame(
        trait_name = trait_cols,
        trait_noise_prob_mean = colMeans(posterior$trait_noise_prob),
        trait_noise_prob_sd = apply(posterior$trait_noise_prob, 2, sd),
        trait_ou_prob_mean = colMeans(posterior$trait_ou_prob),
        trait_ou_prob_sd = apply(posterior$trait_ou_prob, 2, sd),
        theta_mixture_mean = colMeans(posterior$theta_mixture_mean),
        theta_mixture_sd = colMeans(posterior$theta_mixture_sd),
        theta_ou_posterior_mean = colMeans(posterior$theta_ou_posterior_mean),
        theta_ou_posterior_sd = colMeans(posterior$theta_ou_posterior_sd),
        theta_noise_posterior_mean = colMeans(posterior$theta_noise_posterior_mean),
        theta_noise_posterior_sd = colMeans(posterior$theta_noise_posterior_sd),
        likely_noise = colMeans(posterior$trait_noise_prob) > 0.5
    )
    
    write.csv(trait_results, trait_results_file, row.names = FALSE)
    log_msg(paste("Trait-level results saved to:", trait_results_file))
} else {
    log_msg("WARNING: Skipping trait results CSV creation due to failed sampling")
}

# Print detailed results summary
log_msg("=== DETAILED RESULTS SUMMARY ===")

# Parameter interpretation
half_life_mean <- mean(posterior$phylogenetic_half_life)
phi_mean <- mean(posterior$phi)

log_msg(paste("EVOLUTIONARY PARAMETERS:"))
log_msg(paste("  Selection strength (alpha):", round(alpha_mean, 3), "± ", round(sd(posterior$alpha_out), 3)))
log_msg(paste("  Diffusion rate (sigma_x):", round(sigma_x_mean, 3), "± ", round(sd(posterior$sigma_x_out), 3)))
log_msg(paste("  Phylogenetic half-life:", round(half_life_mean, 3), "time units"))
log_msg(paste("  Phylogenetic signal (phi):", round(phi_mean, 3)))

log_msg(paste("POPULATION THETA PARAMETERS:"))
log_msg(paste("  Population mean:", round(theta_pop_mean_est, 3), "± ", round(sd(posterior$theta_pop_mean_out), 3)))
log_msg(paste("  Population SD:", round(theta_pop_sd_est, 3), "± ", round(sd(posterior$theta_pop_sd_out), 3)))

log_msg(paste("MIXTURE MODEL RESULTS:"))
log_msg(paste("  Noise probability:", round(p_noise_mean, 3), "± ", round(sd(posterior$p_noise_out), 3)))
log_msg(paste("  Expected # noise traits:", round(n_likely_noise_mean, 1), "out of", n_traits))

# Trait classification
n_likely_noise_actual <- sum(colMeans(posterior$trait_noise_prob) > 0.5)
log_msg(paste("  Traits classified as noise:", n_likely_noise_actual))

if (n_likely_noise_actual > 0) {
    noise_traits <- trait_cols[colMeans(posterior$trait_noise_prob) > 0.5]
    log_msg(paste("  Likely noise traits:", paste(noise_traits, collapse=", ")))
}

# Convergence assessment
log_msg(paste("CONVERGENCE ASSESSMENT:"))
if (main_max_rhat < 1.02) {
    log_msg("  EXCELLENT: Very good convergence achieved")
} else if (main_max_rhat < 1.05) {
    log_msg("  GOOD: Acceptable convergence")
} else {
    log_msg("  WARNING: Poor convergence - consider increasing iterations")
}

if (sum(divergences) == 0) {
    log_msg("  EXCELLENT: No divergences")
} else if (sum(divergences) < 5) {
    log_msg("  GOOD: Very few divergences")
} else {
    log_msg("  WARNING: Many divergences - consider increasing adapt_delta")
}

# Biological interpretation
log_msg(paste("BIOLOGICAL INTERPRETATION:"))
if (half_life_mean < 0.1 * tree_height) {
    log_msg("  -> Very rapid evolution (weak phylogenetic constraint)")
} else if (half_life_mean < 0.5 * tree_height) {
    log_msg("  -> Moderate evolutionary rate")
} else {
    log_msg("  -> Slow evolution (strong phylogenetic constraint)")
}

if (phi_mean < 0.1) {
    log_msg("  -> Weak phylogenetic signal")
} else if (phi_mean < 0.5) {
    log_msg("  -> Moderate phylogenetic signal")
} else {
    log_msg("  -> Strong phylogenetic signal")
}

# Output file summary
log_msg("=== OUTPUT FILES CREATED ===")
log_msg(paste("1. Main results (RDS):", output_file))
if (!is.na(alpha_mean)) {
    log_msg(paste("2. Summary statistics (CSV):", summary_csv_file))
    if (!is.null(posterior$trait_noise_prob)) {
        log_msg(paste("3. Trait-level results (CSV):", trait_results_file))
    }
}
log_msg(paste("4. Log file:", LOG_FILE))

# Usage instructions
log_msg("=== USAGE INSTRUCTIONS ===")
log_msg("To load results in R:")
log_msg(paste("  result <- readRDS('", output_file, "')", sep=""))

if (!is.na(alpha_mean)) {
    log_msg("  posterior <- result$posterior")
    log_msg("  summary_stats <- result$summary")
    log_msg("")
    log_msg("Key posterior samples:")
    log_msg("  result$posterior$alpha_out           # Selection strength")
    log_msg("  result$posterior$sigma_x_out         # Diffusion rate")
    log_msg("  result$posterior$theta_pop_mean_out  # Population theta mean")
    log_msg("  result$posterior$theta_pop_sd_out    # Population theta SD")
    if (!is.null(posterior$trait_noise_prob)) {
        log_msg("  result$posterior$trait_noise_prob    # Trait noise probabilities")
        log_msg("  result$posterior$theta_mixture_mean  # Best theta estimates per trait")
    }
    log_msg("")
    log_msg("Tree scaling information:")
    log_msg("  result$tree_scaling$original_height  # Original tree height")
    log_msg("  result$tree_scaling$scaling_factor   # Factor used for normalization")
    log_msg("")
    log_msg("To convert results back to original time scale:")
    log_msg("  original_alpha <- result$posterior$alpha_out / result$tree_scaling$scaling_factor")
    log_msg("  original_half_life <- result$posterior$phylogenetic_half_life * result$tree_scaling$scaling_factor")
} else {
    log_msg("  # Note: Sampling failed, limited results available")
    log_msg("  result$diagnostics  # Check convergence issues")
    log_msg("  result$stan_data    # Check input data")
}

log_msg("=== REAL DATA POPULATION THETA SCRIPT COMPLETED ===")

# Final cleanup
gc()