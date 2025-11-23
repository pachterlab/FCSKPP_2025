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
LOG_FILE <- paste0("bayesian_ou_gaussian_mixture_poptheta_", Sys.Date(), ".log")
log_msg <- function(msg) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    entry <- paste(timestamp, "-", msg, "\n")
    cat(entry)
    cat(entry, file = LOG_FILE, append = TRUE)
}

log_msg("=== POPULATION THETA SCRIPT STARTED ===")

# Parse arguments
option_list <- list(
    make_option(c("-i", "--input_dir"), default="simulations"),
    make_option(c("--chains"), type="integer", default=NULL),
    make_option(c("--iter"), type="integer", default=1200),
    make_option(c("--warmup"), type="integer", default=600),
    make_option(c("--cores"), type="integer", default=parallel::detectCores()),
    make_option(c("--max_files"), type="integer", default=NULL),
    make_option(c("--parallel_files"), type="logical", default=TRUE),
    make_option(c("--recompile_model"), type="logical", default=FALSE),
    make_option(c("--thin"), type="integer", default=1),
    make_option(c("--log_every"), type="integer", default=50)
)
args <- parse_args(OptionParser(option_list=option_list))

# Optimize chain/core allocation
if (is.null(args$chains)) {
    args$chains <- min(4, args$cores)
}

# Calculate parallelization strategy
total_cores <- args$cores
files_per_batch <- max(1, floor(total_cores / args$chains))
cores_per_file <- min(args$chains, args$cores)

log_msg(paste("=== PARALLELIZATION STRATEGY ==="))
log_msg(paste("Total cores available:", total_cores))
log_msg(paste("Chains per file:", args$chains))
log_msg(paste("Files in parallel:", files_per_batch))
log_msg(paste("Cores per file:", cores_per_file))
log_msg(paste("Thinning:", args$thin))
log_msg(paste("Log every:", args$log_every, "iterations"))
log_msg(paste("Input dir:", args$input_dir))

# Find files
log_msg("Finding simulation files...")
start_time <- Sys.time()
sim_files <- list.files(args$input_dir, pattern="^sim_.*\\.rds$", full.names=TRUE)
if (!is.null(args$max_files)) sim_files <- sim_files[1:min(args$max_files, length(sim_files))]
log_msg(paste("Found", length(sim_files), "files in", round(difftime(Sys.time(), start_time, units="secs"), 2), "seconds"))

# Modified Stan model with population theta parameters (only 2 extra parameters)
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
    // Now uses population theta parameters instead of fixed prior
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
        // For Gaussian noise with theta ~ N(theta_pop_mean, theta_pop_var):
        // Marginal is N(theta_pop_mean, sigma_noise^2 + theta_pop_var + se^2)
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
    
    // NEW: Population theta parameters (only 2 additional parameters!)
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
    
    // NEW: Priors on population theta parameters
    theta_pop_mean ~ normal(theta_init_mean, theta_init_sd * 2);
    theta_pop_sd ~ exponential(2.0);  // Regularizing prior
    
    // Mixture likelihood using population theta parameters
    // (theta is still integrated out analytically for each trait)
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
    
    // NEW: Population theta outputs
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
            
            // OU process theta posterior (using population parameters)
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
            
            // Gaussian noise theta posterior (using population parameters)
            {
                real sum_y = sum(y[k]);
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

# Process single file function - completely self-contained
process_single_file <- function(sim_file_path, chains, iter, warmup, thin, cores_per_file, log_every, compiled_model) {
    # Load required libraries in worker
    library(rstan)
    library(ape)
    library(tools)
    
    # Set up file-specific logging
    file_log <- file.path(dirname(sim_file_path), paste0("mcmc_poptheta_", tools::file_path_sans_ext(basename(sim_file_path)), ".log"))
    file.create(file_log)
    
    file_log_msg <- function(msg) {
        timestamp <- format(Sys.time(), "%H:%M:%S")
        entry <- paste(timestamp, "-", msg, "\n")
        cat(entry, file = file_log, append = TRUE)
    }
    
    file_log_msg("Starting OU-Gaussian mixture model with population theta parameters")
    
    # Use the pre-compiled model passed as argument
    model <- compiled_model
    
    # Load data
    start_time <- Sys.time()
    sim_data <- readRDS(sim_file_path)
    load_time <- difftime(Sys.time(), start_time, units="secs")
    file_log_msg(paste("Data loaded in", round(load_time, 2), "seconds"))
    
    # Prepare data
    start_time <- Sys.time()
    trait_list <- sim_data$trait_list
    me_values <- sim_data$me_values
    tree <- sim_data$tree
    n_traits <- length(trait_list)
    n_tips <- length(tree$tip.label)
    
    file_log_msg(paste("Data dimensions: n_traits =", n_traits, ", n_tips =", n_tips))
    file_log_msg(paste("Parameters to estimate: 6 (alpha, sigma_x, sigma_noise, p_noise, theta_pop_mean, theta_pop_sd)"))
    
    # Build matrices
    y_list <- lapply(1:n_traits, function(k) as.vector(trait_list[[k]]))
    se_list <- lapply(1:n_traits, function(k) as.vector(me_values[[k]]))
    
    vcv_bm <- vcv(tree)
    tree_height <- max(node.depth.edgelength(tree))
    
    # CHANGED: Compute initial theta estimates using phylogenetic means
    file_log_msg("Computing initial theta estimates using phylogenetic means...")
    
    # Calculate phylogenetic mean for each trait
    trait_phylo_means <- sapply(trait_list, function(trait) {
        trait_values <- as.vector(trait)
        compute_phylo_mean(trait_values, tree)
    })
    
    # These are now just initial estimates, not fixed priors
    theta_init_mean <- mean(trait_phylo_means, na.rm = TRUE)
    theta_init_sd <- max(sd(trait_phylo_means, na.rm = TRUE), 0.1)
    
    file_log_msg(paste("Initial theta estimates: mean =", round(theta_init_mean, 3), ", sd =", round(theta_init_sd, 3)))
    file_log_msg(paste("Tree height:", round(tree_height, 3)))
    
    # CHANGED: Update stan_data list
    stan_data <- list(
        N = n_tips, K = n_traits,
        y = y_list, vcv_bm = vcv_bm, se = se_list,
        tree_height = tree_height,
        theta_init_mean = theta_init_mean,  # Changed from theta_prior_mean
        theta_init_sd = theta_init_sd       # Changed from theta_prior_sd
    )
    
    prep_time <- difftime(Sys.time(), start_time, units="secs")
    file_log_msg(paste("Data prep completed in", round(prep_time, 2), "seconds"))
    
    # Fit model
    file_log_msg(paste("Starting MCMC with", chains, "chains,", iter, "iterations"))
    start_time <- Sys.time()
    
    refresh_rate <- max(25, min(75, iter %/% 16))
    
    control_settings <- list(
        adapt_delta = 0.85,
        max_treedepth = 10,
        stepsize = 0.9,
        metric = "diag_e"
    )
    
    # Simple initialization based on data scale
    all_trait_values <- unlist(y_list)
    y_var <- var(all_trait_values)
    y_sd <- sqrt(y_var)
    
    # CHANGED: Update initialization function
    fit <- sampling(
        model, data = stan_data,
        chains = chains, 
        iter = iter, 
        warmup = warmup,
        thin = thin,
        cores = cores_per_file,
        verbose = TRUE,
        refresh = refresh_rate,
        control = control_settings,
        algorithm = "NUTS",
        open_progress = FALSE,
        init = function() {
            list(
                alpha = 2.0 / tree_height,
                sigma_x = y_sd * 0.8,
                sigma_noise = y_sd * 0.3,
                p_noise = 0.3,
                theta_pop_mean = theta_init_mean,  # NEW: Initialize population theta mean
                theta_pop_sd = theta_init_sd       # NEW: Initialize population theta sd
            )
        }
    )
    
    mcmc_time <- difftime(Sys.time(), start_time, units="secs")
    file_log_msg(paste("MCMC completed in", round(mcmc_time, 2), "seconds"))
    
    # Save results
    start_time <- Sys.time()
    # CHANGED: Update output filename
    output_file <- file.path(dirname(sim_file_path), paste0("bayesian_fit_ou_gaussian_poptheta_", tools::file_path_sans_ext(basename(sim_file_path)), ".rds"))
    
    file_log_msg("Extracting posterior samples...")
    posterior <- rstan::extract(fit)
    summary_stats <- summary(fit)$summary
    
    # CHANGED: Update convergence diagnostics
    max_rhat <- max(summary_stats[, "Rhat"], na.rm = TRUE)
    min_ess <- min(summary_stats[, "n_eff"], na.rm = TRUE)
    main_params <- c("alpha", "sigma_x", "sigma_noise", "p_noise", "theta_pop_mean", "theta_pop_sd")
    main_params_rhat <- summary_stats[main_params, "Rhat"]
    main_max_rhat <- max(main_params_rhat, na.rm=TRUE)
    
    file_log_msg(paste("Initial convergence: Max Rhat =", round(max_rhat, 3), ", Min ESS =", round(min_ess, 0)))
    file_log_msg(paste("Main params Max Rhat =", round(main_max_rhat, 3)))
    
    # Extend sampling if convergence is poor
    if (main_max_rhat > 1.05 || min_ess < 200) {
        file_log_msg("Poor convergence detected - extending sampling...")
        
        additional_iter <- 500
        fit_extended <- sampling(
            model, data = stan_data,
            chains = chains,
            iter = iter + additional_iter,
            warmup = warmup,
            thin = thin,
            cores = cores_per_file,
            verbose = FALSE,
            refresh = 0,
            control = control_settings,
            init = rstan::extract(fit, permuted = FALSE, inc_warmup = FALSE)[warmup, , , drop = FALSE]
        )
        
        fit <- fit_extended
        file_log_msg("Extended sampling completed")
        
        # Recompute diagnostics
        posterior <- rstan::extract(fit)
        summary_stats <- summary(fit)$summary
        max_rhat <- max(summary_stats[, "Rhat"], na.rm = TRUE)
        min_ess <- min(summary_stats[, "n_eff"], na.rm = TRUE)
        main_params_rhat <- summary_stats[main_params, "Rhat"]
        main_max_rhat <- max(main_params_rhat, na.rm=TRUE)
    }
    
    divergences <- get_num_divergent(fit)
    
    file_log_msg(paste("Final convergence: Max Rhat =", round(max_rhat, 3), ", Min ESS =", round(min_ess, 0)))
    file_log_msg(paste("Main params Max Rhat =", round(main_max_rhat, 3)))
    file_log_msg(paste("Divergences =", sum(divergences)))
    
    # CHANGED: Update model summary logging
    p_noise_mean <- mean(posterior$p_noise_out)
    n_likely_noise_mean <- mean(posterior$n_likely_noise)
    alpha_mean <- mean(posterior$alpha_out)
    sigma_x_mean <- mean(posterior$sigma_x_out)
    theta_pop_mean_est <- mean(posterior$theta_pop_mean_out)      # NEW
    theta_pop_sd_est <- mean(posterior$theta_pop_sd_out)         # NEW
    
    file_log_msg(paste("Results: alpha =", round(alpha_mean, 3)))
    file_log_msg(paste("Results: sigma_x =", round(sigma_x_mean, 3)))
    file_log_msg(paste("Results: p_noise =", round(p_noise_mean, 3)))
    file_log_msg(paste("Results: theta_pop_mean =", round(theta_pop_mean_est, 3)))  # NEW
    file_log_msg(paste("Results: theta_pop_sd =", round(theta_pop_sd_est, 3)))      # NEW
    file_log_msg(paste("Expected # noise traits =", round(n_likely_noise_mean, 1)))
    
    # Additional diagnostics for alpha collapse
    alpha_q05 <- quantile(posterior$alpha_out, 0.05)
    alpha_q95 <- quantile(posterior$alpha_out, 0.95)
    file_log_msg(paste("Alpha 90% CI: [", round(alpha_q05, 3), ",", round(alpha_q95, 3), "]"))
    
    # Check for alpha near boundary
    if (alpha_mean < 0.1 / tree_height) {
        file_log_msg("WARNING: Alpha near lower boundary - possible identification issue")
    }
    
    # Signal strength diagnostic
    mean_phi_derived <- mean(posterior$phi)
    file_log_msg(paste("Derived phylogenetic signal (phi):", round(mean_phi_derived, 3)))
    
    # CHANGED: Update result object
    result <- list(
        posterior = posterior,
        summary = summary_stats,
        stan_data = stan_data,
        simulation_file = sim_file_path,
        true_parameters = sim_data$parameters,
        theta_init_info = list(                    # Changed from theta_prior_info
            mean = theta_init_mean,
            sd = theta_init_sd,
            trait_phylo_means = trait_phylo_means
        ),
        timing = list(prep = prep_time, mcmc = mcmc_time),
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
            theta_pop_mean = theta_pop_mean_est,     # NEW
            theta_pop_sd = theta_pop_sd_est          # NEW
        ),
        settings = list(
            mixture_model = "OU-Gaussian-PopTheta",  # Updated name
            n_parameters = 6,                        # Updated count
            thinning = thin,
            parameterization = "population_theta"    # Updated description
        ),
        fit = fit
    )
    
    saveRDS(result, output_file)
    file_log_msg(paste("Results saved to:", basename(output_file)))
    file_log_msg("File processing completed")
    
    return(output_file)
}

# CHANGED: Update logging messages
log_msg("=== COMPILING/LOADING STAN OU-GAUSSIAN MIXTURE MODEL (POPULATION THETA) ===")

# Compile model once in main process
compile_or_load_model <- function(input_dir, recompile = FALSE) {
    model_code <- stan_model_code
    model_hash <- digest::digest(model_code, algo = "md5")
    # CHANGED: Update cache filename
    compiled_model_file <- file.path(input_dir, paste0("compiled_ou_gaussian_poptheta_", substr(model_hash, 1, 8), ".rds"))
    
    log_msg(paste("Model cache file:", basename(compiled_model_file)))
    
    if (!recompile && file.exists(compiled_model_file)) {
        log_msg("Found cached model file, attempting to load...")
        tryCatch({
            file_size <- file.info(compiled_model_file)$size
            if (is.na(file_size) || file_size < 1000) {
                stop("Cached file appears corrupted (too small)")
            }
            
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
                log_msg("Removed corrupted cache file")
            }
        })
    } else if (recompile) {
        log_msg("Force recompilation requested")
    } else {
        log_msg("No cached model found")
    }
    
    # Compile new model
    log_msg("Compiling Stan OU-Gaussian mixture model with population theta parameters (60-90 seconds)...")
    start_time <- Sys.time()
    
    old_auto_write <- getOption("rstan_auto_write", FALSE)
    options(rstan_auto_write = FALSE)
    
    tryCatch({
        model <- stan_model(model_code = model_code, verbose = FALSE, auto_write = FALSE)
        compile_time <- difftime(Sys.time(), start_time, units="secs")
        log_msg(paste("Model compiled successfully in", round(compile_time, 2), "seconds"))
        
        tryCatch({
            saveRDS(model, compiled_model_file)
            log_msg(paste("Model cached to:", basename(compiled_model_file)))
        }, error = function(e) {
            log_msg(paste("Warning: Failed to cache model:", e$message))
        })
        
        return(model)
        
    }, error = function(e) {
        log_msg(paste("COMPILATION FAILED:", e$message))
        stop("Stan model compilation failed: ", e$message)
    }, finally = {
        options(rstan_auto_write = old_auto_write)
    })
}

# Compile model once to create cache
model <- compile_or_load_model(args$input_dir, args$recompile_model)

if (args$parallel_files && files_per_batch > 1 && length(sim_files) > 1) {
    log_msg("=== PROCESSING FILES IN PARALLEL ===")
    overall_start <- Sys.time()
    
    # Use mclapply and pass the compiled model to each worker
    results <- parallel::mclapply(sim_files, function(sim_file) {
        tryCatch({
            process_single_file(sim_file, args$chains, args$iter, args$warmup, 
                              args$thin, cores_per_file, args$log_every, model)
        }, error = function(e) {
            cat("Error processing", basename(sim_file), ":", e$message, "\n")
            return(NULL)
        })
    }, mc.cores = files_per_batch, mc.preschedule = FALSE)
    
    total_time <- difftime(Sys.time(), overall_start, units="mins")
    log_msg(paste("All files completed in", round(total_time, 2), "minutes"))
    
} else {
    log_msg("=== PROCESSING FILES SEQUENTIALLY ===")
    overall_start <- Sys.time()
    results <- list()
    
    for (i in 1:length(sim_files)) {
        log_msg(paste("Processing file", i, "of", length(sim_files), ":", basename(sim_files[i])))
        results[[i]] <- tryCatch({
            process_single_file(sim_files[i], args$chains, args$iter, args$warmup, 
                              args$thin, cores_per_file, args$log_every, model)
        }, error = function(e) {
            log_msg(paste("Error processing", basename(sim_files[i]), ":", e$message))
            return(NULL)
        })
    }
    
    total_time <- difftime(Sys.time(), overall_start, units="mins")
    log_msg(paste("All files completed in", round(total_time, 2), "minutes"))
}

# Create summary
successful_results <- results[!sapply(results, is.null)]
log_msg(paste("Successfully processed", length(successful_results), "out of", length(sim_files), "files"))

# CHANGED: Update summary filename
summary_file <- file.path(args$input_dir, "bayesian_fit_ou_gaussian_poptheta_summary.csv")
summary_data <- data.frame()

for (result_file in successful_results) {
    result <- readRDS(result_file)
    
    # CHANGED: Update summary CSV columns
    summary_row <- data.frame(
        file = basename(result$simulation_file),
        n_traits = result$stan_data$K,
        n_tips = result$stan_data$N,
        prep_time = as.numeric(result$timing$prep),
        mcmc_time = as.numeric(result$timing$mcmc),
        # Initial theta estimates (not priors)
        theta_init_mean = result$theta_init_info$mean,
        theta_init_sd = result$theta_init_info$sd,
        # Estimated population theta parameters
        theta_pop_mean = mean(result$posterior$theta_pop_mean_out),     # NEW
        theta_pop_mean_sd = sd(result$posterior$theta_pop_mean_out),    # NEW
        theta_pop_sd = mean(result$posterior$theta_pop_sd_out),         # NEW
        theta_pop_sd_sd = sd(result$posterior$theta_pop_sd_out),        # NEW
        # Main parameters (same as before)
        alpha_mean = mean(result$posterior$alpha_out),
        alpha_sd = sd(result$posterior$alpha_out),
        sigma_x_mean = mean(result$posterior$sigma_x_out),
        sigma_x_sd = sd(result$posterior$sigma_x_out),
        sigma_noise_mean = mean(result$posterior$sigma_noise_out),
        sigma_noise_sd = sd(result$posterior$sigma_noise_out),
        p_noise_mean = result$mixture_summary$p_noise_mean,
        p_noise_sd = sd(result$posterior$p_noise_out),
        # Derived quantities (same as before)
        sigma_sq_mean = mean(result$posterior$sigma_sq),
        half_life_mean = mean(result$posterior$phylogenetic_half_life),
        stationary_var_mean = mean(result$posterior$stationary_variance),
        phi_derived_mean = mean(result$posterior$phi),
        n_likely_noise_mean = result$mixture_summary$n_likely_noise_mean,
        max_noise_prob_mean = mean(result$posterior$max_noise_prob),
        # Convergence diagnostics (same as before)
        max_rhat = result$diagnostics$max_rhat,
        min_ess = result$diagnostics$min_ess,
        main_max_rhat = result$diagnostics$main_max_rhat,
        divergences = result$diagnostics$divergences,
        mixture_model = result$settings$mixture_model,
        parameterization = result$settings$parameterization
    )
    summary_data <- rbind(summary_data, summary_row)
}

write.csv(summary_data, summary_file, row.names = FALSE)
log_msg(paste("OU-Gaussian population theta summary saved to:", summary_file))

# Additional summary statistics
if (length(successful_results) > 0) {
    log_msg("=== OU-GAUSSIAN MIXTURE MODEL SUMMARY (POPULATION THETA) ===")
    
    # Aggregate diagnostics
    all_max_rhat <- sapply(successful_results, function(f) readRDS(f)$diagnostics$max_rhat)
    all_main_rhat <- sapply(successful_results, function(f) readRDS(f)$diagnostics$main_max_rhat)
    all_p_noise_mean <- sapply(successful_results, function(f) readRDS(f)$mixture_summary$p_noise_mean)
    all_divergences <- sapply(successful_results, function(f) readRDS(f)$diagnostics$divergences)
    all_alpha_mean <- sapply(successful_results, function(f) readRDS(f)$mixture_summary$alpha_mean)
    all_sigma_x_mean <- sapply(successful_results, function(f) readRDS(f)$mixture_summary$sigma_x_mean)
    
    log_msg(paste("Convergence summary:"))
    log_msg(paste("  Overall - Max Rhat range:", round(min(all_max_rhat, na.rm=TRUE), 3), "to", round(max(all_max_rhat, na.rm=TRUE), 3)))
    log_msg(paste("  Main params - Max Rhat range:", round(min(all_main_rhat, na.rm=TRUE), 3), "to", round(max(all_main_rhat, na.rm=TRUE), 3)))
    log_msg(paste("  Files with Rhat > 1.1:", sum(all_max_rhat > 1.1, na.rm=TRUE)))
    log_msg(paste("  Files with divergences:", sum(all_divergences > 0, na.rm=TRUE)))
    log_msg(paste("  Total divergences:", sum(all_divergences, na.rm=TRUE)))
    
    # Parameter estimates summary
    log_msg(paste("Parameter estimates summary:"))
    log_msg(paste("  Alpha (selection strength) range:", round(min(all_alpha_mean), 3), "to", round(max(all_alpha_mean), 3)))
    log_msg(paste("  Mean alpha:", round(mean(all_alpha_mean), 3)))
    log_msg(paste("  Sigma_x (diffusion rate) range:", round(min(all_sigma_x_mean), 3), "to", round(max(all_sigma_x_mean), 3)))
    log_msg(paste("  Mean sigma_x:", round(mean(all_sigma_x_mean), 3)))
    
    # CHANGED: Population theta summary
    all_theta_pop_mean <- sapply(successful_results, function(f) readRDS(f)$mixture_summary$theta_pop_mean)
    all_theta_pop_sd <- sapply(successful_results, function(f) readRDS(f)$mixture_summary$theta_pop_sd)
    
    log_msg(paste("Population theta summary:"))
    log_msg(paste("  Estimated theta population mean range:", round(min(all_theta_pop_mean), 3), "to", round(max(all_theta_pop_mean), 3)))
    log_msg(paste("  Mean estimated theta population mean:", round(mean(all_theta_pop_mean), 3)))
    log_msg(paste("  Estimated theta population SD range:", round(min(all_theta_pop_sd), 3), "to", round(max(all_theta_pop_sd), 3)))
    log_msg(paste("  Mean estimated theta population SD:", round(mean(all_theta_pop_sd), 3)))
    
    # Compare with simulation truth if available
    first_result <- readRDS(successful_results[1])
    if (!is.null(first_result$true_parameters$Theta_mean) && !is.null(first_result$true_parameters$Theta_sd)) {
        true_theta_mean <- first_result$true_parameters$Theta_mean
        true_theta_sd <- first_result$true_parameters$Theta_sd
        
        log_msg(paste("Population theta recovery:"))
        log_msg(paste("  True theta mean:", round(true_theta_mean, 3), "vs estimated mean:", round(mean(all_theta_pop_mean), 3)))
        log_msg(paste("  True theta SD:", round(true_theta_sd, 3), "vs estimated mean:", round(mean(all_theta_pop_sd), 3)))
        log_msg(paste("  Theta mean bias:", round(mean(all_theta_pop_mean) - true_theta_mean, 3)))
        log_msg(paste("  Theta SD bias:", round(mean(all_theta_pop_sd) - true_theta_sd, 3)))
    }
    
    # Mixture model insights
    log_msg(paste("Mixture model insights:"))
    log_msg(paste("  Mean noise probability across files:", round(mean(all_p_noise_mean), 3)))
    log_msg(paste("  Noise probability range:", round(min(all_p_noise_mean), 3), "to", round(max(all_p_noise_mean), 3)))
    log_msg(paste("  Files with substantial noise (p > 0.1):", sum(all_p_noise_mean > 0.1)))
    log_msg(paste("  Files with high noise (p > 0.3):", sum(all_p_noise_mean > 0.3)))
    
    # Performance summary
    all_mcmc_times <- sapply(successful_results, function(f) as.numeric(readRDS(f)$timing$mcmc))
    all_n_traits <- sapply(successful_results, function(f) readRDS(f)$stan_data$K)
    
    log_msg(paste("Performance summary:"))
    log_msg(paste("  MCMC time range:", round(min(all_mcmc_times), 1), "to", round(max(all_mcmc_times), 1), "seconds"))
    log_msg(paste("  Mean MCMC time:", round(mean(all_mcmc_times), 1), "seconds"))
    log_msg(paste("  Mean time per trait:", round(mean(all_mcmc_times / all_n_traits), 2), "seconds"))
    
    # CHANGED: Model summary
    log_msg(paste("=== POPULATION THETA MODEL FEATURES ==="))
    log_msg("  1. OU process vs. independent Gaussian noise at tips")
    log_msg("  2. Direct estimation of alpha (selection strength)")
    log_msg("  3. Direct estimation of sigma_x (diffusion rate)")
    log_msg("  4. Estimation of population theta mean and SD")
    log_msg("  5. Analytical theta integration (exact)")
    log_msg("  6. Only 6 parameters: alpha, sigma_x, sigma_noise, p_noise, theta_pop_mean, theta_pop_sd")
    log_msg("  7. No theta prior misspecification")
    
    # Interpretation guide
    log_msg("=== INTERPRETATION GUIDE ===")
    log_msg("  - alpha: OU selection strength (larger = stronger pull to optimum)")
    log_msg("  - sigma_x: OU diffusion rate (evolutionary rate parameter)")
    log_msg("  - sigma_noise: Standard deviation of Gaussian noise at tips")
    log_msg("  - p_noise: Probability that a trait is pure noise")
    log_msg("  - theta_pop_mean: Population mean of trait optima")
    log_msg("  - theta_pop_sd: Population SD of trait optima")
    log_msg("  - phi (derived): Phylogenetic signal = f(alpha, tree_height)")
    log_msg("  - half_life: Time for process to decay by 50% = ln(2)/alpha")
    log_msg("  - stationary_variance: Long-term variance = sigma_x^2/(2*alpha)")
    log_msg("  - trait_noise_prob[k] > 0.5: Trait k likely from noise process")
    
    # Advantages of population theta approach
    log_msg("=== ADVANTAGES OF POPULATION THETA APPROACH ===")
    log_msg("  1. Eliminates theta prior misspecification")
    log_msg("  2. Learns theta distribution from data")
    log_msg("  3. Maintains analytical exactness")
    log_msg("  4. Minimal computational overhead (only 2 extra parameters)")
    log_msg("  5. Should reduce systematic bias in alpha and sigma_noise")
    log_msg("  6. Direct interpretability of population theta parameters")
    
    # Convergence recommendations
    if (mean(all_max_rhat, na.rm=TRUE) < 1.02) {
        log_msg("EXCELLENT: Very good convergence achieved")
    } else if (mean(all_max_rhat, na.rm=TRUE) < 1.05) {
        log_msg("GOOD: Acceptable convergence")
    } else {
        log_msg("RECOMMENDATION: Consider increasing iterations")
    }
    
    if (mean(all_divergences, na.rm=TRUE) == 0) {
        log_msg("EXCELLENT: No divergences")
    } else if (mean(all_divergences, na.rm=TRUE) < 2) {
        log_msg("GOOD: Very few divergences")
    } else {
        log_msg("RECOMMENDATION: Increase adapt_delta")
    }
    
    # Parameter interpretation
    mean_half_life <- mean(sapply(successful_results, function(f) {
        posterior <- readRDS(f)$posterior
        mean(posterior$phylogenetic_half_life)
    }))
    
    log_msg(paste("=== BIOLOGICAL INTERPRETATION ==="))
    log_msg(paste("  Mean phylogenetic half-life:", round(mean_half_life, 3), "time units"))
    log_msg(paste("  Mean population theta mean:", round(mean(all_theta_pop_mean), 3)))
    log_msg(paste("  Mean population theta SD:", round(mean(all_theta_pop_sd), 3)))
    log_msg(paste("  This suggests traits evolve toward optima that vary with"))
    log_msg(paste("  SD =", round(mean(all_theta_pop_sd), 3), "around a mean of", round(mean(all_theta_pop_mean), 3)))
    
    if (mean_half_life < 0.1) {
        log_msg("  -> Very rapid evolution (weak phylogenetic constraint)")
    } else if (mean_half_life < 1.0) {
        log_msg("  -> Moderate evolutionary rate")
    } else {
        log_msg("  -> Slow evolution (strong phylogenetic constraint)")
    }
}

log_msg("=== POPULATION THETA SCRIPT COMPLETED ===")

# Final cleanup
gc()