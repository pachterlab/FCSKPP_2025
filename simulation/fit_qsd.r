#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(ape)
    library(optimx)
    library(numDeriv)
    library(expm)
})

LOG_FILE <- paste0("log_ml_", Sys.Date(), ".log")
log_msg <- function(msg) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    entry <- paste(timestamp, "-", msg, "\n")
    cat(entry)
    cat(entry, file = LOG_FILE, append = TRUE)
}

log_msg("=== 2D POPULATION THETA ML SCRIPT STARTED (NO ME OPTIMIZATION) ===")

option_list <- list(
    make_option(c("-i", "--input_dir"), default="sim_qsd"),
    make_option(c("--max_files"), type="integer", default=NULL),
    make_option(c("--optimizer"), default="Nelder-Mead"),
    make_option(c("--max_iter"), type="integer", default=2500),
    make_option(c("--reltol"), type="numeric", default=1e-8),
    make_option(c("--n_starts"), type="integer", default=30)
)
args <- parse_args(OptionParser(option_list=option_list))

## Change these based on the kind of input files (Simulated under model 1 vs model 2).
sim_files <- list.files(args$input_dir, pattern="^sim_2d_.*\\.rds$", full.names=TRUE)
# sim_files <- list.files(args$input_dir, pattern="^sim_model2_.*\\.rds$", full.names=TRUE)

if (!is.null(args$max_files)) sim_files <- sim_files[1:min(args$max_files, length(sim_files))]

log_msg(paste("Found", length(sim_files), "files"))
log_msg(paste("Using optimizer:", args$optimizer))
log_msg(paste("Number of starting points:", args$n_starts))

make_positive_definite <- function(M) {
    M <- (M + t(M)) / 2
    if (any(!is.finite(M))) stop("Matrix contains NaN or infinite values")
    eigenvals <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
    min_eigenval <- min(eigenvals)
    if (min_eigenval <= 1e-6) {
        regularization <- 1e-4 - min_eigenval
        diag(M) <- diag(M) + regularization
    }
    return(M)
}

# TODO: re-write this function.
# lyapunov_ab <- function(alpha, beta, M) {
#     v22 <- M[2,2] / (2 * beta)
#     v12 <- (M[1,2] + alpha * v22) / (alpha + beta)
#     v11 <- M[1,1] / (2 * alpha) + v12
#     V <- matrix(c(v11, v12, v12, v22), 2, 2)
#     return(V)
# }

# {{((a + e) s2)/(2 e (2 a + e)), (a s2)/(2 e (2 a + e))}, {(a s2)/(
#   2 e (2 a + e)), ((a + e) s2)/(2 e (2 a + e))}}

lyapunov_epsilon <- function(alpha, epsilon, sigma) {
    a <- alpha
    e <- epsilon
    s2 <- sigma**2
    
    v11 <- ((a + e)*s2)/(2*e*(2*a + e))
    v12 <- (a*s2)/(2*e*(2*a + e))
    v21 <- (a*s2)/(2*e*(2*a + e))
    v22 <-((a + e)*s2)/(2*e*(2*a + e))
    
    V <- matrix(c(v11, v12, v21, v22), 2, 2)
    
    return(V)
}

# Build 2x2 OU covariance matrix (unchanged)
build_2d_ou_cov_matrix_general <- function(vcv_bm, alpha, sigma_x, epsilon) {
    if (!is.matrix(vcv_bm)) vcv_bm <- as.matrix(vcv_bm)
    if (nrow(vcv_bm) != ncol(vcv_bm)) stop("vcv_bm must be square")
    if (!is.numeric(vcv_bm)) vcv_bm <- matrix(as.numeric(vcv_bm), nrow(vcv_bm), ncol(vcv_bm))
    
    N <- nrow(vcv_bm)
    H <- matrix(c(alpha + epsilon, -alpha, -alpha, alpha+epsilon), 2, 2, byrow = TRUE)
    Sigma <- diag(c(sigma_x, sigma_x))
    V_stat <- lyapunov_epsilon(alpha, epsilon, sigma_x)
    
    unique_times <- unique(as.vector(vcv_bm))
    exp_cache <- list()
    for (t_val in unique_times) {
        if (!is.finite(t_val) || length(t_val) == 0) next
        exp_cache[[as.character(t_val)]] <- expm(-H*t_val) 
    }
    vcv_ou <- array(0, dim = c(N, N, 2, 2))
    
    for (i in 1:N) {
        t_i <- vcv_bm[i,i]
        if (length(t_i) == 0 || !is.finite(t_i)) stop("Invalid tip time for tip i")
        exp_neg_H_ti <- exp_cache[[as.character(t_i)]]
        if (is.null(exp_neg_H_ti)) exp_neg_H_ti <- expm(-H*t_i)
        vcv_diag <- V_stat - exp_neg_H_ti %*% V_stat %*% t(exp_neg_H_ti)
        vcv_ou[i,i,,] <- vcv_diag
        
        if (i > 1) {
            for (j in 1:(i-1)) {
                t_j <- vcv_bm[j,j]
                t_mrca <- vcv_bm[i,j]
                
                if (length(t_mrca) == 0 || !is.finite(t_mrca)) stop(sprintf("Invalid MRCA time for pair (%d,%d)", i, j))
                if (t_mrca < -1e-12) stop("Negative MRCA time")
                
                exp_neg_H_tj <- exp_cache[[as.character(t_j)]]; if (is.null(exp_neg_H_tj)) exp_neg_H_tj <- expm(-H*t_j) 
                exp_neg_H_tmrca <- exp_cache[[as.character(t_mrca)]]; if (is.null(exp_neg_H_tmrca)) exp_neg_H_tmrca <- expm(-H*t_mrca) 

                I_branch <- V_stat - exp_neg_H_tmrca %*% V_stat %*% t(exp_neg_H_tmrca)
                
                dt_i <- t_i - t_mrca
                dt_j <- t_j - t_mrca
                if (dt_i < -1e-12 || dt_j < -1e-12) stop("Negative delta times")

                exp_neg_H_dt_i <- expm(-H*dt_i) 
                exp_neg_H_dt_j <- expm(-H*dt_j) 

                vcv_offdiag <- exp_neg_H_dt_i %*% I_branch %*% t(exp_neg_H_dt_j)
                vcv_ou[i,j,,] <- vcv_offdiag
                vcv_ou[j,i,,] <- t(vcv_offdiag)
            }
        }
    }
    
    return(vcv_ou)
}

# Precompute covariance structure for OU component
precompute_ou_cov_structure <- function(vcv_bm, alpha, sigma_x,
                                        theta_mu_pop_sd, epsilon) {
    N <- nrow(vcv_bm)
    
    # Build the OU covariance array
    vcv_ou <- build_2d_ou_cov_matrix_general(vcv_bm, alpha, sigma_x, epsilon)
    
    # Build the full 2N x 2N observation covariance matrix (no measurement error)
    Sigma_obs <- matrix(0, 2*N, 2*N)
    for (i in 1:N) {
        for (j in 1:N) {
            row_start <- (i-1)*2 + 1
            col_start <- (j-1)*2 + 1
            Sigma_obs[row_start:(row_start+1), col_start:(col_start+1)] <- vcv_ou[i,j,,]
        }
    }
    
    # Build coefficient matrix for theta means
    H <- matrix(c(alpha + epsilon, -alpha, -alpha, alpha + epsilon), 2, 2, byrow = TRUE)
    mu_theta_coeff <- matrix(0, 2*N, 2)
    for (i in 1:N) {
        t_i <- vcv_bm[i,i]
        exp_neg_H_ti <- expm(-H*t_i)
        mean_coeff <- diag(2) - exp_neg_H_ti
        mu_theta_coeff[(i-1)*2 + 1,] <- mean_coeff[1,]
        mu_theta_coeff[(i-1)*2 + 2,] <- mean_coeff[2,]
    }
    
    # Build theta population covariance
    theta_pop_cov_2d_raw <- matrix(c(
        theta_mu_pop_sd^2 , theta_mu_pop_sd^2,
        theta_mu_pop_sd^2, theta_mu_pop_sd^2
    ), 2, 2)
    theta_pop_cov_2d <- make_positive_definite(theta_pop_cov_2d_raw)
    
    # Build full marginal covariance
    Sigma_marginal_raw <- Sigma_obs + mu_theta_coeff %*% theta_pop_cov_2d %*% t(mu_theta_coeff)
    Sigma_marginal <- make_positive_definite(Sigma_marginal_raw)
    
    # Cholesky decomposition
    tryCatch({
        L <- chol(Sigma_marginal)
        log_det <- 2 * sum(log(diag(L)))
        
        return(list(
            L = L,
            log_det = log_det,
            mu_theta_coeff = mu_theta_coeff,
            N = N,
            success = TRUE
        ))
    }, error = function(e) {
        return(list(success = FALSE))
    })
}

# # Precompute covariance structure for noise component
# precompute_noise_cov_structure <- function(N, theta_mu_pop_sd, theta_gamma_pop_sd,
#                                           sigma_noise1, sigma_noise2) {
#     # Build theta population covariance
#     theta_cov_2d_raw <- matrix(c(
#         theta_mu_pop_sd^2 + theta_gamma_pop_sd^2, theta_gamma_pop_sd^2,
#         theta_gamma_pop_sd^2, theta_gamma_pop_sd^2
#     ), 2, 2)
#     theta_cov_2d <- make_positive_definite(theta_cov_2d_raw)
    
#     # Build noise covariance (block diagonal structure)
#     Sigma_noise <- matrix(0, 2*N, 2*N)
#     for (i in 1:N) {
#         Sigma_noise[(i-1)*2 + 1, (i-1)*2 + 1] <- theta_cov_2d[1,1] + sigma_noise1^2
#         Sigma_noise[(i-1)*2 + 2, (i-1)*2 + 2] <- theta_cov_2d[2,2] + sigma_noise2^2
#         Sigma_noise[(i-1)*2 + 1, (i-1)*2 + 2] <- theta_cov_2d[1,2]
#         Sigma_noise[(i-1)*2 + 2, (i-1)*2 + 1] <- theta_cov_2d[2,1]
#     }
#     Sigma_noise_pd <- make_positive_definite(Sigma_noise)
    
#     # Cholesky decomposition
#     tryCatch({
#         L <- chol(Sigma_noise_pd)
#         log_det <- 2 * sum(log(diag(L)))
        
#         return(list(
#             L = L,
#             log_det = log_det,
#             N = N,
#             success = TRUE
#         ))
#     }, error = function(e) {
#         return(list(success = FALSE))
#     })
# }

# Fixed noise covariance builder: shared theta across tips + per-tip independent noise
precompute_noise_cov_structure <- function(N, theta_mu_pop_sd,
                                                 sigma_noise1, sigma_noise2) {

  # theta population covariance (2x2)
  theta_cov_2d_raw <- matrix(c(
    theta_mu_pop_sd^2 , theta_mu_pop_sd^2,
    theta_mu_pop_sd^2, theta_mu_pop_sd^2
  ), 2, 2)
  theta_cov_2d <- make_positive_definite(theta_cov_2d_raw)

  # Build full covariance: every tip pair (i,j) shares theta_cov_2d,
  # plus measurement/noise variances on the diagonal blocks
  # Sigma_theta_all = J_N ⊗ theta_cov_2d  (J_N = all-ones NxN)
  Sigma_theta_all <- kronecker(matrix(1, nrow = N, ncol = N), theta_cov_2d)

  # Measurement / per-tip noise variance (2N x 2N diag)
  meas_diag <- diag(rep(c(sigma_noise1^2, sigma_noise2^2), N))

  Sigma_noise_full <- Sigma_theta_all + meas_diag

  Sigma_noise_pd <- make_positive_definite(Sigma_noise_full)

  tryCatch({
    L <- chol(Sigma_noise_pd)
    log_det <- 2 * sum(log(diag(L)))
    return(list(
      L = L,
      log_det = log_det,
      N = N,
      success = TRUE
    ))
  }, error = function(e) {
    return(list(success = FALSE))
  })
}


# Fast likelihood for single trait using precomputed structures
fast_single_trait_ou_ll <- function(y_vec, precomp_ou, theta_mu_pop_mean) {
    if (!precomp_ou$success) return(-1e10)
    
    theta_pop_mean_2d <- c(theta_mu_pop_mean , theta_mu_pop_mean)
    mu_marginal <- precomp_ou$mu_theta_coeff %*% theta_pop_mean_2d
    
    residual <- y_vec - mu_marginal
    quad_form <- sum((backsolve(precomp_ou$L, residual, transpose = TRUE))^2)
    
    ll <- -0.5 * (length(y_vec) * log(2*pi) + precomp_ou$log_det + quad_form)
    
    if (!is.finite(ll)) return(-1e10)
    return(ll)
}

# Fast likelihood for single trait using precomputed noise structure
fast_single_trait_noise_ll <- function(y_vec, precomp_noise, theta_mu_pop_mean) {
    if (!precomp_noise$success) return(-1e10)
    
    mu_noise <- rep(c(theta_mu_pop_mean , theta_mu_pop_mean), precomp_noise$N)
    
    residual <- y_vec - mu_noise
    quad_form <- sum((backsolve(precomp_noise$L, residual, transpose = TRUE))^2)
    
    ll <- -0.5 * (length(y_vec) * log(2*pi) + precomp_noise$log_det + quad_form)
    
    if (!is.finite(ll)) return(-1e10)
    return(ll)
}

# Fast mixture likelihood summing over all traits
fast_mixture_likelihood_all_traits <- function(y_array, vcv_bm, alpha, sigma_x,
                                              sigma_noise1, sigma_noise2, p_noise,
                                              theta_mu_pop_mean, theta_mu_pop_sd, epsilon) {
    N <- nrow(y_array)
    n_traits <- ncol(y_array)
    
    # Precompute covariance structures once
    precomp_ou <- precompute_ou_cov_structure(vcv_bm, alpha, sigma_x,
                                              theta_mu_pop_sd, epsilon)
    precomp_noise <- precompute_noise_cov_structure(N, theta_mu_pop_sd,
                                                    sigma_noise1, sigma_noise2)
    
    if (!precomp_ou$success || !precomp_noise$success) return(-1e10)
    
    # Loop through traits and sum likelihoods
    total_ll <- 0
    for (k in 1:n_traits) {
        # Extract trait data as vector
        y_vec <- numeric(2*N)
        for (i in 1:N) {
            y_vec[(i-1)*2 + 1] <- y_array[i, k, 1]
            y_vec[(i-1)*2 + 2] <- y_array[i, k, 2]
        }
        
        # Compute component likelihoods
        ll_ou <- fast_single_trait_ou_ll(y_vec, precomp_ou, theta_mu_pop_mean)
        ll_noise <- fast_single_trait_noise_ll(y_vec, precomp_noise, theta_mu_pop_mean)
        
        # Mixture likelihood
        if (ll_ou > ll_noise) {
            ll_k <- ll_ou + log(1 - p_noise + p_noise * exp(ll_noise - ll_ou))
        } else {
            ll_k <- ll_noise + log(p_noise + (1 - p_noise) * exp(ll_ou - ll_noise))
        }
        
        if (!is.finite(ll_k)) return(-1e10)
        total_ll <- total_ll + ll_k
    }
    
    return(total_ll)
}

compute_phylo_mean_2d <- function(trait_values, tree) {
    tryCatch({
        vcov_matrix <- ape::vcv(tree)
        det_vcov <- det(vcov_matrix)
        if (abs(det_vcov) < 1e-10) {
            warning("VCV matrix near singular, using arithmetic mean")
            return(apply(trait_values, 1, mean, na.rm = TRUE))
        }
        vcov_inv <- solve(vcov_matrix)
        ones <- rep(1, ncol(trait_values))
        phylo_means <- numeric(2)
        for (d in 1:2) {
            trait_d <- trait_values[d, ]
            phylo_means[d] <- as.numeric((t(ones) %*% vcov_inv %*% trait_d) / (t(ones) %*% vcov_inv %*% ones))
        }
        return(phylo_means)
    }, error = function(e) {
        warning("Phylogenetic mean calculation failed, using arithmetic mean: ", e$message)
        return(apply(trait_values, 1, mean, na.rm = TRUE))
    })
}

# NEW: Transform parameters from constrained to unconstrained space
params_to_unconstrained <- function(alpha, sigma_x, sigma_noise1, sigma_noise2,
                                   p_noise, theta_mu_pop_mean, theta_mu_pop_sd) {
    c(log(alpha), log(sigma_x),
      log(sigma_noise1), log(sigma_noise2), qlogis(p_noise),
      theta_mu_pop_mean, log(theta_mu_pop_sd))
}

# NEW: Transform parameters from unconstrained to constrained space
unconstrained_to_params <- function(x) {
    list(
        alpha = exp(x[1]),
        sigma_x = exp(x[2]),
        sigma_noise1 = exp(x[3]),
        sigma_noise2 = exp(x[4]),
        p_noise = plogis(x[5]),
        theta_mu_pop_mean = x[6],
        theta_mu_pop_sd = exp(x[7])
    )
}

process_single_file_ml <- function(sim_file_path, optimizer, max_iter, reltol, n_starts) {
    epsilon<- 0.01
    
    library(ape); library(tools); library(optimx)
    sim_data <- readRDS(sim_file_path)
    trait_list <- sim_data$trait_list
    tree <- sim_data$tree
    n_traits <- length(trait_list)
    n_tips <- length(tree$tip.label)
    
    # Build y_array (no measurement error, so no se_array needed)
    y_array <- array(0, dim = c(n_tips, n_traits, 2))
    for (i in 1:n_tips) {
        for (k in 1:n_traits) {
            y_array[i, k, 1] <- trait_list[[k]][1, i]
            y_array[i, k, 2] <- trait_list[[k]][2, i]
        }
    }
    
    vcv_bm <- vcv(tree)
    tree_height <- max(node.depth.edgelength(tree))
    all_y1_vals <- as.vector(y_array[, , 1])
    all_y2_vals <- as.vector(y_array[, , 2])
    max_rate <- min(100.0, 40.0 / tree_height)
    data_scale_1 <- sqrt(max(0.01, var(all_y1_vals, na.rm = TRUE)))
    data_scale_2 <- sqrt(max(0.01, var(all_y2_vals, na.rm = TRUE)))
    max_sigma <- max(data_scale_1, data_scale_2) * 10.0

    if (!is.null(sim_data$trait_theta_mus) ) {
        theta_mu_init_mean <- mean(sim_data$trait_theta_mus, na.rm = TRUE)
        theta_mu_init_var <- max(var(sim_data$trait_theta_mus, na.rm = TRUE), 0.01)
    } else {
        trait_phylo_means <- array(0, dim=c(n_traits, 2))
        trait_theta_mus <- numeric(n_traits)
        for (k in 1:n_traits) {
            trait_phylo_means[k, ] <- compute_phylo_mean_2d(trait_list[[k]], tree)
            trait_theta_mus[k] <- trait_phylo_means[k, 1]
        }
        theta_mu_init_mean <- mean(trait_theta_mus, na.rm = TRUE)
        theta_mu_init_var <- max(var(trait_theta_mus, na.rm = TRUE), 0.01)
    }

    safe_theta_mu_mean <- if(is.finite(theta_mu_init_mean)) theta_mu_init_mean else 0.0

    # Objective function using fast likelihood with bounds enforcement
    objective <- function(params) {
        p <- unconstrained_to_params(params)
        
        # Enforce bounds (soft constraints via penalty)
        if (p$alpha > max_rate  || 
            p$sigma_x > max_sigma || 
            p$sigma_noise1 > max_sigma || p$sigma_noise2 > max_sigma ||
            p$theta_mu_pop_sd > max_sigma || 
            abs(p$theta_mu_pop_mean) > max_sigma*5 ) {
            return(1e10)
        }

        # Use fast likelihood that precomputes covariance once
        total_ll <- fast_mixture_likelihood_all_traits(
            y_array, vcv_bm, p$alpha,  p$sigma_x,
            p$sigma_noise1, p$sigma_noise2, p$p_noise,
            p$theta_mu_pop_mean, p$theta_mu_pop_sd, epsilon
        )
        
        if (!is.finite(total_ll)) return(1e10)
        return(-total_ll)
    }

    # Run optimization from multiple starting points
    best_result <- NULL
    best_ll <- -Inf
    
    set.seed(123)  # For reproducibility
    
    for (start_idx in 1:n_starts) {
        cat("\n=== Starting point", start_idx, "of", n_starts, "===\n")
        
        if (start_idx == 1) {
            # Use informed starting values for first attempt
            init_params <- params_to_unconstrained(
                alpha = 0.5,
                sigma_x = pmax(0.1, data_scale_1 * 0.5),
                sigma_noise1 = pmax(0.1, data_scale_1 * 0.3),
                sigma_noise2 = pmax(0.1, data_scale_2 * 0.3),
                p_noise = 0.2,
                theta_mu_pop_mean = safe_theta_mu_mean,
                theta_mu_pop_sd = pmax(0.1, sqrt(theta_mu_init_var))
            )
        } else {
            # Use random perturbations for subsequent attempts
            init_params <- params_to_unconstrained(
                alpha = runif(1, 0.1, 2.0),
                sigma_x = runif(1, 0.1, data_scale_1 * 2),
                sigma_noise1 = runif(1, 0.1, data_scale_1),
                sigma_noise2 = runif(1, 0.1, data_scale_2),
                p_noise = runif(1, 0.05, 0.4),
                theta_mu_pop_mean = rnorm(1, safe_theta_mu_mean, sqrt(theta_mu_init_var)),
                theta_mu_pop_sd = runif(1, 0.1, 2 * sqrt(theta_mu_init_var))
            )
        }

        cat("Initial objective value:", objective(init_params), "\n")

        opt_result <- optimx(par = init_params, fn = objective, method = optimizer,
                   control = list(maxit = max_iter, reltol = reltol, trace = 0))
        
        current_ll <- -opt_result$value[1]
        cat("Final log-likelihood:", current_ll, "\n")
        
        if (current_ll > best_ll ) { # && opt_result$convcode[1] == 0, we want the best result even if not converged
            best_ll <- current_ll
            best_result <- opt_result
            cat("*** New best result! ***\n")
        }
    }
    
    if (is.null(best_result)) {
        cat("Warning: No successful optimization runs\n")
        return(NULL)
    }
    
    cat("\n=== Best result across all starting points ===\n")
    cat("Best log-likelihood:", best_ll, "\n")
    print(best_result)

    best_params <- as.numeric(best_result[1, 1:11])
    p <- unconstrained_to_params(best_params)
    
    log_likelihood <- best_ll

    hessian_result <- tryCatch({
        hess <- numDeriv::hessian(objective, best_params)
        se_params <- sqrt(diag(solve(hess)))
        se_params
    }, error = function(e) {
        rep(NA, 11)
    })

    se_params <- hessian_result

    result <- list(
        estimates = list(
            alpha = p$alpha,
            sigma_x = p$sigma_x,
            sigma_noise1 = p$sigma_noise1,
            sigma_noise2 = p$sigma_noise2,
            p_noise = p$p_noise,
            theta_mu_pop_mean = p$theta_mu_pop_mean,
            theta_mu_pop_sd = p$theta_mu_pop_sd
        ),
        standard_errors = list(
            alpha_se = if(!is.na(se_params[1])) sqrt(p$alpha^2 * se_params[1]^2) else NA,
            sigma_x_se = if(!is.na(se_params[2])) sqrt(p$sigma_x^2 * se_params[2]^2) else NA,
            sigma_noise1_se = if(!is.na(se_params[3])) sqrt(p$sigma_noise1^2 * se_params[3]^2) else NA,
            sigma_noise2_se = if(!is.na(se_params[4])) sqrt(p$sigma_noise2^2 * se_params[4]^2) else NA,
            p_noise_se = if(!is.na(se_params[5])) abs(p$p_noise * (1-p$p_noise) * se_params[5]) else NA,
            theta_mu_pop_mean_se = if(!is.na(se_params[6])) se_params[6] else NA,
            theta_mu_pop_sd_se = if(!is.na(se_params[7])) sqrt(p$theta_mu_pop_sd^2 * se_params[7]^2) else NA
            ),
        optimization = list(
            log_likelihood = log_likelihood,
            aic = 2 * 11 - 2 * log_likelihood,
            bic = log(n_traits * n_tips * 2) * 11 - 2 * log_likelihood,
            convergence_code = best_result$convcode,
            n_iterations = best_result$fevals,
            optimizer = optimizer,
            n_starts = n_starts
        ),
        data_info = list(n_traits = n_traits, n_tips = n_tips, tree_height = tree_height, simulation_file = sim_file_path),
        true_parameters = sim_data$parameters,
        theta_init_info = list(theta_mu_mean = theta_mu_init_mean, theta_mu_var = theta_mu_init_var),
        timing = list(optimization = NA),
        settings = list(mixture_model = "2D-OU-Gaussian-ComponentTheta-ML-NoME", n_parameters = 7, 
                       parameterization = "2d_component_theta", measurement_error = FALSE)
    )

    # Compute trait classification using precomputed structures
    precomp_ou <- precompute_ou_cov_structure(vcv_bm, p$alpha, p$sigma_x,
                                              p$theta_mu_pop_sd, epsilon)
    precomp_noise <- precompute_noise_cov_structure(n_tips, p$theta_mu_pop_sd, 
                                                    p$sigma_noise1, p$sigma_noise2)
    
    trait_probs <- matrix(0, n_traits, 2)
    for (k in 1:n_traits) {
        y_vec <- numeric(2*n_tips)
        for (i in 1:n_tips) {
            y_vec[(i-1)*2 + 1] <- y_array[i, k, 1]
            y_vec[(i-1)*2 + 2] <- y_array[i, k, 2]
        }
        
        ll_ou <- suppressWarnings(fast_single_trait_ou_ll(y_vec, precomp_ou, 
                                                          p$theta_mu_pop_mean))
        ll_noise <- suppressWarnings(fast_single_trait_noise_ll(y_vec, precomp_noise,
                                                                p$theta_mu_pop_mean))
        
        log_p_ou_post <- log(1 - p$p_noise) + ll_ou
        log_p_noise_post <- log(p$p_noise) + ll_noise
        if (log_p_ou_post > log_p_noise_post) {
            log_normalizer <- log_p_ou_post + log(1 + exp(log_p_noise_post - log_p_ou_post))
        } else {
            log_normalizer <- log_p_noise_post + log(1 + exp(log_p_ou_post - log_p_noise_post))
        }
        trait_probs[k, 1] <- exp(log_p_ou_post - log_normalizer)
        trait_probs[k, 2] <- exp(log_p_noise_post - log_normalizer)
    }

    result$trait_classification <- list(
        trait_ou_prob = trait_probs[, 1],
        trait_noise_prob = trait_probs[, 2],
        mean_noise_prob = mean(trait_probs[, 2]),
        n_likely_noise = sum(trait_probs[, 2] > 0.5),
        max_noise_prob = max(trait_probs[, 2])
    )

    output_file <- file.path(dirname(sim_file_path), paste0("ml_fit_", tools::file_path_sans_ext(basename(sim_file_path)), ".rds"))
    suppressWarnings(saveRDS(result, output_file))
    return(output_file)
}

# Main loop (no parallel)
cat("Processing files... \n")
overall_start <- Sys.time()
all_results <- lapply(sim_files, function(file) process_single_file_ml(file, args$optimizer, args$max_iter, args$reltol, args$n_starts))
total_time <- difftime(Sys.time(), overall_start, units="mins")
successful_results <- all_results[!sapply(all_results, is.null)]
cat("Done.\n")
cat(paste("Processed", length(successful_results), "/", length(sim_files), "files in", round(total_time, 2), "minutes\n"))

# Save detailed summary
detailed_summary_file <- file.path(args$input_dir, "ml_detailed_summary.rds")
detailed_summary <- list(
    results_files = successful_results,
    processing_info = list(
        total_files = length(sim_files),
        successful_files = length(successful_results),
        total_time_minutes = as.numeric(total_time),
        optimizer_used = args$optimizer,
        n_starts = args$n_starts,
        max_iter = args$max_iter,
        parallel_used = FALSE,
        timestamp = Sys.time()
    )
)
suppressWarnings(saveRDS(detailed_summary, detailed_summary_file))
cat("Detailed summary saved:", detailed_summary_file, "\n")

success_rate <- round(100 * length(successful_results) / max(1, length(sim_files)), 1)
cat("Completed with", success_rate, "% success rate.\n")

log_msg("=== SCRIPT COMPLETED ===")
log_msg(paste("Success rate:", success_rate, "%"))