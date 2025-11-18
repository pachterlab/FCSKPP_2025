#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(rstan)
  library(MASS)
  library(Matrix)
  library(expm)
})

parse_cmd_args <- function() {
  option_list <- list(
    make_option(c("-o", "--output_dir"), type="character", default="sim2_UB_40",
                help="Output directory for simulated data files"),
    make_option(c("-t", "--tree"), type="character", default=NULL,
                help="Newick file with phylogenetic tree (if not provided, a default tree will be used)"),
    make_option(c("--size"), type="integer", default=NULL,
                help="Size of tree to randomly simulate if desired (if not provided, a default tree will be used)"),
    make_option(c("-n", "--ntraits"), type="integer", default=150,
                help="Number of traits to simulate per repeat"),
    make_option(c("-r", "--repeats"), type="integer", default=100,
                help="Number of parameter combinations to simulate"),
    make_option(c("--alpha_min"), type="numeric", default=0.1,
                help="Minimum value for alpha parameter"),
    make_option(c("--alpha_max"), type="numeric", default=40.0,
                help="Maximum value for alpha parameter"),
    make_option(c("--beta_min"), type="numeric", default=0.1,
                help="Minimum value for beta parameter"),
    make_option(c("--beta_max"), type="numeric", default=40.0,
                help="Maximum value for beta parameter"),
    make_option(c("--sigma1_min"), type="numeric", default=0.1,
                help="Minimum value for sigma_x1 parameter"),
    make_option(c("--sigma1_max"), type="numeric", default=10.0,
                help="Maximum value for sigma_x1 parameter"),
    make_option(c("--sigma2_min"), type="numeric", default=0.1,
                help="Minimum value for sigma_x2 parameter"),
    make_option(c("--sigma2_max"), type="numeric", default=10.0,
                help="Maximum value for sigma_x2 parameter"),
    make_option(c("--sigma_noise1_min"), type="numeric", default=0.1,
                help="Minimum value for sigma_noise1 parameter"),
    make_option(c("--sigma_noise1_max"), type="numeric", default=10.0,
                help="Maximum value for sigma_noise1 parameter"),
    make_option(c("--sigma_noise2_min"), type="numeric", default=0.1,
                help="Minimum value for sigma_noise2 parameter"),
    make_option(c("--sigma_noise2_max"), type="numeric", default=10.0,
                help="Maximum value for sigma_noise2 parameter"),
    make_option(c("--p_min"), type="numeric", default=0.01,
                help="Minimum probability of noise traits"),
    make_option(c("--p_max"), type="numeric", default=0.6,
                help="Maximum probability of noise traits"),
    make_option(c("--theta_mu_mean"), type="numeric", default=0,
                help="Mean of theta_mu distribution for population sampling"),
    make_option(c("--theta_mu_sd"), type="numeric", default=1.0,
                help="Standard deviation of theta_mu distribution for population sampling"),
    make_option(c("--theta_b_mean"), type="numeric", default=0,
                help="Mean of theta_b distribution for population sampling"),
    make_option(c("--theta_b_sd"), type="numeric", default=1.0,
                help="Standard deviation of theta_b distribution for population sampling"),
    make_option(c("--me_scale"), type="numeric", default=0.0,
                help="Scale factor for measurement error"),
    make_option(c("--seed"), type="integer", default=NULL,
                help="Random seed for reproducibility"),
    make_option(c("--verbose"), type="logical", default=FALSE,
                help="Verbose output"),
    make_option(c("--quiet"), type="logical", default=FALSE,
                help="Suppress most output")
  )
  
  opt_parser <- OptionParser(option_list=option_list)
  return(parse_args(opt_parser))
}

get_tree <- function(tree_file = NULL, tree_size=NULL, quiet=FALSE) {
  if (!is.null(tree_file) && file.exists(tree_file)) {
    tree <- read.tree(tree_file)
    if (!quiet) cat("Using tree from:", tree_file, "\n")
  } 
  else if (!is.null(tree_size)) { 
    tree <- rcoal(tree_size)
    tree_height <- max(node.depth.edgelength(tree))
    tree$edge.length <- tree$edge.length / tree_height
    if (!quiet) cat("Simulated tree with", tree_size, "tips\n")
  }
  else {
    newick_tree <- "(Xenopus_tropicalis:351.68654000,(Sus_scrofa:94.00000000,((Rattus_norvegicus:11.64917000,Mus_musculus:11.64917000)'14':75.55083000,(Homo_sapiens:28.82000000,Macaca_mulatta:28.82000000)'13':58.38000000)'25':6.80000000)'37':257.68654000);"
    tree <- read.tree(text = newick_tree)
    tree$node.label <- NULL
    
    tree_height <- max(vcv(tree)[1, ])
    tree$edge.length <- tree$edge.length / tree_height
    if (!quiet) cat("Using default 5-species tree\n")
  }
  
  return(tree)
}

simulate_2d_ou <- function(tree, alpha, beta, sigma_x1, sigma_x2, theta, x0) {
  # Load PCMBase
  if (!requireNamespace("PCMBase", quietly = TRUE)) {
    stop("PCMBase must be installed.")
  }
  library(PCMBase)
  
  n_tips <- length(tree$tip.label)
  
  # Define OU model
  modelOU <- PCM("OU__Global_X0__Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
                 k = 2)
  
  # Fill in parameters
  modelOU$X0[1] <- x0[1]
  modelOU$X0[2] <- x0[2]

  modelOU$Theta[1] <- theta[1]
  modelOU$Theta[2] <- theta[2]
  
  # H_matrix changed from model 1.
  modelOU$H[1,1] <- alpha
  modelOU$H[1,2] <- 0
  modelOU$H[2,1] <- -beta
  modelOU$H[2,2] <- beta   
  modelOU$Sigma_x[1] <- sigma_x1
  modelOU$Sigma_x[4] <-  sigma_x2 # PCMBase uses σ, not σ^2

  # Simulate
  sim_result <- PCMSim(tree, model = modelOU, X0 = x0)
  
  # Extract tip values only
  tip_vals <- sim_result[,1:n_tips , drop = FALSE]
  
  # Transpose to 2 x n_tips matrix (same format as your old code)
  trait_matrix <- tip_vals
  
  return(trait_matrix)
}


generate_me_values_2d <- function(n_traits, n_tips, scale=0.1) {
  me_values <- list()
  
  for (i in 1:n_traits) {
    me_matrix <- matrix(abs(rnorm(2*n_tips, 0, scale)), nrow=2, ncol=n_tips)
    me_values[[i]] <- me_matrix
  }
  
  return(me_values)
}

simulate_2d_traits_mixture <- function(tree, n_traits, me_values, params) {
  trait_list <- list()
  trait_types <- character(n_traits)
  trait_theta_mus <- numeric(n_traits)
  trait_theta_bs <- numeric(n_traits)
  trait_thetas <- array(0, dim=c(n_traits, 2))
  
  n_tips <- length(tree$tip.label)
  
  is_noise <- rbinom(n_traits, 1, params$p_noise) == 1
  n_phylo <- sum(!is_noise)
  n_noise <- sum(is_noise)
  
  for (i in 1:n_traits) {
    trait_theta_mu <- rnorm(1, mean = params$theta_mu_pop_mean, sd = params$theta_mu_pop_sd)
    trait_theta_b <- rnorm(1, mean = params$theta_b_pop_mean, sd = params$theta_b_pop_sd)
    
    trait_theta_mus[i] <- trait_theta_mu
    trait_theta_bs[i] <- trait_theta_b
    
    trait_theta <- c(trait_theta_b, trait_theta_b - trait_theta_mu)
    trait_thetas[i, ] <- trait_theta
    
    trait_x0 <- trait_theta
    
    if (is_noise[i]) {
      trait_data <- matrix(0, nrow=2, ncol=n_tips)
      trait_data[1, ] <- rnorm(n_tips, mean = trait_theta[1], sd = params$sigma_noise1)
      trait_data[2, ] <- rnorm(n_tips, mean = trait_theta[2], sd = params$sigma_noise2)
      trait_types[i] <- "noise"
    } else {
      trait_data <- simulate_2d_ou(tree, params$alpha, params$beta, 
                                       params$sigma_x1, params$sigma_x2, 
                                       trait_theta, trait_x0)
      trait_types[i] <- "phylo"
    }
    
    trait_data_with_error <- trait_data + me_values[[i]] * matrix(rnorm(2*n_tips), nrow=2)
    
    trait_list[[i]] <- trait_data_with_error
  }
  
  return(list(
    traits = trait_list, 
    types = trait_types, 
    theta_mus = trait_theta_mus,
    theta_bs = trait_theta_bs,
    thetas = trait_thetas,
    n_phylo = n_phylo, 
    n_noise = n_noise
  ))
}

generate_random_params <- function(args) {
  params <- list()
  
  params$alpha <- runif(1, args$alpha_min, args$alpha_max)
  params$beta <- runif(1, args$beta_min, args$beta_max)
  params$sigma_x1 <- runif(1, args$sigma1_min, args$sigma1_max)
  params$sigma_x2 <- runif(1, args$sigma2_min, args$sigma2_max)
  params$sigma_noise1 <- runif(1, args$sigma_noise1_min, args$sigma_noise1_max)
  params$sigma_noise2 <- runif(1, args$sigma_noise2_min, args$sigma_noise2_max)
  params$p_noise <- runif(1, args$p_min, args$p_max)
  
  params$theta_mu_pop_mean <- args$theta_mu_mean
  params$theta_mu_pop_sd <- args$theta_mu_sd
  params$theta_b_pop_mean <- args$theta_b_mean
  params$theta_b_pop_sd <- args$theta_b_sd

  return(params)
}

# ## Use log uniform sampling distributions for all positive values.
# generate_random_params <- function(args) {
#   log_uniform <- function(min_val, max_val) {
#     exp(runif(1, log(min_val), log(max_val)))
#   }
  
#   params <- list()
  
#   params$alpha <- log_uniform(args$alpha_min, args$alpha_max)
#   params$beta <- log_uniform(args$beta_min, args$beta_max)
#   params$sigma_x1 <- log_uniform(args$sigma1_min, args$sigma1_max)
#   params$sigma_x2 <- log_uniform(args$sigma2_min, args$sigma2_max)
#   params$sigma_noise1 <- log_uniform(args$sigma_noise1_min, args$sigma_noise1_max)
#   params$sigma_noise2 <- log_uniform(args$sigma_noise2_min, args$sigma_noise2_max)
  
#   params$p_noise <- runif(1, args$p_min, args$p_max)
  
#   params$theta_mu_pop_mean <- args$theta_mu_mean
#   params$theta_mu_pop_sd <- args$theta_mu_sd
#   params$theta_b_pop_mean <- args$theta_b_mean
#   params$theta_b_pop_sd <- args$theta_b_sd
  
#   return(params)
# }


create_filename <- function(output_dir, params, repeat_num) {
  alpha_str <- paste0("A", sprintf("%.2f", params$alpha))
  beta_str <- paste0("B", sprintf("%.2f", params$beta))
  sigma1_str <- paste0("S1_", sprintf("%.2f", params$sigma_x1))
  sigma2_str <- paste0("S2_", sprintf("%.2f", params$sigma_x2))
  sn1_str <- paste0("SN1_", sprintf("%.2f", params$sigma_noise1))
  sn2_str <- paste0("SN2_", sprintf("%.2f", params$sigma_noise2))
  p_str <- paste0("P", sprintf("%.3f", params$p_noise))
  
  theta_str <- paste0("ThetaMu_", sprintf("%.1f", params$theta_mu_pop_sd), 
                     "_ThetaB_", sprintf("%.1f", params$theta_b_pop_sd))
  
  filename <- file.path(output_dir, paste0("sim_model2_", repeat_num, "_", alpha_str, "_", beta_str, "_", 
                                          sigma1_str, "_", sigma2_str, "_", sn1_str, "_", sn2_str, "_", 
                                          p_str, "_", theta_str, ".rds"))
  
  return(filename)
}

main <- function() {
  args <- parse_cmd_args()
  
  if (!is.null(args$seed)) {
    set.seed(args$seed)
  }
  
  if (!dir.exists(args$output_dir)) {
    dir.create(args$output_dir, recursive = TRUE)
    if (!args$quiet) cat("Created output directory:", args$output_dir, "\n")
  }
  
  if (!args$quiet) cat("Loading phylogenetic tree... ")
  tree <- get_tree(args$tree, args$size, args$quiet)
  
  summary_df <- data.frame()
  
  if (!args$quiet) {
    cat("\nStarting 2D OU mixture model simulation:\n")
    cat("- Repeats:", args$repeats, "\n")
    cat("- Traits per repeat:", args$ntraits, "\n")
    cat("- Alpha range: [", args$alpha_min, ",", args$alpha_max, "]\n")
    cat("- Beta range: [", args$beta_min, ",", args$beta_max, "]\n")
    cat("- Sigma_x1 range: [", args$sigma1_min, ",", args$sigma1_max, "]\n")
    cat("- Sigma_x2 range: [", args$sigma2_min, ",", args$sigma2_max, "]\n")
    cat("- Sigma_noise1 range: [", args$sigma_noise1_min, ",", args$sigma_noise1_max, "]\n")
    cat("- Sigma_noise2 range: [", args$sigma_noise2_min, ",", args$sigma_noise2_max, "]\n")
    cat("- Noise probability range: [", args$p_min, ",", args$p_max, "]\n")
    cat("- Theta_mu ~ N(", args$theta_mu_mean, ",", args$theta_mu_sd, ") per trait\n")
    cat("- Theta_b ~ N(", args$theta_b_mean, ",", args$theta_b_sd, ") per trait\n\n")
  }
  
  progress_interval <- max(1, args$repeats %/% 10)
  
  for (repeat_num in 1:args$repeats) {
    if (!args$quiet && (repeat_num %% progress_interval == 0 || repeat_num == 1)) {
      cat("Progress:", repeat_num, "/", args$repeats, "\n")
    }
    
    if (args$verbose) cat(paste0("\n=== Repeat ", repeat_num, " of ", args$repeats, " ===\n"))
    
    params <- generate_random_params(args)
    
    if (args$verbose) {
      cat("Selection matrix parameters: alpha =", sprintf("%.3f", params$alpha), 
          "beta =", sprintf("%.3f", params$beta), "\n")
      cat("Diffusion parameters: sigma_x1 =", sprintf("%.3f", params$sigma_x1),
          "sigma_x2 =", sprintf("%.3f", params$sigma_x2), "\n")
      cat("Noise parameters: sigma_noise1 =", sprintf("%.3f", params$sigma_noise1),
          "sigma_noise2 =", sprintf("%.3f", params$sigma_noise2),
          "P_noise =", sprintf("%.3f", params$p_noise), "\n")
    }
    
    if (args$verbose) cat("Generating 2D measurement error matrices...\n")
    me_values <- generate_me_values_2d(args$ntraits, length(tree$tip.label), scale=args$me_scale)
    
    if (args$verbose) cat(paste0("Simulating ", args$ntraits, " 2D traits under mixture model...\n"))
    sim_result <- simulate_2d_traits_mixture(tree, args$ntraits, me_values, params)
    
    if (args$verbose) {
      cat("Simulated", sim_result$n_phylo, "phylogenetic traits and", sim_result$n_noise, "noise traits\n")
    }
    
    result <- list(
      trait_list = sim_result$traits,
      trait_types = sim_result$types,
      trait_theta_mus = sim_result$theta_mus,
      trait_theta_bs = sim_result$theta_bs,
      trait_thetas = sim_result$thetas,
      me_values = me_values,
      tree = tree,
      parameters = list(
        alpha = params$alpha,
        beta = params$beta,
        sigma_x1 = params$sigma_x1,
        sigma_x2 = params$sigma_x2,
        sigma_noise1 = params$sigma_noise1,
        sigma_noise2 = params$sigma_noise2,
        p_noise = params$p_noise,
        theta_mu_pop_mean = params$theta_mu_pop_mean,
        theta_mu_pop_sd = params$theta_mu_pop_sd,
        theta_b_pop_mean = params$theta_b_pop_mean,
        theta_b_pop_sd = params$theta_b_pop_sd
      ),
      n_traits = args$ntraits,
      n_phylo_traits = sim_result$n_phylo,
      n_noise_traits = sim_result$n_noise,
      tip_labels = tree$tip.label,
      model_type = "2D_OU_mixture"
    )
    
    output_file <- create_filename(args$output_dir, params, repeat_num)
    
    if (args$verbose) cat(paste0("Saving results to ", output_file, "...\n"))
    saveRDS(result, file = output_file)
    
    summary_row <- data.frame(
      repeat_num = repeat_num,
      alpha = params$alpha,
      beta = params$beta,
      sigma_x1 = params$sigma_x1,
      sigma_x2 = params$sigma_x2,
      sigma_noise1 = params$sigma_noise1,
      sigma_noise2 = params$sigma_noise2,
      p_noise = params$p_noise,
      theta_mu_pop_mean = params$theta_mu_pop_mean,
      theta_mu_pop_sd = params$theta_mu_pop_sd,
      theta_b_pop_mean = params$theta_b_pop_mean,
      theta_b_pop_sd = params$theta_b_pop_sd,
      n_phylo_traits = sim_result$n_phylo,
      n_noise_traits = sim_result$n_noise,
      filename = basename(output_file)
    )
    summary_df <- rbind(summary_df, summary_row)
  }
  
  summary_file <- file.path(args$output_dir, "simulation_summary_2d.csv")
  write.csv(summary_df, file = summary_file, row.names = FALSE)
  if (!args$quiet) cat("\nSummary saved to:", summary_file, "\n")
  
  if (!args$quiet) {
    cat("\n=== 2D OU Simulation Complete ===\n")
    cat("Total simulations:", args$repeats, "\n")
    cat("Average phylogenetic traits per simulation:", round(mean(summary_df$n_phylo_traits), 1), "\n")
    cat("Average noise traits per simulation:", round(mean(summary_df$n_noise_traits), 1), "\n")
    cat("Results saved to:", args$output_dir, "\n")
  }
}

if (!interactive()) {
  main()
}
