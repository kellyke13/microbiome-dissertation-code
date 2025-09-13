# ==================== Simulator for ZIPG ==============================

simulate_zipg_data <- function(n = 100, K = 30, d = 3, d_star = 2,
                               beta_sd = 0.5, seed = 123) {
  set.seed(seed)
  X <- matrix(rnorm(n * d), nrow = n, ncol = d)
  X_star <- matrix(rnorm(n * d_star), nrow = n, ncol = d_star)
  M <- rlnorm(n, meanlog = log(15000), sdlog = 0.3)
  
  beta0 <- rnorm(K, mean = 0, sd = 1)
  beta <- matrix(rnorm(K * d, mean = 0, sd = beta_sd), nrow = K)
  beta0_star <- rnorm(K, mean = 0, sd = 1)
  beta_star <- matrix(rnorm(K * d_star, mean = 0, sd = beta_sd), nrow = K)
  p_k <- runif(K, min = 0.2, max = 0.5)
  
  W <- matrix(0, nrow = n, ncol = K)
  Lambda <- matrix(0, n, K)
  Theta <- matrix(0, n, K)
  U <- matrix(0, n, K)
  Mu <- matrix(0, n, K)
  
  for (i in 1:n) {
    for (k in 1:K) {
      log_lambda_ik <- beta0[k] + sum(X[i, ] * beta[k, ]) + log(M[i])
      lambda_ik <- exp(log_lambda_ik)
      
      log_theta_ik <- beta0_star[k] + sum(X_star[i, ] * beta_star[k, ])
      theta_ik <- exp(log_theta_ik)
      
      U_ik <- rgamma(1, shape = 1 / theta_ik, rate = 1 / theta_ik)
      mu_ik <- lambda_ik * U_ik
      
      W[i, k] <- if (runif(1) < p_k[k]) 0 else rpois(1, mu_ik)
      
      Lambda[i, k] <- lambda_ik
      Theta[i, k] <- theta_ik
      U[i, k] <- U_ik
      Mu[i, k] <- mu_ik
    }
  }
  
  list(W = W, X = X, X_star = X_star, M = M,
       beta0 = beta0, beta = beta, beta0_star = beta0_star, beta_star = beta_star,
       p_k = p_k, Lambda = Lambda, Theta = Theta, U = U, Mu = Mu)
}

# ========================== Extraction + metrics ===============================

extract_coef_from_summary <- function(summary_df, d = 3, d_star = 2) {
  ests <- summary_df$Estimation
  names(ests) <- rownames(summary_df)
  
  beta <- ests[grep("^beta[0-9]+$", names(ests))]
  beta_star <- ests[grep("^beta[0-9]+\\*$", names(ests))]
  
  beta <- tail(beta, d)
  beta_star <- tail(beta_star, d_star)
  
  list(
    beta0 = ests["beta0"],
    beta = beta,
    beta0_star = ests["beta0*"],
    beta_star = beta_star,
    gamma = ests["gamma"]
  )
}

get_recovery_metrics_detailed <- function(beta0_hat, beta0_true,
                                          beta_hat_mat, beta_true_mat,
                                          beta0_star_hat, beta0_star_true,
                                          beta_star_hat_mat, beta_star_true,
                                          gamma_hat, gamma_true) {
  beta_full_hat  <- cbind(beta0_hat, beta_hat_mat)
  beta_full_true <- cbind(beta0_true, beta_true_mat)
  colnames(beta_full_hat) <- colnames(beta_full_true) <- paste0("beta", 0:(ncol(beta_full_hat) - 1))
  
  beta_metrics <- data.frame(
    Parameter = colnames(beta_full_hat),
    Bias = colMeans(beta_full_hat - beta_full_true),
    RMSE = sqrt(colMeans((beta_full_hat - beta_full_true)^2))
  )
  
  beta_star_full_hat  <- cbind(beta0_star_hat, beta_star_hat_mat)
  beta_star_full_true <- cbind(beta0_star_true, beta_star_true)
  colnames(beta_star_full_hat) <- colnames(beta_star_full_true) <- paste0("beta", 0:(ncol(beta_star_full_hat) - 1), "*")
  
  beta_star_metrics <- data.frame(
    Parameter = colnames(beta_star_full_hat),
    Bias = colMeans(beta_star_full_hat - beta_star_full_true),
    RMSE = sqrt(colMeans((beta_star_full_hat - beta_star_full_true)^2))
  )
  
  gamma_metrics <- data.frame(
    Parameter = "gamma",
    Bias = mean(gamma_hat - gamma_true),
    RMSE = sqrt(mean((gamma_hat - gamma_true)^2))
  )
  
  summary_table <- rbind(beta_metrics, beta_star_metrics, gamma_metrics)
  rownames(summary_table) <- NULL
  summary_table
}


# ======================== Single run for a given n ============================

run_one_recovery <- function(n, K = 30, d = 3, d_star = 2, beta_sd = 0.5,
                             seed = 1, verbose = FALSE) {
  sim <- simulate_zipg_data(n = n, K = K, d = d, d_star = d_star,
                            beta_sd = beta_sd, seed = seed)
  
  X_df  <- as.data.frame(sim$X);      colnames(X_df)  <- paste0("X",  seq_len(ncol(X_df)))
  Xs_df <- as.data.frame(sim$X_star); colnames(Xs_df) <- paste0("Xs", seq_len(ncol(Xs_df)))
  dat   <- cbind(X_df, Xs_df)
  
  form_X  <- as.formula(paste("~",  paste(colnames(X_df),  collapse = " + ")))
  form_Xs <- as.formula(paste("~",  paste(colnames(Xs_df), collapse = " + ")))
  
  suppressPackageStartupMessages(library(ZIPG))
  fits <- vector("list", K)
  for (k in seq_len(K)) {
    if (verbose && k %% 5 == 0) cat(sprintf("  -> Fitting OTU %d/%d\n", k, K))
    wk <- sim$W[, k]
    fits[[k]] <- tryCatch({
      ZIPG_main(data = dat, W = wk, M = sim$M, X = form_X, X_star = form_Xs)
    }, error = function(e) NULL)
  }
  
  summary_list <- lapply(fits, function(f) tryCatch(ZIPG_summary(f), error = function(e) NULL))
  valid_idx <- which(!sapply(summary_list, is.null))
  if (length(valid_idx) == 0) {
    warning(sprintf("No OTU converged for n = %d (seed = %d). Returning NA metrics.", n, seed))
    out_tbl <- data.frame(Parameter = c(paste0("beta", 0:d), paste0("beta", 0:d, "*"), "gamma"),
                          Bias = NA_real_, RMSE = NA_real_, n = n, seed = seed)
    return(list(metrics_tbl = out_tbl, pairs = NULL))
  }
  
  coefs_list <- lapply(summary_list[valid_idx], extract_coef_from_summary, d = d, d_star = d_star)
  
  beta0_hat         <- sapply(coefs_list, function(l) l$beta0)
  beta_hat_mat      <- t(sapply(coefs_list, function(l) l$beta))
  beta0_star_hat    <- sapply(coefs_list, function(l) l$beta0_star)
  beta_star_hat_mat <- t(sapply(coefs_list, function(l) l$beta_star))
  gamma_hat         <- sapply(coefs_list, function(l) l$gamma)
  
  beta0_true       <- sim$beta0[valid_idx]
  beta_true_mat    <- sim$beta[valid_idx, , drop = FALSE]
  beta0_star_true  <- sim$beta0_star[valid_idx]
  beta_star_true   <- sim$beta_star[valid_idx, , drop = FALSE]
  gamma_true       <- log(sim$p_k[valid_idx] / (1 - sim$p_k[valid_idx]))
  
  tbl <- get_recovery_metrics_detailed(
    beta0_hat, beta0_true,
    beta_hat_mat, beta_true_mat,
    beta0_star_hat, beta0_star_true,
    beta_star_hat_mat, beta_star_true,
    gamma_hat, gamma_true
  )
  tbl$n <- n
  tbl$seed <- seed
  
  pairs <- list(
    beta0_true = beta0_true, beta0_hat = beta0_hat,
    beta_true_mat = beta_true_mat, beta_hat_mat = beta_hat_mat,
    beta0_star_true = beta0_star_true, beta0_star_hat = beta0_star_hat,
    beta_star_true = beta_star_true, beta_star_hat_mat = beta_star_hat_mat,
    gamma_true = gamma_true, gamma_hat = gamma_hat
  )
  
  list(metrics_tbl = tbl, pairs = pairs)
}


append_pairs <- function(A, B) {
  if (is.null(A)) return(B)
  if (is.null(B)) return(A)
  A$beta0_true        <- c(A$beta0_true, B$beta0_true)
  A$beta0_hat         <- c(A$beta0_hat,  B$beta0_hat)
  A$beta_true_mat     <- rbind(A$beta_true_mat,     B$beta_true_mat)
  A$beta_hat_mat      <- rbind(A$beta_hat_mat,      B$beta_hat_mat)
  A$beta0_star_true   <- c(A$beta0_star_true, B$beta0_star_true)
  A$beta0_star_hat    <- c(A$beta0_star_hat,  B$beta0_star_hat)
  A$beta_star_true    <- rbind(A$beta_star_true,    B$beta_star_true)
  A$beta_star_hat_mat <- rbind(A$beta_star_hat_mat, B$beta_star_hat_mat)
  A$gamma_true        <- c(A$gamma_true, B$gamma_true)
  A$gamma_hat         <- c(A$gamma_hat,  B$gamma_hat)
  A
}


# ================= Run over n_grid x reps ==================================

run_recovery_over_n <- function(n_grid = c(50, 100, 200, 400),
                                reps = 20,
                                K = 30, d = 3, d_star = 2, beta_sd = 0.5,
                                base_seed = 20240901, verbose_each = FALSE) {
  
  all_runs <- list()
  pairs_by_n <- vector("list", length(n_grid)); names(pairs_by_n) <- as.character(n_grid)
  idx <- 1
  
  for (n in n_grid) {
    cur_pairs <- NULL
    for (r in seq_len(reps)) {
      seed <- base_seed + 1000 * which(n_grid == n) + r
      one <- run_one_recovery(n = n, K = K, d = d, d_star = d_star,
                              beta_sd = beta_sd, seed = seed, verbose = verbose_each)
      all_runs[[idx]] <- one$metrics_tbl
      cur_pairs <- append_pairs(cur_pairs, one$pairs)
      idx <- idx + 1
    }
    pairs_by_n[[as.character(n)]] <- cur_pairs
  }
  
  all_df <- do.call(rbind, all_runs)
  
  # aggregate by n x Parameter: mean Bias/RMSE over reps
  library(dplyr)
  summary_by_n <- all_df %>%
    group_by(n, Parameter) %>%
    summarise(
      Bias = mean(Bias, na.rm = TRUE),
      RMSE = mean(RMSE, na.rm = TRUE),
      .groups = "drop"
    )
  
  # order params for nice printing
  reorder_params <- function(df, d, d_star) {
    beta_names  <- paste0("beta", 0:d)
    betaS_names <- paste0("beta", 0:d, "*")
    all_levels  <- c(beta_names, betaS_names, "gamma")
    df$Parameter <- factor(df$Parameter, levels = all_levels)
    df[order(df$Parameter), ]
  }
  
  per_n_tables <- lapply(split(summary_by_n, summary_by_n$n),
                         function(df) reorder_params(df[, c("Parameter","Bias","RMSE")], d, d_star))
  
  list(per_n = per_n_tables, raw = all_df, summary = summary_by_n, pairs_by_n = pairs_by_n,
       dims = list(d = d, d_star = d_star))
}

# run the result

res <- run_recovery_over_n(
  n_grid = c(50, 100, 200, 400),
  reps = 20,
  K = 30,
  d = 3, d_star = 2,
  beta_sd = 0.5,
  base_seed = 20240901,
  verbose_each = FALSE
)

for (n_val in names(res$per_n)) {
  df <- res$per_n[[n_val]]
  cat("\n\n=====================\n")
  cat(sprintf(" Parameter recovery (ZIPG), n = %s\n", n_val))
  cat("=====================\n")
  print(df, row.names = FALSE, digits = 4)
}


# =================================== plotting ======================================
plot_recovery <- function(true, est, name) {
  plot(true, est,
       main = paste("Recovery of", name),
       xlab = paste("True", name),
       ylab = paste("Estimated", name),
       pch = 19, col = "#1f77b4", cex = 1.1)
  abline(0, 1, col = "red", lty = 2)
  grid()
}


plot_all_recoveries_by_n <- function(pairs_by_n, d, d_star) {
  for (n_val in names(pairs_by_n)) {
    pairs <- pairs_by_n[[n_val]]
    if (is.null(pairs)) {
      message(sprintf("No pairs available for n = %s (likely no successful fits).", n_val))
      next
    }
    cat(sprintf("\n--- Plotting recovery for n = %s ---\n", n_val))
    
    ## β0
    plot_recovery(pairs$beta0_true, pairs$beta0_hat, paste0("β0 (n=", n_val, ")"))
    
    ## β slopes
    for (j in seq_len(d)) {
      plot_recovery(pairs$beta_true_mat[, j], pairs$beta_hat_mat[, j],
                    paste0("β", j, " (n=", n_val, ")"))
    }
    
    ## β0*
    plot_recovery(pairs$beta0_star_true, pairs$beta0_star_hat, paste0("β0* (n=", n_val, ")"))
    
    ## β* slopes
    for (j in seq_len(d_star)) {
      plot_recovery(pairs$beta_star_true[, j], pairs$beta_star_hat_mat[, j],
                    paste0("β", j, "* (n=", n_val, ")"))
    }
    
    ## γ
    plot_recovery(pairs$gamma_true, pairs$gamma_hat, paste0("γ (n=", n_val, ")"))
  }
}

# plots per n
plot_all_recoveries_by_n(res$pairs_by_n, d = res$dims$d, d_star = res$dims$d_star)
