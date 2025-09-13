library(ZIPG)
library(dplyr)
library(ggplot2)
library(patchwork)

# ======================== Scenario grid ==========================================
p_levels      <- c(0.1, 0.3, 0.5)       # zero-inflation
disp_shifts   <- c(-1, 0, 1)            # log-dispersion intercept shift
beta_scales   <- c(0.5, 1.0, 1.5)       # effect size multipliers

scenario_grid <- expand.grid(
  p = p_levels,
  disp_shift = disp_shifts,
  beta_scale = beta_scales
)

scenario_grid$scenario_id <- seq_len(nrow(scenario_grid))
print(scenario_grid)

# ==================== generate count matrix ==================================

simulate_zipg_ext <- function(n = 200, K = 30, 
                              p = 0.3, disp_shift = 0, beta_scale = 1,
                              m_mean = 1.5e4, m_sd = 2e3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Sequencing depth
  M <- pmax(round(rnorm(n, m_mean, m_sd)), 1000)
  
  # Covariates: binary + 2 continuous
  X <- cbind(
    x1 = rbinom(n, 1, 0.5),
    x2 = rnorm(n),
    x3 = rnorm(n)
  )
  p_cov <- ncol(X)
  
  # Storage for counts
  W <- matrix(0L, n, K)
  
  # Baseline dispersion
  theta_base <- 0.5
  theta <- exp(log(theta_base) + disp_shift)  # apply shift
  
  for (k in 1:K) {
    # Generate abundance coefficients
    beta0 <- rnorm(1, mean = -2, sd = 0.7)
    beta  <- rnorm(p_cov, mean = 0, sd = 0.3) * beta_scale
    
    # Compute log lambda and lambda
    log_lambda <- beta0 + X %*% beta + log(M)
    lambda <- as.vector(exp(log_lambda))
    
    # Gamma noise
    U <- rgamma(n, shape = 1/theta, scale = theta)  # mean 1, var = theta
    
    # Poisson counts
    Y <- rpois(n, lambda * U)
    
    # Apply zero-inflation
    ZI <- rbinom(n, 1, p)  # 1 = extra zero
    W[, k] <- ifelse(ZI == 1, 0, Y)
  }
  
  return(W)
}

evaluate_counts <- function(W) {
  n <- nrow(W); K <- ncol(W)
  
  # Sparsity
  overall_zero_prop <- mean(W == 0)
  
  # Mean–variance relationship
  taxon_means <- colMeans(W)
  taxon_vars  <- apply(W, 2, var)
  dispersion_index <- taxon_vars / taxon_means  # v/m
  
  # log-log regression slope
  valid <- which(taxon_means > 0)
  fit <- lm(log(taxon_vars[valid]) ~ log(taxon_means[valid]))
  slope <- coef(fit)[2]
  
  return(list(
    overall_zero = overall_zero_prop,
    mean_var = data.frame(mean = taxon_means, var = taxon_vars, disp = dispersion_index),
    slope = slope,
    avg_dispersion_index = mean(dispersion_index[taxon_means > 0], na.rm = TRUE)
  ))
}


results <- list()

for (i in seq_len(nrow(scenario_grid))) {
  pars <- scenario_grid[i, ]
  cat("Running scenario", pars$scenario_id, "...\n")
  
  # simulate
  W <- simulate_zipg_ext(n = 200, K = 30,
                         p = pars$p,
                         disp_shift = pars$disp_shift,
                         beta_scale = pars$beta_scale,
                         seed = 123 + i)
  
  # evaluate
  eval_res <- evaluate_counts(W)
  
  results[[i]] <- list(
    scenario_id = pars$scenario_id,
    params = pars,
    eval = eval_res
  )
}

# Extract summary table: overall zero proportion and slope
summary_df <- do.call(rbind, lapply(results, function(r) {
  data.frame(
    scenario_id = r$scenario_id,
    p = r$params$p,
    disp_shift = r$params$disp_shift,
    beta_scale = r$params$beta_scale,
    overall_zero = r$eval$overall_zero,
    slope = r$eval$slope,
    avg_dispersion_index = r$eval$avg_dispersion_index
  )
}))

print(summary_df)


make_main_effect_df <- function(df, factor_col) {
  df %>%
    group_by(.data[[factor_col]]) %>%
    summarise(
      mean_zero = mean(overall_zero),
      se_zero   = sd(overall_zero) / sqrt(n()),
      lower     = mean_zero - 1.96 * se_zero,
      upper     = mean_zero + 1.96 * se_zero,
      .groups = "drop"
    ) %>%
    rename(level = {{factor_col}})
}

main_p    <- make_main_effect_df(summary_df, "p")
main_disp <- make_main_effect_df(summary_df, "disp_shift")
main_beta <- make_main_effect_df(summary_df, "beta_scale")

print(main_p)
print(main_disp)
print(main_beta)


# Plot for zero-inflation p
a <- gg_p <- ggplot(main_p, aes(x = as.factor(level), y = mean_zero)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Zero-inflation probability (p)", y = "Overall zero proportion") +
  ggtitle("Effect of zero-inflation (p)") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  theme_minimal()

# Plot for dispersion shift
b <- gg_disp <- ggplot(main_disp, aes(x = as.factor(level), y = mean_zero)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Dispersion shift (log θ)", y = "Overall zero proportion") +
  ggtitle("Effect of dispersion on sparsity") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  theme_minimal()

# Plot for effect size
c <- gg_beta <- ggplot(main_beta, aes(x = as.factor(level), y = mean_zero)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Effect size multiplier (β)", y = "Overall zero proportion") +
  ggtitle("Effect of effect size on sparsity") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  theme_minimal()


final_plot <- (a | b | c) +
  plot_annotation(
    title = "Main Effects on Sparsity Across ZIPG Parameters",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )

# Print or save
print(final_plot)

# ==================== Mean-Variance Relationship ===============================

# Compute per-taxon DI and average DI for one simulated matrix W (n x K)
compute_dispersion_index <- function(W) {
  mu  <- colMeans(W)
  v   <- apply(W, 2, var)
  ok  <- which(mu > 0 & is.finite(mu) & is.finite(v))
  if (length(ok) == 0L) {
    return(list(DI_k = rep(NA_real_, ncol(W)), avg_DI = NA_real_))
  }
  DI_k <- v[ok] / mu[ok]
  avg_DI <- mean(DI_k, na.rm = TRUE)
  # Return DI aligned to taxa (NAs where mean==0)
  full_DI <- rep(NA_real_, ncol(W))
  full_DI[ok] <- DI_k
  list(DI_k = full_DI, avg_DI = avg_DI)
}

# Parameters
n_reps <- 20    # number of replicates per scenario
n <- 200; K <- 30  # fixed as per your design

# Storage for replicate-level averages
di_reps <- list()

for (i in seq_len(nrow(scenario_grid))) {
  pars <- scenario_grid[i, ]
  cat("Scenario", pars$scenario_id, " (p=", pars$p,
      ", disp_shift=", pars$disp_shift,
      ", beta_scale=", pars$beta_scale, ")\n", sep = "")
  
  for (r in 1:n_reps) {
    W <- simulate_zipg_ext(
      n = n, K = K,
      p = pars$p,
      disp_shift = pars$disp_shift,
      beta_scale = pars$beta_scale,
      seed = 5000 + 100*i + r
    )
    di <- compute_dispersion_index(W)
    di_reps[[length(di_reps) + 1L]] <- data.frame(
      scenario_id = pars$scenario_id,
      p = pars$p,
      disp_shift = pars$disp_shift,
      beta_scale = pars$beta_scale,
      replicate = r,
      avg_DI = di$avg_DI
    )
  }
}

di_reps_df <- do.call(rbind, di_reps)
head(di_reps_df)



summarise_main_effect <- function(df, factor_col) {
  df %>%
    group_by(.data[[factor_col]]) %>%
    summarise(
      mean_avg_DI = mean(avg_DI, na.rm = TRUE),
      se          = sd(avg_DI,  na.rm = TRUE) / sqrt(sum(!is.na(avg_DI))),
      lower       = mean_avg_DI - 1.96 * se,
      upper       = mean_avg_DI + 1.96 * se,
      .groups = "drop"
    ) %>%
    rename(level = {{factor_col}})
}

di_p_summary    <- summarise_main_effect(di_reps_df, "p")
di_disp_summary <- summarise_main_effect(di_reps_df, "disp_shift")
di_beta_summary <- summarise_main_effect(di_reps_df, "beta_scale")

di_p_summary
di_disp_summary
di_beta_summary


theme_set(theme_minimal(base_size = 13))

plot_main_ci <- function(df, xlab, ylab = "Average dispersion index (var/mean)", 
                         ylimits = NULL) {
  p <- ggplot(df, aes(x = factor(level), y = mean_avg_DI)) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.12) +
    labs(x = xlab, y = ylab) +
    theme_minimal()
  
  if (!is.null(ylimits)) {
    p <- p + ylim(ylimits)
  }
  return(p)
}
# Common y-limits
ylims <- range(
  c(di_p_summary$lower, di_p_summary$upper,
    di_disp_summary$lower, di_disp_summary$upper,
    di_beta_summary$lower, di_beta_summary$upper),
  na.rm = TRUE
)

# Individual plots
gg_di_p    <- plot_main_ci(di_p_summary,    "Zero-inflation probability (p)", ylimits = ylims)
gg_di_disp <- plot_main_ci(di_disp_summary, "Log-dispersion shift",           ylimits = ylims)
gg_di_beta <- plot_main_ci(di_beta_summary, "Effect size scale",              ylimits = ylims)

combined_di <- (gg_di_p | gg_di_disp | gg_di_beta) +
  plot_annotation(
    title = "Main effects on dispersion index",
    subtitle = "Zero inflation, dispersion shift, and covariate effect size"
  )

combined_di


# Compute mean and 95% CI
summarise_main_effect <- function(df, factor_col) {
  df %>%
    group_by(.data[[factor_col]]) %>%
    summarise(
      mean_DI = mean(avg_DI, na.rm = TRUE),
      se      = sd(avg_DI,  na.rm = TRUE) / sqrt(sum(!is.na(avg_DI))),
      lower   = mean_DI - 1.96 * se,
      upper   = mean_DI + 1.96 * se,
      .groups = "drop"
    ) %>%
    rename(level = {{factor_col}})
}

# Main effects
di_p_summary    <- summarise_main_effect(di_reps_df, "p")
di_disp_summary <- summarise_main_effect(di_reps_df, "disp_shift")
di_beta_summary <- summarise_main_effect(di_reps_df, "beta_scale")

print(di_p_summary)
print(di_disp_summary)
print(di_beta_summary)
