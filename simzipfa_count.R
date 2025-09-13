library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

# Parameters
n <- 200     # fixed sample size
p <- 100     # number of taxa
n_reps <- 20 # per scenario

tau_values   <- c(0.5, 0.7, 1.0)
sigma_values <- c(0.5, 0.7, 1.0)
K_values     <- c(2, 3, 5)

# compute dispersion index for one count matrix
compute_dispersion_index <- function(W) {
  means <- rowMeans(W)
  vars  <- apply(W, 1, var)
  di <- vars / means
  di[is.nan(di) | is.infinite(di)] <- NA
  return(mean(di, na.rm = TRUE))
}

# ====================== simulator ===============================

simulate_zipfa_data <- function(n = 200, p = 100, K = 3, tau = 0.7,
                                m_mean = 1.5e4, m_sd = 2e3, 
                                sigma_V = 0.7,
                                seed = 1) {
  set.seed(seed)
  
  # Sample latent factors and loadings
  U <- matrix(rnorm(n * K), nrow = n, ncol = K)
  V <- matrix(rnorm(p * K, sd = sigma_V), nrow = p, ncol = K)
  
  # Compute log lambda and Poisson mean
  log_Lambda <- U %*% t(V)
  Lambda <- exp(log_Lambda)
  
  # Compute structural zero probability
  P <- 1 / (1 + Lambda^tau)
  
  # Sample sequencing depths
  M <- pmax(round(rnorm(n, mean = m_mean, sd = m_sd)), 1000)
  
  # Simulate observed counts
  W <- matrix(0L, nrow = n, ncol = p)
  for (i in 1:n) {
    for (j in 1:p) {
      if (runif(1) > P[i, j]) {
        W[i, j] <- rpois(1, Lambda[i, j] * M[i])
      }
    }
  }
  
  return(list(W = W, L = log_Lambda, M = M, P = P))
}


results <- list()

grid <- expand.grid(tau = tau_values, sigma = sigma_values, K = K_values)

for (i in seq_len(nrow(grid))) {
  tau_i   <- grid$tau[i]
  sigma_i <- grid$sigma[i]
  K_i     <- grid$K[i]
  
  cat(sprintf("Running ZIPFA sim: tau = %.1f | sigma = %.1f | K = %d\n", tau_i, sigma_i, K_i))
  
  sparsity_vec <- numeric(n_reps)
  dispindex_vec <- numeric(n_reps)
  
  for (r in 1:n_reps) {
    sim <- simulate_zipfa_data(n = n, p = p, K = K_i, tau = tau_i, sigma_V = sigma_i, seed = r)
    W <- sim$W
    
    sparsity_vec[r]   <- mean(W == 0)
    dispindex_vec[r]  <- compute_dispersion_index(t(W))  # taxon-wise
  }
  
  results[[i]] <- data.frame(
    tau   = tau_i,
    sigma = sigma_i,
    K     = K_i,
    mean_zeros = mean(sparsity_vec),
    se_zeros   = sd(sparsity_vec) / sqrt(n_reps),
    mean_di    = mean(dispindex_vec),
    se_di      = sd(dispindex_vec) / sqrt(n_reps)
  )
}


zipfa_summary <- bind_rows(results)
zipfa_summary <- zipfa_summary %>%
  mutate(
    ci_lower_zeros = mean_zeros - 1.96 * se_zeros,
    ci_upper_zeros = mean_zeros + 1.96 * se_zeros,
    ci_lower_di    = mean_di    - 1.96 * se_di,
    ci_upper_di    = mean_di    + 1.96 * se_di
  )

# View result
print(zipfa_summary)


# Summarise by tau
zipfa_tau_plot <- zipfa_summary %>%
  group_by(tau) %>%
  summarise(
    mean_zeros = mean(mean_zeros),
    se_zeros = mean(se_zeros),
    ci_lower_zeros = mean_zeros - 1.96 * se_zeros,
    ci_upper_zeros = mean_zeros + 1.96 * se_zeros,
    mean_di = mean(mean_di),
    se_di = mean(se_di),
    ci_lower_di = mean_di - 1.96 * se_di,
    ci_upper_di = mean_di + 1.96 * se_di
  )

# Plot zero proportion vs τ
ggplot(zipfa_tau_plot, aes(x = factor(tau), y = mean_zeros)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower_zeros, ymax = ci_upper_zeros), width = 0.15) +
  labs(
    title = "Effect of zero-inflation parameter (τ) on sparsity",
    x = "Zero-inflation parameter (τ)",
    y = "Overall zero proportion"
  ) +
  theme_bw(base_size = 13)

zipfa_sigma_plot <- zipfa_summary %>%
  group_by(sigma) %>%
  summarise(
    mean_zeros = mean(mean_zeros),
    se_zeros = mean(se_zeros),
    ci_lower_zeros = mean_zeros - 1.96 * se_zeros,
    ci_upper_zeros = mean_zeros + 1.96 * se_zeros,
    mean_di = mean(mean_di),
    se_di = mean(se_di),
    ci_lower_di = mean_di - 1.96 * se_di,
    ci_upper_di = mean_di + 1.96 * se_di
  )

# Plot DI vs σ
ggplot(zipfa_sigma_plot, aes(x = factor(sigma), y = mean_di)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower_di, ymax = ci_upper_di), width = 0.15) +
  labs(
    title = "Effect of loading variance (σ) on dispersion",
    x = "Loading variance (σ)",
    y = "Average dispersion index"
  ) +
  theme_bw(base_size = 13)

zipfa_K_plot <- zipfa_summary %>%
  group_by(K) %>%
  summarise(
    mean_zeros = mean(mean_zeros),
    se_zeros = mean(se_zeros),
    ci_lower_zeros = mean_zeros - 1.96 * se_zeros,
    ci_upper_zeros = mean_zeros + 1.96 * se_zeros,
    mean_di = mean(mean_di),
    se_di = mean(se_di),
    ci_lower_di = mean_di - 1.96 * se_di,
    ci_upper_di = mean_di + 1.96 * se_di
  )

# Plot DI vs K
ggplot(zipfa_K_plot, aes(x = factor(K), y = mean_di)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower_di, ymax = ci_upper_di), width = 0.15) +
  labs(
    title = "Effect of number of latent factors (K) on dispersion",
    x = "Number of latent factors (K)",
    y = "Average dispersion index"
  ) +
  theme_bw(base_size = 13)

