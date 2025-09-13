library(ZIPFA)
library(dplyr)
library(ggplot2)
library(tibble)

# ===================== Simulator for ZIPFA =====================================

simulate_zipfa_data <- function(n=200, p=30, K=3, tau=0.7,
                                m_mean=1.5e4, m_sd=2e3, sV=0.7, seed=1){
  set.seed(seed)
  U <- matrix(rnorm(n*K), n, K)
  V <- matrix(rnorm(p*K, sd=sV), p, K)
  L <- U %*% t(V)                     # log Lambda
  Lambda <- exp(L)
  M <- pmax(round(rnorm(n, m_mean, m_sd)), 1000)
  mu <- Lambda * matrix(M, n, p)
  P  <- 1 / (1 + Lambda^tau)
  
  W <- matrix(0L, n, p)
  for(i in 1:n){
    for(k in 1:p){
      if (runif(1) < P[i,k]) next
      W[i,k] <- rpois(1, mu[i,k])
    }
  }
  list(W=W, U=U, V=V, L=L, tau=tau, M=M)
}

# ===================== Parameter recovery ============================

n_values <- c(50, 100, 200, 400)
n_reps   <- 50
K <- 3
p <- 30

results   <- list()

for (n in n_values) {
  cat(sprintf("=== Sample size n = %d ===\n", n))
  metrics_list <- vector("list", n_reps)
  seeds <- sample(1e6, n_reps)
  
  for (rep_i in seq_len(n_reps)) {
    sd_i <- seeds[rep_i]
    sim <- simulate_zipfa_data(n = n, p = p, K = K, seed = sd_i)
    log_lambda_true <- sim$L
    W <- sim$W
    
    fit <- tryCatch({
      ZIPFA(W, k = K, tau = 0.7)
    }, error = function(e) {
      warning(sprintf("ZIPFA failed (n=%d, rep=%d, seed=%d): %s", n, rep_i, sd_i, conditionMessage(e)))
      return(NULL)
    })
    
    if (is.null(fit)) {
      metrics_list[[rep_i]] <- data.frame(FrobError = NA_real_, MSE = NA_real_)
      next
    }
    
    U_hat <- fit$Ufit[[length(fit$Ufit)]]
    V_hat <- fit$Vfit[[length(fit$Vfit)]]
    log_lambda_hat <- U_hat %*% t(V_hat)
    
    diff_mat <- log_lambda_hat - log_lambda_true
    mse_val  <- mean((diff_mat)^2)
    
    metrics_list[[rep_i]] <- data.frame(MSE = mse_val)
  }
  
  res_df <- dplyr::bind_rows(metrics_list)
  res_df$n <- n
  results[[as.character(n)]] <- res_df
}

metrics_all <- dplyr::bind_rows(results)

summary_df <- metrics_all %>%
  filter(!is.na(MSE)) %>%
  group_by(n) %>%
  summarise(
    MSE_mean  = mean(MSE),
    MSE_se    = sd(MSE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    MSE_lower  = MSE_mean - 1.96 * MSE_se,
    MSE_upper  = MSE_mean + 1.96 * MSE_se
  )

print(summary_df)



ggplot(summary_df, aes(x = factor(n), y = MSE_mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = MSE_lower, ymax = MSE_upper), width = 0.2) +
  labs(x = "Sample size (n)", y = "MSE of log Λ estimate",
       title = "MSE of log Λ recovery by sample size") +
  theme_minimal()


summary_df_1 <- tibble::tibble(
   n = c(50, 100, 200, 400),
   MSE_mean = c(96.5, 92.8, 91.4, 90.9),
   MSE_se   = c(1.71, 2.02, 1.92, 1.79),
   MSE_lower = c(93.2, 88.5, 87.7, 87.4),
   MSE_upper = c(99.8, 96.8, 95.2, 94.4)
 )
ggplot(summary_df_1, aes(x = factor(n), y = MSE_mean)) +
   geom_errorbar(aes(ymin = MSE_lower, ymax = MSE_upper),
                 width = 0.18, size = 0.5) +
   geom_point(size = 3) +
   labs(title = "MSE of log Λ recovery by sample size",
        x = "Sample size (n)", y = "MSE of log Λ estimate") +
   coord_cartesian(ylim = c(85, 100)) +
   theme_minimal(base_size = 13)
 