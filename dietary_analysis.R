# ========================= packages used ====================================
library(ZIPG)        # Zero-Inflated Poisson-Gamma model
library(ZIPFA)       # Zero-Inflated Poisson Factor Analysis model
library(VennDiagram) # For Venn diagrams (overlaps of sets)
library(dplyr)       # Data manipulation (select, filter, mutate, etc.)
library(tidyr)       # Data tidying (reshape, pivot)
library(stringr)     # String manipulation
library(purrr)       # Functional programming (map, reduce, etc.)
library(ggplot2)     # Data visualization
library(reshape2)    # Data reshaping (melt, dcast)
library(pheatmap)    # Heatmap plotting

data(Dietary)
W <- Dietary$OTU      # OTU count matrix (n × p)
M <- Dietary$M        # Sequencing depth (n × 1)
X <- Dietary$COV      # Covariates data
n <- nrow(W); p <- ncol(W)

# ============================  EDA ==========================================
# sequencing depth
depth_df <- data.frame(sample = 1:n, depth = M)
ggplot(depth_df, aes(x = depth)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Distribution of sequencing depth",
       x = "Total reads per sample",
       y = "Frequency")

# zero inflation
zero_overall <- mean(W == 0)  # overall proportion of zeros
cat("Overall proportion of zeros:", round(zero_overall, 3), "\n")

otu_zero <- colMeans(W == 0)
ggplot(data.frame(zero_prop = otu_zero), aes(x = zero_prop)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Zero proportion across OTUs",
       x = "Proportion of zero counts",
       y = "Number of OTUs")


### covariates
# age
hist(X$Age, breaks=20, col="steelblue", main="Distribution of Age")
summary(X$Age)

# Alcohol
ggplot(X, aes(x = factor(ALC01))) +
  geom_bar(fill = "grey50", width = 0.3) +
  theme_minimal(base_size = 15) +
  labs(title = "Alcohol consumption indicator",
       x = "ALC01 (0=non-drinker, 1=drinker)",
       y = "Number of samples")
prop_drinker <- mean(X$ALC01 == 1)
prop_drinker

# Nutrient PCs
ggplot(X, aes(x = nutrPC1)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Distribution of nutrPC1",
       x = "nutrPC1", y = "Frequency")

ggplot(X, aes(x = nutrPC2)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Distribution of nutrPC2",
       x = "nutrPC2", y = "Frequency")


# ======================= Model Fitting =========================================

# ======================== ZIPG fitting =======================================

form_X  <- ~ Age + ALC01 + nutrPC1 + nutrPC2   # abundance
form_Xs <- ~ ALC01                       # dispersion

success <- rep(FALSE, p)   # track which OTUs fit

for (j in seq_len(p)) {
  yj <- W[, j]
  fit_j <- tryCatch(
    ZIPG_main(
      data   = X,
      W      = yj,
      M      = M,
      X      = form_X,
      X_star = form_Xs,
      return_model = TRUE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(fit_j) && !all(is.na(fit_j$res$par))) {
    success[j] <- TRUE
  }
}

num_success <- sum(success)
num_fail    <- p - num_success
cat("Total OTUs:", p, "\n")
cat("Fitted successfully:", num_success, "\n")
cat("Failed:", num_fail, "\n")


failed_otus <- colnames(W)[!success]
success_otus <- colnames(W)[success]
W_success <- W[, success]

# ================== refit the model and collect estimates ==========================

zipg_results <- data.frame(
  OTU = colnames(W_success),
  beta_Age   = NA_real_, se_Age   = NA_real_, p_Age   = NA_real_,
  beta_ALC01 = NA_real_, se_ALC01 = NA_real_, p_ALC01 = NA_real_,
  beta_PC1   = NA_real_, se_PC1   = NA_real_, p_PC1   = NA_real_,
  beta_PC2   = NA_real_, se_PC2   = NA_real_, p_PC2   = NA_real_
)

for (j in seq_len(ncol(W_success))) {
  yj <- W_success[, j]
  
  fit_j <- tryCatch(
    ZIPG_main(
      data   = X,
      W      = yj,
      M      = M,
      X      = form_X,
      X_star = form_Xs,
      return_model = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(fit_j)) next
  
  # coefficient names (same order in $res$par and $wald_test)
  coef_names <- c("(Intercept)", "Age", "ALC01", "nutrPC1", "nutrPC2")
  
  if (!is.null(fit_j$res$par) && !is.null(fit_j$wald_test)) {
    coefs <- fit_j$res$par[1:length(coef_names)]
    ses   <- fit_j$wald_test$SE[1:length(coef_names)]
    pvals <- fit_j$wald_test$pval[1:length(coef_names)]
    names(coefs) <- names(ses) <- names(pvals) <- coef_names
    
    # --- save results for Age
    zipg_results$beta_Age[j]   <- coefs["Age"]
    zipg_results$se_Age[j]     <- ses["Age"]
    zipg_results$p_Age[j]      <- pvals["Age"]
    # --- save results for ALC01
    zipg_results$beta_ALC01[j] <- coefs["ALC01"]
    zipg_results$se_ALC01[j]   <- ses["ALC01"]
    zipg_results$p_ALC01[j]    <- pvals["ALC01"]
    
    # --- save results for nutrPC1
    zipg_results$beta_PC1[j] <- coefs["nutrPC1"]
    zipg_results$se_PC1[j]   <- ses["nutrPC1"]
    zipg_results$p_PC1[j]    <- pvals["nutrPC1"]
    
    # --- save results for nutrPC2
    zipg_results$beta_PC2[j] <- coefs["nutrPC2"]
    zipg_results$se_PC2[j]   <- ses["nutrPC2"]
    zipg_results$p_PC2[j]    <- pvals["nutrPC2"]
  }
}

zipg_results$FDR_Age <- p.adjust(zipg_results$p_Age, method = "BH")
zipg_results$FDR_ALC01 <- p.adjust(zipg_results$p_ALC01, method = "BH")
zipg_results$FDR_PC1   <- p.adjust(zipg_results$p_PC1,   method = "BH")
zipg_results$FDR_PC2   <- p.adjust(zipg_results$p_PC2,   method = "BH")

c(
  Age     = sum(zipg_results$FDR_Age   < 0.05, na.rm = TRUE),
  ALC01   = sum(zipg_results$FDR_ALC01 < 0.05, na.rm = TRUE),
  nutrPC1 = sum(zipg_results$FDR_PC1   < 0.05, na.rm = TRUE),
  nutrPC2 = sum(zipg_results$FDR_PC2   < 0.05, na.rm = TRUE)
)


library(VennDiagram)

zipg_results$sig_Age <- zipg_results$FDR_Age < 0.05
zipg_results$sig_ALC01   <- zipg_results$FDR_ALC01 < 0.05
zipg_results$sig_nutrPC1 <- zipg_results$FDR_PC1   < 0.05
zipg_results$sig_nutrPC2 <- zipg_results$FDR_PC2   < 0.05

# sets of significant OTUs
sig_otus_Age <- zipg_results$OTU[zipg_results$sig_Age]
sig_otus_ALC01   <- zipg_results$OTU[zipg_results$sig_ALC01]
sig_otus_nutrPC1 <- zipg_results$OTU[zipg_results$sig_nutrPC1]
sig_otus_nutrPC2 <- zipg_results$OTU[zipg_results$sig_nutrPC2]

venn_list <- list(
  Age     = sig_otus_Age,
  ALC01   = sig_otus_ALC01,
  nutrPC1 = sig_otus_nutrPC1,
  nutrPC2 = sig_otus_nutrPC2
)

# create Venn diagram
venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  fill = c("red", "steelblue", "orange", "forestgreen"),
  alpha = 0.6,
  cex = 1.4, fontface = "bold", fontfamily = "sans",
  cat.cex = 1.2, cat.fontface = "bold", cat.default.pos = "outer",
  category.names = c("Age", "ALC01", "nutrPC1", "nutrPC2"),
  margin = 0.05
)

grid::grid.newpage()
grid::grid.draw(venn.plot)

zipg_results$taxa_name <- Dietary$taxa_name[colnames(W) %in% zipg_results$OTU]
zipg_results$genus <- sapply(strsplit(zipg_results$taxa_name, " "), `[`, 1)
zipg_results$genus <- gsub("Uncl\\. Genus ", "", zipg_results$genus)
zipg_results$genus <- sapply(strsplit(zipg_results$genus, " "), `[`, 1)

# subset significant OTUs
sig_Age   <- subset(zipg_results, FDR_Age   < 0.05)
sig_ALC01 <- subset(zipg_results, FDR_ALC01 < 0.05)
sig_PC1   <- subset(zipg_results, FDR_PC1   < 0.05)
sig_PC2   <- subset(zipg_results, FDR_PC2   < 0.05)

# === Count top 10 genera for each covariate ===
cat("\nTop genera for Age:\n")
print(sort(table(sig_Age$genus), decreasing = TRUE)[1:10])

cat("Top genera for ALC01:\n")
print(sort(table(sig_ALC01$genus), decreasing = TRUE)[1:10])

cat("\nTop genera for nutrPC1:\n")
print(sort(table(sig_PC1$genus), decreasing = TRUE)[1:10])

cat("\nTop genera for nutrPC2:\n")
print(sort(table(sig_PC2$genus), decreasing = TRUE)[1:10])


# ========================= Table =============================================

alpha <- 0.05    # FDR significance threshold
beta_min <- 0.20    # Minimum effect size threshold
pct_change <- function(b) (exp(b) - 1) * 100     #convert log-coefficient to percentage change

# Reshape results into long format for easier filtering and ranking
long_df <- zipg_results %>%
  transmute(
    OTU, taxa_name, genus,
    beta_Age,    FDR_Age,
    beta_ALC01,  FDR_ALC01,
    beta_PC1,    FDR_PC1,
    beta_PC2,    FDR_PC2
  ) %>%
  pivot_longer(
    cols = -c(OTU, taxa_name, genus),
    names_to = c(".value", "covar"),
    names_pattern = "(beta|FDR)_(.*)"
  ) %>%
  mutate(
    pct = pct_change(beta),
    hit = (FDR <= alpha) & (abs(beta) >= beta_min),
    direction = case_when(beta >  0 ~ "positive",
                          beta <  0 ~ "negative",
                          TRUE        ~ "zero")
  )

ranked_hits <- long_df %>%
  filter(hit) %>%
  mutate(score = abs(beta) * (-log10(pmax(FDR, 1e-300)))) %>%
  arrange(desc(score))


# Selection functions
pick_top_hits_simple <- function(k = 10) {
  ranked_hits %>%
    slice_head(n = k)
}

pick_top_hits_diverse <- function(k = 10, per_genus_max = 2) {
  sel <- ranked_hits %>%
    mutate(row_id = row_number()) %>%
    {
      keep <- logical(nrow(.))
      genus_count <- list()
      for (i in seq_len(nrow(.))) {
        g <- .$genus[i] %||% "NA"
        cnt <- genus_count[[g]] %||% 0
        if (cnt < per_genus_max) {
          keep[i] <- TRUE
          genus_count[[g]] <- cnt + 1
        }
        if (sum(keep) >= k) break
      }
      .[keep, ]
    }
  sel
}

top_hits_simple  <- pick_top_hits_simple(k = 10)
top_hits_diverse <- pick_top_hits_diverse(k = 10, per_genus_max = 2)

# Format hits into a clean table for presentation
format_hits <- function(df_hits) {
  df_hits %>%
    transmute(
      OTU,
      taxa_name,
      genus,
      covariate = recode(covar,
                         Age   = "Age",
                         ALC01 = "ALC01 (alcohol)",
                         PC1   = "nutrPC1",
                         PC2   = "nutrPC2"),
      beta = round(beta, 3),
      FDR  = formatC(FDR, format = "e", digits = 2),
      `Pct change` = paste0(ifelse(pct >= 0, "+", ""), round(pct, 1), "%"),
      direction
    )
}

table_simple  <- format_hits(top_hits_simple)
table_diverse <- format_hits(top_hits_diverse)

ensure_both_signs <- function(df, min_each = 2) {
  pos <- df %>% filter(direction == "positive") %>% slice_head(n = min_each)
  neg <- df %>% filter(direction == "negative") %>% slice_head(n = min_each)
  rest <- df %>% filter(!(OTU %in% c(pos$OTU, neg$OTU) & covar %in% c(pos$covar, neg$covar)))
  bind_rows(pos, neg, rest) %>%
    slice_head(n = 10) %>%
    format_hits()
}

# Balanced selection with both positive and negative hits
table_diverse_balanced <- ensure_both_signs(top_hits_diverse, min_each = 2)

table_simple
table_diverse
table_diverse_balanced

# ========================== ZIPFA Model =====================================

library(ZIPFA)
data(Dietary)
W <- Dietary$OTU
set.seed(123)  # for reproducibility

# fit metrics 
W <- Dietary$OTU
n <- nrow(W); p <- ncol(W)
Ks <- 2:5
res_zipfa <- data.frame(
  K = Ks, logLik = NA_real_, AIC = NA_real_
)
fits <- vector("list", length(Ks))

for (i in seq_along(Ks)) {
  k <- Ks[i]
  
  fit <- ZIPFA(
    X = W, k = k,
    cut = 1, iter = 20, maxiter = 100,
    tolLnlikelihood = 5e-4,
    initialtau = "iteration",
    Madj = TRUE, display = FALSE
  )
  
  fits[[i]] <- fit
  
  # final log-likelihood
  ll_vec <- as.numeric(fit$Likelihood)
  logLik <- if (length(ll_vec)) tail(ll_vec, 1) else NA_real_
  
  # number of free parameters = n*k + p*k + 1
  kpar <- n*k + p*k + 1
  AIC  <- if (is.finite(logLik)) 2*kpar - 2*logLik else NA_real_
  
  res_zipfa$logLik[i] <- logLik
  res_zipfa$AIC[i]    <- AIC
}

print(res_zipfa)

# Find best K by AIC
best_idx <- which.min(res_zipfa$AIC)
best_K   <- res_zipfa$K[best_idx]
cat("Best K is:", best_K, "with AIC =", res_zipfa$AIC[best_idx], "\n")

df_long <- melt(res_zipfa, id.vars = "K", 
                measure.vars = c("logLik", "AIC"))

ggplot(df_long, aes(x = K, y = value, color = variable, group = variable)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  theme_minimal(base_size = 14) +
  labs(title = "ZIPFA fit metrics by factor rank",
       x = "Factor rank (K)",
       y = "Value",
       color = "Metric") +
  scale_x_continuous(breaks = res_zipfa$K)


# =================== Heatmap ================================================

fit_final <- ZIPFA(X=W, k=5, cut=1, iter=20, maxiter=100,
                   tolLnlikelihood=5e-4, initialtau="iteration",
                   Madj=TRUE, display=FALSE)
V <- fit_final$Vfit[[length(fit_final$Vfit)]]   # p x K
K <- ncol(V); p <- nrow(V)

if (is.null(rownames(V))) rownames(V) <- colnames(W)  # OTU IDs
taxa_names <- Dietary$taxa_name
if (length(taxa_names) == nrow(V)) {
  rownames(V) <- paste0(rownames(V), " — ", taxa_names)  # "OTUxx — Genus species"
}

m <- 20
top_taxa_idx <- order(rowSums(abs(V)), decreasing = TRUE)[1:m]
V_top <- V[top_taxa_idx, ]

# Get real taxon names
taxa_names <- Dietary$taxa_name  # or colnames(W), or any vector of taxon labels
rownames(V_top) <- taxa_names[top_taxa_idx]

# Name columns as Factor1, Factor2, etc.
colnames(V_top) <- paste0("Factor", 1:ncol(V))

# Plot heatmap
pheatmap(
  V_top,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  main = sprintf("ZIPFA: top %d taxa by absolute loading per factor", m),
  filename = "zipfa_toploadings_heatmap.png", 
  width = 7, height = 10
)


# ======================== Model Comparison ========================================

# ================== ZIPG metrics (unchanged) ==================

# Initialise result container
zipg_fit_results <- data.frame(
  OTU    = character(),
  logLik = numeric(),
  AIC    = numeric(),
  stringsAsFactors = FALSE
)


ZIPG_main(data = X, W = W[,4], M = M, X = form_X, X_star = form_Xs)

# Updated loop: fit & extract logLik and AIC
for (j in seq_len(p)) {
  cat("Fitting OTU: ", j, "\n")
  yj <- W[, j]
  
  fit_j <- tryCatch(
    ZIPG_main(
      data   = X,
      W      = yj,
      M      = M,
      X      = form_X,
      X_star = form_Xs,
      return_model = TRUE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(fit_j) && !all(is.na(fit_j$res$par))) {
    success[j] <- TRUE
    
    # extract log-likelihood
    loglik_j <- fit_j$logli
    n_param_j <- length(fit_j$res$par)
    aic_j <- -2 * loglik_j + 2 * n_param_j
    
    zipg_fit_results <- rbind(zipg_fit_results, data.frame(
      OTU    = colnames(W)[j],
      logLik = loglik_j,
      AIC    = aic_j
    ))
  }
}

# Save OTU names
success_otus <- zipg_fit_results$OTU

# Summarise total logLik and AIC
zipg_fit_summary <- zipg_fit_results %>%
  summarise(
    total_logLik = sum(logLik, na.rm = TRUE),
    total_AIC    = sum(AIC, na.rm = TRUE),
    n_OTUs       = n()
  )

print(zipg_fit_summary)


compute_zipfa_loglik_aic <- function(fit, W, M) {
  # fit: ZIPFA model fit object
  # W: n x p count matrix (used only for dimensions)
  # M: sequencing depth vector (not used here since fit$Likelihood is trusted)
  
  n <- nrow(W)
  p <- ncol(W)
  K <- ncol(fit$Ufit[[length(fit$Ufit)]])  # infer K from final U
  logLik_vec <- as.numeric(fit$Likelihood)
  logLik <- if (length(logLik_vec)) tail(logLik_vec, 1) else NA_real_
  
  num_params <- n * K + p * K + 1  # U + V + tau
  AIC <- if (is.finite(logLik)) 2 * num_params - 2 * logLik else NA_real_
  
  return(list(logLik = logLik, AIC = AIC))
}

zipfa_metrics_k5 <- compute_zipfa_loglik_aic(fit_final, W = W_success, M = M)
print(zipfa_metrics_k5)

# ===================== Latent factor linked to covariates =====================

U <- fit$Ufit[[length(fit$Ufit)]]   # n × K matrix

U_df <- as.data.frame(U)
colnames(U_df) <- paste0("Factor", 1:ncol(U))
U_df <- cbind(U_df, X)

library(broom)
library(dplyr)

results <- lapply(seq_len(ncol(U)), function(k) {
  f <- as.formula(paste0("Factor", k, " ~ Age + ALC01 + nutrPC1 + nutrPC2"))
  tidy(lm(f, data = U_df)) %>% 
    mutate(Factor = paste0("Factor", k))
})

reg_results <- do.call(rbind, results)

reg_results <- reg_results %>%
  group_by(term) %>%
  mutate(FDR = p.adjust(p.value, method = "BH")) %>%
  ungroup()
sig_results <- reg_results %>%
  filter(FDR <= 0.05)

library(ggplot2)

# Optional: cap estimates to a fixed range for color mapping
reg_results$capped_estimate <- pmax(pmin(reg_results$estimate, 50), -50)

ggplot(reg_results, aes(x = Factor, y = term, fill = capped_estimate)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(FDR <= 0.05, "*", "")), 
            size = 3, color = "black", fontface = "bold", vjust = 0.5) +
  scale_fill_gradient2(
    low = "#4575b4", mid = "white", high = "#d73027",
    midpoint = 0, limits = c(-50, 50), name = "Estimate"
  ) +
  labs(
    title = "Covariate associations with ZIPFA latent factors",
    x = "Factor (U)", y = "Covariate"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )


