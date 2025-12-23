# ============================================================
# General Linear Test (GLT) / ANOVA via Hat Matrices
# Author: Mary Acheampong
# Purpose:
#   Compare a full linear model vs a reduced (nested) model using
#   the general linear F-test:
#     H0: Reduced model is adequate
#     Ha: Full model provides a significantly better fit
# ============================================================

#' Hat matrix: H = X (t(X)*X)^(-1)*t(X)
#' @param X numeric matrix (n x p)
#' @return numeric matrix (n x n)
hat_matrix <- function(X) {
  X <- as.matrix(X)
  X %*% solve(t(X) %*% X) %*% t(X)
}

#' General Linear Test (nested models) using hat matrices
#'
#' @param y numeric vector length n
#' @param X_full numeric matrix (n x p_full)
#' @param X_reduced numeric matrix (n x p_reduced), must be nested in X_full
#' @return list with ANOVA-style table + test stats
glt_anova <- function(y, X_full, X_reduced) {
  y <- as.numeric(y)
  X_full <- as.matrix(X_full)
  X_reduced <- as.matrix(X_reduced)
  
  n <- length(y)
  if (nrow(X_full) != n) stop("nrow(X_full) must equal length(y).")
  if (nrow(X_reduced) != n) stop("nrow(X_reduced) must equal length(y).")
  
  # Hat matrices
  Hf <- hat_matrix(X_full)
  Hr <- hat_matrix(X_reduced)
  
  # Degrees of freedom
  rank_full <- qr(X_full)$rank
  rank_reduced <- qr(X_reduced)$rank
  
  df_model <- rank_full - rank_reduced
  df_error <- n - rank_full
  if (df_model <= 0) stop("Reduced model must be nested with fewer parameters than full model.")
  
  # Sums of squares
  # SSD = extra sum of squares explained by full model beyond reduced
  SSD <- as.numeric(t(y) %*% (Hf - Hr) %*% y)
  
  # SSE_full = residual SS under full model
  SSE_full <- as.numeric(t(y) %*% (diag(n) - Hf) %*% y)
  
  MS_model <- SSD / df_model
  MS_error <- SSE_full / df_error
  
  F_value <- MS_model / MS_error
  p_value <- stats::pf(F_value, df_model, df_error, lower.tail = FALSE)
  
  table <- data.frame(
    Source = c("Model (Full vs Reduced)", "Error"),
    Df     = c(df_model, df_error),
    SumSq  = c(SSD, SSE_full),
    MeanSq = c(MS_model, MS_error),
    F      = c(F_value, NA_real_),
    p_value = c(p_value, NA_real_)
  )
  
  list(
    table = table,
    F_value = F_value,
    df1 = df_model,
    df2 = df_error,
    p_value = p_value
  )
}

# -----------------------------
# Example (safe, simulated data)
# -----------------------------
if (interactive()) {
  set.seed(42)
  n <- 60
  x <- rnorm(n)
  g <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  
  # Full model: intercept + x + group
  X_full <- model.matrix(~ x + g)
  
  # Reduced model: intercept + x
  X_reduced <- model.matrix(~ x)
  
  # Generate y with a group effect
  beta <- c(2, 1.2, 0.5, -0.3)  # (Intercept, x, gB, gC)
  y <- as.vector(X_full %*% beta + rnorm(n, sd = 1))
  
  out <- glt_anova(y, X_full, X_reduced)
  print(out$table)
  
  # Compare to base R anova()
  print(anova(lm(y ~ x), lm(y ~ x + g)))
}
