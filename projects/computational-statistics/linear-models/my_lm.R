# ============================================================
# Linear Models via Matrix Algebra (from scratch)
# Author: Mary Acheampong
# Purpose:
#   - Fit OLS using matrix formulas
#   - Produce coefficient table: estimate, SE, t, p-value, CI
# Notes:
#   - Designed for learning/portfolio use (base R, no dependencies)
# ============================================================

#' Summary for one coefficient using a t-test
#' @param est numeric, coefficient estimate
#' @param se numeric, standard error
#' @param df integer, degrees of freedom
#' @param alpha numeric, significance level (default 0.05)
t_test_summary <- function(est, se, df, alpha = 0.05) {
  t0 <- est / se
  p_value <- 2 * stats::pt(-abs(t0), df = df)
  t_crit <- stats::qt(1 - alpha / 2, df = df)
  lower <- est - t_crit * se
  upper <- est + t_crit * se
  
  c(est = est, se = se, t0 = t0, p.value = p_value, lower = lower, upper = upper)
}

#' Fit OLS regression using matrix algebra
#'
#' @param y numeric vector (n x 1)
#' @param X numeric matrix (n x p). Include an intercept column if desired.
#' @param col_names optional character vector of column names for X
#' @param alpha numeric, significance level for CI (default 0.05)
#' @param digits integer, rounding for printed table
#' @return list with coefficients table, beta_hat, vcov, sigma2, df_resid, fitted, resid
my_lm <- function(y, X, col_names = NULL, alpha = 0.05, digits = 5) {
  # Basic checks
  if (!is.numeric(y)) stop("y must be numeric.")
  y <- as.numeric(y)
  
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be numeric.")
  if (nrow(X) != length(y)) stop("nrow(X) must equal length(y).")
  
  n <- nrow(X)
  
  # Use QR to get rank (more stable than relying on p = ncol(X))
  qrX <- qr(X)
  p_rank <- qrX$rank
  if (p_rank < ncol(X)) {
    warning("X appears rank-deficient; results may be unstable. Consider removing collinear columns.")
  }
  
  XtX_inv <- solve(t(X) %*% X)
  beta_hat <- as.vector(XtX_inv %*% t(X) %*% y)
  
  fitted <- as.vector(X %*% beta_hat)
  resid <- y - fitted
  
  df_resid <- n - p_rank
  sigma2 <- as.numeric(t(resid) %*% resid / df_resid)
  
  vcov <- sigma2 * XtX_inv
  se <- sqrt(diag(vcov))
  
  # Build coefficient table
  rows <- Map(t_test_summary, beta_hat, se, MoreArgs = list(df = df_resid, alpha = alpha))
  coef_mat <- do.call(rbind, rows)
  
  # Formatting
  coef_df <- as.data.frame(coef_mat)
  num_cols <- c("est", "se", "t0", "lower", "upper")
  coef_df[num_cols] <- round(coef_df[num_cols], digits)
  coef_df$p.value <- format.pval(coef_df$p.value, eps = 1e-4, digits = 4)
  
  if (!is.null(col_names)) {
    if (length(col_names) != ncol(X)) stop("col_names must have length ncol(X).")
    rownames(coef_df) <- col_names
  }
  
  list(
    coefficients = coef_df,
    beta_hat = beta_hat,
    vcov = vcov,
    sigma2 = sigma2,
    df_resid = df_resid,
    fitted = fitted,
    resid = resid
  )
}

# -----------------------------
# Example (simulated data)
# -----------------------------
if (interactive()) {
  set.seed(123)
  
  n <- 80
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  X <- cbind(Intercept = 1, x1 = x1, x2 = x2)
  
  beta_true <- c(2, 1.5, -0.8)
  y <- as.vector(X %*% beta_true + rnorm(n, sd = 1))
  
  fit <- my_lm(y, X, col_names = colnames(X))
  print(fit$coefficients)
  
  # Compare to base R lm for sanity
  lm_fit <- lm(y ~ x1 + x2)
  print(summary(lm_fit)$coefficients)
}
