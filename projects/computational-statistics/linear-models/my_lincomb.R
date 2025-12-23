# ============================================================
# Linear Combinations / Contrasts for Linear Regression
# Author: Mary Acheampong
# Purpose:
#   Inference for linear combinations of coefficients:
#     H0: a'β = c
#   using t-tests and confidence intervals.
# ============================================================

#' Inference for linear combinations a' beta
#'
#' @param beta_hat numeric vector of estimated coefficients (p x 1)
#' @param vcov_mat variance-covariance matrix of beta_hat (p x p)
#' @param a matrix or vector defining linear combinations (k x p) or (p,)
#' @param c numeric scalar or vector length k (hypothesized values), default 0
#' @param df integer degrees of freedom for t reference
#' @param alpha numeric significance level for CI (default 0.05)
#' @param row_names optional names for each linear combination (length k)
#' @return data.frame with estimate, SE, t, p-value, CI
lincomb_test <- function(beta_hat, vcov_mat, a, c = 0, df, alpha = 0.05, row_names = NULL) {
  beta_hat <- as.numeric(beta_hat)
  vcov_mat <- as.matrix(vcov_mat)
  
  p <- length(beta_hat)
  if (nrow(vcov_mat) != p || ncol(vcov_mat) != p) stop("vcov_mat must be p x p.")
  
  # Ensure 'a' is a matrix (k x p)
  if (is.vector(a)) a <- matrix(a, nrow = 1)
  a <- as.matrix(a)
  
  if (ncol(a) != p) stop("a must have p columns (same length as beta_hat).")
  
  k <- nrow(a)
  
  # Handle c
  if (length(c) == 1) c <- rep(c, k)
  if (length(c) != k) stop("c must be scalar or length equal to number of rows in a.")
  
  est <- as.vector(a %*% beta_hat)
  se <- sqrt(diag(a %*% vcov_mat %*% t(a)))
  
  t0 <- (est - c) / se
  p_value <- 2 * stats::pt(-abs(t0), df = df)
  
  t_crit <- stats::qt(1 - alpha / 2, df = df)
  lower <- est - t_crit * se
  upper <- est + t_crit * se
  
  out <- data.frame(
    estimate = est,
    se = se,
    t = t0,
    p_value = p_value,
    lower = lower,
    upper = upper
  )
  
  if (!is.null(row_names)) {
    if (length(row_names) != k) stop("row_names must have length k.")
    rownames(out) <- row_names
  }
  
  out
}

# -----------------------------
# Example (safe, simulated data)
# -----------------------------
if (interactive()) {
  set.seed(7)
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  X <- cbind(Intercept = 1, x1 = x1, x2 = x2)
  
  beta_true <- c(1, 2, -1)
  y <- as.vector(X %*% beta_true + rnorm(n, sd = 1.5))
  
  # Use the my_lm() function you created in my_lm.R
  # If needed, source it:
  # source("projects/computational-statistics/linear-models/my_lm.R")
  
  fit <- my_lm(y, X, col_names = colnames(X))
  beta_hat <- fit$beta_hat
  vc <- fit$vcov
  df <- fit$df_resid
  
  # Example 1: test H0: (x1 + x2) = 0
  a1 <- c(0, 1, 1)
  res1 <- lincomb_test(beta_hat, vc, a = a1, c = 0, df = df, row_names = "x1 + x2 = 0")
  print(res1)
  
  # Example 2: two contrasts at once
  A <- rbind(
    c(0, 1, 0),  # x1 = 0
    c(0, 0, 1)   # x2 = 0
  )
  res2 <- lincomb_test(beta_hat, vc, a = A, c = 0, df = df, row_names = c("x1 = 0", "x2 = 0"))
  print(res2)
}
