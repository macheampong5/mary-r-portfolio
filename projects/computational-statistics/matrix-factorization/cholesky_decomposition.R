# ============================================================
# Cholesky Decomposition Utilities
# Author: Mary Acheampong
# Purpose:
#   - Compute matrix inverse square root using Cholesky factor
#   - Use it for whitening / multivariate standardization
# Notes:
#   Requires SPD matrix (e.g., covariance).
# ============================================================

#' Inverse square root of an SPD matrix via Cholesky
#' If S = R^T R (chol gives upper-triangular R), then S^{-1/2} = (R^T)^{-1}
#'
#' @param S symmetric positive definite matrix
#' @return list(S_inv_half, R)
inv_sqrt_spd_chol <- function(S) {
  S <- as.matrix(S)
  if (!isTRUE(all.equal(S, t(S)))) stop("S must be symmetric.")
  
  R <- chol(S)              # S = R^T R (R upper triangular)
  S_inv_half <- solve(t(R)) # (R^T)^{-1}
  
  list(S_inv_half = S_inv_half, R = R)
}

#' Whiten data using Cholesky-based inverse square root
#'
#' @param Y numeric matrix n x p
#' @return list(Z, center, S, S_inv_half)
whiten_data_chol <- function(Y) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(Y)
  
  center <- colMeans(Y)
  Yc <- Y - matrix(center, n, p, byrow = TRUE)
  
  S <- stats::cov(Yc)
  out <- inv_sqrt_spd_chol(S)
  
  Z <- Yc %*% t(out$S_inv_half)
  
  list(Z = Z, center = center, S = S, S_inv_half = out$S_inv_half)
}

# -----------------------------
# Example (built-in data)
# -----------------------------
if (interactive()) {
  Y <- as.matrix(iris[, 1:4])
  
  w <- whiten_data_chol(Y)
  Z <- w$Z
  S <- w$S
  S_inv_half <- w$S_inv_half
  p <- ncol(Y)
  
  # Check 1: (S^{-1/2})^T (S^{-1/2}) = S^{-1}
  check1 <- all.equal(t(S_inv_half) %*% S_inv_half, solve(S), check.attributes = FALSE)
  print(check1)
  
  # Check 2: whitening => S^{-1/2} S (S^{-1/2})^T = I
  check2 <- all.equal(S_inv_half %*% S %*% t(S_inv_half), diag(p), check.attributes = FALSE)
  print(check2)
  
  # Empirical covariance of whitened data ~ I
  print(round(cov(Z), 4))
}
