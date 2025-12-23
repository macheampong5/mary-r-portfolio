# ============================================================
# Spectral (Eigen) Decomposition Utilities
# Author: Mary Acheampong
# Purpose:
#   - Compute matrix inverse square root using eigen-decomposition
#   - Use it for whitening / multivariate standardization
# Notes:
#   Works best for symmetric positive definite (SPD) matrices like covariances.
# ============================================================

#' Inverse square root of an SPD matrix via eigen-decomposition
#' For SPD S: S = V diag(lambda) V^T, then S^{-1/2} = V diag(1/sqrt(lambda)) V^T
#'
#' @param S symmetric matrix (typically covariance)
#' @param tol small tolerance for eigenvalues near zero
#' @return list(S_inv_half, eigenvalues, eigenvectors)
inv_sqrt_spd_spectral <- function(S, tol = 1e-12) {
  S <- as.matrix(S)
  if (!isTRUE(all.equal(S, t(S)))) stop("S must be symmetric.")
  
  er <- eigen(S)
  V <- er$vectors
  lam <- er$values
  
  if (any(lam <= tol)) {
    stop("S must be positive definite (all eigenvalues > tol).")
  }
  
  S_inv_half <- V %*% diag(1 / sqrt(lam)) %*% t(V)
  
  list(S_inv_half = S_inv_half, eigenvalues = lam, eigenvectors = V)
}

#' Whiten data using a covariance matrix inverse square root
#' Returns Z = (Y - mean) %*% t(S^{-1/2})
#'
#' @param Y numeric matrix n x p
#' @return list(Z, center, S, S_inv_half)
whiten_data_spectral <- function(Y) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(Y)
  
  center <- colMeans(Y)
  Yc <- Y - matrix(center, n, p, byrow = TRUE)
  
  S <- stats::cov(Yc)  # covariance
  out <- inv_sqrt_spd_spectral(S)
  
  Z <- Yc %*% t(out$S_inv_half)
  
  list(Z = Z, center = center, S = S, S_inv_half = out$S_inv_half)
}

# -----------------------------
# Example (built-in data)
# -----------------------------
if (interactive()) {
  # Use iris numeric features (built-in)
  Y <- as.matrix(iris[, 1:4])
  
  w <- whiten_data_spectral(Y)
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
