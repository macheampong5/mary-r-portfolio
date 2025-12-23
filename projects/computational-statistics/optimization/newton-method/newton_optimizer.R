# ============================================================
# Newton's Method Optimizer (General Purpose)
# Author: Mary Acheampong
# Purpose:
#   Maximize an objective function by solving:
#     theta_{k+1} = theta_k - H(theta_k)^{-1} g(theta_k)
#   where g is the gradient and H is the Hessian of the objective.
#
# Notes:
#   - This is a general optimizer used later for MLE (Beta, Logistic, etc.)
#   - Assumes you provide g(theta) and H(theta)
# ============================================================

#' Newton optimizer for maximization
#'
#' @param theta0 numeric vector, initial guess
#' @param g function(theta) returning gradient vector
#' @param H function(theta) returning Hessian matrix
#' @param eps numeric, convergence tolerance on ||theta_{k+1} - theta_k||
#' @param max_iter integer, maximum number of iterations
#' @param verbose logical, print iteration diagnostics
#' @return list(estimates, hessian, iterations, convergence_criterion, converged)
my_newton <- function(theta0, g, H, eps = 1e-6, max_iter = 50, verbose = TRUE) {
  theta0 <- as.vector(theta0)
  crit <- eps + 1
  it <- 0
  while (crit > eps && it < max_iter) {
    it <- it + 1
    grad <- g(theta0)
    hess <- H(theta0)
    # Newton step: solve(H) %*% g
    step <- tryCatch(
      solve(hess, grad),
      error = function(e) stop("Hessian is singular / not invertible at current iterate.")
    )
    theta1 <- as.vector(theta0 - step)
    diff <- theta1 - theta0
    crit <- sqrt(drop(t(diff) %*% diff))
    if (verbose) {
      print(c(iteration = it, criterion = crit))
    }
    theta0 <- theta1
  }
  list(
    estimates = theta0,
    hessian = H(theta0),
    iterations = it,
    convergence_criterion = crit,
    converged = (crit <= eps)
  )
}

# -----------------------------
# Example 
# Maximize: f(t) = -(t-2)^2  (maximum at t = 2)
# g(t) = -2(t-2), H(t) = -2
# -----------------------------
if (interactive()) {
  f_demo <- function(t) - (t[1] - 2)^2
  g_demo <- function(t) c(-2 * (t[1] - 2))
  H_demo <- function(t) matrix(-2, 1, 1)
  
  out <- my_newton(theta0 = 0, g = g_demo, H = H_demo, eps = 1e-10, verbose = TRUE)
  print(out$estimates)
  print(out$converged)
}
