source("projects/optimization/newton-method/newton_optimizer.R")


### t.o function
t.o = function(est, se, df, alpha = 0.05) {
  t0 = est / se
  pvalue = 2 * pt(-abs(t0), df)
  tcrit = qt(1 - alpha/2, df)
  lower = est - tcrit * se
  upper = est + tcrit * se
  ## Return
  c(est = est, se = se, t0 = t0, p.value = pvalue, lower = lower, upper = upper)
}

### Example data (simulated Beta sample)
if (interactive()) {
  set.seed(123)
  x <- rbeta(200, shape1 = 1.2, shape2 = 3.4)
}


# Method of Moments starting alpha
xbar <- mean(x)
alpha0 <- 2 * xbar / (1 - xbar)
theta0 <- c(alpha0, 2)
## Summaries
S1 <- sum(log(x))
S2 <- sum(log(1-x))
n <- length(x)
## Gradient (score) — returns a vector c(g_alpha, g_beta)
g_beta <- function(theta) {
  alpha <- theta[1]; beta <- theta[2]
  c(
    n*(digamma(alpha+beta) - digamma(alpha)) + S1,
    n*(digamma(alpha+beta) - digamma(beta)) + S2
  )
}
## Hessian — returns a 2x2 matrix
H_beta <- function(theta) {
  alpha <- theta[1]; beta <- theta[2]
  matrix(
    c(
      n*(trigamma(alpha+beta) - trigamma(alpha)), n*trigamma(alpha+beta),
      n*trigamma(alpha+beta), n*(trigamma(alpha+beta) - trigamma(beta))
    ),
    nrow = 2, byrow = TRUE
  )
}
loglik_beta <- function(theta) {
  alpha <- theta[1]; beta <- theta[2]
  n * (lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta)) +
    (alpha) * S1 + (beta) * S2}

## Run Newton's method
result <- my_newton(
  theta0 = theta0,
  g = g_beta,
  H = H_beta,
  eps = 1e-8
)
result

### post-processing using t.o
se <- sqrt(diag(solve(-result$hessian))) ## standard errors
## Apply t.o to each parameter
final_results <- Map(t.o,
                     est = result$estimates,
                     se = se,
                     MoreArgs = list(df = Inf))
## Convert to nice output format
final_output <-as.data.frame(do.call(rbind, final_results))
final_output
