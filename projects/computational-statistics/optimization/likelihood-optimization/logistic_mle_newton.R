source("projects/optimization/newton-method/newton_optimizer.R")
source("projects/optimization/finite-differences/fd_gradient.R")
source("projects/optimization/finite-differences/fd_hessian.R")
dat <- matrix(c(
  1, 1, 14, 50,
  1, 0, 10, 40,
  0, 1, 14, 100,
  0, 0, 27, 159,
  1, 1, 28, 59,
  1, 0, 7, 35,
  0, 1, 19, 75,
  0, 0, 27,200
), ncol = 4, byrow = TRUE)
x1 <- dat[,1] # theft type
x2 <- dat[,2] # prior arrest
y <- dat[,3] # sent to prison
ny <- dat[,4] # not sent to prison
n <- y + ny # group totals
X <- cbind(1,x1,x2)
lm_y <- y/n
start_vals <- my.lm(lm_y, X)$est
start_vals

loglik_logit <- function(beta) {
  eta <- as.vector(X %*% beta) # linear predictor
  p <- exp(eta) / (1 + exp(eta)) # logistic transform
  sum(y * log(p) + (n - y) * log(1 - p))
}
g_beta <- function(b) g.fd(loglik_logit, b)
H_beta <- function(b) H.fd(loglik_logit, b)
## Newton solve
newton_out <- my_newton(start_vals, g = g_beta, H = H_beta)
## Standard errors
se <- sqrt(diag(solve(-newton_out$hessian)))
## Apply t.o() to each parameter
final_results <- Map(t.o,
                     est = newton_out$estimates,
                     se = se,
                     MoreArgs = list(df = Inf))
## Combine into a data frame
final_output <- as.data.frame(do.call(rbind, final_results))
## Format p-values
final_output$p.value <- format.pval(final_output$p.value, eps = 1e-5, digits = 5)
final_output

loglik_logit(newton_out$estimates)
