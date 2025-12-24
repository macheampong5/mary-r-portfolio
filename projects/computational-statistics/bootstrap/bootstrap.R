# ------------------------------------------------------------
# Resampling and Inference Utilities
# Author: Mary Acheampong
# Purpose: Generalized functions for inference and resampling
# ------------------------------------------------------------
samp.o <- function(theta.star, alpha = 0.05) {
  est <- mean(theta.star)
  se  <- sd(theta.star)
  ci  <- quantile(theta.star, probs = c(alpha/2, 1-alpha/2))
  c(est = est, se = se, lower = ci[1], upper = ci[2])
}

my.boot <- function(x, stat, B = 5000,
                    type = c("nonparametric","parametric"),
                    rdist = NULL, args = list()){
  type <- match.arg(type)
  n <- length(x)
  theta.star <- numeric(B)
  
  if (type == "nonparametric") {
    for (b in 1:B) {
      idx <- sample(1:n, n, replace = TRUE)
      x.star <- x[idx,]
      theta.star[b] <- stat(x.star)
    }
  }
  
  if (type == "parametric") {
    for (b in 1:B) {
      x.star <- do.call(rdist, c(list(n), args))
      theta.star[b] <- stat(x.star)
    }
  }
  
  return(theta.star)
}

# -------------------------------
# Simple examples (usage only)
# -------------------------------
if (interactive()) {
  
  # Example: t.o
  t.o(est = 1.5, se = 0.3, df = 100)
  
  # Example: bootstrap mean
  x <- matrix(rnorm(100), ncol = 1)
  stat <- function(z) mean(z)
  
  theta.star <- my.boot(x, stat, B = 1000)
  samp.o(theta.star)
  
}
