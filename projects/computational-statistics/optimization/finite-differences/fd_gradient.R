### Finite-difference gradient
g.fd <- function(f, t, ee = (1e-16)^(1/4)) {
  ## Get number of parameters
  p <- length(t)
  ## Initialize gradient vector (all zeros)
  g <- numeric(p)
  # Loop through each parameter
  for (i in 1:p) {
    ## Create vector: zero everywhere except small step at position i
    ei <- rep(0, p); ei[i] <- ee
    # Central difference formula: [f(t+e) - f(t-e)] / (2e)
    g[i] <- (f(t + ei) - f(t - ei)) / (2 * ee)
  }
  ## Return the approximate gradient vector
  g
}
### Example (simple quadratic function)
if (interactive()) {
  f <- function(t) (t[1]-3)^2 + (t[2]+1)^2
  t0 <- c(0,0)
  g.fd(f, t0)
}
