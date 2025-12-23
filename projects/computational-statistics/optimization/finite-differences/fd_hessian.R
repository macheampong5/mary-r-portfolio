#### Finite-difference Hessian
H.fd <- function(f, t, ee = (1e-16)^(1/4)) {
  ## Get number of parameters
  p <- length(t)
  ## Initialize Hessian matrix (all zeros)
  H <- matrix(0, p, p)
  # Loop through rows of Hessian
  for (i in 1:p) {
    ## Create vector for parameter i
    ei <- rep(0, p); ei[i] <- ee
    ## Diagonal elements: second derivative using central difference
    ## Formula: [f(t+e) - 2f(t) + f(t-e)] / ^2 ; ee = episilon
    H[i, i] <- (f(t + ei) - 2*f(t) + f(t - ei)) / (ee^2)
    ## Loop through columns of Hessian
    for (j in 1:p) {
      ## Only compute upper triangle (Hessian is symmetric)
      if (i < j) {
        # Create vector for parameter j
        ej <- rep(0, p); ej[j] <- ee
        # Off-diagonal elements: mixed partial derivatives
        ##Formula:[f(t+e_i+e_j)-f(t+e_i-e_j)-f(t-e_i+e_j)+f(t-e_i-e_j)]/(4e^2)
        H[i, j] <- (f(t + ei + ej) - f(t + ei - ej) - f(t - ei + ej) +
                      f(t - ei - ej)) / (4 * ee^2)
        # Use symmetry to fill lower triangle
        H[j, i] <- H[i, j] # symmetry
      }
    }
  }
  # Return the approximate Hessian matrix
  H
}



### Example (simple quadratic function)
if (interactive()) {
  f <- function(t) (t[1]-3)^2 + (t[2]+1)^2
  t0 <- c(0,0)
  H.fd(f, t0)
}

