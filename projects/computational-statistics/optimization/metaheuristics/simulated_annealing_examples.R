# Simulated Annealing (SA) examples using optim(method="SANN")

# Example 1: minimize a simple 1D function
if (interactive()) {
  set.seed(0)
  f1 <- function(t) (t[1] - 2)^2 + sin(5*t[1])^2
  sann1 <- optim(
    par = 0,
    fn = f1,
    method = "SANN",
    control = list(fnscale = 1, maxit = 10000, temp = 100, trace = FALSE, REPORT = 1000)
  )
  c(est = sann1$par, objective = sann1$value)
}

# Example 2: minimize a 2D function
if (interactive()) {
  set.seed(0)
  f2 <- function(t) (t[1] - 1)^2 + (t[2] + 2)^2 + 0.2*sin(6*t[1]) + 0.2*cos(6*t[2])
  sann2 <- optim(
    par = c(0, 0),
    fn = f2,
    method = "SANN",
    control = list(fnscale = 1, maxit = 50000, temp = 1000, trace = FALSE, REPORT = 1000)
  )
  c(sann2$par, sann2$value)
}
