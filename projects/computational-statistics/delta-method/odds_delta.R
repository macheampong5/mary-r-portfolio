# t.o function
t.o <- function(est, se, df, alpha = 0.05) {
  t0 = est / se
  pvalue = 2 * pt(-abs(t0), df)
  tcrit = qt(1 - alpha/2, df)
  lower = est - tcrit * se
  upper = est + tcrit * se
  c(est = est, se = se, t0 = t0, p.value = pvalue, lower = lower, upper = upper)
}
n <- 1016
x<- 416
p.hat <- x/n
odds.hat<- p.hat/(1-p.hat)
# Delta Method SE
se.odds <- sqrt( p.hat / ( n * (1 - p.hat)^3 ) )
# Output
result <- t.o(est = odds.hat, se = se.odds, df = Inf)
result
