# Requires the GA package
# install.packages("GA")
library(GA)

# Source the logistic likelihood setup so we reuse X, y, n, and loglik_logit()
source("projects/optimization/likelihood-optimization/logistic_mle_newton.R")
set.seed(0)
ga_out <- ga(type = "real-valued",
             fitness = loglik_logit,
             lower = c(-5, -1, -1),
             upper = c( 0, 2, 2),
             popSize = 50,
             maxiter = 500,
             pcrossover = 0.9,
             pmutation = 0.05,
             elitism = 2,
             monitor = FALSE)

summary(ga_out)

# Best solution found
ga_out@solution
ga_out@fitnessValue


########  
# Simple tuning demo (small grid)
grid <- expand.grid(
  popSize = c(50, 100),
  maxiter = c(200, 500),
  pcrossover = c(0.7, 0.9),
  pmutation = c(0.05, 0.10)
)

results <- data.frame()

set.seed(0)
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  oo <- ga(type="real-valued", fitness=loglik_logit,
           lower=c(-5,-1,-1), upper=c(0,2,2),
           popSize=g$popSize, maxiter=g$maxiter,
           pcrossover=g$pcrossover, pmutation=g$pmutation,
           elitism=2, monitor=FALSE)
  
  results <- rbind(results, cbind(g, bestFitness=oo@fitnessValue))
}

results[order(-results$bestFitness), ]

