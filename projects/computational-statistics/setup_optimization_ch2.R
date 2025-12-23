# ============================================================
# Setup script for Chapter 2: Optimization
# Author: Mary Acheampong
# Purpose: Create folder + empty file structure
# ============================================================

base_dir <- "projects/optimization"

dirs <- c(
  base_dir,
  file.path(base_dir, "newton-method"),
  file.path(base_dir, "finite-differences"),
  file.path(base_dir, "likelihood-optimization"),
  file.path(base_dir, "metaheuristics"),
  file.path(base_dir, "combinatorial-optimization")
)

files <- c(
  file.path(base_dir, "README.md"),
  
  # Newton + derivatives
  file.path(base_dir, "newton-method", "newton_optimizer.R"),
  file.path(base_dir, "finite-differences", "fd_gradient.R"),
  file.path(base_dir, "finite-differences", "fd_hessian.R"),
  
  # Likelihood
  file.path(base_dir, "likelihood-optimization", "beta_mle_newton.R"),
  file.path(base_dir, "likelihood-optimization", "logistic_mle_newton.R"),
  
  # Metaheuristics
  file.path(base_dir, "metaheuristics", "genetic_algorithm_mle.R"),
  file.path(base_dir, "metaheuristics", "simulated_annealing_examples.R"),
  
  # Combinatorial
  file.path(base_dir, "combinatorial-optimization", "pushup_problem_sa.R")
)

for (d in dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

for (f in files) {
  if (!file.exists(f)) file.create(f)
}

cat("Chapter 2 (Optimization) structure created.\n")
