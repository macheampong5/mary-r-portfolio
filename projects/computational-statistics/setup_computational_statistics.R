# ================================
# Setup script for Computational Statistics portfolio
# Author: Mary Acheampong
# Purpose: Create folder and file structure (no content)
# ================================

base_dir <- "projects/computational-statistics"

dirs <- c(
  base_dir,
  file.path(base_dir, "linear-models"),
  file.path(base_dir, "matrix-factorization"),
  file.path(base_dir, "simulation-and-efficiency")
)

files <- c(
  file.path(base_dir, "README.md"),
  file.path(base_dir, "linear-models", "my_lm.R"),
  file.path(base_dir, "linear-models", "my_anova.R"),
  file.path(base_dir, "linear-models", "my_lincomb.R"),
  file.path(base_dir, "matrix-factorization", "spectral_decomposition.R"),
  file.path(base_dir, "matrix-factorization", "cholesky_decomposition.R"),
  file.path(base_dir, "simulation-and-efficiency", "timing_comparisons.R"),
  file.path(base_dir, "simulation-and-efficiency", "vectorization_tricks.R")
)

# Create directories
for (d in dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Create empty files
for (f in files) {
  if (!file.exists(f)) file.create(f)
}

cat("Computational Statistics structure created successfully.\n")
