# WRS Reference Test Suite Creator
# This script creates baseline outputs from v0.45 for backward compatibility testing

# Load current version
library(WRS)

# Set fixed seed for reproducibility
set.seed(12345)

# Generate standard test data
x <- rnorm(50)
y <- rnorm(50) + 0.5
x2 <- rnorm(30)
y2 <- rnorm(30) + 0.3
x3 <- rnorm(25)
y3 <- rnorm(25) - 0.2

# Matrix data for multivariate tests
m1 <- matrix(rnorm(100), ncol = 2)
m2 <- matrix(rnorm(100) + 0.5, ncol = 2)

cat("Creating reference outputs from WRS v0.45...\n")
cat("Testing", length(x), "observations per group\n\n")

# Create reference test suite
# Test ~100 key functions covering major functionality
reference_tests <- list(
  # Two-sample tests
  yuen = yuen(x, y, tr = 0.2),
  yuenbt = yuenbt(x, y, nboot = 100),
  pb2gen_mean = pb2gen(x, y, est = mean, nboot = 100),
  trimpb = trimpb(x, y, nboot = 100),

  # One-way ANOVA
  t1way = t1way(list(x, y, x2)),
  # t1waybt = t1waybt(list(x, y, x2), nboot = 100),  # Uncomment if needed

  # Location estimators
  hd_median = hd(x, q = 0.5),
  hd_q25 = hd(x, q = 0.25),
  hd_q75 = hd(x, q = 0.75),
  mest = mest(x),
  onestep = onestep(x),
  tmean = tmean(x, tr = 0.2),
  winmean = winmean(x, tr = 0.2),

  # Variance/scale
  winvar = winvar(x, tr = 0.2),
  bivar = bivar(x),

  # Correlation
  pbcor = pbcor(x, y),
  wincor = wincor(x, y, tr = 0.2),

  # Outlier detection
  out = out(x),
  outbox = outbox(x),

  # Utilities
  elimna = elimna(c(x, NA, y)),
  winval = winval(x, tr = 0.2),

  # Confidence intervals
  trimse = trimse(x, tr = 0.2),
  hdci = hdci(x, q = 0.5, nboot = 100),
  mestci = mestci(x, nboot = 100)

  # Add more key functions as needed...
  # Regression, ANCOVA, etc. can be added in subsequent tests
)

cat("✓ Created", length(reference_tests), "reference test results\n")

# Save reference outputs
saveRDS(reference_tests, file = "/home/mando/coding/R-Projects/WRS/reference-outputs.rds")
cat("✓ Saved to reference-outputs.rds\n\n")

# Display sample outputs
cat("Sample reference values:\n")
cat("  yuen$p.value:", reference_tests$yuen$p.value, "\n")
cat("  hd(median):", reference_tests$hd_median, "\n")
cat("  pbcor:", reference_tests$pbcor$cor, "\n\n")

cat("Reference test suite created successfully!\n")
cat("Run this script BEFORE any modifications to establish baseline.\n")
