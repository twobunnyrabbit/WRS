#!/usr/bin/env Rscript

# Run backward compatibility tests for Phase 4
# Sources the new modular code and runs the test suite

cat("\n=== Running Backward Compatibility Tests ===\n\n")

# Source all module files from pkg/R/
cat("Loading modular WRS package...\n")

r_dir <- "/home/mando/coding/R-Projects/WRS/pkg/R"

# List of module files in order
module_files <- c(
  '00-utils-core.R',
  'location.R',
  'outliers.R',
  'bootstrap.R',
  'two-sample.R',
  'anova.R',
  'correlation.R',
  'ancova.R',
  'regression.R',
  'mcp.R',
  'covariance.R',
  'regression-advanced.R',
  'medians.R',
  'plotting.R',
  'effect-size.R',
  'power.R',
  'winsorize.R',
  'classification.R',
  'special.R',
  'zzz-internal.R'
)

# Source each module
for (f in module_files) {
  file_path <- file.path(r_dir, f)
  if (!file.exists(file_path)) {
    stop(sprintf("Module file not found: %s", file_path))
  }
  cat(sprintf("  Loading %s...\n", f))
  tryCatch({
    source(file_path)
  }, error = function(e) {
    cat(sprintf("  ✗ Error loading %s:\n", f))
    cat("   ", conditionMessage(e), "\n")
    stop("Failed to load module")
  })
}

cat("\n✓ All modules loaded successfully!\n\n")

# Load and run the test suite
cat("Loading test suite...\n")
source("/home/mando/coding/R-Projects/WRS/pkg/tests/test-backward-compat.R")

cat("Running tests...\n\n")
test_backward_compatibility()

cat("\n=== Test run complete ===\n")
