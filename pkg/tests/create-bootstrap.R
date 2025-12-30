# Extract bootstrap.R - Bootstrap and Resampling Infrastructure
# Generic bootstrap, permutation, and resampling methods

source("pkg/tests/extract-functions.R")

# Bootstrap and resampling infrastructure
# Note: Many specific bootstrap functions (yuenbt, etc.) go in their
# respective modules (two-sample.R, anova.R, etc.)
bootstrap_funcs <- c(
  # Generic bootstrap functions
  "bootse",        # Bootstrap standard error
  "bootdse",       # Bootstrap dependent SE
  "bootdpci",      # Bootstrap dependent PB CI
  "bootdep",       # Bootstrap dependent
  "bootdepsub",    # Subroutine
  "bootcov",       # Bootstrap covariance
  "boot.TM",       # Bootstrap trimmed mean

  # Regression bootstrap (generic)
  "regboot",       # Regression bootstrap
  "regbootg",      # Regression bootstrap group
  "regbootMC",     # MC version

  # ANCOVA bootstrap (generic)
  "ancboot",       # ANCOVA bootstrap
  "ancbootg",      # ANCOVA bootstrap group

  # BCA (Bias-Corrected Accelerated) methods
  "bca.mean",      # BCA mean (in location.R too?)
  "wmw.bca",       # Wilcoxon-Mann-Whitney BCA
  "wmw.bcapv",     # BCA p-value
  "wmw.bcapv.v2",  # BCA p-value v2
  "wmw.bcav2",     # BCA variant 2
  "corblp.bca.C",  # Correlation BCA C
  "pbcan",         # Percentile bootstrap canonical

  # Permutation tests
  "permg",         # Permutation group
  "permg.t",       # Permutation group t-test
  "perm.rho",      # Permutation correlation
  "perm.rho.sub",  # Subroutine
  "bwdepth.perm",  # Bandwidth depth permutation

  # MC subroutines (generic bootstrap helpers)
  "outproMC.sub",  # Outlier projection MC sub (duplicated from outliers?)
  "outproMC.sub2", # Sub 2
  "linWMWMC.sub",  # Linear WMW MC sub
  "linWMWMC.sub2", # Sub 2
  "ridgeGnullMC.sub", # Ridge G null MC sub

  # Additional bootstrap/resampling utilities
  "btsqrk"         # Bootstrap SQRK
)

cat("Creating bootstrap.R with", length(bootstrap_funcs), "bootstrap/resampling functions\n\n")

# Create header
header <- c(
  "# WRS Package - Bootstrap and Resampling Infrastructure",
  "# Generic bootstrap, permutation, and resampling methods",
  "#",
  "# This module contains:",
  "#   - Generic bootstrap functions: bootse, bootdep, bootcov",
  "#   - BCA (Bias-Corrected Accelerated) methods",
  "#   - Permutation test infrastructure",
  "#   - Bootstrap helpers and subroutines",
  "#",
  "# Note: Specific bootstrap analyses (yuenbt, etc.) are in their",
  "# respective modules (two-sample.R, anova.R, regression.R, etc.)",
  ""
)

# Extract the functions
result <- extract_functions(
  func_names = bootstrap_funcs,
  output_file = "pkg/R-new/bootstrap.R",
  header = header
)

cat("\n=== EXTRACTION SUMMARY ===\n")
cat("Functions extracted:", result$found, "/", length(bootstrap_funcs), "\n")
if (length(result$not_found) > 0) {
  cat("\nNot found (", length(result$not_found), "):\n")
  cat(paste(result$not_found, collapse = ", "), "\n")
}

# Test sourcing
cat("\n=== TESTING EXTRACTION ===\n")
cat("Loading all Week 1 modules...\n")

tryCatch({
  source("pkg/R-new/00-utils-core.R")
  source("pkg/R-new/location.R")
  source("pkg/R-new/outliers.R")
  source("pkg/R-new/bootstrap.R")
  cat("✓ All Week 1 modules sourced successfully!\n")
  cat("✓ Bootstrap module loaded without errors\n")
},
error = function(e) {
  cat("✗ ERROR:\n")
  cat("  ", conditionMessage(e), "\n")
})

cat("\n=== WEEK 1 FOUNDATION MODULES COMPLETE ===\n")
cat("Extracted modules:\n")
cat("  1. 00-utils-core.R - 50 core utilities\n")
cat("  2. location.R - 73 location estimators\n")
cat("  3. outliers.R - 64 outlier detection methods\n")
cat("  4. bootstrap.R - ", result$found, " bootstrap/resampling functions\n")
cat("\n=== NEXT STEPS ===\n")
cat("1. Run backward compatibility tests\n")
cat("2. Proceed to Week 2 modules (two-sample, anova, correlation)\n")
