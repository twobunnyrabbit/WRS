# Extract 00-utils-core.R - Foundation Utilities
# Based on dependency analysis: top 50 most-called functions

source("pkg/tests/extract-functions.R")

# Top 50 core utilities from dependency analysis
# These are called by the most other functions
core_utils <- c(
  # Top tier (called by 100+ functions)
  "elimna",      # 928 calls - CRITICAL
  "listm",       # 327 calls
  "matl",        # 219 calls

  # High priority (50+ calls)
  "near",        # 62 calls
  "hd",          # 61 calls
  "winvar",      # 53 calls
  "lplot",       # 51 calls
  "qest",        # 50 calls

  # Important utilities (20-50 calls)
  "pool.a.list", # 46 calls
  "con2way",     # 43 calls
  "rmul",        # 34 calls
  "regYhat",     # 34 calls
  "pdis",        # 30 calls
  "yuen",        # 29 calls
  "akerd",       # 29 calls
  "chi.int2",    # 27 calls
  "outpro",      # 27 calls
  "idealf",      # 26 calls
  "runmean2g",   # 26 calls
  "smmcrit",     # 25 calls
  "depth",       # 23 calls
  "lplot.pred",  # 23 calls
  "trimci",      # 22 calls
  "standm",      # 22 calls
  "con3way",     # 22 calls
  "binmat",      # 22 calls
  "near3d",      # 21 calls

  # Common utilities (15-20 calls)
  "kron",        # 19 calls
  "fac2list",    # 19 calls
  "smmcrit01",   # 18 calls
  "chi.int",     # 18 calls
  "psi.bt",      # 18 calls
  "rdplot",      # 17 calls
  "binom.conf",  # 17 calls
  "pdisMC",      # 16 calls
  "wincor",      # 15 calls
  "winall",      # 15 calls
  "erho.bt",     # 15 calls
  "yuend",       # 15 calls
  "fdepth",      # 15 calls
  "ols",         # 15 calls
  "rplot",       # 15 calls

  # Additional common utilities (10-15 calls)
  "out",         # 14 calls
  "bmp",         # 14 calls
  "trimse",      # 13 calls
  "ancova",      # 13 calls
  "covmtrim",    # 13 calls
  "olshc4",      # 13 calls
  "con.all.pairs", # 13 calls
  "rmmcppb"      # 13 calls
)

cat("Creating 00-utils-core.R with", length(core_utils), "core utility functions\n\n")

# Create header
header <- c(
  "# WRS Package - Core Utilities",
  "# Foundation functions used throughout the package",
  "#",
  "# This file contains the most commonly called utility functions,",
  "# extracted from Rallfun-v45.R based on dependency analysis.",
  "#",
  "# Top functions by call count:",
  "#   - elimna: Called by 928 functions (47% of package!)",
  "#   - listm: Called by 327 functions",
  "#   - matl: Called by 219 functions",
  "#",
  "# These functions MUST be loaded before other modules.",
  ""
)

# Extract the functions
result <- extract_functions(
  func_names = core_utils,
  output_file = "pkg/R-new/00-utils-core.R",
  header = header
)

cat("\n=== EXTRACTION SUMMARY ===\n")
cat("Functions extracted:", result$found, "/", length(core_utils), "\n")
cat("Output file:", result$output_file, "\n")

# Test that the file can be sourced
cat("\n=== TESTING EXTRACTION ===\n")
cat("Attempting to source 00-utils-core.R...\n")

tryCatch({
  source("pkg/R-new/00-utils-core.R")
  cat("✓ File sourced successfully!\n")
  cat("✓ All", result$found, "functions loaded without errors\n")
},
error = function(e) {
  cat("✗ ERROR sourcing file:\n")
  cat("  ", conditionMessage(e), "\n")
  cat("\nThis may indicate missing dependencies that need to be extracted first.\n")
})

cat("\n=== NEXT STEPS ===\n")
cat("1. Review pkg/R-new/00-utils-core.R\n")
cat("2. Test individual functions\n")
cat("3. Proceed to extract location.R module\n")
