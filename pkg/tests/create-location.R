# Extract location.R - Robust Location Estimators
# Functions for estimating central tendency using robust methods

source("pkg/tests/extract-functions.R")

# Core location estimators and related functions
# Based on robust statistics literature and WRS documentation
location_funcs <- c(
  # Primary robust location estimators
  "mest",          # M-estimator
  "mestci",        # M-estimator with CI
  "mestse",        # M-estimator standard error
  "mestseb",       # M-estimator SE bootstrap
  "mom",           # Modified one-step M-estimator
  "momci",         # mom with CI
  "onestep",       # One-step M-estimator
  "tmean",         # Trimmed mean (if not in core)
  "tvar",          # Trimmed variance (if not in core)

  # Harrell-Davis estimator and variants
  # Note: hd is in core utils, but hd-related functions go here
  "hdci",          # HD confidence interval
  "hdpb",          # HD percentile bootstrap
  "hdseb",         # HD standard error bootstrap
  "hdno",          # HD for different purposes
  "hdmq",          # HD mid-quantile
  "hdep",          # HD related
  "IQRhd",         # IQR using HD
  "medhd2g",       # Median HD for 2 groups
  "Dqcihd",        # Dependent quantile CI HD
  "Dqcomhd",       # Dependent quantile comparison HD
  "cbmhd",         # Compare bootstrap medians HD
  "cbmhdMC",       # MC version
  "Dcbmhd",        # Dependent version

  # Trimmed means and variants
  "trimci",        # Already in core, but related funcs here
  # "trimse",      # Already in core utils
  "ghtrim",        # Generalized Huber trim
  "ftrim",         # F trim
  "btrim",         # Bootstrap trim
  "bbtrim",        # Double bootstrap trim
  "bwtrim",        # Bootstrap Winsorized trim
  "bbwtrim",       # Double bootstrap Winsorized trim
  "bwtrimbt",      # Bootstrap Winsorized trim bootstrap-t
  "bbwtrimbt",     # Double bootstrap version
  "bwwtrim",       # Bootstrap double Winsorized trim
  "bwwtrimbt",     # Bootstrap version
  "dtrimpb",       # Dependent trim percentile bootstrap
  "dtrimQS",       # Dependent trim quantile shift
  "dlintrim",      # Dependent linear trim
  "ghtrim",        # Generalized Huber trim
  "bbtrimQS",      # Double bootstrap trim QS

  # Bivariate location
  "biloc",         # Bivariate location
  "dep.loc.summary", # Dependent location summary

  # Location difference/comparison functions
  "loc2dif",       # Location difference for 2 groups
  "loc2difpb",     # with percentile bootstrap
  "loc2gmulpb",    # 2-group multiple PB
  "loc2plot",      # Plot for 2-group location
  "mdifloc",       # Multiple difference location
  "M2m.loc",       # Matrix to matrix location
  "mul.loc2g",     # Multiple location 2 groups
  "Dmul.loc2g",    # Dependent multiple location
  "loc.dif.summary", # Location difference summary

  # Specialized location functions
  "L1median",      # L1 (spatial) median
  "bmean",         # Bootstrap mean
  "dmean",         # Dependent mean
  "mmean",         # Matrix mean
  "ghmean",        # Generalized Huber mean
  "mgvmean",       # Another variant
  "harmonic.mean", # Harmonic mean
  "lognormal.mean", # Log-normal mean
  "lognormal.mom",  # Log-normal MOM
  "bca.mean",      # BCa mean
  "bptdmean",      # Bootstrap t-dependent mean

  # Other location-related
  "center.m",      # Center matrix
  "funloc",        # Function location
  "funlocpb",      # with PB
  "locpre",        # Location prediction
  "locpres1",      # Location prediction variant
  "locreg",        # Location regression
  "locvar",        # Location variance
  "locvarsm",      # Location variance smooth
  "locCV",         # Location cross-validation
  "llocv2",        # Local location variant
  "covloc",        # Covariance location
  "ogk.center",    # OGK center
  "ghdist",        # Generalized Hausdorff distance
  "meancr.cord.oph" # Mean related
)

cat("Creating location.R with", length(location_funcs), "location estimator functions\n\n")

# Create header
header <- c(
  "# WRS Package - Robust Location Estimators",
  "# Methods for estimating central tendency and location",
  "#",
  "# This module contains:",
  "#   - M-estimators: mest, mom, onestep",
  "#   - Harrell-Davis quantile estimators: hd variants",
  "#   - Trimmed means and variants",
  "#   - Bivariate and multivariate location estimators",
  "#   - Location comparison and difference functions",
  "#",
  "# See: Wilcox (2022), Introduction to Robust Estimation",
  ""
)

# Extract the functions
result <- extract_functions(
  func_names = location_funcs,
  output_file = "pkg/R-new/location.R",
  header = header
)

cat("\n=== EXTRACTION SUMMARY ===\n")
cat("Functions extracted:", result$found, "/", length(location_funcs), "\n")
if (length(result$not_found) > 0) {
  cat("\nNot found (", length(result$not_found), "):\n")
  cat(paste(result$not_found, collapse = ", "), "\n")
}

# Test sourcing (requires core utils first)
cat("\n=== TESTING EXTRACTION ===\n")
cat("Loading core utils first, then location.R...\n")

tryCatch({
  source("pkg/R-new/00-utils-core.R")
  source("pkg/R-new/location.R")
  cat("✓ Both files sourced successfully!\n")
  cat("✓ Location module loaded without errors\n")
},
error = function(e) {
  cat("✗ ERROR:\n")
  cat("  ", conditionMessage(e), "\n")
})

cat("\n=== NEXT STEPS ===\n")
cat("1. Review pkg/R-new/location.R\n")
cat("2. Proceed to extract outliers.R\n")
