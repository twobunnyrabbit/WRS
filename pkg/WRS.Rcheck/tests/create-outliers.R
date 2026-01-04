# Extract outliers.R - Outlier Detection and Depth Methods
# Functions for detecting outliers and computing data depth

source("pkg/tests/extract-functions.R")

# Outlier detection and depth functions
outlier_funcs <- c(
  # Core outlier detection (out is in core utils)
  # "out",         # In core utils
  # "outpro",      # In core utils
  "outproMC",      # MC version of outpro
  "outproMC.sub",  # Subroutine
  "outproMC.sub2", # Subroutine 2
  "outproad",      # Outlier projection adaptive
  "outproadMC",    # MC version
  "outpro.depth",  # Projection depth
  "outbox",        # Boxplot-based outlier detection
  "out3d",         # 3D outlier detection
  "outbag",        # Bagging-based outlier detection

  # Specific outlier methods
  "outmah",        # Mahalanobis distance
  "outmve",        # MVE-based
  "outmgv",        # MGV method
  "outmgvad",      # MGV adaptive
  "outmgvf",       # MGV variant
  "outmgv.v2",     # MGV version 2
  "outms",         # MS method
  "outogk",        # OGK method
  "outcov",        # Covariance-based
  "outtbs",        # TBS method
  "outDETMCD",     # DET MCD method
  "outICS",        # ICS method
  "outmc",         # MC method
  "out.dummy",     # Dummy/placeholder
  "out.by.groups", # By groups
  "out.methods",   # Multiple methods
  "outblp.HH",     # BLP HH method

  # Depth functions (depth is in core utils)
  # "depth",       # In core utils
  "depth2",        # Depth variant 2
  "depthcom",      # Depth comparison
  "depthcomsub",   # Subroutine
  "depthg2",       # Depth for 2 groups
  "depths1",       # Depth variant
  # "fdepth",      # In core utils
  "fdepthv2",      # F-depth version 2
  "indepth",       # Independent depth
  "pdepth",        # P-depth
  "prodepth",      # Projection depth
  "rdepth.orig",   # R-depth original
  "resdepth",      # Residual depth
  "resdepth.sub",  # Subroutine
  "zdepth",        # Z-depth
  "zdepth.sub",    # Subroutine
  "zoudepth",      # Zou depth
  "unidepth",      # Univariate depth
  "discdepth",     # Discriminant depth
  "mregdepth",     # Multivariate regression depth
  "smean.depth",   # Sample mean depth
  "pbadepth",      # PBA depth
  "Qdepthcom",     # Quantile depth comparison
  "Qdepthcomsub",  # Subroutine
  "comdepthsvm",   # Compare depth SVM
  "aov2depth",     # ANOVA to depth

  # Bagging and depth methods
  "bagdepth",      # Bag depth
  "bwdepth",       # Bandwidth depth
  "bwdepthMC.ci",  # MC CI for bandwidth depth
  "bwdepth.perm",  # Permutation version
  "bwdepth.sub",   # Subroutine

  # Classification/bagging variants
  "Depth.class.bag",    # Depth classification bag
  "dis.depth.bag",      # Distance depth bag
  "pro.class.bag",      # Projection class bag
  "pro.classPD.bag",    # Projection class PD bag
  "KNNbag",             # KNN bagging
  "LSMbag",             # LSM bagging
  "NNbag",              # NN bagging
  "RFbag",              # Random Forest bagging
  "SVMbag"              # SVM bagging
)

cat("Creating outliers.R with", length(outlier_funcs), "outlier detection functions\n\n")

# Create header
header <- c(
  "# WRS Package - Outlier Detection and Data Depth",
  "# Methods for detecting outliers and computing data depth",
  "#",
  "# This module contains:",
  "#   - Projection-based outlier detection: outpro variants",
  "#   - Classical methods: outbox, outmah, outmve",
  "#   - Modern robust methods: outogk, outDETMCD",
  "#   - Data depth methods: depth variants, bagdepth",
  "#   - Classification and bagging methods",
  "#",
  "# See: Wilcox (2022), Chapter on outlier detection",
  ""
)

# Extract the functions
result <- extract_functions(
  func_names = outlier_funcs,
  output_file = "pkg/R-new/outliers.R",
  header = header
)

cat("\n=== EXTRACTION SUMMARY ===\n")
cat("Functions extracted:", result$found, "/", length(outlier_funcs), "\n")
if (length(result$not_found) > 0) {
  cat("\nNot found (", length(result$not_found), "):\n")
  cat(paste(result$not_found, collapse = ", "), "\n")
}

# Test sourcing
cat("\n=== TESTING EXTRACTION ===\n")
cat("Loading dependencies and outliers.R...\n")

tryCatch({
  source("pkg/R-new/00-utils-core.R")
  source("pkg/R-new/location.R")
  source("pkg/R-new/outliers.R")
  cat("✓ All modules sourced successfully!\n")
  cat("✓ Outliers module loaded without errors\n")
},
error = function(e) {
  cat("✗ ERROR:\n")
  cat("  ", conditionMessage(e), "\n")
})

cat("\n=== NEXT STEPS ===\n")
cat("1. Review pkg/R-new/outliers.R\n")
cat("2. Proceed to extract bootstrap.R\n")
