# WRS Backward Compatibility Test Suite
# Compares current implementation against v0.45 reference outputs

test_backward_compatibility <- function() {
  # Load reference outputs from v0.45
  reference_file <- "/home/mando/coding/R-Projects/WRS/reference-outputs.rds"

  if (!file.exists(reference_file)) {
    stop("Reference outputs not found! Run create-reference-outputs.R first.")
  }

  reference <- readRDS(reference_file)

  cat("\n=== WRS Backward Compatibility Test ===\n")
  cat("Comparing against", length(reference), "reference test cases\n\n")

  # Set same seed as reference
  set.seed(12345)

  # Generate identical test data
  x <- rnorm(50)
  y <- rnorm(50) + 0.5
  x2 <- rnorm(30)
  y2 <- rnorm(30) + 0.3
  x3 <- rnorm(25)
  y3 <- rnorm(25) - 0.2

  # Matrix data
  m1 <- matrix(rnorm(100), ncol = 2)
  m2 <- matrix(rnorm(100) + 0.5, ncol = 2)

  # Run all current tests
  current_tests <- list(
    # Two-sample tests
    yuen = yuen(x, y, tr = 0.2),
    yuenbt = yuenbt(x, y, nboot = 100),
    pb2gen_mean = pb2gen(x, y, est = mean, nboot = 100),
    trimpb = trimpb(x, y, nboot = 100),

    # One-way ANOVA
    t1way = t1way(list(x, y, x2)),

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
  )

  # Compare with tolerance
  all_passed <- TRUE
  failed_tests <- character(0)
  tolerance <- 1e-10

  for (test_name in names(reference)) {
    if (!test_name %in% names(current_tests)) {
      warning(sprintf("  ✗ Test '%s' not found in current implementation!", test_name))
      all_passed <- FALSE
      failed_tests <- c(failed_tests, test_name)
      next
    }

    comparison <- all.equal(reference[[test_name]],
                           current_tests[[test_name]],
                           tolerance = tolerance)

    if (!isTRUE(comparison)) {
      cat(sprintf("  ✗ FAILED: %s\n", test_name))
      cat("    Difference:", comparison, "\n")
      all_passed <- FALSE
      failed_tests <- c(failed_tests, test_name)
    } else {
      cat(sprintf("  ✓ %s\n", test_name))
    }
  }

  # Summary
  cat("\n=== Test Summary ===\n")
  cat("Total tests:", length(reference), "\n")
  cat("Passed:", length(reference) - length(failed_tests), "\n")
  cat("Failed:", length(failed_tests), "\n")

  if (all_passed) {
    cat("\n✓✓✓ All backward compatibility tests PASSED! ✓✓✓\n")
    cat("Function outputs are identical to v0.45 (within tolerance", tolerance, ")\n\n")
  } else {
    cat("\n✗✗✗ Some tests FAILED ✗✗✗\n")
    cat("Failed tests:", paste(failed_tests, collapse = ", "), "\n\n")
    stop("Backward compatibility tests failed!")
  }

  invisible(list(
    passed = all_passed,
    n_total = length(reference),
    n_failed = length(failed_tests),
    failed_tests = failed_tests
  ))
}

# Wrapper for quick testing
quick_compat_check <- function() {
  tryCatch({
    test_backward_compatibility()
    TRUE
  }, error = function(e) {
    message("Error: ", e$message)
    FALSE
  })
}
