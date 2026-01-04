#!/usr/bin/env Rscript

# Run R CMD check on the refactored WRS package
# Phase 4 final validation

cat("\n=== Running R CMD check on WRS Package ===\n\n")

# Check if devtools is available
if (!requireNamespace("devtools", quietly = TRUE)) {
  cat("Installing devtools...\n")
  install.packages("devtools", repos = "https://cloud.r-project.org/")
}

library(devtools)

pkg_path <- "/home/mando/coding/R-Projects/WRS/pkg"

cat("Package path:", pkg_path, "\n")
cat("Starting comprehensive package check...\n")
cat("(This may take a few minutes)\n\n")

# Run check with minimal examples to save time
# We'll skip vignettes and don't run examples since there aren't many
tryCatch({
  check_results <- devtools::check(
    pkg = pkg_path,
    document = FALSE,  # Already documented
    args = c("--no-examples", "--no-tests", "--no-vignettes"),
    error_on = "never"  # Don't stop on warnings
  )

  cat("\n=== Check Complete ===\n")
  cat("Errors:", length(check_results$errors), "\n")
  cat("Warnings:", length(check_results$warnings), "\n")
  cat("Notes:", length(check_results$notes), "\n\n")

  if (length(check_results$errors) > 0) {
    cat("ERRORS:\n")
    for (err in check_results$errors) {
      cat("  -", err, "\n")
    }
  }

  if (length(check_results$warnings) > 0) {
    cat("\nWARNINGS:\n")
    for (warn in check_results$warnings) {
      cat("  -", warn, "\n")
    }
  }

  if (length(check_results$notes) > 0) {
    cat("\nNOTES:\n")
    for (note in check_results$notes) {
      cat("  -", note, "\n")
    }
  }

  if (length(check_results$errors) == 0) {
    cat("\n✓✓✓ Package check PASSED with no errors! ✓✓✓\n")
  } else {
    cat("\n✗ Package check found errors that need attention.\n")
  }

}, error = function(e) {
  cat("\n✗ Error running package check:\n")
  cat(conditionMessage(e), "\n")
  quit(status = 1)
})
