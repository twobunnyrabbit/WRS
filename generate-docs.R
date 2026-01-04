#!/usr/bin/env Rscript

# Generate package documentation with roxygen2
# Phase 4 of WRS refactoring

cat("=== Phase 4: Generating Package Documentation ===\n\n")

# Check if roxygen2 is installed
if (!requireNamespace("roxygen2", quietly = TRUE)) {
  cat("Installing roxygen2...\n")
  install.packages("roxygen2", repos = "https://cloud.r-project.org/")
}

library(roxygen2)

# Set working directory to package root
pkg_path <- "/home/mando/coding/R-Projects/WRS/pkg"
setwd(pkg_path)

cat("Package path:", pkg_path, "\n")
cat("Generating documentation with roxygen2::roxygenize()...\n\n")

# Generate documentation
tryCatch({
  roxygen2::roxygenize(package.dir = pkg_path,
                       roclets = c("namespace", "rd"),
                       load_code = FALSE)
  cat("\n✓ Documentation generated successfully!\n")
}, error = function(e) {
  cat("\n✗ Error generating documentation:\n")
  cat(conditionMessage(e), "\n")
  quit(status = 1)
})

# Count generated .Rd files
rd_files <- list.files(file.path(pkg_path, "man"), pattern = "\\.Rd$")
cat("\nGenerated", length(rd_files), ".Rd documentation files\n")

# Show first 10 files as examples
if (length(rd_files) > 0) {
  cat("\nFirst 10 documentation files:\n")
  print(head(rd_files, 10))
}

cat("\n=== Documentation generation complete! ===\n")
