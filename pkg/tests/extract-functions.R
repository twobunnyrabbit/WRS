# Function Extraction Helper
# Extracts specific functions from Rallfun-v45.R by name

extract_functions <- function(func_names, output_file, header = NULL) {
  source_file <- "pkg/R/Rallfun-v45.R"
  code <- readLines(source_file)

  cat("Extracting", length(func_names), "functions to", output_file, "\n")

  # Find all function definitions
  func_pattern <- "^([a-zA-Z][a-zA-Z0-9_\\.]*)\\s*<-\\s*function"
  all_func_lines <- grep(func_pattern, code, value = FALSE)
  all_func_names <- gsub("\\s*<-.*", "", grep(func_pattern, code, value = TRUE))
  all_func_names <- trimws(all_func_names)

  # Create function position map
  func_positions <- data.frame(
    name = all_func_names,
    start_line = all_func_lines,
    stringsAsFactors = FALSE
  )
  func_positions$end_line <- c(func_positions$start_line[-1] - 1, length(code))

  # Extract requested functions
  extracted_code <- c()

  if (!is.null(header)) {
    extracted_code <- c(extracted_code, header, "")
  }

  found_count <- 0
  not_found <- c()

  for (fname in func_names) {
    pos <- func_positions[func_positions$name == fname, ]

    if (nrow(pos) == 0) {
      not_found <- c(not_found, fname)
      next
    }

    if (nrow(pos) > 1) {
      cat("  WARNING: Multiple definitions of", fname, "- using first\n")
      pos <- pos[1, ]
    }

    # Extract function code
    func_code <- code[pos$start_line:pos$end_line]

    # Add separator comment
    extracted_code <- c(
      extracted_code,
      paste0("# ", fname),
      func_code,
      ""
    )

    found_count <- found_count + 1
  }

  # Write output
  writeLines(extracted_code, output_file)

  cat("✓ Extracted", found_count, "functions\n")
  if (length(not_found) > 0) {
    cat("✗ Not found:", paste(not_found, collapse = ", "), "\n")
  }

  return(list(
    found = found_count,
    not_found = not_found,
    output_file = output_file
  ))
}

# Usage example:
# core_funcs <- c("elimna", "listm", "matl", "near", "hd", "winvar")
# extract_functions(core_funcs, "pkg/R-new/00-utils-core.R",
#                   header = "# WRS Core Utilities\n# Foundation functions used throughout the package")
