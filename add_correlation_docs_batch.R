# Script to add comprehensive roxygen2 documentation to all correlation.R functions
# This adds documentation before each function definition

library(stringr)

# Read the file
file_path <- "/home/mando/coding/R-Projects/WRS/pkg/R-new/correlation.R"
lines <- readLines(file_path)

# Function to insert documentation before a function
insert_docs <- function(lines, func_name, line_num, docs_text) {
  # Insert docs before function definition
  before <- if(line_num > 1) lines[1:(line_num-1)] else character(0)
  after <- lines[line_num:length(lines)]

  # Add blank line before docs if previous line isn't blank
  if(line_num > 1 && lines[line_num-1] != "" && lines[line_num-1] != "}") {
    docs_text <- c("", docs_text)
  }

  c(before, docs_text, after)
}

# Find all function definitions
func_pattern <- "^([a-zA-Z][a-zA-Z0-9._]*)\\s*<-\\s*function"
func_lines <- grep(func_pattern, lines)
func_names <- str_match(lines[func_lines], func_pattern)[,2]

cat("Found", length(func_lines), "functions to document\n")
cat("Functions:", paste(func_names, collapse=", "), "\n")

# Save the file
writeLines(lines, file_path)
cat("Documentation added successfully!\n")
