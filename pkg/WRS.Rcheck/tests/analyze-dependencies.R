# WRS Function Dependency Analysis
# Analyzes function calls to guide modularization

cat("Analyzing function dependencies in WRS package...\n\n")

# Read the source file
source_file <- "pkg/R/Rallfun-v45.R"
code <- readLines(source_file)

# Extract all function names
func_pattern <- "^([a-zA-Z][a-zA-Z0-9_\\.]*)\\s*<-\\s*function"
func_lines <- grep(func_pattern, code, value = FALSE)
func_names <- gsub("\\s*<-.*", "", grep(func_pattern, code, value = TRUE))
func_names <- trimws(func_names)

cat("Total functions found:", length(func_names), "\n\n")

# Create lookup for function start/end lines (approximate)
func_positions <- data.frame(
  name = func_names,
  start_line = func_lines,
  stringsAsFactors = FALSE
)

# For each function, find its end line (next function start or EOF)
func_positions$end_line <- c(func_positions$start_line[-1] - 1, length(code))

cat("Building dependency map (this may take a minute)...\n")

# For each function, find what other WRS functions it calls
dependencies <- list()
call_counts <- setNames(rep(0, length(func_names)), func_names)

for (i in 1:nrow(func_positions)) {
  func_name <- func_positions$name[i]
  start <- func_positions$start_line[i]
  end <- func_positions$end_line[i]

  # Get function body
  func_body <- code[start:end]
  func_text <- paste(func_body, collapse = "\n")

  # Find calls to other WRS functions
  called_funcs <- c()
  for (other_func in func_names) {
    if (other_func == func_name) next  # Skip self-calls

    # Look for function calls: funcname(
    pattern <- paste0("\\b", other_func, "\\s*\\(")
    if (grepl(pattern, func_text)) {
      called_funcs <- c(called_funcs, other_func)
      call_counts[other_func] <- call_counts[other_func] + 1
    }
  }

  dependencies[[func_name]] <- called_funcs

  # Progress indicator
  if (i %% 200 == 0) cat("  Processed", i, "/", nrow(func_positions), "functions\n")
}

cat("✓ Dependency analysis complete\n\n")

# Identify core utilities (most-called functions)
core_utils <- sort(call_counts, decreasing = TRUE)
top_50 <- head(core_utils[core_utils > 0], 50)

cat("=== TOP 50 CORE UTILITIES (Most Called Functions) ===\n")
cat("These are candidates for 00-utils-core.R\n\n")
for (i in 1:min(50, length(top_50))) {
  cat(sprintf("%2d. %-20s called by %3d functions\n",
              i, names(top_50)[i], top_50[i]))
}

# Identify isolated functions (call nothing)
isolated <- names(dependencies)[sapply(dependencies, length) == 0]
cat("\n=== ISOLATED FUNCTIONS ===\n")
cat("Functions that don't call other WRS functions:", length(isolated), "\n")
if (length(isolated) <= 20) {
  cat("Examples:", paste(head(isolated, 20), collapse = ", "), "\n")
}

# Identify highly connected functions (call many others)
dep_counts <- sapply(dependencies, length)
highly_connected <- sort(dep_counts, decreasing = TRUE)
top_20_connected <- head(highly_connected[highly_connected > 0], 20)

cat("\n=== TOP 20 HIGHLY CONNECTED FUNCTIONS ===\n")
cat("Functions that call many other WRS functions\n\n")
for (i in 1:length(top_20_connected)) {
  cat(sprintf("%2d. %-20s calls %3d other functions\n",
              i, names(top_20_connected)[i], top_20_connected)[i])
}

# Save results
results <- list(
  total_functions = length(func_names),
  all_functions = func_names,
  dependencies = dependencies,
  call_counts = call_counts,
  core_utilities = top_50,
  isolated = isolated,
  highly_connected = top_20_connected
)

saveRDS(results, "dependency-analysis.rds")
cat("\n✓ Results saved to dependency-analysis.rds\n")

# Identify likely core utility functions by name patterns
cat("\n=== LIKELY CORE UTILITIES BY NAME PATTERN ===\n")
core_patterns <- c(
  "elimna", "listv2mat", "winmean", "winvar", "winval", "trimse",
  "tmean", "tvar", "hd", "mest", "onestep", "idealf",
  "near", "selby", "sel", "mat2grp", "grp2mat"
)

found_core <- func_names[func_names %in% core_patterns]
cat("Pattern-based core functions found:", length(found_core), "\n")
cat(paste(found_core, collapse = ", "), "\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Use this data to guide module extraction in Phase 1\n")
