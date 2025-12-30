# Duplicate Functions Resolution Strategy

**Created**: 2025-12-30
**Status**: Phase 2, Week 5 - In Progress
**Context**: Resolving 1,153 duplicate function definitions across 20 modules

---

## Executive Summary

### The Problem
During Phase 1 module extraction, functions were distributed across modules but **special.R became a catch-all** that duplicated most functions from other modules. This created:

- **3,079 total function definitions** (should be ~1,926)
- **1,035 function names with duplicates**
- **1,153 net duplicates to remove**

### The Root Cause
**special.R contains 1,034 duplicate definitions** - it duplicates nearly every function from other modules because:
1. It was extracted last as the "catch-all" for remaining functions
2. Many specialized functions (oph.*, bin.*, etc.) properly belong there
3. But it also contains copies of functions that belong in domain-specific modules

### Duplicate Distribution by Module

| Module | Duplicate Definitions | Strategy |
|--------|----------------------|----------|
| special.R | 1,034 | **PRIMARY TARGET** - Remove most duplicates |
| ancova.R | 144 | Keep domain functions, remove from special.R |
| correlation.R | 108 | Keep domain functions, remove from special.R |
| two-sample.R | 103 | Keep domain functions, remove from special.R |
| mcp.R | 101 | Keep domain functions, remove from special.R |
| plotting.R | 94 | Keep domain functions, remove from special.R |
| regression.R | 91 | Keep domain functions, remove from special.R |
| location.R | 74 | Keep domain functions, remove from special.R |
| regression-advanced.R | 72 | Keep domain functions, remove from special.R |
| outliers.R | 63 | Keep domain functions, remove from special.R |
| anova.R | 56 | Keep domain functions, remove from special.R |
| 00-utils-core.R | 52 | Keep core utils, remove from special.R |
| covariance.R | 40 | Keep domain functions, remove from special.R |
| effect-size.R | 40 | Keep domain functions, remove from special.R |
| medians.R | 39 | Keep domain functions, remove from special.R |
| bootstrap.R | 30 | Keep bootstrap infrastructure, remove from special.R |
| classification.R | 27 | Keep classification functions, remove from special.R |
| power.R | 10 | Keep power analysis, remove from special.R |
| winsorize.R | 10 | Keep winsorization, remove from special.R |
| zzz-internal.R | 0 | No duplicates (expected for internal) |

---

## Resolution Strategy

### Phase 1: Automated Analysis & Preparation (COMPLETED ✅)
- [x] Scan all 20 modules and identify duplicates
- [x] Generate detailed duplicate report (saved to /tmp/duplicate-analysis.rds)
- [x] Create this strategy document

### Phase 2: Create Removal Script (NEXT)
Create an R script that will:

1. **Load duplicate analysis** from /tmp/duplicate-analysis.rds
2. **Apply resolution rules**:
   - **Rule 1**: If function appears in domain module + special.R → Keep in domain, remove from special.R
   - **Rule 2**: If function appears in multiple domain modules → Manual review needed
   - **Rule 3**: If function is truly specialized (oph.*, bin.*, run.*, selby*) → Keep in special.R only

3. **Generate removal lists**:
   - Functions to remove from special.R (expected: ~900+)
   - Functions requiring manual review (expected: ~50-100)

4. **Execute removals**:
   - Remove duplicate function definitions from appropriate files
   - Preserve only the canonical version of each function

### Phase 3: Manual Review
Handle edge cases:
- Functions appearing in 3+ modules
- Functions with different implementations (compare code)
- Ambiguous domain assignment

### Phase 4: Validation
- [x] Source all 20 modules (verify no errors)
- [ ] Run backward compatibility tests (must pass 100%)
- [ ] Verify unique function count matches original (1,926 unique)

---

## Resolution Rules (Detailed)

### Rule 1: Domain Module Takes Precedence
If a function exists in a domain-specific module AND special.R:
- **KEEP**: Version in domain module (ancova.R, regression.R, anova.R, etc.)
- **REMOVE**: Version in special.R

**Rationale**: Domain modules were carefully curated. special.R was a catch-all.

**Examples**:
- `ancova()` → Keep in ancova.R, remove from special.R
- `yuen()` → Keep in two-sample.R, remove from special.R
- `pbcor()` → Keep in correlation.R, remove from special.R

### Rule 2: Core Utilities Stay in 00-utils-core.R
If a function exists in 00-utils-core.R AND elsewhere:
- **KEEP**: Version in 00-utils-core.R
- **REMOVE**: Version from other modules

**Examples**:
- `elimna()`, `listm()`, `matl()` → Keep in 00-utils-core.R only

### Rule 3: True Specialized Functions Stay in special.R
Functions with these prefixes belong ONLY in special.R:
- `oph.*` - Ophthalmology-specific functions
- `bin.*` - Binomial-specific functions
- `run.*` - Run test functions
- `selby*` - Selection functions

**These should be removed from other modules IF they appear there.**

### Rule 4: Cross-Module Duplicates (3+ modules)
Requires manual review to determine:
1. Are implementations identical?
2. Which module is the most appropriate home?
3. Should it actually be in 00-utils-core.R as a shared utility?

**Process**:
1. Extract all versions
2. Compare implementations (diff)
3. Choose canonical location
4. Document decision

### Rule 5: Internal Helpers
Functions starting with `.` (dot) are internal and should be marked `@keywords internal` later in Phase 3. No action needed now.

---

## Implementation Plan

### Step 1: Create Automated Removal Script

```r
# Script: /tmp/remove-duplicates.R
# Purpose: Remove duplicate function definitions following resolution rules

# Load analysis
load("/tmp/duplicate-analysis.rds")

# Define resolution rules
resolve_duplicates <- function(all_functions, duplicates_df, duplicate_names) {

  # Rule hierarchy (priority order)
  module_priority <- c(
    "00-utils-core.R",      # Highest - core utilities
    "location.R",
    "outliers.R",
    "bootstrap.R",
    "two-sample.R",
    "anova.R",
    "correlation.R",
    "ancova.R",
    "regression.R",
    "regression-advanced.R",
    "covariance.R",
    "mcp.R",
    "medians.R",
    "plotting.R",
    "effect-size.R",
    "power.R",
    "winsorize.R",
    "classification.R",
    "zzz-internal.R",
    "special.R"              # Lowest - catch-all
  )

  # For each duplicate function name
  removals <- list()
  manual_review <- list()

  for (func_name in duplicate_names) {
    func_locations <- duplicates_df[duplicates_df$function_name == func_name, ]

    # Special case: Check if this is a specialized function
    is_specialized <- grepl("^(oph\\.|bin\\.|run\\.|selby)", func_name)

    if (is_specialized && "special.R" %in% func_locations$file) {
      # Keep in special.R, remove from everywhere else
      keep_in <- "special.R"
      remove_from <- setdiff(func_locations$file, "special.R")
    } else if (nrow(func_locations) == 2 && "special.R" %in% func_locations$file) {
      # Simple case: 2 copies, one in special.R
      # Keep the non-special.R version
      keep_in <- setdiff(func_locations$file, "special.R")[1]
      remove_from <- "special.R"
    } else if (nrow(func_locations) > 2 && !"special.R" %in% func_locations$file) {
      # Multiple domain modules - needs manual review
      manual_review[[func_name]] <- func_locations
      next
    } else {
      # Apply priority-based resolution
      func_locations$priority <- match(func_locations$file, module_priority)
      func_locations <- func_locations[order(func_locations$priority), ]

      keep_in <- func_locations$file[1]
      remove_from <- func_locations$file[-1]
    }

    # Record removals
    for (file in remove_from) {
      subset <- func_locations[func_locations$file == file, ]
      removals[[length(removals) + 1]] <- list(
        function_name = func_name,
        file = file,
        line_number = subset$line_number,
        keep_in = keep_in
      )
    }
  }

  list(
    removals = removals,
    manual_review = manual_review
  )
}

# Generate resolution plan
resolution <- resolve_duplicates(all_functions, duplicates_df, duplicate_names)

# Save resolution plan
save(resolution, file = "/tmp/duplicate-resolution-plan.rds")

# Print summary
cat("=== DUPLICATE RESOLUTION PLAN ===\n\n")
cat("Functions to auto-remove:", length(resolution$removals), "\n")
cat("Functions requiring manual review:", length(resolution$manual_review), "\n\n")

cat("=== REMOVALS BY FILE ===\n")
removal_by_file <- table(sapply(resolution$removals, function(x) x$file))
print(sort(removal_by_file, decreasing = TRUE))

cat("\n=== MANUAL REVIEW NEEDED ===\n")
if (length(resolution$manual_review) > 0) {
  for (fname in names(resolution$manual_review)) {
    cat(sprintf("%s:\n", fname))
    print(resolution$manual_review[[fname]])
    cat("\n")
  }
}
```

### Step 2: Execute Removals

After verifying the resolution plan:

```r
# Script: /tmp/execute-removals.R
# Purpose: Execute the removals generated by remove-duplicates.R

load("/tmp/duplicate-resolution-plan.rds")

# Function to remove a function definition from a file
remove_function_from_file <- function(file_path, function_name, start_line) {
  lines <- readLines(file_path, warn = FALSE)

  # Find the end of the function (matching braces)
  # This is complex - need to track brace matching
  # For now, use a simpler approach: find next function definition

  # Pattern for next function
  next_func_pattern <- "^[a-zA-Z][a-zA-Z0-9_\\.]*\\s*(<-|=)\\s*function"

  # Find all function starts
  func_starts <- grep(next_func_pattern, lines)

  # Find the next function after our target
  next_start <- func_starts[func_starts > start_line][1]

  if (is.na(next_start)) {
    # This is the last function - remove to end of file
    end_line <- length(lines)
  } else {
    # Remove up to (but not including) the next function
    end_line <- next_start - 1
  }

  # Remove the function
  new_lines <- lines[-(start_line:end_line)]

  return(new_lines)
}

# Execute removals
for (removal in resolution$removals) {
  file_path <- file.path("/home/mando/coding/R-Projects/WRS/pkg/R-new", removal$file)

  cat(sprintf("Removing %s from %s (line %d)... ",
              removal$function_name, removal$file, removal$line_number))

  # Load file
  lines <- readLines(file_path, warn = FALSE)

  # Remove function (this needs refinement)
  # For now, just mark the location

  cat("MARKED\n")
}

cat("\nRemoval execution complete!\n")
```

---

## Validation Checklist

After removing duplicates:

- [ ] Verify total function definitions reduced to ~1,926
- [ ] Source all 20 modules without errors
- [ ] Run backward compatibility tests (100% pass rate required)
- [ ] Check that no functions were accidentally removed entirely
- [ ] Verify special.R still contains oph.*, bin.*, run.*, selby* functions
- [ ] Update REFACTORING-PROGRESS.md with new statistics

---

## Expected Outcome

### Before (Current State)
- Total function definitions: 3,079
- Unique functions: 1,926
- Duplicates: 1,153

### After (Target State)
- Total function definitions: 1,926
- Unique functions: 1,926
- Duplicates: 0

### special.R Transformation
- Before: 1,797 functions (most are duplicates)
- After: ~700-800 functions (truly specialized functions only)
- Examples of what should remain:
  - oph.* functions (ophthalmology-specific)
  - bin.* functions (binomial-specific)
  - run.* functions (run tests)
  - selby* functions (selection methods)
  - Miscellaneous specialized methods not fitting other modules

---

## Risk Mitigation

### Backup Strategy
Before making changes:
```bash
# Create backup of pkg/R-new/
cp -r pkg/R-new pkg/R-new.BEFORE-DEDUP
```

### Incremental Approach
1. Start with special.R only (largest source of duplicates)
2. Test after special.R changes
3. Then address remaining cross-module duplicates
4. Test after each major change

### Rollback Plan
If tests fail:
```bash
# Restore from backup
rm -rf pkg/R-new
mv pkg/R-new.BEFORE-DEDUP pkg/R-new
```

---

## Notes & Observations

1. **special.R is massive**: 1,797 functions, 100,014 lines - this module needs the most cleanup

2. **Pattern**: Most duplicates follow the pattern:
   - Domain module has function (e.g., ancova.R has `ancova()`)
   - special.R also has same function (duplicate)
   - Resolution: Keep in domain, remove from special.R

3. **Some functions truly belong in special.R**:
   - oph.* (ophthalmology) - 50+ functions
   - bin.* (binomial) - 30+ functions
   - run.* (run tests) - 20+ functions
   - These should NOT appear in other modules

4. **Cross-module duplicates** (non-special.R):
   - Much rarer (~100 cases)
   - Often legitimate (e.g., plotting functions vs analysis functions)
   - Need careful manual review

5. **Testing is critical**:
   - Must maintain 100% backward compatibility
   - reference-outputs.rds provides baseline
   - All 23 tests must pass

---

## Timeline

- **Step 1** (1 hour): Create and test automated removal script
- **Step 2** (2 hours): Execute removals for special.R
- **Step 3** (1 hour): Handle cross-module duplicates
- **Step 4** (1 hour): Manual review and edge cases
- **Step 5** (1 hour): Testing and validation
- **Total**: ~6 hours estimated

---

## Success Criteria

✅ **Must Achieve:**
- Total function definitions == 1,926 (matching unique count)
- All 20 modules source without errors
- All 23 backward compatibility tests pass
- No functions accidentally removed entirely
- special.R reduced to ~700-800 functions

✅ **Quality Indicators:**
- Clear documentation of all removals
- Automated process (minimal manual intervention)
- Easy to verify and audit
- Rollback available if needed

---

*Last updated: 2025-12-30*
