# Correlation.R Documentation Status

## Overview
- **Total Functions**: 83
- **File Size**: 4,669 lines
- **Documented So Far**: ~17 core functions (20%)
- **Remaining**: ~66 functions (80%)

## Completed Documentation (17 functions)

### 1. Common Parameters Block
Added comprehensive `@inheritParams common-params` documentation block covering:
- x, y, tr, beta, corfun, alpha, nboot, SEED
- plotit, xlab, ylab, outfun, xout, cop, MC, pr, STAND

### 2. Core Correlation Functions (Documented)
1. **wincor.sub** - Winsorized correlation subroutine (@keywords internal)
2. **tau** - Kendall's tau with CI (@export) ✓
3. **tauall** - Kendall's tau for all pairs (@export) ✓
4. **pbcor** - Percentage bend correlation (@export) ✓
5. **corb** - Bootstrap CI for correlation (@export) ✓
6. **runcor** - Running correlation (@export) ✓
7. **pcorb** - Bootstrap CI for Pearson correlation (@export) ✓
8. **correg.sub** - Correlation regression subroutine (@keywords internal)
9. **correg** - Correlation-based regression (@export) ✓
10. **tauloc** - Tau measure of location (@export) ✓
11. **tauvar** - Tau measure of scale (@export) ✓
12. **taulc** - Tau location and scale (@export) ✓
13. **ecor** - Explanatory correlation with outlier detection (@export) ✓
14. **ocor** - Correlation with multivariate outlier detection (@export) ✓
15. **cori** - Correlation at specified value of third variable (@export) ✓
16. **pcor** - Pearson correlation (@export) ✓
17. **mscor** - Multiple skipped correlation (@export) ✓

## Remaining Functions Needing Documentation (66)

### High Priority User-Facing Functions (@export)
1. **rhom** - Robust heteroscedastic regression correlation
2. **tbscor** - Two-sample bootstrap correlation comparison
3. **scor** - Skipped correlation (CRITICAL)
4. **scorci** - Skipped correlation with CI (CRITICAL)
5. **scorciMC** - Skipped correlation CI (parallel)
6. **ogkcor** - OGK correlation
7. **pcorhc4** - Pearson correlation with HC4 standard errors
8. **smcorcom** - Smooth correlation comparison
9. **pcorbv4** - Pearson correlation bootstrap variant 4
10. **corbMC** - Bootstrap correlation (parallel)
11. **qcorp1** - Quantile correlation projection method 1
12. **qcor** - Quantile correlation
13. **tauci** - Kendall's tau with bootstrap CI
14. **tscor** - Theil-Sen correlation
15. **tscorci** - Theil-Sen correlation with CI
16. **wincorci** - Winsorized correlation with CI
17. **corCOMmcp** - Multiple comparison procedure for correlations
18. **regcor** - Regression-based correlation
19. **mscorpb** - Multiple skipped correlation with bootstrap
20. **mscorpbMC** - Multiple skipped correlation bootstrap (parallel)
21. **scorv2** - Skipped correlation version 2
22. **scorreg** - Skipped correlation regression
23. **mscorci** - Multiple skipped correlation with CIs
24. **scorregciH** - Skipped correlation regression CI (Hochberg)
25. **mscorciH** - Multiple skipped correlation CI (Hochberg)
26. **scorregci** - Skipped correlation regression with CI
27. **scorall** - Skipped correlation for all pairs
28. **corxy** - Correlation between two sets of variables
29. **rhohc4bt** - Pearson correlation HC4 bootstrap test
30. **corCOM.DVvsIV** - Multiple comparison: DV vs IVs
31. **corREGorder** - Correlation regression order test
32. **bicor** - Biweight midcorrelation
33. **bicorM** - Biweight midcorrelation matrix
34. **corregci** - Correlation regression with CI
35. **mcd.cor** - MCD correlation
36. **corCOM.PMDPCD** - Correlation comparison PMD vs PCD
37. **corREG.best** - Best subset correlation regression
38. **PcorREG.best.DO** - Pearson correlation regression best subset (dominance)
39. **corREG.DO** - Correlation regression dominance analysis
40. **cor.skip.com** - Skipped correlation comparison
41. **corskip.comPV** - Skipped correlation comparison p-values
42. **rmdif.scores** - RM difference scores
43. **smeancr.cord.oph** - Smooth mean correlation omnibus hypothesis
44. **smeancr.cord** - Smooth mean correlation
45. **part.cor** - Partial correlation
46. **corblp.EP** - Correlation based on linear predictor (EP)
47. **corblp.ci** - Correlation linear predictor with CI
48. **corblppb** - Correlation linear predictor percentile bootstrap
49. **qcorp1.ci** - Quantile correlation P1 with CI
50. **qcor.ci** - Quantile correlation with CI
51. **qcor.ep.ci** - Quantile correlation EP with CI
52. **qcor.R** - Quantile correlation range
53. **qcor.EP** - Quantile correlation expected value
54. **qcor.ep** - Quantile correlation EP method

### Internal/Helper Functions (@keywords internal)
55. **pcorhc4sub** - Pearson correlation HC4 subroutine
56. **scorsubMC** - Skipped correlation MC subroutine
57. **mscorci.cr** - Multiple skipped correlation critical values
58. **mscorci.cr.sub** - Multiple skipped correlation critical values sub
59. **scorci.sub** - Skipped correlation CI subroutine
60. **scorreg.sub** - Skipped correlation regression subroutine (DUPLICATE)
61. **scorreg.cr** - Skipped correlation regression critical values
62. **scorreg.cr.sub** - Skipped correlation regression critical values sub
63. **corREGorder.crit** - Correlation regression order critical values
64. **corCOM.DVvsIV.crit** - Correlation DV vs IV critical values
65. **corCOM.PMDPCD.sub** - Correlation PMD vs PCD subroutine
66. **Note**: There is a duplicate scorreg.sub function (line 1830)

## Documentation Pattern Used

All documented functions follow this comprehensive pattern:

```r
#' Function Title
#'
#' Brief description (one sentence).
#'
#' @inheritParams common-params
#' @param specific_param Description of specific parameter.
#'
#' @return Description of return value with \\item{} for list components.
#'
#' @details
#' Extended description of how the function works, what methods it uses,
#' special considerations, etc.
#'
#' @seealso \\code{\\link{related_function1}}, \\code{\\link{related_function2}}
#'
#' @references
#' Author, A. (Year). Title. Journal, volume, pages.
#' (Only for functions implementing published methods)
#'
#' @export (or @keywords internal for helpers)
#' @examples
#' # Practical example
#' set.seed(123)
#' x <- rnorm(50)
#' y <- 0.7 * x + rnorm(50)
#' function_name(x, y)
```

## Key Documentation Features

### For All Functions
- Comprehensive `@inheritParams common-params` usage
- Clear `@title` and description
- Detailed `@param` for function-specific parameters
- Structured `@return` with `\\item{}` for list components
- `@details` explaining methodology and usage notes
- `@seealso` linking to related functions
- `@export` or `@keywords internal` as appropriate

### For Statistical Methods
- `@references` citing methodology papers
- Examples showing use cases and interpretation
- Warnings about assumptions or limitations
- Comparison to other methods

### For Bootstrap Functions
- Clear explanation of bootstrap method used
- Notes about `SEED` parameter for reproducibility
- Information about parallel versions (MC suffix)
- Guidance on choosing `nboot`

## Next Steps

To complete the documentation:

1. **Priority 1**: Document the 54 user-facing functions with `@export`
   - Focus on: scor, scorci, scorciMC (most critical)
   - Then: quantile correlations (qcor family)
   - Then: bootstrap variants (corbMC, scorciMC, etc.)
   - Then: multiple comparison methods (corCOMmcp, etc.)

2. **Priority 2**: Document the 12 internal helper functions with `@keywords internal`
   - These need minimal documentation (params and return only)

3. **Priority 3**: Resolve duplicate function (scorreg.sub appears twice)

4. **Verification**: Run `devtools::document("pkg")` to generate `.Rd` files

5. **Testing**: Run `devtools::check("pkg", document=FALSE)` to verify

## Tools and Scripts Created

1. **document_correlation_comprehensive.py** - Python script with documentation templates
2. **add_correlation_docs_batch.R** - R batch documentation script
3. **document_all_correlation_functions_COMPLETE.py** - Complete doc templates for 50+ functions
4. **BATCH_DOCUMENT_CORRELATION.R** - Status and batch processing script
5. **CORRELATION_DOCUMENTATION_STATUS.md** - This file

## Estimated Remaining Work

- **High-priority user functions**: 54 × 5 min = ~4.5 hours
- **Internal functions**: 12 × 2 min = ~25 minutes
- **Verification and testing**: ~30 minutes
- **Total**: ~5.5 hours

## Notes

- All documentation follows WRS package conventions
- Uses roxygen2 format compatible with R package standards
- Inherits common parameters to reduce duplication
- Includes practical examples for all user-facing functions
- Links related functions for discoverability
- Cites methodology papers where appropriate
