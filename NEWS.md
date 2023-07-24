
# DMRnet 0.3.4.9002

- work started on group selection with no partition selection
- updated `.onAttach` with citation info

# DMRnet 0.3.4.9001

- hierarchical clustering based on segmentation of 1D points with hclust1d package for PDMR and GLAMER 

# DMRnet 0.3.4.9000

- hierarchical clustering based on segmentation of 1D points for PDMR and GLAMER 

# DMRnet 0.3.4

- PDMR added

# DMRnet 0.3.3

- Updating documentation to reflect that for inference the matrix must be provided without the intercept column
- Fixing incorrect inference for columns not in first-factors-then-numerics sequence
- Candidate fix for warnings related to [issue #33](https://github.com/SzymonNowakowski/DMRnet/issues/33)
- Setting a default `nlambda` value to 100 in `cv.DMRnet()` ([issue #41](https://github.com/SzymonNowakowski/DMRnet/issues/41))
- Fixing problems when few or no models are available in cross validation
- Fixing problems when few or no models are available after `grpreg` ([issue #39](https://github.com/SzymonNowakowski/DMRnet/issues/39))
- Refactoring the code to distribute the no-model fixes from gaussian `DMRnet` into binomial family `DMRnet` and both `SOSnet`s and `GLAMER`s
- Fixing a few other minor long-standing issues in GIC-indexed cross validation
- Fixing df.1se in GIC-indexed cross validation for binomial GLAMER
- Fixing incorrect loglik calculation for the first (largest) model in binomial family, for both `GLAMER` (`-Inf`) and `DMRnet` (incorrect)

# DMRnet 0.3.2

- Improved readability of a getting-started vignette
- Fixed a bug in model-indexed cross validation related to folds with different model sizes
- Added df.1se to GIC-indexed cross validation
- Improved CV plots with `df.1se` model
- Improved readability of README on CRAN ([issue #32](https://github.com/SzymonNowakowski/DMRnet/issues/32))
- Welcome message on package load added

# DMRnet 0.3.1

## Major changes

- GLAMER (Group LAsso MERger) algorithm added
- Rebuild of the existing cross validation routine (models indexed by model dimension)
- New cross validation routine added (models indexed by GIC) and made a new default
- Handling of mismatched factor levels restructured

## Stability improvements

- Adding pivots in SOSnet for gaussian families
- Adding a `try`-`catch` clause to recalculate QR decomposition for pivoted columns within rank in case of failure
- Adding the regularization in gaussian family in `DMR` direct calls and calls via `DMRnet`
- fixing NA values resulting in attempting non-penalized regression in case of non-full rank design matrices

## Other improvements

- Fixing how the parameters get passed in `plot` family of functions
- Returning only non-zero coefficients in `coef`
- Fixing `DRMnet` and cross validation in handling data with a single two-factor column
- Fixing negative Wald statistics in binomial case
- Workaround in GLAMER for [cases of `grpreg` not observing a group constraint](https://github.com/pbreheny/grpreg/issues/54)
