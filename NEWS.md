
# DMRnet 0.3.1.9001

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
