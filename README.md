# DMRnet
This is a fork of the CRAN R package repository for DMRnet — Delete or Merge Regressors Algorithms for Linear and Logistic Model Selection and High-Dimensional Data. The purpose of the fork is to maintain and develop new package versions.

## Installing `DMRnet` package

To install the development package version (currently: 0.3.0.9000) please execute
```
library(devtools)
devtools::install_github("SzymonNowakowski/DMRnet")
```

Alternatively, to install the current stable CRAN version (currently: 0.2.0) please execute

```
install.packages("DMRnet")
```

After that, you can load the installed package into memory with a call to `library(DMRnet)`.


## Changes in DMRnet v. 0.3.0

### GLAMER added
GLAMER stands for Group LAsso MERger and it is a new (simplified in relation to DMRnet) algorithm for which we prove partition selection consistency. It is the first result of that kind for high dimensional scenario. The relevant paper with algorithm description is the following: [Szymon Nowakowski, Piotr Pokarowski and Wojciech Rejchel, 2021. “Group Lasso Merger for Sparse Prediction with High-Dimensional Categorical Data.” arXiv:2112.11114](https://arxiv.org/abs/2112.11114)

To use GLAMER pass `algorithm="glamer"` in a call to `DMRnet` or cross validation routine. GLAMER is not supported in `DMR`.

### Rebuild of the existing cross validation routine (models indexed by model dimension)
The following changes were introduced to improve the existing model dimension-indexed cross validation:
1. The net of lambda values is first calculated from the full data set and then this net is used for particular cross validation folds. 
2. Apart from `df.min` which is the number of parameters of the model with minimal cross-validated error, the routine now returns `df.1se` which is the number of parameters of the smallest model falling under the upper curve of a prediction error plus one standard deviation. It can be used in `predict` for inference by passing `md="df.1se"` instead of the default `md="df.min"`.
3. Corrected handling of mismatched factor levels (see Section [Handling of mismatched factor levels](#handling-of-mismatched-factor-levels)).

To use the existing model dimension-indexed cross validation routine pass `indexation.mode="dimension"` in a call to `cv.DMRnet`.
### New cross validation routine (models indexed by GIC)
A new alternative cross validation routine was introduced to improve the existing model quality. It indexes models by GIC. The method was proposed and first implemented for gaussian family by Piotr Pokarowski. The new method additional features are:
1. The net of lambda values is first calculated from the full data set and then this net is used for particular cross validation folds. 
2. Correct way of handling the mismatched factor levels (see Section [Handling of mismatched factor levels](#handling-of-mismatched-factor-levels)).

To use the new GIC-indexed cross validation routine pass `indexation.mode=GIC"` in a call to `cv.DMRnet`.

The new cross validation routine (models indexed by GIC) is the default starting from version >=0.3.0.
### Handling of mismatched factor levels
The new treatment of factors in cross validation/`predict` and in `DMRnet`/`predict` pairs is based on the following analysis:

Let us assume that
- `Xtr` is training data in cross validation or in a regular call via `DMRnet`->`model`
- `Xte` is test data in cross validation or in a regular call via `model`->`predict`

Without loss of generality, let us consider `Xtr` and `Xte` to be one column only, with factors.

Let us also consider the following definitions:
- `A` is a true set of all factor levels in `Xtr`
- `B` is a true set of all factor levels in `Xte`
- `C=levels(Xtr)` is a set of factor levels in original data that `Xtr` originates from, but it is still assigned to `Xtr` via the `levels()` function. As a rule, when taking subsets, `R` does not eliminate redundant factors, so let us note that `C` is a superset of `A`.

There are 4 classes of problems:
1. `C` is a strict superset of `A`.

   Then if treated naively, `DMRnet(...)` when constructing a model would throw an error,
   because we would end up with `NaN` values in a column dedicated to this superflous factor level (to be exact, it would happen when a columns gets normalized).

   The solution to that is very simple. Before the model gets constructed in `DMRnet` we recalculate the factor level set, `C_new`. Then `C_new=A`.

   **SOLVED**
2. `B` does not contain a level(s) present in `A`.

   (sample case: we did sample to `Xtr` the single Dutch national from the [Insurance data set](https://www.kaggle.com/c/prudential-life-insurance-assessment/data), and he is not present in `Xte`,
   because there is only one instance of Dutch national in the whole Insurance data set).
   As a result `predict(...)` would throw an error, because expanded model-matrix dimensions would be conflicting.

   The solution is simple here, too: in constructing a model make a note about true `A` set (it is stored in `levels.listed` variable in a model)
   and then in `predict(...)` assign the levels of `Xte` to be equal to `A`. Only then create the model-matrix.

   **SOLVED**
3. `B` contains a factor level(s) not present in `A`, AND we are doing CV, so we have access to `Xtr`.

   The solution is to remove the rows with levels that are going to cause problems later in `predict(...)` from `Xte` before the prediction.
   The other solution would be to predict using unknown.factor.levels="NA" flag and then eliminate the `NAs` from comparisons (this solution is NOT used at present)

   **SOLVED**

4. `B` contains a factor level(s) not present in `A`, AND we are NOT doing CV, so we have no access to `Xtr`.

   This case is problematic because this situation gets identified too late - we are already in `predict(...)`.
   At this point, only the model created by `DMRnet(...)` function and passed to `predict(...)` is known.
   We cannot perform inference and we cannot perform any imputation for the problematic data point, either (we don't know `Xtr` and have no access to it).
   
   All that remains is to throw an error (when `unknown.factor.levels="error"`, the default) OR
   eliminate the problematic rows, predict, and then replenish the result with `NAs` in place of problematic values (when `unknown.factor.levels="NA"`).

   And this solution is not fully satisfactory, thus this case remains **PROBLEMATIC**.

### Stability improvements


Generally speaking, matrix rank in real world scenarios is more a numerical concept than a mathematical concept and its value may differ depending on a treshold. Thus various kinds of problems result from data matrices close to singular:
1. The pivots were added in SOSnet for gaussian families and as a consequence a non-full rank data matrix was handled correctly to some extent (see the next point).
2. Numerical instability of QR decomposition in SOSnet for gaussian families (when rank is *much* lower than columns provided) may have resulted in crashes in QR decomposition, thus a `try`-`catch` clause was used to recalculate QR decomposition but only for pivoted columns within rank in case of failure of the original QR attempt.
3. Adding the regularizations in gaussian family in `DMR` (note that this execution path is used in `DMRnet` too) in a form of adding an infinitesimally small diagonal matrix to `R` after QR decomposition calculation.
4. Fixing NA values resulting in attempting non-penalized regression in case of non-full rank design matrices.

Other improvements are the following:
1. Fixing how the parameters get passed in `plot` family of functions.
2. Updating how `coef` presents parameters: only non-zero coefficients get returned.
4. Fixing the error in `DRMnet` and cross validation in handling data with a single two-factor column, by adding `drop=FALSE` statement.
5. Fixing negative values in binomial case coming from `w_stats` from numerical instability when `Kan` was very close to 0 and `Var` is not symmetric, but `w_stats` assumes symmetric `Var` (problem observed in `DMRnet` for binomial family in [Insurance data set](https://www.kaggle.com/c/prudential-life-insurance-assessment/data) - see hard_case_DMRnet_insurance.R` test file in `testing_branch`).
6. There have been [cases of `grpreg` not observing a group constraint](https://github.com/pbreheny/grpreg/issues/54) (i.e. a condition that either all betas are zero, or all betas are non-zero within a group) in [Promoter data set](https://archive.ics.uci.edu/ml/datasets/Molecular+Biology+%28Promoter+Gene+Sequences%29) - see `hard_case_DMRnet_promoter.R` test file in `testing_branch`. Some betas that belonged to groups > 0 were not strictly > 0. It was problematic in GLAMER only, as DMRnet recalculated the t-statistics and was not constrained by initial beta values. It was fixed in GLAMER by adding a small constant to all betas in groups with at least one non-zero beta.


### Weight parameterization
This remains to be introduced to GLAMER and DMRnet algorithms in future versions >0.3.0.
