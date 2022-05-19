# DMRnet
This is a fork of the CRAN R package repository for DMRnet — Delete or Merge Regressors Algorithms for Linear and Logistic Model Selection and High-Dimensional Data. The purpose of the fork is to maintain and develop new package versions.

## Changes in DMRnet v. 0.3.0

### GLAMER added
GLAMER stands for Group LAsso MERger and it is a new (simplified in relation to DMRnet) algorithm for which we prove partition selection consistency. It is the first result of that kind for high dimensional scenario. The relevant paper with algorithm description is the following: [Szymon Nowakowski, Piotr Pokarowski and Wojciech Rejchel. 2021. “Group Lasso Merger for Sparse Prediction with High-Dimensional Categorical Data.” arXiv:2112.11114](https://arxiv.org/abs/2112.11114)

To use GLAMER pass `algorithm="glamer"` in a call to `DMRnet` or cross validation routine. GLAMER is not supported in `DMR`.

### Rebuilt the existing cross validation routine (models indexed by model dimension)
The following changes were introduced to improve the existing model dimension-indexed cross validation:
1. The lambda net is first calculated from the full data and than this net is used for particular cross validation folds. 
2. Apart from `df.min` which is the number of parameters of the model with minimal cross-validated error, the routine now returns `df.1se` which is the number of parameters of the smallest model falling under the upper curve of a prediction error plus one standard deviation. It can be used in `predict` for inference by passing `md="df.1se"` instead of the default `md="df.min"`.

To use the existing model dimension-indexed cross validation routine pass `indexation.mode="dimension"` in a call to `cv.DMRnet`.
### New cross validation routine  (models indexed by GIC)
A new alternative cross validation routine was introduced to improve the existing model quality. It indexes models by GIC. The method was proposed and first implemented for gaussian family by Piotr Pokarowski. The new method has the additional feature, that the lambda net is first calculated from the full data and than this net is used for particular cross validation folds. 

To use the new GIC-indexed cross validation routine pass `indexation.mode=GIC"` in a call to `cv.DMRnet`.
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

### Weights parametrized
This remains to be done
