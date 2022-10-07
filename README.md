<!-- badges: start -->
[![GitHub version](https://img.shields.io/github/r-package/v/SzymonNowakowski/DMRnet?color=yellowgreen&label=GitHub&logo=github)](https://github.com/SzymonNowakowski/DMRnet)
[![CRAN version](https://img.shields.io/cran/v/DMRnet?logo=R)](https://cran.r-project.org/package=DMRnet)
[![downloads](https://cranlogs.r-pkg.org/badges/DMRnet)](https://cran.r-project.org/package=DMRnet)
<!-- badges: end -->


# DMRnet

DMRnet (Delete or Merge Regressors) is a suit of algorithms for linear and logistic model selection with high-dimensional data (i.e. the number of regressors may exceed the number of observations). The predictors can be continuous or categorical. The selected model consists of a subset of numerical regressors and partitions of levels of factors. 

For information on how to get started using DMRnet, see our [getting started vignette](https://cran.r-project.org/package=DMRnet/vignettes/getting-started.html).

## Installing `DMRnet` package

To install the development package version please execute
```
library(devtools)
devtools::install_github("SzymonNowakowski/DMRnet")
```

Alternatively, to install the current stable CRAN version please execute

```
install.packages("DMRnet")
```

After that, you can load the installed package into memory with a call to `library(DMRnet)`.


## Features

### GLAMER added

GLAMER was added in 0.3.1 version of the package. GLAMER stands for Group LAsso MERger and it is a new (simplified in relation to DMRnet) algorithm for which we prove partition selection consistency. It is the first result of that kind for high dimensional scenario. The relevant paper with algorithm description is the following: [Szymon Nowakowski, Piotr Pokarowski and Wojciech Rejchel, 2021. “Group Lasso Merger for Sparse Prediction with High-Dimensional Categorical Data.” arXiv:2112.11114](https://arxiv.org/abs/2112.11114)

To use GLAMER pass `algorithm="glamer"` in a call to `DMRnet()` or cross validation routine `cv.DMRnet()`. GLAMER is not supported in `DMR`.

### Two cross validation routines

A new cross validation routine was introduced to improve the computed model quality. It indexes models by GIC. The method was proposed and first implemented for *gaussian* family by Piotr Pokarowski. Since 0.3.1 version of the package it has been built into `DMRnet` for both *gaussian* and *binomial* families. 

All in all, the cross validation features in the package are the following:

1. Models can be indexed by GIC or by model dimension. The relevant setting is selected with, respectively, `indexation.mode="GIC"` or `indexation.mode="dimension"` parameter in a call to `cv.DMRnet()`. The setting that indexes models by GIC has been the default since 0.3.1 version of the package.

2. The net of lambda values is first calculated from the full data set and then this net is used for particular cross validation folds. The motivation behind this change is to stabilize the results.

3. Apart from `df.min`, which is the model with minimal cross-validated error, the routines now return `df.1se` which is the smallest model falling under the upper curve of a prediction error plus one standard deviation. It can be used in `predict()` for inference by passing `md="df.1se"` instead of the default `md="df.min"`.

4. Cross validation handles the mismatched factor levels in a way that minimizes incorrect behavior (see Section [Handling of mismatched factor levels](#handling-of-mismatched-factor-levels)).

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

   Then, if treated naively, `DMRnet(...)` when constructing a model would throw an error,
   because we would end up with `NaN` values in a column dedicated to this superfluous factor level (to be exact, it would happen when columns get normalized).

   The solution to that is straightforward. Before the model gets constructed in `DMRnet` we recalculate the factor level set, `C_new`. Then `C_new=A`.

   **SOLVED**
   
1. `B` does not contain a level(s) present in `A`.

   (sample case: we did sample to `Xtr` the single Dutch national from the [Insurance data set](https://www.kaggle.com/c/prudential-life-insurance-assessment/data), and he is not present in `Xte`,
   because there is only one instance of Dutch national in the whole Insurance data set).
   As a result `predict(...)` would throw an error, because expanded model-matrix dimensions would be conflicting.

   The solution is simple here, too: in constructing a model make a note about the true `A` set (technically, it gets stored into `levels.listed` variable in a model)
   and then in `predict(...)` assign the levels of `Xte` to be equal to `A`. Only then create the model-matrix.

   **SOLVED**
   
1. `B` contains a factor level(s) not present in `A`, AND we are doing CV, so we have access to `Xtr`.

   The solution is to remove the rows with levels that are going to cause problems later in `predict(...)` from `Xte` before the prediction.
   The other solution would be to predict using unknown.factor.levels="NA" flag and then eliminate the `NAs` from comparisons (this solution is NOT used at present)

   **SOLVED**

1. `B` contains a factor level(s) not present in `A`, AND we are NOT doing CV, so we have no access to `Xtr`.

   This case is problematic because this situation gets identified too late - we are already in `predict(...)`.
   At this point, only the model created by `DMRnet(...)` function
   (which got passed into `predict(...)` function) is known.
   We cannot perform inference and we cannot perform any imputation for the problematic data point, either 
   (we don't know `Xtr` and have no access to it).
   
   All that remains is to throw an error (when `unknown.factor.levels="error"`, the default) OR
   eliminate the problematic rows, predict, and then replenish the result with `NAs` in place 
   of problematic values (when `unknown.factor.levels="NA"`).

   None of this solutions is fully satisfactory, thus this case remains **PROBLEMATIC**.

### Stability improvements

Generally speaking, matrix rank in real world scenarios is more a numerical concept than a mathematical concept and its value may differ depending on a threshold. Thus various kinds of problems result from data matrices close to singular. Since 0.3.1 version of the package, the work has been devoted to improve stability of computations with such ill-defined matrices. See `NEWS.md` for more information on detailed stability improvements.


### Weight parameterization

This remains to be introduced to GLAMER and DMRnet algorithms in future versions.
