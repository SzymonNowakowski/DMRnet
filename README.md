<!-- badges: start -->
[![GitHub version](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/SzymonNowakowski/DMRnet/master/.version.json&style=flat&logo=github)](https://github.com/SzymonNowakowski/DMRnet)
[![CRAN version](https://img.shields.io/cran/v/DMRnet?logo=R)](https://cran.r-project.org/package=DMRnet)
[![downloads](https://cranlogs.r-pkg.org/badges/DMRnet)](https://cran.r-project.org/package=DMRnet)
<!-- badges: end -->


# DMRnet

DMRnet (Delete or Merge Regressors) is a suit of algorithms for linear and logistic model selection with high-dimensional data (i.e. the number of regressors may exceed the number of observations). The predictors can be continuous or categorical. The selected model consists of a subset of numerical regressors and partitions of levels of factors. 

For information on how to get started using DMRnet, see our [getting started vignette](https://cran.r-project.org/web/packages/DMRnet/vignettes/getting-started.html).

## Installing `DMRnet` package

To install the development package version (currently: 0.3.1.9002) please execute
```
library(devtools)
devtools::install_github("SzymonNowakowski/DMRnet")
```

Alternatively, to install the current stable CRAN version (currently: 0.3.1) please execute

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

## Purpose of the `testing_branch`

The `testing_branch` GitHub branch serves double purpose.

### Consistency with previous package versions

The first purpose is testing agreement of the results obtained for the new versions and the older version v. 0.2.0.

1. In the summer of 2021 massive experiments were performed with DMRnet v. 0.2.0 
   and with new implementation of GLAMER and two new external cross validation routines, 
   which got later implemented into DMRnet. In the scope of the experiments 
   model dimension and prediction error was calculated in the 4 data sets in 200 independent runs.

   The two new cross validation routines were
   
   - `cvg` which stands for "(c)ross (v)alidation with models indexed with (G)IC" 
      and since 0.3.1 version of the package is available in DMRnet by simply passing `indexation.mode=GIC"` 
      in a call to `cv.DMRnet`.
      
   - `cv+sd` which is a regular cross validation with models indexed by model dimension, but
      returning (instead of the best model) the smallest model falling under the upper 
      curve of a prediction error plus one standard deviation. Since 0.3.1 version of the package it is now 
      available in `predict` for inference, use it by passing `md="df.1se"` instead of the default `md="df.min"`.
   
   Throughout the experiments, the models were also selected with a regular GIC
   criterion, both in DMRnet and in GLAMER. The GIC criterion was available in v. 0.2.0 for DMRnet and since 0.3.1 version of the package for both DMRnet and for GLAMER. 
   
   The above possibilities give rise to many combinations of which 4 were selected for the final consideration:
   
   - `cv+sd.GLAMER` - GLAMER run with `cv+sd`.
   
   - `cvg.DMRnet` - DMRnet run with `cvg`.
   
   - `gic.GLAMER` - GLAMER run with GIC criterion to select best model.
   
   - `gic.DMRnet` - DMRnet run with GIC criterion to select best model.
   
   The resulting model dimension and prediction error data was retained. 
   It creates opportunity to compare the expected distibutions of model dimension and prediction error 
   with the actual distributions for new package versions
   (we would expect to get the same results with the new 
   versions of the package as were obtained in the summer of 2021 with v. 0.2.0).
   
   Care was taken to reproduce exactly the same results, as an example: the new calculations in new package versions
   were started from the same seed values. However, calculations in v. 0.2.0 were much less stable and required
   `try ... catch` block to complete. Thus, from the first failed (and restarted) calculation in v. 0.2.0 the 
   new and old results start to diverge. It does not eliminate possibility to compare the 
   resulting distributions.
   
   The comparison charts are presented below:
   
   - [Adult data set](https://archive.ics.uci.edu/ml/datasets/Adult) - binomial response
   
     ![Adult](https://github.com/SzymonNowakowski/DMRnet/blob/testing_branch/result_adult.svg)
     
   - [Promoter data set](https://archive.ics.uci.edu/ml/datasets/Molecular+Biology+%28Promoter+Gene+Sequences%29) - binomial response
   
     ![Promoter](https://github.com/SzymonNowakowski/DMRnet/blob/testing_branch/result_promoter.svg)
     
   - [Insurance data set](https://www.kaggle.com/c/prudential-life-insurance-assessment/data) - gaussian response
   
     ![Insurance](https://github.com/SzymonNowakowski/DMRnet/blob/testing_branch/result_insurance.svg)
     
     **Results compared only for 40 initial results (out of 200 available) due to 
     [`LAPACK` bug in `dgesdd` routine](https://github.com/Reference-LAPACK/lapack/issues/672) 
     thus higher variance is to be expected**
     
   - Antigua data set - Averages by block of yields for the Antigua Corn 
     data - available 
     in [DAAG package](https://cran.rstudio.com/web/packages/DAAG/index.html) - gaussian response
     
     ![Antigua](https://github.com/SzymonNowakowski/DMRnet/blob/testing_branch/result_antigua.svg)
   
2. In the beginning of 2022 massive experiments were performed with DMRnet v. 0.2.0 
   and with new implementation of GLAMER and two new external cross validation routines, which got later
   implemented into DMRnet since 0.3.1 version of the package. In the scope of the experiments model dimension and 
   prediction error was calculated in 252 synthetic experiments, each consisting of 200 independent runs.
   
   The possible choices of an algorithm and a cross validation routine (discussed in point 1. above) give rise 
   to many combinations of which 3 were selected for the final consideration:
   
   - `cv+sd.GLAMER` - GLAMER run with `cv+sd`.
   
   - `cvg.GLAMER` - GLAMER run with `cvg`.
   
   - `cvg.DMRnet` - DMRnet run with `cvg`.
   
   
   The resulting model dimension and prediction error data was retained. 
   It creates opportunity to compare the expected distibutions of model dimension and prediction error 
   with the actual distributions for new package versions
   (we would expect to get the same results with the new 
   versions of the package as were obtained in the beginning of 2022 with v. 0.2.0).
   
   ![High Dimensional Simulations](https://github.com/SzymonNowakowski/DMRnet/blob/testing_branch/result_high_dimensional_simulation.svg)
       
### Torture testing 

The second purpose is testing that all previously identified very hard cases pass in new versions of the package.
Many instabilities in v. 0.2.0 were removed *en bloc* by adding new strategy for handling missing factors 
and adding regularization. However, not all bugs were removed that way and the remaining 
bugs have some additional test files dedicated to testing that the particular bug remains fixed.

To this end the following test cases were identified:

- `hard_case_DMRnet_insurance.R` - a compilation of tests for various previously 
 identified bugs in DMRnet with data coming from 
 [Insurance data set](https://www.kaggle.com/c/prudential-life-insurance-assessment/data). 
 Please consult the comments in the file for more detailed information.

- `hard_case_LAPACK_SVD_insurance.R` - the outstanding (not fixed) cases of GLAMER and DMRnet computation 
 failure resulting from 
 [`LAPACK` bug in `dgesdd` routine](https://github.com/Reference-LAPACK/lapack/issues/672) with
 data isolated from 
 [Insurance data set](https://www.kaggle.com/c/prudential-life-insurance-assessment/data). 

- `hard_case_NA_insurance.R` - cases in which NA values were returned for GLAMER, DMRnet and DMR for [Insurance data set](https://www.kaggle.com/c/prudential-life-insurance-assessment/data) because design matrices involved were not of full rank. [PR#27](https://github.com/SzymonNowakowski/DMRnet/pull/27) and [PR#28](https://github.com/SzymonNowakowski/DMRnet/pull/28) fixed those issues.

- `hard_case_GLAMER_promoter.R` - a problem with group constrained identified with GLAMER in
 data isolated from 
 [Promoter data set](https://nbviewer.org/github/SzymonNowakowski/DMRnet/blob/testing_branch/result_promoter.pdf).
 The cause was that `grpreg` didn't observe the so called *group constraint* for one of the columns 
 and returned two non-zero betas and one beta equal to zero.

- `hard_case_SOSnet.R` - a problem in SOSnet with artificial data set 
 (in which a matrix is close to - but not exactly - singular) resulting in 
 numerical instabilities with QR calculations.

- `hard_case_CV_airbnb.R` - a problem in model-dimension indexed CV resulting in mismatching model dimension in rare cases of different number of predictors in-between folds during the cross validation. The side effect of this bug has been larger than optimal models returned in case of this mismatch. In the extreme (even rarer) cases the returned model dimension exceeded the number of a computed `dmr.fit` models, which in effect led to incorrect `predict()` behavior. Fixed in [PR#34](https://github.com/SzymonNowakowski/DMRnet/pull/34) and [PR#35](https://github.com/SzymonNowakowski/DMRnet/pull/35).
