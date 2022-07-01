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
GLAMER stands for Group LAsso MERger and it is a new (simplified in relation to DMRnet) algorithm for which we prove partition selection consistency. It is the first result of that kind for high dimensional scenario. The relevant paper with algorithm description is the following: [Szymon Nowakowski, Piotr Pokarowski and Wojciech Rejchel. 2021. “Group Lasso Merger for Sparse Prediction with High-Dimensional Categorical Data.” arXiv:2112.11114](https://arxiv.org/abs/2112.11114)

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
1. The pivots were added in SOSnet for gaussian families and as a consequence a non-full rank data matrix was handled correctly to some extent (see the next point)
2. Numerical instability of QR decomposition in SOSnet for gaussian families (when rank is *much* lower than columns provided) may have resulted in crashes in QR decomposition, thus a `try`-`catch` clause was used to recalculate QR decomposition but only for pivoted columns within rank in case of failure of the original QR attempt.
3. Adding the regularizations in gaussian family in `DMR` (note that this execution path is used in `DMRnet` too) in a form of adding an infinitesimally small diagonal matrix to `R` after QR decomposition calculation

Other improvements are the following:
1. Fixing how the parameters get passed in `plot` family of functions.
2. Updating how `coef` presents parameters: only non-zero coefficients get returned.
4. Fixing the error in `DRMnet` and cross validation in handling data with a single two-factor column, by adding `drop=FALSE` statement.
5. Fixing negative values in binomial case coming from `w_stats` from numerical instability when `Kan` was very close to 0 and `Var` is not symmetric, but `w_stats` assumes symmetric `Var` (problem observed in `DMRnet` for binomial family in [Insurance data set](https://www.kaggle.com/c/prudential-life-insurance-assessment/data) - see hard_case_DMRnet_insurance.R` test file in `testing_branch`).
6. There have been cases of `grpreg` not observing a group constraint (i.e. a condition that either all betas are zero, or all betas are non-zero within a group) in [Promoter data set](https://archive.ics.uci.edu/ml/datasets/Molecular+Biology+%28Promoter+Gene+Sequences%29) - see `hard_case_DMRnet_promoter.R` test file in `testing_branch`. Some betas that belonged to groups > 0 were not strictly > 0. It was problematic in GLAMER only, as DMRnet recalculated the t-statistics and was not constrained by initial beta values. It was fixed by adding a small constant to all betas in groups with at least one non-zero beta in GLAMER.


### Weight parameterization
This remains to be introduced to GLAMER and DMRnet algorithms in future versions >0.3.0.

### Remaining issues

The only outstanding (not fixed) cases of DMRnet computation failure known to me at present are gathered in `hard_case_GLAMER_insurance.R`  in  `testing_branch` in [Insurance data set](https://www.kaggle.com/c/prudential-life-insurance-assessment/data). 

## Purpose of the `testing_branch`

The `testing_branch` GitHub branch serves double purpose.

### Consistency with previous package versions

The first purpose is testing agreement of the results obtained for the new version v. 0.3.0 and the older version v. 0.2.0.

1. In the summer of 2021 massive experiments were performed with DMRnet v. 0.2.0 
   and with new implementation of GLAMER and two new external cross validation routines, 
   which got later implemented into DMRnet v. 0.3.0. In the scope of the experiments 
   model dimension and prediction error was calculated in the 4 data sets in 200 independent runs.

   The two new cross validation routines were
   - `cvg` which stands for "(c)ross (v)alidation with models indexed with (G)IC" 
      and in the v. 0.3.0 is available in DMRnet by simply passing `indexation.mode=GIC"` 
      in a call to `cv.DMRnet`.
   - `cv+sd` which is a regular cross validation with models indexed by model dimension, but
      returning (instead of the best model) the smallest model falling under the upper 
      curve of a prediction error plus one standard deviation. In v. 0.3.0 it is now 
      available in `predict` for inference, use it by passing `md="df.1se"` instead of the default `md="df.min"`.
   
   Throughout the experiments, the models were also selected with a regular GIC
   criterion, both in DMRnet and in GLAMER. The GIC criterion was available in v. 0.2.0 for DMRnet and in v. 0.3.0 both for DMRnet and for GLAMER. 
   
   The above possibilities give rise to many combinations of which 4 were selected for the final consideration:
   - `cv+sd.GLAMER` - GLAMER run with `cv+sd`.
   - `cvg.DMRnet` - DMRnet run with `cvg`.
   - `gic.GLAMER` - GLAMER run with GIC criterion to select best model.
   - `gic.DMRnet` - DMRnet run with GIC criterion to select best model.
   
   The resulting model dimension and prediction error data was retained. 
   It creates opportunity to compare the expected distibutions of model dimension and prediction error 
   with the actual distributions we get with a new package
   (we would expect to get the same results with the new 
   version v. 0.3.0 as were obtained in the summer of 2021 with v. 0.2.0).
   
   Care was taken to reproduce exactly the same results, as an example: the new calculations in v. 0.3.0
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
   
     **Results unavailable due to 
     [`LAPACK` bug in `dgesdd` routine](https://github.com/Reference-LAPACK/lapack/issues/672)**
   - Antigua data set - Averages by block of yields for the Antigua Corn 
     data - available 
     in [DAAG package](https://cran.rstudio.com/web/packages/DAAG/index.html) - gaussian response
     
     ![Antigua](https://github.com/SzymonNowakowski/DMRnet/blob/testing_branch/result_antigua.svg)
   
2. In the beginning of 2022 massive experiments were performed with DMRnet v. 0.2.0 
   and with new implementation of GLAMER and two new external cross validation routines, which got later
   implemented into DMRnet v. 0.3.0. In the scope of the experiments model dimension and 
   prediction error was calculated in 252 synthetic experiments, each consisting of 200 independent runs.
   
   The possible choices of an algorithm and a cross validation routine (discussed in point 1. above) give rise 
   to many combinations of which 3 were selected for the final consideration:
   - `cv+sd.GLAMER` - GLAMER run with `cv+sd`.
   - `cvg.GLAMER` - GLAMER run with `cvg`.
   - `cvg.DMRnet` - DMRnet run with `cvg`.
   
   
   The resulting model dimension and prediction error data was retained. 
   It creates opportunity to compare the expected distibutions of model dimension and prediction error 
   with the actual distributions we get with a new package
   (we would expect to get the same results with the new 
   version v. 0.3.0 as were obtained in the beginning of 2022 with v. 0.2.0).
   
   **Results unavailable *yet* because of the scale of computations required to complete calculations**
       
### Very hard cases

The second purpose is testing that all previously identified crux cases pass in v. 0.3.0. 
Many instabilities in v. 0.2.0 were removed *en bloc* with adding new strategy for handling missing factors 
and adding regularization. However, not all bugs were removed that way and the remaining 
bugs have some additional test files dedicated to testing that the particular bug remains fixed.

To this end the following test cases were identified:
- `hard_case_DMRnet_insurance.R` - a compilation of tests for various previously 
 identified bugs in DMRnet with data coming from 
 [Insurance data set](https://www.kaggle.com/c/prudential-life-insurance-assessment/data). 
 Please consult the comments in the file for more detailed information.
- `hard_case_GLAMER_insurance.R` - the only outstanding (not fixed) cases of GLAMER computation 
 failure known to me at present. One of them results from 
 [`LAPACK` bug in `dgesdd` routine](https://github.com/Reference-LAPACK/lapack/issues/672) with
 data isolated from 
 [Insurance data set](https://www.kaggle.com/c/prudential-life-insurance-assessment/data). 
- `hard_case_GLAMER_promoter.R` - a problem with group constrained identified with GLAMER in
 data isolated from 
 [Promoter data set](https://nbviewer.org/github/SzymonNowakowski/DMRnet/blob/testing_branch/result_promoter.pdf).
 The cause was that `grpreg` didn't observe the so called *group constraint* for one of the columns 
 and returned two non-zero betas and one beta equal to zero.
- `hard_case_SOSnet.R` - a problem in SOSnet with artificial data set 
 (in which a matrix is close to - but not exactly - singular) resulting in 
 numerical instabilities with QR calculations fixed in v. 0.3.0.

