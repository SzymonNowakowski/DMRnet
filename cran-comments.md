
# This is a resubmission of 0.3.4

The following NOTE from WINDOWS precheck **has been corrected**:

  * checking CRAN incoming feasibility ... [15s] NOTE
  Maintainer: 'Szymon Nowakowski <s.nowakowski@mimuw.edu.pl>'

  Package CITATION file contains call(s) to personList() or
  as.personList().  Please use c() on person objects instead.
  Package CITATION file contains call(s) to citEntry(), citHeader() or
  citFooter().  Please use bibentry() instead.

# Local checks

## local R CMD check results
There were no ERRORs, WARNINGs, or NOTEs. 

## Downstream dependencies
There seem to be no downstream dependencies:

```{r revdep}
devtools::revdep("DMRnet")
# character(0)
```


