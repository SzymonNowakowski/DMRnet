
# This is a 2nd resubmission of 0.3.4

- The following NOTE from WINDOWS precheck **has been corrected**:

  * checking CRAN incoming feasibility ... [15s] NOTE
  Maintainer: 'Szymon Nowakowski <s.nowakowski@mimuw.edu.pl>'

  Package CITATION file contains call(s) to personList() or
  as.personList().  Please use c() on person objects instead.
  Package CITATION file contains call(s) to citEntry(), citHeader() or
  citFooter().  Please use bibentry() instead.
  
- The folowing NOTE from precheck **has been corrected** by reducing examples by roughly 25%:

  Flavor: r-devel-windows-x86_64
  Check: examples, Result: NOTE
    Examples with CPU (user + system) or elapsed time > 10s
           user system elapsed
    DMRnet   10   0.17   10.21

# Local checks

## local R CMD check results
There were no ERRORs, WARNINGs, or NOTEs. 

## Downstream dependencies
There seem to be no downstream dependencies:

```{r revdep}
devtools::revdep("DMRnet")
# character(0)
```


