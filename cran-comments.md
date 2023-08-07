
# `check_win_devel()` and `check_rhub()` checks

- All words from this NOTE are spelled correctly (they are the author name/conference title/publisher of Springer publication)

  * checking CRAN incoming feasibility ... [11s] NOTE
    Maintainer: 'Szymon Nowakowski <s.nowakowski@mimuw.edu.pl>'
    
    Possibly misspelled words in DESCRIPTION:
      Cham (9:552)
      ICCS (9:473)
      Nowakowski (9:306)
      Springer (9:542)
      Szymon (9:299)

# Local checks

## local R CMD check results
There were no ERRORs, WARNINGs, or NOTEs. 

## Downstream dependencies
There seem to be no downstream dependencies:

```{r revdep}
devtools::revdep("DMRnet")
# character(0)
```


