#!/bin/bash
max=8
for beta in 1 2 3 4 5 6; do
  for snr in 1 2 3 4 5 6 7; do
    while [ `ps axu | grep R/bin/exec/R |wc -l` -gt $max ] ; do
      sleep 60
      echo `ps axu | grep R/bin/exec/R |wc -l`
    done
    echo $beta $snr
    Rscript massive_simulations.R  $beta $snr &> ms$beta$snr &
    sleep 60   #to give the Rscript time to get invoked
  done
done
