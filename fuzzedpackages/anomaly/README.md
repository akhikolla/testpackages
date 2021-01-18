# anomaly
Fast anomaly detection in R

## In Brief
This R package implements **CAPA** (**C**ollective **A**nd **P**oint **A**nomalies) introduced by [Fisch, Eckley and Fearnhead (2018)](https://arxiv.org/abs/1806.01947). The package is available on [CRAN](https://CRAN.R-project.org/package=anomaly) and contains lightcurve data from the Kepler telescope to illustrate the algorithm.

## About CAPA
CAPA detects and distinguishes between collective and point anomalies. The algorithm's runtime scales linearly at best and quadratically at worst in the number of datapoints. It is coded in C and can process 10000 datapoints almost instantly.  
