# treeducken
<!-- badges: start -->
[![R build status](https://github.com/wadedismukes/rtreeducken/workflows/R-CMD-check/badge.svg)](https://github.com/wadedismukes/treeducken/actions)
<!-- badges: end -->
`treeducken` is a phylogenetic simulation package with the primary purpose of simulating host and symbiont phylogenies simultaneously
using the machinery of birth-death models. This extends work on simulating cophylogenetic systems by allowing for scenarios
where speciation in the host does not always imply speciation in the symbiont (Keller-Schmidt et al. 2011). 
This simulation also accounts for extinction in its simulation of host-symbiont systems. 
`treeducken` also can simulate under the three-tree model in a manner similar to that used in the [SimPhy program](https://github.com/adamallo/SimPhy) (Mallo et al. 2016).

## Brief usage example

Below is a quick example of how to use `treeducken` to simulate a set of ten
host and symbiont tree sets. To do this we set a host birth rate, a host death 
rate, a symbiont birth rate, a symbiont death rate, a cospeciation rate, and a 
host expansion rate. These parameters are further described in the manuscript 
(in prep.) and the vignette (`cophylogenetic_sim.Rmd`). Finally, we have to set 
the number of sets we would like and the time to simulate until. Note that since
this is a time-based birth-death simulation some of the host and symbiont trees
will have 3 or less tips and be fairly uninteresting (in my opinion). 

In R:
```
library(treeducken)
h_lambda <- 0.5
h_mu <- 0.3
c_lambda <- 0.5

s_lambda <- 1.0
s_mu <- 0.3
s_her <- 0.0

host_symb_sets <- sim_cophylo_bdp(hbr = h_lambda,
                hdr = h_mu,
                sbr = s_lambda,
                cosp_rate = c_lambda,
                sdr = s_mu,
                host_exp_rate = s_her,
                timeToSimTo = 2.0,
                numbsim = 10)
```


## Installation 

For now you may install the development version using the `devtools` package.

In R:
```
library(devtools)
install_github("wadedismukes/rtreeducken")
```
