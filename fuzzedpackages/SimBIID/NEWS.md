# v0.2.0

Added bisection method for choosing tolerances when multiple
data points are being used in the ABC-SMC algorithm.

Added option to set minimum tolerances for ABC-SMC.

Changed way that random seeds are passed to `mclapply()` to
aid reproducibility when setting seeds and using parallelisation.
Runs will only be reproducible if using the same number of cores
each time (which can be specified using `mc.cores` argument to
various functions).

# v0.1.4

Patch release to fix minor bug in the `predict.PMCMC` method.

# v0.1.3

Patch release to fix the order in which `tspan` objects are
updated in code produced from `mparseRcpp`.

# v0.1.2

Patch release to fix a parsing error when compiling models
with only a single transition term.

# v0.1.1

Patch release to fix a compilation error on Solaris, and also
to fix a memory leak error.

# SimBIID 0.1.0

First release of SimBIID to CRAN.
