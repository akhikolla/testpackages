# bama v1.0.1

## Bug Fix
* `lm1` in the C++ was not being initialized. This has been fixed.

# bama v1.0.0

## New Features
* `bama` now contains `fdr.bama`, a `bama` variant that controls for the false
    discovery rate by estimating the null PIP distribution via the permutation
    test. It contains options to run in parallel (ie using the `parallel` package)
    as well, since the permutation test can be time consuming. `fdr.bama` also
    contains a summary function, `summary.fdr.bama` to help summarize results.

## Major Changes
* `bama` (and `fdr.bama`) now contain a weights parameter, `weights`.
* Change k, lm0, lm1, l from hard coded parameters in the C++ code to user
    set parameters (with defaults).
* Revert change that combined C1 and C2.

## Minor Changes
* Partial revamp `bama` documentation.
* Edit DESCRIPTION, and README / vignette introduction.
* Update DOI with (online) version of record.
* Redo `bama` type checking.

# bama v0.9.1

## Major Changes

* The most notable feature in this release is the renaming of `hdbm` to `bama`
    to prevent confusion with regards to its relationship with `hdmm`.
* bama now only accepts a single extra covariate matrix C, instead of C1 and C2.
## Minor Changes

* Added a summary function to `bama` to calculate the PIP, estimate, and
    credible interval of each mediator and optionally rank them.
* Added data quality check to ensure the matrix column norms in M, C1 and C2
    are not too small. If they are too small, `NaNs`
could be generated because of overflow problems when dividing by the norm.
* Improved documentation
* Added checks to C, the matrix of extra covariates, and alpha.a, the initial
alpha.a values.

# hdbm v0.9.0

* Initial version released on CRAN.
