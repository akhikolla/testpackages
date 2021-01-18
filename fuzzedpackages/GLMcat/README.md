
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GLMcat

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ylleonv/GLMcat")
```

We introduce and illustrate the utility of GLMcat, the R package we
developed to estimate generalized linear models implemented under the
unified specification \((r, F, Z)\), where \(r\) represents the ratio of
probabilities (reference, cumulative, adjacent, or sequential), \(F\)
the cumulative distribution function for the linkage, and \(Z\) the
design matrix. We present the properties of the four families of models,
which must be investigated when selecting the components \(r\), \(F\),
and \(Z\). The functions are user-friendly and fairly intuitive;
offering the possibility to choose from a large range of models through
a combination \((r, F, Z)\).
