# pseudorank 1.0.1

[![CRANstatus](https://www.r-pkg.org/badges/version/pseudorank)](https://cran.r-project.org/package=pseudorank/)
[![](https://cranlogs.r-pkg.org/badges/pseudorank)](https://cran.r-project.org/package=pseudorank/)
[![Travis-CI Build Status](https://travis-ci.org/happma/pseudorank.svg?branch=master)](https://travis-ci.org/happma/pseudorank/)
[![Build status](https://ci.appveyor.com/api/projects/status/queq8aa7cpct3j16?svg=true)](https://ci.appveyor.com/project/happma/pseudorank/)

This R package provides a function written in C++ to calculate pseudo-ranks in R and some rank statistics which can opionally use pseudo-ranks instead of ranks.
For a definition and discussion of pseudo-ranks, see for example 

Brunner, E., Bathke, A.C., and Konietschke, F. (2018a). Rank- and Pseudo-Rank Procedures for Independent Observations in Factorial Designs - Using R and SAS. Springer Series in Statistics, Springer, Heidelberg. ISBN: 978-3-030-02912-8.

Happ M, Zimmermann G, Brunner E, Bathke AC (2020). Pseudo-Ranks: How to Calculate Them Efficiently in R. Journal of Statistical Software, Code Snippets, *95*(1), 1-22. doi: 10.18637/jss.v095.c01 (URL:https://doi.org/10.18637/jss.v095.c01).


To install the current development version from github:

``` r
## install devtools package
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
# install package
devtools::install_github("happma/pseudorank")
library(pseudorank)
```

## Calculating Pseudo-Ranks

The function 'pseudorank' can either be used with data.frames or with vectors. Please note that when using a data.frame only one grouping factor can be used.

``` r
# some example data
df <- data.frame(data =rnorm(100),group = c(rep(1,40),rep(2,40),rep(3,20)))
df$group <- as.factor(df$group)

# two ways to calculate pseudo-ranks
# Variant 1: use a vector for the data and a group vector
pseudorank(df$data, df$group)

# Variant 2: use a formual object, Note that only one group factor can be used
# that is, in data~group*group2 only 'group' will be used
pseudorank(data ~ group, df)
```
Similarly to the 'rank' function from base R, there are several different methods which can be used to handle ties in the data. Most nonparametric tests rely on so-called mid-(pseudo)-ranks, these are just the average of minimum and maximum (pseudo)-ranks.

``` r
# some example data
df <- data.frame(data =rnorm(100),group = c(rep(1,40),rep(2,40),rep(3,20)))
df$group <- as.factor(df$group)

# min and max pseudo-ranks
m <- pseudorank(df$data ,df$group, ties.method = "min")
M <- pseudorank(df$data ,df$group, ties.method = "max")

mid <- pseudorank(df$data ,df$group, ties.method = "average") # identical to (m+M)/2

```

In case of missing values there are three options to choose from. These are the
same as for the function 'rank' or 'sort' from base R. It is recommended to use
'na.last = NA' to remove the NAs. If the NAs are kept, they can either be put at the beginning
or the end of your data, then the pseudo-ranks from those NAs depend on the order they appear in the data.
The order does not matter only if the groups containing missing values have the same sample size.
See the following R Code for an example of this problem where observation 1 and 4 are interchanged.
Here, the pseudo-ranks for those two observations are different, all other pseudo-ranks remain
unchanged.

``` r
data = c(NA,7,1,NA,3,3,5.5,6,NA, 3, 1)
group =  as.factor(c(1,1,1,2,2,2,3,3,3,1,3))

# Variant 1: keep NAs and put them last
pseudorank(data, group, na.last = TRUE)

# we change the order of observation 1 and 4 (both NAs)
group =  as.factor(c(2,1,1,1,2,2,3,3,3,1,3))
pseudorank(data, group, na.last = TRUE)

# Variant 2: keep NAs and put them first
pseudorank(data, group, na.last = FALSE)

# Variant 3: remove Nas (recommended)
pseudorank(data, group, na.last = NA)

```

## Hettmansperger-Norton Test for Patterned Alternatives in k-Sample Problems

The test implemented in this package uses pseudo-ranks instead of ranks. This is mainly due to paradoxical results caused by ranks. See 

Brunner, E., Konietschke, F., Bathke, A. C., & Pauly, M. (2018). Ranks and Pseudo-Ranks-Paradoxical Results of Rank Tests. arXiv preprint arXiv:1802.05650.

for a discussion of this problem.

``` r
# create some data, please note that the group factor needs to be ordered
df <- data.frame(data = c(rnorm(40, 3, 1), rnorm(40, 2, 1), rnorm(20, 1, 1)),
  group = c(rep(1,40),rep(2,40),rep(3,20)))
df$group <- factor(df$group, ordered = TRUE)

# you can either test for a decreasing, increasing or custom trend
hettmansperger_norton_test(df$data, df$group, alternative="decreasing")
hettmansperger_norton_test(df$data, df$group, alternative="increasing")
hettmansperger_norton_test(df$data, df$group, alternative="custom", trend = c(1, 3, 2))

```

## Kruskal-Wallis Test with Pseudo-Ranks

The Kruskal-Wallis test implemented in this package can use pseudo-ranks, if the argument 'pseudoranks = TRUE' is used.

``` r
# create some artificial data
x = c(1, 1, 1, 1, 2, 3, 4, 5, 6)
grp = as.factor(c('A','A','B','B','B','D','D','D','D'))
df = data.frame(x = x, grp = grp)

kruskal_wallis_test(x~grp, data=df, pseudoranks=TRUE)
```
