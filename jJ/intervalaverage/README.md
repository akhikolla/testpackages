
<!-- README.md is generated from README.Rmd. Please edit that file -->

# intervalaverage

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/kaufman-lab/intervalaverage.svg?branch=main)](https://travis-ci.com/kaufman-lab/intervalaverage)
<!-- badges: end -->

Perform fast and memory efficient time-weighted averaging of values
measured over intervals into new arbitrary intervals.

This package is useful in the context of data measured or represented as
constant values over intervals on a one-dimensional discrete axis
(e.g. time-integrated averages of a curve over defined periods).

The intervalaverage package was written specifically to deal with air
pollution data recorded or predicted as averages over sampling periods.
Data in this format often needs to be shifted to non-aligned periods or
averaged up to longer durations (e.g. averaging data measured over
sequential non-overlapping periods to calendar years).

This package makes careful use of optimized data.table syntax to reduce
copying so it is well suited to deal with large datasets. The function
inputs and return values are always data.tables (with values, interval
starts and interval ends stored as columns) which provide a
straightforward representation of the data facilitating direct
inspection and manipulation via data.table’s standard interface.
Functions also accept optional grouping variables as arguments; these
group-by operations get passed to data.table which avoids the overhead
of looping directly in R.

## Installation

You can install the package from CRAN (pending):

    install.packages("intervalaverage")

Or you can install the current development version of the from
[GitHub](https://github.com/):

    # install.packages("devtools")
    library(devtools)
    install_github("kaufman-lab/intervalaverage","main",build_vignettes=TRUE)

or you can manually clone this repository from github and install from
source.

## Example

``` r
library(intervalaverage)
#> Loading required package: data.table
```

Consider some PM2.5 data measured over periods occurring over a
repeating 7-day schedule:

``` r
(x <- data.table(start=seq(1L,by=7L,length=6),
                 end=seq(7L,by=7L,length=6),
                 pm25=c(10,12,8,14,22,18)))
#>    start end pm25
#> 1:     1   7   10
#> 2:     8  14   12
#> 3:    15  21    8
#> 4:    22  28   14
#> 5:    29  35   22
#> 6:    36  42   18
```

Note that the start and end columns define closed intervals
(i.e. inclusive).

Imagine that we would like these PM2.5 measurements to be represented
over a different set of intervals, such as the following 7-day schedule
which does not align cleanly with the intervals in x:

``` r
(y <- data.table(start=seq(3L,by=7L,length=6),
                 end=seq(9L,by=7L,length=6)))
#>    start end
#> 1:     3   9
#> 2:    10  16
#> 3:    17  23
#> 4:    24  30
#> 5:    31  37
#> 6:    38  44
```

Representing the PM25 measurements over the exact intervals in y is not
strictly possible, since each y interval contains two x intervals which
only partially overlap.

While there are several ways to compromise and approximate a
representation of x over the periods in y, one reasonable approach is to
average all the values of PM2.5 in x occuring over the interval of y.
But since the two overlap durations of the two x periods into a single y
period are unequal, the PM2.5 values should be averaged in a way that
takes the duration of overlap into account. That is to say, we need a
time-weighted average of values of x into intervals in y.

This is the exact purpose of the intervalaverage function:

``` r
(z <- intervalaverage(x,y,interval_vars=c("start","end"),
                value_vars=c("pm25")))
#>    start end      pm25 yduration xduration nobs_pm25 xminstart xmaxend
#> 1:     3   9 10.571429         7         7         7         3       9
#> 2:    10  16 10.857143         7         7         7        10      16
#> 3:    17  23  9.714286         7         7         7        17      23
#> 4:    24  30 16.285714         7         7         7        24      30
#> 5:    31  37 20.857143         7         7         7        31      37
#> 6:    38  44        NA         7         5         5        38      42
```

(Note that the interval average package always assumes specified
intervals are closed aka inclusive)

The pm25 value of 10.571429 is equal to `(5/7)*10 + (2/7)*12`, where
`5/7` and `2/7` are the respective time-weights of the the PM2.5 values
based on their respective durations of overlap with the period `[3,9]`.

``` r
identical(z[1,pm25], (5/7)*10 + (2/7)*12)
#> [1] TRUE
```

Note that in z the average value of PM2.5 is missing (NA) for the final
period in y. This is because there is insufficient PM25 to calculate an
average. By using a looser completeness requirements we can get an
average for this period in y (in general, this will just be the mean of
non-missing values in x):

``` r
(z <- intervalaverage(x,y,interval_vars=c("start","end"),
                value_vars=c("pm25"),required_percentage = 70))
#>    start end      pm25 yduration xduration nobs_pm25 xminstart xmaxend
#> 1:     3   9 10.571429         7         7         7         3       9
#> 2:    10  16 10.857143         7         7         7        10      16
#> 3:    17  23  9.714286         7         7         7        17      23
#> 4:    24  30 16.285714         7         7         7        24      30
#> 5:    31  37 20.857143         7         7         7        31      37
#> 6:    38  44 18.000000         7         5         5        38      42
```

There is much more described in the vignette. For example, it is
possible to calculate averages within a grouping variable. For example
if PM2.5 were measured at multiple locations, you could take the average
separately at each location in one call to the `intervalaverage`
function. Averaging intervals need not be the same for each locations.
See the vignette for examples.

## Vignette

Read the intro vignette for an in-depth demonstration of the package
functions.

    library(intervalaverage)
    vignette("intervalaverage-intro")

For an overview of the internals of the `intervalaverage` function, see
the technical vignette:

    vignette("intervalaverage-technicaloverview")
