
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RcppUUID

<!-- badges: start -->

[![GitLab CI Build
Status](https://gitlab.com/artemklevtsov/rcppuuid/badges/master/pipeline.svg)](https://gitlab.com/artemklevtsov/rcppuuid/pipelines)
[![AppVeyor Build
status](https://ci.appveyor.com/api/projects/status/if9qot73i61ts59y?svg=true)](https://ci.appveyor.com/project/artemklevtsov/rcppuuid)
[![Codecov Code
Coverage](https://codecov.io/gl/artemklevtsov/rcppuuid/branch/master/graph/badge.svg)](https://codecov.io/gl/artemklevtsov/rcppuuid)
[![CRAN
Status](http://www.r-pkg.org/badges/version/RcppUUID)](https://cran.r-project.org/package=RcppUUID)
[![License: GPL
v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

<!-- badges: end -->

R package to generate Universally Unique Identifiers (UUIDs) version 4
and 5 using Boost C++ library.

## Installation

To install the package from the CRAN run the following command:

``` r
install.packages("RcppUUID", repos = "https://cloud.r-project.org/")
```

Also you could install the dev-version with the `install_gitlab()`
function from the `remotes` package:

``` r
remotes::install_gitlab("artemklevtsov/rcppuuid")
```

This package contains the compiled code, therefore you have to use the
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) to install it
on Windows.

## Usage

### Generate version 4 UUIDs

Generate single UUID:

``` r
RcppUUID::uuid_generate_random()
#> [1] "41b45c2e-9745-4362-8277-95fe998f2ef9"
```

Generate multiple UUIDs:

``` r
RcppUUID::uuid_generate_random(5)
#> [1] "2f669c10-2ad7-461b-a531-103bafed0b15" "fbe4959e-c723-4670-a155-f440c753c509" "074c5472-b314-4e4f-9aab-88bfc72f030a"
#> [4] "ad9ebe8c-4b85-4819-b1bd-6091a2305024" "d5287bae-666c-448b-bcb9-2e41fc173540"
```

Check uniques for the uuids:

``` r
unique_n <- function(x) length(unique(x))
n <- 1000000
unique_n(RcppUUID::uuid_generate_random(n)) == n
#> [1] TRUE
```

Benchmarking:

Single UUID:

``` r
microbenchmark::microbenchmark(
  uuid = uuid::UUIDgenerate(FALSE),
  RcppUUID = RcppUUID::uuid_generate_random()
)
#> Unit: microseconds
#>      expr    min      lq     mean  median      uq      max neval
#>      uuid 15.532 16.0670 44.47824 16.2875 16.5985 2794.521   100
#>  RcppUUID  8.166  8.6705  9.36236  9.4465  9.7080   19.991   100
```

Multiple UUIDs:

``` r
n <- 10000
microbenchmark::microbenchmark(
  uuid = uuid::UUIDgenerate(FALSE, n),
  RcppUUID = RcppUUID::uuid_generate_random(n)
)
#> Unit: milliseconds
#>      expr       min        lq      mean    median        uq       max neval
#>      uuid 74.147780 77.636038 87.449953 79.408061 81.066626 304.77218   100
#>  RcppUUID  7.956943  8.725376  9.133505  9.101491  9.452271  11.24261   100
```

### Generate version 5 UUIDs

Generate version UUIDs based on the text input:

``` r
RcppUUID::uuid_generate_name(letters[1:5])
#> [1] "54a0a790-c611-5b5b-b50e-ff01490ecdfa" "d5080e36-1ba4-5cb3-861c-34b25868f7db" "33ed51b6-a330-5830-bda9-2bac09e15753"
#> [4] "b74b2afe-06d5-5fea-99cc-a7de0b492704" "8535136c-b0d3-5373-aa79-ab67d33a2a8e"
```

For the each unique input will be generated unique UUID. Check
uniqueness:

``` r
uuids <- replicate(10, RcppUUID::uuid_generate_name(letters))
length(unique(as.vector(uuids))) == length(letters)
#> [1] TRUE
```

### Validate UUIDs

``` r
RcppUUID::uuid_validate(NA_character_)
#> [1] FALSE
RcppUUID::uuid_validate("")
#> [1] FALSE
RcppUUID::uuid_validate("not uuid")
#> [1] FALSE
RcppUUID::uuid_validate(RcppUUID::uuid_generate_random(5))
#> [1] TRUE TRUE TRUE TRUE TRUE
RcppUUID::uuid_validate(RcppUUID::uuid_generate_nil(5))
#> [1] TRUE TRUE TRUE TRUE TRUE
RcppUUID::uuid_validate(RcppUUID::uuid_generate_name(letters[1:5]))
#> [1] TRUE TRUE TRUE TRUE TRUE
```

## Bug reports

Use the following command to go to the page for bug report submissions:

``` r
bug.report(package = "RcppUUID")
```

Before reporting a bug or submitting an issue, please do the following:

  - Make sure that you error or issue was not reported or discussed
    earlier. Please, use the search;
  - Check the news list of the current version. Some errors could be
    caused by the package changes. It could be done with `news(package =
    "RcppUUID", Version == packageVersion("RcppUUID"))` command;
  - Make a minimal reproducible example of the code that consistently
    causes the error;
  - Make sure that the error occurs during the execution of a function
    from the `RcppUUID` package, not from other packages;
  - Try to reproduce the error with the last development version of the
    package from the git repository.

Please attach traceback() and sessionInfo() output to bug report. It may
save a lot of time.

## License

The `RcppUUID` package is distributed under
[GPLv2](http://www.gnu.org/licenses/gpl-2.0.html) license.
