
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ced

<!-- badges: start -->

[![GitLab CI Build
Status](https://gitlab.com/artemklevtsov/ced/badges/master/pipeline.svg)](https://gitlab.com/artemklevtsov/ced/pipelines)
[![AppVeyor Build
status](https://ci.appveyor.com/api/projects/status/p3513kcbhh8rp4hk?svg=true)](https://ci.appveyor.com/project/artemklevtsov/ced)
[![Codecov Code
Coverage](https://codecov.io/gl/artemklevtsov/ced/branch/master/graph/badge.svg)](https://codecov.io/gl/artemklevtsov/ced)
[![CRAN
Status](http://www.r-pkg.org/badges/version/ced)](https://cran.r-project.org/package=ced)
[![License: GPL
v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

<!-- badges: end -->

R bindings for the [`Google Compact Encoding
Detection`](https://github.com/google/compact_enc_det) library. Key
features:

  - process character vector;
  - process raw vector;

## Installation

To install the package from the CRAN run the following command:

``` r
install.packages("ced", repos = "https://cloud.r-project.org/")
```

Also you could install the dev-version with the `install_gitlab()`
function from the `remotes` package:

``` r
remotes::install_gitlab("artemklevtsov/ced")
```

This package contains the compiled code, therefore you have to use the
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) to install it
on Windows.

## Example

``` r
# load packages
library(ced)

# detect string encoding
ascii <- "Hello, useR!"
print(ascii)
#> [1] "Hello, useR!"
ced_enc_detect(ascii)
#> [1] "US-ASCII"
utf8 <- "\u4e0b\u5348\u597d"
print(utf8)
#> [1] "下午好"
ced_enc_detect(utf8)
#> [1] "UTF-8"

# detect raw vector encoding
ced_enc_detect(charToRaw(ascii))
#> [1] "US-ASCII"
ced_enc_detect(charToRaw(utf8))
#> [1] "UTF-8"
```

## Bug reports

Use the following command to go to the page for bug report submissions:

``` r
bug.report(package = "ced")
```

Before reporting a bug or submitting an issue, please do the following:

  - Make sure that you error or issue was not reported or discussed
    earlier. Please, use the search;
  - Check the news list of the current version. Some errors could be
    caused by the package changes. It could be done with `news(package =
    "ced", Version == packageVersion("ced"))` command;
  - Make a minimal reproducible example of the code that consistently
    causes the error;
  - Make sure that the error occurs during the execution of a function
    from the `ced` package, not from other packages;
  - Try to reproduce the error with the last development version of the
    package from the git repository.

Please attach traceback() and sessionInfo() output to bug report. It may
save a lot of time.

## License

The `ced` package is distributed under
[GPLv2](http://www.gnu.org/licenses/gpl-2.0.html) license.
