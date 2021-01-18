## RcppGetconf [![Build Status](https://travis-ci.org/eddelbuettel/rcppgetconf.svg)](https://travis-ci.org/eddelbuettel/rcppgetconf) [![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html) [![CRAN](http://www.r-pkg.org/badges/version/RcppGetconf)](https://cran.r-project.org/package=RcppGetconf) [![Dependencies](https://tinyverse.netlify.com/badge/RcppGetconf)](https://cran.r-project.org/package=RcppGetconf) [![Downloads](https://cranlogs.r-pkg.org/badges/RcppGetconf?color=brightgreen)](https://www.r-pkg.org/pkg/RcppGetconf)

Rcpp Read Access to System Configuration Settings

### What is this?

Modern POSIX systems have a binary `getconf` which can access the system
calls `sysconf`, `pathconf` and `confstr`.  This package brings the
values back to R.

### Requirements

This package requires access to these system calls, and definitions of its
data structures in the system header files.  It works on Linux and OS X and
_should_ work on other POSIX-compliant OSs. Contributions would be very
welcome.

### Installation

The package is on [CRAN](https://cran.r-project.org) and can be installed via
a standard

```r
R> install.packages("RcppGetconf")
```

command.

### Status

It contains two useful functions right now.  It currently builds cleanly on
Linux and OS X; reports from other builds would (and PRs where needed) would
be greatly appreciated.

### Author

Dirk Eddelbuettel

### License

GPL (>= 2)
