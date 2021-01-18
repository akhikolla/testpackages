---
title: "NEWS about the PCMBaseCpp R-package"
author: "Venelin Mitov"
date: "23 March, 2020"
output: html_document
---

# PCMBaseCpp 0.1.9
Fixing an R-process crash on CRAN fedora-clang.

# PCMBaseCpp 0.1.8
* Minor bug and documentation fixes.  

# PCMBaseCpp 0.1.7

* Disabled OPENMP compilation in an attempt to fix the failures on "Fedora clang" and "Solaris". Note that these errors could not be reproduced using devtools::check_rhub().
* Generated project web-page: https://venelin.github.io/PCMBaseCpp .

# PCMBaseCpp 0.1.6

* Fixed an crash situation on Windows due to exceptions bypassing the 
`try(..., silent=TRUE)` construct. See
https://github.com/venelin/PCMBaseCpp/issues/1#issue-507813590 and 
https://github.com/RcppCore/Rcpp/issues/1001#issue-509047134 for details.
* Added a `NEWS.md` file to track changes to the package.
