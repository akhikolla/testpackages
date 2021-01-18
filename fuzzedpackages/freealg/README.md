The Free Algebra in R
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![Build
Status](https://travis-ci.org/RobinHankin/freealg.svg?branch=master)](https://travis-ci.org/RobinHankin/freealg)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/freealg)](https://cran.r-project.org/package=freealg)
[![Rdoc](http://www.rdocumentation.org/badges/version/freealg)](http://www.rdocumentation.org/packages/freealg)
[![Codecov test
coverage](https://codecov.io/gh/RobinHankin/freealg/branch/master/graph/badge.svg)](https://codecov.io/gh/RobinHankin/freealg/branch/master)
<!-- badges: end -->

# Overview

The free algebra is an interesting and useful object. Here I present the
`freealg` package which provides some functionality for free algebra.

The package uses `C++`’s STL `map` class for efficiency, which has the
downside that the order of the terms is undefined. This does not matter
as the mathematical value is unaffected by reordering; and the print
method does a good job in producing human-readable output.

# Installation

You can install the released version of `freealg` from
[CRAN](https://CRAN.R-project.org) with:

``` r
# install.packages("freealg")  # uncomment this to install the package
library("freealg")
```

# The free algebra

The free algebra is the free R-module with a basis consisting of all
words over an alphabet of symbols with multiplication of words defined
as concatenation. Thus, with an alphabet of
![\\{x,y,z\\}](https://latex.codecogs.com/png.latex?%5C%7Bx%2Cy%2Cz%5C%7D
"\\{x,y,z\\}") and

  
![
A=\\alpha x^2yx + \\beta zy
](https://latex.codecogs.com/png.latex?%0AA%3D%5Calpha%20x%5E2yx%20%2B%20%5Cbeta%20zy%0A
"
A=\\alpha x^2yx + \\beta zy
")  

and

  
![
B=\\gamma z + \\delta y^4
](https://latex.codecogs.com/png.latex?%0AB%3D%5Cgamma%20z%20%2B%20%5Cdelta%20y%5E4%0A
"
B=\\gamma z + \\delta y^4
")  

we would have

  
![
A\\cdot B=\\left(\\alpha x^2yx+\\beta zy\\right)\\cdot\\left(\\gamma
z+\\delta y^4\\right)=\\alpha\\gamma x^2yxz+\\alpha\\delta
x^2yxy^4+\\beta\\gamma zyz+\\beta\\delta zy^5
](https://latex.codecogs.com/png.latex?%0AA%5Ccdot%20B%3D%5Cleft%28%5Calpha%20x%5E2yx%2B%5Cbeta%20zy%5Cright%29%5Ccdot%5Cleft%28%5Cgamma%20z%2B%5Cdelta%20y%5E4%5Cright%29%3D%5Calpha%5Cgamma%20x%5E2yxz%2B%5Calpha%5Cdelta%20x%5E2yxy%5E4%2B%5Cbeta%5Cgamma%20zyz%2B%5Cbeta%5Cdelta%20zy%5E5%0A
"
A\\cdot B=\\left(\\alpha x^2yx+\\beta zy\\right)\\cdot\\left(\\gamma z+\\delta y^4\\right)=\\alpha\\gamma x^2yxz+\\alpha\\delta x^2yxy^4+\\beta\\gamma zyz+\\beta\\delta zy^5
")  

Note that multiplication is not commutative, but it is associative. A
natural and easily implemented extension is to use upper-case symbols to
represent multiplicative inverses of the lower-case equivalents. Thus if

  
![
C=\\epsilon
X^2](https://latex.codecogs.com/png.latex?%0AC%3D%5Cepsilon%20X%5E2 "
C=\\epsilon X^2")  

we would have

  
![
A\\cdot C=\\left(\\alpha x^2yx+\\beta zy\\right)\\cdot\\epsilon X^2=
\\alpha\\epsilon x^2yX + \\beta\\epsilon zyX^2
](https://latex.codecogs.com/png.latex?%0AA%5Ccdot%20C%3D%5Cleft%28%5Calpha%20x%5E2yx%2B%5Cbeta%20zy%5Cright%29%5Ccdot%5Cepsilon%20X%5E2%3D%0A%5Calpha%5Cepsilon%20x%5E2yX%20%2B%20%5Cbeta%5Cepsilon%20zyX%5E2%0A
"
A\\cdot C=\\left(\\alpha x^2yx+\\beta zy\\right)\\cdot\\epsilon X^2=
\\alpha\\epsilon x^2yX + \\beta\\epsilon zyX^2
")  

and

  
![
C\\cdot A=\\epsilon X^2\\cdot\\left(\\alpha x^2yx+\\beta zy\\right)=
\\alpha\\epsilon yx + \\beta\\epsilon X^2zy.
](https://latex.codecogs.com/png.latex?%0AC%5Ccdot%20A%3D%5Cepsilon%20X%5E2%5Ccdot%5Cleft%28%5Calpha%20x%5E2yx%2B%5Cbeta%20zy%5Cright%29%3D%0A%5Calpha%5Cepsilon%20yx%20%2B%20%5Cbeta%5Cepsilon%20X%5E2zy.%0A
"
C\\cdot A=\\epsilon X^2\\cdot\\left(\\alpha x^2yx+\\beta zy\\right)=
\\alpha\\epsilon yx + \\beta\\epsilon X^2zy.
")  

The system inherits power associativity from distributivity and
associativity of concatenation, but is not commutative.

# The `freealg` package in use

Creating a free algebra object is straightforward. We can coerce from a
character string with natural idiom:

``` r
X <- as.freealg("1 + 3a + 5b + 5abba")
X
#> free algebra element algebraically equal to
#>  + 1 + 3*a + 5*abba + 5*b
```

or use a more formal method:

``` r
freealg(sapply(1:5,seq_len),1:5)
#> free algebra element algebraically equal to
#>  + 1*a + 2*ab + 3*abc + 4*abcd + 5*abcde
```

``` r
Y <- as.freealg("6 - 4a +2aaab")
X+Y
#> free algebra element algebraically equal to
#>  + 7 - 1*a + 2*aaab + 5*abba + 5*b
X*Y
#> free algebra element algebraically equal to
#>  + 6 + 14*a - 12*aa + 6*aaaab + 2*aaab + 30*abba - 20*abbaa + 10*abbaaaab + 30*b - 20*ba + 10*baaab
X^2
#> free algebra element algebraically equal to
#>  + 1 + 6*a + 9*aa + 15*aabba + 15*ab + 10*abba + 15*abbaa + 25*abbaabba + 25*abbab + 10*b + 15*ba + 25*babba + 25*bb
```

We can demonstrate associativity (which is non-trivial):

``` r
set.seed(0)
(x1 <- rfalg(inc=TRUE))
#> free algebra element algebraically equal to
#>  + 7*C + 6*Ca + 4*B + 3*BC + 1*a + 5*aCBB + 2*bc
(x2 <- rfalg(inc=TRUE))
#> free algebra element algebraically equal to
#>  + 6 + 1*CAAA + 2*Ca + 3*Cbcb + 7*aaCA + 4*b + 5*c
(x3 <- rfalg(inc=TRUE))
#> free algebra element algebraically equal to
#>  + 3*C + 5*CbAc + 1*BACB + 2*a + 10*b + 7*cb
```

(function `rfalg()` generates random `freealg` objects). Then

``` r
x1*(x2*x3) == (x1*x2)*x3
#> [1] TRUE
```

# Further information

For more detail, see the package vignette

`vignette("freealg")`
