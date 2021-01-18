**m2r** – Macaulay2 in R
========================

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/m2r)](https://cran.r-project.org/package=m2r)
[![Travis build
status](https://travis-ci.org/dkahle/m2r.svg?branch=master)](https://travis-ci.org/dkahle/m2r)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/dkahle/m2r?branch=master&svg=true)](https://ci.appveyor.com/project/dkahle/m2r)

Overview
--------

**m2r** is a new R package that provides a persistent connection between
[R](https://www.r-project.org) and [Macaulay2
(M2)](http://www.math.uiuc.edu/Macaulay2/).

The package grew out of a collaboration at the [2016 Mathematics
Research
Community](http://www.ams.org/programs/research-communities/mrc-16) on
algebraic statistics, funded by the [National Science
Foundation](http://www.nsf.gov) through the [American Mathematical
Society](http://www.ams.org/home/page).

If you have a feature request, please file an issue!

Getting started
---------------

**m2r** is loaded like any other R package:

``` r
library(m2r)
# Loading required package: mpoly
#   Please cite m2r! See citation("m2r") for details.
#   M2 found in /Applications/Macaulay2-1.10/bin
```

When loaded, **m2r** initializes a persistent connection to a back-end
Macaulay2 session. The basic function in R that accesses this connection
is `m2()`, which simply accepts a character string that is run by the
Macaulay2 session.

``` r
m2("1 + 1")
# Starting M2... done.
# [1] "2"
```

You can see the persistence by setting variables and accessing them
across different `m2()` calls:

``` r
m2("a = 1")
# [1] "1"
m2("a")
# [1] "1"
```

You can check the variables defined in the M2 session with `m2_ls()`:

``` r
m2_ls()
# [1] "a"
```

You can also check if variables exist with `m2_exists()`:

``` r
m2_exists("a")
# [1] TRUE
m2_exists(c("a","b"))
# [1]  TRUE FALSE
```

Apart from the basic connection to M2, **m2r** has basic data structures
and methods to reference and manipulate the M2 objects within R. For
more on this, see the **m2r** internals section below.

Rings, ideals, and Grobner bases
--------------------------------

**m2r** currently has basic support for
[rings](https://en.wikipedia.org/wiki/Ring_(mathematics)) (think:
[polynomial rings](https://en.wikipedia.org/wiki/Polynomial_ring)):

``` r
(QQtxyz <- ring("t", "x", "y", "z", coefring = "QQ"))
# M2 Ring: QQ[t,x,y,z], grevlex order
```

and [ideals](https://en.wikipedia.org/wiki/Ideal_(ring_theory)) of
rings:

``` r
(I <- ideal("t^4 - x", "t^3 - y", "t^2 - z"))
# M2 Ideal of ring QQ[t,x,y,z] (grevlex) with generators : 
# < t^4  -  x,  t^3  -  y,  t^2  -  z >
```

You can compute [Grobner
bases](https://en.wikipedia.org/wiki/Gröbner_basis) as well. The basic
function to do this is `gb()`:

``` r
gb(I)
# z^2  -  x
# z t  -  y
# -1 z x  +  y^2
# -1 x  +  t y
# -1 z y  +  x t
# -1 z  +  t^2
```

Perhaps an easier way to do this is just to list off the polynomials as
character strings:

``` r
gb("t^4 - x", "t^3 - y", "t^2 - z")
# z^2  -  x
# z t  -  y
# -1 z x  +  y^2
# -1 x  +  t y
# -1 z y  +  x t
# -1 z  +  t^2
```

The result is an `mpolyList` object, from the [**mpoly**
package](https://github.com/dkahle/mpoly). You can see the M2 code by
adding `code = TRUE`:

``` r
gb("t^4 - x", "t^3 - y", "t^2 - z", code = TRUE)
# m2rintgb00000003 = gb(m2rintideal00000003); gens m2rintgb00000003
```

You can compute the basis respective of different [monomial
orders](https://en.wikipedia.org/wiki/Monomial_order) as well. The
default ordering is the one in the respective ring, which defaults to
`grevlex`; however, changing the order is as simple as changing the
ring.

``` r
ring("x", "y", "t", "z", coefring = "QQ", order = "lex")
# M2 Ring: QQ[x,y,t,z], lex order
gb("t^4 - x", "t^3 - y", "t^2 - z")
# t^2  -  z
# -1 t z  +  y
# -1 z^2  +  x
```

On a technical level, `ring()`, `ideal()`, and `gb()` use [nonstandard
evaluation
rules](http://adv-r.had.co.nz/Computing-on-the-language.html). A more
stable way to use these functions is to use their standard evaluation
versions `ring_()`, `ideal_()`, and `gb_()`. Each accepts first a data
structure describing the relevant object of interest first as its own
object. For example, at a basic level this simply changes the previous
syntax to

``` r
use_ring(QQtxyz)
poly_chars <- c("t^4 - x", "t^3 - y", "t^2 - z")
gb_(poly_chars)
# z^2  -  x
# z t  -  y
# -1 z x  +  y^2
# -1 x  +  t y
# -1 z y  +  x t
# -1 z  +  t^2
```

`gb_()` is significantly easier to code with than `gb()` in the sense
that its inputs and outputs are more predictable, so we strongly
recommend that you use `gb_()`, especially inside of other functions and
packages.

As far as other kinds of computations are concerned, we present a
potpurri of examples below.

Ideal saturation:

``` r
ring("x", coefring = "QQ")
# M2 Ring: QQ[x], grevlex order
I <- ideal("(x-1) x (x+1)")
saturate(I, "x") # = (x-1) (x+1)
# M2 Ideal of ring QQ[x] (grevlex) with generator : 
# < x^2  -  1 >
```

Radicalization:

``` r
I <- ideal("x^2")
radical(I)
# M2 Ideal of ring QQ[x] (grevlex) with generator : 
# < x >
```

Primary decomposition:

``` r
ring("x", "y", "z", coefring = "QQ")
# M2 Ring: QQ[x,y,z], grevlex order
I <- ideal("x z", "y z")
primary_decomposition(I)
# M2 List of ideals of QQ[x,y,z] (grevlex) : 
# < z >
# < x,  y >
```

Dimension:

``` r
ring("x", "y", coefring = "QQ")
# M2 Ring: QQ[x,y], grevlex order
I <- ideal("y - (x+1)") 
dimension(I)
# [1] 1
```

Factoring integers and polynomials
----------------------------------

You can compute [prime
decompositions](https://en.wikipedia.org/wiki/Integer_factorization) of
integers with `factor_n()`:

``` r
(x <- 2^5 * 3^4 * 5^3 * 7^2 * 11^1)
# [1] 174636000
factor_n(x)
# $prime
# [1]  2  3  5  7 11
# 
# $power
# [1] 5 4 3 2 1
```

You can also [factor
polynomials](https://en.wikipedia.org/wiki/Factorization) over rings
using `factor_poly()`:

``` r
factor_poly("x^4 - y^4")
# $factor
# x  -  y
# x  +  y
# x^2  +  y^2
# 
# $power
# [1] 1 1 1
```

Smith normal form of a matrix
-----------------------------

The Smith normal form of a matrix *M* here refers to the decomposition
of an integer matrix *D = PMQ*, where *D*, *P*, and *Q* are integer
matrices and *D* is diagonal. *P* and *Q* are unimodular matrices (their
determinants are -1 or 1), so they are invertible. This is somewhat like
a singular value decomposition for integer matrices.

``` r
M <- matrix(c(
   2,  4,   4,
  -6,  6,  12,
  10, -4, -16
), nrow = 3, byrow = TRUE)

(mats <- snf(M))
# $D
#      [,1] [,2] [,3]
# [1,]   12    0    0
# [2,]    0    6    0
# [3,]    0    0    2
# M2 Matrix over ZZ[]
# $P
#      [,1] [,2] [,3]
# [1,]    1    0    1
# [2,]    0    1    0
# [3,]    0    0    1
# M2 Matrix over ZZ[]
# $Q
#      [,1] [,2] [,3]
# [1,]    4   -2   -1
# [2,]   -2    3    1
# [3,]    3   -2   -1
# M2 Matrix over ZZ[]
P <- mats$P; D <- mats$D; Q <- mats$Q

P %*% M %*% Q                # = D
#      [,1] [,2] [,3]
# [1,]   12    0    0
# [2,]    0    6    0
# [3,]    0    0    2
solve(P) %*% D %*% solve(Q)  # = M
#      [,1] [,2] [,3]
# [1,]    2    4    4
# [2,]   -6    6   12
# [3,]   10   -4  -16

det(P)
# [1] 1
det(Q)
# [1] -1
```

**m2r** internals: pointers, reference and value functions, and `m2` objects
----------------------------------------------------------------------------

At a basic level, **m2r** works by passing strings between R and M2.
Originating at the R side, these strings are properly formated M2 code
constructed from the inputs to the R functions. That code goes to M2, is
evaluated there, and then “exported” with M2’s function
`toExternalString()`. The resulting string often, but not always,
produces the M2 code needed to recreate the object resulting from the
evaluation, and in that sense is M2’s version of R’s `dput()`. That
string is passed back into R and parsed there into R-style data
structures, typically [S3-classed
lists](http://adv-r.had.co.nz/OO-essentials.html#s3).

The R-side parsing of the external string from M2 is an expensive
process because it is currently implemented in R. Consequently (and for
other reasons, too!), in some cases you’ll want to do a M2 computation
from R, but leave the output in M2. Since you will ultimately want
something in R referencing the result, nearly every **m2r** function
that performs M2 computations has a pointer version. As a simple naming
convention, the name of the function that returns the pointer, called
the reference function, is determined by the name of the ordinary
function, called the value function, by appending a `.`.

For example, we’ve seen that `factor_n()` computes the prime
decomposition of a number. The corresponding reference function is
`factor_n.()`:

``` r
(x <- 2^5 * 3^4 * 5^3 * 7^2 * 11^1)
# [1] 174636000
factor_n.(x)
# M2 Pointer Object
#   ExternalString : new Product from {new Power from {2,5},new Power fro...
#          M2 Name : m2o460
#         M2 Class : Product (WrapperType)
```

All value functions simply wrap reference functions and parse the output
with `m2_parse()`, a general M2 parser, often followed by a little more
parsing. `m2_parse()` typically creates an object of class `m2` so that
R knows what kind of thing it is. For example:

``` r
class(factor_n.(x))
# [1] "m2_pointer" "m2"
```

Even more, `m2_parse()` often creates objects that have an inheritance
structure that references `m2` somewhere in the middle of its class
structure, with specific structure preceding and general structure
succeeding (examples below). Apart from its class, the general principle
we follow here for the object itself is this: if the M2 object has a
direct analogue in R, it is parsed into that kind of R object and
additional M2 properties are kept as metadata (attributes); if there is
no direct analogue in R, the object is an `NA` with metadata.

Perhaps the easiest way to see this is with a matrix. `m2_matrix()`
creates a matrix on the M2 side from input on the R side. In the
following, to make things more clear we use [**magrittr**’s pipe
operator](https://github.com/tidyverse/magrittr), with which the
following calls are semantically equivalent: `g(f(x))` and
`x %>% f %>% g`.

``` r
library(magrittr)
mat <- matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2)
mat %>% m2_matrix.   # = m2_matrix.(mat)
# M2 Pointer Object
#   ExternalString : map((ZZ)^3,(ZZ)^2,{{1, 4}, {2, 5}, {3, 6}})
#          M2 Name : m2rintmatrix00000001
#         M2 Class : Matrix (Type)
mat %>% m2_matrix. %>% m2_parse
#      [,1] [,2]
# [1,]    1    4
# [2,]    2    5
# [3,]    3    6
# M2 Matrix over ZZ[]
mat %>% m2_matrix. %>% m2_parse %>% str
#  'm2_matrix' int [1:3, 1:2] 1 2 3 4 5 6
#  - attr(*, "m2_name")= chr "m2rintmatrix00000003"
#  - attr(*, "m2_meta")=List of 1
#   ..$ ring: 'm2_polynomialring' logi NA
#   .. ..- attr(*, "m2_name")= chr "ZZ"
#   .. ..- attr(*, "m2_meta")=List of 3
#   .. .. ..$ vars    : NULL
#   .. .. ..$ coefring: chr "ZZ"
#   .. .. ..$ order   : chr "grevlex"
mat %>% m2_matrix    # = m2_parse(m2_matrix.(mat))
#      [,1] [,2]
# [1,]    1    4
# [2,]    2    5
# [3,]    3    6
# M2 Matrix over ZZ[]
```

It may be helpful to think of every `m2` object as being a missing value
(`NA`, a `logical(1)`) with two M2 attributes: their name (`m2_name`)
and a capture-all named list (`m2_meta`). These can be accessed with
`m2_name()` and `m2_meta()`. For example, a ring, having no analogous
object in R, is an `NA` with attributes:

``` r
r <- ring("x", "y", coefring = "QQ")
str(r)
#  'm2_polynomialring' logi NA
#  - attr(*, "m2_name")= chr "m2rintring00000006"
#  - attr(*, "m2_meta")=List of 3
#   ..$ vars    :List of 2
#   .. ..$ : chr "x"
#   .. ..$ : chr "y"
#   ..$ coefring: chr "QQ"
#   ..$ order   : chr "grevlex"
class(r)
# [1] "m2_polynomialring" "m2"
m2_name(r)
# [1] "m2rintring00000006"
m2_meta(r)
# $vars
# $vars[[1]]
# [1] "x"
# 
# $vars[[2]]
# [1] "y"
# 
# 
# $coefring
# [1] "QQ"
# 
# $order
# [1] "grevlex"
```

But a matrix of integers isn’t:

``` r
mat <- m2_matrix(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2))
str(mat)
#  'm2_matrix' num [1:3, 1:2] 1 2 3 4 5 6
#  - attr(*, "m2_name")= chr "m2rintmatrix00000005"
#  - attr(*, "m2_meta")=List of 1
#   ..$ ring: 'm2_polynomialring' logi NA
#   .. ..- attr(*, "m2_name")= chr "ZZ"
#   .. ..- attr(*, "m2_meta")=List of 3
#   .. .. ..$ vars    : NULL
#   .. .. ..$ coefring: chr "ZZ"
#   .. .. ..$ order   : chr "grevlex"
class(mat)
# [1] "m2_matrix" "m2"        "matrix"
m2_name(mat)
# [1] "m2rintmatrix00000005"
m2_meta(mat)
# $ring
# M2 Ring: ZZ[], grevlex order
```

Since a matrix of integers is an object in R, it’s represented as one,
and consequently we can compute with it directly as it if it were a
matrix; it is. On the other hand, since a ring is not, it’s an `NA`.
When dealing with M2, object like rings, that is to say objects without
R analogues, are more common than those like integer matrices.

Creating your own **m2r** wrapper
---------------------------------

We’ve already wrapped a number of Macaulay2 functions; for a list of
functions in **m2r**, check out `ls("package:m2r")`. But the list is
very far from exhaustive. To create your own wrapper function of a
Macaulay2 command, you’ll need to create an R file that looks like the
one below. This will create both value (e.g. `f`) and reference/pointer
(e.g. `f.`) versions of the function. As a good example of these at
work, see the scripts for
[`factor_n()`](https://github.com/musicman3320/m2r/blob/master/R/factor_n.R)
or
[`factor_poly()`](https://github.com/musicman3320/m2r/blob/master/R/factor_poly.R).

``` r
#' Function documentation header
#'
#' Function header explanation, can run several lines. Function
#' header explanation, can run several lines. Function header
#' explanation, can run several lines.
#'
#' @param esntl_parm_1 esntl_parm_1 description
#' @param esntl_parm_2 esntl_parm_2 description
#' @param code return only the M2 code? (default: \code{FALSE})
#' @param parse_parm_1 parse_parm_1 description
#' @param parse_parm_2 parse_parm_2 description
#' @param ... ...
#' @name f
#' @return (value version) parsed output or (reference/dot version)
#'   \code{m2_pointer}
#' @examples
#'
#' \dontrun{ requires Macaulay2 be installed
#'
#' # put examples here
#' 1 + 1
#'
#' }
#'





# value version of f (standard user version)
#' @rdname f
#' @export
f <- function(esntl_parm_1, esntl_parm_2, code = FALSE, parse_parm_1, parse_parm_2, ...) {

  # run m2
  args <- as.list(match.call())[-1]
  eargs <- lapply(args, eval, envir = parent.frame())
  pointer <- do.call(f., eargs)
  if(code) return(invisible(pointer))

  # parse output
  parsed_out <- m2_parse(pointer)

  # more parsing, like changing classes and such
  TRUE

  # return
  TRUE

}




# reference version of f (returns pointer to m2 object)
#' @rdname f
#' @export
f. <- function(esntl_parm_1, esntl_parm_2, code = FALSE, ...) {

  # basic arg checking
  TRUE

  # create essential parameters to pass to m2 this step regularizes input to m2, so it
  # is the one that deals with pointers, chars, rings, ideals, mpolyLists, etc.
  TRUE

  # construct m2_code from regularized essential parameters
  TRUE

  # message
  if(code) { message(m2_code); return(invisible(m2_code)) }

  # run m2 and return pointer
  m2.(m2_code)

}
```

Acknowledgements
----------------

This material is based upon work supported by the National Science
Foundation under Grant Nos.
[1321794](https://nsf.gov/awardsearch/showAward?AWD_ID=1321794) and
[1622449](https://nsf.gov/awardsearch/showAward?AWD_ID=1622449).

Installation
------------

Here’s how you can install the current *developmental* version of
**m2r**. Remember you need to have
[Macaulay2](http://www.math.uiuc.edu/Macaulay2/) downloaded; **m2r**
will look for it in your path variable (in the terminal, `echo $PATH`)
as set by `~/.bash_profile` or, if nonexistent, then `~/.bashrc`, then
`~/.profile`.

``` r
# install.packages("devtools")
devtools::install_github("coneill-math/m2r")
```
