# Version 3.1.0 (2020-10-24)

Speed improvement.


# Version 3.0.0 (2019-10-26)

The package now uses `Rcpp` to evaluate the hypergeometric function in all cases, 
for real or complex values of the arguments.


# Version 2.0.0 (2019-10-09)

The package now uses `Rcpp` to evaluate the hypergeometric function when none 
of its arguments is a complex number.


# Version 1.0.1 (2019-09-26)

- Fixed `lmvgamma` for complex values

- Allows complex values `z` with `Re(z)<0` in `mvgamma`

- Added more unit tests


# Version 1.0.0 (2019-09-16)

First release.
