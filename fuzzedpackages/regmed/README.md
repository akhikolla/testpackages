[![CRAM_Status_Badge](http://www.r-pkg.org/badges/version/regmed)](https://CRAN.R-project.org/package=regmed)
[![Downloads](http://cranlogs.r-pkg.org/badges/regmed)](https://CRAN.R-project.org/package=regmed)
[![Total-Downloads](https://cranlogs.r-pkg.org/badges/grand-total/regmed)](https://CRAN.R-project.org/package=regmed)

# The `regmed` Package
Mediation analysis for multiple mediators by penalized structural equation 
models using sparse group lasso. The penalty considers the natural groupings 
of parameters that determine mediation, as well as encourages sparseness of 
the model parameters. 

# The `regmed.grid()` Function

`regmed.grid()` is a function that fits regularized mediation models over a vector 
grid of lambda penalty values. Structural equation models for analysis of multiple mediators
are extended by creating a sparse group lasso penalized model such that
the penalty considers the natural groupings of the pair of parameters
that determine mediation, as well as encourages sparseness of the model
parameters. 

# The `regmed.fit()` Function

`regmed.fit()`  fits a regularized mediation model for a specified lambda penalty value. We provide summary and plot methods implemented on the S3 class created by the function.
