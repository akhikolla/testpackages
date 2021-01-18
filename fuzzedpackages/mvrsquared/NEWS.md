# mvrsquared v0.1.1
This patches an error being thrown during testing on some Linux operating systems.
The root cause seems to be an imprecise calculation introduced in parallel computing.
See the note under `help(calc_rsquared)`.

# mvrsquared v0.1.0 
This version introduces parallel processing at the C++ level using RcppThread.

To calculate R-squared in parallel, set the `threads` argument to a number 
greater than 1 when calling `calc_rsquared`.

# mvrsquared v0.0.3
This version makes some changes to documentation including the README

# mvrsquared v0.0.2
This version includes

* An arXiv citation to the working paper deriving this method
* Changes to examples requested by CRAN

# mvrsquared v0.0.1
This version is the first release!

