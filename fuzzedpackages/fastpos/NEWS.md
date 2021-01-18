# fastpos 0.4.1

* tests now show more information if they fail, also the tolerance for an
individual deviation is higher to avoid false alarms
* progress bars are only shown if R is run interactively (this is required as
per the manual for writing R extensions)

# fastpos 0.4.0

* multiple cores are now supported via the future package, use parameter 
`n_cores`, gains are only visible if `n_studies` is high (overhead)
* tests for multiple cores were included
* `replace` parameter was added to the externally visible functions, and tests
for `replace = false`
* minor improvements in layout/documentation
* test added for old pop_create function
* test coverage is now somewhat higher since the untestable proportion of code
reduced

# fastpos 0.3.0

* Parameter for number of studies is now always n_studies. Some lower level
functions used number_of_studies.
* create_pop now creates a bivariate normal distribution with an exact rho, this
is slower than previously, but it is probably worth it.
* Problems with C++ progressbar (RcppProgress) were fixed. The code should also
be cleaner now. 
* Make C++ and Rstudio interruption similar: depending on the
timing, either Rstudio interrupts by itself and stops everything or C++ returns
a value. In the latter case, R now also stops quietly.
* Since the interrupts cannot be properly tested, code coverage is only 90%.

# fastpos 0.2.0

* If the corridor of stability is not reached, NA is returned by an internal
function that *simulate_pos* build upon. *n_not_breached* is simply the number of
NA values. The maximum sample size is still used for the calculation of the
quantiles if the corridor was not reached, which should be better than ignoring
the specific study altogether (i.e. treating it as an NA value). 

* *simulate_pos* now returns an IntegerVector (instead of NumericVector) to
better handle NA values.

* A simple test for relative precision was added. Another test was added that
unloads the package. Code coverage is now 100%.

* The vignette and readme were improved substantially. Word choice and grammar
were checked by an editor. A section on the speed of *fastpos* in comparison to
*corEvol* was added including a theoretical argument and an empirical test.

# fastpos 0.1.0

* First release