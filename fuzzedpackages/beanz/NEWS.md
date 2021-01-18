# beanz 2.0

* Added a `NEWS.md` file to track changes to the package.

* Added unit test
* Made it consistent throughout the software and the manuscript that the
  half-normal prior is for $\omega$ rather than $\omega^{2}$
* Added initial step-size as an option of the sampler
* Added adapt-delta as an option in the software and set the default value to 0.95
* Added explicit report of Rhat in the results and warnings for problematic
  convergence based on Rhat summaries
* Removed in shiny the option HMC vs. Fixed-params
* The Stan models are optimized by vectorization and non-centered parameterization
* Updated the GUI to allow the user to enter numerical values directly for priors
* Added looic as the measure of goodness of fit and the basis for model comparison 
* Renamed the parameters in lst.par.pri to match with the model instruction page 
* Added line for generating results for no subgroup effect model in relevant example code
* Added model information to the output from function r.rpt.tbl 


# beanz 2.2

* Minor fix in bzSummary() and bzSummaryComp(). Instead of returning matrix,
  they now return data frames.

