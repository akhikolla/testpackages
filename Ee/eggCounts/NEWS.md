## eggCounts 2.3

* No longer requiring integer correction factors


## eggCounts 2.2, 2.2-1:

* Updated makevars files

* Fixed a bug in the unpaired model (thanks to Paul)

* Improved documentation

* Improved examples

* Improved printing of warning messages


## eggCounts 2.1, 2.1-1, 2.1-2:

* Improved usability of fecr_stanExtra() 

* Fixed a minor bug in stan2mcmc() for paired model

* Updated citation information, URL

* Updated Stan version requirement


## eggCounts 2.0:

* Added individual efficacy model for paired design

* Added vignette

* Added functions for selecting priors in small sample size (<10) settings

* Added simplified model for small sample sizes in paired design

* Added fecr_stanExtra() to import external models hosted on gitHub

* Improved function documentations

* Modified default arguments: nsamples, nburnin, nchain, zeroInflation

* Fixed parameter sampling issue in using the unpaired model without zero inflation

* Added simple MCMC convergence diagnostic output based on split chain potential scale reduction factors

* Added possibility to not save all model parameters in fec_stan() and fecr_stan()

* Added possibility to simulate with different efficacy in simData2s()

* Re-named function fecr_probability() to fecr_probs()

* Added 'plot' argument in fecr_probs() to plot posterior density of reduction 

* Added warnings in using fecrtCI() when all of post-treatment counts are zero

* Added warnings in using fec_stan() and fecr_stan() when divergent transitions occur

* Cleaned Stan models code

* Removed un-used dataset


## eggCounts 1.4: 

* Modified default arguments for fec_stan and fecr_stan

* Optimized package onLoad messages

* Added plotCounts function

* Corrected fec_stan help file examples (thanks to Ivailo)

* Corrected coding error with non-default prior for ziunpaired model (thanks to Anja)


## eggCounts 1.3:

* Fixed bug in printing MCMC summaries for meanEPG

* Fixed implementation bug of the zero-inflation models that was in v1.2 

* Fixed minor bug in printing warnings

* Added function to compute probability of reduction based on posterior density

* Added flexibility to un-constrain reduction between 0 and 1

* Updated citation information 


## eggCounts 1.2:

* Allowed direct extraction of summary output

* Fixed minor bug in validating correction factors

* Fixed minor bug in converting multiple mcmc chains

* Updated citation information


## eggCounts 1.1, 1.1-1, 1.1-2:  

* Updated documentation

* Removed redundant functions 

* Updated LinkingTo version requirements

* Updated Stan program syntax 

* Fixed minor bug in the mcmc output of the 1-sample zero-inflation model


## eggCounts 1.0:  

* MCMC inference migrated from Metropolis-Hastlings algorithm and Gibbs sampler to Stan modelling language.

* Further improved bug fix of version 0.4.

* Improved models.

* Full-support on zero-inflation models.


## eggCounts 0.4, 0.4-1:  

* Bug fixes introduced since version 0.1.

* Further improved bug fix of version 0.3.

* Alerts the user in case of atypical data through 'NOTE'.

* Added verboselevel functionality for detailed debugging. 

* Improved help and added tests.


## eggCounts 0.2, 0.3:  

* stable release on CRAN

* Updated package structure: NEWS, ChangeLog and CITATION files.

* Bug fix: all "zero" draws from rgamma are captured and
    set to `minrgamma=.Machine$double.eps`.
    Thus numerical results are likely to change.
