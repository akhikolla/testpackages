# bayesDP 1.3.3
* New package maintainer (Graeme L. Hickey) since package was orphaned
* Updates to README, DESCRIPTION, NAMESPACE
* Added `stop` break to `discount_logit` for `method = mc`

# bayesDP 1.3.2
## Bug fixes and minor improvements
* Minor `bdplm` vignette typo fixes

# bayesDP 1.3.1
## Major new features
* Changes to inputs for `bdpsurvival`
  + Current and (optional) historical data are specified in separate data frames
* Updated normal approximation used for `method = "mc"` of the `bdpbinomial` and `bdpnormal` functions

## Bug fixes and minor improvements
* Summary method for `bdplm` now exists and mimics `lm`
* Removed `bdpbinomial` vignette language around success (vs event) 
* Reported one-arm sample size for `bdpsurvival` print method adjusted to current data only

# bayesDP 1.3.0
## Major new features
* Addition of the bdplm function for two-arm trials
* Users can now choose between 3 discount functions via the discount_function input:
  + Weibull CDF
  + Scaled Weibull CDF - scales the Weibull CDF so that the max possible value is 1
  + Identity - sets the discount weight to the posterior probability
* Removal of bdpregression

## Bug fixes and minor improvements
* Removed two-sided and one-sided function inputs to avoid confusion
* Posterior probabilities for `method = "mc"` switched from `pshisq` to pnorm 
* Updated vignettes to reflect new features

# bayesDP 1.2.0
## Major new features
* Supports one-arm regression analysis
* Two additional modular functions
* Implementation of Monte Carlo-based estimation of alpha discount

## Bug fixes and minor improvements
* Fixes to class slots
* Added `print` input to plot method

# bayesDP 1.1.0
## Major new features
* Supports two-arm survival analysis via hazard rate comparisons
* Completely revamped summary and print methods to produce better formatted results
* Plot method allows users to specify a `type`
* Added vignettes for each of `bdpbinomial`, `bdpnormal`, and `bdpsurvival`
* Implemented the `fix_alpha` input which allows users to set the historical data weight at `alpha_max`

## Bug fixes and minor improvements
* Fixed error with two-arm analysis where models did not fit if either the current or historical control data were not input
* Changed `two_side` input to logical
* Consolidated several internal functions into a single function for computational efficiency gains

# bayesDP 1.0.3
* README update
* Added plot types
* Added Vignettes
* Added logo
* Improved documentation
* Updated `print`, `summary`, `plot` methods
* Refactored `bdpnormal` / `bdpbinomial`

# bayesDP 1.0.2
* Crucial bugfixes

# bayesDP 1.0.1
* User requested bugfixes

# bayesDP 1.0.0
* Initial CRAN release with normal, binomial and survival functions
