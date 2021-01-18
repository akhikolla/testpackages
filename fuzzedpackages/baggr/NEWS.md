# baggr 0.4.0 (February 2020)

### New features

* Covariates can now be used in all baggr() models: in "rubin" model they give meta-regression
  (group-level covariates), while in "full" and "logit" models they can be used for "regular"
  regression (individual-level covariates)
* Priors for covariates are set through the argument prior_beta
* You can work with regression coefficients for covariates 
    + you can access and summarise coefficients through fixed_effects(),
    + you will also see them when printing baggr objects; 
    + when using forest_plot() you can request `show = "covariates"`
* Prototype of pp_check() function now works for Rubin model (thanks to Brice Green)
  you can apply it to generate new datasets according to posterior distribution of treatment effect
  and contrast them with the observed quantities as part of model checking
* baggr_compare() function now has standard output which you can print() or plot(), 
  thanks to Brice Green
* Vignettes and documentation were updated to better describe binary data analysis
* We now give more warnings when plugging in stupid inputs.
  
### Bug fixes

* Messages for setting priors were accidentally given when e.g. running full pooling models
* All models were re-written to standardise our approach and syntax. 
  + "Full" model might now behave differently.
  + "Mutau" model will be re-written and generalised for next release.
  + LOO CV is also disabled for some models. Prompts will be given.

# baggr 0.3.0

### New features

* Binary data models for both summary and individual-level data.
* New vignette for working with binary data; see `vignette("baggr_binary")`.
* Expanded helper functions (esp. `prepare_ma`), esp. for prepping binary data.
* Added forest plot functionality for all types of models.
* Various outputs can now be transformed (main use case is `exp`, but any transform is allowed).
* Reworked vignette sections for pooling and cross-validation.
* Pooling statistics are now calculated for the whole model and better documented.
* More consistent theming, similar to bayesplot (thanks to Brice Green)
* Comparison of leave-one-out cross-validations with `loo_compare` (thanks to Brice Green)
  
### Bug fixes

* Re-enabled missing Cauchy priors



# baggr 0.2.0

### New features

* Users can now define priors in `baggr()` using a syntax similar to `rstanarm`.
  Extra priors are available
* `baggr()` outputs prior predictive distributions; they can be examined using
  `baggr_compare` and `effect_plot`, `effect_draw` -- 2 new functions
* More types of model comparisons are possible
* LOO CV function has been reworked
* Full pooling and no pooling models have been reworked to avoid divergent 
  transitions.



# baggr 0.1.0

First package version for CRAN. 
