# SurvBoost
SurvBoost

High dimensional variable selection method for stratified proportional hazards model. Implementing an extension of gradient boosting methods applied to survival data, incorporating stratification to relax the proportional hazards assumption in many cases.


### Installation

In order to install the package, several other R packages must be installed. 
The code relies on Rcpp, RcppArmadillo, and RcppParallel in order to improve computational speed. 
Additionally the survival package is used for simulation and post selection inference and will be required for installation. 

The following line of R code installs the package from CRAN: 
```
install.packages("SurvBoost")

# Or you can also install the package from github: 
install.packages("devtools")
devtools::install_github("EmilyLMorris/survBoost")
```

### Example 

Simple example of variable selection:
```{r, eval = FALSE}
# Using fixed number of iterations: 
boosting_core(Surv(time,delta) ~ strata(strata_idx) + V1 + V2 + V3 + V4 + V5, 
              data, rate=0.1, control=500) 

# Specifying the number of variables to select: 
boosting_core(formula, data, rate=0.1, control_method="num_selected", control_parameter=5)
```

Check whether it is appropriate to stratify with a certain variable: 
```{r, eval = FALSE}
strata.boosting(data$strata_variable, data$time)
```

Plot the coefficient paths: 
```{r, eval = FALSE}
plot.boosting(boosting_core.output, type="coefficients")
```

### Methods description

[SurvBoost: An R Package for High-Dimensional Variable Selection in the Stratified Proportional Hazards Model via Gradient Boosting](https://arxiv.org/abs/1803.07715)

[![Travis build status](https://travis-ci.org/EmilyLMorris/survBoost.svg?branch=master)](https://travis-ci.org/EmilyLMorris/survBoost)