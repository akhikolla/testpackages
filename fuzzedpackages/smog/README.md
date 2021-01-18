# smog
[![Travis-CI Build Status](https://travis-ci.com/chongma8903/smog.svg?branch=master)](https://travis-ci.com/github/chongma8903/smog)

Structural Modeling by using Overlapped Group Penalty

## Installation
```r
# Install smog from CRAN
install.packages("smog")

# or install the source type package from GitHub:
# install.packages("devtools")
devtools::install_github("chongma8903/smog")
```

## Features
* fits a linear non-penalized phenotype (demographic) variables and penalized groups of prognostic effect and predictive effect.
* satisfies such hierarchy structures that if a predictive effect exists, its prognostic effect must also exist.
* can deal with continuous, binomial or multinomial, and survival response variables.
* incorporates the iterative shrinkage-thresholding algorithm (ISTA) and the alternating direction method of multipliers algorithms (ADMM).

## Usage
Create a new S3 class of `smog`, and the kernal function `smog.default` (or `smog.formula`) returns an object of the S3 class `smog`. The kernel functions include:

* `smog.default`: input the data and parameters to yield a model of the class `smog`.
* `smog.formula`: can accept `formula` to fit the model for the data.
* `predict.smog`: produces the predicted response values for new data, provided a fitted model. 
* `cv.smog`: provides cross-validation analysis based on the data.  
* `cv.cglasso`: cross-validation for conditional group lasso approach  


## Examples
```r
sim = sim_rct_biomarker(n = 100, p = 20, p_prog = 2, p_pred = 2, p_both = 2)
y = sim$Y
x = sim$M
d = 20
g = c(d+1, rep(1:d,2))
v = c(rep(0,1), rep(1,2*d))
label = c("trt", rep(c("prog","pred"), c(d,d)))

sfit1 = cv.smog(x,y,g,v,label,family = "gaussian", type = "AIC")
plot(sfit1)

sfit2 = cv.cglasso(x,y,g,v,label,family = "gaussian", nlambda.max = 20)
plot(sfit2)

```
