# spsurv
An R package for semi-parametric survival analysis.

The *spsurv* package was designed to contribute with a flexible set of semi-parametric survival regression modelings, including proportional hazards (PH), proportional odds (PO), and accelerated failure time (AFT) models for right-censored data.

- Install and load the *spsurv* package using the devtools package.

```r
install.packages("devtools")

library(devtools)
install_github("rvpanaro/spsurv")
```
- Check out the main fitter function examples.

```r
library("KMsurv")
data("larynx")

library(spsurv)

## Maximum Likelihood
fit <- spbp(Surv(time, delta)~age+factor(stage),
                    approach = "mle",  data = larynx)
summary(fit)      

## NUTS sampling (Bayesian)
fit2 <- spbp(Surv(time, delta)~age+factor(stage),
                     approach = "bayes",  data = larynx,
                     iter = 2000, chains = 1, warmup = 1000)
summary(fit2)
```

The *spsurv* already provides:
- Integration with Stan software.
- Estimates either in Bayesian or Frequentist (point estimate) inferential approaches.
- Three survival regression classes: PH, PO and AFT.
- Six distinct prior specifications in a Bayesian analysis.
