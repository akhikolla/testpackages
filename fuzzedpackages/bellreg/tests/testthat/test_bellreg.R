
library(bellreg)

# ML approach:
mle <- bellreg(nf ~ lroll, data = faults, approach = "mle")
summary(mle)
coef(mle)
vcov(mle)
confint(mle)

# Bayesian approach:
bayes <- bellreg(nf ~ lroll, data = faults, approach = "bayes")
summary(bayes)
