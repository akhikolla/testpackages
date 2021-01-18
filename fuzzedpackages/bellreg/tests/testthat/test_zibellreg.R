
library(bellreg)

# ML approach:
mle <- zibellreg(cells ~ smoker+gender|smoker+gender, data = cells, approach = "mle")
summary(mle)
coef(mle)
vcov(mle)
confint(mle)


# Bayesian approach:
bayes <- zibellreg(cells ~ 1|smoker+gender, data = cells, approach = "bayes", chains = 2, iter = 100)
summary(bayes)


