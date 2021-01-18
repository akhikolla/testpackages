library(Countr)
library(dplyr)
library(xtable)

data(fertility)
form <- children ~ german + years_school + voc_train + university + religion +
                   year_birth + rural + age_marriage
gam <- renewalCount(formula = form, data = fertility, dist = "gamma",
                    computeHessian = TRUE, 
                    control = renewal.control(trace = 0, method = "nlminb")
                    )
v1 <- gam$vcov
v2 <- vcov(gam)
all(v1 == v2)

se <- se.coef(gam, type = "asymptotic")
se

ci <- confint(gam, type = "asymptotic")
ci

gam <- addBootSampleObject(gam, R = 400, parallel = "multicore", ncpus = 14)

gam$vcov <- matrix()
varCovar <- vcov(gam, method = "boot")

capboot <- "Bootstrap variance-covariance matrix of model \texttt{gam}."
print(xtable(varCovar, digits = -1, caption = capboot, label = "tab:varCovar"), 
      rotate.colnames = TRUE, floating.environment = "sidewaystable" )

se_boot  <- se.coef(gam, type = "boot")
se_boot

ci_boot  <- confint(gam, level = 0.95, type = "boot", bootType = "norm")
ci_boot

save.image()
