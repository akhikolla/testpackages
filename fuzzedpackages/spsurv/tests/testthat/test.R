
library(spsurv)
data("veteran")

# ML approach:
fit <- spbp(Surv(time, status)~karno+factor(celltype),
            approach = "mle",  data = veteran)
summary(fit)
vcov(fit)
coef(fit)
model.matrix(fit)
# stan_dens(fit)
# traceplot(fit)
# extract(fit)

# Bayesian approach:
fit2 <- spbp(Surv(time, status) ~ karno + factor(celltype),
             approach = "bayes",  data = veteran, chains = 1, iter = 1000, cores = 1)
summary(fit2)
# vcov(fit2)
# coef(fit2)
model.matrix(fit2)
stan_dens(fit2)
traceplot(fit2)
extract(fit2)

## PH model
bpph(Surv(time, status)~karno+factor(celltype),  data = veteran)## PO model
bpph(Surv(time, status)~karno+factor(celltype), approach = "bayes",  data = veteran,
     iter = 1000, chains = 1, cores = 1)## PH model

## PO model
bppo(Surv(time, status)~karno+factor(celltype),  data = veteran)## PO model
bppo(Surv(time, status)~karno+factor(celltype),  approach = "bayes",  data = veteran,
     iter = 1000, chains = 1, cores = 1)## PO model

## AFT model
bpaft(Surv(time, status)~karno+factor(celltype),  data = veteran)## AFT model
bpaft(Surv(time, status)~karno+factor(celltype),  approach = "bayes",  data = veteran,
     iter = 10, chains = 1, cores = 1)## PO model
