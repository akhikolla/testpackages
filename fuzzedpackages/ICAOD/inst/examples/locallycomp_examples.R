## Here we produce the results of Table 2 in in McGree and Eccleston (2008)
# For D- and P-efficiency see, ?leff and ?peff

p <- c(1, -2, 1, -1)
prior4.4 <- uniform(p -1.5, p + 1.5)
formula4.4 <- ~exp(b0+b1*x1+b2*x2+b3*x1*x2)/(1+exp(b0+b1*x1+b2*x2+b3*x1*x2))
prob4.4 <- ~1-1/(1+exp(b0 + b1 * x1 + b2 * x2 + b3 * x1 * x2))
predvars4.4 <-  c("x1", "x2")
parvars4.4 <- c("b0", "b1", "b2", "b3")
lb <- c(-1, -1)
ub <- c(1, 1)


# set checkfreq = Inf to ask for equivalence theorem at final step.
res.0 <- locallycomp(formula = formula4.4, predvars = predvars4.4, parvars = parvars4.4,
                     family = binomial(), prob = prob4.4, lx = lb, ux = ub,
                     alpha = 0, k = 1, inipars = p, iter = 10,
                     ICA.control = ICA.control(checkfreq = Inf))

\dontrun{
res.25 <- locallycomp(formula = formula4.4, predvars = predvars4.4, parvars = parvars4.4,
                      family = binomial(), prob = prob4.4, lx = lb, ux = ub,
                      alpha = .25, k = 4, inipars = p, iter = 350,
                      ICA.control = ICA.control(checkfreq = Inf))

res.5 <- locallycomp(formula = formula4.4, predvars = predvars4.4, parvars = parvars4.4,
                     family = binomial(), prob = prob4.4, lx = lb, ux = ub,
                     alpha = .5, k = 4, inipars = p, iter = 350,
                     ICA.control = ICA.control(checkfreq = Inf))
res.75 <- locallycomp(formula = formula4.4, predvars = predvars4.4, parvars = parvars4.4,
                      family = binomial(), prob = prob4.4, lx = lb, ux = ub,
                      alpha = .75, k = 4, inipars = p, iter = 350,
                      ICA.control = ICA.control(checkfreq = Inf))

res.1 <- locallycomp(formula = formula4.4, predvars = predvars4.4, parvars = parvars4.4,
                     family = binomial(), prob = prob4.4, lx = lb, ux = ub,
                     alpha = 1, k = 4, inipars = p, iter = 350,
                     ICA.control = ICA.control(checkfreq = Inf))

#### computing the D-efficiency
# locally D-optimal design is locally DP-optimal design when alpha = 1.

leff(formula = formula4.4, predvars = predvars4.4, parvars = parvars4.4, family = binomial(),
     x1 = res.0$evol[[10]]$x, w1 = res.0$evol[[10]]$w,
     inipars = p,
     x2 = res.1$evol[[350]]$x, w2 = res.1$evol[[350]]$w)

leff(formula = formula4.4, predvars = predvars4.4, parvars = parvars4.4, family = binomial(),
     x1 = res.25$evol[[350]]$x, w1 = res.25$evol[[350]]$w,
     inipars = p,
     x2 = res.1$evol[[350]]$x, w2 = res.1$evol[[350]]$w)

leff(formula = formula4.4, predvars = predvars4.4, parvars = parvars4.4, family = binomial(),
     x1 = res.5$evol[[350]]$x, w1 = res.5$evol[[350]]$w,
     inipars = p,
     x2 = res.1$evol[[350]]$x, w2 = res.1$evol[[350]]$w)


leff(formula = formula4.4, predvars = predvars4.4, parvars = parvars4.4, family = binomial(),
     x1 = res.75$evol[[350]]$x, w1 = res.75$evol[[350]]$w,
     inipars = p,
     x2 = res.1$evol[[350]]$x, w2 = res.1$evol[[350]]$w)



#### computing the P-efficiency
# locally p-optimal design is locally DP-optimal design when alpha = 0.

leff(formula = formula4.4, predvars = predvars4.4, parvars = parvars4.4, family = binomial(),
     x2 = res.0$evol[[10]]$x, w2 = res.0$evol[[10]]$w,
     prob = prob4.4,
     type = "PA",
     inipars = p,
     x1 = res.25$evol[[350]]$x, w1 = res.25$evol[[350]]$w)

leff(formula = formula4.4, predvars = predvars4.4, parvars = parvars4.4, family = binomial(),
     x2 = res.0$evol[[10]]$x, w2 = res.0$evol[[10]]$w,
     prob = prob4.4,
     inipars = p,
     type = "PA",
     x1 = res.5$evol[[350]]$x, w1 = res.5$evol[[350]]$w)

leff(formula = formula4.4, predvars = predvars4.4, parvars = parvars4.4, family = binomial(),
     x2 = res.0$evol[[10]]$x, w2 = res.0$evol[[10]]$w,
     prob = prob4.4,
     inipars = p,
     type = "PA",
     x1 = res.75$evol[[350]]$x, w1 = res.75$evol[[350]]$w)


leff(formula = formula4.4, predvars = predvars4.4, parvars = parvars4.4, family = binomial(),
     x2 = res.0$evol[[10]]$x, w2 = res.1$evol[[10]]$w,
     prob = prob4.4,
     type = "PA",
     inipars = p,
     x1 = res.1$evol[[350]]$x, w1 = res.1$evol[[350]]$w)


}
