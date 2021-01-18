
p <- c(1, -2, 1, -1)
prior4.4 <- uniform(p -1.5, p + 1.5)
formula4.4 <- ~exp(b0+b1*x1+b2*x2+b3*x1*x2)/(1+exp(b0+b1*x1+b2*x2+b3*x1*x2))
prob4.4 <- ~1-1/(1+exp(b0 + b1 * x1 + b2 * x2 + b3 * x1 * x2))
predvars4.4 <-  c("x1", "x2")
parvars4.4 <- c("b0", "b1", "b2", "b3")
lb <- c(-1, -1)
ub <- c(1, 1)

## That is the optimal design when alpha = .25, see ?locallycomp on how to find it
xopt <- c(-1, -0.389, 1, 0.802, -1, 1, -1, 1)
wopt <- c(0.198, 0.618, 0.084, 0.1)

# We want to verfiy the optimality of the optimal design by the general equivalence theorem.

senslocallycomp(formula = formula4.4, predvars = predvars4.4, parvars = parvars4.4,
                family = binomial(), prob = prob4.4, lx = lb, ux = ub,
                alpha = .25, inipars = p, x = xopt, w = wopt)

\dontrun{
# is this design also optimal when alpha = .3

senslocallycomp(formula = formula4.4, predvars = predvars4.4, parvars = parvars4.4,
                family = binomial(), prob = prob4.4, lx = lb, ux = ub,
                alpha = .3, inipars = p, x = xopt, w = wopt)

# when alpha = .3
senslocallycomp(formula = formula4.4, predvars = predvars4.4, parvars = parvars4.4,
                family = binomial(), prob = prob4.4, lx = lb, ux = ub,
                alpha = .5, inipars = p, x = xopt, w = wopt)
# when alpha = .8
senslocallycomp(formula = formula4.4, predvars = predvars4.4, parvars = parvars4.4,
                family = binomial(), prob = prob4.4, lx = lb, ux = ub,
                alpha = .8, inipars = p, x = xopt, w = wopt)


# when alpha = .9
senslocallycomp(formula = formula4.4, predvars = predvars4.4, parvars = parvars4.4,
                family = binomial(), prob = prob4.4, lx = lb, ux = ub,
                alpha = .9, inipars = p, x = xopt, w = wopt)

## As can be seen, the design looses efficiency as alpha increases.
}
