# Verifying a robust design for the two-parameter logistic model
sensrobust(formula = ~1/(1 + exp(-b *(x - a))),
           predvars = c("x"),
           parvars = c("a", "b"),
           family = binomial(),
           prob = rep(1/4, 4),
           parset = matrix(c(0.5, 1.5, 0.5, 1.5, 4.0, 4.0, 5.0, 5.0), 4, 2),
           x = c(0.260, 1, 1.739), w = c(0.275, 0.449, 0.275),
           lx = -5, ux = 5)


###################################
# user-defined optimality criterion
##################################
# When the model is defined by the formula interface
# Checking the A-optimality  for the 2PL model.
# the criterion function must have argument x, w fimfunc and the parameters defined in 'parvars'.
# use 'fimfunc' as a function of the design points x,  design weights w and
#  the 'parvars' parameters whenever needed.
Aopt <-function(x, w, a, b, fimfunc){
  sum(diag(solve(fimfunc(x = x, w = w, a = a, b = b))))
}
## the sensitivtiy function
# xi_x is a design that put all its mass on x in the definition of the sensitivity function
# x is a vector of design points
Aopt_sens <- function(xi_x, x, w, a, b, fimfunc){
  fim <- fimfunc(x = x, w = w, a = a, b = b)
  M_inv <- solve(fim)
  M_x <- fimfunc(x = xi_x, w = 1, a  = a, b = b)
  sum(diag(M_inv %*% M_x %*%  M_inv)) - sum(diag(M_inv))
}

sensrobust(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
           parvars = c("a", "b"), family = "binomial",
           crtfunc = Aopt,
           sensfunc = Aopt_sens,
           lx = -3, ux = 3,
           prob = c(.25, .5, .25),
           parset = matrix(c(-2, 0, 2, 1.25, 1.25, 1.25), 3, 2),
           x = c(-2.469, 0, 2.469), w = c(.317, .365, .317))
# not optimal. the optimal design has four points. see the last example in ?robust
