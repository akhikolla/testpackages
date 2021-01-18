# Finding a robust design for the two-parameter logistic model
# See how we set a stopping rule.
# The ELB is computed every checkfreq = 30 iterations
# The optimization stops when the ELB is larger than stoptol = .95
res1 <- robust(formula = ~1/(1 + exp(-b *(x - a))),
               predvars = c("x"), parvars = c("a", "b"),
               family = binomial(),
               lx = -5, ux = 5, prob = rep(1/4, 4),
               parset = matrix(c(0.5, 1.5, 0.5, 1.5, 4.0, 4.0, 5.0, 5.0), 4, 2),
               iter = 1, k =3,
               ICA.control = list(stop_rule = "equivalence",
                                  stoptol = .95, checkfreq = 30))

\dontrun{
  res1 <- update(res1, 100)
  # stops at iteration 51
}


\dontrun{
  res1.1 <- robust(formula = ~1/(1 + exp(-b *(x - a))),
                   predvars = c("x"), parvars = c("a", "b"),
                   family = binomial(),
                   lx = -5, ux = 5, prob = rep(1/4, 4),
                   parset = matrix(c(0.5, 1.5, 0.5, 1.5, 4.0, 4.0, 5.0, 5.0), 4, 2),
                   x = c(-3, 0, 3),
                   iter = 150, k =3)
  plot(res1.1)
  # not optimal
}


###################################
# user-defined optimality criterion
##################################
# When the model is defined by the formula interface
# A-optimal design for the 2PL model.
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

res2 <- robust(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
               parvars = c("a", "b"), family = "binomial",
               lx = -3, ux = 3,
               iter = 1, k = 4,
               crtfunc = Aopt,
               sensfunc = Aopt_sens,
               prob = c(.25, .5, .25),
               parset = matrix(c(-2, 0, 2, 1.25, 1.25, 1.25), 3, 2),
               ICA.control = list(checkfreq = 50, stoptol = .999,
                                  stop_rule = "equivalence",
                                  rseed = 1))
\dontrun{
  res2 <- update(res2, 500)
}





# robust c-optimal design
# example from Chaloner and Larntz (1989), Figure 3, but robust design
c_opt <-function(x, w, a, b, fimfunc){
  gam <- log(.95/(1-.95))
  M <- fimfunc(x = x, w = w, a = a, b = b)
  c <- matrix(c(1, -gam * b^(-2)), nrow = 1)
  B <- t(c) %*% c
  sum(diag(B %*% solve(M)))
}

c_sens <- function(xi_x, x, w, a, b, fimfunc){
  gam <- log(.95/(1-.95))
  M <- fimfunc(x = x, w = w, a = a, b = b)
  M_inv <- solve(M)
  M_x <- fimfunc(x = xi_x, w = 1, a = a, b = b)
  c <- matrix(c(1, -gam * b^(-2)), nrow = 1)
  B <- t(c) %*% c
  sum(diag(B %*% M_inv %*% M_x %*%  M_inv)) - sum(diag(B %*% M_inv))
}


res3 <- robust(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
               parvars = c("a", "b"), family = "binomial",
               lx = -1, ux = 1,
               parset = matrix(c(0, 7, .2, 6.5), 2, 2, byrow = TRUE),
               prob = c(.5, .5),
               iter = 1, k = 3,
               crtfunc = c_opt, sensfunc = c_sens,
               ICA.control = list(rseed = 1, checkfreq = Inf))

\dontrun{
  res3 <- update(res3, 300)
}

