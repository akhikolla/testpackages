#################################
# Exponential growth model
################################
# See how we set stopping rule by adjusting 'stop_rule', 'checkfreq' and 'stoptol'
# It calls the 'senslocally' function every checkfreq = 50 iterations to
# calculate the ELB. if ELB is greater than stoptol = .95, then the algoithm stops.

# initializing by one iteration
res1 <- locally(formula = ~a + exp(-b*x), predvars = "x", parvars = c("a", "b"),
                lx = 0, ux = 1, inipars = c(1, 10),
                iter = 1, k = 2,
                ICA.control= ICA.control(rseed = 100,
                                         stop_rule = "equivalence",
                                         checkfreq = 20, stoptol = .95))
\dontrun{
  # update the algorithm
  res1 <- update(res1, 150)
  #stops at iteration 21 because ELB is greater than .95
}

### fixed x, lx and ux are only required for equivalence theorem
\dontrun{
  res1.1 <- locally(formula = ~a + exp(-b*x), predvars = "x", parvars = c("a", "b"),
                    lx = 0, ux = 1, inipars = c(1, 10),
                    iter = 100,
                    x = c(.25, .5, .75),
                    ICA.control= ICA.control(rseed = 100))
  plot(res1.1)
  # we can not have an optimal design using this x
}

################################
## two parameter logistic model
################################
res2 <- locally(formula = ~1/(1 + exp(-b *(x - a))),
                predvars = "x", parvars = c("a", "b"),
                family = binomial(), lx = -3, ux = 3,
                inipars = c(1, 3), iter = 1, k = 2,
                ICA.control= list(rseed = 100, stop_rule = "equivalence",
                                  checkfreq = 50, stoptol = .95))
\dontrun{
  res2 <- update(res2, 100)
  # stops at iteration 51
}




################################
# A model with two predictors
################################
# mixed inhibition model
\dontrun{
  res3 <- locally(formula =  ~ V*S/(Km * (1 + I/Kic)+ S * (1 + I/Kiu)),
                  predvars = c("S", "I"),
                  parvars = c("V", "Km", "Kic", "Kiu"),
                  family = gaussian(),
                  lx = c(0, 0), ux = c(30, 60),
                  k = 4,
                  iter = 300,
                  inipars = c(1.5, 5.2, 3.4, 5.6),
                  ICA.control= list(rseed = 100, stop_rule = "equivalence",
                                    checkfreq = 50, stoptol = .95))
  # stops at iteration 100
}


\dontrun{
  # fixed x
  res3.1 <- locally(formula =  ~ V*S/(Km * (1 + I/Kic)+ S * (1 + I/Kiu)),
                    predvars = c("S", "I"),
                    parvars = c("V", "Km", "Kic", "Kiu"),
                    family = gaussian(),
                    lx = c(0, 0), ux = c(30, 60),
                    iter = 100,
                    x = c(20, 4, 20, 4, 10,  0, 0, 30, 3, 2),
                    inipars = c(1.5, 5.2, 3.4, 5.6),
                    ICA.control= list(rseed = 100))
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

res4 <- locally(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
                parvars = c("a", "b"), family = "binomial",
                lx = -3, ux = 3, inipars = c(1, 1.25),
                iter = 1, k = 2,
                crtfunc = Aopt,
                sensfunc = Aopt_sens,
                ICA.control = list(checkfreq = Inf))
\dontrun{
  res4 <- update(res4, 50)
}

# When the FIM of the model is defined directly via the argument 'fimfunc'
# the criterion function must have argument x, w fimfunc and param.
# use 'fimfunc' as a function of the design points x,  design weights w
# and param whenever needed.
Aopt2 <-function(x, w, param, fimfunc){
  sum(diag(solve(fimfunc(x = x, w = w, param = param))))
}
## the sensitivtiy function
# xi_x is a design that put all its mass on x in the definition of the sensitivity function
# x is a vector of design points
Aopt_sens2 <- function(xi_x, x, w, param, fimfunc){
  fim <- fimfunc(x = x, w = w, param = param)
  M_inv <- solve(fim)
  M_x <- fimfunc(x = xi_x, w = 1, param = param)
  sum(diag(M_inv %*% M_x %*%  M_inv)) - sum(diag(M_inv))
}

res4.1 <- locally(fimfunc = FIM_logistic,
                  lx = -3, ux = 3, inipars = c(1, 1.25),
                  iter = 1, k = 2,
                  crtfunc = Aopt2,
                  sensfunc = Aopt_sens2,
                  ICA.control = list(checkfreq = Inf))
\dontrun{
  res4.1 <- update(res4.1, 50)
  plot(res4.1)
}


# locally c-optimal design
# example from Chaloner and Larntz (1989) Figure 3
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


res4.2 <- locally(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
                  parvars = c("a", "b"), family = "binomial",
                  lx = -1, ux = 1, inipars = c(0, 7),
                  iter = 1, k = 2,
                  crtfunc = c_opt, sensfunc = c_sens,
                  ICA.control = list(rseed = 1, checkfreq = Inf))
\dontrun{
res4.2 <- update(res4.2, 100)
}

