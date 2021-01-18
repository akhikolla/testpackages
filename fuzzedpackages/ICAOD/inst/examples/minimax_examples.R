########################################
# Two-parameter exponential growth model
########################################
res1 <- minimax (formula = ~a + exp(-b*x), predvars = "x", parvars = c("a", "b"),
                 lx = 0, ux = 1, lp = c(1, 1), up = c(1, 10),
                 iter = 1, k = 4,
                 ICA.control= ICA.control(rseed = 100),
                 crt.minimax.control = list(optslist = list(maxeval = 100)))
# The optimal design has 3 points, but we set k = 4 for illustration purpose to
#   show how the algorithm modifies the design by adjusting the weights
# The value of maxeval is changed to reduce the CPU time
\dontrun{
  res1 <- update(res1, 150)
  # iterating the algorithm up to 150 more iterations
}

res1 # print method
plot(res1) # Veryfying the general equivalence theorem

\dontrun{
  ## fixed x
  res1.1 <- minimax (formula = ~a + exp(-b*x), predvars = "x", parvars = c("a", "b"),
                     lx = 0, ux = 1, lp = c(1, 1), up = c(1, 10),
                     x = c(0, .5, 1),
                     iter = 150, k = 3, ICA.control= ICA.control(rseed = 100))
  # not optimal
}

########################################
# Two-parameter logistic model.
########################################
# A little playing with the tuning parameters
# The value of maxeval is reduced to 200 to increase the speed
cont1 <- crt.minimax.control(optslist = list(maxeval = 200))
cont2 <- ICA.control(rseed = 100, checkfreq = Inf, ncount = 60)

\dontrun{
  res2 <- minimax (formula = ~1/(1 + exp(-b *(x - a))), predvars = "x",
                   parvars = c("a", "b"),
                   family = binomial(), lx = -3, ux = 3,
                   lp = c(0, 1), up = c(1, 2.5), iter = 200, k = 3,
                   ICA.control= cont2, crt.minimax.control = cont1)
  print(res2)
  plot(res2)
}

############################################
# An example of a model with two predictors
############################################
# Mixed inhibition model
lower <- c(1, 4, 2, 4)
upper <- c(1, 5, 3, 5)
cont <- crt.minimax.control(optslist = list(maxeval = 100)) # to be faster
\dontrun{
  res3 <- minimax(formula =  ~ V*S/(Km * (1 + I/Kic)+ S * (1 + I/Kiu)),
                  predvars = c("S", "I"),
                  parvars = c("V", "Km", "Kic", "Kiu"),
                  lx = c(0, 0), ux = c(30, 60), k = 4,
                  iter = 100, lp = lower, up = upper,
                  ICA.control= list(rseed = 100),
                  crt.minimax.control = cont)

  res3 <- update(res3, 100)
  print(res3)
  plot(res3) # sensitivity plot
  res3$arg$time
}

# Now consider grid points instead of assuming continuous parameter space
# set n.grid to 5
\dontrun{
  res4 <- minimax(formula =  ~ V*S/(Km * (1 + I/Kic)+ S * (1 + I/Kiu)),
                  predvars = c("S", "I"),
                  parvars = c("V", "Km", "Kic", "Kiu"),
                  lx = c(0, 0), ux = c(30, 60),
                  k = 4, iter = 130, n.grid = 5, lp = lower, up = upper,
                  ICA.control= list(rseed = 100, checkfreq = Inf),
                  crt.minimax.control = cont)
  print(res4)
  plot(res4) # sensitivity plot
}

############################################
# Standardized maximin D-optimal designs
############################################
# Assume the purpose is finding STANDARDIZED designs
# We know from literature that the locally D-optimal design (LDOD)
# for this model has an analytical solution.
# The follwoing function takes the parameter as input and returns
# the design points and weights of LDOD.
# x and w are exactly similar to the arguments of 'fimfunc'.
# x is a vector and returns the design points 'dimension-wise'.
# see explanation of the arguments of 'fimfunc' in 'Details'.

LDOD <- function(V, Km, Kic, Kiu){
  #first dimention is for S and the second one is for I.
  S_min <- 0
  S_max <- 30
  I_min <- 0
  I_max <- 60
  s2 <- max(S_min, S_max*Km*Kiu*(Kic+I_min)/
              (S_max*Kic*I_min+S_max*Kic*Kiu+2*Km*Kiu*I_min+2*Km*Kiu*Kic))
  i3 <- min((2*S_max*Kic*I_min + S_max*Kic*Kiu+2*Km*Kiu*I_min+Km*Kiu*Kic)/
              (Km*Kiu+S_max*Kic), I_max)
  i4 <- min(I_min + (sqrt((Kic+I_min)*(Km*Kic*Kiu+Km*Kiu*I_min+
                                         S_max*Kic*Kiu+S_max*Kic*I_min)/
                            (Km*Kiu+S_max*Kic))), I_max )
  s4 <- max(-Km*Kiu*(Kic+2*I_min-i4)/(Kic*(Kiu+2*I_min-i4)), S_min)
  x <- c(S_max, s2, S_max, s4, I_min, I_min, i3, i4)
  return(list(x = x, w =rep(1/4, 4)))

}
formalArgs(LDOD)
\dontrun{
  minimax(formula =  ~ V*S/(Km * (1 + I/Kic)+ S * (1 + I/Kiu)),
          predvars = c("S", "I"),
          parvars = c("V", "Km", "Kic", "Kiu"),
          lx = c(0, 0), ux = c(30, 60),
          k = 4, iter = 300,
          lp = lower, up = upper,
          ICA.control= list(rseed = 100, checkfreq = Inf),
          crt.minimax.control = cont,
          standardized = TRUE,
          localdes = LDOD)
}


################################################################
# Not necessary!
# The rest of the examples here are only for professional uses.
################################################################
# Imagine you have written your own FIM, say in Rcpp that is faster than
# the FIM created by the formula interface above.

###########################################
# An example of a model with two predictors
###########################################
# For example, th cpp FIM function for the mixed inhibition model is named:
formalArgs(FIM_mixed_inhibition)

# We should reparamterize the arguments to match the standard of the
# argument 'fimfunc' (see 'Details').
myfim <- function(x, w, param){
  npoint <- length(x)/2
  S <- x[1:npoint]
  I <- x[(npoint+1):(npoint*2)]
  out <- FIM_mixed_inhibition(S = S, I = I, w = w, param = param)
  return(out)
}
formalArgs(myfim)

# Finds minimax optimal design, exactly as before, but NOT using the
# formula interface.
\dontrun{
  res5 <- minimax(fimfunc = myfim,
                  lx = c(0, 0), ux = c(30, 60), k = 4,
                  iter = 100, lp = lower, up = upper,
                  ICA.control= list(rseed = 100),
                  crt.minimax.control = cont)
  print(res5)
  plot(res5) # sensitivity plot
}
#########################################
# Standardized maximin D-optimal designs
#########################################
# To match the argument 'localdes' when no formula inteface is used,
# we should reparameterize LDOD.
# The input must be 'param' same as the argument of 'fimfunc'
LDOD2 <- function(param)
  LDOD(V = param[1], Km = param[2], Kic = param[3], Kiu = param[4])

# compare these two:
formalArgs(LDOD)
formalArgs(LDOD2)
\dontrun{
  res6 <- minimax(fimfunc = myfim,
                  lx = c(0, 0), ux = c(30, 60), k = 4,
                  iter = 300, lp = lower, up = upper,
                  ICA.control= list(rseed = 100, checkfreq = Inf),
                  crt.minimax.control = cont,
                  standardized = TRUE,
                  localdes = LDOD2)
  res6
  plot(res6)
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
\dontrun{
res7 <- minimax(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
                parvars = c("a", "b"), family = "binomial",
                lx = -2, ux = 2,
                lp = c(-2, 1), up = c(2, 1.5),
                iter = 400, k = 3,
                crtfunc = Aopt,
                sensfunc = Aopt_sens,
                crt.minimax.control = list(optslist = list(maxeval = 200)),
                ICA.control = list(rseed = 1))
  plot(res7)
}
# with grid points
res7.1 <- minimax(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
                  parvars = c("a", "b"), family = "binomial",
                  lx = -2, ux = 2,
                  lp = c(-2, 1), up = c(2, 1.5),
                  iter = 1, k = 3,
                  crtfunc = Aopt,
                  sensfunc = Aopt_sens,
                  n.grid = 9,
                  ICA.control = list(rseed = 1))
\dontrun{
  res7.1 <- update(res7.1, 400)
  plot(res7.1)
}

# When the FIM of the model is defined directly via the argument 'fimfunc'
# the criterion function must have argument x, w fimfunc and param.
# use 'fimfunc' as a function of the design points x,  design weights w and
#  the 'parvars' parameters whenever needed.
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
\dontrun{
res7.2 <- minimax(fimfunc = FIM_logistic,
                  lx = -2, ux = 2,
                  lp = c(-2, 1), up = c(2, 1.5),
                  iter = 1, k = 3,
                  crtfunc = Aopt2,
                  sensfunc = Aopt_sens2,
                  crt.minimax.control = list(optslist = list(maxeval = 200)),
                  ICA.control = list(rseed = 1))
  res7.2 <- update(res7.2, 200)
  plot(res7.2)
}
# with grid points
res7.3 <- minimax(fimfunc = FIM_logistic,
                  lx = -2, ux = 2,
                  lp = c(-2, 1), up = c(2, 1.5),
                  iter = 1, k = 3,
                  crtfunc = Aopt2,
                  sensfunc = Aopt_sens2,
                  n.grid = 9,
                  ICA.control = list(rseed = 1))
\dontrun{
  res7.3 <- update(res7.2, 200)
  plot(res7.3)
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

\dontrun{
res8 <- minimax(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
                parvars = c("a", "b"), family = "binomial",
                lx = -1, ux = 1,
                lp = c(-.3, 6), up = c(.3, 8),
                iter = 500, k = 3,
                crtfunc = c_opt, sensfunc = c_sens,
                ICA.control = list(rseed = 1, ncount = 100),
                n.grid = 12)
  plot(res8)
}



