##########################
# Power logistic model
##########################
# verifying the minimax D-optimality of a design with points x0 and weights w0
x0 <- c(-4.5515, 0.2130, 2.8075)
w0 <- c(0.4100, 0.3723, 0.2177)
# Power logistic model when s = .2
sensminimax(formula =  ~ (1/(1 + exp(-b * (x-a))))^.2,
            predvars = "x",
            parvars = c("a", "b"),
            family = binomial(),
            x = x0, w = w0,
            lx = -5, ux = 5,
            lp = c(0, 1), up = c(3, 1.5))

##############################
# A model with two predictors
##############################
# Verifying the minimax D-optimality of a design for a model with two predictors
# The model is the mixed inhibition model.
# X0 is the vector of four design points that are:
# (3.4614, 0) (4.2801, 3.1426) (30, 0) (30, 4.0373)
x0 <- c(3.4614, 4.2801, 30, 30, 0, 3.1426, 0, 4.0373)
w0 <- rep(1/4, 4)
sensminimax(formula =  ~ V*S/(Km * (1 + I/Kic)+ S * (1 + I/Kiu)),
            predvars = c("S", "I"),
            parvars = c("V", "Km", "Kic", "Kiu"),
            family = "gaussian",
            x = x0, w = w0,
            lx = c(0, 0), ux = c(30, 60),
            lp = c(1, 4, 2, 4), up = c(1, 5, 3, 5))

##########################################
# Standardized maximin D-optimal designs
##########################################
# Verifying the standardized maximin D-optimality of a design for
# the loglinear model
# First we should define the function for 'localdes' argument
# The function LDOD takes the parameters and returns the points and
# weights of the locally D-optimal design
LDOD <- function(theta0, theta1, theta2){
  ## param is the vector of theta = (theta0, theta1, theta2)
  lx <- 0 # lower bound of the design space
  ux <- 150 # upper bound of the design space
  param <- c()
  param[1] <- theta0
  param[2] <- theta1
  param[3] <- theta2
  xstar <- (ux+param[3]) * (lx + param[3]) *
    (log(ux + param[3]) - log(lx + param[3]))/(ux - lx) - param[3]
  return(list(x = c(lx, xstar, ux) , w = rep(1/3, 3)))
}
x0 <- c(0, 4.2494, 17.0324, 149.9090)
w0 <- c(0.3204, 0.1207, 0.2293, 0.3296)
\dontrun{
  sensminimax(formula = ~theta0 + theta1* log(x + theta2),
              predvars = c("x"),
              parvars = c("theta0", "theta1", "theta2"),
              x = x0, w = w0,
              lx = 0, ux = 150,
              lp = c(2, 2, 1), up = c(2, 2, 15),
              localdes = LDOD,
              standardized = TRUE,
              sens.minimax.control = list(n_seg = 10))
}
################################################################
# Not necessary!
# The rest of the examples here are only for professional uses.
################################################################
# Imagine you have written your own FIM, say in Rcpp that is faster than
# the FIM created by the formula interface here.

##########################
# Power logistic model
##########################
# For example, th cpp FIM function for the power logistic model is named:
FIM_power_logistic
args(FIM_power_logistic)
# The arguments do not match the standard of the argument 'fimfunc'
# in 'sensminimax'
# So we reparameterize it:
myfim1 <- function(x, w, param)
  FIM_power_logistic(x = x, w = w, param =param, s = .2)

args(myfim1)
\dontrun{
  # Verify minimax D-optimality of a design
  sensminimax(fimfunc = myfim1,
              x = c(-4.5515, 0.2130, 2.8075),
              w = c(0.4100, 0.3723, 0.2177),
              lx = -5, ux = 5,
              lp = c(0, 1), up = c(3, 1.5))
}
##############################
# A model with two predictors
##############################
# An example of a  model with two-predictors: mixed inhibition model
# Fisher information matrix:
FIM_mixed_inhibition
args(FIM_mixed_inhibition)

# We should first reparameterize the FIM to match the standard of the
# argument 'fimfunc'
myfim2 <- function(x, w, param){
  npoint <- length(x)/2
  S <- x[1:npoint]
  I <- x[(npoint+1):(npoint*2)]
  out <- FIM_mixed_inhibition(S = S, I = I, w = w, param = param)
  return(out)
}
args(myfim2)
\dontrun{
  # Verifyng minimax D-optimality of a design
  sensminimax(fimfunc = myfim2,
              x = c(3.4614, 4.2801, 30, 30, 0, 3.1426, 0, 4.0373),
              w = rep(1/4, 4),
              lx = c(0, 0), ux = c(30, 60),
              lp = c(1, 4, 2, 4), up = c(1, 5, 3, 5))
}

#########################################
# Standardized maximin D-optimal designs
#########################################
# An example of a user-written FIM function:
help(FIM_loglin)
# An example of verfying standardaized maximin D-optimality for a design
# Look how we re-define the function LDOD above
LDOD2 <- function(param){
  ## param is the vector of theta = (theta0, theta1, theta2)
  lx <- 0 # lower bound of the design space
  ux <- 150 # upper bound of the design space
  xstar <- (ux + param[3]) * (lx + param[3]) *
    (log(ux + param[3]) - log(lx + param[3]))/(ux - lx) - param[3]
  return(list(x = c(lx, xstar, ux) , w = rep(1/3, 3)))
}

args(LDOD2)

sensminimax(fimfunc = FIM_loglin,
            x = x0,
            w = w0,
            lx = 0, ux = 150,
            lp = c(2, 2, 1), up = c(2, 2, 15),
            localdes = LDOD2,
            standardized = TRUE)



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

sensminimax(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
            parvars = c("a", "b"), family = "binomial",
            lp = c(-2, 1), up = c(2, 1.5),
            crtfunc = Aopt,
            lx = -2, ux = 2,
            sensfunc = Aopt_sens,
            x = c(-2, .0033, 2), w = c(.274, .452, .274))

