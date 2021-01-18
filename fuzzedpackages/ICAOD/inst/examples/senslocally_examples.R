############################
# Exponential growth model
############################
# Verifying optimailty of a locally D-optimal design
senslocally(formula = ~a + exp(-b*x),
            predvars = "x", parvars = c("a", "b"),
            x = c(.1, 1), w = c(.5, .5),
            lx = 0, ux = 1, inipars = c(1, 10))


##############################
# A model with two predictors
##############################
x0 <- c(30, 3.861406, 30, 4.600633, 0, 0, 5.111376, 4.168798)
w0 <- rep(.25, 4)
senslocally(formula =  ~ V*S/(Km * (1 + I/Kic)+ S * (1 + I/Kiu)),
            predvars = c("S", "I"),
            parvars = c("V", "Km", "Kic", "Kiu"),
            x = x0, w = w0,
            lx = c(0, 0), ux = c(30, 60),
            inipars = c(1.5, 5.2, 3.4, 5.6))
\dontrun{
  # using package rgl for 3d plot:
  res<- senslocally(formula =  ~ V*S/(Km * (1 + I/Kic)+ S * (1 + I/Kiu)),
                    predvars = c("S", "I"),
                    parvars = c("V", "Km", "Kic", "Kiu"),
                    x = x0, w = w0,
                    lx = c(0, 0), ux = c(30, 60),
                    inipars = c(1.5, 5.2, 3.4, 5.6),
                    plot_3d = "rgl")

}

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

senslocally(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
            parvars = c("a", "b"), family = "binomial",
            inipars = c(0, 1.5),
            crtfunc = Aopt,
            lx = -2, ux = 2,
            sensfunc = Aopt_sens,
            x = c(-1,  1), w = c(.5, .5))
# not optimal
