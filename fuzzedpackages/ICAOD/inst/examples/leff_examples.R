p <- c(1, -2, 1, -1)
prior4.4 <- uniform(p -1.5, p + 1.5)
formula4.4 <- ~exp(b0+b1*x1+b2*x2+b3*x1*x2)/(1+exp(b0+b1*x1+b2*x2+b3*x1*x2))
prob4.4 <- ~1-1/(1+exp(b0 + b1 * x1 + b2 * x2 + b3 * x1 * x2))
predvars4.4 <-  c("x1", "x2")
parvars4.4 <- c("b0", "b1", "b2", "b3")


# Locally D-optimal design is as follows:
## weight and point of D-optimal design
# Point1     Point2     Point3     Point4
# /1.00000 \ /-1.00000\ /0.06801 \ /1.00000 \
# \-1.00000/ \-1.00000/ \1.00000 / \1.00000 /
#   Weight1    Weight2    Weight3    Weight4
# 0.250      0.250      0.250      0.250

xopt_D <- c(1, -1, .0680, 1, -1, -1, 1, 1)
wopt_D <- rep(.25, 4)

# Let see if we use only three of the design points, what is the relative efficiency.
leff(formula = formula4.4, predvars = predvars4.4, parvars = parvars4.4, family = binomial(),
     x1 = c(1, -1, .0680,  -1, -1, 1), w1 = c(.33, .33, .33),
     inipars = p,
     x2 = xopt_D, w2 = wopt_D)
# Wow, it heavily drops!


# Locally P-optimal design has only one support point and is -1 and 1
xopt_P <- c(-1, 1)
wopt_P <- 1

# What is the relative P-efficiency of the D-optimal design with respect to P-optimal design?
leff(formula = formula4.4, predvars = predvars4.4, parvars = parvars4.4, family = binomial(),
     x1 = xopt_D, w1 = wopt_D,
     inipars = p,
     type = "PA",
     prob = prob4.4,
     x2 = xopt_P, w2 = wopt_P)
# .535


