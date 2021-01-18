## ----  message=FALSE, warning=FALSE-------------------------------------------
library(OptCirClust)

O = rnorm(100)
K = 5
Circumference = 6

## ----  message=FALSE, warning=FALSE-------------------------------------------
# Our recommended method is the fast and optimal FOCC:
result_FOCC <- CirClust(O, K, Circumference, method = "FOCC")

# The slow and optimal BOCC:
result_BOCC <- CirClust(O, K, Circumference, method = "BOCC")

# The slow and heuristic HEUC:
result_HEUC <- CirClust(O, K, Circumference, method = "HEUC")

## ----  message=FALSE, warning=FALSE-------------------------------------------
opar <- par(mar=c(0,0,2,0))

plot(result_FOCC, main = "FOCC: fast and optimal\n***Recommended***")

plot(result_BOCC, main = "BOCC: quadratic time\nalways optimal")

plot(result_HEUC, main = "HEUC: heuristic\nnot always optimal")

par(opar)

