
library(BuyseTest)
library(data.table)


## * Settings
n.patients <- c(60,65)
BuyseTest.options(check = TRUE,
                  keep.pairScore = TRUE,
                  method.inference = "none",
                  trace = 0)

## * Simulated data
set.seed(10)
dt.sim <- simBuyseTest(n.T = n.patients[1],
                       n.C = n.patients[2],
                       argsBin = list(p.T = c(0.5,0.75)),
                       argsCont = list(mu.T = 1:3, sigma.T = rep(1,3)),
                       argsTTE = list(rates.T = 1/(1:3), rates.Censoring.T = rep(1,3)))

## * Compare
## iFormula <- treatment ~ tte(eventtime1,status1,1) + tte(eventtime2,status2,1) + tte(eventtime1,status1,0.25) 
iFormula <- treatment ~ tte(eventtime1,status1,1) + tte(eventtime2,status2,1) + tte(eventtime1,status1,0.25) + tte(eventtime2,status2,0)

BuyseTest.options(engine = "GPC_cpp")
e.BT1 <- BuyseTest(iFormula, data = dt.sim,
                   method.inference = "u-statistic", scoring.rule = "Peron")

