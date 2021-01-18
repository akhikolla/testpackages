## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  echo = TRUE,
  message = FALSE,
  eval = TRUE,
  comment = "#>",
  fig = TRUE,
  fig.width = 7,
  fig.height = 5,
  fig.align = 'center'
)

## ----toy-example-data, results='asis'-----------------------------------------
library(knitr)
X1 <- c(4, 2, 2, 14, 6, 9, 4, 0, 1)
X2 <- c(0, 0, 1, 3, 2, 1, 2, 2, 2)
N1 <- rep(148, 9)
N2 <- rep(132, 9)
Y1 <- N1 - X1
Y2 <- N2 - X2
df <- data.frame(X1, Y1, X2, Y2)
kable(df, caption = "Toy Example")

## ----toy-example-5------------------------------------------------------------
library("DiscreteFDR")
DBH.sd.fast <- fast.Discrete(df, alternative = "two.sided", direction = "sd")
print(DBH.sd.fast)

summary(DBH.sd.fast)

## ----toy-example-summary------------------------------------------------------
DBH.sd.fast.summary <- summary(DBH.sd.fast)
DBH.sd.fast.summary$Table

## ----tow-example-2------------------------------------------------------------
p <- fisher.pvalues.support(df, alternative = "two.sided")
raw.pvalues <- p$raw

# or:
raw.pvalues2 <- DBH.sd.fast$Data$raw.pvalues

all(raw.pvalues == raw.pvalues2)

p.adjust(raw.pvalues, method = "BH")

## ----toy-example-crit---------------------------------------------------------
p <- fisher.pvalues.support(df, alternative = "two.sided")
raw.pvalues <- p$raw
pCDFlist <- p$support

DBH.sd.crit <- DBH(raw.pvalues, pCDFlist, 0.05, "sd", TRUE)
crit.vals.BH.disc <- DBH.sd.crit$Critical.values
crit.vals.BH.cont <- 1:9 * 0.05/9
cbind(sort(raw.pvalues), crit.vals.BH.disc, crit.vals.BH.cont)

## ----toy-example-plot-1-------------------------------------------------------
plot(DBH.sd.crit, col = c("red", "blue", "green"), pch = c(4, 2, 19), lwd = 2, type.crit = 'o',
     legend = "topleft", cex = 1.3)

## ----toy-example-plot-2-------------------------------------------------------
plot(DBH.sd.crit, col = c("red", "blue", "green"), pch = c(4, 2, 19), lwd = 2, type.crit = 'o',
     cex = 1.3, ylim = c(0, 0.25), main = "Comparison of discrete and continuous BH procedures")
points(crit.vals.BH.cont, pch = 19, cex = 1.3, lwd = 2)
legend("topright", c("Rejected", "Accepted", "Critical Values (disc.)", "Critical Values (cont.)"),
       col = c("red", "blue", "green", "black"), pch = c(4, 2, 19, 19), lwd = 2, lty = 0)

## ----toy-example-3------------------------------------------------------------
p$support[c(1,5)]

## ----toy-example-4------------------------------------------------------------
pCDFlist <- p$support
stepf <- lapply(pCDFlist, function(x) stepfun(x, c(0, x)))
par(mfcol = c(1, 3), mai = c(1, 0.5, 0.3, 0.1))
plot(stepf[[1]], xlim = c(0, 1), ylim = c(0, 1), do.points = FALSE, lwd = 1, lty = 1, ylab = "F(x)", 
     main = "(a)")
for(i in (2:9)){
  plot(stepf[[i]], add = TRUE, do.points = FALSE, lwd = 1, col = i)
}
segments(0, 0, 1, 1, col = "grey", lty = 2)

#   Plot xi
support <- sort(unique(unlist(pCDFlist)))
components <- lapply(stepf, function(s){s(support) / (1 - s(support))}) 
xi.values <- 1/9 * Reduce('+', components)
xi <- stepfun(support, c(0, xi.values))
plot(xi, xlim = c(0, 0.10), ylim = c(0, 0.10), do.points = FALSE, ylab = expression(xi), main = "(b)")
segments(0, 0, 0.1, 0.1, col = "grey", lty = 2)

#   Plot discrete critical values as well a BH constants and raw p-values
DBH.sd <- DBH(raw.pvalues, pCDFlist, direction = "sd", ret.crit.consts = TRUE)
plot(DBH.sd, col = c("black", "black", "red"), pch = c(4, 4, 19), type.crit = 'p', ylim = c(0, 0.15),
     cex = 1.3, main = "(c)", ylab = "Critical Values")
points(1:9, 0.05 * (1:9) / 9, col = "green", pch = 19, cex = 1.3)

mtext("Figure 1", 1, outer = TRUE, line = -2)

## ----format-am----------------------------------------------------------------
data(amnesia)
amnesia.formatted <- fisher.pvalues.support(amnesia[, 2:3], input = "HG2011")
raw.pvalues <- amnesia.formatted$raw
pCDFlist <- amnesia.formatted$support

## ----DFDR-Pharmacovigilance---------------------------------------------------
DBH.su  <-  DBH(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)
DBH.sd  <-  DBH(raw.pvalues, pCDFlist, direction = "sd", ret.crit.consts = TRUE)
ADBH.su <- ADBH(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)
ADBH.sd <- ADBH(raw.pvalues, pCDFlist, direction = "sd", ret.crit.consts = TRUE)
DBR.su  <-  DBR(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)

## ----Hist-Pharmaco------------------------------------------------------------
hist(DBH.sd)

## ----Plot-Pharmaco------------------------------------------------------------
m <- length(raw.pvalues)
crit.values.BH <- 0.05 * (1:m) / m
scale.points <- 0.7

plot(DBH.su, col = c("black", "black", "orange"), pch = NA, type.crit = 'p', xlim = c(1, 100),
     ylim = c(0, DBH.su$Critical.values[100]), ylab = "critical values", cex = scale.points, main = "")

points(crit.values.BH[1:105],          col = "green",  pch = 19, cex = scale.points)
points(DBH.sd$Critical.values[1:105],  col = "red",    pch = 19, cex = scale.points)
points(ADBH.su$Critical.values[1:105], col = "blue",   pch = 19, cex = scale.points)
points(ADBH.sd$Critical.values[1:105], col = "purple", pch = 19, cex = scale.points)
points(DBR.su$Critical.values[1:105],  col = "yellow", pch = 19, cex = scale.points)
points(sort(raw.pvalues),                              pch = 4,  cex = scale.points)
mtext("Figure 2", 1, outer = TRUE, line = -1)

## ----Reject-Pharmaco----------------------------------------------------------
rej.BH <- length(which(p.adjust(raw.pvalues, method = "BH") <= 0.05))
rej.DBH.su <- length(DBH.su$Indices)
rej.DBH.sd <- length(DBH.sd$Indices)
rej.ADBH.su <- length(ADBH.su$Indices)
rej.ADBH.sd <- length(ADBH.sd$Indices)
rej.DBR.su <- length(DBR.su$Indices)
c(rej.BH, rej.DBH.su, rej.DBH.sd, rej.ADBH.su, rej.ADBH.sd, rej.DBR.su)

## ----Poisson-Setup------------------------------------------------------------
lambda.0 <- c(0.6, 1.2, 0.7, 1.3, 1.0, 0.2, 0.8, 1.3, 0.9)
lambda.vector <- lambda.0
observations <- c(3, 3, 1, 2, 3, 3, 1, 2, 4)
configuration <- cbind(observations, lambda.0)
alpha <- 0.05
m <- length(observations)
print(configuration)

## ----Poisson-RawPValues-------------------------------------------------------
raw.pvalues <- sapply(1:m,function(i){ppois(observations[i] - 1, lambda.vector[i], lower.tail = FALSE)})
print(raw.pvalues)

## ----Poisson-tauMin-----------------------------------------------------------
y.min <- alpha/m * (1 + alpha/m)^(-1)
n.max <- sapply(1:m, function(w){qpois(y.min,        lambda.vector[w], lower.tail = FALSE)}) + 1
t.min <- sapply(1:m, function(w){ppois(n.max[w] - 1, lambda.vector[w], lower.tail = FALSE)})
s.min <- min(t.min)
print(s.min)

## ----Poisson-Supports---------------------------------------------------------
supports <- lapply(1:m, function(w){sort(ppois(0:n.max[w] - 1, lambda.vector[w], lower.tail = FALSE))})
DBH.sd <- DBH(raw.pvalues, supports, direction = "sd", ret.crit.consts = TRUE)

## ----Poisson-BH---------------------------------------------------------------
p.adjust(raw.pvalues, method = "BH")

## ----Poisson-DBH--------------------------------------------------------------
DBH.sd$Adjusted

## ----Poisson-Print------------------------------------------------------------
print(DBH.sd)

## ----Poisson-Rejections-------------------------------------------------------
stepf <- lapply(supports, function(x) stepfun(x, c(0, x)))
par(mfcol = c(1, 3), mai = c(1, 0.5, 0.3, 0.1))

plot(stepf[[1]], xlim = c(0,1), ylim = c(0,1), do.points = FALSE, lwd = 1, lty = 1, ylab = "F(x)", 
     main = "(a)")
for(i in (2:9)){
  plot(stepf[[i]], add = TRUE, do.points = FALSE, lwd = 1, col = i)
}
segments(0, 0, 1, 1, col = "grey", lty = 2)

#   Plot xi
support <- sort(unique(unlist(supports)))
components <- lapply(stepf, function(s){s(support) / (1 - s(support))}) 
xi.values <- 1/9 * Reduce('+', components)
xi <- stepfun(support, c(0, xi.values))
plot(xi, xlim = c(0, 0.10), ylim = c(0, 0.10), do.points = FALSE, ylab = expression(xi), main = "(b)")
segments(0, 0, 0.1, 0.1, col = "grey", lty = 2)

#   Plot discrete critical values as well a BH constants
DBH.sd <- DBH(raw.pvalues, supports, direction = "sd", ret.crit.consts = TRUE)
plot(DBH.sd, col = c("black", "black", "red"), pch = c(4, 4, 19), type.crit = 'p', ylim = c(0, 0.15),
     cex = 1.3, main = "(c)", ylab = "Critical Values")
points(1:9, 0.05 * (1:9) / 9, col = "green", pch = 19, cex = 1.3)

mtext("Figure 3", 1, outer = TRUE, line = -2)

