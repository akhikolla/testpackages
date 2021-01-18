## ---- warning=FALSE, message=FALSE, eval=FALSE---------------------------
#  install.packages("motmot")

## ---- warning=FALSE, message=FALSE---------------------------------------
library(motmot)

## ------------------------------------------------------------------------
data(anolis.tree)
data(anolis.data)
attach(anolis.data)
anolis.tree

## ------------------------------------------------------------------------
sortedData <- sortTraitData(phy = anolis.tree, y = anolis.data, 
  data.name = "Male_SVL", pass.ultrametric = TRUE)
phy <- sortedData$phy
male.length <- sortedData$trait

## ----plot1, fig.cap = "Figure 1. TraitData showing the realtive male snout-vent length at the tips", echo = T, fig.height = 5, fig.width = 5,dpi=200----
traitData.plot(y = male.length, phy = phy, lwd.traits = 2, 
  col.label = "#00008050", tck = -0.01, mgp = c(0, 0.2, 0), 
  cex.axis = 0.5, show.tips = FALSE)

## ------------------------------------------------------------------------
## uncomment to view the tree
# plot(phy, show.tip.label=FALSE, no.margin=TRUE, edge.col="grey20")
# nodelabels(182, 182, bg="black", col="white")
phy.clade <- extract.clade(phy, 182)
temp.mat <- male.length[match(phy.clade$tip.label, rownames(male.length)), ]
male.length.clade <- as.matrix(temp.mat)

## ------------------------------------------------------------------------
bm.ml <- transformPhylo.ML(phy = phy.clade, y = male.length.clade, model = "bm")
bm.ml

## ------------------------------------------------------------------------
lambda.ml <- transformPhylo.ML(phy = phy.clade, y = male.length.clade, 
  model = "lambda")
lambda.ml

## ----plot2, fig.cap = "Figure 2. Profile plot of ML estimation for Pagel's lambda", echo = T, fig.height = 5, fig.width = 5,dpi=200----
lambda.ml <- transformPhylo.ML(phy = phy.clade, y = male.length.clade, 
  model = "lambda", profilePlot = TRUE)

## ------------------------------------------------------------------------
p.value <- 1 - pchisq(lambda.ml$MaximumLikelihood - bm.ml$logLikelihood, 1)
p.value

## ------------------------------------------------------------------------
bm.ml$AICc - lambda.ml$AICc

## ------------------------------------------------------------------------
delta.ml <- transformPhylo.ML(ph = phy.clade, y = male.length.clade,
  model = "delta")
delta.ml

## ----plot3, fig.cap = "Figure 3. Comparison of BM and Kappa transformed trees.", echo = T, fig.height = 5, fig.width = 5, ,dpi=200----
kappa.ml <- transformPhylo.ML(phy = phy.clade, y = male.length.clade, 
  model = "kappa", profilePlot = FALSE, returnPhy = TRUE)
par(mfrow = c(1, 2))
plot.phylo(phy.clade, show.tip.label = FALSE, no.margin = TRUE)
mtext(text = "Original phylogeny", side = 3, cex = 0.7, line = -1)
plot.phylo(kappa.ml$kappaPhy, show.tip.label = FALSE, no.margin = TRUE)
mtext(text = "Kappa model phylogeny", side = 3, cex = 0.7, line = -1)
mtext(text = "Kappa = 1e-8", side = 3, cex = 0.7, line = -2)

## ----plot4, fig.cap = "Figure 4. Profile plot to estimate alpha", echo = T, fig.height = 5, fig.width = 5, ,dpi=200----
ou.ml <- transformPhylo.ML(phy = phy.clade, y = male.length.clade,
  model = "OU", profilePlot = TRUE, upperBound = 2)
ou.ml

## ------------------------------------------------------------------------
p.value <- 1 - pchisq(ou.ml$MaximumLikelihood - bm.ml$logLikelihood, 1)
p.value
bm.ml$AICc - ou.ml$AICc

## ----plot5, fig.cap = "Figure 5. Profile plot to estimate the ACDC parameter", echo = T, fig.height = 5, fig.width = 5,dpi=200----
acdc.ml <- transformPhylo.ML(phy = phy.clade, y = male.length.clade, 
  model = "ACDC", profilePlot = TRUE)
acdc.ml

## ------------------------------------------------------------------------
p.value.2 <- 1 - pchisq(acdc.ml$MaximumLikelihood - bm.ml$logLikelihood , 1)
p.value.2

## ------------------------------------------------------------------------
transformPhylo.ML(phy = phy.clade, y = male.length.clade, model = "ACDC", 
  profilePlot = FALSE, upperBound = -1e-6, print.warning = FALSE)

## ----plot6, fig.cap = "Figure 6. Profile plot to estimate the psi parameter", echo = TRUE, fig.height = 5, fig.width = 5,dpi=200----
psi.ml <- transformPhylo.ML(phy = phy.clade, y = male.length.clade, 
  model = "psi", profilePlot = TRUE)
psi.ml

## ------------------------------------------------------------------------
p.value.psi <- 1 - pchisq(psi.ml$MaximumLikelihood - bm.ml$logLikelihood , 1)
p.value.psi

## ------------------------------------------------------------------------
psi_ext.est <- transformPhylo.ML(phy = phy.clade, y = male.length.clade,
  model = "psi", profilePlot = FALSE, hiddenSpeciation = TRUE, full.phy = phy)
all.equal(psi.ml, psi_ext.est)

## ----plot7, fig.cap = "Figure 7. Two clades used in the multipsi model", echo = T, fig.height = 5, fig.width = 5,dpi=200----
plot(phy.clade, no.margin = TRUE, cex = 0.8)
two.clade.labels <- c(rep("a", 17), rep("b", 37))
edgelabels(two.clade.labels, col=c(rep("blue", 17), rep("red", 37)),
  bg = "white")

## ------------------------------------------------------------------------
transformPhylo.ML(phy = phy.clade, y = male.length.clade, model = "multipsi", 
  branchLabels = c(rep("a", 17), rep("b", 37)), hiddenSpeciation = TRUE,   
  full.phy = phy)

## ------------------------------------------------------------------------
acdc.ml.lambda <- transformPhylo.ML(phy = phy.clade, y = male.length.clade,
  model = "ACDC", lambdaEst = TRUE)
# original ACDC model
acdc.ml
# ACDC model plus lambda
acdc.ml.lambda

## ------------------------------------------------------------------------
# p value of the ACDC and ACDC+lambda models. No significant improvement
1 - pchisq(acdc.ml.lambda$MaximumLikelihood - acdc.ml$MaximumLikelihood , df = 1)
# p value of the BM and ACDC+lambda model comparison. No significant improvement
1 - pchisq(acdc.ml.lambda$MaximumLikelihood - bm.ml$logLikelihood, df = 2)

## ----plot8, fig.cap = "Figure 8. Lineages with different rates of evolution", echo = T, fig.height = 5, fig.width = 5, ,dpi=200----
plot(phy.clade, show.tip.label = FALSE, no.margin = TRUE,
  edge.col = "grey20")
nodelabels(c(32, 49), c(32, 49), bg = "black", col = "white")

## ------------------------------------------------------------------------
cladeRate.ml <- transformPhylo.ML(phy = phy.clade, y = male.length.clade, 
  model = "clade", nodeIDs = c(32, 49))
cladeRate.ml

## ------------------------------------------------------------------------
# tm1 algorithm not run
# tm1.ml <- transformPhylo.ML(y = male.length.clade, phy = phy.clade, 
#   model = "tm1", minCladeSize = 2, nSplits = 3)
# trait.medusa.tm1.summary <- summary(tm1.ml, cutoff = 2, AICc = TRUE)
# tm2 model
tm2.ml <- transformPhylo.ML(y = male.length.clade, phy = phy.clade,
  model = "tm2", minCladeSize = 5, nSplits = 2)

## ----plot9, fig.cap = "Figure 9. The subset of the tree showing the rate heterogeneity estimated from the traitMedusa model", echo = T, fig.height = 5, fig.width = 5,dpi=200----
trait.medusa.tm2.summary <- summary(tm2.ml, cutoff = 2, AICc = TRUE)
trait.medusa.tm2.summary
colour_motmot <- plot(x = trait.medusa.tm2.summary, reconType = "rates",
  type = "fan", cex=0.5, edge.width=2)

## ------------------------------------------------------------------------
## uncomment to run
# set.seed(203);
# calcCutOff(phy.clade, n = 1000, model = "tm2", minCladeSize = 5, nSplits = 1);
##      95% 
## 5.698198 

## ------------------------------------------------------------------------
summary(tm2.ml, cutoff = 5.698198, AICc = TRUE)$Rates

## ------------------------------------------------------------------------
timeSlice.10.ml <- transformPhylo.ML(y = male.length.clade, phy = phy.clade, 
  model = "timeSlice", splitTime = 10)

## ----plot10, fig.cap = "Figure 10. TimeSlice plot with a split at 10 Ma", echo = T, fig.height = 5, fig.width = 5,dpi=200----
outputSummary <- plot(timeSlice.10.ml, cutoff = 0.001, cex = 0.55,
  edge.width = 2, cex.plot = 0.8, colour.ramp = c("blue", "red"),
  label.offset = 0.5)

## ------------------------------------------------------------------------
outputSummary$RatesCI

## ------------------------------------------------------------------------
timeSlice.ml <- transformPhylo.ML(y = male.length.clade, phy = phy.clade,
  model = "timeSlice", nSplits = 1, boundaryAge = 8)

## ----plot11, fig.cap = "Figure 11. TimeSlice plot with Maximum likelihood estimation of split time", echo = T, fig.height = 5, fig.width = 5,dpi=200----
outputSummary <- plot(timeSlice.ml, cutoff = 1, cex = 0.2, edge.width = 2,
  cex.plot = 0.8, colour.ramp = c("blue", "red"), label.offset = 0.5)

## ------------------------------------------------------------------------
modeSlice.ml <- transformPhylo.ML(y = male.length.clade, phy = phy.clade, 
  model = "modeSlice", splitTime = c(40, 30), mode.order = c("ACDC", "OU", "BM"),
  rate.var = TRUE, acdcScalar = TRUE)
modeSlice.ml$AICc
bm.ml$AICc

## ------------------------------------------------------------------------
bm.model <- transformPhylo.ML(male.length.clade, phy = phy.clade, model = "bm")
nested.acdc <- transformPhylo.ML(male.length.clade, phy = phy.clade, 
  model = "ACDC", nodeIDs = 44)
nested.ou <- transformPhylo.ML(male.length.clade, phy = phy.clade, model = "OU", 
  nodeIDs = 44)

1 - pchisq(nested.acdc$MaximumLikelihood - bm.model$logLikelihood, 1)
1 - pchisq(nested.ou$MaximumLikelihood - bm.model$logLikelihood, 1)

## ---- results="hide"-----------------------------------------------------
set.seed(12) # set seed so run will be identical - for example use only
lambda.mcmc <- transformPhylo.MCMC(y = male.length.clade, phy = phy.clade, 
  model = "lambda", mcmc.iteration = 2000, burn.in = 0.25, random.start = FALSE, 
  sample.every = 1)

## ------------------------------------------------------------------------
lambda.mcmc[1:4]

## ----plot12, fig.cap = "Figure 12. MCMC trace for Pagel's lambda", echo = T, fig.height = 5, fig.width = 5, dpi = 200----
mcmc.plot(lambda.mcmc)

## ------------------------------------------------------------------------
data(finches)
emp.tree <- finch.tree
emp.data <- finch.data
param.simulation <- chr.disp.param(emp.tree, n.sim = 100, max.sigma = 8, 
  max.a = 8, ntraits=1, mc.cores = 1)

## ------------------------------------------------------------------------
# Data and phylogeny
data(anolis.tree)
anolis.tree$node.label <- NULL
set.seed(3492)
lm.data <- transformPhylo.sim(phy = anolis.tree, n = 2, model = "bm")
dat <- data.frame(x = lm.data[, 1], y = lm.data[, 2], names = anolis.tree$tip, 
  row.names = anolis.tree$tip)
# pgls from CAPER with matrix inversion
library(caper)
comp.dat <- comparative.data(anolis.tree, dat, names)
time.now <- Sys.time()
matrix.inv.caper <- pgls(y ~ x, data = comp.dat, lambda = "ML")
pgls.time <- Sys.time() - time.now
pgls.time
time.now <- Sys.time()
picModel <- pic.pgls(formula = y ~ x, phy = anolis.tree, y = dat,
  lambda = "ML", return.intercept.stat = FALSE)
pic.time <- Sys.time() - time.now
pic.time

## ------------------------------------------------------------------------
# from caper
summary(matrix.inv.caper)
# from MOTMOT
picModel

