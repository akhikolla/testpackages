## ---- align='center', fig.width=7---------------------------------------------
library('TreeTools', quietly = TRUE, warn.conflicts = FALSE)

tree1 <- BalancedTree(10)
tree2 <- PectinateTree(10)

origPar <- par(mfrow = 1:2, mar = rep(0.2, 4)) # Set up plotting area
plot(tree1)
plot(tree2)

## -----------------------------------------------------------------------------
library('TBRDist')

## -----------------------------------------------------------------------------
USPRDist(tree1, tree2)

## -----------------------------------------------------------------------------
ReplugDist(tree1, tree2)

## -----------------------------------------------------------------------------
TBRDist(tree1, tree2, exact = TRUE)

## -----------------------------------------------------------------------------
TBRDist(tree1, tree2, exact = TRUE, maf = TRUE)

## -----------------------------------------------------------------------------
TBRDist(tree1, tree2, exact = FALSE)

## -----------------------------------------------------------------------------
TBRDist(tree1, list(tree1, tree2), exact = TRUE, maf = TRUE)

## -----------------------------------------------------------------------------
USPRDist(list(tree1, tree2, tree2), list(tree2, tree1, tree2))

## ----echo=FALSE---------------------------------------------------------------
par(origPar)

